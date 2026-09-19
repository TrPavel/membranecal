"""PRE-OUTCOME acquisition only: no local-accuracy endpoint is imported or computed."""

from __future__ import annotations

# Source metadata explanations remain uninterrupted strings.
# ruff: noqa: E501
import argparse
import concurrent.futures
import gzip
import hashlib
import json
import re
import threading
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import UTC, datetime
from pathlib import Path

import gemmi
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "workspace/v02"
V01 = ROOT / "data/frozen/membranecal_real_pilot_v0.1"
LOCK = threading.Lock()
SOURCE_LOCKS = {}


def read(p):
    return json.loads(Path(p).read_text(encoding="utf-8"))


def write(p, d):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(
        json.dumps(d, sort_keys=True, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


CODE_SHA256 = sha(__file__)


def canonical_identity(uniprot, forbidden):
    accession = uniprot["primaryAccession"]
    aliases = {accession, *uniprot.get("secondaryAccessions", [])}
    if aliases & forbidden:
        raise ValueError("V01_CANONICAL_OR_ALIAS_OVERLAP")
    return accession


def confidence_values(qa, metricdefs, length):
    ids = {
        metricdefs["id"][i]
        for i, n in enumerate(metricdefs.get("name", []))
        if n.lower() == "plddt" and metricdefs["mode"][i] == "local"
    }
    indices = [i for i, mid in enumerate(qa.get("metric_id", [])) if mid in ids]
    keys = [
        (qa["model_id"][i], qa["label_asym_id"][i], int(qa["label_seq_id"][i])) for i in indices
    ]
    if len(ids) != 1 or len(keys) != len(set(keys)) or len({(k[0], k[1]) for k in keys}) != 1:
        raise ValueError("AFDB_CONFIDENCE_AMBIGUOUS")
    result = {int(qa["label_seq_id"][i]): float(qa["metric_value"][i]) / 100 for i in indices}
    if set(result) != set(range(1, length + 1)) or any(not 0 <= v <= 1 for v in result.values()):
        raise ValueError("AFDB_CONFIDENCE_QC")
    return result


def get(kind, identifier, url, suffix):
    with LOCK:
        lock = SOURCE_LOCKS.setdefault((kind, identifier, suffix), threading.Lock())
    with lock:
        return _get(kind, identifier, url, suffix)


def _get(kind, identifier, url, suffix):
    path = WORK / "raw" / kind / (identifier + suffix)
    if not path.resolve().is_relative_to(WORK.resolve()):
        raise ValueError("source path escapes v0.2 workspace")
    meta = path.with_name(path.name + ".source.json")
    if path.exists() and meta.exists():
        rec = read(meta)
        if sha(path) != rec["sha256"]:
            raise ValueError("source checksum mismatch")
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(3):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "MembraneCal-research-pilot/0.2"}
            )
            with urllib.request.urlopen(req, timeout=45) as response:
                data = response.read()
            tmp = path.with_suffix(path.suffix + ".part")
            tmp.write_bytes(data)
            tmp.replace(path)
            write(
                meta,
                {
                    "schema_version": "0.2.0",
                    "source": kind,
                    "identifier": identifier,
                    "url": url,
                    "retrieved_utc": datetime.now(UTC).isoformat(),
                    "sha256": sha(path),
                    "bytes": len(data),
                    "raw_path": path.relative_to(ROOT).as_posix(),
                    "transformation": "unmodified HTTP response bytes",
                },
            )
            return path
        except urllib.error.HTTPError as e:
            if e.code in {400, 404, 410}:
                raise
            if attempt == 2:
                raise
        except (OSError, TimeoutError):
            if attempt == 2:
                raise
        time.sleep(1 + attempt)


def discover(limit):
    path = get(
        "opm_index",
        "primary_structures",
        "https://opm-back.cc.lehigh.edu/opm-backend/primary_structures?pageSize=10000",
        ".json",
    )
    payload = read(path)
    if len(payload["objects"]) != payload["total_objects"]:
        raise ValueError("incomplete OPM discovery; pagination needed")
    old = read(V01 / "proteins.json")["records"]
    oldentries = {r["entry_id"].lower() for r in old}
    buckets = defaultdict(list)
    excluded = Counter()
    for r in payload["objects"]:
        if r["pdbid"].lower() in oldentries:
            excluded["V01_STRUCTURE_OVERLAP"] += 1
            continue
        match = re.match(r"[0-9.]+", str(r.get("resolution") or ""))
        if r["type_id"] != 1 or r["classtype_id"] != 1 or not match:
            excluded["NOT_ALPHA_HELICAL_TRANSMEMBRANE_OR_NO_RESOLUTION"] += 1
            continue
        if float(match.group()) > 3.5 or not r.get("thickness"):
            excluded["RESOLUTION_OR_BOUNDARY"] += 1
            continue
        buckets[r["family_id"]].append(r)
    # One metadata representative per named UniProt-code signature; family-round-robin
    # prevents the candidate cap being consumed by thousands of related receptor structures.
    for family, rs in buckets.items():
        seen, unique = set(), []
        for r in sorted(
            rs,
            key=lambda x: (
                float(re.match(r"[0-9.]+", str(x["resolution"])).group()),
                hashlib.sha256(x["pdbid"].encode()).hexdigest(),
            ),
        ):
            signature = tuple(sorted(r.get("uniprotcodes") or [r["pdbid"]]))
            if signature not in seen:
                unique.append(r)
                seen.add(signature)
        buckets[family] = unique
    candidates = []
    while len(candidates) < limit:
        batch = [buckets[k].pop(0) for k in sorted(buckets) if buckets[k]]
        if not batch:
            break
        candidates.extend(batch[: limit - len(candidates)])
    write(
        WORK / "discovery.json",
        {
            "schema_version": "0.2.0",
            "index_sha256": sha(path),
            "total_opm_records": len(payload["objects"]),
            "prefilter_exclusions": dict(excluded),
            "policy": "alpha-helical transmembrane OPM; resolution <=3.5; source boundary; no old PDB; unique code signatures; per-family rounds, best resolution then SHA256(PDB ID)",
            "limit": limit,
            "candidates": candidates,
        },
    )
    print(
        json.dumps({"discovered": len(payload["objects"]), "selected": len(candidates)}), flush=True
    )


def ss_labels(block, asym):
    labels = {}
    for cat, label in [("_struct_conf.", "HELIX"), ("_struct_sheet_range.", "SHEET")]:
        d = block.get_mmcif_category(cat)
        for i, a in enumerate(d.get("beg_label_asym_id", [])):
            if a != asym or d["end_label_asym_id"][i] != asym:
                continue
            if label == "HELIX" and not str(d["conf_type_id"][i]).startswith("HELX"):
                continue
            try:
                lo, hi = int(d["beg_label_seq_id"][i]), int(d["end_label_seq_id"][i])
            except (ValueError, TypeError):
                continue
            for p in range(lo, hi + 1):
                labels[p] = label if p not in labels or labels[p] == label else "CONFLICT"
    return labels


def sifts_mapping(path, author_chain):
    root = ET.fromstring(gzip.decompress(path.read_bytes()))
    result = {}
    for residue in root.iter():
        if residue.tag.rsplit("}", 1)[-1] != "residue":
            continue
        cross = [e.attrib for e in residue if e.tag.rsplit("}", 1)[-1] == "crossRefDb"]
        pdbs = [
            x for x in cross if x.get("dbSource") == "PDB" and x.get("dbChainId") == author_chain
        ]
        ups = [x for x in cross if x.get("dbSource") == "UniProt"]
        if len(pdbs) == 1 and len(ups) == 1:
            key = pdbs[0].get("dbResNum")
            if key and key not in {"null", "None"}:
                val = (ups[0]["dbAccessionId"], int(ups[0]["dbResNum"]), ups[0].get("dbResName"))
                if key in result and result[key] != val:
                    result[key] = None
                else:
                    result[key] = val
    return result


def family_evidence(uniprot):
    phrases = []
    for c in uniprot.get("comments", []):
        if c["commentType"] == "SIMILARITY":
            phrases.extend(t["value"] for t in c.get("texts", []))
    pfam = sorted(
        r["id"] for r in uniprot.get("uniProtKBCrossReferences", []) if r["database"] == "Pfam"
    )
    return phrases, pfam


def prepare_entry(candidate):
    entry = candidate["pdbid"].lower()
    output = WORK / "prepared" / (entry + ".json")
    if output.exists() and read(output).get("acquisition_code_sha256") == CODE_SHA256:
        return read(output)["status"]
    try:
        cp = get("pdb", entry, f"https://files.rcsb.org/download/{entry}.cif.gz", ".cif.gz")
        sp = get(
            "sifts",
            entry,
            f"https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{entry}.xml.gz",
            ".xml.gz",
        )
        op = get("opm", entry, f"https://opm-assets.storage.googleapis.com/pdb/{entry}.pdb", ".pdb")
        block = gemmi.cif.read(str(cp)).sole_block()
        methods = block.get_mmcif_category("_exptl.").get("method", [])
        if not any(m in {"X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"} for m in methods):
            raise ValueError("EXPERIMENTAL_METHOD")
        resolutions = [
            float(v)
            for cat, key in [
                ("_refine.", "ls_d_res_high"),
                ("_em_3d_reconstruction.", "resolution"),
            ]
            for v in block.get_mmcif_category(cat).get(key, [])
            if v not in {None, False, "?", "."}
        ]
        if not resolutions or max(resolutions) > 3.5 or min(resolutions) <= 0:
            raise ValueError("EXPERIMENTAL_RESOLUTION_QC")
        refstructure = gemmi.read_structure(str(cp))
        refstructure.setup_entities()
        opstructure = gemmi.read_structure(str(op))
        # OPM legacy coordinates and numbering are used only for experimental-frame transfer.
        opca = {}
        for chain in opstructure[0]:
            for residue in chain:
                atoms = [a for a in residue if a.name == "CA"]
                if atoms:
                    a = max(atoms, key=lambda a: a.occ)
                    opca[(chain.name, str(residue.seqid))] = (
                        np.array([a.pos.x, a.pos.y, a.pos.z]),
                        residue.name,
                    )
        entity_for = {s: e.name for e in refstructure.entities for s in e.subchains}
        poly = block.get_mmcif_category("_entity_poly.")
        sequences = {
            e: re.sub(r"\s+", "", poly["pdbx_seq_one_letter_code_can"][i])
            for i, e in enumerate(poly.get("entity_id", []))
        }
        old = read(V01 / "proteins.json")["records"]
        oldacc = {r["accession"] for r in old}
        development = read(WORK / "development_sequences.json")["records"]
        oldacc.update(
            a
            for r in development
            for a in [r["accession"], r["current_primary_accession"], *r["secondary_accessions"]]
        )
        oldconstruct = {r["construct_sequence_sha256"] for r in old} | {
            hashlib.sha256(r["construct_sequence"].encode()).hexdigest() for r in old
        }
        available, rejected = [], []
        used_entities = set()
        for chain in refstructure[0]:
            # One label-asym per entity, chosen by deposited order before any prediction outcome.
            for sub in chain.subchains():
                asym = sub.subchain_id()
                ent = entity_for.get(asym)
                if not ent or ent in used_entities or ent not in sequences:
                    continue
                used_entities.add(ent)
                try:
                    seq = sequences[ent]
                    construct_hash = hashlib.sha256(seq.encode()).hexdigest()
                    if construct_hash in oldconstruct:
                        raise ValueError("V01_CONSTRUCT_OVERLAP")
                    mapping = sifts_mapping(sp, chain.name)
                    coords = []
                    for residue in sub:
                        ca = sorted(
                            [a for a in residue if a.name == "CA" and a.occ > 0],
                            key=lambda a: (-a.occ, str(a.altloc)),
                        )
                        if ca and residue.label_seq is not None:
                            a = ca[0]
                            coords.append(
                                {
                                    "label": int(residue.label_seq),
                                    "author": str(residue.seqid),
                                    "name": residue.name,
                                    "xyz": [a.pos.x, a.pos.y, a.pos.z],
                                    "altloc": str(a.altloc).replace("\x00", ""),
                                    "occupancy": float(a.occ),
                                }
                            )
                    hits = [
                        r
                        for r in coords
                        if (chain.name, r["author"]) in opca
                        and opca[(chain.name, r["author"])][1] == r["name"]
                    ]
                    if len(hits) < 30 or len(hits) / max(1, len(coords)) < 0.8:
                        raise ValueError("OPM_TARGET_CHAIN_COVERAGE")
                    moving = np.array([r["xyz"] for r in hits])
                    target = np.array([opca[(chain.name, r["author"])][0] for r in hits])
                    mc, tc = moving.mean(axis=0), target.mean(axis=0)
                    u, _, vt = np.linalg.svd((moving - mc).T @ (target - tc))
                    rot = u @ np.diag([1, 1, np.linalg.det(u @ vt)]) @ vt
                    rmsd = float(
                        np.sqrt(np.mean(np.sum(((moving - mc) @ rot + tc - target) ** 2, axis=1)))
                    )
                    if rmsd > 0.5:
                        raise ValueError("OPM_REFERENCE_TRANSFORM_RMSD")
                    annotated = [mapping.get(r["author"]) for r in coords]
                    accessions = {r[0] for r in annotated if r}
                    if len(accessions) != 1:
                        raise ValueError("SIFTS_MULTI_OR_NO_ACCESSION")
                    acc = next(iter(accessions))
                    if acc.split("-")[0] in oldacc:
                        raise ValueError("V01_CANONICAL_OVERLAP")
                    if "-" in acc:
                        raise ValueError("NONCANONICAL_ISOFORM")
                    up = get(
                        "uniprot", acc, f"https://rest.uniprot.org/uniprotkb/{acc}.json", ".json"
                    )
                    uni = read(up)
                    canonical = uni["sequence"]["value"]
                    sifts_acc = acc
                    acc = canonical_identity(uni, oldacc)
                    tm = [f for f in uni.get("features", []) if f["type"] == "Transmembrane"]
                    if not tm:
                        raise ValueError("NO_UNIPROT_TRANSMEMBRANE_SUPPORT")
                    if not 100 <= len(canonical) <= 1500 or len(seq) < 0.5 * len(canonical):
                        raise ValueError("SEQUENCE_LENGTH_OR_CONSTRUCT_COVERAGE")
                    mapped = [(r, mapping[r["author"]]) for r in coords if mapping.get(r["author"])]
                    if (
                        len(mapped) < 100
                        or len(mapped) / len(coords) < 0.95
                        or len(coords) / len(seq) < 0.7
                    ):
                        raise ValueError("MAPPING_OR_OBSERVABILITY_QC")
                    positions = [m[1] for _, m in mapped]
                    if (
                        len(positions) != len(set(positions))
                        or min(positions) < 1
                        or max(positions) > len(canonical)
                    ):
                        raise ValueError("NONUNIQUE_OR_OUT_OF_RANGE_MAPPING")
                    identities = [
                        gemmi.find_tabulated_residue(r["name"]).one_letter_code.upper()
                        == canonical[m[1] - 1]
                        for r, m in mapped
                    ]
                    if sum(identities) / len(identities) < 0.98:
                        raise ValueError("CONSTRUCT_IDENTITY_LT_98")
                    api = get(
                        "afdb_api",
                        acc,
                        f"https://alphafold.ebi.ac.uk/api/prediction/{acc}",
                        ".json",
                    )
                    models = [
                        m
                        for m in read(api)
                        if m.get("uniprotSequence") == canonical
                        and m.get("uniprotStart") == 1
                        and m.get("uniprotEnd") == len(canonical)
                    ]
                    if len(models) != 1:
                        raise ValueError("AFDB_CANONICAL_PRODUCT_AMBIGUOUS")
                    model = models[0]
                    pred = get("afdb_cif", acc, model["cifUrl"], ".cif")
                    pb = gemmi.cif.read(str(pred)).sole_block()
                    predicted_sequences = pb.get_mmcif_category("_entity_poly.").get(
                        "pdbx_seq_one_letter_code_can", []
                    )
                    if (
                        len(predicted_sequences) != 1
                        or re.sub(r"\s+", "", predicted_sequences[0]) != canonical
                    ):
                        raise ValueError("AFDB_CIF_SEQUENCE_MISMATCH")
                    # Firewall: only local metric IDs/values and sequence identifiers are read.
                    # Prediction Cartesian coordinates are not inspected here.
                    qa = pb.get_mmcif_category("_ma_qa_metric_local.")
                    metricdefs = pb.get_mmcif_category("_ma_qa_metric.")
                    confidence = confidence_values(qa, metricdefs, len(canonical))
                    ss = ss_labels(block, asym)
                    half = float(candidate["thickness"]) / 2
                    boundary = re.search(r"1/2 of bilayer thickness:\s*([0-9.]+)", op.read_text())
                    if not boundary or abs(float(boundary.group(1)) - half) > 0.1:
                        raise ValueError("OPM_BOUNDARY_DISAGREEMENT")
                    residue_rows = []
                    for (r, m), identity in zip(mapped, identities, strict=True):
                        depth = abs(float(((np.array(r["xyz"]) - mc) @ rot + tc)[2]))
                        ratio = depth / half
                        region = (
                            "TM_CORE"
                            if ratio <= 0.8
                            else "MEMBRANE_INTERFACE"
                            if ratio <= 1.2
                            else "STRUCTURED_EXTRAMEMBRANE"
                            if ratio > 1.6
                            else "UNKNOWN_OR_EXCLUDED"
                        )
                        secondary = ss.get(r["label"], "UNAVAILABLE")
                        if secondary not in {"HELIX", "SHEET"}:
                            region = "UNKNOWN_OR_EXCLUDED"
                        residue_rows.append(
                            {
                                "canonical_position": m[1],
                                "label_seq_id": r["label"],
                                "author_residue": r["author"],
                                "reference_ca": r["xyz"],
                                "confidence": confidence[m[1]],
                                "secondary_structure": secondary,
                                "region": region,
                                "depth": depth,
                                "depth_ratio": ratio,
                                "identity_match": identity,
                                "altloc": r["altloc"],
                                "occupancy": r["occupancy"],
                            }
                        )
                    if sum(r["region"] == "TM_CORE" for r in residue_rows) < 20:
                        raise ValueError("TARGET_TM_CORE_SUPPORT_LT_20")
                    phrases, pfam = family_evidence(uni)
                    rec = {
                        "case_id": entry + ":" + ent + ":" + asym,
                        "entry_id": entry,
                        "entity_id": ent,
                        "label_asym_id": asym,
                        "author_chain_id": chain.name,
                        "model_id": str(refstructure[0].num),
                        "accession": acc,
                        "canonical_sequence": canonical,
                        "construct_sequence": seq,
                        "construct_sequence_sha256": construct_hash,
                        "exact_construct_prediction": seq == canonical,
                        "sequence_identity": sum(identities) / len(identities),
                        "observed_fraction": len(coords) / len(seq),
                        "mapping_coverage": len(mapped) / len(coords),
                        "uniprot_family_phrases": phrases,
                        "pfam_ids": pfam,
                        "opm_family_id": candidate["family_id"],
                        "opm_family_name": candidate["family_name_cache"],
                        "protein_name": candidate["name"],
                        "methods": methods,
                        "resolution": float(
                            re.match(r"[0-9.]+", str(candidate["resolution"])).group()
                        ),
                        "deposition_date": block.find_value(
                            "_pdbx_database_status.recvd_initial_deposition_date"
                        ),
                        "release_date": block.find_values(
                            "_pdbx_audit_revision_history.revision_date"
                        )[0]
                        if block.find_values("_pdbx_audit_revision_history.revision_date")
                        else None,
                        "opm_transform_rmsd": rmsd,
                        "half_thickness": half,
                        "opm_matched_atoms": len(hits),
                        "afdb_model_id": model.get("entryId"),
                        "afdb_version": model.get("latestVersion"),
                        "afdb_model_date": model.get("modelCreatedDate"),
                        "prediction_path": pred.relative_to(ROOT).as_posix(),
                        "prediction_sha256": sha(pred),
                        "ligands": sorted(set(block.find_values("_pdbx_entity_nonpoly.comp_id"))),
                        "assembly_details": list(
                            block.find_values("_pdbx_struct_assembly.details")
                        ),
                        "residues": residue_rows,
                    }
                    available.append(rec)
                    rec["sifts_accession"] = sifts_acc
                    rec["opm_reported_resolution"] = rec["resolution"]
                    rec["resolution"] = max(resolutions)
                except (ValueError, KeyError, urllib.error.HTTPError) as e:
                    rejected.append({"entity_id": ent, "label_asym_id": asym, "reason": str(e)})
        result = {
            "schema_version": "0.2.0",
            "entry_id": entry,
            "status": "prepared" if available else "excluded",
            "records": available,
            "exclusions": rejected,
        }
    except Exception as e:
        result = {
            "schema_version": "0.2.0",
            "entry_id": entry,
            "status": "source_or_parser_failure",
            "error_type": type(e).__name__,
            "error": str(e),
            "source_url": getattr(e, "url", None),
            "records": [],
        }
    result["acquisition_code_sha256"] = CODE_SHA256
    write(output, result)
    return result["status"]


def acquire(workers):
    candidates = read(WORK / "discovery.json")["candidates"]
    counts = Counter()
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for i, status in enumerate(pool.map(prepare_entry, candidates), 1):
            counts[status] += 1
            if i % 10 == 0:
                print(
                    json.dumps({"completed": i, "total": len(candidates), "status": dict(counts)}),
                    flush=True,
                )
    write(WORK / "acquisition_summary.json", dict(counts))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=["discover", "acquire"])
    p.add_argument("--limit", type=int, default=800)
    p.add_argument("--workers", type=int, default=4)
    a = p.parse_args()
    discover(a.limit) if a.command == "discover" else acquire(a.workers)
