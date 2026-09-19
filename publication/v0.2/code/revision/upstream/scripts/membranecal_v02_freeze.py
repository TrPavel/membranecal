"""Build outcome-blind v0.2 release. Does not import or evaluate an accuracy metric."""

from __future__ import annotations

# Machine-readable protocol explanations are intentionally uninterrupted strings.
# ruff: noqa: E501
import argparse
import concurrent.futures
import hashlib
import json
from collections import Counter, defaultdict
from itertools import combinations

import pyarrow as pa
import pyarrow.parquet as pq

from scripts.membranecal_v02_acquire import ROOT, V01, WORK, family_evidence, get, read, sha, write
from scripts.membranecal_v02_methods import (
    CONTRASTS,
    component_ids,
    family_names,
    grouping_edge,
    match_protein,
    pfam_complement,
    sequence_evidence,
)

RELEASE = ROOT / "data/frozen/membranecal_independent_validation_v0.2"


def prepare_development():
    """Fresh sequence metadata for development bridges; never old raw checkout files."""
    old = read(V01 / "proteins.json")["records"]

    def one(r):
        acc = r["accession"]
        uni = read(get("uniprot", acc, f"https://rest.uniprot.org/uniprotkb/{acc}.json", ".json"))
        phrases, pfam = family_evidence(uni)
        return {
            **r,
            "canonical_sequence": uni["sequence"]["value"],
            "current_primary_accession": uni["primaryAccession"],
            "secondary_accessions": uni.get("secondaryAccessions", []),
            "uniprot_family_phrases": phrases,
            "pfam_ids": pfam,
        }

    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        enriched = list(pool.map(one, old))
    write(
        WORK / "development_sequences.json",
        {
            "records": enriched,
            "note": "Current canonical sequences retrieved independently for all development proteins; identities linked to committed v0.1 accessions. No new or old accuracy used.",
        },
    )
    print(f"development sequence bridges: {len(enriched)}", flush=True)


def validate_groups():
    raw = read(WORK / "known_family_sequences.json")
    records = {}
    for acc, uni in raw.items():
        phrases, pfam = family_evidence(uni)
        records[acc] = {
            "accession": acc,
            "canonical_sequence": uni["sequence"]["value"],
            "uniprot_family_phrases": phrases,
            "pfam_ids": pfam,
        }
    results = []
    for a, b, expected in [
        ("P29972", "P41181", True),
        ("P08588", "P07550", True),
        ("P41181", "P07550", False),
    ]:
        ev = sequence_evidence(records[a]["canonical_sequence"], records[b]["canonical_sequence"])
        edge = grouping_edge(records[a], records[b], ev)
        results.append(
            {
                "a": a,
                "b": b,
                "expected_edge": expected,
                "edge": edge,
                "evidence": ev,
                "passed": bool(edge) == expected,
            }
        )
    seq = records["P29972"]["canonical_sequence"]
    fragment = sequence_evidence(seq, seq[: int(len(seq) * 0.4)])
    results.append(
        {
            "case": "40_percent_fragment",
            "passed": grouping_edge({}, {}, fragment) is None,
            "evidence": fragment,
        }
    )
    ev = {"identity": 0.8, "coverage_a": 0.6, "coverage_b": 0.6}
    a = {"pfam_ids": ["A", "B"], "uniprot_family_phrases": ["Belongs to the example family."]}
    b = {"pfam_ids": ["B", "C"], "uniprot_family_phrases": a["uniprot_family_phrases"]}
    results.append({"case": "A+B_vs_C+B_shared_domain", "passed": grouping_edge(a, b, ev) is None})
    # Explicitly record an unresolved limitation rather than falsely passing an
    # architecture test that the available identifiers cannot support.
    limitations = [
        "Pfam ID complements lack order and repeated-hit counts; A+B vs B+A and A+B vs A+B+B cannot be resolved from these metadata alone. Full-length sequence evidence and broad-family sensitivity reduce but do not eliminate this uncertainty."
    ]
    if not all(r["passed"] for r in results):
        raise ValueError("known-family validation failed")
    return {"records": records, "tests": results, "limitations": limitations}


def build_groups(records, old):
    all_records = {r["accession"]: r for r in old + records}
    evidence_cache_path = WORK / "pair_evidence.json"
    cache = read(evidence_cache_path) if evidence_cache_path.exists() else {}
    edges = []
    strict = []
    for i, (a, b) in enumerate(combinations(sorted(all_records), 2)):
        x, y = all_records[a], all_records[b]
        sx, sy = x["canonical_sequence"], y["canonical_sequence"]
        rescue_possible = (
            bool(family_names(x) & family_names(y))
            and bool(pfam_complement(x))
            and pfam_complement(x) == pfam_complement(y)
        )
        if min(len(sx), len(sy)) / max(len(sx), len(sy)) < (0.5 if rescue_possible else 0.8):
            continue
        key = hashlib.sha256(
            ("bio1.85-blosum62-local-10-.5-v1:" + sx + ":" + sy).encode()
        ).hexdigest()
        if key not in cache:
            cache[key] = sequence_evidence(sx, sy)
        ev = cache[key]
        kind = grouping_edge(x, y, ev)
        if kind:
            edges.append({"a": a, "b": b, "method": kind, **ev})
            if kind == "FULL_LENGTH_SEQUENCE":
                strict.append((a, b))
        if i and i % 10000 == 0:
            print(f"grouping pairs visited {i}; edges {len(edges)}", flush=True)
            write(evidence_cache_path, cache)
    write(evidence_cache_path, cache)
    raw_components = component_ids(all_records, [(e["a"], e["b"]) for e in edges])
    incident = {e[k] for e in edges for k in ("a", "b")}
    resolved = {a: bool(family_names(r)) or a in incident for a, r in all_records.items()}
    groups = {
        a: raw_components[a] if resolved[a] else "family:UNRESOLVED_POOLED" for a in all_records
    }
    curated_categories = defaultdict(list)
    for a, r in sorted(all_records.items()):
        for name in sorted(family_names(r)):
            curated_categories[name].append(a)
    curated_edges = [
        (members[0], a) for members in curated_categories.values() for a in members[1:]
    ]
    homology = component_ids(
        all_records, [(e["a"], e["b"]) for e in edges] + curated_edges, "homology:"
    )
    old_groups = {homology[r["accession"]] for r in old}
    novel = {
        r["accession"]
        for r in records
        if resolved[r["accession"]] and homology[r["accession"]] not in old_groups
    }
    # Merge primary groups sharing an experimental entry, including nonprimary
    # cohort proteins as bridges. Multiple chains never become independent draws.
    entry_groups = defaultdict(set)
    for r in records:
        entry_groups[r["entry_id"]].add(groups[r["accession"]])
    block_edges = [pair for fs in entry_groups.values() for pair in combinations(sorted(fs), 2)]
    block_edges.extend((groups[a], groups[b]) for a, b in curated_edges)
    block_by = component_ids(set(groups.values()), block_edges, "block:")
    broad_edges = [(e["a"], e["b"]) for e in edges]
    categories = defaultdict(list)
    for a, r in all_records.items():
        for f in family_names(r):
            categories["curated:" + f].append(a)
        if r.get("opm_family_id") is not None:
            categories["opm:" + str(r["opm_family_id"])].append(a)
    for members in categories.values():
        broad_edges.extend((members[0], a) for a in members[1:])
    broad = component_ids(all_records, broad_edges, "broad:")
    strict_groups = component_ids(all_records, strict, "strict:")
    return {
        "schema_version": "0.2.0",
        "family_by_accession": groups,
        "raw_component_by_accession": raw_components,
        "resolved_by_accession": resolved,
        "block_by_family": block_by,
        "curated_dependence_edges": curated_edges,
        "curated_homology_by_accession": homology,
        "family_novel_accessions": sorted(novel),
        "development_accessions": sorted(r["accession"] for r in old),
        "edges": edges,
        "broad_family_by_accession": broad,
        "strict_family_by_accession": strict_groups,
        "limitation": "Operational groups validated on a bounded panel, not established evolutionary-independent units. Curated singleton evidence can be incomplete; unresolved cases pooled; domain order and repeated-domain counts unavailable.",
    }


def build():
    if (RELEASE / "PRE_OUTCOME_FREEZE.json").exists():
        raise ValueError("freeze already sealed; refuse overwrite")
    discovery = read(WORK / "discovery.json")
    prepared = [
        read(WORK / "prepared" / (c["pdbid"].lower() + ".json")) for c in discovery["candidates"]
    ]
    current_sha = sha(ROOT / "scripts/membranecal_v02_acquire.py")
    if any(p["acquisition_code_sha256"] != current_sha for p in prepared):
        raise ValueError("acquisition implementation changed: complete pre-outcome reprepare first")
    old = read(WORK / "development_sequences.json")["records"]
    exclusions = [
        {
            "entry_id": p["entry_id"],
            "status": p["status"],
            "error": p.get("error"),
            "entities": p.get("exclusions", []),
        }
        for p in prepared
    ]
    candidates = [r for p in prepared for r in p["records"]]
    records, seen, constructs = [], set(), set()
    for r in sorted(
        candidates,
        key=lambda r: (
            r["resolution"],
            -len(r["residues"]),
            hashlib.sha256(r["case_id"].encode()).hexdigest(),
        ),
    ):
        reason = (
            "V02_CANONICAL_DUPLICATE"
            if r["accession"] in seen
            else "V02_CONSTRUCT_DUPLICATE"
            if r["construct_sequence_sha256"] in constructs
            else None
        )
        if reason:
            exclusions.append({"case_id": r["case_id"], "reason": reason})
            continue
        records.append(r)
        seen.add(r["accession"])
        constructs.add(r["construct_sequence_sha256"])
    records.sort(key=lambda r: r["accession"])
    old_aliases = {
        a
        for r in old
        for a in [r["accession"], r["current_primary_accession"], *r["secondary_accessions"]]
    }
    assert not seen & old_aliases
    assert not {r["entry_id"] for r in records} & {r["entry_id"] for r in old}
    assert not constructs & {r["construct_sequence_sha256"] for r in old}
    validations = validate_groups()
    groups = build_groups(records, old)
    matching, diagnostics, feasibility = {}, {}, {}
    for name, width, compatible in [
        ("primary", 2, False),
        ("bins_1", 1, False),
        ("bins_5", 5, False),
        ("v01_compatible", 2, True),
    ]:
        sets, diag = [], []
        for r in records:
            m, d = match_protein(r, width=width, compatible=compatible)
            sets.extend(m)
            diag.extend(d)
        matching[name], diagnostics[name] = sets, diag
        fs = {}
        for left, right in CONTRASTS:
            contrast = left + "-" + right
            accessions = {m["accession"] for m in sets if m["contrast"] == contrast}
            family_ids = {
                groups["family_by_accession"][a]
                for a in accessions
                if groups["resolved_by_accession"][a]
            }
            blocks = {groups["block_by_family"][f] for f in family_ids}
            selected = [d for d in diag if d["contrast"] == contrast]
            available = sum(d["available_left"] + d["available_right"] for d in selected)
            fs[contrast] = {
                "proteins": len(accessions),
                "resolved_families": len(family_ids),
                "resolved_blocks": len(blocks),
                "useful": len(accessions) >= 30 and len(family_ids) >= 20 and len(blocks) >= 20,
                "raw_cell_retention": sum(d["retained_residues"] for d in selected) / available
                if available
                else 0,
                "effective_mass_retention": sum(
                    2 * d["eligible_mass"] for d in selected if d["included"]
                )
                / available
                if available
                else 0,
                "max_absolute_confidence_difference": max(
                    (abs(m["confidence_difference"]) for m in sets if m["contrast"] == contrast),
                    default=None,
                ),
            }
        feasibility[name] = fs
    protocol = {
        "schema_version": "0.2.0",
        "stage": "PRE_OUTCOME",
        "parent_commit": "51e8feb0b29fc089f40c238a4b4a72657eef7ae4",
        "cohort": {
            "candidate_cap": 800,
            "scope": "alpha-helical transmembrane OPM, resolution <=3.5 A, X-ray or EM",
            "independence": [
                "canonical accession absent from all 174 v0.1 accepted",
                "experimental entry absent from all v0.1 accepted",
                "exact construct SHA256 absent from all v0.1 accepted",
            ],
            "canonical_length": [100, 1500],
            "construct_canonical_length_ratio_min": 0.5,
            "experimental_observed_fraction_min": 0.7,
            "mapped_observed_fraction_min": 0.95,
            "mapped_residues_min": 100,
            "identity_min": 0.98,
            "primary_core_residues_min": 20,
            "opm_chain_fraction_min": 0.8,
            "opm_rmsd_max_angstrom": 0.5,
            "opm_half_thickness_disagreement_max_angstrom": 0.1,
            "representative": "one deposited label-asym per entity; one canonical/construct across entries sorted resolution, descending mapped residues, SHA256(case ID)",
        },
        "mapping": {
            "source": "SIFTS XML author-chain/residue to UniProt canonical; positive unique positions; reference highest-occupancy CA with altloc lexical tie break; no prediction alignment fitted",
            "noncanonical_isoforms": "excluded",
            "mutations": "<=2% retained and separately sensitivity-excluded",
            "missing_coordinates": "no imputation; endpoint uses paired observed CA only; zero neighbours undefined",
            "missing_frozen_matched_endpoint": "abort primary computation rather than silently rematch or drop outcomes",
            "experimental_model": "first deposited model; selected before outcome access",
        },
        "regions": {
            "source": "OPM oriented experimental coordinates via Kabsch transfer",
            "depth_ratio": "abs(z)/hydrophobic_half_thickness",
            "core_max": 0.8,
            "interface_max": 1.2,
            "extramembrane_min_exclusive": 1.6,
            "primary_secondary_structure": ["HELIX", "SHEET"],
            "other": "UNKNOWN_OR_EXCLUDED; unavailable does not mean disorder",
        },
        "groups": {
            "implementation": "scripts/membranecal_v02_methods.py",
            "aligner": "Biopython 1.85 local affine BLOSUM62; open -10, extend -0.5; no E-value claim",
            "sequence_identity_min": 0.3,
            "full_length_paired_coverage_each_min": 0.8,
            "rescue_paired_coverage_each_min": 0.5,
            "rescue_additional": "same curated family and identical nonempty Pfam complement; shared Pfam alone never sufficient",
            "multidomain": "whole-sequence paired coverage; different complements prohibit rescue; order/repeat uncertainty explicitly unresolved by available metadata",
            "unresolved": "single pooled family; omitted from family-novel; exclusion sensitivity",
            "development_bridges": True,
            "curated_dependence": "same exact normalized explicit whole-family membership links primary bootstrap blocks regardless of sequence threshold; normalize capitalization/whitespace/leading articles; no generic function/domain-only labels",
            "family_novel": "combined sequence/rescue plus curated homology graph excludes components touching development; unresolved excluded; structure-sharing alone does not define homology",
        },
        "matching": {
            "method": "within-protein coarsened exact matching",
            "plddt_bin_width_points": 2,
            "bin_edges": "[0,2),...,[98,100]",
            "bin_arithmetic": "pLDDT points rounded to 6 decimals to avoid binary representation edge drift",
            "secondary_structure": "exact HELIX/SHEET",
            "min_cell_each_region": 3,
            "min_total_common_mass": 20,
            "cell_weight": "min(n_left,n_right), identical on both sides",
            "residue_weight": "cell mass / regional cell residue count; normalize by protein common mass",
        },
        "primary_contrasts": [a + "-" + b for a, b in CONTRASTS[:2]],
        "secondary_contrasts": [CONTRASTS[2][0] + "-" + CONTRASTS[2][1]],
        "effect": "regional (local_CA_agreement - pLDDT/100) difference; equal proteins per operational sequence/rescue group then equal group means; not equal curated-family weights",
        "metric": {
            "name": "superposition-free CA local-distance agreement",
            "reference_neighbor_radius_angstrom": 15,
            "radius_inclusive": True,
            "thresholds_angstrom": [0.5, 1, 2, 4],
            "thresholds_inclusive": True,
            "sequence_separation_exclusion": 0,
            "neighbors": "nonself paired reference CA; finite coordinates required",
            "isolated_residue": "undefined, no zero imputation",
            "secondary_implementation": "independent scalar pair loop on first five accession-sorted proteins after unseal",
        },
        "multiplicity": {
            "method": "Bonferroni two primary contrasts",
            "primary_interval_level": 0.975,
            "secondary_interval_level": 0.95,
            "quantiles_primary": [0.0125, 0.9875],
        },
        "uncertainty": {
            "bootstrap_draws": 20000,
            "seed": 20909,
            "unit": "connected operational-group plus explicit curated-whole-family plus shared experimental-entry blocks; retain development/nonprimary bridges; equal-operational-group ratio estimator",
            "low_block_intervals": "descriptive only; cannot establish primary conclusion",
        },
        "planning": {
            "expected_candidates": 800,
            "expected_attrition_fraction": [0.4, 0.6],
            "expected_common_support_retention": [0.25, 0.6],
            "min_paired_proteins": 30,
            "min_resolved_families": 20,
            "min_resolved_blocks": 20,
            "aim_groups": 40,
            "material_difference": 0.05,
            "equivalence_test": False,
        },
        "sensitivities": [
            "1-point bins",
            "5-point bins",
            "v0.1-compatible SS mask and no SS matching",
            "family-novel relative to combined development graph, exclude unresolved",
            "exact full construct sequence",
            "unresolved excluded",
            "broad curated/OPM family merges",
            "equal curated-homology-component weighting",
            "strict full-length sequence groups",
            "identity 100% only",
            "release after 2018-04-30 (date-based only; no known training membership claim)",
        ],
        "kill_criteria": [
            "any development overlap",
            "outcomes accessed before dedicated freeze commit",
            "source checksum mismatch",
            "no defensible membrane frame/mapping",
            "known-family validation failure",
            "either primary below 30 proteins/20 resolved families/20 blocks: independent cohort not feasible for full question",
        ],
        "decision_rule": "research/membranecal_v02_design.md: symmetric material margin .05; feasibility first, then predeclared method-sensitivity flag, then adjusted-interval material support, otherwise no material effect established",
        "protocol_deviations": [],
        "outcomes_computed": False,
    }
    RELEASE.mkdir(parents=True, exist_ok=True)
    residues = []
    proteins = []
    for r in records:
        proteins.append({k: v for k, v in r.items() if k != "residues"})
        residues.extend({"accession": r["accession"], **row} for row in r["residues"])
    pq.write_table(
        pa.Table.from_pylist(residues), RELEASE / "residues_preoutcome.parquet", compression="zstd"
    )
    write(RELEASE / "protocol.json", protocol)
    write(RELEASE / "proteins.json", {"schema_version": "0.2.0", "records": proteins})
    write(RELEASE / "groups.json", groups)
    write(RELEASE / "group_validation.json", validations)
    write(RELEASE / "development_sequences.json", read(WORK / "development_sequences.json"))
    write(RELEASE / "matched_sets.json", matching)
    write(RELEASE / "common_support.json", diagnostics)
    write(RELEASE / "feasibility.json", feasibility)
    write(RELEASE / "exclusions.json", exclusions)
    write(RELEASE / "discovery.json", discovery)
    source_manifest = [read(p) for p in sorted((WORK / "raw").rglob("*.source.json"))]
    write(RELEASE / "sources.json", {"schema_version": "0.2.0", "records": source_manifest})
    summary = {
        "candidate_entries": len(prepared),
        "source_status": dict(Counter(p["status"] for p in prepared)),
        "accepted_proteins": len(records),
        "accepted_entries": len({r["entry_id"] for r in records}),
        "residues": len(residues),
        "resolved_families": len(
            {
                groups["family_by_accession"][r["accession"]]
                for r in records
                if groups["resolved_by_accession"][r["accession"]]
            }
        ),
        "unresolved_proteins": sum(
            not groups["resolved_by_accession"][r["accession"]] for r in records
        ),
        "family_novel_proteins": len(groups["family_novel_accessions"]),
        "exact_construct_proteins": sum(r["exact_construct_prediction"] for r in records),
        "development_overlap": {"canonical": 0, "entry": 0, "construct": 0},
    }
    write(RELEASE / "cohort_summary.json", summary)
    files = list(RELEASE.glob("*")) + [
        ROOT / "scripts/membranecal_v02_acquire.py",
        ROOT / "scripts/membranecal_v02_methods.py",
        ROOT / "scripts/membranecal_v02_freeze.py",
        ROOT / "scripts/membranecal_v02_analysis.py",
        ROOT / "research/membranecal_v02_design.md",
        ROOT / "data/licences/membranecal_independent_validation_v0.2.json",
    ]
    seal = {
        "schema_version": "0.2.0",
        "identity": "PRE_OUTCOME_FREEZE",
        "outcomes_computed": False,
        "files": {p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(files) if p.is_file()},
    }
    write(RELEASE / "PRE_OUTCOME_FREEZE.json", seal)
    print(json.dumps({"cohort": summary, "feasibility": feasibility}, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["development", "freeze"])
    args = parser.parse_args()
    prepare_development() if args.command == "development" else build()
