import math
from pathlib import Path
from collections import defaultdict
from itertools import combinations
import numpy as np
import pyarrow.parquet as pq
from .common import ROOT, RELEASE, REPORT, read, sha, write, verify_freeze
from .methods import CONTRASTS, REGIONS, cluster_interval, component_ids
def matrix_scores(reference, prediction):
    positions = sorted(set(reference) & set(prediction))
    if not positions:
        return {}
    ref = np.array([reference[k] for k in positions], dtype=float)
    pred = np.array([prediction[k] for k in positions], dtype=float)
    if not np.isfinite(ref).all() or not np.isfinite(pred).all():
        raise ValueError("nonfinite coordinates")
    rd = np.sqrt(np.sum((ref[:, None] - ref[None, :]) ** 2, axis=2))
    pd = np.sqrt(np.sum((pred[:, None] - pred[None, :]) ** 2, axis=2))
    mask = (rd <= 15.0) & ~np.eye(len(positions), dtype=bool)
    counts = mask.sum(axis=1)
    hits = sum(
        ((np.abs(rd - pd) <= threshold) & mask).sum(axis=1) for threshold in (0.5, 1.0, 2.0, 4.0)
    )
    return {
        p: {
            "local_accuracy": float(hits[i] / (4 * counts[i])) if counts[i] else None,
            "neighbors": int(counts[i]),
        }
        for i, p in enumerate(positions)
    }

def scalar_scores(reference, prediction):
    """Independent scalar arithmetic, no NumPy/matrix helper reused."""
    result = {}
    for p in sorted(set(reference) & set(prediction)):
        hits = neighbors = 0
        for q in sorted(set(reference) & set(prediction)):
            if p == q:
                continue
            rd = math.dist(reference[p], reference[q])
            if rd > 15:
                continue
            neighbors += 1
            delta = abs(rd - math.dist(prediction[p], prediction[q]))
            hits += sum(delta <= threshold for threshold in (0.5, 1, 2, 4))
        result[p] = {
            "local_accuracy": hits / (4 * neighbors) if neighbors else None,
            "neighbors": neighbors,
        }
    return result

def load_inputs():
    proteins = read(RELEASE / "proteins.json")["records"]
    rows = pq.read_table(RELEASE / "residues_preoutcome.parquet").to_pylist()
    grouped = defaultdict(list)
    for r in rows:
        grouped[r["accession"]].append(r)
    return proteins, grouped

def verify_score_manifest():
    manifest = read(RELEASE / "score_manifest.json")
    if manifest["freeze_commit"] != read(RELEASE / "unseal_receipt.json")["freeze_commit"]:
        raise ValueError("score freeze identity mismatch")
    for name, expected in manifest["files"].items():
        if sha(RELEASE / name) != expected:
            raise ValueError("score artefact checksum mismatch: " + name)

def contrast_effects(matches, errors, permitted):
    results = {}
    for match in matches:
        acc = match["accession"]
        if acc not in permitted:
            continue
        value = 0.0
        for cell in match["cells"]:
            left = [errors[(acc, p)] for p in cell["left_positions"]]
            right = [errors[(acc, p)] for p in cell["right_positions"]]
            if any(x is None for x in left + right):
                raise ValueError(
                    "Frozen matched endpoint unavailable; do not silently rematch or drop outcome-dependent cases"
                )
            value += cell["mass"] * (float(np.mean(left)) - float(np.mean(right)))
        results[acc] = value / match["mass"]
    return results

def block_map(proteins, family_by, curated_edges=()):
    entries = defaultdict(set)
    for p in proteins:
        entries[p["entry_id"]].add(family_by[p["accession"]])
    edges = [edge for fs in entries.values() for edge in combinations(sorted(fs), 2)]
    edges.extend((family_by[a], family_by[b]) for a, b in curated_edges)
    return component_ids(set(family_by.values()), edges, "block:")

def mean_or_none(values):
    return float(np.mean(values)) if len(values) else None

def analyze(output):
    receipt = read(RELEASE / "unseal_receipt.json")
    verify_freeze(receipt["freeze_commit"])
    verify_score_manifest()
    proteins, grouped = load_inputs()
    protein_by = {p["accession"]: p for p in proteins}
    scores = {
        (r["accession"], r["canonical_position"]): r
        for r in pq.read_table(RELEASE / "scores.parquet").to_pylist()
    }
    errors = {}
    pooled = []
    for acc, rs in grouped.items():
        for r in rs:
            local = scores[(acc, r["canonical_position"])]["local_accuracy"]
            error = local - r["confidence"] if local is not None else None
            errors[(acc, r["canonical_position"])] = error
            if local is not None:
                pooled.append({**r, "local_accuracy": local, "error": error})
    grouping = read(RELEASE / "groups.json")
    matches = read(RELEASE / "matched_sets.json")
    family = grouping["family_by_accession"]
    all_acc = set(protein_by)
    resolved = {a for a in all_acc if grouping["resolved_by_accession"][a]}
    scenarios = {name: (name, all_acc, family) for name in matches}
    scenarios.update(
        {
            "family_novel": ("primary", set(grouping["family_novel_accessions"]), family),
            "exact_construct": (
                "primary",
                {p["accession"] for p in proteins if p["exact_construct_prediction"]},
                family,
            ),
            "unresolved_excluded": ("primary", resolved, family),
            "broad_family": (
                "primary",
                all_acc,
                {
                    a: (
                        grouping["broad_family_by_accession"][a]
                        if grouping["resolved_by_accession"][a]
                        else "broad:UNRESOLVED_POOLED"
                    )
                    for a in family
                },
            ),
            "strict_sequence_groups": (
                "primary",
                all_acc,
                {
                    a: (
                        grouping["strict_family_by_accession"][a]
                        if grouping["resolved_by_accession"][a]
                        else "strict:UNRESOLVED_POOLED"
                    )
                    for a in family
                },
            ),
            "identity_100": (
                "primary",
                {p["accession"] for p in proteins if p["sequence_identity"] == 1},
                family,
            ),
            "curated_homology": (
                "primary",
                all_acc,
                {
                    a: (
                        grouping["curated_homology_by_accession"][a]
                        if grouping["resolved_by_accession"][a]
                        else "homology:UNRESOLVED_POOLED"
                    )
                    for a in family
                },
            ),
            "post_2018_release": (
                "primary",
                {
                    p["accession"]
                    for p in proteins
                    if p["release_date"] and p["release_date"] > "2018-04-30"
                },
                family,
            ),
        }
    )
    results = {}
    for scenario, (match_name, permitted, family_by) in scenarios.items():
        bs = block_map(proteins, family_by, grouping["curated_dependence_edges"])
        contrasts = {}
        for i, (left, right) in enumerate(CONTRASTS):
            name = left + "-" + right
            selected = [m for m in matches[match_name] if m["contrast"] == name]
            effects = contrast_effects(selected, errors, permitted)
            stats = cluster_interval(effects, family_by, bs, primary=i < 2)
            resolved_families = {family_by[a] for a in effects if a in resolved}
            resolved_blocks = {bs[f] for f in resolved_families}
            stats.update(
                protein_effects=effects,
                resolved_families=len(resolved_families),
                resolved_blocks=len(resolved_blocks),
                useful=len(effects) >= 30
                and len(resolved_families) >= 20
                and len(resolved_blocks) >= 20,
            )
            stats["interval_interpretation"] = (
                "primary approximate clustered interval"
                if scenario == "primary" and stats["useful"] and i < 2
                else "descriptive sensitivity/secondary or insufficient support"
            )
            contrasts[name] = stats
        results[scenario] = contrasts
    descriptive = {}
    for region in ["ALL"] + list(REGIONS):
        rs = pooled if region == "ALL" else [r for r in pooled if r["region"] == region]
        byprotein = defaultdict(list)
        for r in rs:
            byprotein[r["accession"]].append(r["error"])
        byfamily = defaultdict(list)
        for acc, values in byprotein.items():
            byfamily[family[acc]].append(float(np.mean(values)))
        descriptive[region] = {
            "residues": len(rs),
            "proteins": len(byprotein),
            "families": len(byfamily),
            "pooled_error": mean_or_none([r["error"] for r in rs]),
            "protein_balanced_error": mean_or_none([np.mean(v) for v in byprotein.values()]),
            "family_balanced_error": mean_or_none([np.mean(v) for v in byfamily.values()]),
            "mean_confidence": mean_or_none([r["confidence"] for r in rs]),
            "mean_accuracy": mean_or_none([r["local_accuracy"] for r in rs]),
        }
    bins = []
    for lo in range(0, 100, 10):
        rs = [
            r
            for r in pooled
            if lo <= r["confidence"] * 100 < lo + 10 or (lo == 90 and r["confidence"] == 1)
        ]
        bins.append(
            {
                "lower_points": lo,
                "upper_points": lo + 10,
                "residues": len(rs),
                "confidence": mean_or_none([r["confidence"] for r in rs]),
                "accuracy": mean_or_none([r["local_accuracy"] for r in rs]),
            }
        )
    overall_corr = (
        float(
            np.corrcoef([r["confidence"] for r in pooled], [r["local_accuracy"] for r in pooled])[
                0, 1
            ]
        )
        if len(pooled) > 1
        else None
    )
    failure_cases = []
    for acc, rs in grouped.items():
        values = [
            errors[(acc, r["canonical_position"])]
            for r in rs
            if errors[(acc, r["canonical_position"])] is not None
        ]
        p = protein_by[acc]
        failure_cases.append(
            {
                "accession": acc,
                "entry_id": p["entry_id"],
                "mean_error": mean_or_none(values),
                "construct_mismatch": "sequence differs from canonical"
                if not p["exact_construct_prediction"]
                else "full sequence identical",
                "substitution_fraction": 1 - p["sequence_identity"],
                "state_conformation_mismatch": "unresolved: coordinate discrepancy alone cannot establish state",
                "mapping_failure": "no failure detected by frozen QC; residual mapping error remains possible",
                "assembly_context": p["assembly_details"],
                "ligands": p["ligands"],
                "missing_experimental_fraction": 1 - p["observed_fraction"],
                "membrane_annotation_uncertainty": {
                    "opm_transfer_rmsd": p["opm_transform_rmsd"],
                    "half_thickness": p["half_thickness"],
                    "caveat": "computed orientation; no independent bilayer observation",
                },
                "confidence_failure": "positive confidence-minus-agreement discrepancy is operational overconfidence; biological mechanism unassigned",
            }
        )
    failure_cases.sort(key=lambda r: (-abs(r["mean_error"] or 0), r["accession"]))
    primary = results["primary"]
    primary_names = [a + "-" + b for a, b in CONTRASTS[:2]]

    def material(s):
        ci = s["interval"]
        return 1 if ci and ci[0] > 0.05 else -1 if ci and ci[1] < -0.05 else 0

    sensitivity_flags = []
    for name in primary_names:
        p = primary[name]
        for scenario in (
            "bins_1",
            "bins_5",
            "family_novel",
            "unresolved_excluded",
            "broad_family",
            "curated_homology",
            "v01_compatible",
        ):
            s = results[scenario][name]
            if not s["useful"] or p["estimate"] is None:
                continue
            shifted_crossing = (
                abs(s["estimate"] - p["estimate"]) >= 0.05 and s["estimate"] * p["estimate"] < 0
            )
            reversed_material = material(p) and s["estimate"] * material(p) < 0
            if shifted_crossing or reversed_material:
                sensitivity_flags.append(
                    {
                        "contrast": name,
                        "sensitivity": scenario,
                        "shift": s["estimate"] - p["estimate"],
                    }
                )
    if not all(primary[n]["useful"] for n in primary_names):
        decision = "MEMBRANECAL_V0_2_INDEPENDENT_COHORT_NOT_FEASIBLE"
    elif sensitivity_flags:
        decision = "MEMBRANECAL_V0_2_METHOD_SENSITIVITY_REMAINS_DOMINANT"
    elif any(material(primary[n]) for n in primary_names):
        decision = "MEMBRANECAL_V0_2_REGIONAL_EFFECT_SUPPORTED"
    else:
        decision = "MEMBRANECAL_V0_2_NO_MATERIAL_REGIONAL_EFFECT_ESTABLISHED"
    result = {
        "schema_version": "0.2.0",
        "milestone": "MEMBRANECAL_INDEPENDENT_VALIDATION_V0_2_COMPLETE",
        "decision": decision,
        "freeze_commit": receipt["freeze_commit"],
        "protocol_sha256": sha(RELEASE / "protocol.json"),
        "cohort": read(RELEASE / "cohort_summary.json"),
        "scored_residues": len(pooled),
        "undefined_residues": len(scores) - len(pooled),
        "contrasts": results,
        "descriptive": descriptive,
        "confidence_bins": bins,
        "pooled_confidence_accuracy_pearson": overall_corr,
        "sensitivity_flags": sensitivity_flags,
        "failure_cases": failure_cases[:20],
        "temporal_caveat": "2018-04-30 is the published original AlphaFold training structure cutoff used as a date-based diagnostic. AFDB model version/date is recorded individually. Post-cutoff does not establish absence of later model/template/sequence exposure; no individual training membership is claimed.",
        "protocol_deviations": [],
    }
    output = Path(output)
    write(output / "results.json", result)
    cohort = result["cohort"]
    lines = [
        "# MembraneCal independent validation v0.2",
        "",
        result["milestone"],
        "",
        decision,
        "",
        f"PRE_OUTCOME_FREEZE: `{receipt['freeze_commit']}`; protocol SHA256 `{result['protocol_sha256']}`.",
        "",
        "## Cohort",
        "",
        f"{cohort['candidate_entries']} candidate entries; {cohort['accepted_proteins']} accepted canonical proteins in {cohort['accepted_entries']} experimental entries; {cohort['residues']:,} mapped residues. {cohort['resolved_families']} resolved operational groups; {cohort['unresolved_proteins']} unresolved proteins. Zero development canonical/alias, entry or construct overlaps. Family-novel subset: {cohort['family_novel_proteins']} proteins; full-sequence exact-construct subset: {cohort['exact_construct_proteins']} proteins.",
        "",
        "## Matched primary results",
        "",
        "Effects and intervals are in percentage points. Equal proteins within operational sequence/rescue groups and equal group means; blocks additionally link explicit curated whole-family annotations and shared entries, 20,000 bootstrap draws, Bonferroni 97.5% intervals. Equal curated-homology weighting is a separate sensitivity.",
        "",
        "| Contrast | Effect | Interval | Proteins | Resolved groups | Blocks | Useful |",
        "|---|---:|---|---:|---:|---:|---|",
    ]
    for name in primary_names:
        s = primary[name]
        effect = f"{100 * s['estimate']:.3f}" if s["estimate"] is not None else "undefined"
        interval = (
            ", ".join(f"{100 * v:.3f}" for v in s["interval"]) if s["interval"] else "undefined"
        )
        lines.append(
            f"| {name} | {effect} | [{interval}] | {s['proteins']} | {s['resolved_families']} | {s['resolved_blocks']} | {s['useful']} |"
        )
    lines += [
        "",
        "## Sensitivity estimates",
        "",
        "Intervals below remain descriptive sensitivity intervals. N/F/B are paired proteins, resolved groups and resolved dependence blocks.",
        "",
        "| Scenario | Contrast | Effect, points | Interval, points | N/F/B | Useful support |",
        "|---|---|---:|---|---|---|",
    ]
    for scenario in results:
        if scenario == "primary":
            continue
        for name in primary_names:
            s = results[scenario][name]
            effect = f"{100 * s['estimate']:.3f}" if s["estimate"] is not None else "undefined"
            ci = (
                ", ".join(f"{100 * v:.3f}" for v in s["interval"]) if s["interval"] else "undefined"
            )
            lines.append(
                f"| {scenario} | {name} | {effect} | [{ci}] | {s['proteins']}/{s['resolved_families']}/{s['resolved_blocks']} | {s['useful']} |"
            )
    lines += [
        "",
        "## Common support",
        "",
        "| Contrast | Raw retained-cell fraction | Effective mass fraction | Max absolute residual pLDDT difference, points |",
        "|---|---:|---:|---:|",
    ]
    for name, s in read(RELEASE / "feasibility.json")["primary"].items():
        imbalance = (
            100 * s["max_absolute_confidence_difference"]
            if s["max_absolute_confidence_difference"] is not None
            else 0
        )
        lines.append(
            f"| {name} | {s['raw_cell_retention']:.3f} | {s['effective_mass_retention']:.3f} | {imbalance:.3f} |"
        )
    lines += [
        "",
        "## Descriptive calibration",
        "",
        f"Pooled confidence–accuracy Pearson correlation: {overall_corr:.4f}. This residue-level description is not an independent-replicate inferential test.",
        "",
        "| Region | Residues | Pooled error, points | Protein-balanced error | Family-balanced error |",
        "|---|---:|---:|---:|---:|",
    ]
    for region, s in descriptive.items():
        values = [
            f"{100 * s[k]:.3f}" if s[k] is not None else "undefined"
            for k in ["pooled_error", "protein_balanced_error", "family_balanced_error"]
        ]
        lines.append(f"| {region} | {s['residues']} | {' | '.join(values)} |")
    lines += ["", "## Family influence", ""]
    for name in primary_names:
        s = primary[name]
        influential = sorted(
            s["leave_one_family_out"].items(), key=lambda x: (-abs(x[1] - s["estimate"]), x[0])
        )[:3]
        lines.append(
            name
            + ": "
            + "; ".join(
                f"omit {f}: {100 * v:.3f} points ({100 * (v - s['estimate']):+.3f} shift)"
                for f, v in influential
            )
            + "."
        )
    lines += [
        "",
        "## Largest protein-average discrepancies",
        "",
        "These are descriptive case investigations; none was excluded after outcomes were seen.",
        "",
        "| Protein | PDB | Mean error, points | Construct | Missing experimental fraction |",
        "|---|---|---:|---|---:|",
    ]
    for case in failure_cases[:8]:
        lines.append(
            f"| {case['accession']} | {case['entry_id']} | {100 * case['mean_error']:.3f} | {case['construct_mismatch']} | {case['missing_experimental_fraction']:.3f} |"
        )
    lines += [
        "",
        "## Matching and grouping",
        "",
        "Frozen 2-point pLDDT bins × HELIX/SHEET; >=3 residues per side/cell and >=20 common mass. Raw retained-cell fraction differs from effective mass retention: see frozen feasibility.json/common_support.json. Continuous confidence equality is not claimed.",
        "",
        "Combined development/validation sequence graph, full-length >=30% identity and >=80% paired coverage; curated-family plus complete Pfam-complement plus >=50% coverage rescue. Shared Pfam alone never creates an edge. Primary bootstrap blocks additionally join exact normalized whole-family annotations regardless of sequence threshold. Family novelty uses those homology links, not shared-entry links. Bounded known-family validation passes; ordered architecture/repeat uncertainty remains. Unresolved proteins are conservatively pooled. These are operational groups, not proven evolutionary-independent units.",
        "",
        "## Sensitivities and leakage",
        "",
        "results.json includes all frozen 1/5-bin, v0.1-compatible, family-novel, exact-construct, unresolved-excluded, broad/strict grouping, 100%-identity and temporal sensitivities, plus family means and leave-one-family-out estimates.",
        "",
        result["temporal_caveat"],
        "",
        "## Interpretation",
        "",
        "Material support requires an adjusted interval entirely beyond ±5 points. No material effect established is not equivalence or proof of a zero effect. Scope is protein-independent alpha-helical transmembrane canonical-AFDB benchmarking with observed structural support. Experimental state/construct and computational orientation limitations remain.",
        "",
        "## Failure modes",
        "",
        "Twenty largest absolute protein-average errors are recorded with construct, mutation, assembly, missingness and membrane-frame evidence. State mismatches and mechanistic explanations remain unresolved without direct evidence; no post-outcome exclusion or retuning.",
        "",
        "## Reproduction",
        "",
        "Run `python -m membranecal_release.verify` using environment/requirements.lock. This replays compact-coordinate scores and all statistical outputs offline. Raw source archives are omitted; source terms and attribution are documented in DATA_LICENSE.md.",
        "",
    ]
    output.mkdir(parents=True, exist_ok=True)
    (output / "report.md").write_text("\n".join(lines), encoding="utf-8", newline="\n")
    write(
        output / "checksums.json",
        {
            p.name: sha(p)
            for p in sorted(output.glob("*"))
            if p.is_file() and p.name != "checksums.json"
        },
    )
    print(decision, flush=True)

def verify_scores():
    verify_score_manifest()
    proteins, grouped = load_inputs()
    scores = defaultdict(list)
    for r in pq.read_table(RELEASE / "scores.parquet").to_pylist():
        scores[r["accession"]].append(r)
    for p in proteins:
        acc = p["accession"]
        ref = {r["canonical_position"]: r["reference_ca"] for r in grouped[acc]}
        pred = {
            r["canonical_position"]: r["prediction_ca"]
            for r in scores[acc]
            if r["prediction_ca"] is not None
        }
        computed = matrix_scores(ref, pred)
        for r in scores[acc]:
            if (
                computed.get(r["canonical_position"], {}).get("local_accuracy")
                != r["local_accuracy"]
            ):
                raise ValueError("compact metric replay mismatch")
    print(f"compact coordinate metric replay verified: {len(proteins)} proteins", flush=True)
