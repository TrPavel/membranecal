"""Independent closure audit. Does not import the frozen study's scientific functions."""

from __future__ import annotations

# Narrative audit qualifications are deliberately uninterrupted strings.
# ruff: noqa: E501
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

SOURCE = "d9641f28fdf50d42057223e3787f999e1a127198"
FREEZE = "5deb8f9abc345d811e0e98680b6122b4046a0d0d"
DATA = Path("data")
RESULT = Path("results/results.json")
REGIONS = ["TM_CORE", "MEMBRANE_INTERFACE", "STRUCTURED_EXTRAMEMBRANE"]
PAIRS = [(REGIONS[1], REGIONS[0]), (REGIONS[2], REGIONS[0]), (REGIONS[2], REGIONS[1])]
TOLERANCE = 1e-12


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def close(a, b):
    if not math.isclose(a, b, rel_tol=0, abs_tol=TOLERANCE):
        raise AssertionError(f"independent reproduction mismatch: {a} versus {b}")


def components(nodes, edges, prefix):
    """Breadth-first graph traversal, independently of the frozen union-find code."""
    adjacency = {n: set() for n in nodes}
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    result = {}
    for start in sorted(adjacency):
        if start in result:
            continue
        stack, visited = [start], {start}
        while stack:
            for neighbor in adjacency[stack.pop()]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        result.update(dict.fromkeys(visited, prefix + min(visited)))
    return result


def names(record):
    labels = [t[8:] for t in record.get("family_tokens", []) if t.startswith("uniprot:")]
    for paragraph in record.get("uniprot_family_phrases", []):
        labels += re.findall(r"Belongs to (.*? family)(?:\.|;|$)", paragraph)
    return {re.sub(r"^(the|an|a) ", "", " ".join(x.lower().split())) for x in labels}


def region(row, compatible):
    z = row["depth_ratio"]
    value = REGIONS[0] if z <= 0.8 else REGIONS[1] if z <= 1.2 else REGIONS[2] if z > 1.6 else None
    if row["secondary_structure"] not in {"HELIX", "SHEET"} and (
        not compatible or value == REGIONS[2]
    ):
        value = None
    return value


def make_matches(rows, errors, width, compatible):
    """Reconstruct all cells from covariates; fsum arithmetic for means/effects."""
    strata = defaultdict(lambda: defaultdict(list))
    confidence = {}
    for r in rows:
        confidence[r["canonical_position"]] = r["confidence"]
        reg = region(r, compatible)
        if reg is None:
            continue
        bin_id = min(int(round(r["confidence"] * 100, 6) // width), 100 // width - 1)
        category = "UNCONTROLLED" if compatible else r["secondary_structure"]
        strata[(bin_id, category)][reg].append(r["canonical_position"])
    result = {}
    for left, right in PAIRS:
        cells = []
        available = sum(len(cell[left]) + len(cell[right]) for cell in strata.values())
        for (bin_id, secondary), cell in sorted(strata.items()):
            x, y = sorted(cell[left]), sorted(cell[right])
            mass = min(len(x), len(y))
            if mass < 3:
                continue
            dc = math.fsum(confidence[p] for p in x) / len(x) - math.fsum(
                confidence[p] for p in y
            ) / len(y)
            de = math.fsum(errors[p] for p in x) / len(x) - math.fsum(errors[p] for p in y) / len(y)
            cells.append(
                {
                    "bin": bin_id,
                    "secondary_structure": secondary,
                    "left_positions": x,
                    "right_positions": y,
                    "mass": mass,
                    "confidence_difference": dc,
                    "error_difference": de,
                }
            )
        total = sum(c["mass"] for c in cells)
        included = total >= 20
        result[left + "-" + right] = {
            "cells": cells,
            "mass": total,
            "included": included,
            "available": available,
            "retained": sum(len(c["left_positions"]) + len(c["right_positions"]) for c in cells)
            if included
            else 0,
            "effective": 2 * total if included else 0,
            "effect": math.fsum(c["mass"] * c["error_difference"] for c in cells) / total
            if included
            else None,
            "confidence_difference": math.fsum(
                c["mass"] * c["confidence_difference"] for c in cells
            )
            / total
            if included
            else None,
        }
    return result


def bootstrap(effects, family, block, draws=20000, seed=20909):
    """Independent frequency-weighted block resampling, not indexed-value sums."""
    values = defaultdict(list)
    for acc in sorted(effects):
        values[family[acc]].append(effects[acc])
    means = {f: math.fsum(v) / len(v) for f, v in sorted(values.items())}
    members = defaultdict(list)
    for f in means:
        members[block[f]].append(f)
    ids = sorted(members)
    totals = np.asarray([math.fsum(means[f] for f in members[b]) for b in ids])
    denominators = np.asarray([len(members[b]) for b in ids])
    samples = np.random.Generator(np.random.PCG64(seed)).integers(0, len(ids), (draws, len(ids)))
    # Represent each draw by the number of occurrences of each entire dependence block.
    distribution = np.asarray(
        [
            float(np.bincount(row, minlength=len(ids)) @ totals)
            / float(np.bincount(row, minlength=len(ids)) @ denominators)
            for row in samples
        ]
    )
    ordered = np.sort(distribution)

    def percentile(q):
        location = q * (draws - 1)
        lo = math.floor(location)
        hi = math.ceil(location)
        return float(ordered[lo] + (ordered[hi] - ordered[lo]) * (location - lo))

    return {
        "estimate": math.fsum(means.values()) / len(means),
        "interval": [percentile(0.0125), percentile(0.9875)],
        "proteins": len(effects),
        "groups": len(means),
        "blocks": len(ids),
        "family_means": means,
        "quantiles": [0.0125, 0.9875],
        "draws": draws,
        "seed": seed,
    }, distribution


def material_direction(interval, margin=0.05):
    return 1 if interval[0] > margin else -1 if interval[1] < -margin else 0


def decision(contrasts):
    primary = contrasts["primary"]
    primaries = [a + "-" + b for a, b in PAIRS[:2]]
    if any(not primary[c]["useful"] for c in primaries):
        return "MEMBRANECAL_V0_2_INDEPENDENT_COHORT_NOT_FEASIBLE", []
    flags = []
    considered = [
        "bins_1",
        "bins_5",
        "family_novel",
        "unresolved_excluded",
        "broad_family",
        "curated_homology",
        "v01_compatible",
    ]
    for c in primaries:
        original = primary[c]
        direction = material_direction(original["interval"])
        for name in considered:
            s = contrasts[name][c]
            if not s["useful"]:
                continue
            shift = s["estimate"] - original["estimate"]
            if (abs(shift) >= 0.05 and s["estimate"] * original["estimate"] < 0) or (
                direction and s["estimate"] * direction < 0
            ):
                flags.append((c, name))
    if flags:
        return "MEMBRANECAL_V0_2_METHOD_SENSITIVITY_REMAINS_DOMINANT", flags
    if any(material_direction(primary[c]["interval"]) for c in primaries):
        return "MEMBRANECAL_V0_2_REGIONAL_EFFECT_SUPPORTED", []
    return "MEMBRANECAL_V0_2_NO_MATERIAL_REGIONAL_EFFECT_ESTABLISHED", []


def run(study, output):
    study = Path(study).resolve()
    output = Path(output)
    root = study / DATA

    from .verify import verify_hashes
    bindings = verify_hashes(study)
    receipt = read(root / "unseal_receipt.json")
    assert receipt["freeze_commit"] == FREEZE and receipt["protocol_sha256"] == digest(
        root / "protocol.json"
    )
    for name, expected in read(root / "score_manifest.json")["files"].items():
        assert digest(root / name) == expected
    protocol = read(root / "protocol.json")
    assert protocol["planning"]["material_difference"] == 0.05
    assert protocol["multiplicity"]["quantiles_primary"] == [0.0125, 0.9875]
    assert protocol["uncertainty"]["bootstrap_draws"] == 20000
    proteins = read(root / "proteins.json")["records"]
    development = read(root / "development_sequences.json")["records"]
    assert len(development) == 174
    old_aliases = {
        a
        for p in development
        for a in [p["accession"], p["current_primary_accession"], *p["secondary_accessions"]]
    }
    assert not old_aliases & {p["accession"] for p in proteins}
    assert not {p["entry_id"] for p in development} & {p["entry_id"] for p in proteins}
    old_hashes = {p["construct_sequence_sha256"] for p in development} | {
        hashlib.sha256(p["construct_sequence"].encode()).hexdigest() for p in development
    }
    assert not old_hashes & {p["construct_sequence_sha256"] for p in proteins}
    groups = read(root / "groups.json")
    all_records = {p["accession"]: p for p in development + proteins}
    edge_pairs = [(e["a"], e["b"]) for e in groups["edges"]]
    operational = components(all_records, edge_pairs, "family:")
    assert operational == groups["raw_component_by_accession"]
    incident = {a for pair in edge_pairs for a in pair}
    resolved = {a: bool(names(p)) or a in incident for a, p in all_records.items()}
    assert resolved == groups["resolved_by_accession"]
    family = {a: operational[a] if resolved[a] else "family:UNRESOLVED_POOLED" for a in all_records}
    assert family == groups["family_by_accession"]
    curated = defaultdict(list)
    for a, p in sorted(all_records.items()):
        for label in sorted(names(p)):
            curated[label].append(a)
    curated_edges = [(members[0], a) for members in curated.values() for a in members[1:]]
    assert {tuple(p) for p in groups["curated_dependence_edges"]} == set(curated_edges)
    homology = components(all_records, edge_pairs + curated_edges, "homology:")
    assert homology == groups["curated_homology_by_accession"]
    old_homology = {homology[p["accession"]] for p in development}
    novel = {
        p["accession"]
        for p in proteins
        if resolved[p["accession"]] and homology[p["accession"]] not in old_homology
    }
    assert novel == set(groups["family_novel_accessions"])
    entry = defaultdict(list)
    for p in proteins:
        entry[p["entry_id"]].append(family[p["accession"]])
    block_edges = [(family[a], family[b]) for a, b in curated_edges] + [
        (fs[0], f) for fs in entry.values() for f in fs[1:]
    ]
    blocks = components(set(family.values()), block_edges, "block:")
    assert blocks == groups["block_by_family"]
    rows = defaultdict(list)
    for r in pq.read_table(root / "residues_preoutcome.parquet").to_pylist():
        rows[r["accession"]].append(r)
        assert (region(r, False) or "UNKNOWN_OR_EXCLUDED") == r["region"]
    scores = {
        (r["accession"], r["canonical_position"]): r["local_accuracy"]
        for r in pq.read_table(root / "scores.parquet").to_pylist()
    }
    assert all(s is not None and 0 <= s <= 1 for s in scores.values())
    matches = read(root / "matched_sets.json")
    feasibility = read(root / "feasibility.json")
    reference = read(study / RESULT)
    scalar_checks = read(root / "metric_validation.json")["checks"]
    raw_checks = read(
        study / "provenance/source_investigation.json"
    )["records"]
    assert len(scalar_checks) == 5 and sum(s["residues"] for s in scalar_checks) == 2341
    assert all(s["matrix_scalar_exact"] for s in scalar_checks)
    assert len(raw_checks) == 13
    assert sum(s["independent_column_coordinate_check"]["residues"] for s in raw_checks) == 4716
    assert all(s["independent_column_coordinate_check"]["passed"] for s in raw_checks)
    supported = {}
    all_effects = {}
    rebuilt_by_scenario = {}
    for scenario, width, compatible in [
        ("primary", 2, False),
        ("bins_1", 1, False),
        ("bins_5", 5, False),
        ("v01_compatible", 2, True),
    ]:
        expected = {(m["accession"], m["contrast"]): m for m in matches[scenario]}
        sums = defaultdict(Counter)
        imbalances = defaultdict(dict)
        effects = defaultdict(dict)
        included_pairs = set()
        for acc, rs in rows.items():
            errors = {
                r["canonical_position"]: scores[(acc, r["canonical_position"])] - r["confidence"]
                for r in rs
            }
            rebuilt = make_matches(rs, errors, width, compatible)
            for contrast, m in rebuilt.items():
                sums[contrast].update({k: m[k] for k in ["available", "retained", "effective"]})
                if not m["included"]:
                    assert (acc, contrast) not in expected
                    continue
                included_pairs.add((acc, contrast))
                frozen = expected[(acc, contrast)]
                assert m["mass"] == frozen["mass"] and len(m["cells"]) == len(frozen["cells"])
                for x, y in zip(m["cells"], frozen["cells"], strict=True):
                    for key in [
                        "bin",
                        "secondary_structure",
                        "left_positions",
                        "right_positions",
                        "mass",
                    ]:
                        assert x[key] == y[key]
                    close(x["confidence_difference"], y["confidence_difference"])
                close(m["confidence_difference"], frozen["confidence_difference"])
                close(
                    m["effect"], reference["contrasts"][scenario][contrast]["protein_effects"][acc]
                )
                effects[contrast][acc] = m["effect"]
                imbalances[contrast][acc] = m["confidence_difference"]
        assert included_pairs == set(expected)
        supported[scenario] = {}
        for contrast, s in sums.items():
            raw = s["retained"] / s["available"]
            effective = s["effective"] / s["available"]
            close(raw, feasibility[scenario][contrast]["raw_cell_retention"])
            close(effective, feasibility[scenario][contrast]["effective_mass_retention"])
            bygroup = defaultdict(list)
            for acc, value in imbalances[contrast].items():
                bygroup[family[acc]].append(value)
            balance = math.fsum(math.fsum(v) / len(v) for v in bygroup.values()) / len(bygroup)
            supported[scenario][contrast] = {
                **dict(s),
                "raw_retention": raw,
                "effective_retention": effective,
                "paired_proteins": len(effects[contrast]),
                "group_balanced_confidence_difference": balance,
                "maximum_absolute_protein_confidence_difference": max(
                    map(abs, imbalances[contrast].values())
                ),
            }
        all_effects[scenario] = effects
        rebuilt_by_scenario[scenario] = len(included_pairs)
    independent = {}
    bootstrap_arrays = {}
    for i, (left, right) in enumerate(PAIRS[:2]):
        name = left + "-" + right
        calculated, distribution = bootstrap(all_effects["primary"][name], family, blocks)
        expected = reference["contrasts"]["primary"][name]
        close(calculated["estimate"], expected["estimate"])
        for x, y in zip(calculated["interval"], expected["interval"], strict=True):
            close(x, y)
        assert (
            calculated["proteins"] == expected["proteins"]
            and calculated["groups"] == expected["families"]
            and calculated["blocks"] == expected["blocks"]
        )
        for f, v in calculated["family_means"].items():
            close(v, expected["family_means"][f])
        calculated["matches_reference_within_absolute_tolerance"] = TOLERANCE
        independent[name] = calculated
        bootstrap_arrays[f"contrast_{i + 1}"] = distribution
    expected_scenarios = {
        "primary",
        "bins_1",
        "bins_5",
        "v01_compatible",
        "family_novel",
        "exact_construct",
        "unresolved_excluded",
        "broad_family",
        "strict_sequence_groups",
        "identity_100",
        "curated_homology",
        "post_2018_release",
    }
    assert set(reference["contrasts"]) == expected_scenarios
    outcome, flags = decision(reference["contrasts"])
    assert outcome == reference["decision"] and flags == reference["sensitivity_flags"] == []
    # Independently recompute all scenario point estimates, including equal-curated weighting.
    sensitivity_checks = {}
    accepted = {p["accession"] for p in proteins}
    permitted = {
        "family_novel": novel,
        "exact_construct": {
            p["accession"] for p in proteins if p["canonical_sequence"] == p["construct_sequence"]
        },
        "unresolved_excluded": {a for a in accepted if resolved[a]},
        "identity_100": {p["accession"] for p in proteins if p["sequence_identity"] == 1},
        "post_2018_release": {
            p["accession"]
            for p in proteins
            if p["release_date"] and p["release_date"] > "2018-04-30"
        },
    }
    assert permitted["exact_construct"] == {
        p["accession"] for p in proteins if p["exact_construct_prediction"]
    }
    for scenario, contrasts in reference["contrasts"].items():
        if scenario == "curated_homology":
            mapping = {
                a: homology[a] if resolved[a] else "homology:UNRESOLVED_POOLED" for a in all_records
            }
        elif scenario == "broad_family":
            mapping = {
                a: groups["broad_family_by_accession"][a]
                if resolved[a]
                else "broad:UNRESOLVED_POOLED"
                for a in all_records
            }
        elif scenario == "strict_sequence_groups":
            mapping = {
                a: groups["strict_family_by_accession"][a]
                if resolved[a]
                else "strict:UNRESOLVED_POOLED"
                for a in all_records
            }
        else:
            mapping = family
        scenario_entries = defaultdict(list)
        for p in proteins:
            scenario_entries[p["entry_id"]].append(mapping[p["accession"]])
        scenario_edges = [(mapping[a], mapping[b]) for a, b in curated_edges] + [
            (fs[0], f) for fs in scenario_entries.values() for f in fs[1:]
        ]
        scenario_blocks = components(set(mapping.values()), scenario_edges, "block:")
        for name, s in contrasts.items():
            group_values = defaultdict(list)
            base = scenario if scenario in all_effects else "primary"
            assert set(s["protein_effects"]) == set(all_effects[base][name]) & permitted.get(
                scenario, accepted
            )
            for acc, value in s["protein_effects"].items():
                close(value, all_effects[base][name][acc])
                group_values[mapping[acc]].append(all_effects[base][name][acc])
            estimate = math.fsum(math.fsum(v) / len(v) for v in group_values.values()) / len(
                group_values
            )
            close(estimate, s["estimate"])
            assert len(group_values) == s["families"]
            resolved_groups = {mapping[a] for a in s["protein_effects"] if resolved[a]}
            assert len(resolved_groups) == s["resolved_families"]
            assert len({scenario_blocks[f] for f in resolved_groups}) == s["resolved_blocks"]
            assert len(s["protein_effects"]) == s["proteins"]
            assert len({scenario_blocks[f] for f in group_values}) == s["blocks"]
            assert s["useful"] == (
                s["proteins"] >= 30 and s["resolved_families"] >= 20 and s["resolved_blocks"] >= 20
            )
        sensitivity_checks[scenario] = (
            "all three point estimates independently reproduced; published intervals retained except independent primary resampling"
        )
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / "independent_bootstrap_draws.npz", **bootstrap_arrays)
    audit = {
        "source_commit": SOURCE,
        "freeze_commit": FREEZE,
        "bound_artifact_count": len(bindings),
        "source_bindings": bindings,
        "development_proteins_checked": len(development),
        "overlap": {"canonical_and_aliases": 0, "entries": 0, "constructs": 0},
        "cohort": reference["cohort"],
        "independence_of_implementation": "No imports of frozen scientific functions; graph traversal, matching reconstruction, fsum estimates, frequency-weighted PCG64 block resampling and manual linear quantiles implemented separately. Same specified RNG stream used for numerical comparability.",
        "primary": independent,
        "support": supported,
        "rebuilt_matched_protein_contrasts": rebuilt_by_scenario,
        "sensitivity_checks": sensitivity_checks,
        "decision": outcome,
        "sensitivity_flags": flags,
        "genuine_numerical_discrepancies": [],
        "tolerance": TOLERANCE,
        "metric_verification_scope": {
            "all_protein_compact_replay": "existing study verifier, run separately; not independent raw verification",
            "independent_scalar": "existing five proteins / 2341 residues",
            "independent_raw_column": "existing 13 proteins / 4716 residues; no all-554 raw verification claimed",
        },
        "grouping_levels_reproduced": [
            "operational components",
            "curated homology components",
            "dependence blocks",
        ],
        "bootstrap_draws_sha256": digest(output / "independent_bootstrap_draws.npz"),
    }
    audit["group_reconstruction_scope"] = (
        "Components independently reconstructed from frozen sequence edges; curated annotation bridges rebuilt. Sequence alignments/edge completeness and broad/strict mappings not independently rederived."
    )
    write(output / "independent_reproduction.json", audit)
    print(
        "Independent estimates, primary resampling, matching, graph levels and decision reproduced.",
        flush=True,
    )
    return audit
