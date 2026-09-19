"""Outcome-blind grouping and matching; deterministic numerical analysis utilities."""

from __future__ import annotations

import math
import re
from collections import defaultdict

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices

REGIONS = ("TM_CORE", "MEMBRANE_INTERFACE", "STRUCTURED_EXTRAMEMBRANE")
CONTRASTS = ((REGIONS[1], REGIONS[0]), (REGIONS[2], REGIONS[0]), (REGIONS[2], REGIONS[1]))


def component_ids(nodes, edges, prefix="family:"):
    parent = {n: n for n in sorted(nodes)}

    def root(n):
        while parent[n] != n:
            parent[n] = parent[parent[n]]
            n = parent[n]
        return n

    for a, b in sorted(edges):
        a, b = sorted((root(a), root(b)))
        parent[b] = a
    return {n: prefix + root(n) for n in sorted(nodes)}


def family_names(record):
    result = set()
    for phrase in record.get("uniprot_family_phrases", []):
        for name in re.findall(r"Belongs to (.*? family)(?:\.|;|$)", phrase):
            result.add(name.lower())
    # Development records already contain provenance-backed family sentences.
    result.update(
        t.removeprefix("uniprot:").lower()
        for t in record.get("family_tokens", [])
        if t.startswith("uniprot:")
    )
    return {re.sub(r"^(?:the|a|an)\s+", "", " ".join(name.split())) for name in result}


def pfam_complement(record):
    return tuple(
        sorted(
            record.get(
                "pfam_ids",
                [
                    t.removeprefix("pfam:")
                    for t in record.get("family_tokens", [])
                    if t.startswith("pfam:")
                ],
            )
        )
    )


def sequence_evidence(a, b):
    """Actual identity; paired-residue coverage, not alignment span including long gaps."""
    if min(len(a), len(b)) / max(len(a), len(b)) < 0.5:
        return {"identity": 0.0, "coverage_a": 0.0, "coverage_b": 0.0, "paired": 0}
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(
        a.replace("U", "X").replace("O", "X"), b.replace("U", "X").replace("O", "X")
    )[0]
    paired = identical = 0
    for (lo_a, hi_a), (lo_b, hi_b) in zip(*alignment.aligned, strict=True):
        paired += int(hi_a - lo_a)
        identical += sum(x == y for x, y in zip(a[lo_a:hi_a], b[lo_b:hi_b], strict=True))
    return {
        "identity": identical / alignment.length if alignment.length else 0.0,
        "coverage_a": paired / len(a),
        "coverage_b": paired / len(b),
        "paired": paired,
    }


def grouping_edge(a, b, evidence):
    coverage = min(evidence["coverage_a"], evidence["coverage_b"])
    if evidence["identity"] >= 0.3 and coverage >= 0.8:
        return "FULL_LENGTH_SEQUENCE"
    # Shared Pfam alone never creates an edge. This conservative remote-family rescue
    # requires curated whole-family evidence, identical complete ID complement AND
    # a substantial sequence alignment. ID complement is not called domain architecture.
    if (
        evidence["identity"] >= 0.3
        and coverage >= 0.5
        and family_names(a) & family_names(b)
        and pfam_complement(a)
        and pfam_complement(a) == pfam_complement(b)
    ):
        return "CURATED_FAMILY_PLUS_COMPLEMENT_PLUS_SEQUENCE"
    return None


def assign_region(ratio, secondary, compatible=False):
    if not math.isfinite(ratio) or ratio < 0:
        return "UNKNOWN_OR_EXCLUDED"
    region = (
        REGIONS[0]
        if ratio <= 0.8
        else REGIONS[1]
        if ratio <= 1.2
        else REGIONS[2]
        if ratio > 1.6
        else "UNKNOWN_OR_EXCLUDED"
    )
    if secondary not in {"HELIX", "SHEET"} and (not compatible or region == REGIONS[2]):
        return "UNKNOWN_OR_EXCLUDED"
    return region


def match_protein(record, width=2, min_cell=3, min_mass=20, compatible=False):
    cells = defaultdict(lambda: defaultdict(list))
    for row in record["residues"]:
        region = assign_region(row["depth_ratio"], row["secondary_structure"], compatible)
        if region not in REGIONS:
            continue
        confidence = round(row["confidence"] * 100, 6)
        bin_id = min(int(confidence // width), int(100 / width) - 1)
        secondary = "UNCONTROLLED" if compatible else row["secondary_structure"]
        cells[(bin_id, secondary)][region].append(row)
    matches, diagnostics = [], []
    for left, right in CONTRASTS:
        kept = []
        for (bin_id, secondary), regions in sorted(cells.items()):
            x, y = regions[left], regions[right]
            if min(len(x), len(y)) < min_cell:
                continue
            kept.append(
                {
                    "bin": bin_id,
                    "secondary_structure": secondary,
                    "left_positions": sorted(r["canonical_position"] for r in x),
                    "right_positions": sorted(r["canonical_position"] for r in y),
                    "mass": min(len(x), len(y)),
                    "confidence_difference": float(
                        np.mean([r["confidence"] for r in x])
                        - np.mean([r["confidence"] for r in y])
                    ),
                }
            )
        mass = sum(c["mass"] for c in kept)
        nleft = sum(len(r[left]) for r in cells.values())
        nright = sum(len(r[right]) for r in cells.values())
        retained = (
            sum(len(c["left_positions"]) + len(c["right_positions"]) for c in kept)
            if mass >= min_mass
            else 0
        )
        diagnostic = {
            "accession": record["accession"],
            "contrast": left + "-" + right,
            "available_left": nleft,
            "available_right": nright,
            "eligible_mass": mass,
            "included": mass >= min_mass,
            "retained_residues": retained,
            "raw_cell_support_fraction": retained / (nleft + nright) if nleft + nright else 0.0,
            "effective_support_fraction": 2 * mass / (nleft + nright)
            if nleft + nright and mass >= min_mass
            else 0.0,
            "retained_left": sum(len(c["left_positions"]) for c in kept) if mass >= min_mass else 0,
            "retained_right": sum(len(c["right_positions"]) for c in kept)
            if mass >= min_mass
            else 0,
        }
        diagnostics.append(diagnostic)
        if mass >= min_mass:
            matches.append(
                {
                    "accession": record["accession"],
                    "contrast": left + "-" + right,
                    "mass": mass,
                    "confidence_difference": sum(
                        c["mass"] * c["confidence_difference"] for c in kept
                    )
                    / mass,
                    "cells": kept,
                }
            )
    return matches, diagnostics


def cluster_interval(effects, family_by, block_by, draws=20000, seed=20909, primary=True):
    """Equal proteins per family, equal families; resample family+entry blocks intact."""
    groups = defaultdict(list)
    for accession, value in sorted(effects.items()):
        groups[family_by[accession]].append(value)
    means = {f: float(np.mean(values)) for f, values in sorted(groups.items())}
    if not means:
        return {"estimate": None, "interval": None, "proteins": 0, "families": 0, "blocks": 0}
    blocks = defaultdict(list)
    for f, value in means.items():
        blocks[block_by[f]].append(value)
    packed = [blocks[k] for k in sorted(blocks)]
    total = np.array([sum(v) for v in packed])
    counts = np.array([len(v) for v in packed])
    rng = np.random.default_rng(seed)
    sampled = rng.integers(0, len(packed), (draws, len(packed)))
    bootstrap = total[sampled].sum(axis=1) / counts[sampled].sum(axis=1)
    alpha = 0.025 if primary else 0.05
    interval = (
        np.quantile(bootstrap, [alpha / 2, 1 - alpha / 2]).tolist() if len(packed) >= 2 else None
    )
    return {
        "estimate": float(np.mean(list(means.values()))),
        "interval": interval,
        "confidence_level": 1 - alpha,
        "proteins": len(effects),
        "families": len(means),
        "blocks": len(packed),
        "draws": draws,
        "seed": seed,
        "family_means": means,
        "protein_balanced": float(np.mean(list(effects.values()))),
        "leave_one_family_out": {
            f: float(np.mean([v for other, v in means.items() if other != f]))
            for f in means
            if len(means) > 1
        },
    }
