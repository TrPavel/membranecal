import csv
import os
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from . import audit
CONTRASTS = [a + "-" + b for a, b in audit.PAIRS[:2]]
LABELS = ["Interface − core", "Structured extra − core"]
SCENARIOS = [
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
]
SCENARIO_LABELS = [
    "Primary",
    "1-point bins",
    "5-point bins",
    "v0.1-compatible regions",
    "Family-novel subset",
    "Exact full sequence",
    "Unresolved excluded",
    "Broad grouping",
    "Strict sequence grouping",
    "Mapped identity 100%",
    "Equal curated homology",
    "Released after 2018-04-30",
]
COLORS = ["#176B87", "#AC5129"]


def text(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(value.rstrip() + "\n", encoding="utf-8", newline="\n")

def table(path, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

def markdown(rows):
    keys = list(rows[0])
    return "\n".join(
        ["| " + " | ".join(keys) + " |", "| " + " | ".join("---" for _ in keys) + " |"]
        + [
            "| " + " | ".join(str(r[k]).replace("|", ";").replace("\n", " ") for k in keys) + " |"
            for r in rows
        ]
    )

def effect(s):
    return (
        f"{100 * s['estimate']:+.3f} [{100 * s['interval'][0]:+.3f}, {100 * s['interval'][1]:+.3f}]"
    )

def tables(study, output, verified):
    r = audit.read(study / audit.RESULT)
    primary = []
    sensitivity = []
    for scenario, label in zip(SCENARIOS, SCENARIO_LABELS, strict=True):
        for contrast in CONTRASTS:
            s = r["contrasts"][scenario][contrast]
            row = {
                "scenario": scenario,
                "label": label,
                "contrast": contrast,
                "estimate_points": s["estimate"] * 100,
                "lower_points": s["interval"][0] * 100,
                "upper_points": s["interval"][1] * 100,
                "marginal_confidence": s["confidence_level"],
                "proteins": s["proteins"],
                "total_groups": s["families"],
                "resolved_groups": s["resolved_families"],
                "blocks": s["blocks"],
                "resolved_blocks": s["resolved_blocks"],
                "useful": s["useful"],
            }
            sensitivity.append(row)
            if scenario == "primary":
                primary.append(row)
    support = [
        {
            "contrast": c,
            **s,
            "denominator_definition": "eligible HELIX/SHEET residues in the two regions before cell/protein support filtering; not all scored residues",
        }
        for c, s in verified["support"]["primary"].items()
        if c in CONTRASTS
    ]
    discovery = audit.read(study / audit.DATA / "discovery.json")
    exclusions = audit.read(study / audit.DATA / "exclusions.json")
    duplicates = sum(x.get("reason") == "V02_CANONICAL_DUPLICATE" for x in exclusions)
    # The acquisition audit explicitly distinguishes candidate entries from protein cases.
    assert duplicates == 68
    cohort = r["cohort"]
    flow = [
        {
            "stage": "OPM index records",
            "count": discovery["total_opm_records"],
            "unit": "index records",
            "note": "Frozen source index",
        },
        *[
            {
                "stage": k,
                "count": v,
                "unit": "index records",
                "note": "Sequential metadata prefilter exclusion",
            }
            for k, v in discovery["prefilter_exclusions"].items()
        ],
        {
            "stage": "Candidates",
            "count": cohort["candidate_entries"],
            "unit": "entries",
            "note": "Unique-code representatives and per-family rounds after prefilter; cap 800",
        },
        *[
            {
                "stage": k,
                "count": v,
                "unit": "candidate entries",
                "note": "Mutually exclusive preparation status",
            }
            for k, v in cohort["source_status"].items()
        ],
        {
            "stage": "Eligible protein cases before canonical selection",
            "count": cohort["accepted_proteins"] + duplicates,
            "unit": "protein cases",
            "note": "Multiple eligible chains can arise from an entry",
        },
        {
            "stage": "Duplicate canonical cases removed",
            "count": duplicates,
            "unit": "protein cases",
            "note": "Outcome-blind representative selection",
        },
        {
            "stage": "Accepted canonical proteins",
            "count": cohort["accepted_proteins"],
            "unit": "proteins",
            "note": f"{cohort['accepted_entries']} experimental entries; {cohort['residues']} mapped scored residues",
        },
        {
            "stage": "Development proteins screened",
            "count": 174,
            "unit": "proteins",
            "note": "Final canonical/alias, PDB and construct overlaps all zero",
        },
    ]
    source = audit.read(
        study / "provenance/source_investigation.json"
    )
    source_by = {x["accession"]: x for x in source["records"]}
    cases = []
    for case in r["failure_cases"][:8]:
        s = source_by[case["accession"]]
        cases.append(
            {
                "accession": case["accession"],
                "entry_id": case["entry_id"],
                "mean_error_points": 100 * case["mean_error"],
                "experimental_title": "; ".join(s["experimental_title"]),
                "observed_residues": s["observed_residues"],
                "canonical_length": s["canonical_length"],
                "construct_length": s["construct_length"],
                "mapped_identity": s["sequence_identity"],
                "exact_full_sequence": s["exact_construct"],
                "source_dois": "; ".join(x for x in s["publication_dois"] if x),
                "raw_pdb_sha256": s["raw_pdb_sha256"],
                "raw_afdb_sha256": s["raw_afdb_sha256"],
                "interpretation": "post-outcome descriptive all-observed-residue mean; no causal state attribution",
            }
        )
    for name, rows in [
        ("primary_effects", primary),
        ("sensitivity_effects", sensitivity),
        ("support_balance", support),
        ("cohort_flow", flow),
        ("descriptive_cases", cases),
    ]:
        table(output / "tables" / (name + ".csv"), rows)
    report = (study / audit.RESULT.parent / "report.md").read_text(encoding="utf-8")
    for contrast in CONTRASTS:
        s = r["contrasts"]["primary"][contrast]
        assert (
            f"| {contrast} | {s['estimate'] * 100:.3f} | [{s['interval'][0] * 100:.3f}, {s['interval'][1] * 100:.3f}] | {s['proteins']} | {s['resolved_families']} | {s['resolved_blocks']} | {s['useful']} |"
            in report
        )
    return r, primary, sensitivity, support, flow, cases

def save(fig, folder, name):
    for suffix in ("svg", "png", "tiff"):
        fig.savefig(
            folder / f"{name}.{suffix}",
            dpi=300,
            **(
                {"metadata": {"Date": None}}
                if suffix == "svg"
                else {"pil_kwargs": {"compression": "tiff_lzw"}}
                if suffix == "tiff"
                else {"metadata": {"Software": "MembraneCal v0.2 review; Matplotlib"}}
            ),
        )
        if suffix == "svg":
            path = folder / f"{name}.svg"
            text(
                path,
                "\n".join(line.rstrip() for line in path.read_text(encoding="utf-8").splitlines()),
            )
        elif suffix == "tiff":
            path = folder / f"{name}.tiff"
            with Image.open(path) as raster:
                rgb = raster.convert("RGB")
            rgb.save(path, compression="tiff_lzw", dpi=(300, 300))
    plt.close(fig)

def figures(output, primary, sensitivity, support, flow, cases):
    plt.rcParams.update(
        {
            "font.family": "Arial" if os.name == "nt" else "Liberation Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.hashsalt": "membranecal-v02-closure",
            "axes.titleweight": "bold",
            "axes.labelcolor": "#27323A",
            "text.color": "#27323A",
            "figure.facecolor": "white",
        }
    )
    folder = output / "figures"
    fig, ax = plt.subplots(figsize=(7.5, 6.8))
    ax.set(xlim=(0, 10), ylim=(0, 10))
    ax.axis("off")
    boxes = [
        (5, 9.2, "8,915 OPM index records"),
        (5, 7.35, "800 metadata-selected candidate entries"),
        (5, 5.5, "509 entries with eligible protein cases"),
        (5, 3.65, "622 eligible protein cases"),
        (5, 1.7, "554 canonical proteins · 473 PDB entries\n223,714 mapped, scored residues"),
    ]
    # Counts are checked against their source table before drawing explicit flow labels.
    assert [flow[0]["count"], flow[4]["count"]] == [8915, 800]
    for x, y, label in boxes:
        ax.text(
            x,
            y,
            label,
            ha="center",
            va="center",
            fontsize=10,
            bbox={"boxstyle": "round,pad=0.65", "fc": "#EDF5F7", "ec": COLORS[0]},
        )
    for (_, y, _), (_, next_y, _) in zip(boxes[:-1], boxes[1:], strict=True):
        ax.annotate(
            "",
            (5, next_y + 0.48),
            (5, y - 0.45),
            arrowprops={"arrowstyle": "->", "color": COLORS[0]},
        )
    ax.text(
        0.1,
        8.1,
        "Prefilter excludes 3,540 type/resolution availability;\n1,497 resolution/boundary; 65 development overlaps.\nThen representatives and cap select 800 from 3,813.",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1},
    )
    ax.text(
        3.0,
        6.4,
        "286 entry exclusions; 5 source/QC failures",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1},
    )
    ax.text(
        3.0,
        4.6,
        "Multiple eligible protein chains per entry",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1},
    )
    ax.text(
        3.0,
        2.65,
        "68 duplicate canonical cases removed",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1},
    )
    ax.text(
        5,
        0.25,
        "Final overlap with all 174 development proteins: 0 canonical/alias · 0 PDB · 0 construct\nNew cohort independence does not establish absence of model training exposure.",
        ha="center",
        fontsize=9,
    )
    ax.set_title("Figure 1  |  Cohort flow and development independence", loc="left", pad=18)
    save(fig, folder, "01_cohort_flow")

    fig, axes = plt.subplots(1, 2, figsize=(7.5, 4.0), gridspec_kw={"width_ratios": [1, 1.3]})
    for ax, bounds, title in zip(
        axes,
        [(-6, 6), (-1.25, 1.45)],
        ["Investigator-chosen ±5 points", "Detail of the same estimates"],
        strict=True,
    ):
        ax.axvline(0, color="#87929B", lw=1)
        for i, s in enumerate(primary):
            ax.errorbar(
                s["estimate_points"],
                1 - i,
                xerr=[
                    [s["estimate_points"] - s["lower_points"]],
                    [s["upper_points"] - s["estimate_points"]],
                ],
                fmt="o",
                color=COLORS[i],
                capsize=4,
                lw=2,
            )
        ax.set(
            xlim=bounds,
            ylim=(-0.55, 1.5),
            yticks=[1, 0],
            yticklabels=["Interface\n− core", "Structured extra\n− core"] if ax == axes[0] else [],
            xlabel="Agreement − confidence\ncontrast (points)",
            title=title,
        )
    for margin in (-5, 5):
        axes[0].axvline(margin, ls="--", color="#6C7277", lw=1)
    for i, s in enumerate(primary):
        axes[1].text(
            -1.17,
            0.68 - i,
            f"{s['estimate_points']:+.3f} [{s['lower_points']:+.3f}, {s['upper_points']:+.3f}]\nn={s['proteins']}; blocks={s['blocks']}",
            fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1},
        )
    fig.suptitle(
        "Figure 2  |  Primary paired effects with 97.5% marginal intervals",
        x=0.01,
        ha="left",
        fontweight="bold",
    )
    fig.tight_layout()
    save(fig, folder, "02_primary_effects")

    by = {s["contrast"]: s for s in support}
    ordered = [by[c] for c in CONTRASTS]
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 4.3))
    x = np.arange(2)
    axes[0].bar(
        x - 0.18,
        [100 * s["raw_retention"] for s in ordered],
        0.34,
        label="Residues in retained cells",
        color=COLORS[0],
    )
    axes[0].bar(
        x + 0.18,
        [100 * s["effective_retention"] for s in ordered],
        0.34,
        label="Effective common mass",
        color="#91BAC8",
    )
    axes[0].set(
        ylim=(0, 105),
        ylabel="% of pre-match eligible regional residues",
        xticks=x,
        xticklabels=["Interface–core", "Structured\nextra–core"],
    )
    axes[0].legend(fontsize=8, loc="upper right")
    for bars in axes[0].containers:
        axes[0].bar_label(bars, fmt="%.2f%%", fontsize=8, padding=3)
    for i, s in enumerate(ordered):
        axes[0].text(
            i,
            5,
            f"denominator\n{s['available']:,}",
            ha="center",
            fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1},
        )
    axes[1].axhline(0, lw=1, color="#9A9FA4")
    axes[1].bar(
        x - 0.18,
        [100 * s["group_balanced_confidence_difference"] for s in ordered],
        0.34,
        label="Signed group-balanced",
        color=COLORS[0],
    )
    axes[1].bar(
        x + 0.18,
        [100 * s["maximum_absolute_protein_confidence_difference"] for s in ordered],
        0.34,
        label="Maximum absolute protein",
        color=COLORS[1],
    )
    axes[1].set(
        ylim=(-0.22, 1.03),
        ylabel="Residual confidence difference\n(pLDDT points)",
        xticks=x,
        xticklabels=["Interface–core", "Structured\nextra–core"],
    )
    axes[1].legend(fontsize=8, loc="upper left")
    fig.suptitle(
        "Figure 3  |  Matched support and residual confidence balance",
        x=0.01,
        ha="left",
        fontweight="bold",
    )
    fig.tight_layout()
    save(fig, folder, "03_support_balance")

    fig, axes = plt.subplots(1, 2, figsize=(7.5, 7.1), sharey=True)
    for i, (ax, contrast) in enumerate(zip(axes, CONTRASTS, strict=True)):
        ax.axvline(0, color="#87929B", lw=1)
        for j, scenario in enumerate(SCENARIOS):
            s = next(
                s for s in sensitivity if s["scenario"] == scenario and s["contrast"] == contrast
            )
            ax.errorbar(
                s["estimate_points"],
                j,
                xerr=[
                    [s["estimate_points"] - s["lower_points"]],
                    [s["upper_points"] - s["estimate_points"]],
                ],
                fmt="o" if j else "D",
                color=COLORS[i],
                capsize=3,
            )
            ax.text(
                1.03,
                j,
                str(s["proteins"]),
                transform=ax.get_yaxis_transform(),
                va="center",
                ha="left",
                fontsize=8,
                clip_on=False,
            )
        ax.set(
            xlim=(-1.45, 2.5),
            ylim=(11.6, -0.7),
            yticks=range(12),
            yticklabels=SCENARIO_LABELS,
            xlabel="Contrast (points)",
            title=LABELS[i],
        )
        ax.text(1.03, -0.4, "n", transform=ax.get_yaxis_transform(), ha="left", fontsize=8)
    fig.suptitle(
        "Figure 4  |  Prespecified sensitivities (descriptive intervals)",
        x=0.01,
        ha="left",
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.008,
        "97.5% marginal intervals; no across-scenario guarantee. n = matched proteins.\nInvestigator-chosen ±5-point margin lies outside this zoom.",
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.025, 1, 1))
    save(fig, folder, "04_sensitivities")

    fig, ax = plt.subplots(figsize=(7.5, 5.7))
    labels = [f"{s['accession']} / {s['entry_id']}" for s in cases]
    vals = [s["mean_error_points"] for s in cases]
    ax.barh(range(8), vals, color=[COLORS[1] if v < 0 else COLORS[0] for v in vals])
    ax.axvline(0, lw=1, color="#87929B")
    ax.set(
        yticks=range(8),
        yticklabels=labels,
        ylim=(7.7, -0.7),
        xlim=(-45, 23),
        xlabel="All-observed-residue mean (agreement − confidence)\npercentage points",
    )
    for i, value in enumerate(vals):
        ax.text(
            value + (0.65 if value > 0 else -0.65),
            i,
            f"{value:+.2f}",
            ha="left" if value > 0 else "right",
            va="center",
            fontsize=9,
        )
    ax.set_title(
        "Figure 5  |  Eight largest absolute protein mean discrepancies", loc="left", pad=17
    )
    fig.text(
        0.5,
        0.01,
        "Post-outcome unmatched means; no inferential intervals or causal state attribution.\nNegative = confidence exceeds agreement under this endpoint.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.025, 1, 1))
    save(fig, folder, "05_descriptive_cases")

def verify_tables(output, reference):
    """Re-read written source tables so rounding, naming and serialization are checked."""
    checked = 0
    for name in ("primary_effects", "sensitivity_effects"):
        with (output / "tables" / (name + ".csv")).open(encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream):
                expected = reference["contrasts"][row["scenario"]][row["contrast"]]
                for field, value in [
                    ("estimate_points", expected["estimate"]),
                    ("lower_points", expected["interval"][0]),
                    ("upper_points", expected["interval"][1]),
                ]:
                    audit.close(float(row[field]) / 100, value)
                assert int(row["proteins"]) == expected["proteins"]
                checked += 1
    return checked

def claims(r, verified):
    primary = r["contrasts"]["primary"]
    common = "C-alpha local agreement minus pLDDT/100; cell-mass weights, equal proteins within equal operational groups"
    prefix = audit.RESULT.as_posix() + "#/"
    return [
        {
            "ID": "C1",
            "Claim": "Small nonzero interface shift",
            "Exact evidence": prefix + "contrasts/primary/" + CONTRASTS[0],
            "Population": "437 matched proteins; 260 resolved + 1 pooled group; 186 blocks",
            "Endpoint and weighting": common,
            "Uncertainty": effect(primary[CONTRASTS[0]]) + " points; 97.5% marginal",
            "Qualification": "Excludes zero; does not establish a >=5-point departure",
        },
        {
            "ID": "C2",
            "Claim": "No material regional departure established",
            "Exact evidence": prefix + "decision; protocol.json#/planning/material_difference",
            "Population": "Interface 437; structured extra 229; separately matched",
            "Endpoint and weighting": common,
            "Uncertainty": effect(primary[CONTRASTS[1]]) + " points for extra-core",
            "Qualification": "Neither entire adjusted interval beyond +/-5; no formal equivalence claim or universal absence",
        },
        {
            "ID": "C3",
            "Claim": "Protein-independent validation relative to v0.1",
            "Exact evidence": "development_sequences.json; proteins.json; audit/independent_reproduction.json#/overlap",
            "Population": "All 554 versus all 174 accepted development proteins",
            "Endpoint and weighting": "Identity/alias, experimental-entry and construct exclusion audit",
            "Uncertainty": "Zero observed overlaps under frozen identifiers",
            "Qualification": "Not proof of absence of AFDB training/template/homolog exposure; family-novel is an operational subset",
        },
        {
            "ID": "C4",
            "Claim": "Primary inference is on common support",
            "Exact evidence": "tables/support_balance.csv; matched_sets.json; feasibility.json",
            "Population": "Contrast-specific HELIX/SHEET support",
            "Endpoint and weighting": common,
            "Uncertainty": "Raw retention 79.31% / 38.46%; effective mass 45.27% / 23.81%",
            "Qualification": "Denominator excludes non-structured and depth-gap residues before matching; no extrapolation to discarded residues",
        },
        {
            "ID": "C5",
            "Claim": "Operational groups and dependence blocks serve distinct roles",
            "Exact evidence": "groups.json; audit/independent_reproduction.json#/group_reconstruction_scope",
            "Population": "Union of 174 development and 554 validation records for grouping",
            "Endpoint and weighting": "Equal operational groups in point estimate; whole dependence blocks resampled; curated-homology sensitivity changes weights",
            "Uncertainty": "20,000 PCG64 block draws; ratio of sampled group sums to counts",
            "Qualification": "Not proven evolutionary families; resolved includes annotation-or-edge; broad/strict mappings and sequence edge derivation not independently rerun",
        },
        {
            "ID": "C6",
            "Claim": "Exact sequence does not establish matching biological state",
            "Exact evidence": "tables/descriptive_cases.csv; proteins.json#/records; source_investigation.json",
            "Population": "198 whole-cohort exact-sequence cases; 152/75 matched",
            "Endpoint and weighting": "Canonical sequence equals construct; descriptive case means separate",
            "Uncertainty": "Exact-subset estimates in tables/sensitivity_effects.csv",
            "Qualification": "Assembly, ligand, membrane, conformation and observability can differ despite identical sequence",
        },
        {
            "ID": "C7",
            "Claim": "Large discrepancies motivate descriptive investigations",
            "Exact evidence": "tables/descriptive_cases.csv; provenance/source_investigation.json",
            "Population": "Eight largest absolute all-residue protein mean errors selected after outcomes",
            "Endpoint and weighting": "Unmatched equal-residue protein means",
            "Uncertainty": "Descriptive; no intervals; 13 raw-column checks total including five scalar cases",
            "Qualification": "Experimental titles establish context, not an AFDB alternative state or causal mechanism",
        },
        {
            "ID": "C8",
            "Claim": "Prospectively frozen analysis protocol",
            "Exact evidence": "provenance/original_PRE_OUTCOME_FREEZE.json; unseal_receipt.json; provenance/registration_search.json",
            "Population": "Existing v0.2 protocol, cohort, covariates, matching and scientific code",
            "Endpoint and weighting": "Git freeze " + audit.FREEZE,
            "Uncertainty": "Git/source hash evidence; public searches 2026-09-09 found no registration",
            "Qualification": "Local Git chronology is not independently timestamped external preregistration; private/unindexed records cannot be excluded",
        },
        {
            "ID": "C9",
            "Claim": "Independent statistical reproduction agrees",
            "Exact evidence": "audit/independent_reproduction.json; audit/independent_bootstrap_draws.npz",
            "Population": "Two primary contrasts; four rebuilt matching scenarios; 36 scenario point estimates",
            "Endpoint and weighting": "Separate implementation with same frozen inputs, specified RNG and estimand",
            "Uncertainty": "Absolute numerical tolerance 1e-12",
            "Qualification": "All-554 compact metric replay uses existing vectorized metric; scalar metric check only 5/2341 and stored raw-column evidence 13/4716",
        },
        {
            "ID": "C10",
            "Claim": "No frozen method-sensitivity flag fired",
            "Exact evidence": prefix + "sensitivity_flags; tables/sensitivity_effects.csv",
            "Population": "11 prespecified sensitivities; seven feed the material-decision rule",
            "Endpoint and weighting": "Scenario-specific support/weights",
            "Uncertainty": "Retained 97.5% marginal intervals for the two contrasts; no across-scenario adjustment",
            "Qualification": "Flag requires >=5-point shift with sign crossing or reversal of established material direction; absence of flag is not method interchangeability",
        },
        {
            "ID": "C11",
            "Claim": "Five-point materiality is investigator chosen",
            "Exact evidence": "protocol.json#/planning; provenance/original_design.md",
            "Population": "Frozen matched endpoint",
            "Endpoint and weighting": "Symmetric margin of 0.05 in unit-interval score",
            "Uncertainty": "Rough planning precision illustration only",
            "Qualification": "No validated biological, measurement-error or utility threshold documented; no retrospective rationale added",
        },
    ]

def supplemental_tables(study, output, results):
    protocol = audit.read(study / audit.DATA / "protocol.json")
    eligibility = [
        {"Frozen field": key, "Criterion": value} for key, value in protocol["cohort"].items()
    ]
    eligibility += [
        {"Frozen field": "mapping." + key, "Criterion": value}
        for key, value in protocol["mapping"].items()
    ]
    table(output / "tables/eligibility_criteria.csv", eligibility)
    influence = []
    for contrast, record in results["contrasts"]["primary"].items():
        for group, value in record["leave_one_family_out"].items():
            influence.append(
                {
                    "contrast": contrast,
                    "omitted_operational_group": group,
                    "estimate_points": value * 100,
                    "shift_points": (value - record["estimate"]) * 100,
                    "scope": "one operational group; not one dependence block",
                }
            )
    table(output / "tables/influence_diagnostics.csv", influence)
    all_effects = [
        {
            "scenario": scenario,
            "contrast": contrast,
            "estimate_points": record["estimate"] * 100,
            "lower_points": record["interval"][0] * 100,
            "upper_points": record["interval"][1] * 100,
            "marginal_confidence": record["confidence_level"],
            "proteins": record["proteins"],
            "resolved_groups": record["resolved_families"],
            "blocks": record["blocks"],
        }
        for scenario, contrasts in results["contrasts"].items()
        for contrast, record in contrasts.items()
    ]
    assert len(all_effects) == 36
    table(output / "tables/all_scenario_effects.csv", all_effects)
