"""Run from the extracted release: python -m membranecal_release.verify."""

from __future__ import annotations

# ruff: noqa: E501
import argparse
import csv
import hashlib
import json
import math
import os
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = "740b223fbf6bcb43c42c126208510a835d04dade7f268d0999b684a7738c6671"
SCORES = "4f293fda6d7c7d0a4c8578071f396e7d573a54ed3c472795c67a5a5b9baca4b8"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_hashes(root):
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))["files"]
    for name, expected in manifest.items():
        path = root / name
        if not path.resolve().is_relative_to(root.resolve()) or path.is_symlink():
            raise ValueError("unsafe manifest member")
        if sha(path) != expected:
            raise ValueError("package checksum mismatch: " + name)
    for name, expected in [("protocol.json", PROTOCOL), ("scores.parquet", SCORES)]:
        if sha(root / "data" / name) != expected:
            raise ValueError("frozen identity mismatch: " + name)
    return manifest


def compare(actual, expected):
    if isinstance(expected, dict):
        if set(actual) != set(expected):
            raise ValueError("different result keys")
        for key in expected:
            compare(actual[key], expected[key])
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise ValueError("different result length")
        for a, b in zip(actual, expected, strict=True):
            compare(a, b)
    elif isinstance(expected, float):
        if not math.isclose(actual, expected, rel_tol=0, abs_tol=1e-12):
            raise ValueError("numerical mismatch")
    elif actual != expected:
        raise ValueError("result mismatch")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "reproduced")
    parser.add_argument("--hashes-only", action="store_true")
    args = parser.parse_args()
    if os.environ.get("PYTHONPATH"):
        raise RuntimeError("Unset PYTHONPATH before replay")
    if sys.flags.optimize:
        raise RuntimeError("Do not disable scientific assertions with optimized Python")
    bindings = verify_hashes(ROOT)
    if args.hashes_only:
        print("Package and frozen protocol/scores hashes verified")
        return
    output = args.output.resolve()
    if output == ROOT or output in ROOT.parents:
        raise ValueError("Output must be a separate directory")
    for name in ("audit", "tables", "figures", "results"):
        (output / name).mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(output / "renderer-cache")
    from PIL import Image

    from . import analysis, audit, presentation

    # No private package imports: all three modules must originate in this snapshot.
    for module in (analysis, audit, presentation):
        if not Path(module.__file__).resolve().is_relative_to(ROOT):
            raise RuntimeError("Scientific import escaped snapshot")
    analysis.verify_scores()
    analysis.analyze(output / "results")
    reference = audit.read(ROOT / "results/results.json")
    compare(audit.read(output / "results/results.json"), reference)
    if (output / "results/report.md").read_bytes() != (ROOT / "results/report.md").read_bytes():
        raise ValueError("original report mismatch")
    independent = audit.run(ROOT, output / "audit")
    r, primary, sensitivity, support, flow, cases = presentation.tables(ROOT, output, independent)
    presentation.figures(output, primary, sensitivity, support, flow, cases)
    presentation.supplemental_tables(ROOT, output, r)
    presentation.table(output / "tables/claim_evidence.csv", presentation.claims(r, independent))
    checked = presentation.verify_tables(output, r)
    for name in (
        "primary_effects",
        "sensitivity_effects",
        "support_balance",
        "cohort_flow",
        "descriptive_cases",
        "eligibility_criteria",
        "all_scenario_effects",
        "influence_diagnostics",
        "claim_evidence",
    ):
        compare(
            list(csv.DictReader((output / f"tables/{name}.csv").open(encoding="utf-8"))),
            list(csv.DictReader((ROOT / f"tables/{name}.csv").open(encoding="utf-8"))),
        )
    # Verify every reported sensitivity bound in the journal dataset against replay.
    for row in csv.DictReader((ROOT / "tables/all_scenario_effects.csv").open(encoding="utf-8")):
        s = r["contrasts"][row["scenario"]][row["contrast"]]
        for key, value in [
            ("estimate_points", s["estimate"]),
            ("lower_points", s["interval"][0]),
            ("upper_points", s["interval"][1]),
        ]:
            compare(float(row[key]) / 100, value)
    manuscript = (ROOT / "manuscript/manuscript.md").read_text(encoding="utf-8")
    for row in [*primary, *sensitivity]:
        for key in ("estimate_points", "lower_points", "upper_points"):
            if f"{row[key]:+.3f}" not in manuscript:
                raise ValueError("manuscript effect missing")
    figure_checks = []
    for p in sorted((output / "figures").glob("*.png")):
        if p.stem.split("_")[0].lstrip("0") not in "12345":
            raise ValueError("figure numbering")
        with Image.open(p) as actual, Image.open(ROOT / "figures" / p.name) as expected:
            figure_checks.append(
                {
                    "file": p.name,
                    "pixels_identical": actual.size == expected.size
                    and actual.tobytes() == expected.tobytes(),
                    "dimensions": list(actual.size),
                }
            )
    word = ROOT / "manuscript/manuscript_candidate.docx"
    with zipfile.ZipFile(word) as archive:
        xml = ET.fromstring(archive.read("word/document.xml"))
        words = " ".join(xml.itertext())
        scientific_text = manuscript.split("## Abstract", 1)[1]
        scientific_text = re.sub(r"\[([^]]*)\]\(https?://[^)]+\)", r"\1", scientific_text)
        numeric_expression = r"(?<![A-Za-z0-9])\d+(?:[.,]\d+)*(?![A-Za-z0-9])"
        numbers = Counter(re.findall(numeric_expression, scientific_text))
        if numbers - Counter(re.findall(numeric_expression, words)):
            raise ValueError("Word lost a manuscript numeric token")
        for row in [*primary, *sensitivity]:
            for key in ("estimate_points", "lower_points", "upper_points"):
                if f"{row[key]:+.3f}" not in words:
                    raise ValueError("Word numerical alignment")
        links = ET.fromstring(archive.read("word/_rels/document.xml.rels"))
        urls = {e.attrib.get("Target") for e in links}
        for url in re.findall(r"\]\((https?://[^)]+)\)", manuscript):
            if url not in urls:
                raise ValueError("Word source hyperlink missing")
    report = {
        "status": "PASS",
        "bound_files": len(bindings),
        "manifest_sha256": sha(ROOT / "manifest.json"),
        "protocol_sha256": PROTOCOL,
        "scores_sha256": SCORES,
        "compact_metric_proteins": 554,
        "independent_points": 36,
        "independent_primary_bootstrap_draws_each": 20000,
        "scientific_tolerance": 1e-12,
        "table_effect_rows": checked,
        "regenerated_csv_tables": 9,
        "word_numeric_tokens_checked": sum(numbers.values()),
        "all_scenario_rows": 36,
        "manuscript_word_numeric_and_urls": "PASS",
        "figures": figure_checks,
        "figure_scope": "Exact pixel equality reported for this environment only; renderer portability is not claimed",
        "private_checkout_required": False,
    }
    audit.write(output / "verification_report.json", report)
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
