"""Offline public numerical replay. Run from the extracted package directory."""

# ruff: noqa: E501
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def hashes():
    manifest = json.loads((ROOT / "MANIFEST.json").read_text())["files"]
    for rel, digest in manifest.items():
        p = ROOT / rel
        if not p.resolve().is_relative_to(ROOT) or p.is_symlink():
            raise ValueError("Unsafe path")
        if hashlib.sha256(p.read_bytes()).hexdigest() != digest:
            raise ValueError("Checksum mismatch: " + rel)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--hashes-only", action="store_true")
    a = parser.parse_args()
    if os.environ.get("PYTHONPATH"):
        raise RuntimeError("Unset PYTHONPATH")
    if sys.flags.optimize:
        raise RuntimeError("Assertions must remain enabled")
    manifest = hashes()
    if a.hashes_only:
        print("All payload hashes verified")
        return
    out = a.output.resolve()
    if out == ROOT or ROOT.is_relative_to(out) or out.is_relative_to(ROOT):
        raise ValueError("Output must be outside package")
    if out.exists() and any(out.iterdir()):
        raise ValueError("Output must be new or empty")
    out.mkdir(parents=True, exist_ok=True)
    for name in ["audit", "results", "tables", "figures"]:
        (out / name).mkdir()
    os.environ["MPLCONFIGDIR"] = str(out / "renderer-cache")
    core = ROOT / "code/numerical"
    sys.path.insert(0, str(core))
    from membranecal_release import analysis, audit, presentation
    from membranecal_release.verify import compare

    for module in (analysis, audit, presentation):
        assert Path(module.__file__).resolve().is_relative_to(core)
    analysis.verify_scores()
    analysis.analyze(out / "results")
    reference = audit.read(core / "results/results.json")
    compare(audit.read(out / "results/results.json"), reference)
    assert (out / "results/report.md").read_bytes() == (core / "results/report.md").read_bytes()
    independent = audit.run(core, out / "audit")
    r, primary, sensitivity, support, flow, cases = presentation.tables(core, out, independent)
    presentation.supplemental_tables(core, out, r)
    presentation.table(out / "tables/claim_evidence.csv", presentation.claims(r, independent))
    presentation.verify_tables(out, r)
    for p in sorted((core / "tables").glob("*.csv")):
        compare(
            list(csv.DictReader((out / "tables" / p.name).open())), list(csv.DictReader(p.open()))
        )
    # Final publication figures use frozen CSVs, whose original result tables were checked above.
    art = out / "publication_figures"
    (art / "tables").mkdir(parents=True)
    (art / "figures").mkdir()
    for p in (ROOT / "data/result_tables").glob("*.csv"):
        shutil.copyfile(p, art / "tables" / p.name)
    (art / "evidence").mkdir()
    spec = importlib.util.spec_from_file_location(
        "publication_art", ROOT / "code/professional_figures.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.ROOT = art
    module.main()
    from PIL import Image

    comparisons = []
    for expected in sorted((ROOT / "figures").glob("*.png")):
        with Image.open(expected) as e, Image.open(art / "figures" / expected.name) as actual:
            equal = e.size == actual.size and e.tobytes() == actual.tobytes()
            comparisons.append({"file": expected.name, "pixels_equal": equal})
    assert len(comparisons) == 5 and all(row["pixels_equal"] for row in comparisons), (
        "Figure replay differs in the recorded runtime; inspect renderer/fonts before claiming identity"
    )
    report = {
        "status": "PASS",
        "bound_files": len(manifest),
        "private_checkout_required": False,
        "raw_upstream_cache_present": False,
        "primary_and_sensitivity_results": "unchanged",
        "independent_points": 36,
        "primary_bootstrap_draws_each": 20000,
        "csv_tables": 9,
        "figures": comparisons,
        "figure_portability": "Pixel identity applies only to recorded runtime/fonts; numerical coordinates and source CSVs are invariant",
    }
    (out / "replay_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
