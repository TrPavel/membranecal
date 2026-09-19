"""Dispatch unchanged frozen numerical code; no new scientific computation policy."""

from __future__ import annotations

import argparse
import json
import os
import runpy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run(mode, output):
    from verify import validate

    if os.environ.get("PYTHONPATH") or sys.flags.optimize:
        raise RuntimeError("Unset PYTHONPATH and keep assertions enabled")
    validate()
    output = output.resolve()
    if output.is_relative_to(ROOT) or ROOT.is_relative_to(output):
        raise ValueError("Output must be separate from repository")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output must be new or empty")
    payload = ROOT / "publication/v0.2"
    if mode == "full":
        sys.argv = [str(payload / "replay.py"), "--output", str(output)]
        runpy.run_path(str(payload / "replay.py"), run_name="__main__")
        return
    core = payload / "code/numerical"
    sys.path.insert(0, str(core))
    from membranecal_release import analysis

    if not Path(analysis.__file__).resolve().is_relative_to(core):
        raise RuntimeError("Scientific import escaped payload")
    analysis.verify_scores()
    output.mkdir(parents=True, exist_ok=True)
    report = {
        "status": "PASS",
        "mode": "compact",
        "compact_score_verification": "all 554 frozen protein records",
        "full_bootstrap_tables_figures": "not run in compact mode",
    }
    (output / "replay_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["compact", "full"], default="compact")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.mode, args.output)
