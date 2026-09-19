"""Reject edits to baseline scientific pins or pinned files in a pull request."""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--base", required=True)
args = parser.parse_args()
baseline = subprocess.check_output(
    ["git", "show", args.base + ":release/frozen-files.json"], cwd=ROOT
)
if baseline != (ROOT / "release/frozen-files.json").read_bytes():
    raise ValueError("Frozen scientific pin set changed")
for name, digest in json.loads(baseline)["files"].items():
    if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
        raise ValueError("Frozen content changed: " + name)
print("Baseline scientific pins and contents unchanged")
