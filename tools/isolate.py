"""Fresh-environment Python audit guard for extracted-checkout replay."""

from __future__ import annotations

import argparse
import json
import os
import runpy
import sys
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("root", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument("--mode", choices=["compact", "full"], default="full")
args = parser.parse_args()
root = args.root.resolve()
output = args.output.resolve()
if not sys.flags.isolated or sys.prefix == sys.base_prefix or os.environ.get("PYTHONPATH"):
    raise RuntimeError("Use fresh-venv python -I and unset PYTHONPATH")
if output.is_relative_to(root) or root.is_relative_to(output):
    raise ValueError("Separate output required")
sys.dont_write_bytecode = True
sys.path.insert(0, str(root / "tools"))
allowed = [root, output, Path(sys.prefix).resolve(), Path(sys.base_prefix).resolve()]
if os.name == "nt":
    allowed.append(Path(os.environ["WINDIR"]) / "Fonts")
else:
    allowed.extend([Path("/usr/share/fonts"), Path("/usr/local/share/fonts"), Path("/etc/fonts")])
counts = {"allowed_python_file_opens": 0, "denied_parent_probe": 0}


def guard(event, arguments):
    if event in {"socket.connect", "subprocess.Popen", "os.system"}:
        raise RuntimeError("Offline replay prohibits network/subprocess rescue")
    if event == "open" and isinstance(arguments[0], (str, bytes, os.PathLike)):
        p = Path(os.fsdecode(arguments[0])).resolve()
        if not any(p.is_relative_to(base) for base in allowed):
            raise PermissionError("Python file access outside approved replay boundaries")
        counts["allowed_python_file_opens"] += 1


sys.addaudithook(guard)
try:
    (root.parent / "forbidden-parent-probe").read_bytes()
except PermissionError:
    counts["denied_parent_probe"] += 1
else:
    raise RuntimeError("Filesystem guard failed")
sys.argv = [str(root / "tools/replay.py"), "--mode", args.mode, "--output", str(output)]
runpy.run_path(str(root / "tools/replay.py"), run_name="__main__")
report = {
    "status": "PASS",
    "mode": args.mode,
    "isolated_python": True,
    "pythonpath_empty": True,
    "guard_scope": "Python audit-hook file/network/subprocess checks; not an OS sandbox",
    **counts,
}
(output / "isolation_report.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report))
