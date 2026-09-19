from pathlib import Path
import hashlib
import json
ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / "data"
REPORT = ROOT / "results"
def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))
def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")
def verify_freeze(commit):
    from .verify import verify_hashes
    if commit != "5deb8f9abc345d811e0e98680b6122b4046a0d0d":
        raise ValueError("freeze identity mismatch")
    verify_hashes(ROOT)
