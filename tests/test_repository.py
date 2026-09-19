import hashlib
import json
from pathlib import Path

import pytest

from tools.verify import safe_member, validate, verify_files

ROOT = Path(__file__).resolve().parents[1]


def test_repository_integrity():
    assert validate()["status"] == "PASS"


@pytest.mark.parametrize("name", ["../secret", "/absolute", "C:/secret", "a\\b", ".git/config"])
def test_rejects_unsafe_member(name):
    with pytest.raises(ValueError):
        safe_member(name)


def test_tampered_payload_is_rejected(tmp_path):
    (tmp_path / "sample").write_bytes(b"original")
    pins = {"sample": hashlib.sha256(b"original").hexdigest()}
    verify_files(tmp_path, pins)
    (tmp_path / "sample").write_bytes(b"changed")
    with pytest.raises(ValueError, match="Checksum mismatch"):
        verify_files(tmp_path, pins)


def test_missing_payload_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="Checksum mismatch"):
        verify_files(tmp_path, {"missing": "0" * 64})


def test_derivation_preserves_every_scientific_file():
    frozen = json.loads((ROOT / "release/frozen-files.json").read_text())["files"]
    entries = json.loads((ROOT / "release/derivation-map.json").read_text())["entries"]
    indexed = {entry["public_path"]: entry for entry in entries}
    for name, sha in frozen.items():
        assert indexed[name]["source_sha256"] == sha
        assert indexed[name]["public_sha256"] == sha


def test_main_results_match_frozen_tables():
    import csv

    path = ROOT / "publication/v0.2/data/result_tables/primary_effects.csv"
    rows = list(csv.DictReader(path.open()))
    assert len(rows) == 2
    observed = {
        tuple(
            round(float(row[key]), 3) for key in ("estimate_points", "lower_points", "upper_points")
        )
        for row in rows
    }
    assert observed == {(0.583, 0.234, 0.969), (0.065, -0.809, 0.940)}
