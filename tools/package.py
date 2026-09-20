"""Refresh reviewed metadata or build a deterministic, allowlisted release ZIP."""

from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

from verify import ROOT, SELF_FILES, inventory, read, sha, validate, verify_files


def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n"
    )


def refresh():
    # Never refresh the scientific pins from current bytes.
    verify_files(ROOT, read(ROOT, "release/frozen-files.json")["files"])
    names = read(ROOT, "release/allowlist.json")
    actual = inventory(ROOT)
    if set(actual) - set(names):
        raise ValueError("Review the explicit allowlist before adding files")
    if set(names) - set(actual) - SELF_FILES - {"release/licence-ledger.json"}:
        raise ValueError("Missing allowlisted files")
    ledger = []
    for name in names:
        if name in {"schemas/citation-file-format-1.2.0.json", "schemas/CFF-LICENSE.txt"}:
            category = "third-party citation-file-format schema and licence"
            licence = "CC-BY-4.0; see schemas/CFF-LICENSE.txt and THIRD_PARTY_NOTICES.md"
        elif name.startswith("publication/v0.2/"):
            category = "mixed study payload; component-level rights retained"
            licence = "See publication/v0.2/licence_ledger.json and THIRD_PARTY_NOTICES.md"
        elif name.endswith(".py"):
            category, licence = "author-created code", "Apache-2.0"
        elif name.endswith(".md"):
            category, licence = "author-created documentation", "CC-BY-4.0"
        else:
            category, licence = "authored repository configuration and metadata", "Apache-2.0"
        ledger.append({"path": name, "category": category, "licence": licence})
    dump(ROOT / "release/licence-ledger.json", {"blanket_data_relicensing": False, "files": ledger})
    files = {n: sha((ROOT / n).read_bytes()) for n in names if n not in SELF_FILES}
    dump(ROOT / "release/manifest.json", {"schema_version": "1.0.0", "files": files})
    files["release/manifest.json"] = sha((ROOT / "release/manifest.json").read_bytes())
    (ROOT / "release/checksums.sha256").write_text(
        "".join(f"{files[n]}  {n}\n" for n in sorted(files)), encoding="utf-8", newline="\n"
    )


def build(output):
    validate()
    output = output.resolve()
    if output.is_relative_to(ROOT) or ROOT.is_relative_to(output):
        raise ValueError("Build output must be separate from the checkout")
    output.mkdir(parents=True, exist_ok=True)
    target = output / "membranecal-0.2.0.zip"
    if target.exists():
        raise ValueError("Do not overwrite an existing release archive")
    names = read(ROOT, "release/allowlist.json")
    with zipfile.ZipFile(target, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for n in names:
            info = zipfile.ZipInfo("membranecal-0.2.0/" + n, (2026, 9, 19, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            z.writestr(info, (ROOT / n).read_bytes(), compresslevel=9)
    receipt = {
        "archive": target.name,
        "sha256": sha(target.read_bytes()),
        "bytes": target.stat().st_size,
        "files": len(names),
        "manifest_sha256": sha((ROOT / "release/manifest.json").read_bytes()),
        "package_version": "0.2.0",
        "study_version": "0.2",
        "doi": "10.5281/zenodo.22843162",
        "doi_state": read(ROOT, "release/metadata.json")["doi_state"],
        "scientific_changes": False,
    }
    target.with_suffix(".zip.sha256").write_text(
        receipt["sha256"] + "  " + target.name + "\n", encoding="utf-8", newline="\n"
    )
    dump(output / "archive-identity.json", receipt)
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.refresh:
        refresh()
    if args.output:
        print(json.dumps(build(args.output)))
