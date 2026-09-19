"""Verify the shipped repository without installing a private project package."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import zipfile
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
IGNORED = {".git", ".venv", "__pycache__", ".pytest_cache", ".ruff_cache", "local-evidence", "dist"}
SELF_FILES = {"release/manifest.json", "release/checksums.sha256"}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def read(root, name):
    return json.loads((root / name).read_text(encoding="utf-8"))


def inventory(root):
    return sorted(
        p.relative_to(root).as_posix()
        for p in root.rglob("*")
        if p.is_file() and not any(part in IGNORED for part in p.relative_to(root).parts)
    )


def safe_member(name):
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or "\\" in name or ":" in name:
        raise ValueError("Unsafe archive or manifest path")
    if any(part.lower() == ".git" for part in path.parts):
        raise ValueError("Git objects must not be distributed")


def verify_files(root, files):
    for name, digest in files.items():
        safe_member(name)
        p = root / name
        if p.is_symlink() or not p.resolve().is_relative_to(root.resolve()):
            raise ValueError("Unsafe filesystem member: " + name)
        if not p.is_file() or sha(p.read_bytes()) != digest:
            raise ValueError("Checksum mismatch: " + name)


def privacy(root, names, structured=False):
    patterns = {
        "unrelated_workstream": re.compile(rb"(?i)membrane" + b"states"),
        "private_repository": re.compile(rb"TrPavel/" + b"membranestate-bench"),
        "absolute_private_path": re.compile(
            rb"(?:[A-Z]:[/\\](?:Users|MembraneState|tmp)|/" + rb"home/|/" + rb"root/)"
        ),
        "credential": re.compile(
            rb"(?:\bgh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{35,}|"
            rb"\bAKIA[A-Z0-9]{16}\b|\brpa_[A-Za-z0-9]{20,}\b|"
            rb"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----)"
        ),
    }
    findings = []
    count = 0

    def check(name, data, depth=0):
        nonlocal count
        count += 1
        safe_member(name.split("!")[-1])
        if depth > 8:
            raise ValueError("Archive nesting limit exceeded")
        if data.startswith(b"PK"):
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                for info in z.infolist():
                    safe_member(info.filename)
                    if (
                        info.file_size > 100_000_000
                        or (info.external_attr >> 16) & 0o170000 == 0o120000
                    ):
                        raise ValueError("Oversized or symlink archive member")
                    check(name + "!" + info.filename, z.read(info), depth + 1)
            return
        for kind, pattern in patterns.items():
            if pattern.search(data):
                findings.append({"path": name, "kind": kind})
        if structured and name.endswith(".parquet"):
            import pyarrow.parquet as pq

            table = pq.read_table(io.BytesIO(data))
            for col in table.column_names:
                if "string" in str(table.schema.field(col).type):
                    check(name + "!column", json.dumps(table[col].to_pylist()).encode(), depth + 1)
        if structured and data.startswith(b"%PDF"):
            from pypdf import PdfReader

            pdf = PdfReader(io.BytesIO(data))
            text = str(pdf.metadata) + "\n" + "\n".join(p.extract_text() or "" for p in pdf.pages)
            check(name + "!text", text.encode(), depth + 1)
        if structured and name.lower().endswith((".png", ".tiff", ".jpg")):
            from PIL import Image

            with Image.open(io.BytesIO(data)) as im:
                check(name + "!metadata", str(im.info).encode(), depth + 1)

    for name in names:
        check(name, (root / name).read_bytes())
    if findings:
        raise ValueError(json.dumps(findings))
    return {"status": "PASS", "components": count, "structured": structured, "findings": []}


def validate(root=ROOT, development=False):
    names = inventory(root)
    allow = read(root, "release/allowlist.json")
    if names != allow:
        raise ValueError("Unexpected or missing files: " + str(sorted(set(names) ^ set(allow))))
    manifest = read(root, "release/manifest.json")["files"]
    if set(manifest) != set(names) - SELF_FILES:
        raise ValueError("Repository manifest coverage mismatch")
    verify_files(root, manifest)
    checks = {}
    for line in (root / "release/checksums.sha256").read_text().splitlines():
        digest, name = line.split("  ", 1)
        checks[name] = digest
    if set(checks) != set(names) - {"release/checksums.sha256"}:
        raise ValueError("Checksum coverage mismatch")
    verify_files(root, checks)
    frozen = read(root, "release/frozen-files.json")["files"]
    verify_files(root, frozen)
    payload = root / "publication/v0.2"
    verify_files(payload, read(payload, "MANIFEST.json")["files"])
    entries = read(root, "release/derivation-map.json")["entries"]
    if len(entries) != 141 or len({e["public_path"] for e in entries}) != 141:
        raise ValueError("Imported derivation coverage mismatch")
    for entry in entries:
        verify_files(root, {entry["public_path"]: entry["public_sha256"]})
        if entry["transformation_type"] == "unchanged_copy":
            if entry["source_sha256"] != entry["public_sha256"]:
                raise ValueError("Unchanged copy was altered")
    ledger = read(root, "release/licence-ledger.json")["files"]
    if {r["path"] for r in ledger} != set(names):
        raise ValueError("Licence coverage mismatch")
    meta = read(root, "release/metadata.json")
    cff = (root / "CITATION.cff").read_text()
    if meta["package_version"] != "0.2.0" or meta["study_version"] != "0.2":
        raise ValueError("Version mismatch")
    if meta["doi"] != "10.5281/zenodo.22843162" or meta["doi_state"] != "reserved_draft":
        raise ValueError("DOI state mismatch")
    if "version: 0.2.0" not in cff or "doi: " + meta["doi"] not in cff:
        raise ValueError("Citation metadata mismatch")
    if (payload / "CITATION.cff").read_bytes() != (root / "CITATION.cff").read_bytes():
        raise ValueError("Payload citation differs")
    if "[x]" in (root / "PUBLIC_RELEASE_OWNER_REVIEW.md").read_text().lower():
        raise ValueError("Owner approval must not be synthesized during staging")
    result = {
        "status": "PASS",
        "files": len(names),
        "imported_files": len(entries),
        "frozen_files": len(frozen),
        "licence_classified_files": len(ledger),
        "privacy": privacy(root, names, structured=development),
        "doi_state": meta["doi_state"],
    }
    if development:
        import jsonschema
        import yaml

        schema = read(root, "schemas/release-metadata.schema.json")
        jsonschema.Draft202012Validator.check_schema(schema)
        jsonschema.validate(meta, schema)
        parsed = yaml.safe_load(cff)
        jsonschema.validate(parsed, read(root, "schemas/citation-file-format-1.2.0.json"))
        if parsed["authors"][0]["orcid"] != "https://orcid.org/0009-0007-6030-6019":
            raise ValueError("Author identity mismatch")
        validate_links(root, names)
        result["metadata_schema_and_links"] = "PASS"
    return result


def validate_links(root, names):
    # Historical payload contains inspectable source paths, not website navigation.
    for name in names:
        if name.startswith("publication/v0.2/") or not name.endswith(".md"):
            continue
        for target in re.findall(r"\]\(([^)]+)\)", (root / name).read_text(encoding="utf-8")):
            if "://" in target or target.startswith(("#", "mailto:")):
                continue
            if not ((root / name).parent / target.split("#")[0]).exists():
                raise ValueError("Broken local link: " + name + " -> " + target)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--development", action="store_true")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    report = validate(development=args.development)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))
