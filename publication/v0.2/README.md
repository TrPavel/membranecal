# MembraneCal v0.2 reproducibility package

MembraneCal v0.2 is a frozen scientific study distributed as repository/package release v0.2.0. This package contains the preserved manuscript, five frozen publication figures, S1/S2 supporting information, compact frozen inputs, numerical result tables, standalone numerical replay and post-outcome validation code. Reproducibility package DOI: 10.5281/zenodo.22843162. No journal publication or peer review is claimed.

The DOCX is the corrected PLOS submission version: funding, competing interests and CRediT are supplied through submission metadata; references use verified Vancouver formatting. The unchanged reading PDF/Markdown preserve the complete owner reading version and declarations. Both have the same scientific content and citation mapping. See provenance/compliance_changes.json and provenance/reference_style_audit.json.

## Offline replay

Use Python 3.12 in a fresh virtual environment. Install `python -m pip install --require-hashes -r code/numerical/environment/requirements.lock`. From this extracted top-level directory, unset PYTHONPATH, then run `python replay.py --output ../replayed`. Verify only with `python replay.py --hashes-only --output ../unused`. Do not use optimized Python. Inputs need no raw upstream cache or private checkout. The code regenerates compact scores, all statistical outputs, nine original tables and five publication figures; primary bootstrap checks use 20,000 draws each. Figure pixels depend on the documented Matplotlib/Pillow versions and installed Arial fonts; numerical tables use the original tolerance. The minimum raw-resource acquisition/mapping/grouping source and bounded post-outcome diagnostics are under code/revision. See its upstream/README.md for historical-source reconstruction limits. Full current-source acquisition is not a substitute for this historical freeze.

The original immutable RC1 hash is recorded for provenance, but old manuscripts and owner-planning forms are not the current publication materials. This package uses unchanged scientific modules and a new packaging/replay wrapper. Authored code is Apache-2.0, authored text/figures CC BY 4.0; third-party data remain source-specific. See licence_ledger.json. S2 is post-outcome validation, not a new primary study.

MANIFEST.json binds every payload file except itself and CHECKSUMS.sha256; CHECKSUMS.sha256 additionally binds MANIFEST.json. The external archive identity binds the entire container.

Current administrative availability wording is in ../DATA_AVAILABILITY.md. Preserved manuscript and historical receipts predate DOI injection; their scientific bytes are unchanged.
