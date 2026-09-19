# Export provenance

The repository is a new Git history containing only a reviewed standalone package
and repository support files. It carries no parent research repository Git objects.
The source archive SHA256 is recorded in `release/source-audit.json`.

The captured build input was commit `c63b3b6266a27f43a56734e51d0198a015813ae7`,
tree `803f18bd9a480482d57e0a4fc7b4d33f208b5ca6`; the audit record was at
`5f143a170fb70b800a40b167c105b20cb2dddff1`. Scientific source commit:
`d9641f28fdf50d42057223e3787f999e1a127198`; pre-outcome freeze:
`5deb8f9abc345d811e0e98680b6122b4046a0d0d`.

Each entry in [derivation-map.json](../release/derivation-map.json) supplies source
commit, relative source path, source SHA256, transformation type, destination and
destination SHA256. A generated archive member is explicitly typed as such: it is
not falsely represented as a tracked file in that commit. The source archive hash
binds its bytes. Retained scientific binding and code-derivation records provide
the earlier source lineage.

Allowed changes are unchanged copy, path/wrapper adaptation, metadata update and
privacy redaction. Scientific inputs, implementations, results and figure files
remain pinned byte-for-byte. New repository documentation/tools have a separate
authored-file licence inventory and enter the clean Git history normally.

Current DOI/status metadata supersedes earlier administrative statements embedded
in the preserved manuscript and historical receipts. Those originals remain
unchanged for traceability; see `publication/DATA_AVAILABILITY.md`. A reserved DOI
is not public availability. The original manuscript PDF/DOCX is not a newly
submitted journal manuscript as a consequence of this export.

Checksums demonstrate content identity. They do not independently certify when
private commits were created. The chronology remains a local Git/receipt record,
not externally registered preregistration. No full private history is needed for
the tested public numerical replay.
