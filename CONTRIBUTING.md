# Contributing

Open a focused issue describing the problem, affected version and a minimal
reproducer. Do not attach credentials, private structures or restricted datasets.
Scientific questions should identify the relevant claim, method and evidence.

Separate documentation, engineering and scientific proposals. Frozen v0.2 files
must remain byte-identical. A changed cohort, endpoint, matching rule, weighting or
uncertainty method requires a new scientific identity and an explicit comparability
statement; it cannot enter as a packaging patch.

Install the frozen dependencies and `requirements-dev.txt`, then run:

```bash
python -m ruff check tools tests
python -m ruff format --check tools tests
python -m pytest -q
python tools/verify.py --development
```

Frozen imported code is excluded from reformatting and lint rewrites. Its integrity
and behaviour are checked through pinned hashes and replay. New repository tools
are linted and tested. Deliberate non-scientific file additions require reviewing
`release/allowlist.json`, the derivation record where applicable, and regenerated
manifests via `python tools/package.py --refresh`.

Contributions are accepted under the licence applicable to the changed material:
Apache-2.0 for authored code and CC BY 4.0 for authored prose/figures. Confirm that
you have the right to contribute; retain all third-party attribution. Maintainers
remain responsible for checking AI-assisted contributions and scientific claims.
There is no automatic authorship entitlement from a software contribution.
