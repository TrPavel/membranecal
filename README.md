# MembraneCal

[![CI](https://github.com/TrPavel/membranecal/actions/workflows/ci.yml/badge.svg)](https://github.com/TrPavel/membranecal/actions/workflows/ci.yml)
[![DOI reserved](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.22843162%20reserved-lightgrey)](https://doi.org/10.5281/zenodo.22843162)

**Do matched membrane regions differ in local structural agreement minus model confidence?**
MembraneCal study v0.2 evaluates this question for high-confidence helical residues in
membrane proteins, using existing AlphaFold Database predictions and experimentally
observed structures. It is a frozen scientific study with an offline reproducibility package.

**Status:** private release-candidate preparation for repository version **0.2.0**.
The Zenodo DOI is reserved; its record is a draft, not a published dataset.
The manuscript is prepared for PLOS ONE; no journal publication is claimed.

## What v0.2 evaluates

The capped cohort contains **554 proteins and 223,714 paired-observed residues**.
The primary matched support is high-confidence **HELIX** support, not all residues
in that cohort. Within-protein matching compares membrane interface and structured
extramembrane regions with the transmembrane core. The endpoint is a
superposition-free, intrachain C-alpha local-distance agreement minus pLDDT/100.
See [methods](docs/methods.md) for weighting, uncertainty and the exact frozen protocol.

| Primary contrast | Estimate, percentage points | 97.5% interval |
| --- | ---: | ---: |
| Interface minus core | +0.583 | [+0.234, +0.969] |
| Structured extramembrane minus core | +0.065 | [-0.809, +0.940] |

The interface contrast is small and positive in this matched support. The
extramembrane interval includes zero. Neither contrast establishes an effect outside
the investigator-chosen ±5 percentage-point materiality margin. This is not an
equivalence test or evidence of universal calibration.

## Scope and limitations

This study does not provide a universal protein correctness score, biological
validation, stability or function prediction, or a test of alternative-state recovery.
Construct correspondence, observed residues, exposure controls, membrane-frame
provenance and the C-alpha/intrachain endpoint limit interpretation.
Read the [limitations](docs/limitations.md) before using the numbers.

## Reproduce

Use a fresh **Python 3.12** environment. Dependency installation requires network
access; the numerical replay itself is offline and needs no GPU or account.

```bash
git clone https://github.com/TrPavel/membranecal.git
cd membranecal
python -m venv .venv
# Activate .venv using your shell's activation command.
python -m pip install --require-hashes -r publication/v0.2/code/numerical/environment/requirements.lock
python tools/verify.py
python tools/replay.py --mode compact --output ../membranecal-compact
```

While the repository is private, cloning requires authorized GitHub access.
Unset `PYTHONPATH`; do not use `python -O`. Compact replay verifies the frozen
coordinate-derived scores. For all 36 scenario estimates, bootstrap checks, nine
tables and five figures, use the [full reproduction instructions](docs/reproduction.md).
Full pixel-identical figure replay targets Windows with the documented Arial fonts.

## Frozen identity and availability

- Scientific study: **MembraneCal v0.2**; proposed repository/package release: **v0.2.0**.
- [Release metadata](release/metadata.json), [frozen file identities](release/frozen-files.json),
  [complete repository manifest](release/manifest.json).
- [Standalone study payload](publication/v0.2/README.md),
  [result tables](publication/v0.2/data/result_tables),
  [manuscript](publication/v0.2/manuscript/manuscript.md),
  [supporting information](publication/v0.2/supplementary).
- Reserved package DOI: **10.5281/zenodo.22843162**. Availability is currently private
  staging; see the [Data Availability draft](publication/DATA_AVAILABILITY.md).

Compact inputs support offline result reproduction. Historical raw upstream files
are not redistributed in full, and today's upstream downloads are not an exact
reconstruction of the historical freeze. No tag or GitHub Release has been created.

## Cite the exact version

Pavel Trofimchik. *MembraneCal v0.2: reproducibility package for matched regional
contrasts in structural agreement minus confidence among high-confidence helical
residues in membrane proteins*. Package version 0.2.0 (release candidate).
DOI **10.5281/zenodo.22843162**, reserved, not yet published.

Use [CITATION.cff](CITATION.cff) and record the exact release/commit used.
Author: [Pavel Trofimchik](https://orcid.org/0009-0007-6030-6019), Independent Researcher,
Saint Petersburg, Russian Federation. There is no article DOI yet.

## Rights and provenance

Author-created code is [Apache-2.0](LICENSE); author-created text and figures are
CC BY 4.0. Third-party and derived data retain source-specific terms:
[data licences](DATA_LICENSE.md), [notices](THIRD_PARTY_NOTICES.md),
[machine-readable ledger](release/licence-ledger.json).

The repository starts with a clean export history. Its [derivation map](release/derivation-map.json)
links each imported member to the audited source archive and records any metadata
adaptation. Source commit identities are provenance attestations, not a claim that
private Git history or independently certified preregistration is supplied.
See [provenance](docs/provenance.md).

MembraneCal is an independent scientific line within the wider **MembraneState Bench**
programme. This repository contains only MembraneCal; no wider research checkout or
service is needed for replay. [Contributions](CONTRIBUTING.md) and
[release policy](docs/release-policy.md) preserve the frozen study boundary.
