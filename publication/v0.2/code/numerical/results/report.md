# MembraneCal independent validation v0.2

MEMBRANECAL_INDEPENDENT_VALIDATION_V0_2_COMPLETE

MEMBRANECAL_V0_2_NO_MATERIAL_REGIONAL_EFFECT_ESTABLISHED

PRE_OUTCOME_FREEZE: `5deb8f9abc345d811e0e98680b6122b4046a0d0d`; protocol SHA256 `740b223fbf6bcb43c42c126208510a835d04dade7f268d0999b684a7738c6671`.

## Cohort

800 candidate entries; 554 accepted canonical proteins in 473 experimental entries; 223,714 mapped residues. 335 resolved operational groups; 53 unresolved proteins. Zero development canonical/alias, entry or construct overlaps. Family-novel subset: 366 proteins; full-sequence exact-construct subset: 198 proteins.

## Matched primary results

Effects and intervals are in percentage points. Equal proteins within operational sequence/rescue groups and equal group means; blocks additionally link explicit curated whole-family annotations and shared entries, 20,000 bootstrap draws, Bonferroni 97.5% intervals. Equal curated-homology weighting is a separate sensitivity.

| Contrast | Effect | Interval | Proteins | Resolved groups | Blocks | Useful |
|---|---:|---|---:|---:|---:|---|
| MEMBRANE_INTERFACE-TM_CORE | 0.583 | [0.234, 0.969] | 437 | 260 | 186 | True |
| STRUCTURED_EXTRAMEMBRANE-TM_CORE | 0.065 | [-0.809, 0.940] | 229 | 140 | 95 | True |

## Sensitivity estimates

Intervals below remain descriptive sensitivity intervals. N/F/B are paired proteins, resolved groups and resolved dependence blocks.

| Scenario | Contrast | Effect, points | Interval, points | N/F/B | Useful support |
|---|---|---:|---|---|---|
| bins_1 | MEMBRANE_INTERFACE-TM_CORE | 0.572 | [0.237, 0.955] | 398/240/173 | True |
| bins_1 | STRUCTURED_EXTRAMEMBRANE-TM_CORE | -0.075 | [-0.945, 0.836] | 187/111/76 | True |
| bins_5 | MEMBRANE_INTERFACE-TM_CORE | 0.364 | [0.049, 0.702] | 462/275/195 | True |
| bins_5 | STRUCTURED_EXTRAMEMBRANE-TM_CORE | -0.051 | [-0.898, 0.793] | 263/157/105 | True |
| v01_compatible | MEMBRANE_INTERFACE-TM_CORE | 0.344 | [0.012, 0.701] | 459/276/196 | True |
| v01_compatible | STRUCTURED_EXTRAMEMBRANE-TM_CORE | 0.098 | [-0.765, 0.927] | 268/158/105 | True |
| family_novel | MEMBRANE_INTERFACE-TM_CORE | 0.807 | [0.412, 1.246] | 284/206/165 | True |
| family_novel | STRUCTURED_EXTRAMEMBRANE-TM_CORE | 0.044 | [-1.106, 1.181] | 138/101/79 | True |
| exact_construct | MEMBRANE_INTERFACE-TM_CORE | 0.260 | [-0.350, 0.965] | 152/103/78 | True |
| exact_construct | STRUCTURED_EXTRAMEMBRANE-TM_CORE | 0.695 | [-0.680, 2.177] | 75/54/43 | True |
| unresolved_excluded | MEMBRANE_INTERFACE-TM_CORE | 0.584 | [0.234, 0.971] | 393/260/186 | True |
| unresolved_excluded | STRUCTURED_EXTRAMEMBRANE-TM_CORE | 0.038 | [-0.827, 0.919] | 212/140/95 | True |
| broad_family | MEMBRANE_INTERFACE-TM_CORE | 0.688 | [0.264, 1.174] | 437/170/168 | True |
| broad_family | STRUCTURED_EXTRAMEMBRANE-TM_CORE | 0.257 | [-0.890, 1.356] | 229/90/89 | True |
| strict_sequence_groups | MEMBRANE_INTERFACE-TM_CORE | 0.562 | [0.223, 0.941] | 437/270/186 | True |
| strict_sequence_groups | STRUCTURED_EXTRAMEMBRANE-TM_CORE | 0.083 | [-0.764, 0.955] | 229/145/95 | True |
| identity_100 | MEMBRANE_INTERFACE-TM_CORE | 0.600 | [0.189, 1.070] | 338/212/165 | True |
| identity_100 | STRUCTURED_EXTRAMEMBRANE-TM_CORE | 0.451 | [-0.501, 1.429] | 173/109/81 | True |
| curated_homology | MEMBRANE_INTERFACE-TM_CORE | 0.615 | [0.215, 1.060] | 437/201/186 | True |
| curated_homology | STRUCTURED_EXTRAMEMBRANE-TM_CORE | 0.158 | [-0.869, 1.158] | 229/106/95 | True |
| post_2018_release | MEMBRANE_INTERFACE-TM_CORE | 0.656 | [0.232, 1.116] | 276/178/131 | True |
| post_2018_release | STRUCTURED_EXTRAMEMBRANE-TM_CORE | 0.160 | [-0.830, 1.168] | 178/117/80 | True |

## Common support

| Contrast | Raw retained-cell fraction | Effective mass fraction | Max absolute residual pLDDT difference, points |
|---|---:|---:|---:|
| MEMBRANE_INTERFACE-TM_CORE | 0.793 | 0.453 | 0.634 |
| STRUCTURED_EXTRAMEMBRANE-MEMBRANE_INTERFACE | 0.366 | 0.210 | 0.475 |
| STRUCTURED_EXTRAMEMBRANE-TM_CORE | 0.385 | 0.238 | 0.595 |

## Descriptive calibration

Pooled confidence–accuracy Pearson correlation: 0.6454. This residue-level description is not an independent-replicate inferential test.

| Region | Residues | Pooled error, points | Protein-balanced error | Family-balanced error |
|---|---:|---:|---:|---:|
| ALL | 223714 | -0.093 | -0.094 | -0.402 |
| TM_CORE | 74017 | -0.376 | -0.374 | -0.452 |
| MEMBRANE_INTERFACE | 33079 | 0.227 | 0.175 | -0.019 |
| STRUCTURED_EXTRAMEMBRANE | 46054 | 0.086 | 0.921 | 0.555 |

## Family influence

MEMBRANE_INTERFACE-TM_CORE: omit family:O31539: 0.508 points (-0.075 shift); omit family:G2Q5N0: 0.536 points (-0.046 shift); omit family:Q03S56: 0.538 points (-0.045 shift).
STRUCTURED_EXTRAMEMBRANE-TM_CORE: omit family:G2Q5N0: -0.063 points (-0.127 shift); omit family:P10903: 0.175 points (+0.110 shift); omit family:F1NCD6: 0.170 points (+0.105 shift).

## Largest protein-average discrepancies

These are descriptive case investigations; none was excluded after outcomes were seen.

| Protein | PDB | Mean error, points | Construct | Missing experimental fraction |
|---|---|---:|---|---:|
| P77335 | 6mrt | -38.677 | sequence differs from canonical | 0.120 |
| Q2FZP8 | 6s7v | -22.118 | sequence differs from canonical | 0.000 |
| G2QNH0 | 6gci | -21.119 | sequence differs from canonical | 0.073 |
| Q6LPW0 | 7qha | 15.309 | full sequence identical | 0.100 |
| G3IEF0 | 8djm | -14.277 | sequence differs from canonical | 0.066 |
| P30878 | 7l17 | 13.763 | sequence differs from canonical | 0.066 |
| A5F5Y6 | 8evu | 13.604 | full sequence identical | 0.038 |
| A0KLE1 | 6h2f | -13.366 | sequence differs from canonical | 0.098 |

## Matching and grouping

Frozen 2-point pLDDT bins × HELIX/SHEET; >=3 residues per side/cell and >=20 common mass. Raw retained-cell fraction differs from effective mass retention: see frozen feasibility.json/common_support.json. Continuous confidence equality is not claimed.

Combined development/validation sequence graph, full-length >=30% identity and >=80% paired coverage; curated-family plus complete Pfam-complement plus >=50% coverage rescue. Shared Pfam alone never creates an edge. Primary bootstrap blocks additionally join exact normalized whole-family annotations regardless of sequence threshold. Family novelty uses those homology links, not shared-entry links. Bounded known-family validation passes; ordered architecture/repeat uncertainty remains. Unresolved proteins are conservatively pooled. These are operational groups, not proven evolutionary-independent units.

## Sensitivities and leakage

results.json includes all frozen 1/5-bin, v0.1-compatible, family-novel, exact-construct, unresolved-excluded, broad/strict grouping, 100%-identity and temporal sensitivities, plus family means and leave-one-family-out estimates.

2018-04-30 is the published original AlphaFold training structure cutoff used as a date-based diagnostic. AFDB model version/date is recorded individually. Post-cutoff does not establish absence of later model/template/sequence exposure; no individual training membership is claimed.

## Interpretation

Material support requires an adjusted interval entirely beyond ±5 points. No material effect established is not equivalence or proof of a zero effect. Scope is protein-independent alpha-helical transmembrane canonical-AFDB benchmarking with observed structural support. Experimental state/construct and computational orientation limitations remain.

## Failure modes

Twenty largest absolute protein-average errors are recorded with construct, mutation, assembly, missingness and membrane-frame evidence. State mismatches and mechanistic explanations remain unresolved without direct evidence; no post-outcome exclusion or retuning.

## Reproduction

Run `python -m membranecal_release.verify` using environment/requirements.lock. This replays compact-coordinate scores and all statistical outputs offline. Raw source archives are omitted; source terms and attribution are documented in DATA_LICENSE.md.
