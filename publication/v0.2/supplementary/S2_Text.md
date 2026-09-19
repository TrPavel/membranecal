# S2 Text. Post-outcome validation and revision evidence

Pavel Trofimchik. MembraneCal v0.2 manuscript revision, 15 September 2026.

## Scope and unchanged analysis

This addendum responds to the author-requested internal adversarial audit dated 15 September 2026. It is an AI-assisted project check, not external human peer review. No primary protein, residue, region boundary, matching stratum, group, block, exclusion, score, confidence interval or decision margin has been replaced. The diagnostics below were added after outcomes. They clarify the original analysis rather than selecting a new primary specification.

The immutable original archive is membranecal-v0.2-publication-rc1.zip, SHA256 5682be4d6b88245de84e13d5e64b70819ec3309f7bd9c24fbc071a5f51a7baf4. Its retained 554 proteins, 473 experimental entries and 223,714 paired-observed residues remain the analysis inputs. Original primary interface-core discrepancy is +0.5827439587742238 percentage points with 97.5% interval [+0.2342525902034268, +0.9693052051067319]; structured extra-core is +0.0648016673628679 [-0.8088097517660738, +0.9397327572570777]. Original S1 Text and S1 Dataset remain available in RC1; their historical owner/licence status is superseded only for the authored revision by the owner decisions below. Retain their source-specific qualifications.

## Focused comparison with prior work

Table A. The distinct questions addressed by selected primary studies

| Study | Population and endpoint | Analytical focus and relationship to this study |
| --- | --- | --- |
| Hegedűs et al., 2022 | Transmembrane predictions, human membrane-protein confidence and ABC transporter structures | Regional pLDDT, topology and fold plausibility motivate asking whether confidence tracks structural agreement. The present estimand additionally matches confidence/secondary structure within proteins and weights recorded groups. |
| Dobson et al., TmAlphaFold, 2023 | Predicted alpha-helical membrane proteins; computed membrane placement and structural conflicts | Membrane localization and quality annotation provide context for spatial labels. These annotations do not themselves estimate observed-reference agreement minus confidence in paired common support. |
| Jambrich et al., 2023 | Human transmembrane proteome; experimental/template and predicted structural coverage | Coverage, confidence and sequence relationships motivate explicit selection and exposure limits. This study instead estimates two conditional regional differences on available paired coordinates. |
| Edmunds et al., 2024 | AlphaFold models; empirical quality scores including C-alpha lDDT and confidence-based ranking | Tests how accuracy self-estimates relate to empirical model quality and ranking. Our question uses region-specific matching and an explicit group-weighted discrepancy contrast rather than general model ranking. |
| Xie and Huang, 2024 | 16 membrane transporters with 32 inward/outward-facing experimental structures | Tests capture of alternative conformations. The present discrepancy endpoint cannot diagnose an alternative state or establish a causal state mechanism. |
| Present study | 554 development-disjoint proteins; paired-observed intrachain C-alpha agreement minus normalized confidence | Within-protein common support, then equal protein and operational-group weighting; recorded blocks for uncertainty. Primary support is almost entirely high-confidence HELIX. No broad priority claim or universal calibration claim follows. |

Primary sources (accessed 15 September 2026): https://pmc.ncbi.nlm.nih.gov/articles/PMC8761152/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC9825488/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC10662385/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC11322044/ ; https://pubmed.ncbi.nlm.nih.gov/38564295/ . This is a focused rationale, not a systematic literature review.

## What the fixed weights represent

Within protein p and matching cell c, common mass is min(n_left, n_right). Each cell receives this mass divided by the total retained common mass of that protein. Protein contributions are divided by the number of included proteins in its operational group, and group contributions by the number of included groups. Summing these weights over HELIX/SHEET and confidence bins gives actual estimand composition, not raw residue composition.

Interface-core uses 2,074 HELIX cells and one SHEET cell (P39301): HELIX weight 99.982895457%, SHEET 0.017104543%. Structured extra-core uses 1,116 HELIX cells and no SHEET cells: HELIX weight 100%. Confidence bins at least 80 carry 98.650845% and 99.064756%; bins at least 90 carry 83.807050% and 80.553830%. Effective mass retention is 45.27% and 23.81% of eligible regional residues. These restrictions define a narrow high-confidence helical conditional target.

The unchanged weights give Δ(agreement-confidence) = Δagreement - Δconfidence. In percentage points the interface components are +0.484268472238 and -0.098475486537, yielding +0.582743958774; extra-core components are +0.001758500316 and -0.063043167047, yielding +0.064801667363 (rounding at the final digit). Thus the confidence component contributes about 16.9% of the interface point-estimate magnitude algebraically. This is neither an estimate of bias nor a causal decomposition. The weighted observed-neighbor contrasts are approximately -6.2 and -3.2 neighbors: local reference geometry remains part of the endpoint context. Full-precision calculations are in numerical_audit.json and population_diagnostics.json.

## Predictor metadata and chronology

Table B. Recorded AFDB product model-date metadata

| Recorded value | Protein count | Interpretation |
| --- | --- | --- |
| 2022-06-01 | 135 | Product metadata; not a verified inference date |
| 2025-08-01 | 417 | Product metadata; not a uniform predictor configuration |
| 2025-03-31 | 2 | Product metadata; not proof of training exclusion |

All retrieved products carry database version 6. This does not identify an AlphaFold 6 predictor. No homogeneous inference campaign is claimed. Development-disjointness, exact full sequence and temporal subsets address different recorded proxies and do not eliminate training/template exposure.

The analytical freeze is commit 5deb8f9abc345d811e0e98680b6122b4046a0d0d, recorded 9 September 2026 12:26:49 UTC. Computation receipt starts at 12:27:12 UTC (seconds shown here; original receipt retains its native precision). The original scientific source is d9641f28fdf50d42057223e3787f999e1a127198. Protocol SHA256 is 740b223fbf6bcb43c42c126208510a835d04dade7f268d0999b684a7738c6671; scores SHA256 is 4f293fda6d7c7d0a4c8578071f396e7d573a54ed3c472795c67a5a5b9baca4b8. These identify retained historical artifacts, not external certification.

Before the analytical freeze, AFDB CIFs were already downloaded and parsed for sequence/confidence; coordinate bytes were present. The acquisition path did not extract prediction coordinates for endpoint calculation. No physical blinding, certified absence of inspection, or external preregistration is claimed. The new diagnostics and textual edits are explicitly post-outcome.

## Direct published-function endpoint comparison

The unchanged file alphafold/model/lddt.py was retrieved from Google DeepMind commit dbaafbcdea0cf39fabe502928cab55754c1d5dc7. File SHA256: 2a3bb36fed58cff0beb92da0af40ca74c5b9892947f01ab825e7547691d99d96. Source: https://github.com/google-deepmind/alphafold/blob/dbaafbcdea0cf39fabe502928cab55754c1d5dc7/alphafold/model/lddt.py . The preserved file includes the upstream Apache-2.0 notice.

All 554 proteins and 223,714 residues were evaluated on CPU with Python 3.12, JAX/jaxlib 0.6.2, NumPy 2.3.5, scipy 1.16.1, ml_dtypes 0.5.3 and opt_einsum 3.4.0. Float64 and float32 were run separately. The official function uses strict distance and threshold comparisons plus 1e-10 stabilizers. The frozen score uses inclusive comparisons. Identical paired-observed intrachain coordinates and correspondence were retained; no full-canonical or all-atom target was substituted. Zero-masked padding to multiples of 128 reduced JIT compilation shapes; padding outputs were discarded. This validation uses the pinned public implementation, not a reconstruction of the exact historical training environment. Semantic identity with the precise 2021 training implementation/environment has not been established by this check. No prediction pipeline was executed.

Table C. Direct unchanged JAX function versus frozen endpoint

| Diagnostic | float64 | float32 |
| --- | --- | --- |
| Maximum residue difference, percentage points | 1.40624734e-9 | 1.13636228 |
| Residues with difference greater than 1e-4 on unit scale | 0 | 30 |
| Interface estimate, percentage points | 0.582743958774975 | 0.582738882167135 |
| Interface change, percentage points | 7.51135e-13 | -5.07661e-6 |
| Extra-core estimate, percentage points | 0.064801667364171 | 0.064806934373719 |
| Extra-core change, percentage points | 1.30305e-12 | 5.26701e-6 |

The supplied NumPy port and the direct float32 JAX function differ slightly in primary rounding; the manuscript reports the actual direct JAX values. The direct float64 and float32 runs are retained separately in evidence/direct_endpoint_float64.json and direct_endpoint_float32.json. All original inclusive scores and neighbor counts matched exactly in the supplied audit rerun. This confirms implementation behavior on these retained masks; it does not independently validate every raw mapping or biological state. No alternative endpoint replaced the primary scores.

## Bounded upstream checks

The added raw-case panel was selected using frozen covariates and mapping features without consulting scores. Cases cover the minimum observed fraction (Q6J8X0), minimum mapping coverage (O66905), maximum OPM transform RMSD (Q6J8X0), first accession with alternate location (A0A0H2V8D3), insertion code (A0R612), author/label chain mismatch (A0A0H3CQA2), and mapped identity mismatch (A0A0H2V8D3). Deduplication yields five proteins and 2,725 residues. Selection is recorded before checks execute; it is post-outcome validation, not a newly prespecified primary sample.

The checks parse raw mmCIF atom columns, select chain/entity/model and occupied C-alpha conformers, compare residue identifiers and coordinates, independently parse SIFTS XML correspondence, and recompute the OPM rigid-frame transfer with an independent Kabsch implementation. All five cases passed; maximum depth discrepancy was zero at retained precision. All 15 raw source hashes were checked against the immutable RC1 source records as well as cache metadata. These hashes and per-case counts appear in bounded_upstream.json and raw_source_bindings.json. Raw source objects are not redistributed.

Fifteen sequence pairs were selected from recorded full-length/rescue edges nearest the identity or coverage thresholds and from length-nearest candidates for the first three unresolved accessions. The original exposed alignment and edge functions were rerun and compared with recorded evidence. All selected checks passed. This validates selected transformation replay, not independent alignment methodology, all-pairs graph completeness, domain-order handling or remote evolutionary homology. Complete sequence-edge discovery was not rerun. Earlier checks on 13 proteins/4,716 residues comprised eight outcome-selected discrepancy cases and five accession-sorted cases; they are not a random sample. The new panel broadens inspected features without converting targeted checks into full-source validation.

Largest recorded primary blocks contain 11 of 261 and eight of 141 groups. Maximum leave-one-block changes are approximately 0.0748 and 0.1274 percentage points. These diagnostics concern recorded dependence, not unrecorded shared ancestry. Pooling unresolved proteins is not described as guaranteeing conservative coverage.

## Reproduction and source reconstruction

The original offline verifier reproduces compact-coordinate scores, original scenario results and two independently coded primary bootstraps. It does not rebuild the historical internet inputs. The minimal original acquisition, freeze/mapping and group-generation sources are included byte-for-byte under upstream/scripts with provenance hashes. upstream/scaffold.py creates a separate source workspace from retained metadata after validating RC1 integrity. The source README explains current-source acquisition and the historical index/cache requirements without assuming private Git history.

Source reconstruction requires the identified historical objects or source refresh, appropriate source permissions, and explicit comparison of hashes. The original full candidate annotations are not shipped. The historical freeze CLI expects additional original design/panel artifacts and is not advertised as a one-command independent reconstruction. Selected raw transformations and alignment edges were executed; a fresh complete upstream cohort rebuild was not. These are explicit limits, not silently closed by successful numerical replay.

## Audit-response ledger

| Item | Disposition and evidence | Remaining scope limitation |
| --- | --- | --- |
| R01 | Corrected title, abstract, Methods, Results Table 2 and Discussion; verified actual estimand weights | Low-confidence, SHEET and unmatched populations are not covered by the primary inference |
| R02 | Direct unchanged upstream function run on all 554 proteins in both dtypes; Table C and executable reports | Same observed mask only; not all-atom/full-assembly validation |
| R03 | Main Table 3 reports agreement/confidence decomposition; composite units corrected | Residual matching imbalance is not a causal bias estimate |
| R04 | Continuous estimates lead; five-point rule retained as investigator convention | No validated biological/utility margin or formal equivalence design |
| R05 | Group assumptions clarified, complete generation source exposed, 15 selected pairs replayed | No claim of exhaustive homology discovery or conservative pooling |
| R06 | Five covariate-selected raw cases/2,725 residues checked with source hashes | No claim that all 554 raw mappings were independently checked |
| R07 | Minimal upstream source, scaffold and bounded reconstruction instructions added | Historical source availability and permissions remain external requirements |
| R08 | Analytical barrier and preexisting coordinate bytes stated; hashes moved here | Local Git/receipt chronology is not external preregistration or physical blinding |
| R09 | Mixed product dates tabulated; v6 terminology and exposure limits clarified | Product fields do not establish inference dates or training independence |
| R10 | Same-mask endpoint and neighbor-context contrasts disclosed | Missing residues, interchain context and biological state remain limitations |
| R11 | All 11 sensitivities retained; exact-sequence interval and block influence described | Flag stability is narrower than universal method robustness; no subgroup-difference inference |
| R12 | Unmatched cases labeled descriptive; source titles separated from causal hypotheses | State/assembly mechanisms not demonstrated |
| R13 | Positive rationale and focused prior-study comparison added | No systematic priority or first-calibration claim |
| R14 | Owner identity, funding, CRediT, AI scope and selected licences added | Affiliation/email confirmed by PLOS; scoped numerical OPM permission received; Zenodo URL/DOI and final owner approvals remain outstanding |
| R15 | Prefilter order corrected in text and Figure 1; hashes/internal identifiers moved to SI | Original RC1 historical files intentionally retain their original wording |

Every scientific comment is corrected in text, checked diagnostically, or retained as an explicit scope limit. External owner/source/deposit requirements are not marked complete by this ledger.

## Owner decisions and disclosure

Pavel Trofimchik is the sole author and corresponding author, trofimchikpavel@icloud.com; ORCID https://orcid.org/0009-0007-6030-6019 . Current affiliation is Independent Researcher, Saint Petersburg, Russian Federation. PLOS ONE confirmed this affiliation and the personal corresponding email as acceptable (case 09849353, owner-supplied correspondence).

The author received no specific funding for this work. The author develops MembraneState Bench, the independent research and benchmarking programme within which this study was conducted, and is exploring potential future commercial applications of related benchmarking infrastructure. No commercial entity funded this study. All supplied CRediT roles are included in the manuscript. This disclosure does not invent a company, income or additional conflicts.

AI assistance was substantial in methodological discussion, software development, code review, literature-assisted drafting, manuscript preparation and reproducibility packaging. OpenAI ChatGPT substantially assisted methodological critique, literature-assisted reasoning and manuscript development; OpenAI Codex was used for software development, revision and verification. AI systems are not authors; the human author retains final scientific decisions, verification, interpretation and responsibility. Unprovided historical model identifiers and session histories are not inferred. The author-requested audit and independently coded checks are internal activities, not independent external human peer review.

The author selected Apache-2.0 for authored code and CC BY 4.0 for authored text and figures. The addendum licence notice scopes this choice to authored material; it does not relicense third-party data. A Zenodo deposit is planned, with no final URL or DOI yet. Explicit permission from Andrei Lomize covers numerical OPM values described for this reproducibility package, including required orientation, depth and thickness information, with attribution and transformation provenance retained. Unrelated raw assets and arbitrary database redistribution are outside this permission. The package is a local review revision, not an assertion of completed public access or journal submission readiness. No submission, merge, publication, DOI creation or contact with third parties was performed as part of this revision.
