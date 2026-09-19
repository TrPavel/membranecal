# S1 Text Scientific supplement

This supplement documents the frozen MembraneCal v0.2 study and post-outcome reproduction scope. All effect units are percentage points. It adds reporting and navigation, not an analysis specification. Paths refer to this extracted release; numerical tables are also collected in S1 Dataset. Frozen JSON basenames refer to data/.

## Table A Eligibility and mapping criteria

The following entries are rendered directly from the frozen protocol. Fractions are unit-scale proportions; lengths are residues; distance thresholds are Å. Resolution is at most 3.5 Å for X-ray or EM structures. Noncanonical isoforms are excluded. Source boundary and chronological decisions remain frozen.

| Frozen field | Criterion |
| --- | --- |
| candidate_cap | 800 |
| canonical_length | [100, 1500] |
| construct_canonical_length_ratio_min | 0.5 |
| experimental_observed_fraction_min | 0.7 |
| identity_min | 0.98 |
| independence | ['canonical accession absent from all 174 v0.1 accepted', 'experimental entry absent from all v0.1 accepted', 'exact construct SHA256 absent from all v0.1 accepted'] |
| mapped_observed_fraction_min | 0.95 |
| mapped_residues_min | 100 |
| opm_chain_fraction_min | 0.8 |
| opm_half_thickness_disagreement_max_angstrom | 0.1 |
| opm_rmsd_max_angstrom | 0.5 |
| primary_core_residues_min | 20 |
| representative | one deposited label-asym per entity; one canonical/construct across entries sorted resolution, descending mapped residues, SHA256(case ID) |
| scope | alpha-helical transmembrane OPM, resolution <=3.5 A, X-ray or EM |
| mapping.experimental_model | first deposited model; selected before outcome access |
| mapping.missing_coordinates | no imputation; endpoint uses paired observed CA only; zero neighbours undefined |
| mapping.missing_frozen_matched_endpoint | abort primary computation rather than silently rematch or drop outcomes |
| mapping.mutations | <=2% retained and separately sensitivity-excluded |
| mapping.noncanonical_isoforms | excluded |
| mapping.source | SIFTS XML author-chain/residue to UniProt canonical; positive unique positions; reference highest-occupancy CA with altloc lexical tie break; no prediction alignment fitted |

## Cohort inventory and exclusions

`data/proteins.json` is the full 554-protein inventory, including accessions, constructs, entry identifiers, source versions and QC metadata. `discovery.json` and `exclusions.json` document index selection, bounded candidate sampling, preparation failures and duplicate canonical selection. Figure 1 and `tables/cohort_flow.csv` distinguish 8,915 index records, 800 candidates, 509 prepared entries, 622 eligible cases and 554 retained proteins from 473 entries. `development_sequences.json` supplies all 174 accepted development proteins for canonical/alias, entry and construct exclusions. No overlap was found; this does not exclude training exposure.

## Grouping and dependence

`groups.json` contains the frozen operational groups, sequence/rescue edges, curated homology links and dependence blocks. Primary sequence links use local BLOSUM62 alignment with identity at least 0.30 and paired full-sequence coverage at least 0.80; rescue still requires identity at least 0.30 and paired coverage at least 0.50, the same explicit whole-family annotation and identical nonempty Pfam complement. Shared domains alone are insufficient. Domain order and repeat counts are not resolved by the complement.

The whole cohort has 335 resolved operational groups and 53 unresolved proteins pooled as one group. Equal proteins within equal operational groups determine the point estimate. Resampling whole dependence blocks preserves linked group means and recomputes the sampled group-sum/group-count ratio. Blocks are not equally weighted in the point estimate. Curated links and shared experimental entries define recorded dependence; annotation-or-edge resolution is not proof of evolutionary family membership. Pooling unresolved proteins changes the target and is not guaranteed conservative.

## Support and balance

`matched_sets.json`, `common_support.json`, `feasibility.json` and `tables/support_balance.csv` document exact HELIX/SHEET matching and two-point confidence bins, three residues per region/cell and at least 20 common-mass units per protein. Interface/core retains 437 proteins, 260 resolved plus one pooled group and 186 blocks; structured extra/core retains 229 proteins, 140 resolved plus one pooled group and 95 blocks. Eligible-residue denominators are 107,096 and 120,071, already restricted by structure and region, not all 223,714 scored residues. Raw retention is 79.31%/38.46%, effective mass retention 45.27%/23.81%. Maximum absolute residual protein confidence differences are 0.634/0.595 pLDDT points. Figure 3 displays both mass loss and residual balance.

## Table B Prespecified sensitivities and secondary contrast

`tables/all_scenario_effects.csv` supplies all 36 effects: primary plus eleven frozen sensitivity scenarios, each with three contrasts. `tables/sensitivity_effects.csv` and Figure 4 display the two focal contrasts. The secondary structured-extra/interface contrast has 95% intervals; focal contrasts have 97.5% marginal intervals. No simultaneous across-scenario guarantee is claimed. Only seven prespecified scenarios feed the frozen threshold-and-sign sensitivity flag; no flag is not evidence of universal method invariance.

Exact full-sequence equality applies to 198 whole-cohort proteins and 152/75 matched proteins; it does not establish equal biological state. Family novelty applies to 366 whole-cohort proteins and 284/138 matched proteins under the recorded curated component definition, excluding unresolved proteins and components touching development. The post-2018-release sensitivity is a date proxy, not measured model-training membership. Original `results.json` supplies scenario-specific eligibility and weights; the supplement does not substitute a common population across contrasts.

## Table C Influence diagnostics

`tables/influence_diagnostics.csv` renders every frozen primary leave-one-operational-group-out estimate and its signed change from the full estimate, including the secondary contrast. Largest absolute changes for the two primary contrasts are 0.075 and 0.127 points. These omit operational groups, not whole curated components or bootstrap blocks. No new influence-based exclusions were applied.

## Metric verification and descriptive cases

The endpoint uses nonself observed mapped Cα neighbors within the selected chain and inclusive 15 Å, with inclusive 0.5, 1, 2 and 4 Å distance-error thresholds. Interchain contacts are not evaluated. Undefined scores would not be silently imputed; none occurred. Positive discrepancy means local agreement exceeds pLDDT/100, negative means confidence exceeds agreement under this endpoint.

Three verification scopes remain distinct: original vectorized compact-coordinate replay on all 554 proteins; an independent scalar comparison on five proteins/2,341 residues; and retained independent raw mmCIF-column coordinate checks on 13 proteins/4,716 residues. The offline entry command repeats compact scoring, independently reproduces primary statistics and hash-verifies retained raw-check evidence. It does not redownload raw structures or independently recompute all sequence alignments and alternative grouping maps.

Figure 5 and `tables/descriptive_cases.csv` give the eight largest absolute all-observed-residue protein-average discrepancies selected after outcomes, with source titles, DOIs, raw hashes, construct lengths and mapped identity. These are unmatched descriptive means without inferential intervals. Context such as pore, prepore or inhibited carrier does not establish the AFDB predicted state or prove a causal mechanism. Full investigation: `provenance/source_investigation.json`.

## Reproduction and licences

Run `python -m membranecal_release.verify` from the extracted release after installing environment/requirements.lock. The content manifest binds every included file; the original scientific SHA identifiers are recorded separately. All 36 point estimates and two primary bootstrap intervals are independently checked within 1e-12; the original scientific functions replay all reported sensitivity intervals. Presentation is post-outcome.

DATA_LICENSE.md and provenance/opm_review.json describe attribution, transformations and remaining source-scope qualifications. No authored licence has been selected. manifest.json binds this snapshot; provenance/authority.json records the original scientific identifiers. Numerical reproduction uses the established tolerance, while figure rendering may differ by platform.

## Release navigation

See README.md for the single reproduction command, DATA_LICENSE.md for source-specific notices and provenance/authority.json for external scientific identifiers. S1 Dataset contains the numerical tables.
