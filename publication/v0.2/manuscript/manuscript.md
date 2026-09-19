# Matched regional contrasts in structural agreement minus confidence among high-confidence helical residues in membrane proteins

Pavel Trofimchik

Independent Researcher, Saint Petersburg, Russian Federation

Correspondence: trofimchikpavel@icloud.com

ORCID: https://orcid.org/0009-0007-6030-6019

Short title: Matched membrane agreement confidence contrasts

Keywords: AlphaFold; membrane proteins; local structural agreement; confidence; matched comparison.

## Abstract

**Background.** Regional confidence distributions can reflect structural composition rather than differences in confidence relative to observed accuracy. We asked whether local structural agreement minus normalized AlphaFold confidence differs between membrane regions within proteins after matching secondary structure and confidence strata.

**Methods.** The validation cohort contained 554 canonical proteins from 473 experimental entries, disjoint from 174 development proteins by accession, experimental entry and construct. An analytical protocol was fixed before endpoint computation. Interface and structured-extramembrane residues were separately matched to core residues within proteins using HELIX/SHEET categories and two-point pLDDT bins. Cell weights reflected common support, followed by equal protein and operational-group weighting. Recorded dependence blocks were resampled 20,000 times to obtain 97.5% marginal intervals for two primary contrasts.

**Results.** Interface-minus-core discrepancy was +0.583 percentage points (97.5% interval +0.234 to +0.969; 437 proteins); structured-extramembrane-minus-core was +0.065 (-0.809 to +0.940; 229 proteins). HELIX residues contributed 99.9829% and 100% of the respective estimand weights; pLDDT strata at least 80 contributed 98.6508% and 99.0648%. Effective common-support mass retained 45.27% and 23.81% of eligible residues. Post-outcome validation against the published AlphaFold score function on the same observed-residue mask changed either primary estimate by less than 0.000006 percentage points.

**Conclusions.** A small positive interface discrepancy contrast is supported in predominantly high-confidence helical common support; the extramembrane interval includes zero. These estimates quantify a narrow conditional comparison, not a rise in interface pLDDT, general membrane-protein calibration or a causal membrane effect. The fixed five-point decision convention is retained without assigning it biological significance or claiming formal equivalence.

## Introduction

AlphaFold confidence scores are useful model annotations, but their interpretation depends on the structural endpoint and evaluation population [1,2,3,4]. The lDDT framework quantifies local structural agreement without requiring global superposition [5]. AlphaFold2 pLDDT itself predicts lDDT-Cα [1]; the original lDDT framework also supports all-atom assessment [5]. Our explicitly specified C-alpha implementation and restriction to observed mapped residues are not claimed to reconstruct every detail of AlphaFold's confidence target or training loss. Subtracting normalized pLDDT yields an operational discrepancy in the evaluated structural context, rather than a complete probabilistic calibration assessment.

Membrane proteins provide a demanding setting because experimental constructs, oligomeric assemblies, ligands and conformational contexts can differ from database models. Computed membrane placement supplies useful spatial annotations [6], but does not directly measure a biological bilayer. Comparisons of core, interface and extramembrane residues can also confound region with secondary structure or confidence. A design that compares residues within proteins and restricts both sides to overlapping covariates offers a narrower empirical question.

Previous studies establish why regional comparison requires more than confidence distributions. Hegedűs and colleagues compared transmembrane confidence and topology and assessed ABC transporter folds [7]. TmAlphaFold evaluates predicted membrane placement and structural conflicts [8], and human-proteome analyses relate structural coverage to confidence and template availability [9]. State-oriented transporter benchmarks ask whether alternative conformations can be captured [10]. General model-quality studies compare pLDDT with observed structural scores and model rankings [11,12]. Our question combines these perspectives: within the same protein, does agreement minus confidence differ by region when secondary structure and confidence strata overlap? Matching, explicit support restrictions and group-weighted uncertainty define a conditional estimand that these broader distribution, placement and ranking comparisons do not directly estimate. Table A in S2 Text compares their endpoints, populations and analytical units; it is a focused comparison, not a systematic priority claim.

An exploratory cohort informed the design. The validation cohort excludes its accepted proteins by canonical and alias accessions, experimental entries and identical constructs. Development-disjointness is distinct from independence from training data, templates or homologous sequences; sequence and date proxies cannot resolve every exposure pathway [13]. We report the continuous matched estimates and uncertainty, followed by the fixed five-point decision rule and all prespecified sensitivity scenarios.

## Methods

### Analytical protocol and provenance

The protocol, selected cohort, covariates, groups, matched sets and analysis code were recorded before local-accuracy endpoint computation. The barrier was analytical: complete AFDB CIF files had already been downloaded and parsed for sequence and confidence, so coordinate bytes were present; the acquisition code did not extract prediction coordinates for the endpoint. The recorded Git freeze time was 9 September 2026 at 12:26:49 UTC and the computation receipt began at 12:27:12 UTC. These internal records support procedural traceability, not physical blinding or independent certification of non-inspection. Only a local analytical freeze is documented; no external preregistration is claimed. Full identifiers, receipt precision and file hashes are supplied in S2 Text and the provenance records.

The primary analysis, cohort, exclusions, matching, grouping, confidence bins, resampling and decision margin are unchanged. The validation diagnostics added during manuscript revision are post-outcome checks of implementation and scope. An independently coded calculation within the same AI-assisted project is not independent external peer review or a new experimental validation cohort.

### Cohort and development exclusions

The OPM index contained 8,915 records. The acquisition code first excluded 65 development-structure overlaps, then 3,540 records failing the alpha-helical-transmembrane/resolution-availability filter, and 1,497 under the resolution/source-boundary filter. Unique-code representatives and per-family rounds selected at most 800 candidate entries, ordered by resolution and a deterministic PDB-ID hash. This capped, nonrandom selection is not the total population of eligible membrane structures. Experimental resolution was required to be at most 3.5 Å. Complete mapping, construct, observability and membrane-transfer thresholds appear in Table A in S1 Text.

Of the 800 entries, 509 produced eligible protein cases, 286 were excluded and five had source/parser/QC failures. The latter included three SIFTS 404 responses and two deposited-resolution QC failures. The eligible entries yielded 622 protein cases; removing 68 duplicate canonical cases left 554 proteins in 473 experimental entries (Figure 1). Counts of entries and protein cases are distinct units. There were 279 X-ray and 275 electron-microscopy protein cases. SIFTS provided structural-to-UniProt correspondence [14,15]. All 174 accepted development proteins were screened using original and current primary accessions, secondary aliases, experimental entry identifiers and construct sequences/hashes. Final overlap counts were zero in all three exclusion categories.

The evaluated products were retrieved AFDB canonical-sequence structures carrying database product version 6 [2,3,4,16]. Their model-date metadata were 2022-06-01 for 135 proteins, 2025-08-01 for 417 and 2025-03-31 for two (Table B in S2 Text). These are recorded product fields, not verified inference dates or evidence of a homogeneous model configuration. Database version 6 is not a predictor named AlphaFold 6. No new prediction campaign was run. Of 554 proteins, 198 had identical full canonical and experimental construct sequences; 366 belonged to the operational family-novel subset. Neither condition establishes absence of training or template exposure.

### Endpoint, observability and membrane regions

For each observed and mapped reference C-alpha, neighbors were all other shared mapped C-alphas within the selected protein chain and an inclusive 15 Å reference distance. The endpoint is intrachain; interchain or assembly contacts are not scored. For each neighbor, reference and predicted distances were compared at inclusive absolute-difference thresholds 0.5, 1, 2 and 4 Å. The residue score was the fraction of these neighbor–threshold tests satisfied. There was no sequence-separation exclusion. Missing prediction correspondence or no neighbors would yield an undefined score; none of the 223,714 evaluated residues was undefined. This is a C-alpha local-distance agreement endpoint related to lDDT [5], not the canonical all-atom implementation. It cannot evaluate unobserved experimental residues or establish biological correctness.

The residue discrepancy was e = C-alpha agreement - pLDDT/100. A positive value means agreement exceeds normalized confidence for this observed-residue endpoint; it does not by itself establish absolute underconfidence or explain a modeling failure. All reported discrepancy effects use percentage points, defined as 100 times a difference on the unit scale. Under the same fixed weights, the regional contrast obeys the identity Δ(agreement - confidence) = Δagreement - Δconfidence. Coarse confidence bins permit a nonzero residual confidence term; the contrast is not an increase in pLDDT.

Experimental coordinates were transferred to the frozen OPM orientation. With `d = |z| / OPM half-thickness`, core comprised `d ≤ 0.8`, interface `0.8 < d ≤ 1.2`, and structured extramembrane `d > 1.6`; the intervening gap was excluded. All primary regions required HELIX or SHEET annotations. These boundaries are computed operational labels, not direct measurements of lipid contacts or bilayer position [6].

### Within-protein matching and target population

The primary contrasts were interface minus core and structured extramembrane minus core. Residues were matched within each protein in exact HELIX/SHEET strata and two-point pLDDT bins `[0,2), …, [98,100]`. Confidence multiplied by 100 was rounded to six decimals before assigning bin edges. Each retained cell required at least three residues on each side. Its common mass was `m = min(n_left, n_right)`. Within each side of a cell, residues received equal weight, and the two sides received the same total mass. Protein effects were cell-mass-weighted differences, normalized by that protein's total common mass; a protein required total mass at least 20 for a contrast. This coarsening uses matching mechanics related to coarsened exact matching [17], without claiming causal identification.

The target is the population retained by this matching and support rule. “Raw cell retention” counts the actual regional residues in retained cells of qualifying proteins, divided by all eligible residues in the two regions before cell/protein matching. “Effective mass retention” uses twice the total common mass in the numerator and the same denominator. That denominator already excludes non-HELIX/SHEET and depth-gap residues; it is not all 223,714 scored residues. Protein eligibility differs by contrast, so subtracting the two reported primary effects does not recover a third paired contrast on a common population. Discarded support is not extrapolated.

### Operational groups, curated homology and dependence blocks

The grouping universe was the union of 174 development and 554 validation records. The frozen procedure used full canonical sequences with Biopython 1.85 local alignment, BLOSUM62 and affine gap scores −10 and −0.5 [18]. Primary sequence edges required at least 30% identity, calculated with the alignment-column denominator, and paired coverage of at least 80% of each full sequence. A rescue permitted coverage of at least 50% on both sequences only with the same explicit UniProt whole-family annotation and an identical nonempty Pfam-ID complement. Shared Pfam IDs alone never established an edge. The complement representation does not encode domain order or repeat counts.

Connected components of qualifying edges formed **operational groups**. A protein was considered resolved if it carried an explicit whole-family annotation or participated in a qualifying edge; this does not mean its evolutionary family was independently proven. Unresolved validation proteins (53 in the whole cohort) were pooled into one additional group. Primary estimates first averaged protein effects equally within each group and then averaged groups equally. This weighting limits domination by heavily sampled operational groups, but splitting one evolutionary family into multiple groups can increase its weight. Pooling unresolved proteins changes the estimand and does not guarantee conservative interval coverage.

**Curated homology components** added exact normalized explicit whole-family-annotation links to the sequence-edge graph, retaining development and non-primary bridges. Family novelty required a resolved validation protein whose curated component did not touch development. Experimental-entry links were not used to establish homology. A prespecified sensitivity reweighted equal curated homology components instead of equal operational groups; neither term means a phylogenetically established family.

**Dependence blocks** joined operational groups through curated links and shared experimental entries. Whole blocks were resampled, keeping all associated group means together. The point estimate remained equal-operational-group weighted. A draw was the sum of all sampled group means divided by the number of sampled groups; this is a ratio estimator, not equal-block weighting. Resampling accounts for recorded dependence under the design, not every unknown evolutionary or experimental relationship [19].

### Uncertainty, materiality and sensitivities

The frozen bootstrap used 20,000 draws, NumPy PCG64 seed 20909 and sorted block identifiers. Linear percentile quantiles 0.0125 and 0.9875 gave 97.5% marginal intervals for the two primary contrasts. The Bonferroni allocation targets approximately 95% familywise coverage across those two contrasts, conditional on adequacy of the block bootstrap [19,20,21]. It does not establish exact finite-sample coverage. The secondary contrast used 95% intervals. The two focal contrasts in sensitivity analyses retained 97.5% marginal intervals; these are descriptive and are not simultaneously adjusted across scenarios.

The symmetric five-percentage-point margin was chosen by the investigator before endpoint computation as a planning and decision convention. The design recorded a rough precision illustration, not a formal power calculation or validation of biological importance, downstream utility or measurement error. We therefore interpret the continuous estimates and intervals first and report the unchanged margin rule second.

A material positive direction required the entire adjusted interval above +0.05 on the unit scale (+5 percentage points); a material negative direction required it below −0.05 (-5 percentage points). A contrast first required at least 30 proteins, 20 resolved groups and 20 resolved dependence blocks. The decision then considered the frozen method-sensitivity flag, followed by material-direction support. Seven prespecified sensitivities entered that flag: one- and five-point bins, family-novel, unresolved-excluded, broad grouping, curated-homology weighting and v0.1-compatible regions. A useful sensitivity had to shift by at least five points with a sign crossing, or reverse an established material direction, to trigger it. Thus no flag does not mean all reasonable methods are interchangeable.

Eleven sensitivity scenarios were retained: the preceding seven plus exact full-sequence construct, strict sequence grouping, mapped identity 100%, and experimental release after 30 April 2018. The v0.1-compatible scenario retained the confidence matching, left core/interface secondary structure unmasked and did not match secondary-structure categories, while requiring structured extramembrane residues. The temporal subset reduces one simple exposure proxy; release date alone cannot reconstruct model-training exposure. No formal equivalence test was specified or performed.

### Reproduction and descriptive case investigation

An independently coded statistical implementation within this AI-assisted project rebuilt matching cells, group components from recorded edges, weighted effects and primary block-bootstrap draws. It used compensated sums and independently implemented percentile interpolation, with the same specified random-number stream for numerical comparison. This tests arithmetic on the retained inputs; it does not provide independent data, external peer review or validation of every upstream source transformation.

The independently coded audit verifies both primary resampling calculations and all 36 point estimates: 12 scenarios including the primary specification, each with three contrasts. These are not 36 independent tests. It also checks subsets, group/block counts, support denominators and decision logic. Sensitivity intervals are reproduced by the original analysis; an additional independently coded resampling of all sensitivity intervals is not claimed. The initial graph audit rebuilt components from frozen edges rather than repeating every alignment or discovering all remote homologues. The complete acquisition and group-generation source is included in the revision addendum, and a bounded selected-pair check is described below.

Verification distinguishes compact-coordinate arithmetic for all 554 proteins, the earlier scalar score comparison on five proteins/2,341 residues, and earlier raw-column checks on 13 proteins/4,716 residues selected as eight large-discrepancy cases plus five accession-sorted cases. These raw cases are targeted, not a random validation sample. During revision, five additional covariate-selected cases/2,725 residues were checked against raw mmCIF columns, SIFTS mappings and OPM frame transfer. Selection used observability, mapping coverage, transform RMSD, alternate locations, insertion codes, chain identifiers and identity mismatch, without inspecting scores. Fifteen sequence pairs near recorded edge thresholds or involving unresolved proteins were realigned. The panel checks selected transformations, not complete graph discovery, all-source correctness or evolutionary family completeness (S2 Text).

### Post-outcome comparison with the AlphaFold endpoint

The unmodified published AlphaFold lddt function [22] was run with JAX 0.6.2 on CPU for all 554 proteins/223,714 paired-observed residues, separately in float64 and float32. Its strict distance/threshold comparisons and 1e-10 stabilizers differ from the inclusive study score. Inputs retained exactly the same observed intrachain coordinates and residue correspondence. Zero-masked padding to multiples of 128 reduced compilation shapes; padding outputs were discarded. This validation uses the pinned public implementation, not a reconstruction of the exact historical training environment. It is a direct score-function check, not a run of the AlphaFold prediction pipeline. NumPy ports, source hashes, masks, dtypes and primary-effect comparisons are retained in S2 Text and executable evidence. No alternative scores replaced the primary data.

### AI assistance and author responsibility

AI assistance was substantial in methodological discussion, software development, code review, literature-assisted drafting, manuscript preparation and reproducibility packaging. OpenAI ChatGPT substantially assisted methodological critique, literature-assisted reasoning and manuscript development; OpenAI Codex was used in software development, revision and verification. AI systems are not authors. Final scientific decisions, verification, interpretation and responsibility remain with the human author. The independently coded checks and author-requested adversarial audit are internal project activities, not independent external human peer review. The tools and documented environment are recorded with the revision; historical session-specific model identifiers that have not been provided are not inferred.

## Results

### Primary matched contrasts

The interface contrast retained 437 proteins, 260 resolved groups plus one pooled unresolved group, and 186 dependence blocks. The extramembrane contrast retained 229 proteins, 140 resolved groups plus one pooled group, and 95 blocks. The block counts equal the resolved-block counts here because the unresolved pool shares recorded dependence links with resolved groups; an extra group does not automatically add a separate block.

Table 1. Primary matched discrepancy contrasts

| Contrast | Effect [97.5% interval], percentage points | Proteins | Resolved / total groups | Blocks |
| --- | --- | --- | --- | --- |
| Interface − core | +0.583 [+0.234, +0.969] | 437 | 260 / 261 | 186 |
| Structured extra − core | +0.065 [-0.809, +0.940] | 229 | 140 / 141 | 95 |

The interface estimate was +0.583 percentage points (97.5% interval +0.234 to +0.969), supporting a small positive discrepancy contrast. The structured-extramembrane estimate was +0.065 (-0.809 to +0.940), with an interval spanning zero (Table 1; Figure 2). Both intervals lie within the investigator-chosen ±5-point range; the specified rule did not establish a departure beyond that margin. This was not a formal equivalence test and does not establish a biologically negligible effect.

### Support, balance and sensitivity

For interface-core, raw cell retention was 79.31% and effective mass retention 45.27%; for structured extra-core, these were 38.46% and 23.81%. These denominators refer to eligible regional residues before matching, not all 223,714 scored residues (Figure 3). Post-outcome description of the fixed weights shows that the primary targets are nearly entirely helical and strongly concentrated in high-confidence strata (Table 2). The presence of SHEET residues elsewhere in the cohort does not imply meaningful primary support for beta structure.

Table 2. Actual composition of the primary estimands

| Weight or matched-cell count | Interface-core | Structured extra-core |
| --- | --- | --- |
| HELIX weight, % | 99.9829 | 100.0000 |
| SHEET weight, % | 0.0171 | 0.0000 |
| pLDDT stratum at least 80, weight % | 98.6508 | 99.0648 |
| pLDDT stratum at least 90, weight % | 83.8070 | 80.5538 |
| HELIX matched cells | 2074 | 1116 |
| SHEET matched cells | 1 | 0 |

Weights apply the complete cell, protein and operational-group hierarchy. They are percentages of estimand weight, not percentages of all cohort residues. Strata at least 90 are contained within strata at least 80. These are post-outcome descriptive diagnostics of unchanged primary matching.

The maximum absolute protein-level residual confidence differences were 0.634 and 0.595 pLDDT points. Signed group-balanced differences were -0.098475 and -0.063043 on that scale, respectively. With the same fixed weights, the decomposition in Table 3 explains how this residual confidence term enters the discrepancy contrast. Coarsened matching improves comparability without making continuous confidence identical or controlling unmeasured factors.

Table 3. Additive decomposition using unchanged primary weights

| Component, percentage points | Interface-core | Structured extra-core |
| --- | --- | --- |
| Local agreement contrast | +0.484268 | +0.001759 |
| Normalized confidence contrast | -0.098475 | -0.063043 |
| Agreement minus confidence contrast | +0.582744 | +0.064802 |

The normalized-confidence component in percentage points is numerically equal to its pLDDT-point difference; the composite endpoint is not a pLDDT effect. The residual confidence term accounts algebraically for about 16.9% of the interface point estimate, not 16.9% of statistical bias or a causal mechanism. The weighted observed-neighbor contrasts were approximately -6.2 and -3.2 neighbors, respectively, emphasizing dependence on the observed intrachain reference context.

Table 4. Prespecified focal sensitivity contrasts

| Scenario | Interface-core, percentage points | Extra-core, percentage points |
| --- | --- | --- |
| Primary | +0.583 [+0.234, +0.969] | +0.065 [-0.809, +0.940] |
| 1-point bins | +0.572 [+0.237, +0.955] | -0.075 [-0.945, +0.836] |
| 5-point bins | +0.364 [+0.049, +0.702] | -0.051 [-0.898, +0.793] |
| v0.1-compatible regions | +0.344 [+0.012, +0.701] | +0.098 [-0.765, +0.927] |
| Family-novel subset | +0.807 [+0.412, +1.246] | +0.044 [-1.106, +1.181] |
| Exact full sequence | +0.260 [-0.350, +0.965] | +0.695 [-0.680, +2.177] |
| Unresolved excluded | +0.584 [+0.234, +0.971] | +0.038 [-0.827, +0.919] |
| Broad grouping | +0.688 [+0.264, +1.174] | +0.257 [-0.890, +1.356] |
| Strict sequence grouping | +0.562 [+0.223, +0.941] | +0.083 [-0.764, +0.955] |
| Mapped identity 100% | +0.600 [+0.189, +1.070] | +0.451 [-0.501, +1.429] |
| Equal curated homology | +0.615 [+0.215, +1.060] | +0.158 [-0.869, +1.158] |
| Released after 2018-04-30 | +0.656 [+0.232, +1.116] | +0.160 [-0.830, +1.168] |

All 11 prespecified sensitivities are displayed with the primary specification (Table 4; Figure 4). Family-novel matching retained 284/138 proteins; exact full-sequence matching retained 152/75. The latter gave an interface estimate of +0.260 [-0.350, +0.965] percentage points. Its interval crossing zero does not establish no effect or a difference from the primary population. Equal-curated-homology estimates were +0.615 and +0.158, reflecting a different weighting target. No prespecified threshold-and-sign flag fired; this supports the stability of that particular margin decision, not significance or invariance of the small shift under every method. Sensitivity intervals remain descriptive across scenarios.

Post-outcome leave-one-recorded-block diagnostics found maximum point shifts of 0.0748 and 0.1274 percentage points. The largest recorded blocks contained 11/261 and 8/141 primary groups. This provides no evidence that one recorded block dominates the estimates, but cannot exclude dependence missing from the graph. These block diagnostics are distinct from the earlier leave-one-operational-group calculations.

### Descriptive discrepancies and source contexts

The eight largest absolute unmatched protein-average discrepancies are shown in Table 5 and Figure 5. They are all-observed-residue means, a different estimand from the matched regional contrasts, and are selected descriptively without inferential intervals. Full sequence match in Table 5 means that the full experimental construct sequence equals the canonical sequence underlying the retrieved AFDB prediction; it does not describe a new construct-specific prediction campaign.

Table 5. Descriptive unmatched source cases

| Protein / entry | Mean discrepancy, percentage points | Experimental title | Full sequence match |
| --- | --- | --- | --- |
| P77335 / 6mrt | -38.677 | 12-meric ClyA pore complex | False |
| Q2FZP8 / 6s7v | -22.118 | Lipoteichoic acids flippase LtaA | False |
| G2QNH0 / 6gci | -21.119 | Structure of the bongkrekic acid-inhibited mitochondrial ADP/ATP carrier | False |
| Q6LPW0 / 7qha | +15.309 | Cryo-EM structure of the Tripartite ATP-independent Periplasmic (TRAP) transporter SiaQM from Photobacterium profundum in amphipol | True |
| G3IEF0 / 8djm | -14.277 | HMGCR-UBIAD1 Complex State 1 | False |
| P30878 / 7l17 | +13.763 | Crystal structure of sugar-bound melibiose permease MelB | False |
| A5F5Y6 / 8evu | +13.604 | Cryo EM structure of Vibrio cholerae NQR | True |
| A0KLE1 / 6h2f | -13.366 | Structure of the pre-pore AhlB of the tripartite alpha-pore forming toxin, AHL, from Aeromonas hydrophila. | False |

The experimental titles identify a 12-member ClyA pore for P77335/6mrt, a bongkrekic-acid-inhibited carrier for G2QNH0/6gci, an HMGCR–UBIAD1 complex state for G3IEF0/8djm, and a prepore assembly for A0KLE1/6h2f. These establish source context; they do not establish a particular alternative AFDB conformation or prove that assembly/state caused the discrepancy. Q2FZP8/6s7v has complete construct observability and 100% identity at mapped residues but differs from the full canonical sequence. Missing experimental coordinates therefore do not explain this particular large mean discrepancy. Conversely, Q6LPW0/7qha and A5F5Y6/8evu are exact full-sequence cases with positive discrepancies, so construct mismatch is not a universal explanation either. Raw source hashes, titles and source-publication DOIs are retained in the case table and original investigation.

### Reproduction outcome

The original numerical replay and independently coded primary arithmetic agreed within the specified absolute computational tolerance of 1e-12 on the unit scale. This tolerance describes software agreement, not experimental measurement precision. All 223,714 inclusive residue scores and neighbor counts also matched the author-requested audit exactly. In the direct published-function check, the maximum float64 residue difference was 1.41e-11 on the unit scale and primary-effect changes were below 1.4e-12 percentage points. In float32, 30 residues differed by more than 1e-4 on the unit scale; the maximum residue difference was 1.1364 percentage points. The primary effects became +0.582738882 and +0.064806934 percentage points, changes of -0.000005077 and +0.000005267. This supports numerical comparability for the stated mask without claiming universal per-residue or cross-renderer bitwise identity. The five additional raw-case checks and 15 selected alignment checks passed within their bounded scope.

## Discussion

The main finding is a +0.583-percentage-point interface contrast in agreement minus confidence, with a 97.5% interval from +0.234 to +0.969. The extramembrane contrast is close to zero but its interval includes both negative and positive small effects. These estimates concern predominantly high-confidence helical residues with common support. They constrain a coarse expectation of a large regional discrepancy within this selected population, while leaving low-confidence, beta-structured and unmatched residues largely unaddressed. Their practical meaning comes from the estimates, support and reference context; the fixed five-point convention does not establish biological relevance or irrelevance.

Within-protein matching reduces some differences in secondary structure and confidence composition. It also determines which residues can inform the comparison: more than 98.6% of weight comes from pLDDT strata at least 80 and effectively all weight is HELIX. Equal-group aggregation reduces domination by densely represented operational groups. The substantial loss of common support, especially outside the membrane, prevents extrapolation to all membrane proteins or all their residues. Residual confidence and observed-neighbor differences remain explicit features of the conditional comparison, not evidence of a causal membrane effect.

The grouping analysis makes dependence assumptions visible, rather than establishing evolutionary truth. Primary operational groups, curated homology components and resampling blocks have different functions. Equal curated-component weighting is a useful alternative target and is reported transparently; agreement in scale under this sensitivity is not proof that either grouping is complete. Unknown homology, shared experimental practices and annotation errors may remain.

Large descriptive discrepancies coexist with small matched regional averages. This is not a contradiction: unmatched case means and matched regional contrasts answer different questions. Experimental assembly and ligand metadata can suggest why canonical monomer products merit closer investigation, but causal mechanisms require additional evidence. Prior experimental comparisons emphasize the value and remaining limitations of predictions as structural hypotheses [23]. Neither sequence identity nor successful coordinate parsing establishes a matched biological state, while a large local-distance disagreement alone cannot identify the alternative state predicted by a model.

The analytical freeze provides a recorded separation between design and endpoint computation. Since coordinate-containing prediction files were already present, it cannot certify physical blinding or absence of all outcome access. The post-outcome checks strengthen implementation evidence without changing that historical claim or converting the study into externally preregistered research.

## Limitations

The observed-residue C-alpha endpoint shares its atom selection with the target of AlphaFold2 pLDDT. The direct implementation comparison is reassuring for these fixed inputs, but does not validate a full-canonical, all-atom or full-assembly endpoint. Missing residues, truncations and interchain contacts can change the neighborhood and the relevant structural context. OPM boundaries are computed placements rather than measurements of lipid contacts. The capped available-structure cohort and nearly all-helical high-confidence matching support are not a probability sample of membrane-protein biology.

The cohort excludes development proteins but does not prove absence of training, template or homolog exposure. The family-novel and temporal analyses only address recorded operational proxies. Construct, ligand, oligomeric and conformational differences remain possible. Group annotations, sequence thresholds and unresolved pooling affect weights and dependence assumptions. Percentile block intervals have approximate coverage and do not encompass every source of structural or annotation uncertainty.

The five-point margin has no externally validated biological or utility basis, and no formal equivalence design was used. Sensitivity intervals lack a simultaneous across-scenario guarantee. Earlier raw checks comprised 13 targeted proteins; the added five covariate-selected cases broaden the inspected edge cases but do not validate every source mapping. Selected alignment replay does not establish complete evolutionary relationships. The independent arithmetic shares retained data and graph assumptions with the primary calculation. Some limitations concern packaging and transparency; unobserved biological context, exposure uncertainty and restricted support remain properties of the study design.

## Data and code availability

The original content-addressed RC1 contains the fixed compact inputs, numerical results and executable numerical replay. The revision addendum provides the revised manuscript, post-outcome diagnostics and the minimal acquisition, mapping and group-generation sources. S1 Text and S1 Dataset retain the original scientific definitions and numerical evidence; S2 Text and S2 Dataset document this revision. Numerical replay from compact inputs is offline and does not need private Git history. Reconstructing upstream inputs separately requires the identified source objects, network access or a matching cache, and source-specific permissions. Historical hashes identify required bytes but do not guarantee their continued availability.

The original command python -m membranecal_release.verify reproduces the immutable RC1; it does not by itself verify this revised manuscript or rebuild all upstream resources. Revision-specific checks and source-workspace instructions are included in the addendum. A Zenodo deposit is planned, but no final public or reviewer-access URL or DOI is available yet. Author-created code is designated Apache-2.0 and authored text/figures CC BY 4.0, as selected by the author. Third-party data retain source-specific terms. Explicit permission from Andrei Lomize covers redistribution of the numerical OPM values described for this reproducibility package, including the required membrane orientation, depth and thickness information, with attribution and transformation provenance retained. This permission does not extend to unrelated raw OPM assets or arbitrary database redistribution. Public release and submission access arrangements are therefore not represented as complete.


## Figure legends

**Figure 1.** Cohort selection and primary matching support. Prefilters follow code order: development-entry overlap, type/resolution availability, then resolution/boundary. Entry counts and protein-case counts are distinguished. Primary matching retains 437 interface-core and 229 structured-extramembrane-core proteins, with 45.27% and 23.81% effective mass retention relative to eligible regional residues. This is a capped, development-disjoint sample; absence of training exposure is not established.

**Figure 2.** Equal-operational-group paired effects, with 97.5% marginal block-percentile intervals. Left panel shows the frozen ±5-point scale; right panel enlarges the same estimates. Approximate Bonferroni familywise coverage applies to the two primary contrasts. No equivalence band is asserted. Discrepancy units are percentage points. Source: tables/primary_effects.csv.

**Figure 3.** Common support and residual confidence balance. Panel a shows the percentage of eligible regional residues in retained cells (filled circles) and effective common mass (open diamonds); N is the pre-match eligible-residue denominator. Panel b shows the signed group-balanced confidence difference (filled circles) and maximum absolute protein-level difference (open triangles), in pLDDT points. These are descriptive balance diagnostics, not effect intervals. Source: tables/support_balance.csv.

**Figure 4.** All 11 prespecified sensitivity scenarios alongside primary estimates. Intervals for these two contrasts are 97.5% marginal; sensitivity intervals are descriptive and not simultaneously adjusted across scenarios. The ±5-point margin lies outside the zoomed plotting range. n is matched protein count; contrast units are percentage points. Source: `tables/sensitivity_effects.csv`.

**Figure 5.** Eight post-outcome cases selected by largest absolute all-observed-residue protein mean discrepancy. Markers show descriptive unmatched means in percentage points; connecting lines extend to zero and do not represent uncertainty intervals. Negative values denote confidence exceeding agreement for this endpoint. No causal state-mechanism claim is made. Context and original publication DOIs: tables/descriptive_cases.csv.

## References

1. Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. *Nature*. 2021;596:583–589. [doi:10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2).

2. Varadi M, Anyango S, Deshpande M, et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. *Nucleic Acids Research*. 2022;50:D439–D444. [doi:10.1093/nar/gkab1061](https://doi.org/10.1093/nar/gkab1061).

3. Varadi M, Bertoni D, Magana P, et al. AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. *Nucleic Acids Research*. 2024;52:D368–D375. [doi:10.1093/nar/gkad1011](https://doi.org/10.1093/nar/gkad1011).

4. Bertoni D, Tsenkov M, Magana P, et al. AlphaFold Protein Structure Database 2025: a redesigned interface and updated structural coverage. *Nucleic Acids Research*. 2026;54:D358–D362. Published online 22 November 2025. [doi:10.1093/nar/gkaf1226](https://doi.org/10.1093/nar/gkaf1226).

5. Mariani V, Biasini M, Barbato A, Schwede T. lDDT: a local superposition-free score for comparing protein structures and models using distance difference tests. *Bioinformatics*. 2013;29:2722–2728. [doi:10.1093/bioinformatics/btt473](https://doi.org/10.1093/bioinformatics/btt473).

6. Lomize MA, Pogozheva ID, Joo H, Mosberg HI, Lomize AL. OPM database and PPM web server: resources for positioning of proteins in membranes. *Nucleic Acids Research*. 2012;40:D370–D376. [doi:10.1093/nar/gkr703](https://doi.org/10.1093/nar/gkr703).

7. Hegedűs T, Geisler M, Lukács GL, Farkas B. Ins and outs of AlphaFold2 transmembrane protein structure predictions. *Cellular and Molecular Life Sciences*. 2022;79:73. [doi:10.1007/s00018-021-04112-1](https://doi.org/10.1007/s00018-021-04112-1).

8. Dobson L, Szekeres LI, Gerdán C, Langó T, Zeke A, Tusnády GE. TmAlphaFold database: membrane localization and evaluation of AlphaFold2 predicted alpha-helical transmembrane protein structures. *Nucleic Acids Research*. 2023;51:D517–D522. [doi:10.1093/nar/gkac928](https://doi.org/10.1093/nar/gkac928).

9. Jambrich MA, Tusnády GE, Dobson L. How AlphaFold2 shaped the structural coverage of the human transmembrane proteome. *Scientific Reports*. 2023;13:20283. [doi:10.1038/s41598-023-47204-7](https://doi.org/10.1038/s41598-023-47204-7).

10. Xie T, Huang J. Can Protein Structure Prediction Methods Capture Alternative Conformations of Membrane Transporters? *Journal of Chemical Information and Modeling*. 2024;64:3524–3536. [doi:10.1021/acs.jcim.3c01936](https://doi.org/10.1021/acs.jcim.3c01936).

11. Edmunds NS, Genc AG, McGuffin LJ. Benchmarking of AlphaFold2 accuracy self-estimates as indicators of empirical model quality and ranking: a comparison with independent model quality assessment programmes. *Bioinformatics*. 2024;40:btae491. [doi:10.1093/bioinformatics/btae491](https://doi.org/10.1093/bioinformatics/btae491).

12. Mahtha SK, Venkadesan S, Mohanty D. Comparative evaluation of the prediction accuracy of AlphaFold and ESMFold for monomeric and dimeric proteins. *NAR Genomics and Bioinformatics*. 2026;8:lqag002. [doi:10.1093/nargab/lqag002](https://doi.org/10.1093/nargab/lqag002).

13. AlQuraishi M. ProteinNet: a standardized data set for machine learning of protein structure. *BMC Bioinformatics*. 2019;20:311. [doi:10.1186/s12859-019-2932-0](https://doi.org/10.1186/s12859-019-2932-0).

14. Velankar S, Dana JM, Jacobsen J, et al. SIFTS: Structure Integration with Function, Taxonomy and Sequences resource. *Nucleic Acids Research*. 2013;41:D483–D489. [doi:10.1093/nar/gks1258](https://doi.org/10.1093/nar/gks1258).

15. Dana JM, Gutmanas A, Tyagi N, et al. SIFTS: updated Structure Integration with Function, Taxonomy and Sequences resource allows 40-fold increase in coverage of structure-based annotations for proteins. *Nucleic Acids Research*. 2019;47:D482–D489. [doi:10.1093/nar/gky1114](https://doi.org/10.1093/nar/gky1114).

16. EMBL-EBI. [AlphaFold database release notes](https://www.ebi.ac.uk/pdbe/news/alphafold-database-release-notes), 21 October 2025. Database release v6 terminology; accessed 9 September 2026.

17. Iacus SM, King G, Porro G. Causal Inference without Balance Checking: Coarsened Exact Matching. *Political Analysis*. 2012;20:1–24. [doi:10.1093/pan/mpr013](https://doi.org/10.1093/pan/mpr013). Cited for matching mechanics, not causal identification in this study.

18. Biopython. [Pairwise sequence alignment, version 1.85](https://biopython.org/docs/1.85/Tutorial/chapter_pairwise.html). Accessed 9 September 2026.

19. Field CA, Welsh AH. Bootstrapping clustered data. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*. 2007;69:369–390. [doi:10.1111/j.1467-9868.2007.00593.x](https://doi.org/10.1111/j.1467-9868.2007.00593.x).

20. Dunn OJ. Multiple Comparisons among Means. *Journal of the American Statistical Association*. 1961;56:52–64. [doi:10.1080/01621459.1961.10482090](https://doi.org/10.1080/01621459.1961.10482090).

21. NIST/SEMATECH. [Bonferroni's method](https://www.itl.nist.gov/div898/handbook/prc/section4/prc473.htm), e-Handbook of Statistical Methods. Accessed 9 September 2026.

22. Google DeepMind. AlphaFold lddt score function. Revision dbaafbcdea0cf39fabe502928cab55754c1d5dc7, alphafold/model/lddt.py. [Pinned source](https://github.com/google-deepmind/alphafold/blob/dbaafbcdea0cf39fabe502928cab55754c1d5dc7/alphafold/model/lddt.py). Accessed 15 September 2026.

23. Terwilliger TC, Liebschner D, Croll TI, et al. AlphaFold predictions are valuable hypotheses and accelerate but do not replace experimental structure determination. *Nature Methods*. 2024;21:110–116. [doi:10.1038/s41592-023-02087-4](https://doi.org/10.1038/s41592-023-02087-4).

## Supporting information captions

**S1 Text. Scientific supplement.** Original eligibility, mapping, grouping, common-support and sensitivity definitions. Retained from RC1; revision clarifications are in S2 Text.

**S1 Dataset. Original numerical evidence.** Nine original CSV tables and README, including all 36 scenario contrasts.

**S2 Text. Post-outcome validation and revision evidence.** Prior-study comparison, predictor-product metadata, endpoint comparisons, bounded raw/grouping QA, analytical chronology, source-reconstruction scope and audit-response ledger.

**S2 Dataset. Revision diagnostics and source code.** Machine-readable population/decomposition, endpoint and bounded-QA reports, executable diagnostic code, and minimal acquisition/mapping/group-generation sources with checksums. No raw source archives are included.

## Declarations

**Funding.** The author received no specific funding for this work.

**Competing interests.** The author develops MembraneState Bench, the independent research and benchmarking programme within which this study was conducted, and is exploring potential future commercial applications of related benchmarking infrastructure. No commercial entity funded this study.

**Author contributions.** Pavel Trofimchik: Conceptualization; Data curation; Formal analysis; Investigation; Methodology; Project administration; Resources; Software; Validation; Visualization; Writing - original draft; Writing - review and editing.

**Ethics and data reuse.** This study is a computational reuse of existing protein-structure and annotation resources; no new human-participant or animal experiments were conducted for this analysis. Source attribution and redistribution conditions are recorded separately.
