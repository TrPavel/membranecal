# MembraneCal v0.2: outcome-blind design rationale

Written during acquisition, before any v0.2 local-accuracy computation. Parent v0.1
is development evidence only. Scope is alpha-helical transmembrane proteins with
public canonical AFDB products, not all membrane proteins or all predictions.

## Endpoint and estimands

The primary ecological estimand compares observed experimental chains with canonical
AFDB monomer predictions. Exact full construct/canonical sequence identity defines
the available exact-construct subset; it does not establish matching ligands,
oligomers or state. No new prediction campaign is commissioned.

The inherited C-alpha distance agreement is superposition-free, using only paired
observed residues and reference neighbours within 15 Angstrom, inclusive thresholds
0.5, 1, 2, 4 Angstrom, no sequence-separation exclusion, no imputation of missing
coordinates. This is a C-alpha endpoint, not canonical all-atom lDDT. A separately
written scalar implementation will check a bounded subset after unseal.

## Regions and matching

OPM's hydrophobic slab is a reproducible computed orientation, not experimental
membrane positioning. Retain the development boundaries: |z|/half-thickness <=0.8
core, (0.8,1.2] interface, >1.6 extramembrane; the intervening band is excluded.
There is no pre-outcome external basis here for fitting new boundaries. Restrict all
three primary regions to reference mmCIF HELIX or SHEET, and match that category
exactly. Missing annotation is unavailable, not disorder. The v0.1-compatible
sensitivity restores unmasked core/interface and omits secondary-structure matching.

Within each protein use fixed 2-pLDDT-point bins and exact secondary-structure
categories. Require three residues on each side of a cell and total common mass
20 per contrast. Cell mass is min(n_left,n_right); give both sides that same mass,
with equal residue weights within each side. This controls coarsened support, not
continuous equality; report residual confidence differences. Freeze 1- and 5-point
bins as sensitivities. No extrapolation outside retained cells.

## Family and dependence

Do not equate shared Pfam with evolutionary independence. Compare full canonical
sequences over the union of development and validation proteins, retaining bridges.
Use Biopython 1.85 Smith-Waterman/BLOSUM62, affine gaps -10/-0.5. Primary edges
require >=30% actual identity over alignment columns and >=80% paired-residue
coverage of each entire sequence. The aligner has no calibrated E-value; none is
claimed. Compared with MMseqs2 this slower implementation is locally reproducible
without a new platform-specific binary. The actual-identity/full-coverage principle
is informed by the [MMseqs2 guide](https://mmseqs.com/latest/userguide.pdf).

Remote-family rescue additionally requires >=30% identity, >=50% paired coverage
of both sequences, the same explicit UniProt 'Belongs to ... family' annotation,
and the same nonempty complete Pfam ID complement. An ID complement is not an
ordered domain architecture: repeated/domain order uncertainty is retained as a
limitation. No common single domain alone creates a rescue edge. Validate AQP1/AQP2,
ADRB1/ADRB2, AQP2/ADRB2 negative, fragment and multi-domain adversarial examples.
Use connected components, not assumed pairwise independent families. Isolated
proteins without curated family evidence are unresolved and conservatively pooled
into one primary uncertainty group. Sensitivities omit unresolved proteins and
merge groups sharing a curated family/OPM family to probe residual remote homology.
Primary uncertainty additionally links components with the same normalized,
explicit UniProt whole-family label, retaining development and nonprimary bridges.
Family-novel sensitivity uses sequence/rescue plus these curated homology edges,
excludes any such component touching development and all unresolved cases, and
does not use shared experimental entries as evidence of homology.
'Novel' is relative to this operational grouping, not proof
of previously unseen evolutionary ancestry.

## Inference, sample size and decision

Two primary contrasts: interface-core and extramembrane-core. Extra-interface is
secondary. Equal proteins within operational group, equal operational-group means.
Thus a curated family split into more groups receives more point-estimate weight;
an equal-curated-homology-component sensitivity explicitly measures this choice.
Bootstrap 20,000 fixed-seed samples of connected sequence/rescue-group,
curated-whole-family and experimental-entry blocks, retaining
all proteins within a drawn block. This preserves both family and shared-entry
dependence. Two 97.5% marginal percentile intervals implement Bonferroni nominal
95% familywise coverage; secondary intervals are 95%. Finite-cluster percentile
coverage is approximate, not guaranteed.

Planning pool: 800 OPM representatives from a broader metadata census. Anticipate
40-60% attrition and 25-60% within-contrast support retention; these are planning
ranges, not results or a formal power calculation. Minimum useful sample per
primary contrast: 30 paired proteins, 20 resolved operational families and 20
independent resampling blocks. Aim for 40 groups. At SD 0.10, 20 independent groups
give order-of-magnitude 97.5% half-width 0.05; 40 give 0.035. Precision may be worse.
Use a material difference of 0.05 endpoint units (5 points) in either direction,
chosen symmetrically before outcomes; it is NOT a formal equivalence test margin.

Decision order: if either primary lacks minimum useful support, report independent
cohort not feasible for the full planned question while still reporting frozen
estimates. Otherwise material support requires a primary adjusted interval entirely
above +0.05 or below -0.05. If primary 1/5-bin, family-novel (when useful),
unresolved-excluded, broad-family, curated-homology or v0.1-compatible sensitivity shifts a primary
estimate by >=0.05 and crosses zero, or reverses an established material direction,
method sensitivity dominates. Otherwise report supported material regional effect
or no material regional effect established. No-effect-established is not equivalence.

References materially informing matching and multiplicity:
[Iacus, King & Porro 2012](https://doi.org/10.1093/pan/mpr013),
[NIST Bonferroni guidance](https://www.itl.nist.gov/div898/handbook/prc/section4/prc473.htm).
Confidence is a model output and local agreement a benchmark proxy; their numerical
difference is operational calibration, not a universal probabilistic calibration.

## Acquisition constraints

Raw public files and source-response hashes remain in this worktree's ignored
workspace. Commit compact transformed mappings/metadata only. Inherit the
artefact-specific OPM licence decision from v0.1, retain its evidence hash, and
attribute every source. RCSB CC0 and EMBL-EBI terms rechecked 2026-09-09.
UniProt and AFDB attribution/licence decisions are inherited from the immediately
preceding v0.1 audit. No prohibited bulk source archives enter Git.

Implementation semantics are checked against the pinned
[Biopython 1.85 pairwise-alignment documentation](https://biopython.org/docs/1.85/Tutorial/chapter_pairwise.html).
The temporal diagnostic uses the original 2018-04-30 structure cutoff stated in
the [official AFDB FAQ](https://www.alphafold.ebi.ac.uk/faq). It does not identify
training membership of any specific protein or establish absence of sequence,
template, later-model or homolog exposure. Individual AFDB product versions and
model dates are retained in the cohort.
