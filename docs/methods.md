# Frozen methods

Scientific identity: **MembraneCal study v0.2**. The authoritative executable
definition is [protocol.json](../publication/v0.2/code/numerical/data/protocol.json),
with the [original design](../publication/v0.2/provenance/original_design.md).

The capped cohort includes 554 proteins and 223,714 paired-observed residues.
Construct/canonical mapping uses SIFTS with explicit identity, observed-coverage
and representative-selection rules. OPM-derived membrane frames define the
operational regions. The primary estimand is a regional difference in local
C-alpha structural agreement minus pLDDT/100, on high-confidence HELIX support.

Within-protein coarsened exact matching uses 2-point pLDDT bins and secondary
structure. Matched-cell common mass is shared across regions. Aggregation uses
equal proteins within operational sequence/rescue groups, then equal group means;
this is not equal weighting of curated biological families. Dependence-aware
bootstrap blocks follow the frozen sequence/family rules.

The two primary intervals are 97.5% (Bonferroni adjustment for two contrasts),
with 20,000 draws for each. The ±5 percentage-point margin is investigator-chosen,
not a validated clinical/biophysical threshold or formal equivalence test.
All 36 frozen scenario point estimates are retained; supplementary post-outcome
diagnostics are labelled separately and do not redefine the primary analysis.

The endpoint uses intrachain paired C-alpha positions, a 15 Å reference-neighbour
radius and distance thresholds 0.5, 1, 2 and 4 Å, without superposition. It is not
an all-atom, interchain or functional endpoint. Missing coordinates are not imputed.

No scientific function, frozen input, estimate, interval or figure is changed by
repository packaging. Source/freeze identities and SHA256 bindings are in
[frozen-files.json](../release/frozen-files.json) and the payload provenance.
