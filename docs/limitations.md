# Interpretation boundaries

- Primary matched support consists of high-confidence HELIX residues. The total
  accepted cohort is not the primary matched sample size, and results do not
  generalize automatically to low-confidence, disordered or non-helical residues.
- The C-alpha local-distance endpoint measures intrachain agreement with a selected
  experimentally observed structure. It does not validate atomic chemistry,
  multimeric interfaces, membrane insertion, stability, function or state populations.
- Experimental constructs and canonical predictions differ. Mapping/coverage rules
  reduce mismatch but cannot remove every construct, ligand or state ambiguity.
- Exposure controls are exposure-aware, not proof of absence of training leakage.
  Operational sequence/rescue groups are not universally curated family definitions.
- The capped source population and exclusions limit representativeness. OPM
  orientation is source-derived annotation, not experimental membrane-position truth.
- A small positive interface contrast and an interval spanning zero for the
  extramembrane contrast do not establish universal calibration or equivalence.
- The protocol was prospectively frozen in local Git before recorded unsealing.
  This is not independently timestamped external preregistration.
- Compact inputs enable numerical replay. Historical raw-source reconstruction is
  separate and requires the exact historical objects and their source permissions.
- Figure pixel identity depends on the renderer and fonts; supported numerical
  comparisons use the frozen tolerance. No regulatory use is established.

See the [manuscript](../publication/v0.2/manuscript/manuscript.md) for full limitations.
