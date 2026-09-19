# Source terms and provenance

Current OPM permission is described in ../THIRD_PARTY_NOTICES.md and supersedes any historical permission uncertainty below. These historical source observations do not create a current deposit gate.

# Source-specific data notices

No blanket project-code relicensing of third-party or derived data is asserted.

PDB: reference_ca coordinates, experimental identifiers, constructs, secondary structure and core experimental metadata derive from PDBx/mmCIF; CC0-origin. Policy: https://www.rcsb.org/pages/usage-policy . Attribution to wwPDB/PDB and source structure publications is retained as scientific provenance.

UniProt: canonical sequences, accession aliases, family phrases and Pfam cross-reference identifiers in proteins.json, development_sequences.json and groups.json derive from UniProt. Copyrightable data are CC BY 4.0: https://www.uniprot.org/help/license/ . Attribute the UniProt Consortium; retain original source URLs, retrieval dates and hashes. No separate Pfam archive is distributed.

AlphaFold DB: scores.parquet retains predicted C-alpha coordinates and calculated local agreement; residues_preoutcome.parquet retains pLDDT confidence. Predictions and confidence originate in AFDB models, CC BY 4.0: https://alphafold.ebi.ac.uk/assets/License-Disclaimer.pdf . Attribute the AlphaFold DB team/Google DeepMind/EMBL-EBI and cited model/database publications. Computed local agreement and statistical summaries are study transformations, not unmodified database outputs.

OPM: residues_preoutcome.parquet retains absolute membrane-frame depth in angstroms, continuous normalized membrane depth and derived region labels; proteins.json retains half-thickness, alignment RMSD, matched atom counts and OPM family/experimental metadata. discovery.json contains the selected candidate identifiers and prefilter counts only; broad copied annotations are omitted. Reference C-alpha coordinates remain the PDB reference frame, not an OPM oriented-coordinate archive. Upstream oriented coordinates were used to transfer the membrane frame. The current OPM site footer links CC BY 3.0; attribute The Regents of the University of Michigan and OPM, Lomize and colleagues, https://opm.phar.umich.edu/ and https://doi.org/10.1093/nar/gkr703 . Indicate filtering, coordinate-frame transfer and normalized-depth calculation. Explicit permission from Andrei Lomize covers the described numerical values for this MembraneCal package. The permission note in ../THIRD_PARTY_NOTICES.md defines its limits; no raw OPM archive is included.

SIFTS: mapped canonical positions and mapping QC derive from SIFTS; raw XML is omitted. https://www.ebi.ac.uk/about/terms-of-use/ applies alongside source-specific upstream rights; attribute EMBL-EBI/SIFTS and the cited SIFTS publications. This is narrower source-specific treatment, not an invented blanket CC0 grant.

Source URLs, original raw checksums, dates and transformations are retained in data/sources.json. No raw source archives are redistributed. Replacing all OPM depth values by live retrieval would lose offline replay and depend on mutable upstream availability; keeping the frozen transformed values with transparent qualification preserves auditability. OPM numerical redistribution permission is received; the final public deposit still requires owner approval.
