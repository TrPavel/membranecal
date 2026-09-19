# Acquisition and grouping source

These three scripts are byte-identical copies of the scientific acquisition, freeze and methods sources used in the original study. They contain the download, alias/exclusion, SIFTS mapping, entity/asym selection, occupancy/alternate-location selection, membrane-frame transfer, sequence alignment, edge discovery, matching and aggregation implementations. Their hashes are in ../evidence/source_bindings.json. They import only standard Python, NumPy, Gemmi, PyArrow and Biopython, and these copied scripts. They do not import the private project package or another worktree.

Numerical replay from the immutable RC1 is the tested primary reproduction route. Reconstructing source transformations is a separate task requiring upstream network access, source-specific permissions and historical source bytes. Current upstream responses can differ, disappear or introduce parsing changes. A current rerun must not silently replace the original cohort or be described as the historical analytical freeze.

## Separate source workspace

Use Python 3.12 with the RC1 documented dependencies plus Gemmi 0.7.3. The historical alignment version is Biopython 1.85. From this directory run:

```
python scaffold.py --release-root PATH_TO_ORIGINAL_RC1 --output PATH_TO_NEW_SOURCE_WORKSPACE
```

This validates every RC1 payload hash and prepares only a new self-contained directory. It seeds development identifiers and construct sequences from RC1, retaining the original exclusions. It does not download or run scientific analysis. Change directory to the new source workspace; the original modules can then be inspected or imported using `from scripts import membranecal_v02_acquire, membranecal_v02_freeze`.

For current-source reconstruction, the original acquisition command `python -m scripts.membranecal_v02_acquire discover --limit 800` performs index selection, and `python -m scripts.membranecal_v02_acquire acquire --workers 4` prepares structures. Check `--help` for the exact command arguments. This is not a claim that today's index recreates the original selection. To recreate the original selection, obtain the index with the RC1 discovery SHA256 and the individually recorded source objects; arrange the source cache with its URL/hash metadata as implemented in `_get`. Compare each source to `expected_source_hashes.json` before using it. The original full candidate annotations were intentionally not redistributed; they must come from that identified historical index.

The scaffold uses the retained development metadata. `prepare_development` would retrieve current UniProt records and is therefore a distinct source refresh, not the retained historical input. `build_groups(records, old)` in the copied freeze module is directly callable with `validation_records.json` and `workspace/v02/development_sequences.json`; it generates sequence evidence and grouping from those frozen sequences without scores. It can be expensive because it visits all pairs. The added bounded audit replays 15 covariate-selected pairs; a new exhaustive search for remote homology is not claimed.

The historical `freeze` command also expects a known-family validation panel, historical prepared records and provenance/design files. It is included as inspectable original source, not advertised as a one-command new valid freeze. Its dictionary construction, matching and group-generation functions expose all transformations. No source-reconstruction output can inherit the original private Git SHA or its chronology. New outputs belong in the new workspace and require a separate comparison/correction record.

## Data and rights

The source code is authored project software licensed under the owner-selected Apache-2.0 terms. The owner-selected licence does not grant additional rights to downloaded PDB, SIFTS, UniProt, AFDB or OPM objects. Explicit permission covers the described numerical OPM values for this MembraneCal package; attribution and provenance must be retained. It does not cover unrelated raw assets or arbitrary database redistribution. Raw caches are not shipped in this addendum. Names such as `raw_path` inside the code are relative runtime output fields, not leaked cache contents or a requirement to access the author's local files.
