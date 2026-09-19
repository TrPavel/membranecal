"""Outcome-agnostic case selection for a bounded post-outcome upstream QA panel."""
import argparse
import gzip
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import gemmi
import numpy as np
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'upstream'))
from scripts.membranecal_v02_methods import grouping_edge, sequence_evidence


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--release-root', type=Path, required=True)
    parser.add_argument('--raw-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if sys.flags.optimize:
        raise RuntimeError('Assertions must remain enabled.')
    read = lambda name: json.loads((args.release_root / 'data' / name).read_text())
    frozen_source_hashes = {r['sha256'] for r in read('sources.json')['records']}
    proteins = read('proteins.json')['records']
    byacc = {p['accession']: p for p in proteins}
    residues = pq.read_table(args.release_root / 'data/residues_preoutcome.parquet').to_pylist()
    byres = defaultdict(list)
    for row in residues:
        byres[row['accession']].append(row)
    selection = {}
    for key, reverse in [('observed_fraction', False), ('mapping_coverage', False), ('opm_transform_rmsd', True)]:
        selection[key] = sorted(proteins, key=lambda p: ((-1 if reverse else 1) * p[key], p['accession']))[0]['accession']
    features = {
        'alternate_location': sorted({r['accession'] for r in residues if r['altloc']}),
        'insertion_code': sorted({r['accession'] for r in residues if re.search('[A-Za-z]', r['author_residue'])}),
        'author_label_chain_difference': sorted(p['accession'] for p in proteins if p['author_chain_id'] != p['label_asym_id']),
        'mapped_residue_identity_mismatch': sorted({r['accession'] for r in residues if not r['identity_match']}),
    }
    for label, candidates in features.items():
        selection[label] = candidates[0] if candidates else None
    groups = read('groups.json')
    edges = groups['edges']
    chosen = []
    for method in sorted({e['method'] for e in edges}):
        subset = [e for e in edges if e['method'] == method]
        for key, boundary in [('identity', .3), ('coverage_a', .8 if method == 'FULL_LENGTH_SEQUENCE' else .5)]:
            chosen.extend(sorted(subset, key=lambda e: (abs(e[key] - boundary), e['a'], e['b']))[:3])
    pairs = sorted({(e['a'], e['b']) for e in chosen})
    all_records = {**{p['accession']: p for p in read('development_sequences.json')['records']}, **byacc}
    unknown = sorted(a for a, value in groups['resolved_by_accession'].items() if not value)
    for a in unknown[:3]:
        b = min((b for b in all_records if b != a), key=lambda b: (abs(len(all_records[a]['canonical_sequence']) - len(all_records[b]['canonical_sequence'])), b))
        pairs.append(tuple(sorted((a, b))))
    pairs = sorted(set(pairs))
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'bounded_selection.json').write_text(json.dumps({'selection_basis': 'Frozen covariates and graph evidence only; no scores used', 'raw_cases': selection, 'group_pairs': pairs, 'label': 'post-outcome validation'}, indent=2))
    raw_results = []
    for acc in sorted({a for a in selection.values() if a}):
        protein = byacc[acc]
        entry = protein['entry_id'].lower()
        paths = [args.raw_root / 'pdb' / (entry + '.cif.gz'), args.raw_root / 'sifts' / (entry + '.xml.gz'), args.raw_root / 'opm' / (entry + '.pdb')]
        hashes = {}
        for path in paths:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            meta = json.loads(Path(str(path) + '.source.json').read_text())
            assert digest == meta['sha256']
            assert digest in frozen_source_hashes, (path.name, digest)
            hashes[path.name] = digest
        block = gemmi.cif.read(str(paths[0])).sole_block()
        atoms = block.get_mmcif_category('_atom_site.')
        selected = {}
        names = {}
        for i, atom in enumerate(atoms['label_atom_id']):
            if atom != 'CA' or atoms['label_asym_id'][i] != protein['label_asym_id'] or str(atoms['pdbx_PDB_model_num'][i]) != protein['model_id']:
                continue
            if atoms['label_seq_id'][i] in (False, None, '.', '?') or float(atoms['occupancy'][i]) <= 0:
                continue
            assert atoms['label_entity_id'][i] == protein['entity_id']
            assert atoms['auth_asym_id'][i] == protein['author_chain_id']
            label = int(atoms['label_seq_id'][i])
            ins = atoms['pdbx_PDB_ins_code'][i]
            author = str(atoms['auth_seq_id'][i]) + (str(ins) if ins not in (False, None, '.', '?') else '')
            alt = atoms['label_alt_id'][i] or ''
            rank = (-float(atoms['occupancy'][i]), str(alt))
            xyz = np.array([float(atoms[k][i]) for k in ['Cartn_x', 'Cartn_y', 'Cartn_z']])
            item = (rank, xyz, author, atoms['label_comp_id'][i])
            if label not in selected or rank < selected[label][0]:
                selected[label] = item
        mapping = defaultdict(set)
        for residue in ET.fromstring(gzip.decompress(paths[1].read_bytes())).iter():
            if residue.tag.split('}')[-1] != 'residue':
                continue
            cross = [e.attrib for e in residue if e.tag.split('}')[-1] == 'crossRefDb']
            pdb = [e for e in cross if e.get('dbSource') == 'PDB' and e.get('dbChainId') == protein['author_chain_id']]
            uni = [e for e in cross if e.get('dbSource') == 'UniProt']
            if len(pdb) == len(uni) == 1:
                mapping[pdb[0].get('dbResNum')].add((uni[0]['dbAccessionId'], int(uni[0]['dbResNum'])))
        op = gemmi.read_structure(str(paths[2]))
        opca = {}
        for chain in op[0]:
            if chain.name != protein['author_chain_id']:
                continue
            for residue in chain:
                ca = [a for a in residue if a.name == 'CA']
                if ca:
                    a = max(ca, key=lambda a: a.occ)
                    opca[str(residue.seqid)] = (np.array([a.pos.x, a.pos.y, a.pos.z]), residue.name)
        hits = [item for item in selected.values() if item[2] in opca and item[3] == opca[item[2]][1]]
        moving = np.array([x[1] for x in hits]); target = np.array([opca[x[2]][0] for x in hits])
        mc, tc = moving.mean(0), target.mean(0)
        u, _, vt = np.linalg.svd((moving - mc).T @ (target - tc))
        rotation = u @ np.diag([1, 1, np.linalg.det(u @ vt)]) @ vt
        max_depth = 0
        for row in byres[acc]:
            item = selected[row['label_seq_id']]
            assert item[2] == row['author_residue']
            assert np.array_equal(item[1], row['reference_ca'])
            assert mapping[item[2]] == {(protein.get('sifts_accession', acc), row['canonical_position'])}
            depth = abs(((item[1] - mc) @ rotation + tc)[2])
            max_depth = max(max_depth, abs(depth - row['depth']))
        assert max_depth < 1e-8
        raw_results.append({'accession': acc, 'case_id': protein['case_id'], 'residues': len(byres[acc]), 'raw_hashes': hashes, 'mapping_and_coordinates': 'PASS', 'max_depth_difference_A': float(max_depth), 'matched_opm_atoms': len(hits)})
    edge_by = {(e['a'], e['b']): e for e in edges}
    pair_results = []
    for a, b in pairs:
        evidence = sequence_evidence(all_records[a]['canonical_sequence'], all_records[b]['canonical_sequence'])
        kind = grouping_edge(all_records[a], all_records[b], evidence)
        original = edge_by.get((a, b), edge_by.get((b, a)))
        expected = original['method'] if original else None
        assert kind == expected, (a, b, kind, expected)
        if original:
            assert all(abs(evidence[k] - original[k]) < 1e-12 for k in evidence)
        pair_results.append({'a': a, 'b': b, 'method': kind, **evidence})
    report = {'status': 'PASS', 'raw_cases': raw_results, 'group_pairs': pair_results, 'scope': 'Bounded independent raw-column/SIFTS/frame check and selected alignment replay; no graph completeness, phylogeny, all-source validation or random representativeness claim'}
    (args.output / 'bounded_upstream.json').write_text(json.dumps(report, indent=2) + '\n')
    print('PASS raw cases', len(raw_results), 'residues', sum(r['residues'] for r in raw_results), 'group pairs', len(pair_results))


if __name__ == '__main__':
    main()
