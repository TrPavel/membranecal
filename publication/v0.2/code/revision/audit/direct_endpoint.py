"""Run unchanged upstream lddt on every compact case; post-outcome validation only."""
import argparse
import hashlib
import json
import math
import os
from collections import defaultdict
from pathlib import Path

os.environ['JAX_PLATFORMS'] = 'cpu'
import jax
import jax.numpy as jnp
import numpy as np
from official_lddt import lddt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--release-root', type=Path, required=True)
    parser.add_argument('--audit-output', type=Path, required=True)
    parser.add_argument('--dtype', choices=['float64', 'float32'], required=True)
    args = parser.parse_args()
    jax.config.update('jax_enable_x64', args.dtype == 'float64')
    dtype = getattr(np, args.dtype)
    output = args.audit_output
    r = dict(np.load(output / 'residues_preoutcome.npz', allow_pickle=True))
    s = dict(np.load(output / 'scores.npz', allow_pickle=True))
    n = len(r['accession'])
    values = np.empty(n)
    score = jax.jit(lambda p, t, m: lddt(p, t, m, per_residue=True))
    accessions = sorted(set(r['accession']))
    for k, accession in enumerate(accessions):
        ix = np.flatnonzero(r['accession'] == accession)
        # Zero-mask padding reduces compilation shapes. All real observed points remain;
        # dummy points have zero pair weight and their per-residue outputs are discarded.
        size = math.ceil(len(ix) / 128) * 128
        ref = np.zeros((1, size, 3), dtype=dtype)
        pred = np.zeros_like(ref)
        mask = np.zeros((1, size, 1), dtype=dtype)
        ref[0, :len(ix)] = r['reference_ca'][ix]
        pred[0, :len(ix)] = s['prediction_ca'][ix]
        mask[0, :len(ix), 0] = 1
        values[ix] = np.asarray(score(jnp.asarray(pred), jnp.asarray(ref), jnp.asarray(mask)))[0, :len(ix)]
        if k % 100 == 0:
            print(args.dtype, k, '/', len(accessions), flush=True)
    groups = json.loads((args.release_root / 'data/groups.json').read_text())
    matches = json.loads((args.release_root / 'data/matched_sets.json').read_text())['primary']
    reference = json.loads((args.release_root / 'results/results.json').read_text())['contrasts']['primary']
    index = {(a, int(p)): i for i, (a, p) in enumerate(zip(r['accession'], r['canonical_position']))}
    effects = {}
    for contrast in ['MEMBRANE_INTERFACE-TM_CORE', 'STRUCTURED_EXTRAMEMBRANE-TM_CORE']:
        by = defaultdict(list)
        for row in matches:
            if row['contrast'] != contrast:
                continue
            acc = row['accession']
            terms = []
            for cell in row['cells']:
                left = [index[acc, p] for p in cell['left_positions']]
                right = [index[acc, p] for p in cell['right_positions']]
                v = values - r['confidence']
                terms.append(cell['mass'] * (float(v[left].mean()) - float(v[right].mean())))
            by[groups['family_by_accession'][acc]].append(math.fsum(terms) / row['mass'])
        estimate = math.fsum(math.fsum(v) / len(v) for v in by.values()) / len(by)
        effects[contrast] = {'estimate_pp': 100 * estimate, 'difference_from_frozen_pp': 100 * (estimate - reference[contrast]['estimate'])}
    difference = np.abs(values - s['local_accuracy'])
    report = {
        'source_sha256': hashlib.sha256(Path(__file__).with_name('official_lddt.py').read_bytes()).hexdigest(),
        'source_revision': 'dbaafbcdea0cf39fabe502928cab55754c1d5dc7',
        'source_modified': False, 'jax': jax.__version__, 'dtype': args.dtype,
        'platform': 'Windows CPU', 'execution': 'JAX jit of unchanged function, masked zero padding to multiples of 128',
        'scope': 'All paired-observed intrachain compact coordinates; not full canonical or all-atom truth',
        'proteins': len(accessions), 'residues': n,
        'max_residue_difference_pp': float(100 * difference.max()),
        'residues_difference_gt_1e_4_unit': int((difference > 1e-4).sum()),
        'primary_effects': effects,
    }
    assert all(abs(v['difference_from_frozen_pp']) < 0.001 for v in effects.values()), report
    np.savez_compressed(output / f'direct_endpoint_{args.dtype}.npz', values=values)
    (output / f'direct_endpoint_{args.dtype}.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
