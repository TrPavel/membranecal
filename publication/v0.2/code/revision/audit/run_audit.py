"""Independent audit diagnostics. Read-only with respect to the original RC1."""
import argparse,hashlib,json,os,runpy,sys
from pathlib import Path
import numpy as np

def main():
    p=argparse.ArgumentParser();p.add_argument('--release-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if sys.flags.optimize: raise RuntimeError('Assertions must remain enabled.')
    root=a.release_root.resolve(); out=a.output.resolve()
    if root==out or root in out.parents or out in root.parents: raise ValueError('Use a separate, non-nested audit output directory.')
    if out.exists() and any(out.iterdir()): raise ValueError('Audit output must be fresh or empty.')
    manifest=json.loads((root/'manifest.json').read_text())['files']
    for n,h in manifest.items():
        f=(root/n).resolve()
        if not f.is_relative_to(root) or hashlib.sha256(f.read_bytes()).hexdigest()!=h:raise ValueError(f'Input integrity failed: {n}')
    try: import pyarrow.parquet as pq
    except ImportError as e:raise RuntimeError('Install pyarrow in this audit environment first.') from e
    out.mkdir(parents=True,exist_ok=True)
    for stem in ['residues_preoutcome','scores']:
        tab=pq.read_table(root/'data'/f'{stem}.parquet')
        cols={k:np.asarray(tab[k].to_pylist()) for k in tab.column_names}
        np.savez_compressed(out/f'{stem}.npz',**cols)
    os.environ['MC_RELEASE_ROOT']=str(root);os.environ['MC_AUDIT_OUTPUT']=str(out)
    here=Path(__file__).resolve().parent
    runpy.run_path(str(here/'numerical_audit.py'),run_name='__main__')
    runpy.run_path(str(here/'population_diagnostics.py'),run_name='__main__')
    print(f'Audit diagnostics written to {out}; original release unchanged.')
if __name__=='__main__':main()
