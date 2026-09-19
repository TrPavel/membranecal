"""Post-outcome description of the fixed matched estimand; no study edits."""
import json,os
from pathlib import Path
from collections import Counter
B=Path(os.environ['MC_RELEASE_ROOT']); O=Path(os.environ['MC_AUDIT_OUTPUT'])
m=json.loads((B/'data/matched_sets.json').read_text())['primary']
g=json.loads((B/'data/groups.json').read_text()); family=g['family_by_accession']; result={}
for contrast in ['MEMBRANE_INTERFACE-TM_CORE','STRUCTURED_EXTRAMEMBRANE-TM_CORE']:
    rows=[x for x in m if x['contrast']==contrast]
    counts=Counter(family[x['accession']] for x in rows)
    cells=Counter(); mass=Counter(); weighted=Counter(); sheet=[]
    for r in rows:
        for c in r['cells']:
            key=c['secondary_structure'];cells[key]+=1;mass[key]+=c['mass']
            w=c['mass']/r['mass']/counts[family[r['accession']]]/len(counts)
            weighted[key]+=w
            if key=='SHEET':sheet.append({'accession':r['accession'],'bin':c['bin'],'mass':c['mass'],'protein_mass':r['mass']})
            if c['bin']>=45:weighted['plddt_bin_ge90']+=w
            if c['bin']>=40:weighted['plddt_bin_ge80']+=w
    result[contrast]={'cells':dict(cells),'common_mass':dict(mass),'estimand_weight_fraction':dict(weighted),'sheet_cells':sheet}
(O/'population_diagnostics.json').write_text(json.dumps(result,indent=2)+'\n')
