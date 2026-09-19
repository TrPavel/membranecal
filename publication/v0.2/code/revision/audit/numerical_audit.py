import json, math
from pathlib import Path
from collections import defaultdict,Counter
import numpy as np
import os
ROOT=Path(os.environ['MC_RELEASE_ROOT']); OUT=Path(os.environ['MC_AUDIT_OUTPUT'])
def read(n):return json.loads((ROOT/n).read_text())
def conv(o):
 if isinstance(o,np.generic):return o.item()
 if isinstance(o,np.ndarray):return o.tolist()
 raise TypeError(type(o))
r=dict(np.load(OUT/'residues_preoutcome.npz',allow_pickle=True));s=dict(np.load(OUT/'scores.npz',allow_pickle=True))
assert np.array_equal(r['accession'],s['accession']) and np.array_equal(r['canonical_position'],s['canonical_position'])
n=len(s['accession']);bounds=np.r_[0,np.flatnonzero(r['accession'][1:]!=r['accession'][:-1])+1,n]
recomp=np.empty(n);strict=np.empty(n);strict32=np.empty(n);neighbor_counts=np.empty(n,dtype=int)
perprot=[];eq_radius=0;eq_threshold=Counter()
for lo,hi in zip(bounds[:-1],bounds[1:]):
 a=s['accession'][lo];ref=r['reference_ca'][lo:hi];pred=s['prediction_ca'][lo:hi]
 rd=np.sqrt(np.sum((ref[:,None]-ref[None,:])**2,axis=2));pd=np.sqrt(np.sum((pred[:,None]-pred[None,:])**2,axis=2));delta=np.abs(rd-pd)
 mask=(rd<=15)&~np.eye(hi-lo,dtype=bool);count=mask.sum(1)
 assert (count>0).all()
 recomp[lo:hi]=sum(((delta<=t)&mask).sum(1) for t in [.5,1,2,4])/(4*count)
 neighbor_counts[lo:hi]=count
 eq_radius+=int(np.sum(rd==15))
 for t in [.5,1,2,4]:eq_threshold[t]+=int(np.sum((delta==t)&mask))
 # Literal NumPy port of official AF2 lddt.py operations (float64 reference).
 rt=np.sqrt(1e-10+np.sum((ref[:,None]-ref[None,:])**2,axis=-1))
 pt=np.sqrt(1e-10+np.sum((pred[:,None]-pred[None,:])**2,axis=-1))
 ms=(rt<15)*(1-np.eye(hi-lo));ds=np.abs(rt-pt)
 hits=.25*sum((ds<t).astype(np.float32) for t in [.5,1,2,4])
 strict[lo:hi]=(1e-10+(ms*hits).sum(-1))/(1e-10+ms.sum(-1))
 ref32=ref.astype(np.float32);pred32=pred.astype(np.float32)
 rt=np.sqrt(1e-10+np.sum((ref32[:,None]-ref32[None,:])**2,axis=-1));pt=np.sqrt(1e-10+np.sum((pred32[:,None]-pred32[None,:])**2,axis=-1))
 ms=(rt<15).astype(np.float32)*(1-np.eye(hi-lo,dtype=np.float32));ds=np.abs(rt-pt)
 hits=.25*sum((ds<t).astype(np.float32) for t in [.5,1,2,4])
 strict32[lo:hi]=(1e-10+(ms*hits).sum(-1))/(1e-10+ms.sum(-1))
 perprot.append({'accession':a,'n':hi-lo,'max_score_diff':float(np.max(np.abs(recomp[lo:hi]-s['local_accuracy'][lo:hi]))),'strict64_max_pp':float(100*np.max(np.abs(strict[lo:hi]-recomp[lo:hi]))),'strict32_max_pp':float(100*np.max(np.abs(strict32[lo:hi]-recomp[lo:hi])))})
np.savez_compressed(OUT/'audit_rescores.npz',inclusive=recomp,af2_port64=strict,af2_port32=strict32)
proteins=read('data/proteins.json')['records'];pby={p['accession']:p for p in proteins};g=read('data/groups.json');m=read('data/matched_sets.json');stats=read('results/results.json');idx={(a,int(p)):i for i,(a,p) in enumerate(zip(r['accession'],r['canonical_position']))}
family=g['family_by_accession'];blocks=g['block_by_family'];empirical={};errors=recomp-r['confidence']
for contrast in ['MEMBRANE_INTERFACE-TM_CORE','STRUCTURED_EXTRAMEMBRANE-TM_CORE']:
 rows=[a for a in m['primary'] if a['contrast']==contrast];vals=defaultdict(dict)
 for mt in rows:
  a=mt['accession'];total=mt['mass']
  for field,data in [('error',errors),('accuracy',recomp),('confidence',r['confidence']),('af2_port64_error',strict-r['confidence']),('af2_port32_error',strict32-r['confidence']),('neighbors',neighbor_counts)]:
   terms=[]
   for cell in mt['cells']:
    li=np.array([idx[a,p] for p in cell['left_positions']]);ri=np.array([idx[a,p] for p in cell['right_positions']]);terms.append(cell['mass']*(float(data[li].mean())-float(data[ri].mean())))
   vals[field][a]=math.fsum(terms)/total
 gm={}
 for field,v in vals.items():
  by=defaultdict(list)
  for a,val in v.items():by[family[a]].append(val)
  gm[field]={f:math.fsum(v)/len(v) for f,v in by.items()}
 summaries={field:100*math.fsum(by.values())/len(by) for field,by in gm.items() if field!='neighbors'}
 byblock=defaultdict(list)
 for f,val in gm['error'].items():byblock[blocks[f]].append(val)
 sizes=np.array([len(v) for v in byblock.values()]);estimate=math.fsum(gm['error'].values())/len(gm['error'])
 shifts=[]
 for block,v in byblock.items():
  z=(math.fsum(gm['error'].values())-math.fsum(v))/(len(gm['error'])-len(v))
  shifts.append({'block':block,'groups':len(v),'weight':len(v)/len(gm['error']),'omit_effect_pp':z*100,'shift_pp':100*(z-estimate)})
 # Independently construct block sums/counts and same predeclared RNG; no frozen analysis imports.
 ids=sorted(byblock);sumv=np.array([math.fsum(byblock[b]) for b in ids]);cnt=np.array([len(byblock[b]) for b in ids]);draw=np.random.default_rng(20909).integers(len(ids),size=(20000,len(ids)))
 draws=sumv[draw].sum(1)/cnt[draw].sum(1);ci=np.quantile(draws,[.0125,.9875])*100
 original=stats['contrasts']['primary'][contrast]
 # Weights/covariate means across actual retained cells, same estimand.
 pooledbin=Counter();leftidx=[];rightidx=[]
 for mt in rows:
  for c in mt['cells']:pooledbin[c['bin']]+=c['mass'];leftidx.extend(idx[mt['accession'],p] for p in c['left_positions']);rightidx.extend(idx[mt['accession'],p] for p in c['right_positions'])
 empirical[contrast]={'effects_pp':summaries,'interval_pp':ci.tolist(),'original_estimate_diff':estimate-original['estimate'],'original_ci_diff':(ci/100-np.array(original['interval'])).tolist(),'proteins':len(rows),'groups':len(gm['error']),'blocks':len(byblock),'block_size_distribution':dict(Counter(sizes.tolist())),'kish_effective_blocks':float(sizes.sum()**2/(sizes@sizes)),'largest_weight_blocks':sorted(shifts,key=lambda z:-z['weight'])[:10],'largest_omit_block_shifts':sorted(shifts,key=lambda z:-abs(z['shift_pp']))[:10],'common_mass_by_2point_bin':dict(sorted(pooledbin.items())),'paired_ss_cells':dict(Counter(c['secondary_structure'] for mt in rows for c in mt['cells'])),'mean_neighbor_contrast':math.fsum(gm['neighbors'].values())/len(gm['neighbors'])}
summary={'n':n,'proteins':len(perprot),'inclusive_score_exact_count':int((recomp==s['local_accuracy']).sum()),'neighbor_exact_count':int((neighbor_counts==s['neighbors']).sum()),'strict64_absdiff_max_pp':float(100*np.max(np.abs(strict-recomp))),'strict64_absdiff_gt_1e_9_count':int((np.abs(strict-recomp)>1e-9).sum()),'strict32_absdiff_max_pp':float(100*np.max(np.abs(strict32-recomp))),'strict32_gt_1e_4_count':int((np.abs(strict32-recomp)>1e-4).sum()),'exact_radius_hits_ordered':eq_radius,'exact_threshold_hits_ordered':dict(eq_threshold),'primary_diagnostics':empirical,'metadata':{'afdb_model_dates':dict(Counter(p['afdb_model_date'] for p in proteins)),'ss_counts':dict(Counter(r['secondary_structure'])),'region_counts':dict(Counter(r['region'])),'mean_confidence':float(r['confidence'].mean()),'mean_accuracy':float(recomp.mean())},'top_strict32_delta_proteins':sorted(perprot,key=lambda p:-p['strict32_max_pp'])[:15]}
(OUT/'numerical_audit.json').write_text(json.dumps(summary,indent=2,default=conv)+'\n')
print(json.dumps(summary,indent=2,default=conv)[:28000])
