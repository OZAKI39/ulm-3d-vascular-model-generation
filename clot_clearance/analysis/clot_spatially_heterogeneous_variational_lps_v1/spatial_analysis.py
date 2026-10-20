from common import *
import gzip,csv
from scipy.stats import spearmanr

def analyze_spatial(ref,results,topology):
 field=load_npz(O/'HETEROGENEITY_FIELD.npz');b=ref['bonds'];nb=len(b);bulk=(ref['bond_types']==1)&~ref['initial_cut'];mid=ref['X'][b].mean(1);allrows=[];out={}
 for mode in CASES:
  case=O/'cases'/mode;soft=np.full(nb,np.nan);fail=soft.copy();final=load_npz(case/'keys/final.npz');d=final['d'];sc=field['bond_sc'] if mode=='heterogeneous' else SC[ref['bond_types']]
  with gzip.open(case/'damage_events.jsonl.gz','rt') as f:
   for line in f:
    row=json.loads(line);ids=np.array(row['original_bond_ids'])-1;values=np.array(row['d']);new=ids[np.isnan(soft[ids])&(values>0)];soft[new]=row['time_s'];new=ids[np.isnan(fail[ids])&(values==1)];fail[new]=row['time_s']
  np.savez_compressed(case/'bond_event_times.npz',first_softening_s=soft,first_failure_s=fail,original_bond_ids=np.arange(1,nb+1),final_damage=d)
  for i in np.flatnonzero(~ref['initial_cut']):
   allrows.append(dict(case=mode,original_bond_id=i+1,atom_i=int(b[i,0]+1),atom_j=int(b[i,1]+1),type=int(ref['bond_types'][i]),reference_x_m=mid[i,0],reference_y_m=mid[i,1],reference_z_m=mid[i,2],reference_q_ij=field['bond_q'][i],actual_sc=sc[i],first_softening_s=None if np.isnan(soft[i]) else soft[i],first_failure_s=None if np.isnan(fail[i]) else fail[i],final_damage=d[i]))
  q=field['bond_q'];quartiles=np.quantile(q[bulk],[0,.25,.5,.75,1]);groups=[]
  for j in range(4):
   sel=bulk&(q>=quartiles[j])&(q<=quartiles[j+1] if j==3 else q<quartiles[j+1]);groups.append(dict(q_quartile=j+1,q_range=quartiles[j:j+2],original_eligible_bonds=int(sel.sum()),softened=int(np.isfinite(soft[sel]).sum()),failed=int(np.isfinite(fail[sel]).sum()),failed_fraction=float(np.mean(d[sel]==1)),mean_final_damage=float(np.mean(d[sel])),mean_sc=float(np.mean(sc[sel]))))
  hit=bulk&np.isfinite(soft);coef=spearmanr(q[hit],soft[hit]).statistic if hit.sum()>2 else None
  fd=bulk&(d==1);positions=mid[fd];centroid=positions.mean(0) if len(positions) else None;spread=np.sqrt(np.mean(np.sum((positions-centroid)**2,axis=1))) if len(positions) else None
  # Fixed reference octants: spatial coverage, not a crack reconstruction.
  center=np.mean(ref['X'][:int(ref['clot_count'])],axis=0);bins=(mid>center).astype(int)@np.array([1,2,4]);counts=np.bincount(bins[fd],minlength=8);concentration=float(np.sum((counts/counts.sum())**2)) if counts.sum() else None
  reg=json.loads((case/'COMPONENT_REGISTRY.json').read_text());composition=[]
  for r in topology[mode]['fragment_lifetimes']:
   ids=np.array(r['atom_ids'])-1;composition.append(dict(component_id=r['component_id'],nodes=len(ids),mean_q=float(field['q'][ids].mean()),q_range=[float(field['q'][ids].min()),float(field['q'][ids].max())],birth_s=r['birth_s'],lifetime_s=r['lifetime_s'],right_censored=r['right_censored']))
  secondary=results[mode]['native']['secondary_events'];cohesive=[e for e in secondary if reg[str(e['parent_component_id'])]['nodes']>=20]
  out[mode]=dict(q_quartile_damage=groups,q_vs_first_softening_spearman=float(coef) if coef is not None and np.isfinite(coef) else None,q_vs_final_damage_spearman=float(spearmanr(q[bulk],d[bulk]).statistic),failed_bond_midpoint_centroid_m=centroid,failed_bond_midpoint_rms_spread_m=spread,failed_bonds_per_reference_octant=counts,octant_Herfindahl=concentration,fragment_property_composition=composition,secondary_all_count=len(secondary),secondary_parent_at_least20_count=len(cohesive),secondary_parent_at_least20_events=cohesive,interpretation='Descriptive single-realization associations; geometry, original notch, load and nonlocal coupling confound spatial correlations. Uniform case uses q only as a fixed reference covariate, not a material property. No p-value or independent sample claim.')
 U=load_npz(O/'cases/uniform/keys/final.npz')['d'];H=load_npz(O/'cases/heterogeneous/keys/final.npz')['d'];fu=(U==1)&bulk;fh=(H==1)&bulk
 out['difference']=dict(final_bulk_failed_jaccard=float(np.sum(fu&fh)/np.sum(fu|fh)) if np.any(fu|fh) else None,uniform_only_failed=int(np.sum(fu&~fh)),heterogeneous_only_failed=int(np.sum(fh&~fu)),mean_absolute_bulk_damage_difference=float(np.mean(abs(U[bulk]-H[bulk]))),first_softening_same_original_bonds=results['uniform']['native']['events'].get('first_softening',{}).get('original_bond_ids')==results['heterogeneous']['native']['events'].get('first_softening',{}).get('original_bond_ids'),first_failure_same_original_bonds=results['uniform']['native']['events'].get('first_failure',{}).get('original_bond_ids')==results['heterogeneous']['native']['events'].get('first_failure',{}).get('original_bond_ids'))
 with (O/'NATIVE_EVENT_HISTORY.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=allrows[0]);w.writeheader();w.writerows(allrows)
 dump(O/'SPATIAL_DAMAGE_ANALYSIS.json',out);return out
