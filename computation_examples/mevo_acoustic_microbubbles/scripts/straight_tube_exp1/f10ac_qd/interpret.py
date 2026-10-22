"""Conservative mode classification, seed robustness and frequency-scoped admission."""
import numpy as np
from .common import *

def main():
 val=read(BASE/'VALIDATION.json')
 if (BASE/'OUTSIDE_VALIDATION.json').exists():
  od=read(BASE/'OUTSIDE_VALIDATION.json')
  if od.get('status')=='COMPLETE':val['groups'].append(od['group']);val['comparisons']+=od['comparisons']
 write(BASE/'VALIDATION_COMBINED.json',val)
 rows=read(BASE/'ALL_RUNS.json');byid={r['id']:r for r in rows};modes=[]
 for g in val['groups']:
  ref=byid[g['replicates'][0]['id']];noise=np.max([r['noise_bound'] for r in g['replicates']],axis=0);d=np.array(ref['descriptor']) if ref['descriptor'] is not None else np.zeros(2)
  # Categories reflect measured signs, never claim every occupied cell is a new mode.
  signs=np.where(abs(d)>np.maximum(1e-4,5*noise),np.sign(d),0).astype(int)
  principal=np.array(ref['final']['principal_rms_lengths_3d_m'])/ref['initial']['principal_rms_lengths_3d_m']-1
  cc=[c for c in val['comparisons'] if c['reference']==ref['id'] and 'errors' in c];posnoise=max([c['errors']['fixed_id_m'] for c in cc],default=np.inf)
  extent=np.array(ref['final']['principal_rms_lengths_3d_m'])-ref['initial']['principal_rms_lengths_3d_m']
  absolute=bool(np.any(extent>max(1e-9,5*posnoise)))
  modes.append(dict(elite=g['elite'],source=g['source'],frequency=g['frequency'],descriptor=ref['descriptor'],category=signs.tolist(),classification=('Rg contraction' if signs[0]<0 else 'Rg expansion' if signs[0]>0 else 'Rg unresolved')+'; '+('aspect increase' if signs[1]>0 else 'aspect decrease' if signs[1]<0 else 'aspect unresolved'),repeated=g['signed_descriptor_repeated'],connected_all=g['connected_all'],principal_fraction=principal,principal_absolute_change_m=extent,absolute_extension_resolved=absolute,noise_bound=noise))
 confirmed=[m for m in modes if m['repeated'] and any(m['category'])];categories=sorted(set(tuple(m['category']) for m in confirmed));historical=[]
 for p in AB.glob('discovery/*/SUMMARY.json'):
  r=read(p)
  if r['state_id']=='S2_r0' and r['status']=='PASS':
   a,b=r['initial_shape'],r['final_shape'];d=np.log([b['Rg_3d_m']/a['Rg_3d_m'],b['aspect_3d']/a['aspect_3d']]);historical.append(dict(id=r['run_id'],descriptor=d,category=np.where(abs(d)>1e-4,np.sign(d),0).astype(int)))
 historical_categories=set(tuple(h['category']) for h in historical);novel=[c for c in categories if c not in historical_categories]
 # Global response diversity is not a local control-Jacobian rank.
 matrix=np.array([m['descriptor'] for m in confirmed]);sv=np.linalg.svd(matrix,compute_uv=False) if len(matrix) else np.array([]);floor=max([np.linalg.norm(m['noise_bound'])*5 for m in confirmed],default=np.inf);effective=int(np.count_nonzero(sv>max(1e-4,floor)))
 comparisons=read(AD/'LONG_CONVERGENCE.json');long=read(AD/'LONG.json');frequencies={}
 for f in [2400000.,2500000.,2600000.]:
  rr=[r for r in long if r['task']['frequency']==f];items=[]
  for r in rr:
   cc=[c for c in comparisons if c['reference']==r['id']];items.append(dict(id=r['id'],action=r['id'].split('_')[-1],complete=r['status']=='PASS',original_numerical_pass=len(cc)==2 and all(c['admitted'] for c in cc),failure=r.get('reasons',[]),scope='S2_r0; N64; 2ms; amplitudes 700 Pa-eq; listed action only'))
  frequencies[str(int(f))]=dict(physics=read(AD/str(int(f))/'PHYSICS_GATE.json')['PASS'],dynamic_actions=items,scoped_dynamic_admitted=any(i['original_numerical_pass'] for i in items),all_tested_long_actions_pass=bool(items) and all(i['original_numerical_pass'] for i in items))
 robust=read(AD/'ROBUSTNESS.json');shell_rows=[]
 for shell in [.8,1.,1.2]:
  rs=[r for r in (long if shell==1 else robust) if r['status']=='PASS' and r['task'].get('shell',1)==shell and ('X180' in r['id']) and r['task'].get('variant','DT1')=='DT1' and not r['task'].get('offset_m')]
  rs=sorted(rs,key=lambda r:r['task']['frequency'])
  shell_rows.append(dict(shell=shell,frequency_Hz=[r['task']['frequency'] for r in rs],delta_log_aspect=[r['descriptor'][1] for r in rs],ranking_high_to_low=[r['task']['frequency'] for r in sorted(rs,key=lambda r:r['descriptor'][1],reverse=True)]))
 complete_shell=all(len(r['frequency_Hz'])==3 for r in shell_rows)
 rankings=[r['ranking_high_to_low'] for r in shell_rows]
 radial=[]
 for shell in [.8,1.,1.2]:
  data=[read(AD/str(f)/f'RESPONSE_{shell:.1f}.json') for f in [2400000,2500000,2600000]]
  yy=[next(n for n in d['nodes'] if n['pressure_Pa']==1500)['metrics']['max_R_over_R0'] for d in data]
  radial.append(dict(shell=shell,max_R_over_R0_at_1500Pa=yy,largest_tested_response_frequency_Hz=[2400000,2500000,2600000][int(np.argmax(yy))],global_resonance_identified=False))
 write(AD/'SHELL_SENSITIVITY_SUMMARY.json',dict(trajectories=shell_rows,complete=complete_shell,frequency_shape_order_changes=complete_shell and any(r!=rankings[0] for r in rankings),radial=radial,scope='Three sampled frequencies do not locate the global resonance'))
 archives=[]
 for path in list((BASE/'search_v2').glob('*/RESULT.json'))+list((BASE/'extensions').glob('*/RESULT.json')):
  record=read(path);f=record.get('frequency',2500000.)
  for category in ['archive','connected_archive']:
   if category not in record:continue
   ar=record[category]['elites']
   for i,solution in enumerate(ar['solution']):
    matches=[]
    for r in rows:
     t=r['task']
     if r['status']!='PASS' or t['state']!='S2_r0' or t.get('frequency',2500000.)!=f or t.get('variant','DT1')!='DT1':continue
     a=t['amplitudes'];p=t['phases'];x=np.array([(a[2]-200)/900,(np.degrees(p[3])-150)/90,a[0]/1100,p[1]/(2*np.pi)])
     if np.allclose(x,solution,rtol=0,atol=1e-12):matches.append(r)
    if not matches:raise RuntimeError('Archive elite has no actual trajectory: '+str(path))
    own_ids={r['id'] for r in record['rows']};r=next((r for r in matches if r['id'] in own_ids),matches[0])
    if not np.allclose(r['descriptor'],ar['measures'][i],rtol=0,atol=1e-10):raise RuntimeError('Elite measure does not match actual stored trajectory')
    archives.append(dict(archive=str(path.relative_to(BASE)),task_type=category,cell=ar['index'][i],solution=solution,objective=ar['objective'][i],measures=ar['measures'][i],run_id=r['id'],trajectory_path=r['trajectory_path'],trajectory_sha256=r['trajectory_sha256'],particle_ids='trajectory.h5:/particle_id',forces='trajectory.h5:/force_samples',contact_reactions='trajectory.h5:/force_samples/pair_contact_reaction_N',final_shape=r['final'],original_solver_quality=r['max_metrics'],scope='Original full-step gates passed; dt/cache/seed admission only for explicitly validated representatives'))
 write(BASE/'ELITE_MANIFEST.json',archives)
 write(BASE/'INTERPRETATION.json',dict(modes=modes,confirmed_sign_categories=categories,distinct_sign_category_count=len(categories),new_vs_historical_sign_categories=novel,historical=historical,descriptor_response_singular_values=sv,noise_threshold=max(1e-4,floor),resolved_descriptor_span_rank=effective,rank_scope='Finite responses across different actions/frequencies, not a local controllability Jacobian and not proof of arbitrary shape control',frequency_admission=frequencies,identity=identity()))
if __name__=='__main__':main()
