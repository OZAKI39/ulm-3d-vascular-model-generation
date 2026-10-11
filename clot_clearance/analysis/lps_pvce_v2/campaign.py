from analysis.lps_pvce_v2.common import *
import subprocess,sys,traceback

def continuum():
 results={};start=time.perf_counter()
 for level in ['coarse','refined']:
  rows=[]
  for case in ['A0','A1','A2']:
   out=BASE/'continuum'/level/case;out.mkdir(parents=True,exist_ok=True)
   if (out/'RESULT.json').exists():
    # Resume only a completed stored result after stdout serialization failure; never re-solve it.
    returncode=0
   else:
    with (out/'stdout.log').open('w') as so,(out/'stderr.log').open('w') as se:p=subprocess.run([sys.executable,'-m','analysis.lps_pvce_v2.fem',level,case],stdout=so,stderr=se)
    returncode=p.returncode
   if returncode!=0 or not (out/'RESULT.json').exists():save(BASE/'CONTINUUM_REFERENCE_AUDIT.json',dict(status='FAILED_REFERENCE',failed=str(out),levels=results,wall_seconds=time.perf_counter()-start));return
   r=json.loads((out/'RESULT.json').read_text());rows.append(r);print('FEM',level,case,r['accepted'],r['U_J'],flush=True)
   if not r['accepted']:save(BASE/'CONTINUUM_REFERENCE_AUDIT.json',dict(status='FAILED_REFERENCE',levels={**results,level:rows},wall_seconds=time.perf_counter()-start));return
  U=np.array([r['U_J'] for r in rows]);results[level]=dict(states=rows,G_two_J_m2=float((U[0]-U[1])/3e-8),G_three_J_m2=float((3*U[0]-4*U[1]+U[2])/6e-8))
  save(BASE/'CONTINUUM_REFERENCE_PROGRESS.json',results)
 changes={k:results['refined'][k]/results['coarse'][k]-1 for k in ['G_two_J_m2','G_three_J_m2']};mesh_ok=max(abs(v) for v in changes.values())<=.05
 save(BASE/'CONTINUUM_REFERENCE_AUDIT.json',dict(status='VALID_COMPATIBLE_REFERENCE',mesh_scope='ONE_REFINEMENT_5PCT_SCREEN_PASSED' if mesh_ok else 'MESH_UNRESOLVED_WITHIN_BUDGET',levels=results,relative_G_change=changes,wall_seconds=time.perf_counter()-start,mesh_limit='Only one x-z tip-band refinement through full thickness; no asymptotic convergence or y-resolution certification.'))

def onset():
 from analysis.lps_pvce_v2.lps import run
 cfg=json.loads((BASE/'V2_DESIGN_FREEZE.json').read_text())['B'];base=json.loads((BASE/'BASELINE_STATE_AUDIT.json').read_text());old=json.loads((OLD/'A0/EQUILIBRIUM.json').read_text());samples=[dict(name='V1_A0_reused',opening_m=2e-6,accepted=True,risk=base,U_J=old['U_J'],last=old['last'],state_path=str(OLD/'A0/equilibrium.npz'))];lo=samples[0];hi=None;stop=None
 for i in range(6):
  if i==0:amp=cfg['search']['initial_m']
  elif hi is None:amp=samples[-1]['opening_m']*1.25
  else:
   low=lo['opening_m'];high=hi['opening_m'];root=low+(1-lo['risk']['max_utilization'])*(high-low)/(hi['risk']['max_utilization']-lo['risk']['max_utilization']);target=root*(.99 if i%2==0 else 1.01);amp=float(np.clip(target,low+.1*(high-low),high-.1*(high-low)))
  warm=min(samples,key=lambda s:abs(s['opening_m']-amp));r=run(f'A0_search_{i+1:02d}',amp,warm=warm['state_path'],warm_amp=warm['opening_m']);r['state_path']=str(BASE/'lps'/r['name']/'equilibrium.npz');samples.append(r)
  if not r['accepted']:stop='STOPPED_STATIC_OR_CONSISTENCY_GATE';break
  ordered=sorted(samples,key=lambda s:s['opening_m']);util=[s['risk']['max_utilization'] for s in ordered];reaction=[s['last']['top_reaction_N'] for s in ordered]
  if np.any(np.diff(util)<=0) or np.any(np.diff(reaction)<=0):stop='STOPPED_NONMONOTONE_SAMPLED_BRANCH';break
  if r['risk']['max_utilization']<1:lo=r
  else:hi=r
  width=(hi['opening_m']-lo['opening_m'])/lo['opening_m'] if hi else None
  save(BASE/'A0_SEARCH_PROGRESS.json',dict(samples=samples,lower=lo['opening_m'],upper=hi['opening_m'] if hi else None,relative_width=width))
  if hi and width<=.03:stop='RESOLVED_3PCT_BRACKET';break
 if stop is None:stop='BUDGET_REACHED_BRACKET_UNRESOLVED'
 location=hi['risk']['location'] if hi and stop=='RESOLVED_3PCT_BRACKET' else 'UNDETERMINED';lowloc=lo['risk']['location'];admit=stop=='RESOLVED_3PCT_BRACKET' and location==lowloc and location in ['NOTCH_INTERIOR_FIRST','NOTCH_SIDE_SURFACE_FIRST'];decision=dict(admitted=admit,reason='MEANINGFUL_NOTCH_FIRST_AND_STATIC_GATES' if admit else ('BOUNDARY_CONTROLLED_FIRST_FAILURE' if location=='GRIP_FIRST' else 'STATIC_GATE_NOT_MET' if stop=='STOPPED_STATIC_OR_CONSISTENCY_GATE' else 'ADMISSION_GATE_NOT_MET'),high_A1_A2='NOT_RUN',G_high_two_J_m2=None,G_high_three_J_m2=None)
 if admit:
  high=[lo]
  for case in ['A1','A2']:
   r=run(case+'_high',lo['opening_m'],case);high.append(r)
   if not r['accepted']:break
  decision['high_A1_A2']='COMPLETED' if len(high)==3 and all(s['accepted'] for s in high) else 'INCOMPLETE_STATIC_FAILURE'
  if decision['high_A1_A2']=='COMPLETED':
   U=np.array([s['U_J'] for s in high]);rr=dict(status='NEAR_ONSET_SPECIMEN_PVCE_G_DRIVE',opening_m=lo['opening_m'],states=high,G_high_two_J_m2=float((U[0]-U[1])/3e-8),G_high_three_J_m2=float((3*U[0]-4*U[1]+U[2])/6e-8),scope='Frozen-topology hypothetical branch if any virtual state exceeds native threshold; not materialGc');save(BASE/'NEAR_ONSET_PVCE_RESULTS.json',rr);decision.update({k:rr[k] for k in ['G_high_two_J_m2','G_high_three_J_m2']})
 search=dict(status=stop,new_A0_states=len(samples)-1,samples=samples,critical_opening_bracket_m=[lo['opening_m'],hi['opening_m']] if hi else None,relative_bracket_width=(hi['opening_m']-lo['opening_m'])/lo['opening_m'] if hi else None,lower_state=lo['name'],upper_state=hi['name'] if hi else None,first_location=location,admission=decision,critical_scope='First sampled monotone fixed-topology elastic-branch threshold interval, not dynamic rupture or a uniqueness proof; upper state is above criterion and hypothetical.')
 save(BASE/'A0_FIRST_FAILURE_SEARCH.json',search);save(BASE/'FIRST_FAILURE_LOCATION_AUDIT.json',dict(label=location,interpretation=decision['reason'],lower=lo['risk'],upper=hi['risk'] if hi else None,lower_reaction_N=lo['last']['top_reaction_N'],upper_reaction_N=hi['last']['top_reaction_N'] if hi else None,geometry_source='V1 own reference and m; V3 spatial definitions only, not its load/time values.'))
 print('ONSET FINISHED',stop,search['critical_opening_bracket_m'],location,decision,flush=True)
if __name__=='__main__':
 try:globals()[sys.argv[1]]()
 except BaseException:(BASE/(sys.argv[1]+'_ERROR.txt')).write_text(traceback.format_exc());raise
