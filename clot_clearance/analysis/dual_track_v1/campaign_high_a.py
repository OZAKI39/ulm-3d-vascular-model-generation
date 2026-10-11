from analysis.dual_track_v1.common import *
from analysis.dual_track_v1.high_helpers import risk
from analysis.dual_track_v1.high_a import run
import subprocess,time,traceback
assert json.loads((A/'PVCE_ENERGY_RELEASE_RESULTS.json').read_text())['status']=='PASS';assert json.loads((A/'STATIC_SOLVER_RECOVERY.json').read_text())['status']=='PASS'
c=json.loads((A/'PVCE_DESIGN_FREEZE.json').read_text());c['criteria']['finite_difference_force_scaled_absolute']=1e-6;c['protocol']='V1 scaling and L-BFGS-B, at most one verified Y-only exact Hessian-vector Newton correction per wide high candidate, max600CG iterations, unchanged acceptance criteria; all within same900s/8GiB formal equilibrium budget.';c['frozen_unix']=time.time();save(A/'HIGH_NUMERICAL_FREEZE.json',c)
p=ROOT/'analysis/dual_track_v1/hvp_y.cpp';r=subprocess.run(['g++','-O3','-shared','-fPIC',str(p),'-o',str(p.with_suffix('.so'))],capture_output=True,text=True);assert r.returncode==0
z=dict(np.load(A/'A0/equilibrium.npz'));r0=risk(z);amp=.95*.15/r0['max_stretch']*2e-6;records=[];start=time.monotonic();decision=''
try:
 for i in range(2):
  r=run('wide_A0_search_'+str(i+1),amp,warm=A/'A0/equilibrium.npz');records.append(r)
  if not r['accepted']:decision='STOP_STATIC_EQUILIBRIUM_FAILED';break
  if r['risk']['leading_region']!='notch_front_periodic_layer':decision='STOP_NON_NOTCH_COMPETITOR_LEADS';break
  if r['risk']['max_utilization']>=1:decision='THRESHOLD_BRACKET_CANDIDATE';break
  amp=1.05*amp/r['risk']['max_utilization']
 # A high PVCE requires an accepted subcritical notch leader with >2% margin.
 candidate=next((v for v in reversed(records) if v['accepted'] and v['risk']['max_utilization']<1 and v['risk']['leading_region']=='notch_front_periodic_layer'),None)
 if decision=='THRESHOLD_BRACKET_CANDIDATE' and candidate:
  competing=max(v['max_utilization'] or 0 for k,v in candidate['risk']['regions'].items() if k!='notch_front_periodic_layer')
  if competing<.98*candidate['risk']['max_utilization']:
   high=[candidate]
   for case in ['A1','A2']:
    rr=run('high_'+case,candidate['opening_m'],case,warm=A/case/'equilibrium.npz');high.append(rr)
    if not rr['accepted']:break
   if len(high)==3 and all(v['accepted'] for v in high):
    u=np.array([v['U_J'] for v in high]);save(A/'NEAR_ONSET_PVCE_RESULTS.json',dict(identity='SPECIMEN_G_DRIVE_NOT_GC',opening_m=candidate['opening_m'],U_J=u,G_two_J_m2=(u[0]-u[1])/2e-8,G_three_J_m2=(3*u[0]-4*u[1]+u[2])/4e-8))
 save(A/'FIRST_FAILURE_AUDIT.json',dict(decision=decision,low_risk=r0,high_A0_candidates=records,actual_natural_break_run=False,first_failure_identified=False,reason='Frozen-topology static branch diagnostic; no natural dynamic fracture. A competitor lead terminates high PVCE.',wall_s=time.monotonic()-start))
except BaseException as e:save(A/'HIGH_STOP.json',dict(error=repr(e),traceback=traceback.format_exc(),states_completed=len(records),wall_s=time.monotonic()-start));raise
