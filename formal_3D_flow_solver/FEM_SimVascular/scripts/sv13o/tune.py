"""Bounded A/B/(conditional C)/D tuning, at most one ambiguous repeat."""
import json,re,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv13o import selection_action,confirm_selection
R=ROOT/'reports/sv1_3o';C=ROOT/'configs/sv1_3o';S=ROOT/'scripts/sv13o'
def read(n):return json.loads((R/(n+'.json')).read_text())
def write(n,d):(R/(n+'.json')).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def execute(name,mode,candidate):
 assert not (R/(name+'_acceptance.json')).exists(),'Do not repeat completed work'
 print('START '+name,flush=True)
 p=subprocess.run([sys.executable,'-B',S/'run_case.py',name,mode,candidate])
 d=read(name+'_acceptance');assert p.returncode==0 or d['status']=='FAIL'
 return d
while not (R/'PERF_A_acceptance.json').exists():time.sleep(3)
A=read('PERF_A_acceptance');assert A['status']=='PASS'
assert A['start_step']==60 and A['stop_step']==70
window=read('checkpoint_window');window.update(status='PASS',native_reload='PASS: original native binary restart loader; exactly steps61..70',native_reload_run='PERF_A',runtime_checkpoint_sha256=A['initial_checkpoint_sha256'],new_binary_family_sha256=A['solver_sha256']);write('checkpoint_window',window)
log=(ROOT/'logs/sv1_3o/remote/PERF_A.log').read_text()
blocks=[int(x) for x in re.findall(r'total subdomain blocks\s*=\s*(\d+)',log)]
sizes=[tuple(map(int,x)) for x in re.findall(r'rows=(\d+), cols=(\d+)',log)]
one_block=bool(blocks) and set(blocks)=={1} and bool(sizes) and set(sizes)=={(281452,281452)}
write('asm_topology',dict(status='PASS' if one_block else 'NOT_SINGLE_BLOCK',blocks_per_solve=blocks,matrix_sizes=sorted(set(sizes)),MPI_ranks=1,local_rows=281452 if one_block else None,global_rows=281452 if one_block else None,direct_ILU_eligible=one_block,source='PERF_A actual KSP/PC/Mat views; one rank local=global'))
winner=A;winner_s=A['wall_time_s'];trials={'A':A};decisions=[]
def consider(candidate):
 global winner,winner_s
 d=execute('PERF_'+candidate,'window',candidate);trials[candidate]=d
 decision=selection_action(winner_s,d);entry=dict(candidate=candidate,incumbent=winner['candidate'],incumbent_wall_s=winner_s,**decision)
 if decision['action']=='ACCEPT':winner=d;winner_s=d['wall_time_s']
 elif decision['action']=='CONFIRM_ONCE':
  confirm=execute('PERF_'+candidate+'_CONFIRM','window',candidate);result=confirm_selection(winner_s,d,confirm);entry['confirmation']=result
  if result['accepted']:winner=d;winner_s=result['arithmetic_mean_s']
 entry.update(winner_after=winner['candidate'],winner_selection_time_s=winner_s);decisions.append(entry)
 write('tuning_progress',dict(trials={k:{x:v.get(x) for x in ('status','wall_time_s','statistics')} for k,v in trials.items()},decisions=decisions,winner=winner['candidate']))
 print('DECISION '+json.dumps(entry),flush=True)
if one_block:consider('B')
else:write('PERF_B_acceptance',dict(status='NOT_RUN',reason='Actual ASM is not a single full-size block'))
B=read('PERF_B_acceptance')
c_eligible=one_block and (B['status']!='PASS' or 1-B['wall_time_s']/A['wall_time_s']<.10)
if c_eligible:consider('C')
else:write('PERF_C_acceptance',dict(status='NOT_RUN',reason='B already improves A by at least10%, or direct ILU topology prerequisite is absent',conditional_candidate=True))
options=winner['PETSC_OPTIONS'].replace('-ksp_gmres_restart 100','-ksp_gmres_restart 200')
assert options!=winner['PETSC_OPTIONS']
(C/'candidate_D.json').write_text(json.dumps(dict(base_candidate=winner['candidate'],PETSC_OPTIONS=options,only_change='GMRES restart100 ->200'),indent=2)+'\n')
consider('D')
result=dict(status='PASS',candidate=winner['candidate'],PETSC_OPTIONS=winner['PETSC_OPTIONS'],selection_wall_time_s=winner_s,baseline_A_wall_s=A['wall_time_s'],fixed_window_reduction=1-winner_s/A['wall_time_s'],decisions=decisions,window_start=60,window_end=70,checkpoint_sha256=A['initial_checkpoint_sha256'],config_sha256=A['config_sha256'],solver_sha256=A['solver_sha256'],PETSc_library_sha256=A['PETSc_library_sha256'],scientific_equivalence='DEFERRED',criterion='Same10-step window wall time and solver health; not iteration count; no median')
write('winner',result)
print('WINNER '+json.dumps({k:result[k] for k in ('candidate','selection_wall_time_s','fixed_window_reduction')}),flush=True)
# Profiling is limited to the two selected configurations and never used as timing evidence.
profileA=execute('PROFILE_A','profile','A');assert profileA['status']=='PASS'
if winner['candidate']=='A':profileWinner=profileA
else:profileWinner=execute('PROFILE_WINNER','profile',winner['candidate']);assert profileWinner['status']=='PASS'
write('profile_runs',dict(status='PASS',baseline=profileA['name'],winner=profileWinner['name'],runs=1 if winner['candidate']=='A' else 2,excluded_from_winner_selection=True,scope='Fixed10-step windows only; PETSc log_view with GPU event timing, no Nsight or transfer tracing'))
print('Tuning and permitted lightweight profiles finished; ready for the single full steady run.',flush=True)
