"""Freeze Stage P and create only isolated Stage Q work areas."""
import gzip,json,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import inventory,sha256,git_state,now
from sv_validation.sv13n import checkpoint_one_rank
R=ROOT/'reports/sv1_3q';C=ROOT/'configs/sv1_3q';P=ROOT/'reports/sv1_3p'
for kind in ('reports','configs','outputs','logs','patches','benchmarks'):(ROOT/kind/'sv1_3q').mkdir(parents=True,exist_ok=True)
def read(p):return json.loads(p.read_text())
def write(p,d):p.write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')
assert not (R/'reference_freeze.json').exists()
shutil.copyfile('/home/lzy/.codex/attachments/e5e90640-5d3d-4280-8e80-49e636e82e9a/pasted-text.txt',R/'USER_REQUEST.txt')
delivery=read(P/'delivery_manifest.json')
for f in delivery['files']:assert sha256(ROOT/f['path'])==f['sha256'],f['path']
files={};scopes=[]
for kind in ('reports','outputs','logs','configs','patches','benchmarks'):
 scopes += [str(p.relative_to(ROOT)) for p in (ROOT/kind).iterdir() if p.is_dir() and p.name!='sv1_3q']
scopes += ['inputs','external/sv13o','external/sv13p','external/petsc325/svMultiPhysics-compat','external/petsc325_cuda13_hypre','external/petsc325_cuda13_amgx']
for name in scopes:files.update({name+'/'+k:v for k,v in inventory(ROOT/name)['files'].items() if '__pycache__' not in k})
for name in ('configs','src','scripts','tests'):
 for k,v in inventory(ROOT/name)['files'].items():
  if 'sv13q' not in k and 'sv1_3q' not in k and '__pycache__' not in k:files[name+'/'+k]=v
(R/'history_baseline.json.gz').write_bytes(gzip.compress(json.dumps(dict(files=files,scopes=scopes,timestamp=now())).encode(),mtime=0))
shutil.copyfile(P/'old_fem_baseline.json.gz',R/'old_fem_baseline.json.gz')
policy=read(ROOT/'configs/sv1_3p/policy.json');ref=read(P/'reference_freeze.json');full=read(P/'remote/REAL_VASCULAR_GPU_PC_WINNER_execution.json');rows=full['history']['linear_solves']
r10=[r for r in rows if r['step']==10][-1];r20=[r for r in rows if r['step']==20][-1]
early=dict(source='Stage P unique full steady run, native reported elapsed time at final nonlinear rows of steps 10 and 20',start_row=r10,end_row=r20,wall_time_s=r20['reported_elapsed_s']-r10['reported_elapsed_s'],timing_resolution_s=1.,iterations=sum(r['linear_iterations'] for r in rows if 10<r['step']<=20),KSP_solves=sum(10<r['step']<=20 for r in rows),rebuild_count=10,standalone_restart_baseline=False,limitation='Historical contiguous interval excludes restart/startup cost, whereas new early windows include them. Single observed comparison, not formal benchmark; no R1 rerun.')
cp=ROOT/'outputs/sv1_3p/REAL_VASCULAR_GPU_PC_WINNER/1-procs/stFile_010.bin';cpcheck=checkpoint_one_rank(cp,10,policy['production_policy']['dt_s'])
write(R/'early_checkpoint.json',dict(status='PASS',source=str(cp.relative_to(ROOT)),remote_source='sv1_3p/outputs/REAL_VASCULAR_GPU_PC_WINNER/1-procs/stFile_010.bin',sha256=sha256(cp),native=cpcheck,binary_family='Stage P adapter on frozen PETSc 3.25.5/CUDA13.2',rank=1,source_run_acceptance=read(P/'REAL_VASCULAR_GPU_PC_WINNER_acceptance.json')['status']))
write(R/'baseline_r1.json',dict(status='PASS',late=read(P/'P1_WINDOW_acceptance.json'),early=early,full=read(P/'winner_steady_candidate.json'),no_baseline_rerun=True))
policy.update(candidate_order=['R2','R3','R5','RA'],adaptive_threshold=1.5,max_reuse_age=5,I_ref_rule='First healthy fresh-ILU solve; updated only after a healthy fresh solve',interval_anchor='First actual solve after startup/restart; no PC exists in native checkpoint',age_definition='current timestep minus last rebuild timestep; rebuild before solve when age >= interval or RA age >= 5',adaptive_recovery='At most one retry per stale linear solve; original RHS restored, same matrix and zero initial guess; retain both attempts',smoke_start=0,smoke_end=3,early_start=10,early_end=20,early_checkpoint_sha256=sha256(cp),early_checkpoint_remote='sv1_3p/outputs/REAL_VASCULAR_GPU_PC_WINNER/1-procs/stFile_010.bin',early_baseline=early,late_baseline=read(P/'P1_WINDOW_acceptance.json')['wall_time_s'],ranking='Top 2 healthy late windows -> early windows -> minimum early+late observed wall; early window must not exceed historical R1 by >5%',smoke_resource_budget_s=900.,maximum_full_steady_runs=1)
policy.pop('P1_options',None);policy.pop('early_slow_reject',None)
write(C/'policy.json',policy)
write(R/'reference_freeze.json',dict(ref,Stage_P_delivery_sha256=sha256(P/'delivery_manifest.json'),Stage_P_files_verified=len(delivery['files']),GPU_baseline=read(P/'winner_steady_candidate.json'),window_baseline=read(P/'P1_WINDOW_acceptance.json'),svmp=read(P/'remote/svmp_reuse_build.json'),official_source_git=git_state(ROOT/'external/svMultiPhysics')))
source=ROOT/'external/sv13p/svMultiPhysics-reuse';target=ROOT/'external/sv13q/svMultiPhysics-reuse';assert not target.exists();shutil.copytree(source,target)
for n in ('remote.py','runner_remote.py','environment_remote.py','invoke.py','fetch_case.py','sync_remote.py','reload_field.py','flow_parser.py'):
 text=(ROOT/'scripts/sv13p'/n).read_text().replace('sv1_3p','sv1_3q').replace('sv13p','sv13q')
 if n=='environment_remote.py':text=text.replace("'sv1_3o','sv1_3g", "'sv1_3o','sv1_3p','sv1_3g")
 (ROOT/'scripts/sv13q'/n).write_text(text)
write(R/'official_sources.json',dict(PETSc_version='3.25.5',urls=['https://petsc.org/release/manualpages/KSP/KSPSetReusePreconditioner/','https://petsc.org/release/manualpages/KSP/KSPSolve/','https://petsc.org/release/manualpages/KSP/KSPSetErrorIfNotConverged/'],source_audit='itfunc.c: in-place KSPSolve overwrites RHS; diagonal_scale_fix restores scaled matrix on return, including negative convergence reason; in-place solve requires zero initial guess.'))
print(f'Frozen Stage P: {len(delivery["files"])} files; early R1 observed interval = {early["wall_time_s"]} s; native step10 checkpoint valid.')
