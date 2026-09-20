"""Freeze Stage N and the production inputs before any Stage O run."""
import difflib,gzip,json,shutil,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import inventory,sha256,git_state,now
from sv_validation.sv13n import checkpoint_one_rank
R=ROOT/'reports/sv1_3o';C=ROOT/'configs/sv1_3o';N=ROOT/'reports/sv1_3n'
assert not (R/'reference_freeze.json').exists()
for kind in ('reports','configs','outputs','logs','patches','benchmarks'):(ROOT/kind/'sv1_3o').mkdir(parents=True,exist_ok=True)
def read(p):return json.loads(p.read_text())
def write(p,d):p.write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')
shutil.copyfile('/home/lzy/.codex/attachments/459b3735-555b-4882-80ad-a93d71cf4ac2/pasted-text.txt',R/'USER_REQUEST.txt')
delivery=read(N/'delivery_manifest.json')
assert delivery['status']=='PASS'
for f in delivery['files']:assert sha256(ROOT/f['path'])==f['sha256'],f['path']
candidate=read(N/'gpu_steady_candidate.json');assert candidate['status']=='PASS'
production=read(ROOT/'configs/sv1_3/policy.json')
files={};scopes=[]
for kind in ('reports','outputs','logs','configs','patches','benchmarks'):
 for p in (ROOT/kind).iterdir():
  if p.is_dir() and p.name!='sv1_3o':scopes.append(str(p.relative_to(ROOT)))
scopes+=['inputs','external/petsc325/svMultiPhysics-compat']
for name in scopes:files.update({name+'/'+k:v for k,v in inventory(ROOT/name)['files'].items() if '__pycache__' not in k})
for name in ('configs','src','scripts','tests'):
 for k,v in inventory(ROOT/name)['files'].items():
  if 'sv13o' not in k and 'sv1_3o' not in k and '__pycache__' not in k:files[name+'/'+k]=v
(R/'history_baseline.json.gz').write_bytes(gzip.compress(json.dumps(dict(files=files,scopes=scopes,timestamp=now())).encode(),mtime=0))
shutil.copyfile(N/'old_fem_baseline.json.gz',R/'old_fem_baseline.json.gz')
source=ROOT/'external/petsc325/svMultiPhysics-compat';target=ROOT/'external/sv13o/svMultiPhysics'
assert not target.exists();shutil.copytree(source,target)
names=subprocess.check_output(['git','-C',str(source),'ls-files','-z']).decode().split('\0')
before={n:sha256(source/n) for n in names if n and (source/n).is_file()}
name='Code/Source/solver/main.cpp';old=(source/name).read_text()
needle='    if (save_vtu) {\n      vtk_xml::write_vtus(simulation, solutions, /* lAvg = */ false);\n    }'
assert old.count(needle)==1
replacement='''    // Preserve the normal VTU cadence, but save the completed final state
    // whenever the native stop condition is reached. Restart is written above.
    if (save_vtu || (com_mod.saveVTK && reached_stop_time_step)) {
      vtk_xml::write_vtus(simulation, solutions, /* lAvg = */ false);
    }'''
new=old.replace(needle,replacement);(target/name).write_text(new)
patch=''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='a/'+name,tofile='b/'+name))
patchpath=ROOT/'patches/sv1_3o/final_output_on_stop.patch';patchpath.write_text(patch)
after={n:sha256(target/n) for n in before};assert [n for n in before if before[n]!=after[n]]==[name]
write(R/'source_patch.json',dict(status='PASS',base='Stage N adopted lifecycle adapter',base_commit=git_state(source)['head'],source=str(target.relative_to(ROOT)),changed_files=[name],patch=str(patchpath.relative_to(ROOT)),patch_sha256=sha256(patchpath),before=before,after=after,scope='Output condition only after converged timestep; no assembly, solve, time integration, restart layout or lifecycle modification'))
cp=ROOT/'outputs/sv1_3n/REAL_VASCULAR_GPU/1-procs/stFile_060.bin'
checkpoint=checkpoint_one_rank(cp,60,production['dt_s'])
saved=next(s for s in read(N/'gpu_steady_history.json')['states'] if s['step']==60)
assert checkpoint['sha256']==saved['checkpoint']['sha256']
write(R/'checkpoint_window.json',dict(status='STRUCTURALLY_VERIFIED',start_step=60,end_step=70,timesteps=10,checkpoint=checkpoint,source_binary_sha256=candidate['solver_sha256'],PETSc_library_sha256=candidate['PETSc_library_sha256'],rank_count=1,backend='seqcuda',restart_layout_changed=False,native_reload='PENDING: first PERF_A must load this checkpoint and execute exactly steps61..70'))
baseline=dict(classification='DEVELOPMENT_BASELINE_OBSERVED',wall_time_s=candidate['wall_time_s'],steps=candidate['steps'],KSP_solves=candidate['linear_solves'],total_iterations=candidate['total_KSP_iterations'],mean_iterations=candidate['total_KSP_iterations']/candidate['linear_solves'],first_full_steady_step=candidate['first_full_steady_step'],stop_step=candidate['stop_step'],not_a_repeated_benchmark=True)
write(R/'baseline_runtime.json',baseline)
common='-skip_petscrc -ksp_type gmres -ksp_pc_side right -ksp_norm_type unpreconditioned -ksp_rtol 1e-10 -ksp_atol 1e-24 -ksp_max_it 2000 -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_monitor_true_residual -ksp_converged_reason -ksp_view -options_view -options_left -use_gpu_aware_mpi 0 -mat_type aijcusparse -vec_type cuda -vec_view ::ascii_info -mat_view ::ascii_info -ksp_view_mat ::ascii_info -ksp_view_rhs ::ascii_info -ksp_view_solution ::ascii_info'
options={'A':common+' -ksp_gmres_restart 100 -pc_type asm -pc_asm_overlap 2 -sub_ksp_type preonly -sub_pc_type ilu -sub_pc_factor_levels 2',
 'B':common+' -ksp_gmres_restart 100 -pc_type ilu -pc_factor_levels 2',
 'C':common+' -ksp_gmres_restart 100 -pc_type ilu -pc_factor_levels 1'}
policy=dict(mpi_ranks=1,GPUs=1,OMP_NUM_THREADS=1,production_policy=production,production_policy_path='configs/sv1_3/policy.json',production_policy_sha256=sha256(ROOT/'configs/sv1_3/policy.json'),VTU_cadence=10,restart_cadence=10,force_final_VTU=True,force_final_checkpoint=True,window_start=60,window_end=70,checkpoint_sha256=checkpoint['sha256'],candidate_options=options,max_candidates=['A','B','C','D'],D='Current winner PC, GMRES restart200',C_eligibility='B is unhealthy/not run or B improves A by less than10%; optional high-PC-cost branch not needed',direct_accept_reduction=.10,no_gain_threshold=.05,confirmation='Only one extra timing in5%-10% improvement band. Accept if both faster and arithmetic-mean reduction>=5%; otherwise retain incumbent. No median.',profiling='One additional10-step lightweight PETSc log_view for baseline and final winner, shared if winnerA. Profile runs excluded from selection.',scientific_equivalence='DEFERRED',CPU_runs=False,multi_rank_cuda=False,full_steady_runs=1)
write(C/'policy.json',policy)
write(R/'reference_freeze.json',dict(status='PASS',timestamp=now(),Stage_N_delivery_sha256=sha256(N/'delivery_manifest.json'),Stage_N_files_verified=len(delivery['files']),GPU_baseline=candidate,baseline_runtime=baseline,PETSc=read(N/'petsc325_source.json'),PETSc_build=read(N/'petsc_gpu13_build.json'),svmp=read(N/'svmp_gpu_build.json'),CPU_production='CPU_EARLY_STOP_PRODUCTION',science=read(ROOT/'reports/sv1_3/reference_freeze.json'),production_policy_sha256=policy['production_policy_sha256'],scientific_equivalence='DEFERRED',official_source_git=git_state(ROOT/'external/svMultiPhysics')))
print(json.dumps(dict(frozen_history_entries=len(files),Stage_N_files=len(delivery['files']),checkpoint=checkpoint['sha256'],source_patch_files=[name])))
