"""Freeze all successful baselines and create the isolated within-timestep adapter."""
import difflib,gzip,json,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import inventory,sha256,git_state,now
R=ROOT/'reports/sv1_3p';C=ROOT/'configs/sv1_3p';O=ROOT/'reports/sv1_3o'
for kind in ('reports','configs','outputs','logs','patches','benchmarks'):(ROOT/kind/'sv1_3p').mkdir(parents=True,exist_ok=True)
def read(p):return json.loads(p.read_text())
def write(p,d):p.write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')
assert not (R/'reference_freeze.json').exists()
shutil.copyfile('/home/lzy/.codex/attachments/b8af267b-193e-403b-af29-4e13d30f37c6/pasted-text.txt',R/'USER_REQUEST.txt')
delivery=read(O/'delivery_manifest.json')
for f in delivery['files']:assert sha256(ROOT/f['path'])==f['sha256'],f['path']
oldref=read(O/'reference_freeze.json');baseline=read(O/'PERF_A_acceptance.json');candidate=read(O/'optimized_steady_candidate.json')
scopes=[];files={}
for kind in ('reports','outputs','logs','configs','patches','benchmarks'):
 scopes += [str(p.relative_to(ROOT)) for p in (ROOT/kind).iterdir() if p.is_dir() and p.name!='sv1_3p']
scopes+=['inputs','external/sv13o','external/petsc325/svMultiPhysics-compat']
for name in scopes:files.update({name+'/'+k:v for k,v in inventory(ROOT/name)['files'].items() if '__pycache__' not in k})
for name in ('configs','src','scripts','tests'):
 for k,v in inventory(ROOT/name)['files'].items():
  if 'sv13p' not in k and 'sv1_3p' not in k and '__pycache__' not in k:files[name+'/'+k]=v
(R/'history_baseline.json.gz').write_bytes(gzip.compress(json.dumps(dict(files=files,scopes=scopes,timestamp=now())).encode(),mtime=0))
shutil.copyfile(O/'old_fem_baseline.json.gz',R/'old_fem_baseline.json.gz')
source=ROOT/'external/sv13o/svMultiPhysics';target=ROOT/'external/sv13p/svMultiPhysics-reuse'
assert not target.exists();shutil.copytree(source,target)
before=read(O/'source_patch.json')['after'];name='Code/Source/solver/petsc_impl.cpp';old=(source/name).read_text()
audit=dict(status='PASS',timestamp=now(),source=str(source.relative_to(ROOT)),source_sha256=sha256(source/name),
 KSP_created='petsc_create_linearsolver called from equation initialize, once per equation',PC_created='KSPGetPC during solver initialization',
 KSP_PC_persist_across_linear_solves=True,destruction='petsc_destroy_all only at process cleanup; existing Stage N lifecycle unchanged',
 matrix_update='petsc_set_values -> MatZeroEntries + MatSetValues + assembly on existing Mat',
 repeated_setup='Each petsc_solve calls KSPSetOperators, KSPSetUp and KSPSolve',
 target_policy='REUSE_WITHIN_TIMESTEP',timestep_source='ComMod.cTS, existing integer timestep',
 official_API='KSPSetReusePreconditioner(ksp, PETSC_FALSE on first solve of a timestep; PETSC_TRUE subsequently)',
 source_evidence={'KSPCreate':old[:old.index('KSPCreate(')].count('\n')+1,'KSPDestroy':old[:old.index('KSPDestroy(')].count('\n')+1,'solve_method':old[:old.index('void PetscLinearAlgebra::PetscImpl::solve(')].count('\n')+1},
 official_source='https://petsc.org/release/manualpages/KSP/KSPSetReusePreconditioner/')
write(R/'ksp_lifecycle_audit.json',audit)
needle='    consts::PreconditionerType preconditioner_;'
new=old.replace(needle,needle+'''

    // Stage P: only the PC strategy changes; KSP and matrix ownership stay intact.
    PetscBool reuse_within_timestep_ = PETSC_FALSE;
    int last_pc_timestep_ = -1;
    int pc_rebuild_requests_ = 0;
    int pc_solve_requests_ = 0;
''')
needle='      phys, equation.dof, eq_num, com_mod.nEq);'
assert new.count(needle)==1
new=new.replace(needle,needle+'''
  PetscCallAbort(MPI_COMM_WORLD, PetscOptionsGetBool(nullptr, nullptr,
      "-sv_pc_reuse_within_timestep", &reuse_within_timestep_, nullptr));
''')
needle='  petsc_solve(&lEq.FSILS.RI.fNorm'
assert new.count(needle)==1
addition='''  if (reuse_within_timestep_) {
    const PetscBool reuse = last_pc_timestep_ == com_mod.cTS ? PETSC_TRUE : PETSC_FALSE;
    PetscCallAbort(MPI_COMM_WORLD, KSPSetReusePreconditioner(psol[com_mod.cEq].ksp, reuse));
    last_pc_timestep_ = com_mod.cTS;
    ++pc_solve_requests_;
    if (!reuse) ++pc_rebuild_requests_;
    PC pc = nullptr;
    PetscCallAbort(MPI_COMM_WORLD, KSPGetPC(psol[com_mod.cEq].ksp, &pc));
    PetscCallAbort(MPI_COMM_WORLD, PetscPrintf(MPI_COMM_WORLD,
        "SV13P_REUSE step=%d equation=%d reuse=%d rebuild_requests=%d solve_requests=%d ksp=%p pc=%p\\n",
        com_mod.cTS, com_mod.cEq, static_cast<int>(reuse), pc_rebuild_requests_, pc_solve_requests_,
        static_cast<void*>(psol[com_mod.cEq].ksp), static_cast<void*>(pc)));
  }

'''
new=new.replace(needle,addition+needle);(target/name).write_text(new)
patch=''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='a/'+name,tofile='b/'+name));p=ROOT/'patches/sv1_3p/reuse_within_timestep.patch';p.write_text(patch)
after={n:sha256(target/n) for n in before};assert [n for n in before if before[n]!=after[n]]==[name]
write(R/'source_patch.json',dict(status='PASS',source=str(target.relative_to(ROOT)),base='Frozen Stage O',before=before,after=after,changed_files=[name],patch=str(p.relative_to(ROOT)),patch_sha256=sha256(p),scientific_operations_unchanged=True))
policy=read(ROOT/'configs/sv1_3o/policy.json');base_options=baseline['PETSC_OPTIONS']
write(C/'policy.json',dict(MPI_ranks=1,GPUs=1,OMP_NUM_THREADS=1,window_start=60,window_end=70,checkpoint_sha256=policy['checkpoint_sha256'],production_policy=policy['production_policy'],production_policy_path=policy['production_policy_path'],production_policy_sha256=policy['production_policy_sha256'],VTU_cadence=10,restart_cadence=10,baseline_options=base_options,candidate_order=['P1','P2','P3','P4','P5'],P1_options=base_options+' -sv_pc_reuse_within_timestep true',profiling='Plain PETSc log_view collected in each run; no extra profiling CFD, no GPU synchronization timing flag. Event times inclusive.',early_slow_reject=dict(min_completed_steps=3,projected_wall_ratio=1.5,no_meaningful_iteration_improvement_ratio=.95),winner_threshold=.10,maximum_full_steady_runs=1,scientific_equivalence='DEFERRED',no_repeats=True))
write(R/'reference_freeze.json',dict(status='PASS',timestamp=now(),Stage_O_delivery_sha256=sha256(O/'delivery_manifest.json'),Stage_O_files_verified=len(delivery['files']),GPU_baseline=candidate,window_baseline=baseline,baseline_profile=read(O/'profile_summary.json'),PETSc=oldref['PETSc'],PETSc_build=oldref['PETSc_build'],svmp=read(O/'svmp_build.json'),science=oldref['science'],CPU_production='CPU_EARLY_STOP_PRODUCTION',production_policy_sha256=policy['production_policy_sha256'],scientific_equivalence='DEFERRED',official_source_git=git_state(ROOT/'external/svMultiPhysics')))
write(R/'checkpoint_window.json',read(O/'checkpoint_window.json'))
print(f'Stage O frozen: {len(delivery["files"])} files. Lifecycle audited; one-file P1 adapter prepared.')
