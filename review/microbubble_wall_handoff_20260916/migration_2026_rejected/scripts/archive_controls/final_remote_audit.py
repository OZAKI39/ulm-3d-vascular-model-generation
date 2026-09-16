from pathlib import Path
import json,hashlib,subprocess,shutil,os
R=Path('/workspace/microbubble_lammps/results/lammps2026_bubble_interaction_20260916_091244');W=Path('/workspace/lammps_migration/new_20260916_091244')
assert json.loads((R/'MIGRATION_EXECUTION_STATE_V2.json').read_text())['state']=='REGRESSIONS_PASS_SOURCE_GATE_FAIL'
py='/workspace/microbubble_lammps/work/palabos_lammps_coupling_v0_20260915_224019/venv/bin/python'
p=subprocess.run([py,'-B',str(R/'scripts/finalize_migration_candidate.py'),'--root',str(R),'--output-dir',str(R/'validation/remote_independent')],capture_output=True,text=True,timeout=180);print(p.stdout);print(p.stderr);assert p.returncode==0
final=json.loads((R/'validation/remote_independent/INDEPENDENT_MIGRATION_EVALUATION.json').read_text());assert final['status']=='FAIL' and all(v['status']=='PASS'for v in final['regressions'].values())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
custom={}
for stage,prefix in [('coupling','coupled'),('passive','passive')]:
 custom[stage]={'entrypoint_and_math_source_sha256':{str(p.relative_to(R/'migration'/stage)):sha(p)for p in (R/'migration'/stage/'src').glob('*') if p.is_file() and p.suffix in ['.cpp','.hpp']},'binaries':{p.name:sha(p)for p in (W/stage/'build').glob('*')if p.is_file() and(p.name.startswith(prefix+'_lmp_')or p.name=='libfrozen_flow_coupling.so')},'linked_lammps':{kind:sha(W/('build_'+kind)/'liblammps.a')for kind in ['cpu','gpu']}}
(R/'provenance/CUSTOM_BUILD_FINAL.json').write_text(json.dumps(custom,indent=2)+'\n')
oldroot=Path('/workspace/microbubble_lammps/work/lammps_particle_engine_20260915_214759');before=json.loads((R/'provenance/OLD_ACTIVE_BEFORE.json').read_text());assert all(sha(Path(a['path']))==a['sha256']for a in before)
ret={'status':'NOT_EXECUTED_PHASE_A_SOURCE_GATE_FAILED','removal_authorized_by_gate':False,'files_deleted':0,'entries':[],'old_release':'22Jul2025 Update 6','old_commit':'9c5ab448c78a14fd534619622162ba418d6a1fb1','reason':'User section 19 forbids deleting old installation when migration fails.'}
(R/'LAMMPS_2025_RETIREMENT_MANIFEST.json').write_text(json.dumps(ret,indent=2)+'\n')
a={'status':'OLD_VALIDATED_BASELINE_RETAINED_NO_PROMOTION','ACTIVE_LAMMPS_COUNT':1,'ACTIVE_LAMMPS_RELEASE':'22Jul2025 Update 6','ACTIVE_LAMMPS_COMMIT':'9c5ab448c78a14fd534619622162ba418d6a1fb1','active_definition':'installation used by previously validated scientific pipeline; no new active symlink/PATH binding created','active_root':str(oldroot),'ACTIVE_CPU_BINARY_SHA256':sha(oldroot/'build_cpu/lmp'),'ACTIVE_GPU_BINARY_SHA256':sha(oldroot/'build_gpu/lmp'),'OLD_ACTIVE_LAMMPS_PRESENT':'YES','NEW_ACTIVE_LAMMPS_PRESENT':'NO','OLD_SOURCE_PRESENT_IN_ACTIVE_WORKSPACE':'YES','OLD_BUILD_PRESENT_IN_ACTIVE_WORKSPACE':'YES','OLD_BINARY_PRESENT_IN_ACTIVE_WORKSPACE':'YES','OLD_RUNTIME_LIBRARY_REFERENCES':'RETAINED_NOT_RETIRED','HISTORICAL_OLD_RESULTS_PRESERVED':'YES','LAMMPS_VERSION_COEXISTENCE':'NO_TWO_ACTIVE_INSTALLATIONS','isolated_nonactive_candidate':str(W),'candidate_promoted':False,'active_installation_unchanged':True,'which_lmp':shutil.which('lmp'),'PATH_changed':False,'old_source_binary_fingerprints_match':True,'candidate_quarantine_marker':'NOT_PROMOTED.json'}
(R/'LAMMPS_ACTIVE_INSTALL_AUDIT.json').write_text(json.dumps(a,indent=2)+'\n')
(R/'PHASE_C_NOT_STARTED.json').write_text(json.dumps({'status':'BLOCKED_BY_PHASE_A_SOURCE_CORRECTION_GATE','production_interaction_implemented':False,'bubble_interaction_case_runs':0,'phase_C_numerical_results':None,'phase_C_visualizations':'NOT_GENERATED_NO_NUMERICAL_RESULTS','reference_finalizer_for_phase_C':'NOT_IMPLEMENTED','human_visual_review':'PENDING','no_false_PASS':True},indent=2)+'\n')
for f in ['build_worker.stdout','build_worker.stderr','migration.stdout','migration.stderr','requalification.stdout','requalification.stderr','supervisor.log','supervisord.conf']:
 shutil.copy2(W/f,R/'provenance'/f)
print(subprocess.run(['supervisorctl','-c',str(W/'supervisord.conf'),'status'],capture_output=True,text=True).stdout)
subprocess.run(['supervisorctl','-c',str(W/'supervisord.conf'),'shutdown'],check=True)
print('ALL_TASK_RUNNERS_TERMINAL; OWN_SUPERVISOR_SHUTDOWN')
print(subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,process_name','--format=csv,noheader'],text=True))
