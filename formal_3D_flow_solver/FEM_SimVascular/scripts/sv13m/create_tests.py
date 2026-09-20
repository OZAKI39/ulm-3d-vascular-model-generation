"""Install permanent Stage M acceptance/regression tests once."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
header='from sv13m_support import *\nfrom sv_validation.sv13m import *\n'
cases={
'reference_freeze':'''def test_references_unchanged():
 from sv_validation.provenance import sha256
 d=accepted('reference_manifest')
 assert d['stack']['svmp_gpu_build']['source_unmodified']
 for f in d['files']:assert sha256(ROOT/f['path'])==f['sha256']
 assert d['MPI_prefix']==d['stack']['compatibility_winner']['MPI_prefix']
''',
'ghost_failure_reproduction':'''def test_once_actual_cuda_and_ghost_error():
 d=accepted('baseline_ghost_failure');assert d['run_count']==1 and d['reproduced']
 assert 'seqcuda' in d['vector_types'] and 'seqaijcusparse' in d['matrix_types']
 assert any('VecGhostUpdateBegin' in l for l in d['stack'])
def test_unpatched_cpu_pass_cuda_fail():
 d=read('ghost_probe_before')
 for r in d['runs']:
  if r['variant'].startswith('CPU'):ghost_run_gate(r)
  else:
   with pytest.raises(GateError):ghost_run_gate(r)
''',
'upstream_diff':'''def test_official_history_is_hashed():
 from sv_validation.provenance import sha256
 d=accepted('upstream_audit');assert d['official_repository']=='https://gitlab.com/petsc/petsc.git'
 assert d['first_fixed_commit']=='8ba8e283c7d3ba8c2bd73cc674b79590a073ab5a' and d['first_fixed_tag']=='v3.25.0'
 for c in d['commands']:assert c['exit_code']==0 and sha256(ROOT/c['artifact'])==c['sha256']
def test_local_adaptation_not_misrepresented_as_upstream():
 d=read('repair_03');assert d['upstream_exact_commit'] is None and 'NOT claimed' in d['provenance']
''',
'patch_scope':'''def test_actual_patch_scope():
 d=accepted('patch_integrity');patch_scope_gate(d['modified_files'],d['expected_files'],d['solver_changes'])
@pytest.mark.parametrize('extra,solver',[(['src/ksp/ksp/interface/itfunc.c'],[]),([],['petsc_impl.cpp'])])
def test_unexpected_change_rejected(extra,solver):
 expected=['src/vec/vec/impls/mpi/commonmpvec.c']
 with pytest.raises(GateError):patch_scope_gate(expected+extra,expected,solver)
''',
'patch_provenance':'''def test_provenance_complete_and_bounded():
 d=accepted('repair_iterations');assert len(d['iterations'])<=3
 for r in d['iterations']:
  for k in ('trigger','root_cause','modified_files','modified_functions','added_lines','deleted_lines','tests_added','before_result','after_result'):assert k in r
 text=(ROOT/'patches/sv1_3m/PATCH_PROVENANCE.md').read_text();assert 'NONE — compatibility only' in text and 'Revert:' in text
''',
'petsc_clean_build':'''def test_fresh_source_and_original_configuration():
 d=accepted('petsc_build');assert all(d['selftest_observed'].values())
 assert len(d['configuration_diff_from_L'])==2
 assert all(x['after'].startswith(('--prefix=','PETSC_ARCH=')) for x in d['configuration_diff_from_L'])
 for name in ('configure','make','install','selftest'):assert d[name]['exit_code']==0
 integrity=accepted('patch_integrity');assert integrity['fresh_source'] and not integrity['previous_objects_reused']
''',
'ghost_forward':'''def test_all_forward_data():
 ghost_gate(read('ghost_probe_after'))
def test_silent_corruption_rejected():
 d=ghost();d['values'][0][-2]='999'
 with pytest.raises(GateError):ghost_run_gate(d)
def test_falsified_expected_values_rejected():
 d=ghost();d['values'][0][-2:]=['999','999']
 with pytest.raises(GateError):ghost_run_gate(d)
''',
'ghost_reverse':'''def test_reverse_including_nonzero_remote_contributions():
 d=ghost();ghost_run_gate(d)
 assert any(ph=='reverse' and float(a)>100 for ph,rank,i,a,e in d['values'])
def test_forward_pass_reverse_fail_rejected():
 d=ghost();d['checks'][0][2]='1'
 with pytest.raises(GateError):ghost_run_gate(d)
def test_repair02_error_suppression_not_accepted():
 d=next(r for r in read('ghost_probe_repair_02')['runs'] if r['variant']=='CUDA_MPI')
 with pytest.raises(GateError):ghost_run_gate(d)
''',
'ghost_local_form':'''def test_actual_device_arithmetic_and_local_writeback():
 d=accepted('ghost_coherence');assert len(d['runs'])==6 and all(r['accepted'] for r in d['runs'])
def test_wrong_layout_rejected():
 d=ghost();d['locals'][0][2]='4'
 with pytest.raises(GateError):ghost_run_gate(d)
''',
'petsc_gpu_regression':'''def test_three_same_gpu_sparse_solves():
 d=accepted('petsc_gpu_smoke');old=read('reference_manifest')['stack']['petsc_gpu_smoke']['runs'][0]
 assert len(d['runs'])==3
 for r in d['runs']:
  assert r['accepted'] and r['mat_type']=='seqaijcusparse' and r['vec_type']=='seqcuda'
  assert r['iterations']==old['iterations']
  assert abs(r['true_relative_residual']-old['true_relative_residual'])<1e-14
  assert abs(r['solution_error_inf']-old['solution_error_inf'])<1e-14
''',
'svmp_link':'''def test_correct_link_and_clean_solver_source():
 d=accepted('svmp_gpu_link');assert d['resolved_PETSc']==d['expected_PETSc'] and 'ghostfix' in d['resolved_PETSc']
 assert d['resolved_MPI']==d['expected_MPI'] and d['source_unmodified']
''',
'official_gpu_smoke':'''def test_official_complete_flow():flow_gate(accepted('svmp_gpu_smoke'))
def test_missing_vtu_rejected():
 d=good_flow();d['VTU_count']=0
 with pytest.raises(GateError):flow_gate(d)
''',
'petcs_error_precedence':'''def test_later_ksp_convergence_never_overrides_petsc_error():
 d=good_flow();d.update(PETSc_error_detected=True,petsc_reason='CONVERGED_RTOL')
 with pytest.raises(GateError,match='PETSC_HARD_ERROR'):flow_gate(d)
def test_original_failure_contains_misleading_convergence():
 d=read('baseline_ghost_failure');assert d['history']['petsc_reasons'] and d['status']=='PASS' and d['exit_code']!=0
''',
'gpu_proof':'''def test_real_twenty_steps():flow_gate(accepted('GPU_PROOF_20_acceptance'),steps=20)
def test_nan_fields_rejected():
 d=good_flow();d['velocity_finite']=False
 with pytest.raises(GateError):flow_gate(d,steps=20)
''',
'cpu_proof':'''def test_matching_cpu_proof():
 d=accepted('CPU_PROOF_20_1R_acceptance');flow_gate(d,gpu=False,steps=20)
 gpu=accepted('GPU_PROOF_20_acceptance');assert d['xml_sha256']==gpu['xml_sha256']
 assert d['command'][-2]==gpu['command'][-2]
''',
'science_equivalence':'''def test_actual_science():science_gate(accepted('science_equivalence'))
@pytest.mark.parametrize('field',list(SCIENCE_LIMITS))
def test_outside_frozen_tolerance_rejected(field):
 d=dict(pressure_shift_applied=False,errors={k:0 for k in SCIENCE_LIMITS});d['errors'][field]=1.01*SCIENCE_LIMITS[field]
 with pytest.raises(GateError):science_gate(d)
def test_pressure_alignment_forbidden():
 with pytest.raises(GateError):science_gate(dict(pressure_shift_applied=True,errors={k:0 for k in SCIENCE_LIMITS}))
''',
'gpu_residency':'''def test_actual_profile_separate_from_benchmark():
 d=accepted('gpu_residency');assert d['profiling_case'].startswith('PROFILE_') and d['evidence']
 assert d['peak_VRAM_MiB']>0
''',
'ghost_transfer':'''def test_transfer_attribution_explicit():
 d=accepted('ghost_transfer');assert all(k in d for k in ('ghost','MatMult','PC','growth_relation','limitations'))
''',
'benchmark':'''def test_observed_repeatability():
 d=accepted('benchmark')
 for runs in d['runs'].values():timing_gate(runs)
 assert d['speedup']==min(d['medians']['CPU1'],d['medians']['CPU4'])/d['medians']['GPU1']
def test_single_sample_no_conclusion():
 with pytest.raises(GateError):timing_gate([dict(wall_time_s=1)])
def test_profile_timing_rejected():
 row=dict(accepted=True,steps_completed=20,initial_state='t=0',wall_time_s=1,profiling=True)
 with pytest.raises(GateError):timing_gate([row,row])
def test_variability_requires_third_and_median_uses_all():
 row=dict(accepted=True,steps_completed=20,initial_state='t=0',profiling=False)
 with pytest.raises(GateError):timing_gate([dict(row,wall_time_s=10),dict(row,wall_time_s=12)])
 assert timing_gate([dict(row,wall_time_s=10),dict(row,wall_time_s=12),dict(row,wall_time_s=11)])==11
''',
'visuals':'''def test_required_artifacts_and_not_run_visibility():
 d=accepted('visuals')
 from PIL import Image
 for f in d['figures']:
  image=Image.open(ROOT/f['path']);image.verify()
  assert f['title'] and f['data_source']
 assert len(d['figures'])==(10 if read('svmp_gpu_smoke')['status']=='PASS' else 5)
''',
'history_preservation':'''def test_history_including_stage_l_immutable():
 d=accepted('preservation_audit');assert d['historical']['entries']>=2000 and not d['historical']['changes']
 assert not d['reference_changed'] and d['old_FEM']['status']=='PASS'
def test_remote_historical_stacks_immutable():accepted('remote_preservation')
'''
}
for label,variant in [('cpu_seq_ghost','CPU_seq'),('cuda_seq_ghost','CUDA_seq'),('cpu_mpi_ghost','CPU_MPI'),('cuda_mpi_ghost','CUDA_MPI')]:
 cases[label]=f'''def test_three_actual_runs():
 rows=[r for r in read('ghost_probe_after')['runs'] if r['variant']=={variant!r}]
 assert len(rows)==3
 for r in rows:ghost_run_gate(r)
'''
 if variant.startswith('CPU'):cases[label]+='''def test_cpu_regression_must_be_rejected():
 d=ghost('CPU_MPI');d['values'][0][-2]='-100'
 with pytest.raises(GateError):ghost_run_gate(d)
'''
for name,body in cases.items():
 p=ROOT/'tests'/('test_sv13m_'+name+'.py');assert not p.exists();p.write_text(header+body)
print('Created',len(cases),'permanent test files')
