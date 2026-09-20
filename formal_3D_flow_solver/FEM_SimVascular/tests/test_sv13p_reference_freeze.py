from sv13p_support import *
def test_frozen_production_and_stage_o():
 d=accepted('reference_freeze');assert d['CPU_production']=='CPU_EARLY_STOP_PRODUCTION'
 assert d['window_baseline']['wall_time_s']==pytest.approx(260.169647604)
 assert d['window_baseline']['statistics']['total_iterations']==6603
 assert d['GPU_baseline']['wall_time_s']==pytest.approx(2430.023319)
 assert d['PETSc_build']['CUDA_version']=='13.2' and d['PETSc_build']['version']=='3.25.5'
 assert d['production_policy_sha256']==hashlib.sha256((ROOT/'configs/sv1_3/policy.json').read_bytes()).hexdigest()
 assert d['scientific_equivalence']=='DEFERRED'
def test_only_pc_adapter_changes_source():
 d=accepted('source_patch');assert d['changed_files']==['Code/Source/solver/petsc_impl.cpp']
 assert [n for n in d['before'] if d['before'][n]!=d['after'][n]]==d['changed_files']
 for n,h in d['after'].items():assert hashlib.sha256((ROOT/d['source']/n).read_bytes()).hexdigest()==h
