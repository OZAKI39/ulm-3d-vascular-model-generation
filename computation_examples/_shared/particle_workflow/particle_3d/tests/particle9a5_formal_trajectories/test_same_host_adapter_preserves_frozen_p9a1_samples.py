"""Exact parity against the original runner on the same host, not cross-BLAS bytes."""
from pathlib import Path
import importlib.util,sys
import numpy as np

def test_same_host_adapter_preserves_frozen_p9a1_samples(tmp_path,root,core,real_env,parallel_subset,monkeypatch):
 path=root/'particle_3d/reports/network_derived_flow_mb_validation_v1/scripts/runner.py'
 spec=importlib.util.spec_from_file_location('p9a5_test_frozen_runner',path)
 runner=importlib.util.module_from_spec(spec);sys.modules[spec.name]=runner;spec.loader.exec_module(runner)
 e=core['events'][0];sample=real_env.field.sample(e['birth_center_m'])
 runner.ENV=real_env;runner.ROOT=root;runner.IDENTITY={'role':'FROZEN_REFERENCE_TEST_ONLY'}
 runner.INITIAL={str(e['particle_id']):dict(position_m=e['birth_center_m'],q=e['q'],velocity_m_s=sample.velocity_m_s.tolist(),omega_s_inv=(.5*sample.vorticity_s_inv).tolist())}
 from particle_3d.formal_cohort_p9a5 import DT
 from particle_3d.particle6_stepper import bind_query_dependency
 import particle_3d.particle82a_integration as legacy
 reference=bind_query_dependency(legacy.integrate_admitted,{})
 reference.__defaults__=(DT,None,False)
 monkeypatch.setattr(legacy,'integrate_admitted',reference)
 runner.mb_job((e,str(tmp_path)))
 old=np.load(tmp_path/'trajectories/mb_000001.npz')['samples']
 new=np.load(parallel_subset[0]['target']/'1/trajectory.npz')['samples']
 assert np.array_equal(old,new)
