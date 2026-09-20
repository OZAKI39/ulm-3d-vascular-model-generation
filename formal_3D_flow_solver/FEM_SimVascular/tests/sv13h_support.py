"""Read captured evidence; tests never launch solvers or rewrite historical artifacts."""
import copy, json
from pathlib import Path
import pytest
from sv_validation.provenance import sha256
from sv_validation.sv13h import *
ROOT = Path(__file__).resolve().parents[1]
R = ROOT / 'reports/sv1_3h'
def load(name): return json.loads((R / (name + '.json')).read_text())
def actual(name):
    d = load(name)
    if d['status'] == 'NOT_RUN':
        assert d['executed'] is False and d['measurements'] is None
        pytest.skip(d['reason'])
    assert d['status'] == 'PASS', d.get('reason', d)
    return d
def policy(): return json.loads((ROOT / 'configs/sv1_3h/policy.json').read_text())
def run_record(t):
    return dict(wall_time_s=t, profiling=False, steps=20, initial_state='t=0', science_pass=True)
def smoke_fixture():
    return dict(exit_code=0, mat_type='seqaijcusparse', vec_type='seqcuda', converged_reason=2,
                linear_failures=0, true_relative_residual=1e-13, tolerance=1e-10)
def proof_fixture():
    return dict(executed=True, steps=20, linear_failures=0, nonlinear_failures=0,
                velocity_finite=True, pressure_finite=True, Qin=3e-15,
                Qout={'OUTLET_01':1e-15,'OUTLET_02':1e-15,'OUTLET_03':1e-15},
                mass_error=0., wall_noslip=True, reload_pass=True)
def equivalent_fixture():
    return dict(velocity_relative_L2=1e-8, pressure_relative_L2=1e-8, Qin_relative=1e-9,
                Qout_relative={'OUTLET_01':1e-8,'OUTLET_02':1e-8,'OUTLET_03':1e-8},
                mass_error_difference=1e-9, max_velocity_relative=1e-8)

