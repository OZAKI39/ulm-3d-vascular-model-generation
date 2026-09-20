import pytest
from sv_validation.sv13 import stop_gate,frozen_policy
def fixture():
    return ([{'E_u':1e-8,'E_Q':1e-9} for _ in range(5)],
        {'epsilon_Q':0.,'epsilon_mass':0.,'velocity_finite':True,'pressure_finite':True,'wall_noslip_pass':True})
def test_only_four_intervals_never_stop():
    i,s=fixture();assert not stop_gate(i[:4],s,0,0,frozen_policy())
def test_five_joint_good_intervals_stop():
    i,s=fixture();assert stop_gate(i,s,0,0,frozen_policy())
@pytest.mark.parametrize('key',['epsilon_Q','epsilon_mass'])
def test_steady_with_bad_mass_never_stops(key):
    i,s=fixture();s[key]=2e-6;assert not stop_gate(i,s,0,0,frozen_policy())
@pytest.mark.parametrize('counts',[(1,0),(0,1)])
def test_solver_failure_not_success(counts):
    i,s=fixture();assert not stop_gate(i,s,*counts,frozen_policy())
@pytest.mark.parametrize('key',['velocity_finite','pressure_finite','wall_noslip_pass'])
def test_bad_field_never_stops(key):
    i,s=fixture();s[key]=False;assert not stop_gate(i,s,0,0,frozen_policy())
def test_nonfinite_and_broken_joint_history_rejected():
    i,s=fixture();i[-3]['E_Q']=float('nan');assert not stop_gate(i,s,0,0,frozen_policy())
    i[-3]['E_Q']=2e-6;assert not stop_gate(i,s,0,0,frozen_policy())
