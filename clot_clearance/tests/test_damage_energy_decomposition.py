import json
from pathlib import Path
import numpy as np
import pytest
from pd_clot.regularization.equilibrium import make_fixture
from pd_clot.regularization.damage_energy_audit import audit_state,decompose,counterfactuals

ROOT=Path(__file__).resolve().parents[1]


def historical():
    p=ROOT/'results/streaming_erosion_equilibrium_v1/dt_erosion_0p01'
    c=json.loads((p/'CONFIG.json').read_text());z=np.load(p/'damage_jump_001/state.npz')
    return c,make_fixture(c),z,json.loads((p/'damage_jump_001/AUDIT.json').read_text())


@pytest.mark.parametrize('factor',[0.,.5,1.,2.])
def test_exact_decomposition_with_stabilization_sensitivity(factor):
    c,cl,z,_=historical();c['material']['stabilization_alpha']*=factor
    a=audit_state(cl,z['x'],c['material'],c['safety'],D=z['D_before']);b=audit_state(cl,z['x'],c['material'],c['safety'],D=z['D_after'])
    p,s=decompose(cl,a,b)
    assert abs(s['decomposition_residual_J'])<1e-25
    assert np.max(abs(p['decomposition_residual']))<1e-27
    if factor==0:assert s['delta_U_stabilization_J']==0


def test_historical_damage_jump_reproduces_signed_energy():
    c,cl,z,original=historical();a=audit_state(cl,z['x'],c['material'],c['safety'],D=z['D_before']);b=audit_state(cl,z['x'],c['material'],c['safety'],D=z['D_after'])
    _,s=decompose(cl,a,b)
    assert s['delta_U_total_J']==original['delta_U_damage_J'] and s['delta_U_total_J']<0
    assert a.global_energy['U_total_J']==original['U_before_damage_J']
    assert b.global_energy['U_total_J']==original['U_after_damage_J']


def test_freeze_F_counterfactual_equals_softening_only():
    c,cl,z,_=historical();a=audit_state(cl,z['x'],c['material'],c['safety'],D=z['D_before']);b=audit_state(cl,z['x'],c['material'],c['safety'],D=z['D_after'])
    _,s=decompose(cl,a,b);f=counterfactuals(cl,z['x'],1-z['D_before'],1-z['D_after'],c['material'],c['safety'])
    assert 'NOT USED IN PRODUCTION' in f['label']
    np.testing.assert_allclose(f['freeze_F_constitutive_only_delta_U_J'],s['delta_U_softening_J'],rtol=1e-9,atol=1e-27)
