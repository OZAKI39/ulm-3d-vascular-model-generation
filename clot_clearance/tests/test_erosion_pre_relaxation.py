"""Clock separation, freeze invariants and failure gates through real mechanics."""
import copy,json
from pathlib import Path
import numpy as np
import pytest
from pd_clot.regularization.erosion_equilibrium_runner import run
from pd_clot.regularization.equilibrium_adr import equilibrate,FrozenForceEvaluator,make_fixture

ROOT=Path(__file__).resolve().parents[1]


def config():
    c=json.loads((ROOT/'configs/streaming_erosion_equilibrium_guarded_v1.json').read_text())
    c['clot']['cells']=[4,3,3];c['streaming_erosion']['number_of_erosion_steps']=2
    return c


@pytest.fixture(scope='module')
def erosion(tmp_path_factory):
    p=tmp_path_factory.mktemp('equilibrated_erosion')/'run';c=config();s=run(c,p)
    return p,c,s


def test_pre_relaxation_damage_stays_zero(erosion):
    p,c,s=erosion;z=np.load(p/'equilibrium_000/state.npz');np.testing.assert_array_equal(z['bond_D'],0)
    h=np.load(p/'equilibrium_000/history.npz');np.testing.assert_array_equal(h['maximum_damage'],0)
    assert s['pre_relaxation']['converged']


def test_pre_relaxation_exposure_and_cycles_stay_zero(erosion):
    h=np.load(erosion[0]/'equilibrium_000/history.npz')
    np.testing.assert_array_equal(h['erosion_exposure_time_s'],0)
    np.testing.assert_array_equal(h['represented_ultrasound_cycles'],0)
    assert h['equilibrium_pseudo_time'][-1]>0


def test_each_update_uses_equilibrated_driver_and_only_slow_clock(erosion):
    p,c,s=erosion;assert s['status']=='COMPLETED_EQUILIBRATED_EROSION_DEVELOPMENT_RUN'
    z=np.load(p/'states.npz');np.testing.assert_array_equal(z['erosion_exposure_time_s'],[0,.01,.02])
    assert all(z['equilibrium_passed'])
    for j in json.loads((p/'damage_jumps.json').read_text()):
        assert j['equilibrium_passed_before_update'] and j['equilibrium_relative_force_residual']<=1e-4


def test_signed_energy_jump_recomputed_at_identical_coordinates(erosion):
    p,c,s=erosion;cl=make_fixture(c)
    for jp in sorted(p.glob('damage_jump_*')):
        z=np.load(jp/'state.npz');a=json.loads((jp/'AUDIT.json').read_text())
        before=FrozenForceEvaluator(cl,1-z['D_before'],c)(z['x']).stored_energy_J
        after=FrozenForceEvaluator(cl,1-z['D_after'],c)(z['x']).stored_energy_J
        assert before==a['U_before_damage_J'] and after==a['U_after_damage_J']
        assert before-after==a['delta_U_damage_J']


def test_initial_failure_prevents_damage_and_clock(tmp_path):
    c=config();c['mechanical_equilibrium']['maximum_iterations']=1
    p=tmp_path/'failed';s=run(c,p)
    assert s['status']=='INCOMPLETE_MECHANICAL_EQUILIBRIUM'
    assert s['accepted_damage_updates']==0 and s['saved_exposure_s']==0
    z=np.load(p/'states.npz');np.testing.assert_array_equal(z['bond_D'],0)


def test_midrun_failure_prevents_next_damage_update(tmp_path):
    c=config();calls=0
    def fail_second(*args,**kwargs):
        nonlocal calls
        calls+=1;args=list(args);args[4]=copy.deepcopy(args[4])
        if calls==2:args[4]['mechanical_equilibrium']['maximum_iterations']=1
        return equilibrate(*args,**kwargs)
    p=tmp_path/'midrun';s=run(c,p,solver=fail_second)
    assert s['status']=='INCOMPLETE_MECHANICAL_EQUILIBRIUM'
    assert s['accepted_damage_updates']==1 and s['saved_exposure_s']==.01
    z=np.load(p/'states.npz');j=np.load(p/'damage_jump_001/state.npz')
    np.testing.assert_array_equal(z['bond_D'][-1],j['D_after'])


def test_carrier_metadata_does_not_change_erosion(erosion,tmp_path):
    c=copy.deepcopy(erosion[1]);c['damage']['DeltaN']=998877
    c['streaming']['ultrasound_frequency_Hz']=2.5e9;c['simulation']['representative_frequency_Hz']=99000
    p=tmp_path/'metadata';s=run(c,p)
    a=np.load(erosion[0]/'states.npz');b=np.load(p/'states.npz')
    np.testing.assert_array_equal(a['bond_D'],b['bond_D']);np.testing.assert_array_equal(a['x'],b['x'])
