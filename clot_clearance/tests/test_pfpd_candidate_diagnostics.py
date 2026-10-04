"""Prescribed pattern diagnostics, normalization, formula and preservation gates."""
import json,hashlib
from pathlib import Path
import numpy as np
import pytest
from scipy.spatial import cKDTree
from pd_clot.pfpd_candidate import Family,evaluate,Material
from pd_clot.pfpd_candidate.degradation import energetic,kinematic
from pd_clot.pfpd_candidate.normalization import c0,coupon,crack_energy,bulk_planar_ratio
ROOT=Path(__file__).resolve().parents[1];O=ROOT/'verification/pfpd_candidate_v1'


def test_energetic_endpoints():np.testing.assert_array_equal(energetic([0,1]),[1,0])

def test_delayed_h_intact_up_to_threshold():np.testing.assert_array_equal(kinematic([0,.0033,.5,.95]),np.ones(4))

def test_delayed_h_fully_broken_zero():assert kinematic(1.)==0

def test_sc_one_disables_even_fully_broken_kinematics():np.testing.assert_array_equal(kinematic([0,.95,1.],1),np.ones(3))

def test_sc_zero_h_equals_g():
    s=np.linspace(0,1,123);np.testing.assert_array_equal(kinematic(s,0),energetic(s))

@pytest.mark.parametrize('bad',[-.01,1.01,np.nan,np.inf])
def test_invalid_phase_proxy_rejected(bad):
    with pytest.raises(ValueError):energetic(bad)
    with pytest.raises(ValueError):kinematic(bad)

@pytest.mark.parametrize('kind,expected',[('constant',3/8),('linear',.3),('quadratic',.25)])
def test_c0_analytic_matches_quadrature(kind,expected):assert abs(c0(kind)['c0']-expected)<1e-14

def test_quadratic_integral_parts():
    q=c0();assert abs(q['numerator']-1/60)<1e-15 and abs(q['denominator']-1/15)<1e-15


def test_discrete_crack_pair_sum_matches_offset_enumeration():
    n=6;h=.000125;X=(np.indices((n,n,n)).reshape(3,-1).T+.5)*h
    f=Family.from_arrays(X,np.full(len(X),h**3),cKDTree(X).query_pairs(2.01*h,output_type='ndarray'),2.01*h)
    i,j=f.pairs.T;s=((X[i,2]-n*h/2)*(X[j,2]-n*h/2)<0).astype(float)
    explicit=crack_energy(f,s);optimized=coupon(n,2.01,side=n*h)['crack_energy_J']
    assert abs(explicit-optimized)/explicit<1e-12
    # Directed sum equals twice the unordered contribution, not half of it.
    directed=np.sum(f.wbar[f.oriented_bond]*f.volume[f.source]*f.volume[f.target]*s[f.oriented_bond]/(.25*f.horizon))
    assert abs(explicit-directed)/explicit<1e-14


def test_discrete_griffith_converges_with_lattice_resolution():
    r=[bulk_planar_ratio(m) for m in [3.01,8.01,16.01,32.01]]
    assert np.all(np.diff(r)>0) and abs(r[-1]-1)<.001


def test_all_35_old_energy_and_reconstruction_audit_reproduced():
    old=json.loads((ROOT/'verification/damage_energy_audit_v1/actual_run_replay/damage_jump_audit.json').read_text())
    rows=json.loads((O/'historical_compatibility/ALL_35_JUMPS.json').read_text());lookup={(r['case'],r['macro_step']):r for r in old}
    assert len(rows)==len(old)==35
    for r in rows:
        ref=lookup[(r['case'],int(r['jump'].split('_')[-1]))]
        assert r['replay_error_J']==0
        assert r['A_delta_U_J']==-ref['delta_U_total_J']
        assert r['A_reconstruction_delta_U_J']==-ref['delta_U_reconstruction_J']
        assert r['A_max_relative_delta_F']==ref['max_relative_F_change']
        assert r['A_max_relative_delta_K']==ref['max_relative_K_change']


def test_all_35_separation_and_BC_equality():
    paths=list((O/'historical_compatibility').glob('dt_erosion_*/damage_jump_*/comparison.npz'));assert len(paths)==35
    for p in paths:
        with np.load(p) as z:
            assert np.max(z['D_after'])<.95
            for field in ['K','F','Ftilde']:
                np.testing.assert_array_equal(z['B_'+field+'_before'],z['B_'+field+'_after'])
                for end in ['before','after']:np.testing.assert_array_equal(z['B_'+field+'_'+end],z['C_'+field+'_'+end])
            assert np.sum(z['B_bond_energy_after']-z['B_bond_energy_before'])<0


def test_historical_prescribed_state_recomputed_and_read_only():
    from pd_clot.regularization.equilibrium import make_fixture
    c=json.loads((ROOT/'configs/streaming_erosion_equilibrium_guarded_v1.json').read_text());cl=make_fixture(c);f=Family.from_cloud(cl)
    p=ROOT/'results/streaming_erosion_equilibrium_v1/dt_erosion_0p01/damage_jump_001/state.npz'
    with np.load(p) as z:x=z['x'].copy();d0=z['D_before'].copy();d1=z['D_after'].copy()
    for a in [x,d0,d1]:a.setflags(write=False)
    states=[evaluate(f,x,Material('neo_hookean'),g=energetic(d),h=kinematic(d)) for d in [d0,d1]]
    np.testing.assert_array_equal(states[0].Fbar,states[1].Fbar)
    assert states[1].energy<states[0].energy


def test_production_outputs_recomputed_unchanged():
    from pd_clot.regularization.equilibrium import make_fixture
    from pd_clot.regularization.damage_energy_audit import audit_state
    c=json.loads((ROOT/'configs/streaming_erosion_equilibrium_guarded_v1.json').read_text());cl=make_fixture(c)
    p=ROOT/'results/streaming_erosion_equilibrium_v1/dt_erosion_0p01/damage_jump_010'
    prior=json.loads((p/'AUDIT.json').read_text())
    with np.load(p/'state.npz') as z:
        for d,key in [('D_before','U_before_damage_J'),('D_after','U_after_damage_J')]:
            s=audit_state(cl,z['x'],c['material'],c['safety'],D=z[d]);assert s.global_energy['U_total_J']==prior[key]


def test_preservation_all_preexisting_files():
    before=json.loads((O/'provenance/PRE_RUN_MANIFEST.json').read_text())['files']
    for rel,info in before.items():
        p=ROOT/rel;assert p.exists(),rel
        assert hashlib.sha256(p.read_bytes()).hexdigest()==info['sha256'],rel


def test_ambiguous_canonical_evolution_not_exposed():
    gate=json.loads((O/'FORMULA_CONSISTENCY_GATE.json').read_text())
    assert gate['status']=='UNRESOLVED_IN_V1' and gate['canonical_evolution_implemented'] is False
    import pd_clot.pfpd_candidate.degradation as d
    assert not hasattr(d,'evolve') and not hasattr(d,'update_history')
