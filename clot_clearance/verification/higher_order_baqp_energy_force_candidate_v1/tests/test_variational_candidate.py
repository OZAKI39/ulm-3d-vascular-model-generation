"""The twelve requested focused gates, using saved evidence and live small checks."""
import csv,hashlib,json
from pathlib import Path
import numpy as np
import pytest
from pd_clot.pfpd_candidate.kernel import Family
from pd_clot.pfpd_candidate.baqp import evaluate as historical_evaluate
from pd_clot.pfpd_candidate.energy import Material
from pd_clot.pfpd_candidate.higher_order_rk.mechanics_candidate import Reference,GradientOperator,VariationalBAQP,c1_operator
from pd_clot.pfpd_candidate.higher_order_rk.rk_kinematics import gradient
from pd_clot.pfpd_candidate.higher_order_rk.adjoint_accumulation import delta_bond_gradient
from pd_clot.pfpd_candidate.higher_order_rk.variational_energy import bond_state
from pd_clot.pfpd_candidate.higher_order_rk.validation_states import smooth_direction

OUT=Path(__file__).resolve().parents[1];ROOT=OUT.parents[1];OLD=ROOT/'verification/higher_order_baqp_rk_oracle_v1'
def read(n):return json.loads((OUT/n).read_text())
def rows(n):
    with (OUT/n).open() as f:return list(csv.DictReader(f))
CFG=read('CANDIDATE_CONFIGURATION.json')

@pytest.fixture(scope='module')
def small():
    with np.load(OUT/'N7_REFERENCE_GEOMETRY.npz') as a:
        f=Family.from_arrays(a['X'],a['volume'],a['pairs'],CFG['delta_m'],a['omega']);s,t=a['source'],a['target']
    with np.load(OUT/'N7_PD2_REFERENCE_CACHE.npz') as a:g=a['g']
    ref=Reference.from_family(f);p=VariationalBAQP(ref,GradientOperator('PD2',s,t,g))
    c=VariationalBAQP(ref,c1_operator(f,s,t))
    return f,p,c

def test_01_cached_PD2_kinematics_preserved():
    ident=read('PD2_CACHE_IDENTITY.json')
    assert ident['primary_cache_sha256']==CFG['primary_PD2_cache_sha256']
    with np.load(ROOT/CFG['primary_PD2_cache']) as a:
        assert hashlib.sha256(a['g'].tobytes()).hexdigest()==ident['g_array_sha256']
    for name in ['axial','shear','cross']:
        with np.load(OUT/(name.upper()+'_STATE.npz')) as a,np.load(OLD/(name.upper()+'_NODAL_STATE.npz')) as b:
            np.testing.assert_array_equal(a['Fbar'],b['PD2'])
            np.testing.assert_array_equal(a['x'],b['x'])

def test_02_delta_Ftilde_identity_and_constitutive_stress(small):
    f,p,c=small
    with np.load(OUT/'C1_FORCE_REGRESSION.npz') as a:x=a['quadratic_shear_x']
    v=smooth_direction(f.X,CFG['small_L_m'],CFG['direction_seeds'][0]);e=CFG['h_m']*CFG['delta_Ftilde_step_rel']
    dF=gradient(v,p.operator.source,p.operator.target,p.operator.g)
    analytic=delta_bond_gradient(p.reference,v,dF,slice(None))
    plus=bond_state(p.reference,x+e*v,p.nodes(x+e*v),slice(None),p.material)[0]
    minus=bond_state(p.reference,x-e*v,p.nodes(x-e*v),slice(None),p.material)[0]
    error=np.linalg.norm((plus-minus)/(2*e)-analytic,axis=(1,2)).max()/np.linalg.norm(analytic,axis=(1,2)).max()
    assert error<=CFG['delta_Ftilde_gate']
    assert read('CONSTITUTIVE_STRESS_CHECK.json')['max_normalized_error']<=CFG['constitutive_gate']

def test_03_C1_generic_energy(small):
    f,p,c=small
    with np.load(OUT/'C1_FORCE_REGRESSION.npz') as a:
        for name in ['affine','quadratic_shear','seeded_smooth']:
            x=a[name+'_x'];old=historical_evaluate(f,x,p.material)
            assert abs(c.energy(x)-old.energy)/old.energy<=CFG['C1_energy_gate']

def test_04_C1_generic_force(small):
    f,p,c=small
    with np.load(OUT/'C1_FORCE_REGRESSION.npz') as a:
        for name in ['affine','quadratic_shear','seeded_smooth']:
            actual=c.evaluate(a[name+'_x']).force;expected=a[name+'_historical_force']
            e=np.linalg.norm(actual-expected,axis=1);scale=max(4*CFG['h_m']**2,np.linalg.norm(expected,axis=1).max())
            assert e.max()/scale<=CFG['C1_force_max_gate']
            assert np.sqrt(np.average(e*e,weights=f.volume))/scale<=CFG['C1_force_RMS_gate']

def test_05_rigid_translation():
    a=read('OBJECTIVITY_RESULTS.json')['translation']
    assert a['energy_J']<=CFG['rigid_energy_absolute_gate_J'] and a['max_force_N']<=a['force_gate_N']
    assert max(a['node_F_max_error'],a['bond_F_max_error'])<=CFG['rigid_F_max_gate']

def test_06_finite_rotation_and_covariance():
    obj=read('OBJECTIVITY_RESULTS.json');a=obj['rotation']
    assert CFG['rotation_angle_degrees']==37 and a['Green_strain_max']<=CFG['rigid_strain_max_gate']
    assert a['energy_J']<=CFG['rigid_energy_absolute_gate_J'] and a['max_force_N']<=a['force_gate_N']
    assert obj['superposed_rotation']['force_covariance']['max_normalized']<=1e-10

def test_07_global_force_balance():
    for r in rows('GLOBAL_BALANCE_RESULTS.csv'):assert float(r['force_residual'])<=CFG['global_force_gate']
    # Recompute from saved N19 forces rather than relying only on PASS flags.
    for name in ['finite_affine_shear','axial','shear','cross','affine_axial','affine_shear']:
        with np.load(OUT/(name.upper()+'_STATE.npz')) as a:
            force=a['force'];scale=max(4*CFG['h_m']**2,np.linalg.norm(force,axis=1).max())
            assert np.linalg.norm(force.sum(0))/scale<=CFG['global_force_gate']

def test_08_current_configuration_moment():
    for r in rows('GLOBAL_BALANCE_RESULTS.csv'):assert float(r['moment_residual'])<=CFG['global_moment_gate']
    for name in ['finite_affine_shear','axial','shear','cross','affine_axial','affine_shear']:
        with np.load(OUT/(name.upper()+'_STATE.npz')) as a:
            x,force=a['x'],a['force'];scale=max(4*CFG['h_m']**2,np.linalg.norm(force,axis=1).max())
            assert np.linalg.norm(np.cross(x-x.mean(0),force).sum(0))/(scale*CFG['primary_L_m'])<=CFG['global_moment_gate']

def test_09_exactly_six_directional_checks_all_steps():
    rr=rows('ENERGY_FORCE_DIRECTIONAL_DERIVATIVES.csv');assert len(rr)==18
    groups={(r['state'],r['direction']) for r in rr};assert len(groups)==6
    for name,k in groups:
        r=[v for v in rr if (v['state'],v['direction'])==(name,k)]
        assert [float(v['eps_rel']) for v in r]==CFG['FD_eps_rel']
        with np.load(OUT/(name.upper()+'_STATE.npz')) as a,np.load(OUT/f'FD_{name}_{k}_DIRECTION.npz') as d:
            work=-np.sum(a['force']*d['v'])
            assert np.linalg.norm(d['v'].mean(0))<1e-15
        for v in r:
            fd=(float(v['U_plus_J'])-float(v['U_minus_J']))/(2*float(v['eps_m']))
            err=abs(fd-work)/float(v['denominator_N'])
            assert err<=CFG['FD_all_steps_gross_gate']
            if v['middle_step']=='True':assert err<=CFG['FD_middle_gate']

def test_10_affine_energy_and_bondwise_equivalence():
    a=read('AFFINE_ENERGY_COMPARISON.json');assert set(a)=={'affine_axial','affine_shear'}
    for r in a.values():
        assert r['relative_energy_difference']<=CFG['affine_energy_gate']
        assert r['bond_F_difference_over_0p004']<=CFG['affine_bond_F_gate']
        assert r['bond_density_relative_difference']<=CFG['affine_bond_density_gate']

def test_11_quadratic_energy_matches_U1_without_forcing_U2():
    rr=rows('QUADRATIC_ENERGY_DECOMPOSITION.csv');assert len(rr)==3
    for r in rr:
        U=float(r['U_PD2_J']);U1=float(r['U_exactkin_currentweights_J']);U2=float(r['U_continuum_J'])
        assert abs(U-U1)/float(r['energy_scale_J'])<=CFG['quadratic_U1_gate']
        assert float(r['Gauss32_vs64_relative'])<=1e-12
        np.testing.assert_allclose(float(r['remaining_weight_quadrature_difference_J']),U-U2,rtol=1e-12,atol=1e-30)

def test_12_selected_DOF_slow_oracle(small):
    f,p,c=small;rr=rows('SELECTED_DOF_FORCE_ORACLE.csv');assert len(rr)==12
    with np.load(OUT/'N7_SELECTED_STATE.npz') as a:
        force=p.evaluate(a['x']).force
    for r in rr:
        slow=-(float(r['U_plus_J'])-float(r['U_minus_J']))/(2*float(r['eps_m']))
        actual=force[int(r['node']),int(r['component'])]
        assert abs(slow-actual)/float(r['F_ref_N'])<=CFG['selected_DOF_gate']
    assert not p.operator.g.flags.writeable
