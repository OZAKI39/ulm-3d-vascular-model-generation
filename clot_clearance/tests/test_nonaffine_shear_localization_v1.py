"""The eleven requested shear checks; no prior campaigns or equilibrium reruns."""
import json
from pathlib import Path
import numpy as np
from pd_clot.pfpd_candidate.baqp import evaluate
from pd_clot.pfpd_candidate.energy import Material, constitutive
from pd_clot.pfpd_candidate.nonlocal_boundary.geometry import block_family
from pd_clot.pfpd_candidate.nonlocal_boundary.nonaffine_manufactured import QuadraticShear
from pd_clot.pfpd_candidate.nonlocal_boundary.field_bias import bias_samples, bias_summary

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'verification/nonlocal_boundary_nonaffine_shear_localization_v1'


def read(name): return json.loads((OUT/name).read_text())
def arrays(): return np.load(OUT/'NODEWISE_SHEAR_DECOMPOSITION.npz')
def maxnorm(x): return np.linalg.norm(x, axis=-1).max()
def scale(): return read('EXACT_NONAFFINE_SHEAR_REFERENCE.json')['scales']['F_ref_N']
def material(): return Material(model='svk', mu=1000., bulk=666.6666666666666)


def test_01_closed_form_shear_F_E_S_P_W():
    f=QuadraticShear();z=np.array([-.5,-.2,0.,.3,.5])*f.length;r=f.response(z)
    expected=np.broadcast_to(np.eye(3),(5,3,3)).copy();expected[:,0,2]=[.002,.0026,.003,.0036,.004]
    np.testing.assert_allclose(r['F'],expected,rtol=0,atol=1e-18)
    np.testing.assert_allclose(r['E'],.5*(expected.swapaxes(1,2)@expected-np.eye(3)),rtol=0,atol=1e-15)
    np.testing.assert_allclose(r['S'],2*f.mu*r['E'],rtol=0,atol=0)
    W,P=constitutive(expected,material())
    np.testing.assert_allclose(r['P'],P,rtol=2e-12,atol=2e-13)
    np.testing.assert_allclose(r['W'],W,rtol=2e-13,atol=1e-17)
    np.testing.assert_array_equal(np.linalg.det(r['F']),np.ones(5))


def test_02_body_force_is_negative_P_column_derivative():
    f=QuadraticShear();z=np.linspace(-.5,.5,7)*f.length;dz=1e-5*f.length
    derivative=(f.response(z+dz)['P'][:,:,2]-f.response(z-dz)['P'][:,:,2])/(2*dz)
    np.testing.assert_allclose(-f.response(z)['b'],derivative,rtol=1e-10,atol=1e-10)
    assert np.all(f.response(z)['b'][:,2]<0)


def test_03_continuum_force_and_all_six_face_channels():
    f=QuadraticShear();channels={c:f.continuum(c) for c in ['combined','normal','tangential']}
    for c in channels.values():
        assert c['force_closure_scaled']<=1e-10
        assert {v['face'] for v in c['faces']}=={'x-','x+','y-','y+','z-','z+'}
    for key in ['surface_resultant_N','body_resultant_N']:
        np.testing.assert_allclose(np.asarray(channels['normal'][key])+channels['tangential'][key],channels['combined'][key],rtol=1e-13,atol=1e-20)
    xp=next(v for v in channels['combined']['faces'] if v['face']=='x+')
    assert xp['force_N'][0]>0 and xp['force_N'][2]>0


def test_04_continuum_reference_moment_and_work():
    f=QuadraticShear();c=f.continuum()
    assert c['moment_closure_scaled']<=1e-10 and c['current_moment_closure_scaled']<=1e-10
    # Reference-coordinate torque is not zero for nonsymmetric first Piola P.
    expected=f.mu*(f.gamma**3+f.gamma*f.beta_s**2)*f.length**3
    np.testing.assert_allclose(np.asarray(c['surface_moment_Nm'])+c['body_moment_Nm'],[0,expected,0],rtol=1e-8,atol=1e-23)
    assert np.linalg.norm(c['current_moment_Nm'])<1e-23
    assert c['max_work_closure_scaled']<=1e-10
    assert {r['variation'] for r in c['work']} >= {'translation_x','translation_z','A_xz','A_zz','A_zx','quadratic_x'}


def test_05_affine_shear_B1_limit():
    r=read('AFFINE_SHEAR_LIMIT_RECOVERY.json')
    assert r['PASS'] and r['relative_L2_error']<=1e-10
    assert read('HISTORICAL_INTEGRITY.json')['PASS']


def test_06_affine_bulk_zero_and_nonaffine_oracle_validation():
    rows=read('AFFINE_SHEAR_ORACLE.json');h=read('PREDECLARED_CONFIG.json')['h_m']
    assert len(rows)==3 and all(r['beta']==0 for r in rows)
    assert max(np.linalg.norm(r['q_bulk_N_m3'])*h**3/scale() for r in rows)<=1e-10
    assert all(c['FULL_FORCE_STENCIL_CLOSED'] for c in read('BULK_STENCIL_CLOSURE.json'))
    assert read('BULK_SHEAR_ORACLE_VALIDATION.json')['VERIFIED']
    assert max(r['error_scaled'] for r in read('SELECTIVE_VS_FULL_SHEAR.json'))<=1e-10


def test_07_nodewise_decomposition_and_channel_accounting():
    with arrays() as a:
        rhs=(a['q_bulk']+a['b_local'])*a['volume'][:,None]
        assert len(np.unique(a['volume']))==4
        assert maxnorm(a['f_BC_localized']-(a['q_bulk']*a['volume'][:,None]-a['f_FD']))/scale()<=1e-10
        assert maxnorm(a['f_FD']+a['f_BC_localized']+a['f_body_local']-rhs)/scale()<=1e-10
        assert maxnorm(a['f_FD']+a['f_external_NL'])/scale()<=1e-10
        assert maxnorm(a['f_BC_normal']+a['f_BC_tangential']-a['f_BC_localized'])/scale()<=1e-10
    for c in read('SHEAR_CHANNEL_ACCOUNTING.json')['channels'].values():
        assert c['force']['closure_scaled']<=1e-10 and c['moment']['closure_scaled']<=1e-10


def test_08_naive_cancellation_contains_interior_physics():
    with arrays() as a:
        assert maxnorm(a['f_BC_naive']-a['f_BC_localized']+a['f_bulk_target'])/scale()<=1e-10
    center=read('N27_SHEAR_DEPTH_CHECK.json')['points'][-1]
    assert center['naive_scaled']>1e-5 and center['localized_scaled']<=1e-10


def test_09_N27_closed_center_and_independent_scaffold():
    r=read('N27_SHEAR_DEPTH_CHECK.json');p=r['points'][-1]
    assert r['PASS'] and p['fully_closed'] and p['localized_scaled']<=1e-10
    assert r['closure'][-1]['every_D2_node_has_complete_family']
    assert r['closure'][-1]['levels'][3]['nonbulk_volume_node_count']==0
    assert r['independent_reference']['n']==29


def test_10_NL_perturbed_start_recovery():
    r=read('SHEAR_NL_MANUFACTURED_RECOVERY.json')
    assert r['PASS'] and r['force_residual']<=1e-8 and r['F_recovery_error_max']<=1e-5
    with np.load(OUT/'SHEAR_NL_RECOVERED_STATE.npz') as s, arrays() as a:
        assert abs(s['initial'][:,0]-s['target'][:,0]).max()>1e-8
        np.testing.assert_array_equal(s['initial'][s['fixed']],s['target'][s['fixed']])
        np.testing.assert_array_equal(s['external'],a['f_external_NL'])
        assert np.linalg.norm((s['force_internal']+s['external'])[~s['fixed']])/r['solver_force_scale_N']<=1e-8
        assert np.linalg.norm(s['Fbar']-a['Fbar_target'],axis=(1,2)).max()/.004<=1e-5
        assert maxnorm(s['x']-s['target'])<maxnorm(s['initial']-s['target'])*1e-4


def test_11_Fbar_Ftilde_Psi_diagnostics_are_deterministic():
    summary=read('SHEAR_FIELD_BIAS_SUMMARY.json')
    assert all(summary['deterministic'].values())
    # Fresh small geometry: verify both node and bond-midpoint diagnostics.
    f=QuadraticShear();g,family=block_family(7,6*62.5e-6,251.25e-6)
    st=evaluate(family,f.position(family.X),material());i,j=family.pairs.T
    for F,X,w in [(st.Fbar,family.X,family.volume),(st.Ftilde,(family.X[i]+family.X[j])/2,2*family.wbar*family.volume[i]*family.volume[j])]:
        exact=f.response(X[:,2])['F'];one=bias_samples(F,exact,material());two=bias_samples(F,exact,material())
        for key in one: np.testing.assert_array_equal(one[key],two[key])
        np.testing.assert_array_equal(one['F_error'],np.linalg.norm(F-exact,axis=(1,2)))
        q=bias_summary(one,w,.004,summary['Psi_plus_reference_scale_Pa'])
        assert q['F_error_max_scaled']>0 and q['Psi_plus_error_RMS_Pa']>=0
    with np.load(OUT/'SHEAR_FIELD_BIAS_SAMPLES.npz') as a:
        for name in ['Fbar','Ftilde']:
            samples={k:a[name+'_'+k] for k in ['F_error','Psi_plus_error']}
            metrics=bias_summary(samples,a[name+'_weight'],.004,summary['Psi_plus_reference_scale_Pa'])
            for k,v in metrics.items(): assert v==summary['samples'][name][k]
