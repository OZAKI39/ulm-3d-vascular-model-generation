"""Only the nine requested nonaffine localization checks; no scientific reruns."""
import json
from pathlib import Path
import numpy as np
from pd_clot.pfpd_candidate.energy import Material, constitutive
from pd_clot.pfpd_candidate.nonlocal_boundary.nonaffine_manufactured import QuadraticAxial

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'verification/nonlocal_boundary_nonaffine_localization_v1'


def read(name): return json.loads((OUT/name).read_text())
def arrays(): return np.load(OUT/'NODEWISE_RESIDUAL_DECOMPOSITION.npz')
def maxnorm(x): return np.linalg.norm(x, axis=-1).max()
def scale(): return read('EXACT_NONAFFINE_REFERENCE.json')['scales']['F_ref_N']


def test_quadratic_axial_closed_form_reference():
    field=QuadraticAxial();z=np.array([-.5,-.2,0.,.3,.5])*field.length;r=field.response(z)
    W,P=constitutive(r['F'],Material(model='svk',mu=1000.,bulk=666.6666666666666))
    np.testing.assert_allclose(r['P'],P,rtol=2e-13,atol=2e-13)
    np.testing.assert_allclose(r['W'],W,rtol=2e-14,atol=1e-17)
    np.testing.assert_allclose(r['E'],.5*(r['F'].swapaxes(-1,-2)@r['F']-np.eye(3)),rtol=0,atol=1e-15)
    np.testing.assert_allclose(r['S'],2*field.mu*r['E'],rtol=0,atol=0)
    dz=1e-7*field.length
    derivative=(field.response(z+dz)['P'][:,2,2]-field.response(z-dz)['P'][:,2,2])/(2*dz)
    np.testing.assert_allclose(-r['b'][:,2],derivative,rtol=2e-6)
    np.testing.assert_allclose(r['F'][:,2,2],[1.002,1.0026,1.003,1.0036,1.004],rtol=0,atol=3e-16)
    assert (r['J']>0).all()


def test_local_continuum_balance_and_polynomial_work():
    c=QuadraticAxial().continuum()
    assert c['force_closure_scaled']<1e-12
    assert c['max_work_closure_scaled']<1e-12
    assert {r['variation'] for r in c['work']} >= {'translation_x','translation_z','A_zz','A_xz','quadratic_z'}


def test_affine_limit_recovers_previous_B1():
    r=read('AFFINE_LIMIT_RECOVERY.json')
    assert r['PASS'] and r['relative_L2_error']<=1e-10
    assert read('HISTORICAL_INTEGRITY.json')['PASS']


def test_bulk_oracle_zero_affine_and_closed_dependencies():
    rows=read('BULK_ORACLE_RAW.json');h=read('PREDECLARED_CONFIG.json')['h_m']
    affine=[r for r in rows if r['label']=='affine']; assert len(affine)==3
    assert max(np.linalg.norm(r['q_bulk_N_m3'])*h**3/scale() for r in affine)<=1e-10
    assert all(c['FULL_FORCE_STENCIL_CLOSED'] for c in read('BULK_STENCIL_CLOSURE.json'))
    assert read('BULK_ORACLE_VALIDATION.json')['VERIFIED']
    assert max(r['difference_scaled'] for r in read('SELECTIVE_VS_FULL_NONAFFINE.json'))<=1e-10


def test_localized_correction_identity_and_actual_volumes():
    with arrays() as a:
        assert maxnorm(a['f_BC_localized']-(a['q_bulk']*a['volume'][:,None]-a['f_FD']))/scale()<=1e-10
        assert len(np.unique(a['volume']))==4
        assert maxnorm(a['f_bulk_target'])>1e-4*scale()


def test_naive_minus_localized_retains_bulk_force():
    with arrays() as a:
        assert maxnorm(a['f_BC_naive']-a['f_BC_localized']+a['f_bulk_target'])/scale()<=1e-10
    center=read('N27_DEPTH_CHECK.json')['points'][-1]
    assert center['naive_scaled']>1e-5
    assert center['localized_scaled']<=1e-10


def test_nodewise_corrected_residual_equals_bulk_model():
    with arrays() as a:
        rhs=(a['q_bulk']+a['b_local'])*a['volume'][:,None]
        assert maxnorm(a['f_FD']+a['f_BC_localized']+a['f_body_local']-rhs)/scale()<=1e-10
        assert maxnorm(a['f_FD']+a['f_external_NL'])/scale()<=1e-10
    assert read('FORCE_BALANCE.json')['closure_scaled']<=1e-10
    assert read('MOMENT_BALANCE.json')['closure_scaled']<=1e-10


def test_N27_fully_closed_center_boundary_correction_zero():
    r=read('N27_DEPTH_CHECK.json');p=r['points'][-1]
    assert r['PASS'] and p['fully_closed'] and p['localized_scaled']<=1e-10
    assert r['closure'][-1]['every_D2_node_has_complete_family']
    assert r['closure'][-1]['levels'][3]['nonbulk_volume_node_count']==0
    # Independent N29 shifted scaffold avoids simply subtracting an identical call.
    assert r['independent_reference']['n']==29


def test_NL_recovery_from_nonzero_perturbed_state():
    r=read('NL_MANUFACTURED_RECOVERY.json')
    assert r['PASS'] and r['force_residual']<=1e-8 and r['F_recovery_error_max']<=1e-5
    with np.load(OUT/'NL_RECOVERED_STATE.npz') as s, arrays() as a:
        assert maxnorm(s['initial']-s['target'])>1e-8
        residual=s['force_internal']+s['external']
        assert np.linalg.norm(residual[~s['fixed']])/r['solver_force_scale_N']<=1e-8
        err=np.linalg.norm(s['Fbar']-a['Fbar_target'],axis=(1,2)).max()/.004
        assert err<=1e-5
        assert maxnorm(s['x']-s['target'])<maxnorm(s['initial']-s['target'])*1e-4
