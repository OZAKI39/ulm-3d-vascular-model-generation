"""Only the twelve requested kinematic/energy-driver oracle tests."""
import json
from pathlib import Path
import numpy as np
import pytest
from pd_clot.pfpd_candidate.higher_order_rk.basis import basis,selector
from pd_clot.pfpd_candidate.higher_order_rk.moment_matrix import local_system
from pd_clot.pfpd_candidate.higher_order_rk.rk_kinematics import gradient
from pd_clot.pfpd_candidate.higher_order_rk.manufactured_fields import exact
from pd_clot.pfpd_candidate.higher_order_rk.bond_diagnostic import associated_gradient,psi,metrics

OUT=Path(__file__).resolve().parents[1]
def read(n):return json.loads((OUT/n).read_text())

@pytest.fixture(scope='module')
def data():
    with np.load(OUT/'REFERENCE_GEOMETRY.npz') as a,np.load(OUT/'PD2_REFERENCE_CACHE.npz') as b:
        return {**{k:a[k] for k in a.files},**{k:b[k] for k in b.files}}

def assertF(error,V):
    met=metrics(error,V,.004)
    assert met['max_scaled']<=1e-10 and met['weighted_RMS_scaled']<=1e-11

def test_01_dimensionless_basis_physical_selector_and_scaling(data):
    delta=251.25e-6;v=np.array([[.1,-.2,.3]])
    np.testing.assert_allclose(basis(v),[[.1,-.2,.3,.01,.04,.09,-.02,.03,-.06]])
    np.testing.assert_array_equal(selector(delta)[:3]*delta,np.eye(3))
    # Rescale the same physical cloud and horizon by 2: derivative weights halve.
    k=0;sl=slice(data['ptr'][k],data['ptr'][k+1]);j=data['target'][sl];xi=data['X'][j]-data['X'][k]
    omega=(1-np.linalg.norm(xi,axis=1)/delta)**2
    a=local_system(xi,omega,data['volume'][j],delta)[-1]
    b=local_system(2*xi,omega,8*data['volume'][j],2*delta)[-1]
    np.testing.assert_allclose(2*b,a,rtol=1e-11,atol=1e-10)
    np.testing.assert_allclose(xi.T@a,np.eye(3),rtol=0,atol=1e-11)

def test_02_moment_matrix_symmetry_and_direct_local_sum(data):
    M=data['M'];assert np.min(data['rank'])==9
    assert np.linalg.norm(M-M.transpose(0,2,1))/np.linalg.norm(M)<1e-13
    # Re-sum outer products rather than the implementation's matrix product.
    for k in [0,180,3429,6500]:
        sl=slice(data['ptr'][k],data['ptr'][k+1]);j=data['target'][sl];xi=(data['X'][j]-data['X'][k])/251.25e-6
        Q=basis(xi);w=(1-np.linalg.norm(xi,axis=1))**2*data['volume'][j]
        manual=sum(ww*np.outer(q,q) for ww,q in zip(w,Q))
        np.testing.assert_allclose(manual,M[k],rtol=1e-12,atol=1e-27)

def test_03_all_node_linear_polynomial_reproduction(data):
    g=gradient(data['X'],data['source'],data['target'],data['g'])-np.eye(3)
    for k in range(3):
        m=metrics(g[:,k,:],data['volume'],1.)
        assert m['max_scaled']<=1e-10 and m['weighted_RMS_scaled']<=1e-11

def test_04_all_six_quadratic_terms_at_every_node(data):
    X=data['X'];x,y,z=X.T;p=np.column_stack([x*x,y*y,z*z,x*y,x*z,y*z])
    expected=np.zeros((len(X),6,3));expected[:,0,0]=2*x;expected[:,1,1]=2*y;expected[:,2,2]=2*z
    expected[:,3,0]=y;expected[:,3,1]=x;expected[:,4,0]=z;expected[:,4,2]=x;expected[:,5,1]=z;expected[:,5,2]=y
    got=gradient(p,data['source'],data['target'],data['g'])
    for k in range(6):
        m=metrics(got[:,k]-expected[:,k],data['volume'],18*62.5e-6)
        assert m['max_scaled']<=1e-10 and m['weighted_RMS_scaled']<=1e-11

def affine(name,data):
    x,F=exact(name,data['X'],18*62.5e-6)
    got=np.eye(3)+gradient(x-data['X'],data['source'],data['target'],data['g'])
    assertF(got-F,data['volume'])
def test_05_affine_axial(data):affine('affine_axial',data)
def test_06_affine_shear(data):affine('affine_shear',data)

def quadratic(name,data):
    with np.load(OUT/(name.upper()+'_NODAL_STATE.npz')) as a:
        x,F=exact(name,data['X'],18*62.5e-6)
        np.testing.assert_array_equal(x,a['x']);np.testing.assert_array_equal(F,a['F_exact'])
        got=np.eye(3)+gradient(a['x']-data['X'],data['source'],data['target'],data['g'])
        np.testing.assert_array_equal(got,a['PD2']);assertF(got-F,data['volume'])
def test_07_quadratic_axial(data):quadratic('axial',data)
def test_08_quadratic_shear(data):quadratic('shear',data)
def test_09_XZ_cross_term(data):
    quadratic('cross',data)
    # Optional diagonal basis cannot satisfy the required all-P2 gate on this cloud.
    s=read('C1_PD2_SUMMARY.json')
    assert s['cross']['RK2_DIAG']['Fbar']['max_scaled']>1e-10

def test_10_Ftilde_midpoint_exactness_after_nodal_gates(data):
    summary=read('C1_PD2_SUMMARY.json')
    assert all(summary[f]['PD2']['Fbar']['max_scaled']<=1e-10 for f in ['axial','shear','cross'])
    ids=np.arange(0,len(data['pairs']),683);pp=data['pairs'][ids];X=data['X'];mid=.5*(X[pp[:,0]]+X[pp[:,1]])
    for name in ['axial','shear','cross']:
        with np.load(OUT/(name.upper()+'_NODAL_STATE.npz')) as a,np.load(OUT/(name.upper()+'_BOND_DIAGNOSTICS.npz')) as b:
            Fe=exact(name,mid,18*62.5e-6)[1];got=associated_gradient(X,a['x'],pp,a['PD2'])
            assertF(got-Fe,data['bond_weight'][ids]);assertF(b['PD2_F_error'],data['bond_weight'])
            # Exact quadratic midpoint secant (including the cross term).
            xi=X[pp[:,1]]-X[pp[:,0]];eta=a['x'][pp[:,1]]-a['x'][pp[:,0]]
            assert np.linalg.norm(eta-np.einsum('bij,bj->bi',Fe,xi),axis=1).max()/251.25e-6/.004<1e-10

def test_11_Psi_plus_exact_state_unchanged_split(data):
    summary=read('C1_PD2_SUMMARY.json')
    for name in ['axial','shear','cross']:
        with np.load(OUT/(name.upper()+'_NODAL_STATE.npz')) as a,np.load(OUT/(name.upper()+'_BOND_DIAGNOSTICS.npz')) as b:
            ref=summary[name]['Psi_ref_Pa'];e=psi(a['PD2'])-psi(a['F_exact'])
            for err,w in [(e,data['volume']),(b['PD2_Psi_error'],data['bond_weight'])]:
                m=metrics(err,w,ref);assert m['max_scaled']<=1e-8 and m['weighted_RMS_scaled']<=1e-9

def test_12_installed_Peridynamicsjl_representative_nodes():
    path=OUT/'PERIDYNAMICSJL_CROSSCHECK.json'
    if not path.exists():pytest.skip('Package diagnostic unavailable; see audit report')
    a=read(path.name)
    if a['status']=='NOT_AVAILABLE_WITHOUT_PACKAGE_MODIFICATION':pytest.skip(a['status'])
    assert a['status']=='PASS' and a['force_calls']==0
    for r in a['records']:
        assert np.linalg.norm(np.asarray(r['python_F'])-r['julia_F'])/.004<=1e-10
