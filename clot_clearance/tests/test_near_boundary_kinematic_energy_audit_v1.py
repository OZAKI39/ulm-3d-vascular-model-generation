"""Only the ten requested fixed-state causal checks; no equilibrium reruns."""
import json
from pathlib import Path
import numpy as np
from pd_clot.pfpd_candidate.energy import Material
from pd_clot.pfpd_candidate.kernel import radial
from pd_clot.pfpd_candidate.nonlocal_boundary.nonaffine_manufactured import QuadraticShear,QuadraticAxial
from pd_clot.pfpd_candidate.nonlocal_boundary.kinematic_energy_diagnostics import continuous_ball,psi

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'verification/baqp_near_boundary_kinematic_energy_audit_v1'


def read(n):return json.loads((OUT/n).read_text())
def records():return read('CAUSAL_IDENTITY_SUMMARY.json')
def material():return Material(model='svk',mu=1000.,bulk=666.6666666666666)


def test_01_continuous_ball_odd_moments_and_K():
    d=251.25e-6;r=continuous_ball(d)
    np.testing.assert_array_equal(r['M3_analytic'],np.zeros((3,3,3)))
    assert r['M3_numerical_scaled']<=1e-12 and r['K_relative_error']<=1e-12
    # Independent polynomial antiderivative, not a lattice reference.
    integral=d**5*(1/5-2/6+1/7)
    np.testing.assert_allclose(r['K_analytic'],np.eye(3)*(4*np.pi/3)*integral,rtol=1e-14,atol=0)


def test_02_full_lattice_symmetry_and_quadratic_reproduction():
    r=read('FULL_LATTICE_REFERENCE.json')
    assert r['PASS'] and r['M3_scaled']<=1e-10
    assert all(v['max_scaled']<=1e-10 for v in r['full_reconstruction'].values())
    with np.load(OUT/'FULL_LATTICE_RECONSTRUCTION.npz') as a:
        for name in ['shear','axial']:assert np.linalg.norm(a[name],axis=(1,2)).max()/.004<=1e-10


def test_03_exact_quadratic_Fbar_identity_direct_neighbors():
    with np.load(OUT/'MOMENT_TENSORS.npz') as a:
        X=a['X'];V=a['volume'];K=a['K'];M=a['M3'];inv=a['invK'];f=QuadraticShear()
        for idx in [0,180,3429,6500]:
            xi=X-X[idx];r=np.linalg.norm(xi,axis=1);mask=(r>0)&(r<=251.25e-6);xi=xi[mask];w=radial(r[mask]/251.25e-6)*V[mask]
            directK=sum(weight*np.outer(v,v) for weight,v in zip(w,xi))
            directM=sum(weight*np.einsum('c,d,b->cdb',v,v,v) for weight,v in zip(w,xi))
            np.testing.assert_allclose(directK,K[idx],rtol=1e-13,atol=1e-34)
            np.testing.assert_allclose(directM,M[idx],rtol=1e-12,atol=1e-38)
            eta=f.position(X[mask])-f.position(X[idx])
            recovered=np.einsum('i,ia,ib->ab',w,eta,xi)@inv[idx]
            predicted=np.zeros((3,3));predicted[0]=.001/f.length*(M[idx,2,2,:]@inv[idx])
            assert np.linalg.norm(recovered-f.response(X[idx,2])['F']-predicted)/.004<=1e-10
    for row in records():
        assert row['Fbar']['residual_max_scaled']<=1e-10 and row['Fbar']['residual_RMS_scaled']<=1e-11


def specialized(name,row,field):
    with np.load(OUT/'MOMENT_TENSORS.npz') as a, np.load(OUT/(name.upper()+'_NODE_DIAGNOSTICS.npz')) as d:
        pred=.001/field.length*np.einsum('ib,ibn->in',a['M3'][:,2,2,:],a['invK'])
        assert np.linalg.norm(d['actual_delta_F'][:,row,:]-pred,axis=1).max()/.004<=1e-10
        other=d['actual_delta_F'].copy();other[:,row,:]=0
        assert np.linalg.norm(other,axis=(1,2)).max()/.004<=1e-10


def test_04_specialized_axial_formula():specialized('axial',2,QuadraticAxial())


def test_05_specialized_shear_formula():specialized('shear',0,QuadraticShear())


def test_06_support_and_volume_decomposition():
    with np.load(OUT/'SOURCE_CONTRIBUTIONS.npz') as a:
        assert np.linalg.norm(a['support']+a['volume']-a['total'],axis=1).max()/.004<=1e-10
    s=read('SOURCE_DECOMPOSITION.json')
    assert s['variants']['R1_full']['max_scaled']<=1e-10
    assert s['variants']['R2_support_uniform']['max_scaled']>1e-3
    assert read('HISTORICAL_INTEGRITY.json')['PASS']


def test_07_exact_Ftilde_projection_identity():
    for r in records():
        for key in ['Ftilde','Ftilde_moment_chain']:
            assert r[key]['residual_max_scaled']<=1e-10 and r[key]['residual_RMS_scaled']<=1e-11
    for name in ['SHEAR','AXIAL']:
        with np.load(OUT/(name+'_BOND_DIAGNOSTICS.npz')) as b:
            assert b['projection_residual'].max()/.004<=1e-10
            valid=np.isfinite(b['projection_reduction_ratio'])
            assert b['projection_reduction_ratio'][valid].max()<=1+1e-12
            assert b['projection_reduction_ratio'][valid].min()>=0


def test_08_axial_shear_geometric_row_mapping():
    with np.load(OUT/'AXIAL_NODE_DIAGNOSTICS.npz') as a, np.load(OUT/'SHEAR_NODE_DIAGNOSTICS.npz') as b:
        assert np.linalg.norm(a['actual_delta_F'][:,2,:]-b['actual_delta_F'][:,0,:],axis=1).max()/.004<=1e-10
        assert np.linalg.norm(a['predicted_delta_F'][:,2,:]-b['predicted_delta_F'][:,0,:],axis=1).max()/.004<=1e-10
    assert read('AXIAL_SHEAR_CROSSCHECK_SUMMARY.json')['PASS']


def test_09_nonlinear_Psi_causal_reconstruction():
    for r in records():
        assert r['Psi_Fbar']['residual_max_scaled']<=1e-8 and r['Psi_Ftilde']['residual_max_scaled']<=1e-8
    with np.load(OUT/'BOND_GEOMETRY.npz') as g, np.load(OUT/'MOMENT_TENSORS.npz') as m:
        ids=np.linspace(0,len(g['pairs'])-1,41,dtype=int);i,j=g['pairs'][ids].T;X=m['X'];xi=X[j]-X[i]
        n=xi/np.linalg.norm(xi,axis=1)[:,None];proj=np.eye(3)-n[:,:,None]*n[:,None,:]
        for name,field in [('SHEAR',QuadraticShear()),('AXIAL',QuadraticAxial())]:
            with np.load(OUT/(name+'_NODE_DIAGNOSTICS.npz')) as d, np.load(OUT/(name+'_BOND_DIAGNOSTICS.npz')) as b:
                F=field.response(.5*(X[i,2]+X[j,2]))['F']+(.5*(d['predicted_delta_F'][i]+d['predicted_delta_F'][j]))@proj
                np.testing.assert_allclose(psi(F,material()),b['Psi_predicted'][ids],rtol=1e-12,atol=1e-16)


def test_10_global_energy_decomposition_deterministic():
    r=read('ENERGY_SOURCE_DECOMPOSITION.json')
    for name in ['shear','axial']:
        assert r[name]['deterministic'] and max(r[name]['Gauss32_vs64_relative'])<=1e-12
        with np.load(OUT/(name.upper()+'_BOND_DIAGNOSTICS.npz')) as a:
            for channel,prefix in [('intact','U'),('positive','Uplus')]:
                e=r[name][channel];u0=float(a[prefix+'0_bond_J'].sum());u1=float(a[prefix+'1_bond_J'].sum())
                assert u0==e['U0_J'] and u1==e['U1_J']
                assert float(a[prefix+'0_bond_J'].sum())==u0
                assert abs(e['kinematic_bias_J']+e['weight_quadrature_bias_J']-e['total_bias_J'])/e['U2_J']<=1e-12
