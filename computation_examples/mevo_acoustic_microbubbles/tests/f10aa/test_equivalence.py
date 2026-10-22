import numpy as np,pytest,h5py
from scipy import linalg
from mevo_mb.fast_sim import f10g_contact as ref
from mevo_mb.fast_sim import accelerated_contact as opt
from mevo_mb.fast_sim.accelerated_sequential import sequential_step
from scripts.straight_tube_exp1.f10aa.common import BASE,ROOT,R,policy

@pytest.mark.parametrize('n',[32,64,128])
def test_all_pair_products_and_adjoint(n):
    x=np.random.default_rng(n).normal(size=(n,3))*1e-5
    dense,p,g=ref.contact_geometry(x,R);s,pp,gg=opt.contact_geometry(x,R)
    v=np.random.default_rng(8).normal(size=3*n);lam=np.random.default_rng(9).normal(size=len(p))
    np.testing.assert_array_equal(p,pp);np.testing.assert_array_equal(g,gg)
    np.testing.assert_allclose(s@v,dense@v,atol=2e-15,rtol=1e-13)
    np.testing.assert_allclose(s.T@lam,dense.T@lam,atol=2e-14,rtol=1e-12)
    assert abs(np.dot(s@v,lam)-np.dot(v,s.T@lam))<1e-11
    assert s.nbytes<dense.nbytes/10

def test_sqp_algorithm_and_gates_are_literal_reference():
    old=(ROOT/'src/mevo_mb/fast_sim/sequential_hard_sphere.py').read_text()
    old=old.replace('from .f10g_contact import contact_geometry\nfrom .hard_sphere import active_set_qp','from .accelerated_contact import contact_geometry, active_set_qp, factorize')
    old=old.replace('def sequential_step(x, force, radius, zeta, dt, A, policy):','def sequential_step(x, force, radius, zeta, dt, A, policy, *, backend="cpu"):')
    old=old.replace('factor = linalg.cho_factor(A, lower=True)','factor = factorize(A, backend)').replace('linalg.cho_solve(factor, ','factor.solve(')
    old=old.replace('Opt-in bounded true-sphere SQP, fixed original F/A and physical time.','F10AA representation-only SQP, copied from frozen F10Q with identical gates.')
    assert (ROOT/'src/mevo_mb/fast_sim/accelerated_sequential.py').read_text()==old

def test_frozen_fixtures_against_original_dense_kkt():
    p=policy();largest=0
    with h5py.File(BASE/'baseline/REFERENCE.h5') as h:
      for name,q in h.items():
        x,F,A=q['x'][:],q['force_N'][:],q['A'][:];z=q.attrs['zeta'];dt=1e-6
        assembled=opt.resistance(x,R)['A']
        np.testing.assert_allclose(assembled,A,rtol=1e-13,atol=5e-13)
        r=sequential_step(x,F,R,z,dt,assembled,p);assert r['converged']==q['reference'].attrs['converged']
        assert all(r['gates'].values());largest=max(largest,r['active_count'])
        By,_,gap=ref.contact_geometry(r['positions_m'],R)
        lam=r['pair_force_N'];velocity=r['velocity_m_s'].ravel()
        balance=(z*A@velocity-F.ravel()-By.T@lam).reshape(x.shape)
        scale=max(np.linalg.norm(F,axis=1).max(),abs(lam).max(initial=0),1e-18)
        # Independent dense evaluation uses original A/B and physical units.
        assert np.linalg.norm(balance,axis=1).max()<=p['force_residual_atol_N']+p['force_residual_rtol']*scale
        assert gap.min()>=-p['hard_geometry_m']
        assert abs(lam*gap).max(initial=0)/scale<=p['actual_geometry_complementarity_m']
        assert lam.min(initial=0)>=-p['dual_negative_N']
        assert np.linalg.norm(r['positions_m']-q['reference/positions_m'][:],axis=1).max()<1e-14
    assert largest>=500

def test_illegal_inputs_fail_closed():
    x=np.array([[0,0,0],[R,0,0]])
    with pytest.raises(ValueError):sequential_step(x,np.zeros_like(x),R,1.,1e-6,np.eye(6),policy())
    x[1,0]=3*R
    with pytest.raises(np.linalg.LinAlgError):sequential_step(x,np.zeros_like(x),R,1.,1e-6,-np.eye(6),policy())
    p=policy();p['max_sequential_qp']=0
    with pytest.raises(ValueError):sequential_step(x,np.zeros_like(x),R,1.,1e-6,np.eye(6),p)

def test_degenerate_inconsistent_constraints_are_rejected():
    # Engineered infeasible algebraic QP, not a physical bubble trajectory.
    i=np.array([0,0]);j=np.array([1,1]);normals=np.array([[1.,0,0],[-1.,0,0]])
    B=opt.PairMatrix(2,i,j,normals);A=np.eye(6);factor=opt.factorize(A)
    _,_,status,_=opt.active_set_qp(A,factor,np.zeros(6),B,np.array([-1.,-1.]))
    assert status=='active_set_inconsistent_equalities'
