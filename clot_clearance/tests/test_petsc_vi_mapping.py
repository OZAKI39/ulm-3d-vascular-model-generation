"""VI semantics and mapping regressions, with actual PETSc evidence."""
import ast
import numpy as np
from scripts import petsc_vi_audit as audit
from pd_clot.pfpd_candidate import petsc_vi_candidate as vi


def test_mechanical_bounds_unbounded():
    p=audit.Problem(5,.0045);lo,hi=vi.bounds(p)
    assert np.isneginf(lo[:p.nx]).all() and np.isposinf(hi[:p.nx]).all()


def test_phase_bounds_use_history():
    p=audit.Problem(5,.0045);p.lower=np.linspace(0,.1,p.ns);lo,hi=vi.bounds(p)
    assert np.array_equal(lo[p.nx:],p.lower) and np.all(hi[p.nx:]==p.sc)


def test_lower_semantics_actual_petsc():
    r=audit.read(audit.OUT/'mapping/SCALAR_AND_VECTOR.json')['cases'][0]
    assert r['passed'] and r['value']==0 and r['q']>0 and r['reason']>0


def test_upper_semantics_actual_petsc():
    r=audit.read(audit.OUT/'mapping/SCALAR_AND_VECTOR.json')['cases'][1]
    assert r['passed'] and r['value']==1 and r['q']<0 and r['reason']>0


def test_vector_roundtrip_actual_petsc():
    assert audit.read(audit.OUT/'mapping/SCALAR_AND_VECTOR.json')['roundtrip']


def test_residual_is_raw_q_not_projected():
    x,s,eps=audit.initial();p=audit.Problem(5,eps);z=p.pack(x,s);r=vi.residual(p,z)
    _,_,_,qs=audit.prev.q_fields(p,x,s)
    np.testing.assert_array_equal(r[p.nx:],qs)
    assert np.max(abs(qs))>.1 and p.metrics(x,s)['projected_phase_KKT']==0


def test_diagnostic_extension_agrees_inside_domain_and_is_not_default():
    x,s,eps=audit.initial();p=audit.Problem(5,eps);z=p.pack(x,s)
    np.testing.assert_allclose(vi.residual(p,z),vi.residual(p,z,diagnostic_extension=True),rtol=1e-13,atol=1e-14)
    z[p.nx]=-.001
    import pytest
    with pytest.raises(ValueError,match='Early phase'):vi.residual(p,z)


def test_acceptance_does_not_trust_reason(monkeypatch):
    p=audit.Problem(5,.0045);z=p.pack(audit.prev.v2.affine(p),np.zeros(p.ns))
    monkeypatch.setattr(audit,'full_acceptance',lambda *a:dict(valid=False))
    assert not audit.acceptance(p,z,4)['scientifically_accepted']
    monkeypatch.setattr(audit,'full_acceptance',lambda *a:dict(valid=True))
    assert audit.acceptance(p,z,-5)['scientifically_accepted']
