import numpy as np
from scripts import c7_early_branch as a
from scripts import c7_reduced_vi as r


def test_expanded_state_keeps_excluded_exactly_at_lower():
    p=a.problem(5,.0045);ids=np.array([0,7,20]);rp=r.ReducedProblem(p,ids);u=np.zeros(p.nx+3);u[-3:]=[.001,.002,.003];full=rp.expand(u)
    excluded=np.ones(p.ns,bool);excluded[ids]=False
    assert np.array_equal(full[p.nx:][excluded],p.lower[excluded])
    assert np.array_equal(rp.compress(full),u)


def test_pool_contains_all_weak_and_non_strong_lower_variables(monkeypatch):
    p=a.problem(5,.0045);s=np.zeros(p.ns);qs=np.ones(p.ns);qs[[4,15]]=[0,-1e-3]
    monkeypatch.setattr(r.old,'q_fields',lambda *args:(None,None,None,qs))
    ids=r.initial_pool(p,p.family.X,s)
    assert {4,15}<=set(ids)
    excluded=np.setdiff1d(np.arange(p.ns),ids)
    assert np.all(qs[excluded]>a.cfg()['q_tol'])


def test_empty_predictor_pool_is_integer_and_expands_on_full_kkt(monkeypatch):
    p=a.problem(5,.0045);s=np.zeros(p.ns);qs=np.ones(p.ns)
    monkeypatch.setattr(r.old,'q_fields',lambda *args:(None,None,None,qs))
    ids=r.initial_pool(p,p.family.X,s);assert ids.dtype.kind in 'iu' and len(ids)==0
    qs[30]=-2e-8;z=p.pack(p.prescribe(p.family.X),s)
    bad=r.excluded_violations(p,z,ids);assert np.array_equal(bad,[30])
    expanded=r.expand_pool(p,ids,bad);assert 30 in expanded
    assert set(r.neighbors(p.family,[30]))<=set(expanded)


def test_first_ring_is_not_recursive_closure():
    p=a.problem(5,.0045);ids=r.neighbors(p.family,[0]);nodes=np.unique(p.family.pairs[[0]])
    assert np.all(np.any(np.isin(p.family.pairs[ids],nodes),axis=1))
    assert len(ids)<p.ns


def test_positive_history_stays_in_pool(monkeypatch):
    p=a.problem(5,.0045);s=np.zeros(p.ns);s[31]=.02;qs=np.ones(p.ns)
    monkeypatch.setattr(r.old,'q_fields',lambda *args:(None,None,None,qs))
    assert 31 in r.initial_pool(p,p.family.X,s)
