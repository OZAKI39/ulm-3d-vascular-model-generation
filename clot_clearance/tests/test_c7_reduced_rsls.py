import ast
import numpy as np
from scripts import c7_early_branch as a
from scripts import c7_reduced_vi as r


def test_c5_first_positive_root_matches_full_rsls():
    t=a.read(a.OUT/'c5_regression/first_A1/SUMMARY.json');assert t['passed']
    assert t['agreement_errors']['phase_Linf']<=1e-8 and t['agreement_errors']['scaled_state_L2']<=1e-8
    assert a.read(a.OUT/'jacobian/C5_first_A1.json')['passed']


def test_c5_old_failure_root_matches_full_rsls():
    t=a.read(a.OUT/'c5_regression/same_load_as_original/SUMMARY.json');assert t['passed']
    assert t['agreement_errors']['energy_relative']<=1e-9 and t['agreement_errors']['reaction_relative']<=1e-8
    assert a.read(a.OUT/'jacobian/C5_same_load_as_original.json')['passed']


def test_c7_reduced_jacobian_all_directions_and_scales():
    t=a.read(a.OUT/'jacobian/C7_FIRST_WEAK.json');assert t['passed'] and len(t['rows'])==30
    assert t['max_relative_error']<=5e-5
    assert t['shape'][0]==735+t['pool_size'] and t['pool_size']<.25*46059


def test_solver_selection_has_no_oracle_bond_ids():
    # IDs from earlier C5/C7 activation outputs must never be solver constants.
    for path in [a.ROOT/'scripts/c7_early_branch.py',a.ROOT/'scripts/c7_reduced_vi.py']:
        tree=ast.parse(path.read_text())
        assert not any(isinstance(n,ast.Constant) and n.value in (4488,5059) for n in ast.walk(tree))
    branch=a.read(a.OUT/'branch/ALL.json');first=branch['states'][0]
    assert first['no_oracle_initialization']
    initial=a.arrays(a.OUT/first['label']/'INITIAL_STATE.npz')
    assert np.count_nonzero(initial['lower'])==0 and np.count_nonzero(initial['s'])==0


def test_real_regressions_triggered_full_kkt_pool_expansion():
    for s in a.read(a.OUT/'c5_regression/ALL.json')['states']:
        assert s['initial_pool_size']==0 and s['expansion_count']==1
        assert len(s['pool_history'][0]['violated_excluded_IDs'])>0
        assert s['excluded_violation_max']<=1e-8 and s['full_KKT_bond_count']==5392


def test_checkpoint_resume_returns_identical_state_without_solving(monkeypatch):
    x,s,e=a.prior.initial();target=a.prior.targets()[0];p=a.problem(5,target['strain'],s);x=a.old.tracker.displacement_predictor(p,x,e)
    previous=a.hashlib.sha256(x.tobytes()+s.tobytes()).hexdigest()
    def forbidden(*args,**kwargs):raise AssertionError('Resume recomputed the solution')
    monkeypatch.setattr(r,'solve_once',forbidden)
    row,data=a.solve_pool(p,x,s,'c5_regression/first_A1',None,previous)
    stored=a.arrays(a.OUT/'c5_regression/first_A1/STATE.npz')
    assert row['scientifically_accepted'] and np.array_equal(data['x'],stored['x']) and np.array_equal(data['s'],stored['s'])
    altered=x.copy();altered[p.free,0]+=1e-12
    assert a.checkpoint_key(p,x,s,previous,'coupled')!=a.checkpoint_key(p,altered,s,previous,'coupled')
