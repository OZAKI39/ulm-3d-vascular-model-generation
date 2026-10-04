import ast
from scripts import petsc_vi_audit as audit
from pd_clot.pfpd_candidate import petsc_vi_candidate as vi


def test_sparse_jacobian_all_states_directions_steps():
    g=audit.read(audit.OUT/'jacobian/GATE.json')
    assert g['passed'] and len(g['states'])==3
    for s in g['states']:
        assert len(s['checks'])==30 and len({r['direction'] for r in s['checks']})==10
        assert len({r['epsilon'] for r in s['checks']})==3
        assert s['max_error']<=5e-5


def test_phase_diagonal_sparse_storage():
    s=audit.read(audit.OUT/'JACOBIAN_SPARSITY.json')
    assert s['shape']==[5617,5617] and s['block_nnz']['Jss']==5392
    assert s['structural_nnz']==225**2+2*225*5392+5392
    assert s['CSR_bytes']<40_000_000


def test_derivative_source_reuse_has_only_declared_changes():
    old,loop=vi._source_nodes();new=vi.adapted_constructor().body[0]
    old_body=[ast.unparse(n) for n in old.body]
    new_body=[ast.unparse(n) for n in new.body]
    removed=[s for s in old_body if s not in new_body]
    added=[s for s in new_body if s not in old_body]
    assert len(removed)==3 and len(added)==2
    assert 'self.C[self.clipped] = 0' in removed
    assert 'self.diag = self.Qs / self.D' in added
