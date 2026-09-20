from sv13n_support import *
def test_real_solver_finalize_once_before_mpi():lifecycle_gate(accepted('solver_finalize_lifecycle'))
def test_minimal_cuda_object_cleanup():
 d=json.loads((R/'remote/minimal_lifecycle.json').read_text());assert d['status']=='PASS' and d['execution']['exit_code']==0
@pytest.mark.parametrize('state',[[1,0,0],[0,1,2],[1,0,1]])
def test_wrong_order_or_double_finalize_rejected(state):
 with pytest.raises(GateError):lifecycle_gate(dict(normal_exit=True,observed_at_MPI_Finalize=state))
