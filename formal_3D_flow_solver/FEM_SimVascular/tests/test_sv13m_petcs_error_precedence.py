from sv13m_support import *
from sv_validation.sv13m import *
def test_later_ksp_convergence_never_overrides_petsc_error():
 d=good_flow();d.update(PETSc_error_detected=True,petsc_reason='CONVERGED_RTOL')
 with pytest.raises(GateError,match='PETSC_HARD_ERROR'):flow_gate(d)
def test_original_failure_contains_misleading_convergence():
 d=read('baseline_ghost_failure');assert d['history']['petsc_reasons'] and d['status']=='PASS' and d['exit_code']!=0

def test_post_finalize_mpi_error_rejected_even_with_fields():
 d=good_flow();d['MPI_error_detected']=True
 with pytest.raises(GateError):flow_gate(d)

def test_ksp_finalize_reproducer_detects_exit_failure():
 rows=read('finalize_ksp_probe')['runs']
 assert rows[0]['proper_finalize'] is False and rows[0]['exit_code']!=0
 assert rows[1]['proper_finalize'] is True and rows[1]['exit_code']==0
