from sv13_support import artifact,ROOT
from sv_validation.sv12 import checkpoint_audit
def test_native_automatic_stop_with_complete_checkpoint():
    v=artifact('cpu_validation');e=artifact('cpu_execution');ref=artifact('reference_freeze')
    assert v['status']=='PASS' and v['native_exit']==0 and v['checkpoint_complete']
    assert e['linear_failures']==e['nonlinear_failures']==e['ill_conditioned_warnings']==0
    request=e['stop_request'];assert request['kind']=='STEADY'
    assert e['last_step']==request['requested_final_step']
    assert e['steps_executed']==e['last_step']-e['initial_step']
    cp=artifact('cpu_final_checkpoint')
    actual=checkpoint_audit(cp['path'],e['last_step'],ref['dt'])
    assert actual['sha256']==cp['sha256']
    assert (ROOT/'outputs/sv1_3/cpu_early_stop/STOP_SIM').read_text().strip()==str(e['last_step'])
def test_cpu_configuration_is_baseline_identical():
    e=artifact('cpu_execution');ref=artifact('reference_freeze')
    assert e['PETSC_OPTIONS']==ref['PETSc']['PETSC_OPTIONS']
    assert e['mpi_ranks']==ref['PETSc']['mpi_ranks']
    assert (ROOT/'configs/sv1_2/sv_flow.xml').read_bytes()==(ROOT/'outputs/sv1_3/cpu_early_stop/solver.xml').read_bytes()
