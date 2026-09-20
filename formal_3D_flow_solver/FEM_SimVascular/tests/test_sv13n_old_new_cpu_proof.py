from sv13n_support import *
def test_old_cpu_safely_stopped_by_user():
 d=accepted('OLD_PETSC_CPU_PROOF_20_acceptance');assert d['fixed20_requirement_superseded']
 assert d['exit_code']==0 and d['steps_completed']<20 and d['reload_pass']
def test_new_cpu_proof_deferred():accepted('NEW_PETSC_CPU_PROOF_20_acceptance')
