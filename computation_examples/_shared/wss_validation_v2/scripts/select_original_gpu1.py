"""Select original H0 GPU1 backend BEFORE any new-mesh CFD is started.

Preserve all previous unexecuted input variants and hashes. No mesh/BC changes.
"""
import shutil,time
from case_common import *
for name in ['vessel_medium','vessel_fine']:
    case=V/'stage3'/name
    assert not (case/'run/solver.log').exists()
    archive=case/'reports/unexecuted_cpu8_inputs'
    assert not archive.exists()
    archive.mkdir(parents=True)
    for rel in ['policy.json','input_hashes.json','run/PETSC_OPTIONS.txt']:
        target=archive/rel;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(case/rel,target)
    before=json.loads((case/'input_hashes.json').read_text())
    options=case/'run/PETSC_OPTIONS.txt';s=options.read_text()
    assert '-mat_type aij -vec_type standard' in s
    options.write_text(s.replace('-mat_type aij -vec_type standard','-mat_type aijcusparse -vec_type cuda'))
    policy=json.loads((case/'policy.json').read_text())
    policy.update(MPI_ranks=1,linear_algebra_backend='Original H0 GPU PETSc aijcusparse/cuda; one MPI rank; ASM overlap2/ILU2',execution_selection_reason='Same as qualified original H0; current actual GPU1 time-step experiment converges linear solves, GPU8 rejected; unexecuted CPU8 inputs archived')
    dump(case/'policy.json',policy);lock_case(case)
    after=json.loads((case/'input_hashes.json').read_text())
    changes={k:dict(before=before[k],after=after[k]) for k in before if before[k]!=after[k]}
    assert set(changes)=={'policy.json','run/PETSC_OPTIONS.txt'}
    dump(case/'reports/backend_selection_before_solve.json',dict(unix=time.time(),same_mesh_and_XML=True,changes=changes,previous_inputs_executed=False,archive=str(archive),numerical_and_steady_tolerances_unchanged=True))
    print(name,changes)
