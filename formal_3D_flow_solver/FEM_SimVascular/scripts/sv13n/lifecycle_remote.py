"""Read PETSc lifecycle state at the application's MPI_Finalize, without mutation."""
import re,shutil,sys
from runner_remote import *
phase=sys.argv[1];assert phase in ('before','after')
sv=load('svmp_cpu_build');w=load('petsc_cpu_build')
mpi=json.loads((BASE/'configs/baseline_L_mpi_application_gate.json').read_text())
case=BASE/'outputs'/('OFFICIAL_CPU_LIFECYCLE_'+phase.upper())
assert not case.exists();shutil.copytree(BASE/'outputs/official_fluid_gpu_smoke',case)
options=json.loads((BASE/'configs/petsc_options.json').read_text())['PETSC_OPTIONS']+' -use_gpu_aware_mpi 0 -mat_type aij -vec_type standard'
commands=case/'lifecycle.gdb'
commands.write_text('set pagination off\nset confirm off\nset breakpoint pending on\nset $petsc_finalize_calls=0\nbreak PetscFinalize\ncommands\nsilent\nset $petsc_finalize_calls=$petsc_finalize_calls+1\ncontinue\nend\nbreak MPI_Finalize\nrun\nprintf "SV13N_PETSC_STATE initialized=%d finalized=%d count=%d\\n", *(unsigned char*)&PetscInitializeCalled, *(unsigned char*)&PetscFinalizeCalled, $petsc_finalize_calls\nbt\ncontinue\n')
record=run([mpi['working_launcher'],'-n','1',w['candidate_wrapper'],'gdb','--batch','-x',commands,'--args',sv['executable'],'solver.xml'],
 'lifecycle_'+phase,cwd=case,timeout=180,cuda=w['candidate_wrapper'],extra_env={'LD_LIBRARY_PATH':sv['runtime_library_path'],'PETSC_OPTIONS':options})
raw=text(record);m=re.search(r'SV13N_PETSC_STATE initialized=(\d+) finalized=(\d+) count=(\d+)',raw)
observed=tuple(map(int,m.groups())) if m else None
expected=(1,0,0) if phase=='before' else (0,1,1)
d=dict(status='PASS' if okay(record) and observed==expected and 'exited normally' in raw else 'FAIL',phase=phase,observation=observed,expected=expected,
 source='Read-only GDB inspection of exported PETSc state immediately before application MPI_Finalize',execution=record,solver_sha256=sv['executable_sha256'])
write('lifecycle_'+phase+'_gate',d)
assert d['status']=='PASS',d
