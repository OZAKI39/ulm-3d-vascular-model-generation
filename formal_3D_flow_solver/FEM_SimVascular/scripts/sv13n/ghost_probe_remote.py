"""Run deterministic forward/reverse/local form tests, preserving every run."""
import sys,re
from runner_remote import *
phase=sys.argv[1];assert phase in ('old','gpu13','gpu123')
w=json.loads((BASE/'configs/baseline_L_compatibility_winner.json').read_text()) if phase=='old' else load('petsc_'+phase+'_build')
prefix=Path(w['prefix']);mpi=Path(w['MPI_prefix']);wrapper=w['candidate_wrapper']
launcher=json.loads((BASE/'configs/baseline_L_mpi_application_gate.json').read_text())['working_launcher']
binary=B/('petsc_cuda_ghost_probe_'+phase)
extra={'LD_LIBRARY_PATH':str(prefix/'lib'),'PETSC_OPTIONS':'-skip_petscrc -use_gpu_aware_mpi 0'}
c=run([mpi/'bin/mpicc',B/'petsc_cuda_ghost_probe.c','-I'+str(prefix/'include'),'-L'+str(prefix/'lib'),'-Wl,-rpath,'+str(prefix/'lib'),'-lpetsc','-lm','-o',binary],phase+'_ghost_compile',cuda=wrapper,extra_env=extra)
assert okay(c)
runs=[]
for variant,n,t in [('CPU_seq',1,'standard'),('CUDA_seq',1,'cuda'),('CPU_MPI',2,'standard'),('CUDA_MPI',2,'cuda')]:
 for i in range(1,2):
  r=run([launcher,'-n',str(n),wrapper,binary,'-vec_type',t],phase+'_'+variant+'_'+str(i),timeout=60,cuda=wrapper,extra_env=extra)
  txt=text(r)
  r.update(variant=variant,repetition=i,ranks=n,types=re.findall(r'TYPE rank=(\d+) type=(\S+) base=(\S+) owned=(\d+) ghost_blocks=(\d+) ghost_index=(-?\d+)',txt),locals=re.findall(r'LOCAL rank=(\d+) type=(\S+) size=(-?\d+) expected=(\d+) present=(\d+)',txt),values=re.findall(r'VALUE phase=(\w+) rank=(\d+) local=(\d+) actual=(\S+) expected=(\S+)',txt),checks=re.findall(r'CHECK rank=(\d+) forward_error=(\S+) reverse_error=(\S+) local_size=(\d+) expected_size=(\d+) errors=(\d+)',txt),hard_error='PETSC ERROR' in txt,error_stack=[line for line in txt.splitlines() if 'PETSC ERROR' in line])
  r['accepted']=okay(r) and not r['hard_error'] and len(r['checks'])==n and 'GHOST_PROBE PASS' in txt and all(float(a)==float(e) for _,_,_,a,e in r['values']) and all(float(f)==float(b)==0 and ls==es and err=='0' for _,f,b,ls,es,err in r['checks'])
  runs.append(r)
  write('ghost_probe_'+phase,dict(status='RUNNING',runs=runs,compile=c))
write('ghost_probe_'+phase,dict(status='PASS' if all(r['accepted'] for r in runs) else 'FAIL',runs=runs,compile=c,source_sha256=digest(B/'petsc_cuda_ghost_probe.c'),sequential_note='No remote ghosts on one rank, as in actual solver. Nonzero ghost entries and reverse owner contributions are checked on two ranks.'))
