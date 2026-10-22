from common import *
import subprocess,datetime
old=read(IDENTITY/'PRODUCTION_SOURCE_IDENTITY.json');sources=read(IDENTITY/'EXACT_DEPLOYED_SOURCE_MANIFEST.json')
assert all(sha(RETAINED/n)==h for n,h in sources.items());assert sha(old['actual_binary'])==old['actual_binary_sha256'];assert sha(old['PETSc']['path'])==old['PETSc']['sha256']
inputs={}
for d in [TARGET,PREVIOUS]:
 for n,h in read(d/'OUTPUT_SHA256.json').items():assert sha(d/n)==h;inputs[str((d/n).relative_to(BASE))]=h
for n in ['production_workflow/GPU_MANUFACTURED_RUN_IDENTITY.json','production_workflow/FIELD_VALIDATION.json']:
 inputs[n]=sha(BASE/n)
for n in ['run/solver.xml','run/PETSC_OPTIONS.txt','run/solver.log','run/8-procs/stFile_071.bin','run/8-procs/result_071.vtu']:
 inputs[str((P1/n).relative_to(BASE))]=sha(P1/n)
inputs[str((TARGET/'runtime/PIPE_CURVED_ELEVATION.npz').relative_to(BASE))]=sha(TARGET/'runtime/PIPE_CURVED_ELEVATION.npz')
assert inputs[str((TARGET/'runtime/PIPE_CURVED_ELEVATION.npz').relative_to(BASE))]=='dbad7f2be2e4eb4a718f2db275b6ab3bf663915f229088367530db705407d48e'
inputs[str((SOURCE/'sv13q_reuse.h').relative_to(BASE))]=sha(SOURCE/'sv13q_reuse.h')
put('INPUT_SHA256.json',inputs)
res={k:(Path('/sys/fs/cgroup')/k).read_text() for k in ['cpu.max','cpuset.cpus.effective','memory.max','memory.current','memory.stat','memory.events']}
for k,cmd in {'GPU':['nvidia-smi','--query-gpu=name,driver_version,memory.total,memory.used,memory.free','--format=csv'],'other_GPU_processes':['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv'],'disk':['df','-h',str(BASE)],'processes':['ps','-eo','pid,ppid,pcpu,rss,comm','--sort=-pcpu']}.items():res[k]=subprocess.check_output(cmd,text=True)
for k,cmd in {'HEAD':['rev-parse','HEAD'],'branch':['branch','--show-current'],'status':['status','--short']}.items():res[k]=subprocess.check_output(['git','-c',f'safe.directory={BASE}','-C',str(BASE),*cmd],text=True)
res.update(UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),**managed().cgroup_memory());put('RESOURCES.json',res)
install=Path(old['PETSc']['path']).parent.parent;conf=(install/'include/petscconf.h').read_text();(HERE/'PETSC_BUILD_CONF.txt').write_text(conf)
psrc=install.parent/'source_gpu13_configfix';paths=['src/mat/impls/aij/seq/aijfact.c','src/ksp/pc/impls/asm/asm.c','src/ksp/pc/impls/fieldsplit/fieldsplit.c','src/ksp/ksp/utils/schurm/schurm.c']
evidence={}
for n in paths:
 p=psrc/n;t=p.read_text();evidence[n]=dict(sha256=sha(p),path=str(p));(HERE/(Path(n).name+'.source.txt')).write_text('\n'.join(f'{i+1}: {line}' for i,line in enumerate(t.splitlines()) if (n.endswith('aijfact.c') and i<180) or (n.endswith('asm.c') and any(w in line for w in ['MatIncreaseOverlap','MatCreateSubMatrices','ISCreate','PCSetUp_ASM'])) or (n.endswith('fieldsplit.c') and any(w in line for w in ['selfp','SELFP','MatSchur','MatCreateSubMatrix'])) or (n.endswith('schurm.c') and any(w in line for w in ['DIAG','SELFP','MatMat','MatDiagonal','VecReciprocal']))))
put('PETSC_BUILD_AUDIT.json',dict(version='3.25.5',library=old['PETSc'],config_sha256=sha(install/'include/petscconf.h'),int_bytes=4,scalar_bytes=8,HYPRE=False,GAMG=True,external_LU=False,packages='blaslapack:cuda:mathlib:mpi:pthread',source=evidence,LU_symbolic_allocates_numeric_arrays=True,LU_symbolic_evidence='aijfact.c MatLUFactorSymbolic_SeqAIJ PetscShmgetAllocateArray(bdiag[0]+1,sizeof(PetscScalar),&b->a)',production_source_files_verified=len(sources)))
print('Frozen inputs, source and PETSc library verified',len(inputs),flush=True)

cuda=read(TARGET/'CUDA_LIBRARY_IDENTITY.json');verified={}
def cuda_walk(v):
 if isinstance(v,dict):
  if 'path' in v and 'sha256' in v:
   assert sha(v['path'])==v['sha256'];verified[v['path']]=v['sha256']
  for path,digest in v.get('current_libraries_sha256',{}).items():
   assert sha(path)==digest;verified[path]=digest
  for x in v.values():cuda_walk(x)
 elif isinstance(v,list):
  for x in v:cuda_walk(x)
cuda_walk(cuda)
put('CURRENT_SOURCE_LIBRARY_IDENTITY.json',dict(production_binary=old['actual_binary'],production_binary_sha256=old['actual_binary_sha256'],production_source_commit='c3f0bb892b765b718f61069ecd9726dbc6d177fd',source_manifest_sha256=sha(IDENTITY/'EXACT_DEPLOYED_SOURCE_MANIFEST.json'),verified_source_files=len(sources),PETSc=old['PETSc'],CUDA_libraries_verified=verified,CUDA_toolkit=subprocess.check_output(['/usr/local/cuda/bin/nvcc','--version'],text=True),verified_UTC=res['UTC']))
