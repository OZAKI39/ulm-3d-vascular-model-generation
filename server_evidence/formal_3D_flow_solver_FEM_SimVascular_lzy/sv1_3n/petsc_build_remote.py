"""Clean official 3.25.5 build; CPU first, no patches to PETSc."""
import sys,tarfile,hashlib,shutil,re,urllib.request
from runner_remote import *
candidate=sys.argv[1];assert candidate in ('cpu','gpu13','gpu123')
attempt=sys.argv[2] if len(sys.argv)>2 else 'initial';assert attempt in ('initial','configfix')
run_id=candidate if attempt=='initial' else candidate+'_'+attempt
if candidate!='cpu':
 strategy=json.loads((BASE/'configs/convergence_strategy.json').read_text())
 assert strategy['strategy']=='GPU_DEVELOPMENT_ONLY','Explicit GPU-first user strategy required'
manifest=json.loads((BASE/'configs/petsc325_source.json').read_text())
old=json.loads((BASE/'configs/baseline_L_compatibility_winner.json').read_text())
mpi=json.loads((BASE/'configs/baseline_L_mpi_application_gate.json').read_text())
E=BASE/'external/petsc325';E.mkdir(parents=True,exist_ok=True)
archive=E/'petsc-3.25.5.tar.gz'
if not archive.exists():urllib.request.urlretrieve(manifest['archive_url'],archive)
assert digest(archive)==manifest['source_archive_sha256']
source=E/('source_'+run_id);prefix=E/('install_'+candidate);arch='arch-sv13n-'+candidate
assert not source.exists() and not prefix.exists()
original={}
with tarfile.open(archive) as t:
 for m in t:
  if m.isfile():original[str(Path(m.name).relative_to('petsc-3.25.5'))]=hashlib.sha256(t.extractfile(m).read()).hexdigest()
 temp=E/('extract_'+run_id);temp.mkdir();t.extractall(temp,filter='data')
(temp/'petsc-3.25.5').rename(source)
assert all(digest(source/n)==h for n,h in original.items())
write('petsc_'+run_id+'_source_integrity',dict(status='PASS',phase='before_configure',files=original,modifications=[],source=str(source)))
mpiprefix=Path(old['MPI_prefix']);launcher=mpi['working_launcher']
cuda_path='/usr/local/cuda' if candidate=='gpu13' else old['CUDA_prefix'] if candidate=='gpu123' else ''
wrapper=BASE/('use_'+candidate+'_env.sh')
paths=[str(mpiprefix/'bin'),'/usr/bin','/bin']
if cuda_path:paths.insert(0,cuda_path+'/bin')
if candidate=='gpu123':paths.insert(0,str(BASE.parent/'sv1_3j/external/compat_cuda/host-bin'))
cc='/usr/bin/gcc' if candidate=='gpu13' else '/usr/bin/gcc-12'
cxx='/usr/bin/g++' if candidate=='gpu13' else '/usr/bin/g++-12'
wrapper.write_text('#!/bin/sh\nexport PATH='+':'.join(paths)+'\nexport OMPI_CC='+cc+' OMPI_CXX='+cxx+' OMPI_FC=/usr/bin/gfortran-12\nexport LD_LIBRARY_PATH='+str(mpiprefix/'lib')+(':'+cuda_path+'/lib64' if cuda_path else '')+':${LD_LIBRARY_PATH:-}\nexec "$@"\n');wrapper.chmod(0o755)
extra={'LD_LIBRARY_PATH':str(prefix/'lib'),'PETSC_OPTIONS':'-skip_petscrc -use_gpu_aware_mpi 0'}
cmd=['/usr/bin/python3','-B','configure','--prefix='+str(prefix),'PETSC_ARCH='+arch,
 '--with-cc='+str(mpiprefix/'bin/mpicc'),'--with-cxx='+str(mpiprefix/'bin/mpicxx'),'--with-mpiexec='+launcher,
 '--with-fc=0','--with-fortran-bindings=0','--with-debugging=0','--with-shared-libraries=1','--with-scalar-type=real','--with-precision=double','--with-64-bit-indices=0',
 next(a for a in old['configure_command'] if a.startswith('--with-blaslapack-lib=')),
 '--with-x=0','--with-cxx-dialect=C++17','--with-hip=0','--with-opencl=0','--with-cuda='+('1' if cuda_path else '0')]
if cuda_path:
 # Verified in this exact archive, not copied on trust from 3.19 flags.
 interface=(source/'config/BuildSystem/config/packages/CUDA.py').read_text()
 assert 'with-cuda-arch' in interface
 gpu_probe=next(p for p in load('pre_install_environment')['probes'] if '--query-gpu=name,driver_version,memory.total,compute_cap' in p['command'])
 assert gpu_probe['exit_code']==0
 arch_observed=gpu_probe['stdout'].strip().split(',')[-1].strip().replace('.','')
 assert arch_observed=='89' and '4090' in gpu_probe['stdout']
 cmd+=['--with-cuda-dir='+cuda_path,'--with-cudac='+cuda_path+'/bin/nvcc','--with-cuda-arch='+arch_observed,'--with-cuda-dialect=C++17','CUDAFLAGS=-ccbin '+cxx]
d=dict(status='BUILDING',candidate=candidate,version='3.25.5',commit=manifest['commit'],source=str(source),prefix=str(prefix),PETSC_ARCH=arch,MPI_prefix=str(mpiprefix),launcher=launcher,candidate_wrapper=str(wrapper),CUDA_prefix=cuda_path,CUDA_version='13.2' if candidate=='gpu13' else '12.3.2' if candidate=='gpu123' else 'disabled',source_archive_sha256=digest(archive),configure_command=cmd,steps=[])
def execute(args,name,timeout=3600):
 r=run(args,run_id+'_'+name,cwd=source,timeout=timeout,cuda=wrapper,extra_env=extra);d['steps'].append(r);write('petsc_'+candidate+'_build',d);write('petsc_'+run_id+'_'+attempt+'_build',d)
 if not okay(r):d.update(status='FAIL',failed_step=name);write('petsc_'+candidate+'_build',d);write('petsc_'+run_id+'_'+attempt+'_build',d);raise SystemExit(1)
 return r
d.update(CC=cc,CXX=cxx,host_compiler_selection='System default candidate for CUDA13; acceptance requires actual NVCC and PETSc configure tests. No allow-unsupported-compiler flag.')
execute([cc,'--version'],'compiler_c',30);execute([cxx,'--version'],'compiler_cxx',30)
if cuda_path:execute([cuda_path+'/bin/nvcc','--version'],'compiler_cuda',30)
cfg=execute(cmd,'configure');make=execute(['make','-j8','V=1','PETSC_DIR='+str(source),'PETSC_ARCH='+arch,'all'],'make')
install=execute(['make','PETSC_DIR='+str(source),'PETSC_ARCH='+arch,'install'],'install')
check=execute(['make','PETSC_DIR='+str(prefix),'PETSC_ARCH=','check'],'check',600)
observed={k:s in text(check) for k,s in [('CPU','run successfully with 1 MPI process'),('MPI2','run successfully with 2 MPI processes')]}
if cuda_path:observed['CUDA']='run successfully with cuda' in text(check)
assert all(observed.values()) and not re.search(r'Possible (?:error|problem)',text(check)),text(check)
changed=[n for n,h in original.items() if not (source/n).is_file() or digest(source/n)!=h]
assert not changed
conf=source/arch/'lib/petsc/conf';evidence=R/(run_id+'_build_evidence');evidence.mkdir()
for p in (conf/'configure.log',conf/'make.log',conf/'petscvariables',source/arch/'include/petscconf.h'):shutil.copyfile(p,evidence/p.name)
ldd=execute(['ldd',prefix/'lib/libpetsc.so'],'ldd',30);assert 'not found' not in text(ldd)
resolved=re.search(r'libmpi[^\s]*\s+=>\s+(\S+)',text(ldd));assert resolved and Path(resolved[1]).resolve()==(mpiprefix/'lib/libmpi.so').resolve()
d.update(status='PASS',selftest=check,selftest_observed=observed,source_unmodified=True,source_files_verified=len(original),source_modifications=changed,library_sha256=digest(prefix/'lib/libpetsc.so'),MPI_link=resolved[1])
write('petsc_'+candidate+'_build',d)
write('petsc_'+run_id+'_source_integrity',dict(status='PASS',files=original,modifications=changed,source=str(source)))
print(candidate+' PETSc clean build and selftests PASS',flush=True)
