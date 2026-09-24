"""Clean rebuild of the exact J winner, changing only MPI and isolated build paths."""
import hashlib,shutil,tarfile,sys
from runner_remote import *
mpi=load('mpi_application_gate');assert mpi['status']=='PASS'
j=BASE.parent/'sv1_3j';old=json.loads((j/'reports/compatibility_winner.json').read_text())
assert old['status']=='PASS' and old['CUDA_version']=='12.3.2'
E=BASE/'external/compat_cuda';E.mkdir(exist_ok=True)
source=E/'petsc-3.19.6-cuda123-mpif';prefix=E/'petsc-cuda123-mpif';arch='arch-sv13l-cuda123-mpif'
archive=BASE.parent/'sv1_3/petsc-3.19.6.tar.gz'
assert digest(archive)==old['source_archive_sha256'] and not source.exists()
original={}
with tarfile.open(archive) as t:
    for member in t:
        if member.isfile():original[str(Path(member.name).relative_to('petsc-3.19.6'))]=hashlib.sha256(t.extractfile(member).read()).hexdigest()
    temp=E/('source_extract'+('_retry1' if len(sys.argv)>1 else ''));temp.mkdir();t.extractall(temp,filter='data')
(temp/'petsc-3.19.6').rename(source)
assert all(digest(source/n)==h for n,h in original.items())
write('petsc_source_integrity_before',{'status':'PASS','archive_sha256':digest(archive),'source':str(source),'original_files':original,'source_files':len(original),'modifications':[],'fresh_source':True,'previous_objects_reused':False,'PETSC_ARCH':arch})
mpiprefix=Path(mpi['prefix']);wrapper=mpi['cuda_wrapper']
replace={'--prefix=':'--prefix='+str(prefix),'PETSC_ARCH=':'PETSC_ARCH='+arch,
         '--with-cc=':'--with-cc='+str(mpiprefix/'bin/mpicc'),
         '--with-cxx=':'--with-cxx='+str(mpiprefix/'bin/mpicxx'),
         '--with-mpiexec=':'--with-mpiexec='+mpi['working_launcher']}
command=[next((replacement for start,replacement in replace.items() if arg.startswith(start)),arg) for arg in old['configure_command']]
diff=[{'before':a,'after':b} for a,b in zip(old['configure_command'],command) if a!=b]
assert len(diff)==5
extra={'OMPI_CC':'/usr/bin/gcc-12','OMPI_CXX':'/usr/bin/g++-12','PETSC_OPTIONS':'-use_gpu_aware_mpi 0'}
cfg=run(command,'petsc_configure_retry1' if len(sys.argv)>1 else 'petsc_configure',cwd=source,timeout=3600,cuda=wrapper,extra_env=extra)
write('petsc_configure',dict(cfg,status='PASS' if okay(cfg) else 'FAIL',configuration_diff_from_J=diff))
if not okay(cfg):raise SystemExit('PETSC_REBUILD_WITH_MPI_FORTRAN_FAIL')
conf=source/arch/'lib/petsc/conf';variables=(conf/'petscvariables').read_text()
for key in ('CXX_FLAGS','CUDAC_FLAGS'):
    assert '-std=c++17' in next(line for line in variables.splitlines() if line.startswith(key+' ='))
assert '#define PETSC_HAVE_CUDA 1' in (source/arch/'include/petscconf.h').read_text()
make=run(['make','-j8','V=1','PETSC_DIR='+str(source),'PETSC_ARCH='+arch,'all'],'petsc_make',cwd=source,timeout=3600,cuda=wrapper,extra_env=extra)
changed=[n for n,h in original.items() if not (source/n).is_file() or digest(source/n)!=h]
write('petsc_source_integrity_after',{'status':'FAIL' if changed else 'PASS','modifications':changed,'source_files':len(original)})
evidence=R/'petsc_build_evidence';evidence.mkdir()
for p in (conf/'configure.log',conf/'make.log',conf/'petscvariables',source/arch/'include/petscconf.h'):
    if p.is_file():shutil.copyfile(p,evidence/p.name)
d={'status':'PASS' if okay(make) and not changed else 'FAIL','version':'3.19.6','key':'cuda123-mpif','source':str(source),'prefix':str(prefix),'PETSC_ARCH':arch,
   'source_archive_sha256':digest(archive),'source_unmodified':not changed,'source_files_verified':len(original),
   'CUDA_version':'12.3.2','CUDA_prefix':old['CUDA_prefix'],'host_compiler':load('compiler_environment')['compilers']['g++-12']['number'],
   'MPI_prefix':str(mpiprefix),'configure':cfg,'make':make,'configure_command':command,'configuration_diff_from_J':diff,'candidate_wrapper':wrapper}
write('petsc_rebuild',d)
if d['status']!='PASS':raise SystemExit('PETSC_REBUILD_WITH_MPI_FORTRAN_FAIL')
write('compatibility_winner',d)
configs=BASE/'configs';configs.mkdir(exist_ok=True)
(configs/'mpi_resolution.json').write_text(json.dumps(mpi,indent=2))
print('Clean PETSc CUDA12.3 MPI-Fortran build PASS; original source unchanged',flush=True)
