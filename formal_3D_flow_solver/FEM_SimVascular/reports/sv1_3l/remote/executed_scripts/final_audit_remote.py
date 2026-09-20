"""Preserve old stacks byte-for-byte and archive the new native MPI/PETSc/solver."""
import shutil,subprocess,tarfile
from environment_remote import snapshot
from runner_remote import *
before=load('pre_install_environment');after=snapshot();write('final_environment',after)
keys=('driver_files','cuda_default_link','cuda13_prefix','cuda13_files','cuda126_prefix',
      'cuda126_files','default_compilers','ld_configuration','historical_gpu_stacks')
checks={key:before[key]==after[key] for key in keys}
probe=load('compiler_probe')
checks['existing_Fortran_paths']=all(shutil.which(n)==probe[n] for n in ('gfortran','gfortran-13'))
write('environment_preservation',{'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks})
assert all(checks.values()),'PRESERVATION_FAILURE'
w=load('compatibility_winner');original=load('petsc_source_integrity_before');source=Path(w['source'])
changed=[n for n,h in original['original_files'].items() if not (source/n).is_file() or digest(source/n)!=h]
sv=load('svmp_source_before');svsource=Path(sv['source'])
svchanged=[n for n,h in sv['files'].items() if not (svsource/n).is_file() or digest(svsource/n)!=h]
state=subprocess.check_output(['git','-C',str(svsource),'status','--porcelain'],text=True)
write('final_source_integrity',{'status':'PASS' if not changed and not svchanged and not state.strip() else 'FAIL',
    'PETSc_original_files':len(original['original_files']),'PETSc_changes':changed,
    'svMultiPhysics_tracked_files':len(sv['files']),'svMultiPhysics_changes':svchanged,'svMultiPhysics_git_status':state,
    'PETSc_upgraded':False,'source_patches_applied':False})
assert not changed and not svchanged and not state.strip()
evidence=R/'mpi_build_evidence';evidence.mkdir(exist_ok=True)
mpisource=BASE/'external/openmpi-4.1.6'
for p in (mpisource/'config.log',mpisource/'opal/include/opal_config.h',mpisource/'ompi/include/mpi.h'):
    if p.is_file():shutil.copyfile(p,evidence/p.name)
O=BASE/'outputs/native';O.mkdir(exist_ok=True)
paths=[]
for name,prefix in [('gpu_mpi_fortran',load('mpi_fortran_build')['prefix']),('petsc-cuda123-mpif',w['prefix'])]:
    archive=O/(name+'.tar.gz');assert not archive.exists()
    with tarfile.open(archive,'w:gz') as t:t.add(prefix,arcname=name)
    paths.append(archive)
paths += [Path(load('svmp_gpu_build')['executable']),B/'mpi_hello',B/'mpi_fortran_smoke',B/'mpi_datatype_probe',B/'petsc_cuda_smoke']
write('native_artifacts',{'status':'PASS','artifacts':[{'remote_path':str(p),'sha256':digest(p),'bytes':p.stat().st_size} for p in paths],
    'CUDA_installation':'Unchanged Stage J CUDA12.3.2, source recovery manifest retained in J',
    'dependency_note':'Native binaries retain recorded remote paths; archives are not a WSL installation.'})
print('Old stack and source preservation PASS; new MPI/PETSc native artifacts archived',flush=True)
