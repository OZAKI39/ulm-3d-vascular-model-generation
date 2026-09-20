"""Read-only preservation checks and stage-local native artifact archives."""
import json,subprocess,tarfile
from pathlib import Path
from environment_remote import snapshot
from runner_remote import *
before=load('pre_install_environment');after=snapshot()
write('final_environment',after)
keys=('driver_files','cuda_default_link','cuda13_prefix','cuda13_files','cuda126_prefix','cuda126_files','default_compilers','ld_configuration')
checks={k:before[k]==after[k] for k in keys}
mpi=json.loads((BASE/'configs/mpi_resolution.json').read_text());native=json.loads((BASE/'configs/mpi_native_artifacts.json').read_text())
checks['MPI_artifacts']=all(digest(a['realpath'])==a['sha256'] for a in native)
checks['MPI_wrapper']=digest(mpi['working_launcher'])==mpi['wrapper_sha256']
write('environment_preservation',{'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks})
assert all(checks.values()),'PRESERVATION_FAILURE'
w=load('compatibility_winner');original=load(w['key']+'_source_integrity_before');source=Path(w['source'])
changed=[n for n,h in original['original_files'].items() if not (source/n).is_file() or digest(source/n)!=h]
sv=load('svmp_source_before');svsource=Path(sv['source'])
svchanged=[n for n,h in sv['files'].items() if not (svsource/n).is_file() or digest(svsource/n)!=h]
state=subprocess.check_output(['git','-C',str(svsource),'status','--porcelain'],text=True)
write('final_source_integrity',{'status':'PASS' if not changed and not svchanged and not state.strip() else 'FAIL',
    'PETSc_original_files':len(original['original_files']),'PETSc_changes':changed,
    'svMultiPhysics_tracked_files':len(sv['files']),'svMultiPhysics_changes':svchanged,'svMultiPhysics_git_status':state,
    'PETSc_upgraded':False,'source_patches_applied':False})
assert not changed and not svchanged and not state.strip()
O=BASE/'outputs/native';O.mkdir(exist_ok=True)
archive=O/'petsc_cuda123_install.tar.gz'
assert not archive.exists()
with tarfile.open(archive,'w:gz') as t:t.add(w['prefix'],arcname='petsc-cuda123')
svbuild=load('svmp_gpu_build');paths=[archive,Path(svbuild['executable']),B/'cuda123_kernel_smoke',B/'petsc_cuda_smoke',B/'mpi_datatype_probe']
write('native_artifacts',{'status':'PASS','artifacts':[{'remote_path':str(p),'sha256':digest(p),'bytes':p.stat().st_size} for p in paths],
    'CUDA_installer_mirrored':False,'CUDA_installer_manifest':'cuda123_source_manifest.json',
    'MPI_install_already_mirrored_in_stage':'SV1.3G'})
print('Preservation and source audits PASS; native winner artifacts archived',flush=True)
