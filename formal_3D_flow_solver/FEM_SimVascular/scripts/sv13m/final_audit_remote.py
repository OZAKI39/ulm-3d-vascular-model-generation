"""Read-only history verification and preservation of the experimental native stack."""
import shutil,subprocess,tarfile
from environment_remote import snapshot
from runner_remote import *
before=load('pre_install_environment');after=snapshot();write('final_environment',after)
keys=('driver_files','cuda_default_link','cuda13_prefix','cuda13_files','cuda126_prefix','cuda126_files','default_compilers','ld_configuration','historical_gpu_stacks')
checks={k:before[k]==after[k] for k in keys}
write('remote_preservation',{'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'historical_stages':['SV1.3J','SV1.3L','SV1.3G MPI']});assert all(checks.values())
w=load('compatibility_winner');integrity=load(w['repair_iteration']+'_source_integrity');source=Path(w['source']);original=integrity['original_files'];patched=integrity['patched_files']
modified=[n for n,h in original.items() if not (source/n).is_file() or digest(source/n)!=h]
unexpected=[n for n,h in patched.items() if not (source/n).is_file() or digest(source/n)!=h]
# Newly added files under source/include are forbidden; generated objects live in PETSC_ARCH.
added=[str(p.relative_to(source)) for folder in ('src','include') for p in (source/folder).rglob('*') if p.is_file() and str(p.relative_to(source)) not in original]
sv=load('svmp_source_before');svsource=Path(sv['source']);svchanged=[n for n,h in sv['files'].items() if not (svsource/n).is_file() or digest(svsource/n)!=h]
state=subprocess.check_output(['git','-C',str(svsource),'status','--porcelain'],text=True)
assert not unexpected and not added and not svchanged and not state.strip()
assert sorted(modified)==sorted(integrity['expected_changed_files'])
write('patch_integrity',dict(status='PASS',modified_files=modified,expected_files=integrity['expected_changed_files'],unexpected_changes=unexpected,added_source_files=added,solver_changes=svchanged,solver_git_status=state,source_files=len(original),solver_source_files=len(sv['files']),before_tree_sha256=integrity['before_tree_sha256'],after_tree_sha256=integrity['after_tree_sha256'],fresh_source=True,previous_objects_reused=False,original_archive_sha256=digest(integrity['archive']),patch_sha256=integrity['patch_sha256']))
O=BASE/'outputs/native';O.mkdir(exist_ok=True);archive=O/'petsc-cuda123-ghostfix.tar.gz';assert not archive.exists()
with tarfile.open(archive,'w:gz') as t:t.add(w['prefix'],arcname='petsc-cuda123-ghostfix')
paths=[archive,Path(integrity['archive']),Path(load('svmp_gpu_build')['executable']),B/'petsc_cuda_smoke',B/'petsc_cuda_ghost_probe_repair_03',B/'petsc_cuda_ghost_coherence',B/'petsc_finalize_probe']
write('native_artifacts',dict(status='PASS',stack_label=w['stack_label'],artifacts=[dict(remote_path=str(p),sha256=digest(p),bytes=p.stat().st_size) for p in paths],dependency_note='MPI and CUDA remain unchanged in Stage L/J; native binaries retain the recorded remote RPATH. Not a WSL installed stack.'))
print('Remote history, expected patch scope and unchanged solver source PASS; native stack archived.',flush=True)
