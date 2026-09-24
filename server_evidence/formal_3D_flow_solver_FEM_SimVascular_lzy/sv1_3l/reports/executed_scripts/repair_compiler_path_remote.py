"""Archive the failed harness attempt; restore J's GCC12 PATH lookup without source edits."""
import shutil
from runner_remote import *
source=BASE/'external/compat_cuda/petsc-3.19.6-cuda123-mpif'
previous=load('petsc_configure');assert previous['exit_code']==1
detail=(source/'configure.log').read_text(errors='replace')
assert 'unsupported GNU version' in detail
original=load('petsc_source_integrity_before')
assert all(digest(source/n)==h for n,h in original['original_files'].items())
archive=BASE/'outputs/setup_attempt_01';archive.mkdir(parents=True)
shutil.copyfile(source/'configure.log',L/'petsc_configure_attempt1_detail.log')
write('petsc_configure_attempt1',previous)
write('petsc_source_integrity_attempt1',original)
source.rename(archive/'petsc_source')
cuda_wrapper=BASE/'scripts/use_cuda123_mpif_env.sh'
cuda_wrapper.chmod(0o755)
code=B/'nvcc_host_check.cu'
code.write_text('static_assert(__GNUC__==12,"NVCC must use GCC12");\nstatic_assert(__cplusplus==201703L,"C++17 required");\nint probe(){return 42;}\n')
check=run(['nvcc','-std=c++17','-arch=sm_89','-c',code,'-o',B/'nvcc_host_check.o'],'nvcc_host_path_check',cuda=cuda_wrapper)
assert okay(check)
host=run(['gcc','-dumpfullversion'],'nvcc_host_gcc_version',cuda=cuda_wrapper)
assert text(host).strip().split('.')[0]=='12'
write('compiler_path_correction',{'status':'PASS','first_attempt':'petsc_configure_attempt1.json',
    'cause':'New CUDA wrapper initially omitted the J private host-bin path; NVCC discovered system GCC13.',
    'correction':'Restored the existing J host-bin lookup selecting GCC12; no unsupported-compiler flag or source patch.',
    'failed_source_archive':str(archive/'petsc_source'),'fresh_source_required':True,
    'wrapper_sha256':digest(cuda_wrapper),'compile_check':check,'host_gcc':text(host).strip()})
