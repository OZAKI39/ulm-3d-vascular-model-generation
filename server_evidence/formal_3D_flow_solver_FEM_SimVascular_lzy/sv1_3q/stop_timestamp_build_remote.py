from runner_remote import *
sv=load('svmp_reuse_build');w=sv['PETSc_build'];source=Path(sv['source']);extra={'LD_LIBRARY_PATH':sv['runtime_library_path'],'OMPI_CC':w['CC'],'OMPI_CXX':w['CXX']}
write('svmp_reuse_build_before_stop_timestamp',sv)
r=run(['patch','--batch','--fuzz=0','-p1','-i',BASE/'patches/stop_ack_timestamp.patch'],'stop_timestamp_patch',cwd=source);assert okay(r)
r=run(['cmake','--build',sv['build'],'--parallel','4'],'stop_timestamp_make',timeout=3600,cuda=w['candidate_wrapper'],extra_env=extra);assert okay(r)
patch=json.loads((BASE/'configs/source_patch.json').read_text())
for n,h in patch['after'].items():assert digest(source/n)==h,n
sv.update(executable_sha256=digest(sv['executable']),source_manifest_sha256=digest(BASE/'configs/source_patch.json'));sv['steps'].append(r);write('svmp_reuse_build',sv)
print('Stage Q final instrumented binary ready; no CFD yet.',flush=True)
