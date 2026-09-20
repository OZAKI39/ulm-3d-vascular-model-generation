from runner_remote import *
sv=load('svmp_reuse_build');w=sv['PETSc_build']
r=run(['cmake','--build',Path(sv['build'])/'svMultiPhysics-build','--parallel','4'],'stop_timestamp_inner_make',timeout=3600,cuda=w['candidate_wrapper'],extra_env={'LD_LIBRARY_PATH':sv['runtime_library_path'],'OMPI_CC':w['CC'],'OMPI_CXX':w['CXX']});assert okay(r)
assert b'sv13q_stop_ack.json' in Path(sv['executable']).read_bytes()
sv['steps'].append(r);sv['executable_sha256']=digest(sv['executable']);write('svmp_reuse_build',sv)
print('Incremental inner binary includes actual stop acknowledgement instrumentation.',flush=True)
