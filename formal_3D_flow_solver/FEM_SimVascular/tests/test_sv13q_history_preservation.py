from sv13q_support import *
def test_historical_and_native_stack_protected():
 for n in ('preservation_audit','remote/remote_preservation','remote/final_input_integrity','native_artifact_mirror'):assert read(n)['status']=='PASS'
 p=read('source_patch')
 assert set(p['changed_files'])=={'Code/Source/solver/petsc_impl.cpp','Code/Source/solver/main.cpp'}
 assert p['added_files']==['Code/Source/solver/sv13q_reuse.h']
 for n,h in p['after'].items():assert digest(ROOT/p['source']/n)==h
