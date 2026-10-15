"""Preformal ABI fix: rebuild generated creator for changed derived-class size."""
from common import *
import shutil,subprocess,time
assert not (O/'FORMAL_PILOT_LAUNCHED.json').exists()
backup=O/'preflight_registry_build_failure';backup.mkdir(exist_ok=False)
for name in ['BUILD_AUDIT.json','SOURCE_AND_LIBRARY_IDENTITY.json','PREFLIGHT_PROCESS.json']:shutil.copy2(O/name,backup/name)
shutil.copy2(LIB,backup/'failed_liblammps.so');
for name in ['parity.log','initial.npz','perturbed.npz']:
 if (O/'preflight'/name).exists():shutil.copy2(O/'preflight'/name,backup/name)
audit=json.loads((O/'BUILD_AUDIT.json').read_text());objects=list(audit['object_sha256']);cmd=audit['commands'][0];flags=cmd[:cmd.index('-c')]
src=NATIVE/'build/styles/style_bond.cpp';obj=O/'build/style_bond.cpp.o';command=flags+['-c',str(src),'-o',str(obj)];start=time.perf_counter()
with (O/'build/registry_repair.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
matches=[i for i,p in enumerate(objects) if Path(p).name=='style_bond.cpp.o'];assert len(matches)==1;objects[matches[0]]=str(obj)
rsp=O/'build/objects.rsp';rsp.write_text('\n'.join(objects)+'\n')
with (O/'build/registry_repair.log').open('a') as log:subprocess.run(audit['commands'][-1],cwd=audit['link_cwd'],stdout=log,stderr=subprocess.STDOUT,check=True)
audit['compiled_translation_units'].append('style_bond.cpp');audit['commands'].insert(-1,command);audit['object_sha256']={p:sha(p) for p in objects};audit['reused_objects']-=1;audit['registry_repair_wall_s']=time.perf_counter()-start;dump(O/'BUILD_AUDIT.json',audit)
ident=json.loads((O/'SOURCE_AND_LIBRARY_IDENTITY.json').read_text());ident['new_library_sha256']=sha(LIB);ident['generated_registry_recompiled']=True;ident['registry_source_sha256']=sha(src);dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',ident)
dump(O/'REGISTRY_ABI_REPAIR.json',dict(status='REBUILT_PENDING_RECHECK',cause='Generated style_bond.cpp creator retained sizeof the old single-type class; initial force checks passed but destruction aborted. Recompiled unchanged generated creator with the new header.',preformal=True,formal_cases=0,physics_changed=False,previous_library_sha256=sha(backup/'failed_liblammps.so'),new_library_sha256=sha(LIB),backup=str(backup)))
print('REGISTRY_ABI_REBUILT',sha(LIB),flush=True)
