from common import *
import subprocess,time,difflib,shutil,resource
assert not LIB.exists()
old=json.loads((B/'BUILD_AUDIT.json').read_text());objects=list(old['object_sha256']);flags=old['commands'][0][:old['commands'][0].index('-c')];flags[1:1]=['-I'+str(A)]
src=R/'verification/lps_native_softening_event_resolved_integration_v1/20261013T235945Z'
sources={'bond_bpm_peri_variational.cpp':A/'bond_bpm_peri_variational.cpp','fix_lps_variational_dynamic.cpp':A/'fix_lps_variational_dynamic.cpp','force.cpp':src/'source/src/force.cpp','style_bond.cpp':src/'build/styles/style_bond.cpp'}
for p in objects:assert sha(p)==old['object_sha256'][p],p
start=time.perf_counter();commands=[]
with (O/'build/build.log').open('w') as log:
 for name,path in sources.items():
  obj=O/'build'/(name+'.o');cmd=flags+['-c',str(path),'-o',str(obj)];commands.append(cmd);subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
  match=[i for i,p in enumerate(objects) if Path(p).name==name+'.o'];assert len(match)==1;objects[match[0]]=str(obj)
 rsp=O/'build/objects.rsp';rsp.write_text('\n'.join(objects)+'\n');link=old['commands'][-1].copy();link[link.index('-o')+1]=str(O/'build/liblammps.so.0');link=[('@'+str(rsp)) if q.startswith('@') else q for q in link];commands.append(link);subprocess.run(link,cwd=old['link_cwd'],stdout=log,stderr=subprocess.STDOUT,check=True)
LIB.symlink_to('liblammps.so.0')
before=R/'analysis/clot_multifragment_crushing_mechanism_v1';patch=''
for p in A.glob('*'):
 if p.suffix not in ['.cpp','.h']:continue
 patch+=''.join(difflib.unified_diff((before/p.name).read_text().splitlines(True),p.read_text().splitlines(True),fromfile='baseline/'+p.name,tofile='new/'+p.name))
(O/'NATIVE_SOURCE.diff').write_text(patch)
dump(O/'BUILD_AUDIT.json',dict(status='PASS',library_sha256=sha(LIB),baseline_library_sha256=sha(OLDLIB),baseline_unchanged=sha(OLDLIB)==old['library_sha256'],compiled=list(sources),commands=commands,link_cwd=old['link_cwd'],object_sha256={p:sha(p) for p in objects},wall_s=time.perf_counter()-start,upstream_commit='8de817dd79bfe4525d5d39246a212d833e6dee07',source_sha256={p.name:sha(p) for p in A.glob('*') if p.suffix in ['.cpp','.h']},fixed_damage_force_kernel_unchanged=True))
print('ISOLATED_BUILD_PASS',sha(LIB),flush=True)
