from common import *
import subprocess,time,difflib

assert not LIB.exists()
old=json.loads((B/'BUILD_AUDIT.json').read_text());objects=list(old['object_sha256'])
assert sha(OLDLIB)==old['library_sha256']
for p in objects:assert sha(p)==old['object_sha256'][p],p
flags=old['commands'][0][:old['commands'][0].index('-c')];flags[1:1]=['-I'+str(A)]
start=time.perf_counter();commands=[];source=A/'fix_lps_variational_dynamic.cpp';obj=O/'build/fix_lps_variational_dynamic.cpp.o'
with (O/'build/build.log').open('w') as log:
    cmd=flags+['-c',str(source),'-o',str(obj)];commands.append(cmd);subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    matches=[i for i,p in enumerate(objects) if Path(p).name==obj.name];assert len(matches)==1;objects[matches[0]]=str(obj)
    response=O/'build/objects.rsp';response.write_text('\n'.join(objects)+'\n')
    link=old['commands'][-1].copy();link[link.index('-o')+1]=str(O/'build/liblammps.so.0');link=[('@'+str(response)) if p.startswith('@') else p for p in link];commands.append(link)
    subprocess.run(link,cwd=old['link_cwd'],stdout=log,stderr=subprocess.STDOUT,check=True)
LIB.symlink_to('liblammps.so.0')
assert sha(OLDLIB)==old['library_sha256']
patch=''.join(difflib.unified_diff((HA/source.name).read_text().splitlines(True),source.read_text().splitlines(True),fromfile='accepted_single_indenter/'+source.name,tofile='multisite/'+source.name))
(O/'NATIVE_SOURCE.diff').write_text(patch)
dump(O/'BUILD_AUDIT.json',dict(status='PASS',library_sha256=sha(LIB),baseline_library_sha256=sha(OLDLIB),baseline_unchanged=True,compiled=[source.name],commands=commands,link_cwd=old['link_cwd'],object_sha256={p:sha(p) for p in objects},wall_s=time.perf_counter()-start,upstream_commit=old['upstream_commit'],source_sha256={p.name:sha(p) for p in [source,A/'fix_lps_variational_dynamic.h']},material_object_reused_unchanged=True,change_scope='Only admit two named fix-indent sources after physical force conversion. Coordinator layout, damage kernel, original-ID history and restart unchanged.'))
print('ISOLATED_MULTISITE_BUILD_PASS',sha(LIB))
