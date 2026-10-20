"""Old material implementation, ONLY parser node-capacity guard enlarged for Test A."""
from common import *
import subprocess,time,difflib
out=O/'baseline_capacity';out.mkdir(exist_ok=True);assert not (out/'liblammps.so').exists()
old=json.loads((PREVIOUS/'BUILD_AUDIT.json').read_text());objects=list(old['object_sha256'])
for p in objects:assert sha(p)==old['object_sha256'][p]
src=PREVIOUS_A/'bond_bpm_peri_variational.cpp';text=src.read_text();new=text.replace('n>4096','n>1000000');assert new!=text
(out/src.name).write_text(new);(out/'CAPACITY_ONLY.diff').write_text(''.join(difflib.unified_diff(text.splitlines(True),new.splitlines(True),fromfile=str(src),tofile=str(out/src.name))))
flags=old['commands'][0][:old['commands'][0].index('-c')];flags[1:1]=['-I'+str(PREVIOUS_A)]
cmd=flags+['-c',str(out/src.name),'-o',str(out/(src.name+'.o'))];start=time.perf_counter()
with (out/'build.log').open('w') as log:
 subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
 idx=[i for i,p in enumerate(objects) if Path(p).name==src.name+'.o'];assert len(idx)==1;objects[idx[0]]=str(out/(src.name+'.o'))
 rsp=out/'objects.rsp';rsp.write_text('\n'.join(objects)+'\n');link=old['commands'][-1].copy();link[link.index('-o')+1]=str(out/'liblammps.so.0');link=[('@'+str(rsp)) if q.startswith('@') else q for q in link]
 subprocess.run(link,cwd=old['link_cwd'],stdout=log,stderr=subprocess.STDOUT,check=True)
(out/'liblammps.so').symlink_to('liblammps.so.0')
dump(O/'BASELINE_CAPACITY_AUDIT.json',dict(status='PASS',reason='Latest old library limits reference and restart to 4096 nodes; new specimen has 10000. Exactly two parsing guards enlarged; all material, history and force code unchanged.',baseline_sha256=sha(OLDLIB),comparison_sha256=sha(out/'liblammps.so'),source_sha256=sha(out/src.name),commands=[cmd,link],wall_s=time.perf_counter()-start))
print('CAPACITY_BASELINE_PASS',flush=True)
