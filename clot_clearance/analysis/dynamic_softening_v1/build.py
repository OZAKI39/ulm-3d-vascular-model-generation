from common import *
import subprocess,shutil,time,resource,os,difflib
start=time.perf_counter();assert (O/'DYNAMIC_DAMAGE_COMMIT_DESIGN_FREEZE.json').exists()
patch=[];count=0
names=['bond_bpm_peri_variational.h','bond_bpm_peri_variational.cpp','fix_lps_variational_dynamic.h','fix_lps_variational_dynamic.cpp']
for name in names:
 src=A/name;baseline=PREV/'source/src/BPM'/name
 before=baseline.read_text().splitlines(True) if baseline.exists() else []
 shutil.copy2(src,O/'source/src/BPM'/name)
 patch.extend(difflib.unified_diff(before,src.read_text().splitlines(True),fromfile='baseline/src/BPM/'+name,tofile='dynamic/src/BPM/'+name))
book=O/'source/src/.gitignore';base=(PREV/'source/src/.gitignore').read_text();added=['/'+n for n in names if '/'+n not in base.splitlines()]
newbook=base+'\n# Isolated variational research package copies (legacy make)\n'+'\n'.join(added)+'\n';book.write_text(newbook)
patch.extend(difflib.unified_diff(base.splitlines(True),newbook.splitlines(True),fromfile='baseline/src/.gitignore',tofile='dynamic/src/.gitignore'))
(O/'NATIVE_SOURCE_PATCH.diff').write_text(''.join(patch))
cmd=['cmake','-S',str(O/'source/cmake'),'-B',str(O/'build'),'-D','CMAKE_BUILD_TYPE=Release','-D','BUILD_MPI=ON','-D','BUILD_OMP=ON','-D','BUILD_SHARED_LIBS=ON','-D','PKG_BPM=ON','-D','PKG_MOLECULE=ON','-D','MPI_CXX_COMPILER=/workspace/y_clot/deps/openmpi/bin/mpicxx']
records=[]
for name,c in [('configure',cmd),('compile',['cmake','--build',str(O/'build'),'-j','12'])]:
 t=time.perf_counter()
 with (O/(name+'.log')).open('w') as f:
  z=subprocess.Popen(c,stdout=f,stderr=subprocess.STDOUT);peak_rss=0;peak_processes=0;cpu_ticks={};peak_cpu=0;last_ticks=0;last_t=time.perf_counter()
  while z.poll() is None:
   todo=[z.pid];rss=0;seen=[]
   while todo:
    pid=todo.pop()
    try:
     stat=Path(f'/proc/{pid}/stat').read_text().split();rss+=int(stat[23])*os.sysconf('SC_PAGE_SIZE');cpu_ticks[pid]=max(cpu_ticks.get(pid,0),int(stat[13])+int(stat[14]));seen.append(pid);todo.extend(map(int,Path(f'/proc/{pid}/task/{pid}/children').read_text().split()))
    except (FileNotFoundError,ProcessLookupError):pass
   now=time.perf_counter();total=sum(cpu_ticks.values());peak_cpu=max(peak_cpu,(total-last_ticks)/os.sysconf('SC_CLK_TCK')/max(now-last_t,1e-6));last_ticks=total;last_t=now;peak_rss=max(peak_rss,rss);peak_processes=max(peak_processes,len(seen))
   if now-t>900:z.terminate();raise RuntimeError('compile budget exceeded')
   time.sleep(.5)
 records.append(dict(stage=name,cmd=c,wall_s=time.perf_counter()-t,returncode=z.returncode,sampled_peak_process_tree_RSS_bytes=peak_rss,sampled_peak_CPU_cores=peak_cpu,sampled_peak_processes=peak_processes,child_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss));dump(O/'BUILD_RESOURCES.json',dict(stages=records,wall_s=time.perf_counter()-start,compile_jobs=12))
 if z.returncode:print((O/(name+'.log')).read_text()[-10000:]);raise SystemExit(z.returncode)
dump(O/'ISOLATED_BINARY_IDENTITY.json',dict(library=str(LIB),sha256=sha(LIB),size=LIB.stat().st_size,source_patch_sha256=sha(O/'NATIVE_SOURCE_PATCH.diff'),new_sources={p.name:sha(p) for p in (O/'source/src/BPM').glob('*variational*')},original_src_files_verified=count))
print('BUILD_COMPLETE',sha(LIB),time.perf_counter()-start,flush=True)
