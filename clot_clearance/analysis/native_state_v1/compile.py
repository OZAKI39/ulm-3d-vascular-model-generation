from pathlib import Path
import subprocess,time,resource,json,sys
C=json.loads(Path(__file__).with_name('config.json').read_text());O=Path.cwd()/'verification'/C['task']/C['run'];start=time.perf_counter();owner=len(sys.argv)>1
commands=[['cmake','-S',str(O/('source_owner_fix/cmake' if owner else 'source/cmake')),'-B',str(O/('build_owner_fix' if owner else 'build')),'-DCMAKE_BUILD_TYPE=Release','-DBUILD_MPI=ON','-DBUILD_SHARED_LIBS=ON','-DPKG_BPM=ON','-DPKG_MOLECULE=ON','-DDOWNLOAD_POTENTIALS=OFF','-DMPI_CXX_COMPILER=/workspace/y_clot/deps/openmpi/bin/mpicxx'],['cmake','--build',str(O/('build_owner_fix' if owner else 'build')),'--parallel','8']]
code=0
try:
 for cmd in commands:subprocess.run(cmd,check=True)
except subprocess.CalledProcessError as e:code=e.returncode;raise
finally:
 r=resource.getrusage(resource.RUSAGE_CHILDREN);(O/('OWNER_COMPILE_PERFORMANCE.json' if owner else 'COMPILE_PERFORMANCE.json')).write_text(json.dumps({'commands':commands,'exit_code':code,'wall_s':time.perf_counter()-start,'cpu_user_s':r.ru_utime,'cpu_system_s':r.ru_stime,'max_child_RSS_KiB':r.ru_maxrss,'parallelism':8},indent=2))
