from pathlib import Path
import subprocess, json, os, hashlib, time, datetime, traceback
W=Path(__file__).resolve().parents[1]; S=W/'upstream/lammps'; P=W/'build_provenance'
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name,d): (P/name).write_text(json.dumps(d,indent=2)+'\n')
def command(args,name,env):
 with (P/name).open('w') as f:
  f.write('COMMAND_JSON='+json.dumps(args)+'\n');f.flush()
  p=subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,env=env)
 if p.returncode: raise RuntimeError(f'{name} exited {p.returncode}')
try:
 start=time.time(); env=os.environ.copy();env['PATH']='/usr/local/cuda/bin:'+env['PATH'];env['NVCC_WRAPPER_DEFAULT_COMPILER']='/usr/bin/g++'
 base=['-DCMAKE_BUILD_TYPE=Release','-DBUILD_MPI=ON','-DBUILD_OMP=OFF','-DPKG_GRANULAR=ON','-DPKG_LATBOLTZ=OFF','-DBUILD_SHARED_LIBS=OFF','-DCMAKE_C_COMPILER=/usr/bin/gcc']
 configs={'cpu':base+['-DPKG_KOKKOS=OFF','-DCMAKE_CXX_COMPILER=/usr/bin/g++'], 'gpu':base+['-DPKG_KOKKOS=ON','-DKokkos_ENABLE_CUDA=ON','-DKokkos_ENABLE_SERIAL=ON','-DKokkos_ENABLE_OPENMP=OFF','-DKokkos_ARCH_ADA89=ON','-DCMAKE_CXX_COMPILER='+str(S/'lib/kokkos/bin/nvcc_wrapper')]}
 records={}
 for kind,flags in configs.items():
  b=W/('build_'+kind); ts=time.time();save('BUILD_STATUS.json',{'state':'RUNNING','phase':kind,'start_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
  command(['cmake','-S',str(S/'cmake'),'-B',str(b),'-G','Ninja']+flags,kind+'_configure.log',env)
  command(['cmake','--build',str(b),'--parallel','8'],kind+'_build.log',env)
  command([str(b/'lmp'),'-h'],kind+'_lmp_help.txt',env)
  helptext=(P/(kind+'_lmp_help.txt')).read_text();cache=(b/'CMakeCache.txt').read_text(); (P/(kind+'_CMakeCache.txt')).write_text(cache)
  assert 'PKG_LATBOLTZ:BOOL=OFF' in cache
  assert 'sphere' in helptext and 'nve/sphere' in helptext and 'gran/hooke/history' in helptext
  if kind=='gpu': assert all(x in helptext for x in ['KOKKOS','nve/sphere/kk','gran/hooke/history/kk']) and 'Kokkos_ARCH_ADA89:BOOL=ON' in cache
  records[kind]={'binary':str(b/'lmp'),'sha256':sha(b/'lmp'),'cmake_flags':flags,'wall_seconds':time.time()-ts,'status':'PASS'}
 meta=json.loads((P/'OFFICIAL_RELEASE_METADATA.json').read_text()); source=json.loads((P/'LAMMPS_SOURCE_PROVENANCE.json').read_text())
 save('LAMMPS_BUILD_PROVENANCE.json',{'status':'PASS','release':meta,'source':source,'builds':records,'source_common_to_both':True,'CMAKE_VERSION':subprocess.check_output(['cmake','--version'],text=True),'C_COMPILER':subprocess.check_output(['gcc','--version'],text=True),'CXX_COMPILER':subprocess.check_output(['g++','--version'],text=True),'CUDA_VERSION':subprocess.check_output(['/usr/local/cuda/bin/nvcc','--version'],text=True),'MPI_VERSION':subprocess.check_output(['mpirun','--version'],text=True),'GPU_MODEL':'NVIDIA GeForce RTX 4090','GPU_COMPUTE_CAPABILITY':'8.9','KOKKOS_ARCHITECTURE':'ADA89','ENABLED_PACKAGES':{'cpu':['GRANULAR'],'gpu':['GRANULAR','KOKKOS']},'LATBOLTZ_ENABLED':'NO','build_wall_seconds':time.time()-start})
 save('BUILD_STATUS.json',{'state':'PASS','elapsed_seconds':time.time()-start})
except Exception as e:
 save('BUILD_STATUS.json',{'state':'FAIL','error':repr(e),'traceback':traceback.format_exc()});raise
