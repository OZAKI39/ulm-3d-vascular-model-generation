from pathlib import Path
import subprocess,json,os,hashlib,time,datetime,traceback,urllib.request,tarfile,shutil
W=Path('/workspace/lammps_migration/new_20260916_091244');R=Path('/workspace/microbubble_lammps/results/lammps2026_bubble_interaction_20260916_091244');P=R/'provenance';S=W/'upstream/lammps'
def save(n,d):(P/n).write_text(json.dumps(d,indent=2)+'\n')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def cmd(a,n,env):
 with (P/n).open('w') as f:
  f.write('COMMAND_JSON='+json.dumps(a)+'\n');f.flush();p=subprocess.run(a,stdout=f,stderr=subprocess.STDOUT,env=env)
 if p.returncode:raise RuntimeError(n+' failed '+str(p.returncode))
try:
 t=time.time();target=json.loads((P/'TARGET_RELEASE_FROZEN.json').read_text());archive=W/'upstream/source.tar.gz';archive.parent.mkdir()
 save('BUILD_STATUS.json',{'state':'RUNNING','phase':'DOWNLOAD'})
 urllib.request.urlretrieve(target['source_download_url'],archive)
 with tarfile.open(archive) as tf:
  members=tf.getmembers();top=members[0].name.split('/')[0]
  assert all(not m.name.startswith('/') and '..' not in Path(m.name).parts for m in members)
  tf.extractall(archive.parent,filter='data')
 (archive.parent/top).rename(S)
 save('LAMMPS_2026_SOURCE_PROVENANCE.json',{**target,'source_archive_sha256':sha(archive),'download_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_root':str(S),'archive':str(archive),'version_header':(S/'src/version.h').read_text()})
 audit=P/'upstream_lubrication';audit.mkdir()
 for pattern in ['src/COLLOID/pair_lubricate*','doc/src/pair_lubricate*']:
  for f in S.glob(pattern):shutil.copy2(f,audit/f.name)
 env=os.environ.copy();env.update(HWLOC_COMPONENTS='-gl',OMP_NUM_THREADS='1',NVCC_WRAPPER_DEFAULT_COMPILER='/usr/bin/g++');env['PATH']='/usr/local/cuda/bin:'+env['PATH']
 base=['-DCMAKE_BUILD_TYPE=Release','-DBUILD_MPI=ON','-DBUILD_OMP=OFF','-DPKG_GRANULAR=ON','-DPKG_COLLOID=ON','-DPKG_LATBOLTZ=OFF','-DBUILD_SHARED_LIBS=OFF','-DCMAKE_C_COMPILER=/usr/bin/gcc']
 cfg={'cpu':base+['-DPKG_KOKKOS=OFF','-DCMAKE_CXX_COMPILER=/usr/bin/g++'],'gpu':base+['-DPKG_KOKKOS=ON','-DKokkos_ENABLE_CUDA=ON','-DKokkos_ENABLE_SERIAL=ON','-DKokkos_ENABLE_OPENMP=OFF','-DKokkos_ARCH_ADA89=ON','-DCMAKE_CXX_COMPILER='+str(S/'lib/kokkos/bin/nvcc_wrapper')]}
 records={}
 for backend,flags in cfg.items():
  b=W/('build_'+backend);bt=time.time();save('BUILD_STATUS.json',{'state':'RUNNING','phase':backend})
  cmd(['cmake','-S',str(S/'cmake'),'-B',str(b),'-G','Ninja']+flags,backend+'_configure.log',env)
  cmd(['cmake','--build',str(b),'--parallel','6'],backend+'_build.log',env)
  cmd([str(b/'lmp'),'-h'],backend+'_lmp_help.txt',env)
  shutil.copy2(b/'CMakeCache.txt',P/(backend+'_CMakeCache.txt'))
  records[backend]={'binary':str(b/'lmp'),'sha256':sha(b/'lmp'),'flags':flags,'wall_s':time.time()-bt}
  save('LAMMPS_2026_BUILD_PROVENANCE.json',records)
 save('BUILD_STATUS.json',{'state':'PASS','wall_s':time.time()-t})
except Exception as e:save('BUILD_STATUS.json',{'state':'FAIL','error':repr(e),'traceback':traceback.format_exc()});raise
