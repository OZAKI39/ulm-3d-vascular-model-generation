from pathlib import Path
import json,sys,hashlib,datetime
sys.path.insert(0,str(Path(__file__).parent));import remote_exec as rex
s=rex.STATE;r=Path(s['local_result']);p=r/'provenance'
old=Path('/home/lzy/projects/compre_output')
contracts={
 'engine':old/'lammps_particle_engine/20260915_214759/contracts/MINIMAL_LAMMPS_PARTICLE_ENGINE_CONTRACT.json',
 'coupling':old/'palabos_lammps_coupling_v0/20260915_224019/contracts/PALABOS_LAMMPS_COUPLING_V0_CONTRACT.json',
 'passive':old/'passive_transport_v0/20260915_233516/contracts/PASSIVE_MICROBUBBLE_TRANSPORT_V0_CONTRACT.json',
 'passive_timestep':old/'passive_transport_v0/20260915_233516/contracts/PASSIVE_TRANSPORT_TIMESTEP_CONTRACT.json'}
c={'contract_name':'LAMMPS_2026_MIGRATION_CONTRACT','status':'FROZEN_BEFORE_RESULTS','frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'target':json.loads((p/'TARGET_RELEASE_FROZEN.json').read_text()),'phase_order':['A_NEW_BUILD_AND_ALL_REGRESSIONS','B_PROMOTE_AND_RETIRE_OLD_ACTIVE_ONLY','C_INTERACTION_DEVELOPMENT'],'baseline_contracts':{k:{'sha256':hashlib.sha256(v.read_bytes()).hexdigest(),'contents':json.loads(v.read_text())}for k,v in contracts.items()},'required':['source_provenance','cpu_build','gpu_build','source_correction','diameter','ballistic','technical_contact','polydisperse_1000','mpi','kokkos','coupling_uniform_linear_force_invalid','passive_case0_1_2_3_4'],'new_vs_old_case2_position_max_m':1e-10,'same_fields_populations_initial_conditions':True,'no_old_deletion_before_all_pass':True,'history_immutable':True,'math_compatibility_changes':'FORBIDDEN; registration/header API only if necessary; patch recorded','float':'double','C_adv':0.25,'integrator':'RK2_MIDPOINT','failure_policy':'preserve old active installation and stop dependent phases'}
(r/'contracts').mkdir();(r/'scripts').mkdir();(r/'contracts/LAMMPS_2026_MIGRATION_CONTRACT.json').write_text(json.dumps(c,indent=2)+'\n')
worker='''from pathlib import Path
import subprocess,json,os,hashlib,time,datetime,traceback,urllib.request,tarfile,shutil
W=Path(__WORK__);R=Path(__RESULT__);P=R/'provenance';S=W/'upstream/lammps'
def save(n,d):(P/n).write_text(json.dumps(d,indent=2)+'\\n')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def cmd(a,n,env):
 with (P/n).open('w') as f:
  f.write('COMMAND_JSON='+json.dumps(a)+'\\n');f.flush();p=subprocess.run(a,stdout=f,stderr=subprocess.STDOUT,env=env)
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
'''.replace('__WORK__',repr(s['remote_work'])).replace('__RESULT__',repr(s['remote_result']))
(r/'scripts/build_worker.py').write_text(worker)
files={str(f.relative_to(r)):f.read_text() for f in r.rglob('*') if f.is_file()}
remote='''from pathlib import Path
import json,subprocess,os,hashlib
S=__STATE__;W=Path(S['remote_work']);R=Path(S['remote_result'])
assert not W.exists() and not R.exists()
W.mkdir(parents=True);R.mkdir(parents=True)
for name,content in __FILES__.items():
 f=R/name;f.parent.mkdir(parents=True,exist_ok=True);f.write_text(content)
old=Path('/workspace/microbubble_lammps/work/lammps_particle_engine_20260915_214759')
assert '22 Jul 2025' in (old/'upstream/lammps/src/version.h').read_text()
inventory=[]
for root in [old,Path('/workspace/microbubble_lammps/work/palabos_lammps_coupling_v0_20260915_224019'),Path('/workspace/microbubble_lammps/work/passive_transport_v0_20260915_233516')]:
 for f in root.rglob('*'):
  if f.is_file() and (f.name in ['lmp','liblammps.a'] or f.name.startswith(('passive_lmp','coupled_lmp'))):inventory.append({'path':str(f),'size':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
(R/'provenance/OLD_ACTIVE_BEFORE.json').write_text(json.dumps(inventory,indent=2)+'\\n')
conf=f"""[unix_http_server]
file={W}/s.sock
[supervisord]
logfile={W}/supervisor.log
pidfile={W}/supervisor.pid
childlogdir={W}
[rpcinterface:supervisor]
supervisor.rpcinterface_factory=supervisor.rpcinterface:make_main_rpcinterface
[supervisorctl]
serverurl=unix://{W}/s.sock
[program:build]
command=/usr/bin/python3 {R}/scripts/build_worker.py
autostart=true
autorestart=false
startsecs=0
stdout_logfile={W}/build_worker.stdout
stderr_logfile={W}/build_worker.stderr
"""
(W/'supervisord.conf').write_text(conf)
subprocess.run(['supervisord','-c',str(W/'supervisord.conf')],check=True)
print(json.dumps({'state':'BUILD_LAUNCHED','work':str(W),'result':str(R),'old_inventory':inventory},indent=2))
'''.replace('__STATE__',repr(s)).replace('__FILES__',repr(files))
x=rex.run(remote,90);print(x.stdout);print(x.stderr,file=sys.stderr);raise SystemExit(x.returncode)
