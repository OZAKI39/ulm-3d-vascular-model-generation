from common import *
import os,sys,time,resource,subprocess,shutil
start=time.perf_counter()
# Hard gate: accepted immutable mathematical sources and every sealed asset.
manifest=json.loads((PREV/'PAYLOAD_MANIFEST.json').read_text())
assert sha(PREV/'PAYLOAD_MANIFEST.json')=='3265476d2e6b464f0e2b771fc52731310f97bdc746e54287654725a2abff6985'
conflict=[p for p,v in manifest['files'].items() if sha(R/p)!=v['sha256']];assert not conflict,conflict
prior=json.loads((PREV/'SOURCE_AND_LIBRARY_IDENTITY.json').read_text())
assert sha(OLD/'build_owner_fix/liblammps.so')==EXPECTED
assert sha(R/'tools/lammps-native-y-wall-shared-v1-build/liblammps.so')=='d1b2f53d6b942d1ca9277a8b8cd9c195ce75ea6f616681ec21c311e4faf3b41a'
O.mkdir(parents=True,exist_ok=False);F.mkdir(parents=True,exist_ok=False);(O/'provenance').mkdir()
def cmd(c):
 z=subprocess.run(c,cwd=R,capture_output=True,text=True);return dict(code=z.returncode,stdout=z.stdout,stderr=z.stderr)
files={p:dict(sha256=sha(p),bytes=Path(p).stat().st_size) for p in prior['files']}
identity=dict(files=files,pwd=str(R),git_status=cmd(['git','status','--short']),git_branch=cmd(['git','branch','--show-current']),git_head=cmd(['git','rev-parse','HEAD']),compiler=cmd(['c++','--version']),MPI=cmd(['mpirun','--version']),python=sys.version,executable=sys.executable,affinity=sorted(os.sched_getaffinity(0)),cgroup={str(p):p.read_text().strip() for p in [Path('/sys/fs/cgroup/cpu,cpuacct/cpu.cfs_quota_us'),Path('/sys/fs/cgroup/cpu,cpuacct/cpu.cfs_period_us'),Path('/sys/fs/cgroup/memory/memory.limit_in_bytes')]},memory=cmd(['free','-b']),disk=cmd(['df','-B1','/workspace']),jobs=cmd(['ps','-eo','pid,comm,pcpu,pmem,etime','--sort=-pcpu']))
import lammps
identity['lammps_python_binding']=dict(path=lammps.__file__,sha256=sha(lammps.__file__))
dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',identity)
dump(O/'OFFLINE_MODEL_PROVENANCE.json',dict(status='PASS',verified_assets=len(manifest['files']),manifest_sha256=sha(PREV/'PAYLOAD_MANIFEST.json'),formal_model_sha256=sha(R/'analysis/softening_feasibility_v1/model.py'),definition_sha256=sha(PREV/'SOFTENING_FREE_ENERGY_DEFINITION.md'),evolution_sha256=sha(PREV/'DAMAGE_EVOLUTION_DEFINITION.md'),reference_sha256=sha(PREV/'provenance/reference.npz'),frozen_sources=json.loads((PREV/'DESIGN_FREEZE.json').read_text())['sources']))
protected={}
for root in [R/'verification',R/'results']:
 for p in root.rglob('*'):
  if p.is_file() and TASK not in p.parts:protected[str(p.relative_to(R))]=dict(sha256=sha(p),bytes=p.stat().st_size)
dump(O/'HISTORICAL_PROTECTION_BEFORE.json',protected)
for n in ['reference.npz','body.data','postbreak_same_coordinates.npz']:shutil.copyfile(PREV/'provenance'/n,O/'provenance'/n)
z=np.load(O/'provenance/reference.npz');assert z['X'].shape==(341,3) and z['bonds'].shape==(12047,2)
# Recover actual immutable original data-file bond IDs, never current loop indices.
s=(O/'provenance/body.data').read_text().split('Bonds')[1].splitlines();rows=[]
for line in s:
 try:v=list(map(int,line.split()))
 except ValueError:continue
 if len(v)==4:rows.append(v)
by_pair={tuple(sorted((r[2],r[3]))):r for r in rows};assert len(by_pair)==12047
out=['LPSVAR_REFERENCE_V1',f'341 12047 {float(z["V"][0]):.17g}']
for i,m in enumerate(z['m']):out.append(f'{i+1} {m:.17g}')
ids=[]
for b,(ij,r,h) in enumerate(zip(z['bonds'],z['r0'],z['taper'])):
 i,j=map(int,ij+1);row=by_pair[tuple(sorted((i,j)))];ids.append(row[0]);out.append(f'{row[0]} {i} {j} {row[1]} {r:.17g} {h:.17g} 0 0 0')
(O/'provenance/reference.native').write_text('\n'.join(out)+'\n');np.save(O/'provenance/original_bond_ids.npy',np.array(ids,dtype=np.int64))
freeze=dict(run_id=RUN,model='EXACT_FROZEN_ENDPOINT_STAR_VARIATIONAL_COMPLIANCE',model_provenance=json.loads((O/'OFFLINE_MODEL_PROVENANCE.json').read_text()),material=dict(K_Pa=K,G_Pa=G,rho=1000,dx=DX,horizon=DELTA,alpha=0,s_onset=S0,sc=[.18,.24]),interface='new BondStyle bpm/peri/variational plus explicit collective C API for diagnostic import, accepted commit and read-only export; old bpm/peri files unchanged',units='native atom f retains original force-density convention N/m3; exported force N and complete U J; native bond energy U/V',storage='replicated small-reference graph keyed by original bond ID + immutable endpoint tags + type + zero image; reject duplicate endpoints, periodicity, variable volumes, >4096 nodes; no production scaling claim',MPI='collect only owned coordinates by immutable atom tag with Allreduce and unique-owner assertion; identical replicated graph/history across ranks, collective commit once per monotonically increasing token; forces added only to owned atoms; no ghost/local-index history',restart='versioned native bond-style payload stores full original graph, fixed m,V,r,taper,kappa,d,failed flags,cumulative D,event and token counters; missing/wrong version rejected',commit_order='TRIAL compute never mutates history; explicit accepted commit at fixed current coordinates evaluates U_old, updates kappa/d, evaluates U_new, accumulates signed D, removes newly failed native owned bonds and their 1-2 exclusions once; next run0 rebuilds neighbors. Repeat accepted token at unchanged coordinates is no-op; different coordinates rejected.',local_solver='exact diagonal-plus-rank-one closed endpoint solve copied algebraically from frozen model.py; check 6 endpoint stars against dense solve',force_parity_relative=2e-10,force_parity_absolute_N=2e-17,energy_parity_relative=2e-10,energy_parity_absolute_J=1e-20,gradient_relative_global=3e-6,gradient_absolute_N=2e-12,gradient_steps_m=[5e-10,1e-10],history_abs=2e-13,work_balance_relative=2e-8,zero_cross_force_N=0.0,reference_paths={n:sha(PREV/(n+'.npz')) for n in ['instant_reference','soft_sc018','soft_sc024','planar_partial_sc024']},admission=['identity','undamaged','partial','gradient','two sc native controlled paths','complete separation','sorting / 2MPI migration / restart','conditional contact'],budget=dict(compile_jobs=12,compile_max_seconds=900,test_threads=1,MPI_ranks=2,main_nodes=341,formal_loading_paths=2,numerical_bug_repair_budget=1,no_V3=True,no_large_ModeI=True,no_GPU=True))
dump(O/'NATIVE_SOFTENING_DESIGN_FREEZE.json',freeze)
# Whole isolated source copy, no symlink/hardlink mutation of old code.
shutil.copytree(OLD/'source_owner_fix',O/'source',symlinks=True,ignore=shutil.ignore_patterns('.git','__pycache__'))
dump(O/'PREPARE_RESOURCES.json',dict(wall_s=time.perf_counter()-start,RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,protected_files=len(protected)))
print('PREPARED',RUN,len(protected),identity['cgroup'],flush=True)
