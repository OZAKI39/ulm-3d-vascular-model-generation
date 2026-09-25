from pathlib import Path
import subprocess,json,os,hashlib,re
R=Path('/workspace/microbubble_lammps/results/lammps2026_bubble_interaction_20260916_091244');W=Path('/workspace/lammps_migration/new_20260916_091244');env=dict(os.environ,HWLOC_COMPONENTS='-gl')
items={};links={}
for kind in ['cpu','gpu']:
 text=(R/'provenance'/(kind+'_lmp_help.txt')).read_text();tokens=set(text.split());cache=(R/'provenance'/(kind+'_CMakeCache.txt')).read_text();items[kind]={'sphere':'sphere' in tokens,'nve/sphere':'nve/sphere' in tokens,'nve/sphere/kk':'nve/sphere/kk' in tokens,'nve/noforce':'nve/noforce' in tokens,'lubricate/poly':'lubricate/poly' in tokens,'lubricateU/poly':'lubricateU/poly' in tokens,'lubricate/poly/kk':'lubricate/poly/kk' in tokens,'lubricateU/poly/kk':'lubricateU/poly/kk' in tokens,'COLLOID':'COLLOID' in tokens,'GRANULAR':'GRANULAR' in tokens,'MPI':'BUILD_MPI:BOOL=ON' in cache,'KOKKOS':'PKG_KOKKOS:BOOL=ON' in cache,'CUDA':'Kokkos_ENABLE_CUDA:BOOL=ON' in cache,'LATBOLTZ':'PKG_LATBOLTZ:BOOL=ON' in cache,'binary_sha256':hashlib.sha256((W/('build_'+kind)/'lmp').read_bytes()).hexdigest()}
 exe=W/('build_'+kind)/'lmp';links[str(exe)]=subprocess.check_output(['ldd',str(exe)],text=True,env=env)
for stage,prefix in [('coupling','coupled'),('passive','passive')]:
 for kind in ['cpu','gpu']:
  exe=W/stage/'build'/(prefix+'_lmp_'+kind)
  if exe.exists():links[str(exe)]=subprocess.check_output(['ldd',str(exe)],text=True,env=env)
required=['sphere','nve/sphere','nve/noforce','lubricate/poly','lubricateU/poly','COLLOID','GRANULAR','MPI']
ok=all(items[k][t] for k in items for t in required) and items['gpu']['nve/sphere/kk'] and items['gpu']['CUDA'] and not any(v['LATBOLTZ'] for v in items.values())
(R/'provenance/CAPABILITY_AUDIT.json').write_text(json.dumps({'status':'PASS' if ok else 'FAIL','builds':items,'lubricateU_poly_accelerator_support':'NO_KOKKOS_SUFFIX_IN_FROZEN_BUILD; CPU implementation only','GPU_PERFORMANCE_READY':'NO'},indent=2)+'\n')
(R/'provenance/CANDIDATE_LINKAGE_AUDIT.json').write_text(json.dumps({'candidate_only':True,'linkage':links,'candidate_dynamic_old_lammps_references':sum('lammps_particle_engine_20260915_214759' in text for text in links.values()),'static_lammps_linkage':'build command + archive SHA recorded; not visible in ldd','PATH_unchanged':True,'which_lmp':subprocess.run(['bash','-c','command -v lmp'],capture_output=True,text=True,env=env).stdout.strip() or 'NO_EXISTING_PATH_BINDING'},indent=2)+'\n')
print(json.dumps({'status':'PASS'if ok else 'FAIL','builds':items},indent=2))
