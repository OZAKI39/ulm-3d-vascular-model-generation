from pathlib import Path
import subprocess,json,os,shutil
R=Path('/workspace/microbubble_lammps/results/stable_rigid_sphere_near_field_v0_20260916_101617');W=Path('/workspace/microbubble_lammps/work/stable_rigid_sphere_near_field_v0_20260916_101617');B='/workspace/microbubble_lammps/work/lammps_particle_engine_20260915_214759'
env=dict(os.environ,HWLOC_COMPONENTS='-gl',OMP_NUM_THREADS='1')
for label,cmd in [('CONFIGURE',['cmake','-S',str(R),'-B',str(W/'build'),'-G','Ninja','-DCMAKE_BUILD_TYPE=Release','-DENGINE_ROOT='+B]),('BUILD',['cmake','--build',str(W/'build'),'-j4'])]:
 with open(R/'logs'/f'{label}.log','w') as f:r=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
 print(label,r.returncode,flush=True)
 if r.returncode:print((R/'logs'/f'{label}.log').read_text()[-5000:]);raise SystemExit(r.returncode)
(R/'build_math').mkdir(exist_ok=True);shutil.copy2(W/'build/librigid_math.so',R/'build_math/librigid_math.so')
print('READY',flush=True)
