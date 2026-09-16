from pathlib import Path
import subprocess,os,json,hashlib
p=Path('/workspace/microbubble_lammps/results/lammps2026_bubble_interaction_20260916_091244/source_correction_audit');w=Path('/workspace/lammps_migration/new_20260916_091244')
env=dict(os.environ,HWLOC_COMPONENTS='-gl',OMP_NUM_THREADS='1')
with (p/'compile.log').open('w') as f:
 subprocess.run(['g++','-O2','-std=c++17',str(p/'source_coefficient_probe.cpp'),'-o',str(w/'source_coefficient_probe')],check=True,stdout=f,stderr=subprocess.STDOUT)
subprocess.run([str(w/'source_coefficient_probe'),str(p/'pairs.txt'),str(p/'SOURCE_COEFFICIENTS.csv')],check=True)
for d in sorted(p.glob('actual_lammps_*')):
 assert not (d/'RUN_METRICS.json').exists()
 with (d/'stdout.log').open('w') as f:
  a=[str(w/'build_cpu/lmp'),'-in','in.lammps','-log','log.lammps'];v=subprocess.run(a,cwd=d,stdout=f,stderr=subprocess.STDOUT,env=env,timeout=60)
 (d/'RUN_METRICS.json').write_text(json.dumps({'returncode':v.returncode,'command':a,'steps':0},indent=2)+'\n');assert v.returncode==0,d
print('SOURCE_PROBES_COMPLETED')
