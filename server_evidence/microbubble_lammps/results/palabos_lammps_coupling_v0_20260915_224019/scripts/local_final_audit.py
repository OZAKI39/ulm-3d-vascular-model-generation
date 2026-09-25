from pathlib import Path
import json,hashlib,subprocess,sys,time
r=Path('/home/lzy/projects/compre_output/palabos_lammps_coupling_v0/20260915_224019')
assert json.loads((r/'REMOTE_TO_WSL_INTEGRITY.json').read_text())['status']=='PASS'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
checks={}
for label,p in [('SONOVUE_SAMPLER',Path('/home/lzy/projects/sonovue_size_distribution_v0')),('LAMMPS_ENGINE_RESULTS',Path('/home/lzy/projects/compre_output/lammps_particle_engine/20260915_214759')),('PBS_BSA_RESULTS',Path('/home/lzy/projects/compre_output/pure_fluid_new_medium_smoke/20260915_161937'))]:
 fails=[];count=0
 for line in (p/'SHA256SUMS').read_text().splitlines():
  digest,name=line.split('  ',1);count+=1
  if sha(p/name)!=digest:fails.append(name)
 checks[label]={'checked_files':count,'mismatches':fails}
assert all(not x['mismatches'] for x in checks.values()),checks
(r/'LOCAL_FROZEN_BASELINE_AUDIT.json').write_text(json.dumps({'status':'PASS','checks':checks},indent=2)+'\n')
for args,log in [(['src/finalize_palabos_lammps_coupling.py','--root',str(r),'--output',str(r/'LOCAL_INDEPENDENT_FINALIZER.json')],'LOCAL_FINALIZER.log'),(['scripts/audit_geometry_local.py',str(r)],'LOCAL_GEOMETRY_AUDIT.log')]:
 with (r/log).open('w') as f:p=subprocess.run(['/usr/bin/python3','-B']+args,cwd=r,stdout=f,stderr=subprocess.STDOUT)
 print(log,'EXIT',p.returncode,flush=True);p.check_returncode()
print('LOCAL_INDEPENDENT_FINALIZER_AND_EXACT_GEOMETRY=PASS')
