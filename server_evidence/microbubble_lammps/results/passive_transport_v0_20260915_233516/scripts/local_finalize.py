from pathlib import Path
import subprocess,json,hashlib,shutil,sys
c=Path(__file__).parent;r=Path(json.loads((c/'TASK_PATHS.json').read_text())['local_result']);assert json.loads((r/'REMOTE_TO_WSL_INTEGRITY.json').read_text())['status']=='PASS'
checks={}
for label,p in [('COUPLING',Path('/home/lzy/projects/compre_output/palabos_lammps_coupling_v0/20260915_224019')),('SAMPLER',Path('/home/lzy/projects/sonovue_size_distribution_v0')),('PBS_BSA',Path('/home/lzy/projects/compre_output/pure_fluid_new_medium_smoke/20260915_161937'))]:
 out=subprocess.run(['sha256sum','-c','SHA256SUMS'],cwd=p,text=True,capture_output=True);checks[label]={'status':'PASS' if out.returncode==0 else 'FAIL','files':len(out.stdout.splitlines())};assert out.returncode==0
(r/'LOCAL_BASELINE_AUDIT.json').write_text(json.dumps(checks,indent=2)+'\n')
for n in ['local_geometry_audit.py','visualize_passive.py']:shutil.copy2(c/n,r/'scripts'/n)
shutil.copy2(c/'finalize_passive_transport_v0.py',r/'finalize_passive_transport_v0.py')
with (r/'LOCAL_FINALIZER.log').open('w') as f:out=subprocess.run(['/usr/bin/python3','-B',str(r/'finalize_passive_transport_v0.py'),'--root',str(r)],stdout=f,stderr=subprocess.STDOUT)
print('LOCAL_FULL_FINALIZER',out.returncode,flush=True);out.check_returncode()
