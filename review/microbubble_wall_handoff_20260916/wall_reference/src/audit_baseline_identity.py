"""Read-only whole-manifest baseline verification; never import baseline code."""
from pathlib import Path
import hashlib,json,sys,time
R=Path(__file__).resolve().parents[1]
roots=[Path('/home/lzy/projects/sonovue_size_distribution_v0')]+[Path('/home/lzy/projects/compre_output')/p for p in ['pure_fluid_new_medium_smoke/20260915_161937','stable_rigid_sphere_near_field_v0/20260916_101617','palabos_lammps_coupling_v0/20260915_224019','passive_transport_v0/20260915_233516','lammps_particle_engine/20260915_214759']]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
results=[]
for root in roots:
 manifest=root/'SHA256SUMS';bad=[];n=0
 for line in manifest.read_text().splitlines():
  h,name=line.split(None,1);name=name.lstrip('* ');p=root/name;n+=1
  if not p.is_file() or sha(p)!=h:bad.append(name)
 results.append({'root':str(root),'manifest_sha256':sha(manifest),'checked_files':n,'mismatches':bad})
 print(root.name,n,len(bad),flush=True)
status='PASS' if all(not r['mismatches'] for r in results) else 'FAIL'
(R/'provenance'/('BASELINE_IDENTITY_'+sys.argv[1]+'.json')).write_text(json.dumps({'status':status,'baselines':results,'remote_access':'NOT_USED; this task runs entirely locally'},indent=2)+'\n');assert status=='PASS'
