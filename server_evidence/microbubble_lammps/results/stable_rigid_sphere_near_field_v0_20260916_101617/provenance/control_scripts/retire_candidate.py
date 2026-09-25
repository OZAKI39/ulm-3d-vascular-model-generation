from pathlib import Path
import hashlib,json,shutil,os,subprocess
R=Path('/workspace/microbubble_lammps/results/stable_rigid_sphere_near_field_v0_20260916_101617')
m=json.loads((R/'LAMMPS_2SEP2026_CANDIDATE_RETIREMENT_MANIFEST.json').read_text());C=Path(m['retirement_root']);assert str(C)=='/workspace/lammps_migration/new_20260916_091244' and m['status']=='READY_FOR_RETIREMENT'
# Deletion is explicitly authorized by task sections 2 and 62. Manifest already persisted.
shutil.rmtree(C)
P=Path(m['preserved_failure_evidence']);r=subprocess.run(['sha256sum','--check','--quiet','SHA256SUMS'],cwd=P,capture_output=True,text=True);assert r.returncode==0
B=Path('/workspace/microbubble_lammps/work/lammps_particle_engine_20260915_214759');active=Path('/workspace/lammps_active')
assert not active.exists() and not active.is_symlink();active.symlink_to(B,target_is_directory=True)
audit={'status':'PASS','active_release':'22Jul2025 Update 6','active_commit':'9c5ab448c78a14fd534619622162ba418d6a1fb1','ACTIVE_LAMMPS_COUNT':1,'LAMMPS_VERSION_COEXISTENCE':'NO','active_symlink':str(active),'active_target':str(active.resolve()),'candidate_source_present':C.exists(),'candidate_build_present':C.exists(),'candidate_binary_present':C.exists(),'candidate_failure_evidence_preserved':'YES','previous_evidence_sha256_verification':'PASS','custom_binaries':[]}
for name in ['palabos_lammps_coupling_v0_20260915_224019','passive_transport_v0_20260915_233516']:
 d=Path('/workspace/microbubble_lammps/work')/name/'build'
 for p in sorted(d.iterdir()):
  if p.is_file() and (os.access(p,os.X_OK) or p.suffix=='.so'):
   ld=subprocess.run(['ldd',str(p)],capture_output=True,text=True); assert str(C) not in ld.stdout
   audit['custom_binaries'].append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'ldd':ld.stdout,'candidate_linkage':False})
(R/'LAMMPS_ACTIVE_INSTALL_AUDIT.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps({k:v for k,v in audit.items() if k!='custom_binaries'},indent=2))
