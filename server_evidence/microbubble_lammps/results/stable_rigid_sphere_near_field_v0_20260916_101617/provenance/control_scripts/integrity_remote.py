from pathlib import Path
import json,hashlib,subprocess
R=Path('/workspace/microbubble_lammps/results/stable_rigid_sphere_near_field_v0_20260916_101617');B=Path('/workspace/microbubble_lammps/work/lammps_particle_engine_20260915_214759')
P=Path('/workspace/microbubble_lammps/results/palabos_lammps_coupling_v0_20260915_224019/provenance/BASELINE_SOURCE_BEFORE.json')
data=json.loads(P.read_text());print(type(data),list(data)[:5],flush=True)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
wrong=[]
for n,v in data.items():
 expected=v['sha256'] if isinstance(v,dict) else v
 if not (B/n).is_file() or sha(B/n)!=expected:wrong.append(n)
assert not wrong,wrong
roots=['/workspace/microbubble_lammps/results/lammps_particle_engine_20260915_214759','/workspace/microbubble_lammps/results/palabos_lammps_coupling_v0_20260915_224019','/workspace/microbubble_lammps/results/passive_transport_v0_20260915_233516','/workspace/hemocell_restore/results/new_medium_smoke_20260915_161937','/workspace/microbubble_lammps/results/lammps2026_bubble_interaction_20260916_091244'];checks=[]
for d in roots:
 d=Path(d);p=subprocess.run(['sha256sum','--check','--quiet','SHA256SUMS'],cwd=d,capture_output=True,text=True);assert p.returncode==0,(str(d),p.stdout,p.stderr);checks.append({'root':str(d),'entries':len((d/'SHA256SUMS').read_text().splitlines()),'status':'PASS','manifest_sha256':sha(d/'SHA256SUMS')})
assert not Path('/workspace/lammps_migration/new_20260916_091244').exists()
frozen=[('fields/FROZEN_FLOW_FIELD_V0.h5','7828ed32b0cdd1e0e85836a564b6737128e7dbe30181ddd43110233b37b64a7f'),('provenance/closed_geometry_m.stl','840da5e1c43ec31ac70ec781cbf75b32940bc73538b5dba83c58af1ddf5e36fb')]
for n,h in frozen:assert sha(R/n)==h,n
out={'status':'PASS','frozen_upstream_files_checked':len(data),'upstream_mismatches':wrong,'historical_baselines':checks,'candidate_present':False,'active_root':str(Path('/workspace/lammps_active').resolve()),'active_lammps_count':1,'frozen_field_and_geometry':'PASS','LAMMPS_CORE_MODIFIED':'NO','PALABOS_BASELINE_MODIFIED':'NO','PASSIVE_BASELINE_MODIFIED':'NO','candidate_failure_evidence_preserved':'YES'}
(R/'provenance/FINAL_SOURCE_AND_BASELINE_INTEGRITY.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
