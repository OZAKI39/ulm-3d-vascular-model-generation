from pathlib import Path
import hashlib,json,shutil,os,time,subprocess
R=Path('/workspace/microbubble_lammps/results/lammps2026_bubble_interaction_20260916_091244');W=Path('/workspace/lammps_migration/new_20260916_091244')
old=Path('/workspace/microbubble_lammps/work/lammps_particle_engine_20260915_214759')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
roots={'engine':Path('/workspace/microbubble_lammps/results/lammps_particle_engine_20260915_214759'),'coupling':Path('/workspace/microbubble_lammps/results/palabos_lammps_coupling_v0_20260915_224019'),'passive':Path('/workspace/microbubble_lammps/results/passive_transport_v0_20260915_233516'),'PBS_BSA':Path('/workspace/hemocell_restore/results/new_medium_smoke_20260915_161937')}
audit={}
for label,root in roots.items():
 lines=(root/'SHA256SUMS').read_text().splitlines();bad=[]
 for line in lines:
  h,n=line.split('  ',1)
  if sha(root/n)!=h:bad.append(n)
 audit[label]={'files':len(lines),'mismatches':bad,'root':str(root),'manifest_sha256':sha(root/'SHA256SUMS')};assert not bad
expected=json.loads((roots['coupling']/'provenance/BASELINE_SOURCE_BEFORE.json').read_text());bad=[n for n,h in expected.items() if sha(old/n)!=h];assert not bad
audit['old_source_binary']={'files':len(expected),'mismatches':bad}
audit['status']='PASS';audit['epoch']=time.time();(R/'provenance/HISTORICAL_BASELINE_INTEGRITY.json').write_text(json.dumps(audit,indent=2)+'\n')
for name in ['LICENSE','COPYING']:
 p=W/'upstream/lammps'/name
 if p.exists():shutil.copy2(p,R/'provenance'/('UPSTREAM_'+name))
for n in ['pair_lubricate_poly.cpp','pair_lubricateU_poly.cpp']:
 shutil.copy2(old/'upstream/lammps/src/COLLOID'/n,R/'provenance/upstream_lubrication'/('OLD_ACTIVE_'+n))
shutil.copy2(R/'SOURCE_CORRECTION_AUDIT.json',R/'source_correction_audit/SOURCE_CORRECTION_AUDIT.json')
(W/'NOT_PROMOTED.json').write_text(json.dumps({'status':'CANDIDATE_ONLY_NOT_ACTIVE','source_gate':'FAIL','reason':'Documented unequal-radius full correction not present as required; no promotion or deletion authorized by failed gate','old_active':str(old),'phase_B':'NOT_STARTED','phase_C':'NOT_STARTED'},indent=2)+'\n')
print(json.dumps(audit,indent=2))
