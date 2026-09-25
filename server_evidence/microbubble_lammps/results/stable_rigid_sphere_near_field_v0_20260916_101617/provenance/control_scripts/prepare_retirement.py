from pathlib import Path
import hashlib,json,os,subprocess,shutil,time
R=Path('/workspace/microbubble_lammps/results/stable_rigid_sphere_near_field_v0_20260916_101617');W=Path('/workspace/microbubble_lammps/work/stable_rigid_sphere_near_field_v0_20260916_101617')
for p in [R,W,R/'provenance',R/'src',R/'contracts',R/'validation',R/'theory',R/'scripts',R/'logs']:p.mkdir(parents=True,exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(name,d): (R/name).write_text(json.dumps(d,indent=2)+'\n')
B=Path('/workspace/microbubble_lammps/work/lammps_particle_engine_20260915_214759')
expect={'build_cpu/lmp':'188ebc166c967cb004af3df6e75b1903e4f5df7123eca46264899831a845bee4','build_gpu/lmp':'82dd028a0fd648b5e2355c6b955a367d7bfac4847b0ed69640e5d101fed04f65'}
actual={k:sha(B/k) for k in expect};assert actual==expect,actual
P=Path('/workspace/microbubble_lammps/results/lammps2026_bubble_interaction_20260916_091244')
p=subprocess.run(['sha256sum','--check','--quiet','SHA256SUMS'],cwd=P,capture_output=True,text=True);assert p.returncode==0,p.stdout+p.stderr
save('provenance/PREVIOUS_FAILURE_EVIDENCE_PRE_RETIREMENT.json',{'path':str(P),'sha256sums':sha(P/'SHA256SUMS'),'entries':len((P/'SHA256SUMS').read_text().splitlines()),'verification':'PASS','stable_binaries':actual})
C=Path('/workspace/lammps_migration/new_20260916_091244');assert C.is_dir() and not C.is_symlink()
active=[]
for p in Path('/proc').glob('[0-9]*'):
 try:
  exe=os.readlink(p/'exe');cmd=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
  if str(C) in exe:active.append({'pid':p.name,'exe':exe,'cmd':cmd})
 except (OSError,PermissionError):pass
assert not active,active
files=[]
for f in sorted(C.rglob('*')):
 if f.is_symlink():files.append({'path':str(f),'kind':'symlink','target':os.readlink(f)})
 elif f.is_file():files.append({'path':str(f),'kind':'file','size_bytes':f.stat().st_size,'sha256':sha(f)})
# Preserve candidate supervisor/build/run logs that may extend beyond the previous sealed report.
logs=[]
for f in C.iterdir():
 if f.is_file() and f.suffix in ('.log','.json','.conf'):
  target=R/'provenance'/'retired_candidate_control'/f.name;target.parent.mkdir(exist_ok=True);shutil.copy2(f,target);logs.append(str(target.relative_to(R)))
manifest={'status':'READY_FOR_RETIREMENT','candidate_release':'2Sep2026','candidate_commit':'d71abe6102c44577442ba7f03b7378a83166b9fd','retirement_root':str(C),'source_path':str(C/'upstream/lammps'),'build_paths':[str(C/'build_cpu'),str(C/'build_gpu'),str(C/'coupling/build'),str(C/'passive/build')],'binary_paths':[str(C/'build_cpu/lmp'),str(C/'build_gpu/lmp')],'files':files,'file_count':sum(x['kind']=='file' for x in files),'total_bytes':sum(x.get('size_bytes',0) for x in files),'preserved_failure_evidence':str(P),'additional_preserved_logs':logs,'candidate_running_processes':active}
save('LAMMPS_2SEP2026_CANDIDATE_RETIREMENT_MANIFEST.json',manifest)
for n in ['pair_lubricate_poly.cpp','pair_lubricateU_poly.cpp','pair_lubricate.cpp','pair_lubricateU.cpp']:
 shutil.copy2(B/'upstream/lammps/src/COLLOID'/n,R/'provenance'/n)
for n in ['pair_lubricate.rst','pair_lubricateU.rst']:
 shutil.copy2(B/'upstream/lammps/doc/src'/n,R/'provenance'/n)
print(json.dumps({k:v for k,v in manifest.items() if k!='files'},indent=2),flush=True)
