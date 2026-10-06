from pathlib import Path
import subprocess,json,hashlib,os
DST=Path('/home/lzy/projects/computation_examples/_archives/github_sync_brava_20261006');META=DST/'sync_metadata/brava_p2_20261006';BASE='/workspace/brava_flow_roi_18mlmin_20260928'
code=r'''
from pathlib import Path
import hashlib,json
b=Path('/workspace/brava_flow_roi_18mlmin_20260928');selected=[];excluded=[]
rounds=['native_p2p1_review_20261006','common_curved_geometry_review_20261006','target_bc_pipe_review_20261006','p2_linear_algebra_review_20261006']
for name in rounds:
 d=b/'population_v1/analysis'/name
 selected += [p for p in (d/'monitoring').glob('*/*') if p.name in {'job.log','EXECUTION.json'}]
 for p in (d/'runtime').rglob('*'):
  if p.is_symlink() or not p.is_file():continue
  if p.stat().st_size>=25_000_000:excluded.append(dict(path=str(p.relative_to(b)),bytes=p.stat().st_size,reason='Generated runtime matrix/build/duplicate input; retained on server'))
# Frozen indispensable inputs and byte-identical small run configuration.
selected += [b/'population_v1/analysis/target_bc_pipe_review_20261006/runtime/PIPE_CURVED_ELEVATION.npz']
p1=b/'production_workflow/manufactured_pipe/R1_GMSH_V2/cfd_gpu/run/8-procs'
selected += [p1/'stFile_010.bin',p1/'result_010.vtu']
for name in ['p2_linear_algebra_review_20261006','p2_linear_algebra_rerun_20261006']:
 d=b/'population_v1/analysis'/name/'runtime'
 for c in ['p1_baseline','p1_asm0','p1_schur','p1_candidate_fixed_observer','p2_system']:
  for n in ['solver.xml','PETSC_OPTIONS.txt','inlet_exact_velocity.dat','outlet_exact_traction.vtp']:
   p=d/c/'run'/n
   if p.exists():selected.append(p)
files=[]
for p in sorted(set(selected)):
 with p.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
 files.append(dict(path=str(p.relative_to(b)),bytes=p.stat().st_size,sha256=h))
print(json.dumps(dict(files=files,excluded_large_runtime=excluded)))
'''
r=subprocess.run(['ssh','vast4090','python3','-'],input=code,text=True,capture_output=True,check=True);manifest=json.loads(r.stdout);data=manifest['files'];(META/'SERVER_EVIDENCE.json').write_text(json.dumps(manifest,indent=2)+'\n')
paths=Path('/tmp/brava_github_sync_20261006/server_files.txt');paths.write_text(''.join(x['path']+'\n' for x in data))
subprocess.run(['rsync','-a','--files-from='+str(paths),'vast4090:'+BASE+'/',str(DST/'brava_flow_roi_18mlmin')+'/'],check=True)
for x in data:
 p=DST/'brava_flow_roi_18mlmin'/x['path']
 with p.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
 assert h==x['sha256'],x['path']
p=META/'INCLUDED_FILES.json';v=json.loads(p.read_text());existing={x['path'] for x in v['files']}
for x in data:
 path='brava_flow_roi_18mlmin/'+x['path']
 if path not in existing:v['files'].append(dict(path=path,source=BASE+'/'+x['path'],sha256=x['sha256'],bytes=x['bytes'],kind='file',role='Frozen required input / solver configuration / supervised log and receipt'))
p.write_text(json.dumps(v,indent=2)+'\n')
print('Server evidence files',len(data),'bytes',sum(x['bytes'] for x in data),'excluded generated runtime >=25MB',len(manifest['excluded_large_runtime']))
