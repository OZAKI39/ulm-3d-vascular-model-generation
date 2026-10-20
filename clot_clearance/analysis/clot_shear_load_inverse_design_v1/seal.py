from common import *
import subprocess,difflib,shutil,datetime
assert json.loads((O/'DELIVERY_VALIDATION.json').read_text())['status']=='PASS'
assert json.loads((O/'WSL_DOCUMENT_UPDATE.json').read_text())['status']=='PASS'
visual=json.loads((O/'VISUAL_REVIEW.json').read_text());assert visual['status']=='PASS'
for relative,expected in visual['figure_sha256'].items():assert sha(R/relative)==expected,relative
baseline=json.loads((O/'BASELINE_MANIFEST.json').read_text())['files']
for path,entry in baseline.items():assert sha(Path(path))==entry['sha256'],path
freeze=json.loads((O/'SOURCE_FREEZE.json').read_text())
for relative,want in freeze['implementation_sha256'].items():assert sha(R/relative)==want,relative
for relative,want in freeze['input_sha256'].items():assert sha(O/relative)==want,relative
dump(O/'BASELINE_FINAL_IDENTITY.json',dict(status='PASS',files_verified=len(baseline),scope='Reused native library, prior sealed manifest, model/graph dependencies and copied reference evidence remain byte-identical.'))
services={}
for name in ['clot_follower_'+name+'_'+RUN for name in CASES]+['clot_follower_postprocess_'+RUN]:
 result=subprocess.run(['supervisorctl','status',name],capture_output=True,text=True).stdout;assert 'EXITED' in result,result;services[name]=result.strip()
dump(O/'SERVICE_FINAL_STATE.json',services)
deps=[LIB,R/'analysis/v3_softening_pilot_v1/native.py',R/'analysis/softening_feasibility_v1/model.py',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py',HA/'analyze.py']
dump(O/'DEPENDENCY_SHA256.json',{str(p):sha(p) for p in deps})
sources=[p for p in A.iterdir() if p.is_file() and p.suffix in ['.py','.sh','.md']];dump(O/'SOURCE_SHA256.json',{str(p.relative_to(R)):sha(p) for p in sources})
patch=''.join(''.join(difflib.unified_diff([],p.read_text().splitlines(True),fromfile='/dev/null',tofile='b/clot_clearance/'+str(p.relative_to(R)))) for p in sorted(sources));(O/'SOURCE_DIFF.patch').write_text(patch);shutil.copy2(A/'REPRODUCE.md',O/'REPRODUCE.md')
ep=json.loads((O/'ENTRYPOINT_UPDATE.json').read_text());assert ep['status']=='PASS' and sha(R/ep['path'])==ep['after_sha256']
files=[dict(path=ep['path'],kind='file',bytes=(R/ep['path']).stat().st_size,sha256=ep['after_sha256'])]
for root in [A,O,F]:
 for p in sorted(root.rglob('*')):
  if '__pycache__' in p.parts or p.name in ['DELIVERY_MANIFEST.json','LOCAL_READBACK.json']:continue
  if not p.is_file():continue
  assert '.tmp' not in p.name,p
  if p.is_symlink():
   target=p.resolve();assert target.is_relative_to(O);files.append(dict(path=str(p.relative_to(R)),kind='symlink',target=os.readlink(p),target_sha256=sha(target)))
  else:files.append(dict(path=str(p.relative_to(R)),kind='file',bytes=p.stat().st_size,sha256=sha(p)))
dump(O/'DELIVERY_MANIFEST.json',dict(status='SEALED',task=TASK,run=RUN,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),file_count=len(files),regular_file_bytes=sum(q.get('bytes',0) for q in files),files=files,server_dependencies='Reused native binary and upstream actual HDF5 remain on server; source/hash/input/accepted-state evidence is delivered. This is not a server image.'))
print('SEALED',len(files),sum(q.get('bytes',0) for q in files),sha(O/'DELIVERY_MANIFEST.json'))
