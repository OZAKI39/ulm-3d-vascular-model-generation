from common import *
import datetime
assert json.loads((O/'DELIVERY_VALIDATION.json').read_text())['status']=='PASS'
assert json.loads((O/'VIDEO_GENERATION.json').read_text())['status']=='PASS'
for name in ['WSL_DOCUMENT_UPDATE.json','ENTRYPOINT_UPDATE.json']:assert (O/name).exists()
freeze=json.loads((O/'SOURCE_FREEZE.json').read_text())
for p,h in freeze['implementation_sha256'].items():assert sha(R/p)==h,p
for p,h in freeze['input_sha256'].items():assert sha(O/p)==h,p
assert sha(OLDLIB)==freeze['historical_library_sha256']
dump(O/'POSTPROCESS_SOURCE_FREEZE.json',dict(files={str(p.relative_to(R)):sha(p) for p in A.glob('*.py')},runtime_freeze_sha256=sha(O/'SOURCE_FREEZE.json')))
excludes={'DELIVERY_MANIFEST.json','LOCAL_READBACK.json'};files={};dirs={}
for root in [A,O,F]:
 for p in sorted(root.rglob('*')):
  if p.is_symlink() and p.is_dir():dirs[str(p.relative_to(R))]=os.readlink(p);continue
  if not p.is_file() or '__pycache__' in p.parts or p.name in excludes or p.name.endswith('.tmp'):continue
  item=dict(bytes=p.stat().st_size,sha256=sha(p))
  if p.is_symlink():item['symlink_target']=os.readlink(p)
  files[str(p.relative_to(R))]=item
dump(O/'DELIVERY_MANIFEST.json',dict(task=TASK,run=RUN,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=files,directory_symlinks=dirs,total_bytes=sum(q['bytes'] for q in files.values()),file_count=len(files),dependencies='New isolated native library included; unchanged upstream link objects, older libraries and shared scientific dependencies remain on server with hashes in BUILD_AUDIT and SOURCE_FREEZE.'))
print('SEALED',len(files),sum(q['bytes'] for q in files.values()),sha(O/'DELIVERY_MANIFEST.json'),flush=True)
