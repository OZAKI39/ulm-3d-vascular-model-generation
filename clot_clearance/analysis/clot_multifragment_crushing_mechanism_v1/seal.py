from common import *
import subprocess,difflib,shutil,datetime,platform,scipy
assert json.loads((O/'DELIVERY_VALIDATION.json').read_text())['status']=='PASS'
assert json.loads((O/'WSL_DOCUMENT_UPDATE.json').read_text())['status']=='PASS'
for name in ['clot_crushing_uniform_20261017','clot_crushing_heterogeneous_20261017']:
 result=subprocess.run(['supervisorctl','status',name],capture_output=True,text=True).stdout;assert 'EXITED' in result,result
src=R/'verification/lps_native_softening_event_resolved_integration_v1/20261013T235945Z/source/src'
deps=[LIB,OLDLIB,R/'analysis/v3_softening_pilot_v1/native.py',R/'analysis/softening_feasibility_v1/model.py',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py',src/'fix_indent.cpp',src/'version.h']
dump(O/'DEPENDENCY_SHA256.json',{str(p):sha(p) for p in deps});dump(O/'ENVIRONMENT.json',dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,lammps_version_header=(src/'version.h').read_text(),new_dependencies=[],native_build_audit='BUILD_AUDIT.json'))
sources=[p for p in A.iterdir() if p.is_file() and p.suffix in ['.py','.sh','.md','.cpp','.h']];dump(O/'SOURCE_SHA256.json',{str(p.relative_to(R)):sha(p) for p in sources});patch=''
for p in sorted(sources):
 previous=R/'analysis/v3_softening_pilot_v1'/p.name;old=previous.read_text().splitlines(True) if p.suffix in ['.cpp','.h'] else [];patch+=''.join(difflib.unified_diff(old,p.read_text().splitlines(True),fromfile='baseline/'+p.name if old else '/dev/null',tofile='b/clot_clearance/'+str(p.relative_to(R))))
(O/'SOURCE_DIFF.patch').write_text(patch);shutil.copy2(A/'REPRODUCE.md',O/'REPRODUCE.md')
excluded=[];files=[];cache={}
for root in [A,O,F]:
 for p in sorted(root.rglob('*')):
  if '__pycache__' in p.parts or p.name in ['DELIVERY_MANIFEST.json','LOCAL_READBACK.json']:continue
  if not p.is_file():continue
  if p.parent==O/'build' and (p.name.endswith('.o') or p.name.startswith('liblammps.so')):
   excluded.append(dict(path=str(p.relative_to(R)),bytes=p.stat().st_size,sha256=sha(p),reason='Server-only reproducible build object/library; source, command and SHA delivered'));continue
  assert '.tmp' not in p.name,p
  if p.is_symlink():
   target=p.resolve();assert target.is_relative_to(R);cache.setdefault(target,sha(target));files.append(dict(path=str(p.relative_to(R)),kind='symlink',target=os.readlink(p),target_sha256=cache[target]))
  else:files.append(dict(path=str(p.relative_to(R)),kind='file',bytes=p.stat().st_size,sha256=sha(p)))
dump(O/'DELIVERY_MANIFEST.json',dict(status='SEALED',task=TASK,run=RUN,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),file_count=len(files),regular_file_bytes=sum(q.get('bytes',0) for q in files),server_only_builds=excluded,files=files));print('SEALED',len(files),sum(q.get('bytes',0) for q in files),sha(O/'DELIVERY_MANIFEST.json'))
