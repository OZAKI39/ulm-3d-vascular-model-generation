"""Seal only this round after simulations, report, WSL update and validation."""
from common import *
import difflib,subprocess,datetime,shutil

def main():
 assert json.loads((O/'DELIVERY_VALIDATION.json').read_text())['status']=='PASS'
 for service in ['clot_mechanical_bracketing_20261016','clot_mechanical_supplement_20261016','clot_mechanical_continuation_20261016','clot_mechanical_midpoint_20261016']:
  q=subprocess.run(['supervisorctl','status',service],capture_output=True,text=True)
  assert 'RUNNING' not in q.stdout and 'STARTING' not in q.stdout,q.stdout
 deps=[LIB,R/'analysis/v3_softening_pilot_v1/native.py',R/'analysis/v3_softening_pilot_v1/bond_bpm_peri_variational.cpp',R/'analysis/v3_softening_pilot_v1/fix_lps_variational_dynamic.cpp',R/'analysis/softening_feasibility_v1/model.py',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py']
 dump(O/'DEPENDENCY_SHA256.json',{str(p):sha(p) for p in deps})
 sources=[p for p in A.iterdir() if p.is_file() and p.suffix in ['.py','.sh','.md']]
 dump(O/'SOURCE_SHA256.json',{str(p.relative_to(R)):sha(p) for p in sources})
 patch=''.join(''.join(difflib.unified_diff([],p.read_text().splitlines(True),fromfile='/dev/null',tofile='b/clot_clearance/'+str(p.relative_to(R)))) for p in sorted(sources))
 (O/'SOURCE_DIFF.patch').write_text(patch);shutil.copy2(A/'REPRODUCE.md',O/'REPRODUCE.md')
 book=json.loads((O/'LOAD_CASES.json').read_text());entries=[]
 for req in sorted((O/'requests').glob('*.json')):
  q=json.loads(req.read_text());case=Path(q['case']);summary=case/'SUMMARY.json';failed=case/'FAILED_PREFLIGHT.json';assert summary.exists() or failed.exists(),case
  result=json.loads((summary if summary.exists() else failed).read_text());entries.append(dict(request=str(req.relative_to(R)),request_sha256=sha(req),case=str(case.relative_to(R)),parameters=q['params'],status=result['status'],kind='CONTINUATION' if q['params'].get('replay_source') else 'PHYSICAL_CASE_ATTEMPT'))
 book['executed_request_index']=entries;book['completed_unique_physical_conditions']=len({(q['parameters']['mode'],q['parameters']['peak_pa']) for q in entries if q['status']=='COMPLETED'});dump(O/'LOAD_CASES.json',book)
 manifest=O/'DELIVERY_MANIFEST.json';excluded={'DELIVERY_MANIFEST.json','LOCAL_READBACK.json'};files=[];cache={}
 for root in [A,O,F]:
  for p in sorted(root.rglob('*')):
   if '__pycache__' in p.parts or p.name in excluded:continue
   if p.is_symlink():
    target=p.resolve();assert target.is_relative_to(R) and target.is_file(),p
    if target not in cache:cache[target]=sha(target)
    files.append(dict(path=str(p.relative_to(R)),kind='symlink',target=os.readlink(p),target_sha256=cache[target]));continue
   if not p.is_file():continue
   assert '.tmp' not in p.name,p
   files.append(dict(path=str(p.relative_to(R)),kind='file',bytes=p.stat().st_size,sha256=sha(p)))
 dump(manifest,dict(status='SEALED',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),task=TASK,run=RUN,file_count=len(files),regular_file_bytes=sum(q.get('bytes',0) for q in files),native_library_sha256=sha(LIB),files=files))
 print('SEALED',len(files),sum(q.get('bytes',0) for q in files),sha(manifest))
if __name__=='__main__':main()
