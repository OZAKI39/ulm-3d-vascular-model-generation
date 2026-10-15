from common import *
import shutil,time
shutil.copy2(Path(__file__),O/'reproduction_code/manifest.py')
skip={'ARTIFACT_MANIFEST.json','WSL_DELIVERY_VERIFICATION.json','SERVER_DELIVERY_PROOF.json'};files={}
for base in [O,F,R/'analysis/mode_i_v1']:
 for p in sorted(base.rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts and p.name not in skip:files[str(p.relative_to(R))]={'sha256':sha(p),'bytes':p.stat().st_size}
manifest={'task':C['task'],'run_id':C['run'],'created_UTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'root':str(R),'files':files,'count':len(files),'total_bytes':sum(v['bytes'] for v in files.values()),'self_and_delivery_proofs_excluded_to_avoid_hash_cycles':sorted(skip),'historical_data_not_duplicated_except_declared_reference_and_native_input':True};dump(O/'ARTIFACT_MANIFEST.json',manifest)
print(json.dumps({'files':len(files),'bytes':manifest['total_bytes'],'manifest_sha256':sha(O/'ARTIFACT_MANIFEST.json')}))
