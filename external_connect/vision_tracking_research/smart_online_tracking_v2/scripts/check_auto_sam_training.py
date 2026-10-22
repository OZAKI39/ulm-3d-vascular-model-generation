"""Actual checkpoint load/revert and early cancellation; never promotes a test model."""
from pathlib import Path
import sys,json,copy,subprocess
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from sam_auto_training import read_json,atomic_json,RUNTIME,gpu_slot
from sam_worker import SamWorker
OUT=ROOT/'evidence/auto_sam_training_20261016';OUT.mkdir(exist_ok=True)
job=Path(read_json(RUNTIME/'last_attempt.json')['job']);result=read_json(job/'result.json')
req=read_json(ROOT/'evidence/annotation_repair_v4_20261016/reported_job/request.json')
req.update(mode='point_instances',object_class='SINGLE_BUBBLE',software_test=True,use_validated_adaptation=False)
req['positive_yx']=req['positive_yx'][:1]
worker=SamWorker();report={'software_probe_only':True,'human_labels_modified':False,'candidate_promoted':False,'checks':{}}
def probe(name,adaptation):
 folder=OUT/name;folder.mkdir(exist_ok=True);r={**req,'output':str(folder/'candidate.npz'),'token':name}
 atomic_json(folder/'request.json',r)
 # Explicit in-process test adapter bypasses registry resolution only in this harness.
 # No production active_model.json is created or changed.
 worker.in_job=True;worker.adaptation=adaptation
 with gpu_slot(interactive=True):worker.run(r)
 worker.in_job=False
 with np.load(r['output']) as z:return z['raw_target_000'].copy(),json.loads(str(z['details']))
base,meta=probe('original_probe',{})
adapt={'id':'SOFTWARE_CANDIDATE_LOAD_PROBE','checkpoint':result['checkpoint'],'checkpoint_sha256':result['checkpoint_sha256']}
candidate,meta=probe('candidate_probe',adapt)
report['checks']['actual_candidate_loaded']=meta['model_id']==adapt['id'] and meta['adapted_decoder']
restored,meta=probe('restored_probe',{})
report['checks']['original_restored_exactly']=np.array_equal(base,restored) and meta['model_id']=='ORIGINAL_SAM_VIT_B'
broken,meta=probe('corrupt_checksum_probe',{**adapt,'id':'CORRUPT_CHECKSUM_TEST','checkpoint_sha256':'wrong'})
report['checks']['corrupt_checkpoint_falls_back_to_original']=np.array_equal(base,broken) and bool(meta['adaptation_error']) and not meta['adapted_decoder']
report['candidate_changed_predicted_pixels']=int(np.count_nonzero(base!=candidate))
report['checks']['production_registry_still_original']=not read_json(RUNTIME/'active_model.json',{}).get('admitted',False)
# Pre-cancelled real worker must not allocate/train or create a checkpoint.
cancel=OUT/'cancel_probe';cancel.mkdir(exist_ok=True);atomic_json(cancel/'job.json',read_json(job/'job.json'));(cancel/'CANCEL').touch()
with (cancel/'worker.log').open('w') as log:
 proc=subprocess.run([sys.executable,'-B',str(ROOT/'scripts/sam_auto_train.py'),str(cancel/'job.json')],stdout=log,stderr=subprocess.STDOUT,timeout=40)
c=read_json(cancel/'result.json');report['checks']['cancel_without_training_or_checkpoint']=c['status']=='CANCELLED' and c['steps']==0 and not (cancel/'decoder.pt').exists()
report['PASS']=all(report['checks'].values());atomic_json(OUT/'worker_checks.json',report);print(json.dumps(report),flush=True)
if not report['PASS']:raise SystemExit(1)
