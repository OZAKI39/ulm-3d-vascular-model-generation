from common import *
import ast,time
def main():
 start=time.perf_counter();status=json.loads((O/'FINAL_STATUS.json').read_text());assert status['delivery_status']=='READY_FOR_WSL_AND_SEAL'
 assert json.loads((F/'FIGURE_QA.json').read_text())['status']=='PASS'
 assert json.loads((O/'WSL_DOCUMENT_UPDATE.json').read_text())['status']=='PASS'
 verified={}
 for base in [B,PREV]:
  m=json.loads((base/'PAYLOAD_MANIFEST.json').read_text())
  for name,v in m['files'].items():assert sha(R/name)==v['sha256'],name
  verified[str(base.relative_to(R))]=len(m['files'])
 for name,h in json.loads((O/'PROTECTED_INPUTS.json').read_text()).items():assert sha(name)==h,name
 freeze=json.loads((O/'LOAD_DIRECTION_FREEZE.json').read_text());assert sha(O/'LOAD_DIRECTION_FREEZE.json')==(O/'LOAD_DIRECTION_FREEZE.sha256').read_text().strip()
 for name,h in freeze['driver_sha256'].items():assert sha(A/name)==h,name
 for name,h in freeze['input_sha256'].items():assert sha(O/'provenance'/name)==h,name
 assert sha(O/'REFERENCE_NORMAL_FIELD.npz')==freeze['reference_normal_sha256'];assert sha(LIB)==freeze['library_sha256']
 for p in A.glob('*.py'):ast.parse(p.read_text())
 source={str(p.relative_to(R)):sha(p) for p in A.iterdir() if p.is_file()};identity=json.loads((O/'SOURCE_AND_LIBRARY_IDENTITY.json').read_text());identity['final_local_sources']=source;identity['unchanged_historical_manifests']=verified;dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',identity)
 required=['REVIEW_V3_CHUNK_FRAGMENTATION_LOAD_DIRECTION_CONTROL_V1.md','LOAD_DIRECTION_FREEZE.json','REFERENCE_NORMAL_FIELD.npz','NORMAL_MAPPING_AUDIT.json','SOURCE_AND_LIBRARY_IDENTITY.json','EARLY_WINDOW_COMPARISON.json','V3_DIRECTION_CONTROL_DYNAMIC_HISTORY.csv','V3_DIRECTION_CONTROL_TOPOLOGY.json','CHUNK_COHESION_AUDIT.json','BULK_INTERFACE_FAILURE_AUDIT.json','DYNAMIC_ENERGY_LEDGER.json','HISTORICAL_DIRECTION_COMPARISON.json','RUN_RESOURCES.json','REPRODUCE.txt','OPEN_RESULTS.html']
 for name in required:assert (O/name).is_file(),name
 dump(O/'FINAL_INTEGRITY_CHECK.json',dict(status='PASS',historical_manifests_unchanged=verified,frozen_direction_and_inputs_unchanged=True,native_library_unchanged=True,required_outputs_present=required,wall_s=time.perf_counter()-start))
 status['delivery_status']='COMPLETE';dump(O/'FINAL_STATUS.json',status)
 payload={}
 for folder in [A,O,F]:
  for p in sorted(folder.rglob('*')):
   if not p.is_file() or '__pycache__' in p.parts or p.name in ['PAYLOAD_MANIFEST.json','LOCAL_READBACK.json'] or p.name.endswith('.tmp'):continue
   payload[str(p.relative_to(R))]=dict(bytes=p.stat().st_size,sha256=sha(p))
 dump(O/'PAYLOAD_MANIFEST.json',dict(status='SEALED',task=TASK,run=RUN,files=payload,total_files=len(payload),total_bytes=sum(x['bytes'] for x in payload.values()),scope='Isolated Python changes, frozen direction field, actual accepted native evidence, checkpoint/rollback logs, diagnostics, PNG/SVG, actual WSL document backup/diff/readback; original native library and old raw trajectories referenced by SHA, not duplicated.'))
 print('SEALED',len(payload),'files',sum(x['bytes'] for x in payload.values()),flush=True)
if __name__=='__main__':main()
