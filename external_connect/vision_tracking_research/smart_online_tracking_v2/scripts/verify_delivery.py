"""Verify actual delivery, protected inputs and immutable SHA snapshot."""
from pathlib import Path
import os,sys,json,hashlib,datetime,subprocess,shutil
import numpy as np,av
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from annotation_store import load_annotation,sha

def main():
 checks=[]
 def check(name,condition,details=None):
  checks.append({'check':name,'pass':bool(condition),'details':details})
  if not condition:raise AssertionError(name)
 required=['SMART_TRACKING_V2_DEVELOPMENT.md','ANNOTATION_READY.md','MICROBUBBLE_ANNOTATION_PROTOCOL.md',
           'DATASET_SPLIT_AND_PROVENANCE.md','OPENCV_VS_CELLPOSE.md','LIVE_INFERENCE_BENCHMARK.md',
           'GROUP_AND_BUBBLE_ADMISSION.md','FINAL_GATE.json','FINAL_RECOMMENDATION.md']
 for name in required:check('report '+name,(ROOT/name).stat().st_size>100)
 protected=json.loads((ROOT/'evidence/protected_sources_before.json').read_text());after=[]
 for r in protected:
  actual=sha(r['path']);check('protected '+r['path'],actual==r['sha256']);after.append({**r,'after_sha256':actual,'unchanged':True})
 (ROOT/'evidence/protected_sources_after.json').write_text(json.dumps(after,indent=2))
 rows=json.loads((ROOT/'data/candidate_manifest.json').read_text())['rows'];check('84 unique original candidates',len(rows)==len({r['id'] for r in rows})==84)
 confirmed=0
 for r in rows:
  check('candidate '+r['id'],sha(ROOT/r['image'])==r['image_sha256'] and not r['resized'])
  meta,masks=load_annotation(ROOT/'data/annotations',r)
  confirmed+=bool(meta and meta['expert_confirmed'])
 check('original video checksum',sha(rows[0]['source_file'])==rows[0]['source_sha256'])
 blocks={k:[r['source_frame_index'] for r in rows if r['split']==k] for k in ['train','validation','test']}
 check('temporal split guards',max(blocks['train'])<min(blocks['validation'])<max(blocks['validation'])<min(blocks['test']))
 for name in ['annotation_repair_sam_smoke.json','annotation_repair_manual_smoke.json']:
  d=json.loads((ROOT/'evidence'/name).read_text());check(name,d['PASS'])
 check('software contracts passed','Ran 8 tests' in (ROOT/'logs/software_contract_tests_final.log').read_text() and '\nOK\n' in (ROOT/'logs/software_contract_tests_final.log').read_text())
 for kind in ['png','pdf']:
  check('four figures '+kind,len(list((ROOT/'figures').glob('Fig_VisionV2_*.'+kind)))==4)
  for p in (ROOT/'figures').glob('*.'+kind):
   check('figure signature '+p.name,p.read_bytes().startswith(b'%PDF-' if kind=='pdf' else b'\x89PNG'))
 info=json.loads((ROOT/'evidence/video_validation.json').read_text());timing=[]
 live=[json.loads(s) for s in (ROOT/'results/live_opencv_20s/recorded/frames.jsonl').read_text().splitlines()]
 paired=[json.loads(s) for s in (ROOT/'results/paired_comparison/records.jsonl').read_text().splitlines()]
 expected={'Video_VisionV2_live_tracking.mp4':[r['source_time_s'] for r in live],
           'Video_VisionV2_method_comparison.mp4':[r['source_time_s'] for r in paired if r['kind']=='CAUSAL_REPLAY_OFFLINE' and r['route']=='opencv']}
 for name,v in info.items():
  path=ROOT/'videos'/name;check('video SHA '+name,sha(path)==v['sha256'])
  with av.open(str(path)) as c:
   st=c.streams.video[0];times=[float(f.time) for f in c.decode(st)]
  source=np.asarray(expected[name]);source-=source[0]
  check('no frame duplication or omission '+name,len(times)==len(source))
  delta=float(np.max(np.abs(np.array(times)-source)))
  check('source timestamp reproduction '+name,delta<=1/90000+1e-8,delta)
  timing.append({'video':name,'max_pts_difference_s':delta,'source_frames':len(source),'no_interpolation':True})
 (ROOT/'evidence/video_source_timing_verification.json').write_text(json.dumps(timing,indent=2))
 for model,expected_sha in [('models/sam_vit_b_01ec64.pth','ec2df62732614e57411cdcf32a23ffdf28910380d03139ee0f4fcbe91eb8c912'),
                            ('models/cellpose/cyto3','2dc3087a8abd7da46d1ab0ddd5824639933cc3ff63b382af3fa1939a392db93c')]:
  check('model '+model,sha(ROOT/model)==expected_sha)
 for env in ['.venv','.venv_annotation']:
  installed=json.loads(subprocess.check_output([str(ROOT/env/'bin/python'),'-m','pip','list','--format=json'],text=True))
  check('no SAM2 '+env,not any(r['name'].lower().replace('-','').replace('_','')=='sam2' for r in installed))
 gate=json.loads((ROOT/'FINAL_GATE.json').read_text())
 audit=json.loads((ROOT/'results/annotation_audit_current.json').read_text())
 check('truthful annotation snapshot',audit['expert_confirmed_frames']==gate['expert_confirmed_frames'])
 for r in audit['frames']:
  if r['annotation_revision']:
   folder=ROOT/'data/annotations'/r['id']
   check('snapshot annotation revision '+r['id'],sha(folder/(r['annotation_revision']+'.json'))==r['metadata_sha256'] and sha(folder/r['mask_file'])==r['mask_sha256'])
 check('control always disabled',gate['HARDWARE_CONTROL_DISABLED'] and not gate['BUBBLE_DETECTION_VALIDATED'])
 for tag in ['live_opencv_20s','live_cellpose_20s_resumed']:
  mock=json.loads((ROOT/'results'/tag/'mock.json').read_text())
  check('mock controls '+tag,all(not r['control_permission'] and not r['bubble_valid'] for r in mock))
 doc=json.loads((ROOT/'evidence/progress_document/audit.json').read_text())
 check('progress backup',sha(doc['backup'])==doc['before_sha256'])
 check('progress latest',sha('/home/lzy/projects/microbubble_aggregation.md')==doc['after_sha256'])
 compact=ROOT/'evidence/compact_results';compact.mkdir(exist_ok=True)
 for tag in ['live_opencv_20s','live_cellpose_20s_resumed','paired_comparison']:
  for name in ['summary.json','execution_identity.json']:
   shutil.copy2(ROOT/'results'/tag/name,compact/f'{tag}_{name}')
 for name in ['software_contract_tests_final.log','annotation_close_and_navigation_test.log','annotation_repair_final.log','paired_comparison.log']:
  shutil.copy2(ROOT/'logs'/name,compact/name)
 report={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'PASS':True,'checks':checks,'count':len(checks)}
 (ROOT/'evidence/delivery_verification.json').write_text(json.dumps(report,indent=2))
 manifest=write_manifest(report['utc'])
 print(json.dumps({'PASS':True,'checks':len(checks),'hashed_files':len(manifest['files']),'expert_confirmed_frames':confirmed}))

def write_manifest(utc=None):
 # Inputs/artifacts are retained locally; environments, weights, external repos and changing labels excluded.
 omit={'.venv','.venv_annotation','third_party','models','__pycache__','.git'};files=[]
 for current,dirs,names in os.walk(ROOT):
  dirs[:]=sorted(d for d in dirs if d not in omit)
  for name in sorted(names):
   path=Path(current)/name;rel=path.relative_to(ROOT).as_posix()
   if rel.startswith('data/annotations/') or rel.startswith('results/annotation_jobs/') or name.endswith(('.lock','.pyc')):continue
   if rel in ['evidence/delivery_manifest.json','evidence/final_git_identity.json','logs/delivery_verification_final.log']:continue
   if rel.startswith('logs/annotation_') and rel not in ['logs/annotation_repair_final.log','logs/annotation_close_and_navigation_test.log']:continue
   if path.is_file():files.append({'path':rel,'bytes':path.stat().st_size,'sha256':sha(path)})
 manifest={'utc':utc or datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'source, fixed candidates, saved runs, logs and artifacts; dynamic user annotations excluded',
           'files':files,'model_identity_in':'MODEL_AND_ENVIRONMENT.md','annotation_snapshot':'results/annotation_audit_current.json'}
 (ROOT/'evidence/delivery_manifest.json').write_text(json.dumps(manifest,indent=2))
 return manifest
if __name__=='__main__':
 if sys.argv[1:]==['--manifest-only']:
  manifest=write_manifest();print(json.dumps({'manifest_refreshed':True,'hashed_files':len(manifest['files']),'checks_rerun':False}))
 elif not sys.argv[1:]:main()
 else:raise SystemExit('Usage: verify_delivery.py [--manifest-only]')
