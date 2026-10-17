"""Recompute claims from saved runs; never invokes camera or instrument code."""
from pathlib import Path
import json,hashlib,argparse,math,ast,subprocess
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 return h.hexdigest()
def verify():
 checks=[]
 def check(name,value):
  checks.append({'name':name,'passed':bool(value)})
  if not value:raise AssertionError(name)
 for p in (ROOT/'scripts').glob('*.py'):ast.parse(p.read_text())
 check('source_input_sha256',sha(Path(read(ROOT/'evidence/input_identity.json')['path']))==read(ROOT/'evidence/input_identity.json')['sha256'])
 before=read(ROOT/'evidence/protected_sources_before.json')
 check('46_original_sources_unchanged',all(sha(ROOT.parent/r['path'])==r['sha256'] for r in before))
 checks_done=read(ROOT/'evidence/interface_checks.json');check('actual_receiver_and_mock_tests_pass',all(r['passed'] for r in checks_done))
 for folder in sorted((ROOT/'results').iterdir()):
  if not folder.is_dir() or not(folder/'summary.json').exists():continue
  s=read(folder/'summary.json')
  if folder.name=='cotracker_online':
   w=read(folder/'windows.json');check('cotracker_online_bounded_inputs',s['online_only'] and all(len(r['input_frame_ids'])<=16 for r in w));check('cotracker_no_false_fresh_admission',s['fresh_emissions']==sum(r['latest_frame_age_ms']<=100 for r in w));continue
  frames=read(folder/'frames.json');decoder=read(folder/'decoder.json');mock=read(folder/'mock_received.json');fid=[r['frame_id'] for r in frames]
  check(folder.name+'_monotonic_unique_frame_ids',fid==sorted(set(fid)))
  check(folder.name+'_honest_drop_count',s['dropped_latest_frames']==len(decoder)-len(frames) and s['processed_frames']==len(frames))
  check(folder.name+'_honest_valid_rate',math.isclose(s['mock_valid_visual_rate_wall_hz'],sum(r['valid_visual'] for r in mock)/s['wall_s']))
  check(folder.name+'_expired_messages_rejected',all(not r['valid_visual'] for r in mock if r['age_ms']>100))
  check(folder.name+'_no_hardware_admission',all(not r['bubble_control_admitted'] for r in mock))
  for line in (folder/'states.jsonl').read_text().splitlines():
   m=json.loads(line)
   assert m['capture_timestamp'] is None and not m['control_permission']
   assert m['group_shape'] is None
   assert all(not o['visual_valid'] for o in m['targets'] if o['status']=='LOST')
  check(folder.name+'_no_fabricated_capture_or_group_shape',True)
 figchecks=[]
 for p in sorted((ROOT/'figures').glob('*.png')):
  with Image.open(p) as im:im.verify()
  pdf=p.with_suffix('.pdf');check(p.stem+'_png_pdf_pair',pdf.exists() and pdf.read_bytes().startswith(b'%PDF-') and b'%%EOF' in pdf.read_bytes()[-1024:]);figchecks.append(p.stem)
 check('four_real_figure_pairs',len(figchecks)==4)
 for f in ['OPEN_SOURCE_TRACKING_SURVEY.md','EXISTING_INTERFACES_AUDIT.md','REALTIME_TRACKING_BENCHMARK.md','FINAL_RECOMMENDATION.md','INTERFACE_DESIGN.md']:
  check(f+'_exists',(ROOT/f).stat().st_size>1000)
 return checks
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--write-manifest',action='store_true');args=p.parse_args();checks=verify()
 result={'checks':checks,'passed':all(c['passed'] for c in checks),'check_count':len(checks)}
 if args.write_manifest:
  (ROOT/'evidence/delivery_validation.json').write_text(json.dumps(result,indent=2))
  rows=[]
  for f in sorted(ROOT.rglob('*')):
   if not f.is_file() or any(x in f.relative_to(ROOT).parts for x in ['.git','.venv','third_party','__pycache__']):continue
   if f.name in ['INPUT_OUTPUT_SHA256.json','delivery_receipt.json']:continue
   rows.append({'path':f.relative_to(ROOT).as_posix(),'bytes':f.stat().st_size,'sha256':sha(f)})
  (ROOT/'INPUT_OUTPUT_SHA256.json').write_text(json.dumps({'files':rows,'source_video':read(ROOT/'evidence/input_identity.json'),'third_party_source':'evidence/source_execution_identity.json','checkpoint':'evidence/cotracker_weights.json','excluded':'installed environments, third_party trees, .git, manifest itself, delivery receipt'},indent=2))
 else:
  manifest=read(ROOT/'INPUT_OUTPUT_SHA256.json');changed=[v['path'] for v in manifest['files'] if sha(ROOT/v['path'])!=v['sha256']];assert not changed,changed
 print(json.dumps({'passed':True,'checks':len(checks),'manifest': 'written' if args.write_manifest else 'verified'}))
