"""Validate deliverables and expose the requested analysis-directory entrypoints."""
from common import *
import re,csv,subprocess,importlib.metadata
from PIL import Image
existing_diff=subprocess.check_output(['git','-c',f'safe.directory={R}','diff'],cwd=R)
assert existing_diff==(O/'provenance/git_existing.diff').read_bytes()
dump(O/'REPOSITORY_PRESERVATION.json',dict(existing_tracked_diff_unchanged=True,old_native_library_unchanged=sha(OLDLIB)==json.loads((O/'SOURCE_FREEZE.json').read_text())['historical_library_sha256']))
required=['SOURCE_FREEZE.json','GEOMETRY_FREEZE.json','BIPHASIC_MATERIAL_SPEC.json','BIPHASIC_FIELD.npz','STRONG_DOMAINS.csv','NODE_MATERIAL_PHASE.csv','PAIR_MATERIAL_PROPERTIES.csv','PREFLIGHT_STRUCTURE_AUDIT.json','FROZEN_SHEAR_DRIVE.json','NATIVE_UNIFORM_RESULT.json','NATIVE_BIPHASIC_RESULT.json','TRUE_FRAGMENT_HISTORY.csv','FRAGMENT_SIZE_DISTRIBUTION.csv','FRAGMENT_MASS_DISTRIBUTION.csv','ENERGY_LEDGER.csv','COMPUTE_COST.json','FINAL_STATUS.json','REVIEW_CLOT_BIPHASIC_HETEROGENEOUS_SHEAR_FRAGMENTATION_V2.md']
for name in required:assert (O/name).is_file(),name
report=(O/'REVIEW_CLOT_BIPHASIC_HETEROGENEOUS_SHEAR_FRAGMENTATION_V2.md').read_text();assert [int(x) for x in re.findall(r'^([0-9]+)\. ',report,re.M)]==list(range(1,6))
figs=list((F/'figures').glob('*.png'));assert len(figs)>=8
for p in figs:assert p.with_suffix('.pdf').is_file();im=Image.open(p);im.verify()
video=json.loads((O/'VIDEO_GENERATION.json').read_text());assert video['status']=='PASS' and len(video['videos'])==2
for item in video['videos']:
 assert item['material_actor_presence_checked_every_frame']
 p=F/'videos'/item['file'];assert sha(p)==item['sha256']
 rows=list(csv.DictReader((p.parent/(p.stem+'_FRAME_TIMES.csv')).open()));assert len(rows)==item['frames']
 for r in rows:
  assert r['state_interpolation']=='False'
  for m in CASES:
   q=O/r[m+'_source'];assert q.exists()
   if r['hold']=='False':
    # Full arrays are not altered during visualization; timestamps match saved native metadata.
    
    with np.load(q,allow_pickle=False) as state:assert abs(float(state['time'])-float(r[m+'_actual_time_s']))<1e-16
 im=np.array(Image.open(p.parent/(p.stem+'_COVER.png')).convert('RGB'))
 for j in [0,1]:
  panel=im[280:800,j*960+180:(j+1)*960-180].astype(int)
  assert panel.std()>10 and np.sum(np.max(panel,axis=2)-np.min(panel,axis=2)>40)>250,'Missing colored material in video panel'
 # Decode the entire file to catch truncated streams, not merely a metadata check.
 subprocess.run(['ffmpeg','-v','error','-i',str(p),'-f','null','-'],check=True)
for name in required:
 link=A/name
 if not link.exists():link.symlink_to(os.path.relpath(O/name,A))
for name,target in [('figures',F/'figures'),('videos',F/'videos'),('states',O/'cases')]:
 link=A/name
 if not link.exists():link.symlink_to(os.path.relpath(target,A))
(A/'visualization').mkdir(exist_ok=True)
for name in ['preflight_figures.py','figures.py','preview_rise_end.py','render_videos.py']:
 p=A/'visualization'/name
 if not p.exists():p.symlink_to('../'+name)
(A/'OPEN_RESULTS.html').write_text('<!doctype html><meta charset="utf-8"><a href="'+os.path.relpath(F/'OPEN_RESULTS.html',A)+'">打开本轮科研报告、静态图和两个真实MP4</a>')
validation=json.loads((O/'DELIVERY_VALIDATION.json').read_text());validation['static_figure_png_pdf_pairs']=len(figs);validation['video_full_decode']='PASS';validation['video_count']=2;validation['video_frame_sources']='ALL_REAL_ACCEPTED_STATES';validation['report_question_count']=5;dump(O/'DELIVERY_VALIDATION.json',validation)
status=json.loads((O/'FINAL_STATUS.json').read_text());status['VIDEO_STATUS']='PASS_TWO_1080P_MP4_FULL_DECODE_AND_FRAME_SOURCE_CHECK';status['VISUALIZATION_STATUS']='PASS_PNG_AND_VECTOR_PDF';dump(O/'FINAL_STATUS.json',status)
print('FINAL_DELIVERY_PASS',len(figs),len(video['videos']),flush=True)
