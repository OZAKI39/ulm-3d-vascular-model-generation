"""Apply the user's 2026-10-16 correction: frames 1--5 were only partly labeled.

Requires the GUI to be closed. Retains every pixel and immutable revision; only
the active review scope is revoked. Old model candidates are annotated, not erased.
"""
from pathlib import Path
from datetime import datetime, timezone
import json,sys,shutil,subprocess
import numpy as np
from PyQt6.QtCore import QLockFile
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from annotation_store import load_annotation,sha
from annotation_snapshot import revoke_full_frame,dataset_readiness

IDS=('video_000000','video_000018','video_000036','video_000054','video_000072')

def main():
 lock=QLockFile(str(ROOT/'data/annotation_workbench.lock'));lock.setStaleLockTime(0)
 if not lock.tryLock(100):raise RuntimeError('Close the annotation workbench before scope correction')
 try:
  directory=ROOT/'data/annotations';stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
  out=ROOT/'evidence/review_scope_correction_20261016'/stamp;out.mkdir(parents=True,exist_ok=False)
  records={r['id']:r for r in json.loads((ROOT/'data/candidate_manifest.json').read_text())['rows']}
  loaded={key:load_annotation(directory,records[key]) for key in IDS}
  if not all(meta and meta['expert_confirmed'] for meta,_ in loaded.values()):
   raise ValueError('Correction already applied or review states changed; do not overwrite newer work')
  active=ROOT/'data/auto_training/active_model.json'
  if active.exists() and json.loads(active.read_text()).get('admitted'):
   raise ValueError('An admitted model exists; inspect its training snapshot before correction')
  before={str(p.relative_to(directory)):sha(p) for p in sorted(directory.rglob('*')) if p.is_file()}
  report={'user_clarification':'Frames 1--5 contain partial labels; whole-frame review was not completed.',
          'utc':datetime.now(timezone.utc).isoformat(),'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
          'before_sha256':before,'corrections':[],'invalidated_training_runs':[]}
  (out/'report.json').write_text(json.dumps(report,indent=2))
  for key,(meta,masks) in loaded.items():
   backup=out/'before'/key;backup.mkdir(parents=True)
   for name in ('latest.json',meta['mask_file'],meta['revision']+'.json'):shutil.copy2(directory/key/name,backup/name)
   _,new=revoke_full_frame(directory,records[key]);actual,after_masks=load_annotation(directory,records[key])
   unchanged=all(np.array_equal(a,after_masks[k]) for k,a in masks.items())
   if not unchanged:raise AssertionError('Unexpected pixel modification')
   report['corrections'].append({'id':key,'before_revision':meta['revision'],'after_revision':new['revision'],
                                'all_class_pixels_unchanged':unchanged,'saved_counts':actual['instance_counts'],
                                'status':actual['status']})
   (out/'report.json').write_text(json.dumps(report,indent=2))
  for samples_path in sorted((ROOT/'results/sam_auto_training').glob('*/samples.json')):
   samples=json.loads(samples_path.read_text());affected=sorted({s['image_id'] for s in samples.get('train',[]) if s['image_id'] in IDS})
   if not affected:continue
   invalid={'status':'INVALIDATED_BY_USER_REVIEW_SCOPE_CORRECTION','affected_frames':affected,
            'reason':'Partial annotations were marked as whole-frame reviewed; unmarked pixels are not valid background truth.',
            'training_or_generalization_evidence_usable':False,'correction_report':str(out/'report.json'),
            'original_results_and_checkpoint_preserved':True}
   (samples_path.parent/'review_scope_invalidation.json').write_text(json.dumps(invalid,indent=2))
   report['invalidated_training_runs'].append(str(samples_path.parent.relative_to(ROOT)))
  after={str(p.relative_to(directory)):sha(p) for p in sorted(directory.rglob('*')) if p.is_file()}
  expected={key+'/latest.json' for key in IDS}
  modified={name for name,digest in before.items() if after.get(name)!=digest}
  if modified!=expected:raise AssertionError('An unrelated existing annotation file changed')
  report.update(after_sha256=after,only_expected_pointers_changed=True,
                readiness=dataset_readiness(ROOT),all_pixels_preserved=True,
                source_sha256={name:sha(ROOT/name) for name in ['scripts/correct_partial_reviews_20261016.py','src/annotation_snapshot.py','src/annotation_store.py']})
  (out/'report.json').write_text(json.dumps(report,indent=2))
  (out.parent/'latest.json').write_text(json.dumps({'run':str(out)},indent=2))
  print(json.dumps({'output':str(out),'corrected_frames':len(report['corrections']),
                    'pixels_preserved':True,'invalidated_training_runs':len(report['invalidated_training_runs']),
                    'readiness':report['readiness']['split_summary']},ensure_ascii=False))
 finally:lock.unlock()

if __name__=='__main__':main()
