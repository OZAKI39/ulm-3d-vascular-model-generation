"""Data correction / immutable export fixtures; no real labels are modified."""
import json,sys,tempfile,unittest
from pathlib import Path
import numpy as np
from PIL import Image
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from annotation_store import CLASSES,save_annotation,load_annotation,sha
from annotation_snapshot import revoke_full_frame,freeze_reviewed,export_frozen,dataset_readiness
from sam_auto_training import snapshot,snapshot_still_valid


class Snapshots(unittest.TestCase):
 def setUp(self):
  self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name);(self.root/'data').mkdir()
  Image.fromarray(np.full((60,80,3),123,np.uint8)).save(self.root/'image.png')
  self.record={'id':'fixture','image':'image.png','image_sha256':sha(self.root/'image.png'),
               'source_sha256':'fixture','split':'train','resolution_wh':[80,60]}
  (self.root/'data/candidate_manifest.json').write_text(json.dumps({'rows':[self.record]}))
  (self.root/'protocol.json').write_text(json.dumps({'no_random_adjacent_split':True}))
  self.masks={k:np.zeros((60,80),np.uint32) for k in CLASSES};self.masks['SINGLE_BUBBLE'][20:25,30:35]=4
  self.directory=self.root/'data/annotations'
 def tearDown(self):self.temp.cleanup()
 def save(self,region=False):
  return save_annotation(self.directory,self.record,self.masks,'Fixture reviewer',True,{'SINGLE_BUBBLE:4':'a'},'Keep notes',
                         review_scope='REGIONS' if region else 'FULL_FRAME',confirm_region=[10,20,40,60] if region else None)
 def test_revoke_preserves_pixels_notes_tracks_and_history_but_invalidates_training(self):
  old=self.save();snap=snapshot(self.root);old_file=self.directory/'fixture'/old['mask_file'];digest=sha(old_file)
  _,new=revoke_full_frame(self.directory,self.record);meta,masks=load_annotation(self.directory,self.record)
  for key in CLASSES:np.testing.assert_array_equal(masks[key],self.masks[key])
  self.assertEqual(meta['notes'],old['notes']);self.assertEqual(meta['track_ids'],old['track_ids'])
  self.assertEqual(meta['status'],'DRAFT_UNREVIEWED');self.assertEqual(meta['review_scope'],'REGIONS')
  self.assertEqual(sha(old_file),digest);self.assertFalse(snapshot_still_valid(snap,self.root));self.assertEqual(snapshot(self.root)['rows'],[])
 def test_freeze_exports_captured_revision_not_new_live_work(self):
  self.save(region=True);destination=self.root/'frozen';frozen=freeze_reviewed(self.root,destination)
  self.masks['SINGLE_BUBBLE'].fill(0);save_annotation(self.directory,self.record,self.masks)
  result=export_frozen(destination)
  self.assertEqual(len(result['train']),1);self.assertEqual(result['validation'],[])
  self.assertEqual(result['train'][0]['box_yxyx'],[10,20,40,60])
  import tifffile
  a=tifffile.imread(destination/'cellpose_train/fixture_roi000_masks.tif')
  self.assertEqual(np.count_nonzero(a),25)
  with Image.open(destination/'cellpose_train/fixture_roi000.png') as image:self.assertEqual(image.size,(40,30))
  self.assertTrue(all(sha(destination/name)==digest for name,digest in frozen['snapshot_sha256'].items()))
 def test_draft_never_frozen_and_destination_not_overwritten(self):
  save_annotation(self.directory,self.record,self.masks)
  with self.assertRaisesRegex(ValueError,'No reviewed'):freeze_reviewed(self.root,self.root/'draft')
  self.save(region=True);freeze_reviewed(self.root,self.root/'good')
  with self.assertRaises(FileExistsError):freeze_reviewed(self.root,self.root/'good')
 def test_tampered_snapshot_blocks_export(self):
  self.save(region=True);dest=self.root/'frozen';frozen=freeze_reviewed(self.root,dest)
  (dest/frozen['rows'][0]['frozen_mask_file']).write_bytes(b'changed')
  with self.assertRaisesRegex(ValueError,'Frozen input'):export_frozen(dest)
  self.assertFalse((dest/'cellpose_train').exists())
 def test_train_crops_do_not_supply_independent_validation(self):
  self.save(region=True);ready=dataset_readiness(self.root)
  self.assertFalse(ready['independent_evaluation_ready']);self.assertEqual(ready['split_summary']['train']['reviewed_frames'],1)
  self.assertEqual(ready['split_summary']['train']['reviewed_counts']['SINGLE_BUBBLE'],1)

if __name__=='__main__':unittest.main()
