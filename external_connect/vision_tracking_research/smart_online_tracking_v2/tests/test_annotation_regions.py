"""Review-coverage and instance-edit invariants; temporary software fixtures only."""
import unittest,tempfile,json,sys
from pathlib import Path
import numpy as np
from PIL import Image
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from annotation_store import CLASSES,save_annotation,load_annotation,sha,export_cellpose
from annotation_regions import reviewed_boxes,valid_regions
from annotation_batch import merge_proposals,split_instance
from sam_auto_training import snapshot,build_samples,snapshot_still_valid,admission

class Regions(unittest.TestCase):
 def setUp(self):
  self.temp=tempfile.TemporaryDirectory(prefix='REGION_SOFTWARE_FIXTURE_');self.root=Path(self.temp.name)
  self.folder=self.root/'data/annotations';self.folder.mkdir(parents=True)
  image=self.root/'image.png';Image.fromarray(np.full((256,320,3),170,np.uint8)).save(image)
  self.r={'id':'fixture','image':'image.png','image_sha256':sha(image),'source_sha256':'fixture','split':'train','resolution_wh':[320,256]}
  (self.root/'data/candidate_manifest.json').write_text(json.dumps({'rows':[self.r]}))
  self.m={k:np.zeros((256,320),np.uint32) for k in CLASSES}
  self.m['SINGLE_BUBBLE'][80:86,100:106]=1
  self.m['SINGLE_BUBBLE'][80:86,130:136]=2
  self.m['SINGLE_BUBBLE'][190:196,250:256]=3
  self.box=[40,60,160,180]
 def tearDown(self):self.temp.cleanup()
 def save(self,**kwargs):
  return save_annotation(self.folder,self.r,self.m,reviewer='Fixture human',review_scope='REGIONS',**kwargs)
 def test_partial_confirmation_is_not_full_and_draft_outside_survives(self):
  meta=self.save(confirm=True,confirm_region=self.box)
  self.assertFalse(meta['expert_confirmed']);self.assertEqual(meta['status'],'REGIONS_CONFIRMED')
  self.assertEqual(meta['reviewed_instance_counts']['SINGLE_BUBBLE'],2)
  _,m=load_annotation(self.folder,self.r);self.assertTrue(np.array_equal(m['SINGLE_BUBBLE'],self.m['SINGLE_BUBBLE']))
  self.assertEqual(reviewed_boxes(meta,m),[self.box])
 def test_only_edited_region_invalidated(self):
  a=self.save(confirm=True,confirm_region=self.box)
  b=self.save(confirm=True,review_regions=a['review_regions'],confirm_region=[180,220,230,290])
  self.m['SINGLE_BUBBLE'][80,100]=0
  c=self.save(review_regions=b['review_regions'])
  self.assertEqual(len(c['review_regions']),1);self.assertEqual(c['invalidated_region_ids'],[a['review_regions'][0]['id']])
  self.assertFalse(c['expert_confirmed'])
 def test_outside_class_overlap_allowed_inside_rejected(self):
  self.m['AGGREGATE'][190:196,250:256]=1
  self.save(confirm=True,confirm_region=self.box)
  self.m['AGGREGATE'][80:86,100:106]=1
  with self.assertRaises(ValueError):self.save(confirm=True,confirm_region=self.box)
 def test_draft_never_used_and_roi_samples_have_no_outside_context(self):
  self.save();self.assertEqual(snapshot(self.root)['rows'],[])
  a=self.save(confirm=True,confirm_region=self.box);snap=snapshot(self.root)
  samples,excluded=build_samples(snap,'train',self.root)
  self.assertEqual({s['instance'] for s in samples},{1,2});self.assertEqual(excluded,[])
  for sample in samples:
   y0,x0,y1,x1=sample['crop_yxyx'];ry0,rx0,ry1,rx1=self.box
   self.assertTrue(ry0<=y0<y1<=ry1 and rx0<=x0<x1<=rx1)
   yy,xx=sample['point_yx'];self.assertTrue(y0<=yy<y1 and x0<=xx<x1)
   self.assertTrue(np.array_equal(sample['mask'],self.m['SINGLE_BUBBLE'][y0:y1,x0:x1]==sample['instance']))
  self.assertTrue(snapshot_still_valid(snap,self.root))
  self.save(review_regions=[]);self.assertFalse(snapshot_still_valid(snap,self.root))
 def test_overlapping_review_regions_deduplicate_instances(self):
  a=self.save(confirm=True,confirm_region=self.box)
  self.save(confirm=True,review_regions=a['review_regions'],confirm_region=[30,50,170,190])
  samples,_=build_samples(snapshot(self.root),'train',self.root)
  self.assertEqual(len(samples),2)
 def test_truncated_object_excluded(self):
  self.m['SINGLE_BUBBLE'][78:86,175:187]=4
  self.save(confirm=True,confirm_region=self.box)
  samples,excluded=build_samples(snapshot(self.root),'train',self.root)
  self.assertNotIn(4,[s['instance'] for s in samples]);self.assertTrue(any(e['instance']==4 and e['reason']=='target_truncated' for e in excluded))
 def test_export_crops_only_reviewed_area(self):
  self.save(confirm=True,confirm_region=self.box)
  rows=export_cellpose(self.root,[self.r],self.root/'export')
  self.assertEqual(rows[0]['source_id'],'fixture');self.assertEqual(rows[0]['box_yxyx'],self.box)
  im=Image.open(self.root/'export/fixture_roi000.png');self.assertEqual(im.size,(120,120))
  import tifffile
  self.assertEqual(int(tifffile.imread(self.root/'export/fixture_roi000_masks.tif').max()),2)
 def test_ignore_region_is_not_exported_as_background(self):
  self.m['UNCERTAIN'][50:55,70:75]=1;self.save(confirm=True,confirm_region=self.box)
  rows=export_cellpose(self.root,[self.r],self.root/'export')
  self.assertTrue(rows[0]['status'].startswith('EXCLUDED'));self.assertEqual(list((self.root/'export').glob('*.png')),[])
 def test_two_rois_same_frame_not_independent_validation(self):
  self.r['split']='validation';(self.root/'data/candidate_manifest.json').write_text(json.dumps({'rows':[self.r]}))
  a=self.save(confirm=True,confirm_region=self.box)
  self.save(confirm=True,review_regions=a['review_regions'],confirm_region=[180,220,230,290])
  samples,_=build_samples(snapshot(self.root),'validation',self.root)
  frames=len(set(s['image_id'] for s in samples));self.assertEqual(frames,1)
  self.assertFalse(admission({'validation_frames':frames,'validation_objects':100})[0])
 def test_invalid_region_metadata_fails_closed(self):
  meta=self.save(confirm=True,confirm_region=self.box)
  meta['review_regions'][0]['box_yxyx']=[-1,0,30,40]
  self.assertEqual(valid_regions(meta,self.m),[])
 def test_roi_without_reviewer_not_confirmed(self):
  with self.assertRaises(ValueError):save_annotation(self.folder,self.r,self.m,confirm_region=self.box,review_scope='REGIONS')

class InstanceEdits(unittest.TestCase):
 def test_merge_preserves_existing_pixels_and_skips_whole_conflict(self):
  masks={k:np.zeros((20,20),np.uint32) for k in CLASSES};masks['SINGLE_BUBBLE'][5,5]=9
  masks['UNCERTAIN'][9,9]=1
  proposals=np.zeros((10,10),np.uint32);proposals[1:3,1:3]=1;proposals[5:7,5:7]=2;proposals[8:10,8:10]=3
  result,accepted,skipped=merge_proposals(masks,proposals,[0,0,10,10])
  self.assertEqual(skipped,[2,3]);self.assertEqual(accepted,[{'candidate_id':1,'instance_id':10}])
  self.assertEqual(result[5,5],9);self.assertEqual(result[6,6],0);self.assertEqual(result[9,9],0)
 def test_split_exact_union_and_reject_invalid_seeds(self):
  labels=np.zeros((40,40),np.uint32);labels[10:20,8:27]=4
  out,new=split_instance(labels,4,[[15,12],[15,23]])
  self.assertEqual(new,5);self.assertTrue(np.array_equal(out>0,labels>0));self.assertEqual(set(np.unique(out)),{0,4,5})
  for pts in ([[15,12],[15,12]],[[1,1],[15,23]],[[15,12]]):
   with self.assertRaises(ValueError):split_instance(labels,4,pts)
 def test_disconnected_split_refused(self):
  labels=np.zeros((40,40),np.uint32);labels[10:12,10:12]=1;labels[20:22,20:22]=1
  with self.assertRaises(ValueError):split_instance(labels,1,[[10,10],[20,20]])

if __name__=='__main__':unittest.main(verbosity=2)
