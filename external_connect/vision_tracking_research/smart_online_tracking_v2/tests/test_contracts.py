"""Software fixtures only: never writes expert labels to the research dataset."""
from pathlib import Path
import tempfile,sys,unittest,json,time
import numpy as np,cv2
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from annotation_store import CLASSES,save_annotation,load_annotation,export_cellpose,sha
from vision_engine import VisionEngine,region_shape,neighborhood_for_admitted_masks
from runtime import enforce_freshness,MockConsumer,invalid_state

class Contracts(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
  (self.root/'data/annotations').mkdir(parents=True);cv2.imwrite(str(self.root/'input.png'),np.zeros((30,40,3),np.uint8))
  self.record={'id':'synthetic_software_fixture','resolution_wh':[40,30],'image_sha256':sha(self.root/'input.png'),
               'source_sha256':'SYNTHETIC_NOT_RESEARCH','image':'input.png','split':'train'}
  self.masks={k:np.zeros((30,40),np.uint32) for k in CLASSES}
 def tearDown(self):self.tmp.cleanup()
 def test_multiclass_revision_roundtrip(self):
  for i,k in enumerate(CLASSES):self.masks[k][i*5:i*5+3,2:5]=i+1
  a=save_annotation(self.root/'data/annotations',self.record,self.masks)
  b,loaded=load_annotation(self.root/'data/annotations',self.record)
  self.assertFalse(b['expert_confirmed']);self.assertEqual(b['instance_counts'],{k:1 for k in CLASSES})
  for k in CLASSES:np.testing.assert_array_equal(loaded[k],self.masks[k])
  c=save_annotation(self.root/'data/annotations',self.record,self.masks)
  self.assertNotEqual(a['revision'],c['revision'])
 def test_confirmation_gates_and_checksum(self):
  self.masks['SINGLE_BUBBLE'][2:5,2:5]=1;self.masks['AGGREGATE'][3:6,3:6]=2
  with self.assertRaises(ValueError):save_annotation(self.root/'data/annotations',self.record,self.masks,'fixture',True)
  self.masks['AGGREGATE'][:]=0
  with self.assertRaises(ValueError):save_annotation(self.root/'data/annotations',self.record,self.masks,'',True)
  with self.assertRaises(ValueError):save_annotation(self.root/'data/annotations',self.record,self.masks,'SOFTWARE_TEST',True)
  a=save_annotation(self.root/'data/annotations',self.record,self.masks)
  (self.root/'data/annotations'/self.record['id']/a['mask_file']).write_bytes(b'corrupt')
  with self.assertRaises(ValueError):load_annotation(self.root/'data/annotations',self.record)
 def test_export_pending_and_confirmed_fixture(self):
  save_annotation(self.root/'data/annotations',self.record,self.masks)
  self.assertEqual(export_cellpose(self.root,[self.record],self.root/'pending'),[])
  # Artificial reviewer exists only in this TemporaryDirectory, never in research data.
  self.masks['SINGLE_BUBBLE'][2:5,2:5]=1
  save_annotation(self.root/'data/annotations',self.record,self.masks,'UNIT_FIXTURE_ONLY',True)
  rows=export_cellpose(self.root,[self.record],self.root/'export');self.assertEqual(rows[0]['status'],'EXPORTED')
  import tifffile
  self.assertEqual(int(tifffile.imread(self.root/'export/synthetic_software_fixture_masks.tif').max()),1)
  self.masks['UNCERTAIN'][20:22,20:22]=1
  save_annotation(self.root/'data/annotations',self.record,self.masks,'UNIT_FIXTURE_ONLY',True)
  self.assertTrue(export_cellpose(self.root,[self.record],self.root/'ignored')[0]['status'].startswith('EXCLUDED'))
 def test_shape_tensor_translation(self):
  mask=np.zeros((40,60),np.uint8);mask[10:20,20:40]=1
  s=region_shape(mask,(100,200));self.assertEqual(s['area_px2'],200)
  np.testing.assert_allclose(s['centroid_xy_px'],[129.5,214.5])
  np.testing.assert_allclose(s['Q_px2'],[[33.25,0],[0,8.25]])
  self.assertAlmostEqual(s['Rg_px']**2,41.5)
 def test_stale_and_no_frame_erase_coordinates(self):
  m={'arrival_mono':10,'age_upper_bound_at_arrival_ms':20,'targets':[{'xy_px':[1,2]}],
     'candidates':[1],'unreviewed_region_shapes':[1],'group_shape':1,'flow':1}
  self.assertFalse(enforce_freshness(m,now=10.2));self.assertEqual(m['targets'],[]);self.assertIsNone(m['group_shape'])
  m={'arrival_mono':10,'age_upper_bound_at_arrival_ms':1,'image_available':False}
  self.assertFalse(enforce_freshness(m,now=10.001))
 def test_tracker_resets_and_lost_has_no_coordinates(self):
  e=VisionEngine('opencv');image=np.zeros((30,40,3),np.uint8)
  e.detect=lambda f:([{'label':1,'xy_local':[10,10],'area':9,'score':.5,'score_semantics':'fixture'}],image[:,:,0],np.pad(np.ones((3,3),np.int32),((9,18),(9,28))))
  def h(seq,t,session='a'):return {'seq':seq,'source_time_s':t,'session':session,'source':'FIXTURE','source_time_semantics':'test',
                              'arrival_mono':time.monotonic(),'age_upper_bound_at_arrival_ms':0,'roi_xywh':[0,0,40,30]}
  e.process(image,h(1,1));msg,_,_=e.process(image,h(2,1.03));self.assertTrue(msg['targets'])
  e.detect=lambda f:([],image[:,:,0],np.zeros((30,40),np.int32))
  msg,_,_=e.process(image,h(3,1.06));self.assertTrue(all(t['xy_px'] is None for t in msg['targets']))
  e.process(image,h(20,1.5));self.assertEqual(e.reset_events[-1]['reason'],'time_or_frame_gap')
  e.process(image,h(1,1.6,'b'));self.assertEqual(e.reset_events[-1]['reason'],'source_session_changed')
 def test_mock_rejects_duplicates_and_never_controls(self):
  m=MockConsumer()
  try:
   for _ in range(2):
    m.send({'arrival_mono':time.monotonic(),'age_upper_bound_at_arrival_ms':1,'session':'fixture','frame_id':2,
            'visual_candidate_available':True,'status':'NO_VALID_TARGET'});time.sleep(.03)
  finally:m.close()
  self.assertEqual(sum(r['candidate_valid'] for r in m.rows),1)
  self.assertTrue(all(not r['control_permission'] and not r['bubble_valid'] for r in m.rows))
 def test_disconnect_and_stale_entry_clear_mock(self):
  m=MockConsumer()
  try:
   m.send({'arrival_mono':time.monotonic(),'age_upper_bound_at_arrival_ms':1,'session':'old','frame_id':1,
           'visual_candidate_available':True,'status':'NO_VALID_TARGET'});time.sleep(.03)
   m.send(invalid_state('old',2,'STALE_FRAME_BEFORE_INFERENCE'));time.sleep(.03)
   m.send(invalid_state('old',3,'STREAM_ENDED_OR_RUN_STOPPED'));time.sleep(.03)
  finally:m.close()
  self.assertEqual(len(m.rows),3)
  self.assertTrue(m.rows[0]['candidate_valid'])
  self.assertTrue(all(not r['candidate_valid'] for r in m.rows[1:]))

if __name__=='__main__':unittest.main(verbosity=2)
