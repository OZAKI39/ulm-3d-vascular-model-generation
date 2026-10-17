"""Instance semantics, component extraction and non-destructive acceptance."""
import unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sam_prompt_geometry import point_instance_requests,point_component,single_bubble_checks
from sam_batch_review import assemble_candidates

class Instances(unittest.TestCase):
 def test_roi_does_not_merge_points_or_enlarge_boxes(self):
  points=[[100,100],[130,130],[300,300]]
  plans=point_instance_requests((500,500),points,[[95,95],[400,400]],[50,50,450,450],16)
  self.assertEqual(len(plans),3)
  for i,p in enumerate(plans):
   self.assertEqual(p['positive_yx'],[points[i]])
   np.testing.assert_allclose(np.array(p['box_yxyx'])[2:]-p['box_yxyx'][:2],16)
  self.assertEqual(plans[0]['excluded_other_point_indices'],[1])
  self.assertEqual(plans[0]['user_negative_indices'],[0])
  self.assertEqual(plans[2]['user_negative_indices'],[])
 def test_duplicate_points_and_outside_roi(self):
  p=point_instance_requests((50,50),[[20,20],[20,21]],[],None,16)
  self.assertEqual(p[0]['duplicate_point_indices'],[1])
  with self.assertRaises(ValueError):point_instance_requests((50,50),[[20,20]],[],[0,0,10,10],16)
 def test_component_keeps_seed_boundary_and_holes(self):
  raw=np.zeros((128,128),bool);raw[60:68,60:68]=True;raw[63,63]=False;raw[5:20,5:20]=True
  before=raw.copy();component,info=point_component(raw,[61,61])
  self.assertEqual(component.sum(),63);self.assertFalse(component[63,63]);self.assertEqual(info['removed_other_component_area_px'],225)
  np.testing.assert_array_equal(raw,before)
  empty,info=point_component(raw,[40,40]);self.assertFalse(empty.any());self.assertFalse(info['seed_in_raw_mask'])
 def test_connected_background_is_rejected_not_cropped(self):
  raw=np.zeros((128,128),bool);raw[20:110,20:110]=True
  component,_=point_component(raw,[64,64]);np.testing.assert_array_equal(component,raw)
  check=single_bubble_checks(component,[[64,64]],[],[10,10,120,120],16)
  self.assertFalse(check['allowed'])
 def test_overlap_rejects_both_without_order_clipping(self):
  masks=np.ones((1,4,4),bool)
  target=lambda y,x:{'crop_yxyx':[y,x,y+4,x+4],'proposals':[{'allowed':True,'reasons':[]}]}
  labels,items,rejected=assemble_candidates([masks,masks],[target(2,2),target(3,3)],[0,0],[True,True],(20,20))
  self.assertFalse(labels.any());self.assertEqual(set(rejected),{0,1})
  labels,items,rejected=assemble_candidates([masks,masks],[target(2,2),target(10,10)],[0,0],[True,True],(20,20))
  self.assertEqual(set(np.unique(labels)),{0,1,2});self.assertEqual(len(items),2)
  occupied=np.zeros((20,20),bool);occupied[2,2]=True
  labels,items,rejected=assemble_candidates([masks],[target(2,2)],[0],[True],(20,20),occupied)
  self.assertEqual(set(rejected),{0});self.assertFalse(labels.any())

if __name__=='__main__':unittest.main(verbosity=2)
