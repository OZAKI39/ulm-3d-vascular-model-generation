"""Diameter ranges operate in original pixels; predictions are not reviewed truth."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from annotation_sizes import diameter_range,size_plan,collect_candidates

class DiameterRanges(unittest.TestCase):
 def plan(self,lo=4,hi=24,model='cyto3'):
  return size_plan({'model':model,'diameter_range_px':[lo,hi]})
 def combine(self,raw,plan):return collect_candidates(raw,[0,0,100,100],[200,300,300,400],plan)
 def test_invalid_limits_rejected_instead_of_sorted_or_clipped(self):
  for values in ([3,10],[4,101],[24,8],[16,16],[float('nan'),24],[4,float('inf')],[],[1,2,3]):
   with self.subTest(values=values),self.assertRaises(ValueError):diameter_range(values)
  self.assertEqual(diameter_range([4,100]),[4.,100.])
 def test_three_real_inference_scales_not_only_mean(self):
  self.assertEqual(self.plan(8,32)['inference_diameters_px'],[8,16,32])
  self.assertEqual(len(self.plan(4,100)['inference_diameters_px']),3)
  for model in ('microsam_lm','cellpose4'):
   p=self.plan(model=model);self.assertEqual(p['inference_diameters_px'],[None]);self.assertEqual(p['range_applies_to'],'candidate_filter_only')
 def test_scalar_legacy_is_explicit_and_mixed_arguments_refused(self):
  p=size_plan({'model':'cyto3','diameter_px':16})
  self.assertEqual(p['mode'],'LEGACY_NOMINAL_DIAMETER');self.assertEqual(p['inference_diameters_px'],[16])
  with self.assertRaises(ValueError):size_plan({'model':'cyto3','diameter_px':16,'diameter_range_px':[8,24]})
 def test_small_and_large_objects_preserved_at_original_coordinates(self):
  raw=[np.zeros((100,100),np.uint32) for _ in range(3)]
  raw[0][10:16,10:16]=1;raw[2][60:77,60:77]=1
  labels,objects,rejected=self.combine(raw,self.plan(6,24))
  self.assertEqual(len(objects),2);self.assertEqual(rejected,[])
  self.assertTrue(np.array_equal(labels>0,(raw[0]>0)|(raw[2]>0)))
  self.assertTrue(all(6<=o['equivalent_diameter_px']<=24 for o in objects))
  for obj in objects:
   y,x=np.where(labels==obj['id']);self.assertEqual(obj['centroid_yx'],[float(y.mean()+200),float(x.mean()+300)])
 def test_lower_and_upper_filters_apply_to_all_scales(self):
  raw=[np.zeros((100,100),np.uint32) for _ in range(3)]
  raw[0][10:13,10:13]=1;raw[0][30:39,30:39]=2;raw[0][60:85,60:85]=3
  labels,objects,rejected=self.combine(raw,self.plan(8,20))
  self.assertEqual(len(objects),1);self.assertEqual(objects[0]['raw_id'],2)
  self.assertEqual({r['reasons'][0] for r in rejected},{'below_diameter_range','above_diameter_range'})
 def test_duplicate_predictions_count_once_without_editing_boundary(self):
  raw=[np.zeros((100,100),np.uint32) for _ in range(3)]
  for i,a in enumerate(raw):a[30:42,30:42]=i+4
  labels,objects,rejected=self.combine(raw,self.plan(8,24))
  self.assertEqual(len(objects),1);self.assertTrue(np.array_equal(labels>0,raw[0]>0))
  self.assertEqual(len(rejected),2);self.assertTrue(all(r['reasons']==['duplicate_across_scales'] for r in rejected))
 def test_partial_overlap_is_logged_not_cut_into_fake_instances(self):
  raw=[np.zeros((100,100),np.uint32) for _ in range(3)]
  raw[0][20:30,20:30]=1;raw[1][28:38,28:38]=1
  labels,objects,rejected=self.combine(raw,self.plan(6,24))
  self.assertEqual(len(objects),1);self.assertEqual(np.count_nonzero(labels),100)
  self.assertEqual(rejected[0]['reasons'],['cross_scale_overlap_requires_review'])
  self.assertTrue(any(np.array_equal(labels>0,a>0) for a in raw))
 def test_border_rejected_and_nonround_mask_not_changed(self):
  raw=[np.zeros((100,100),np.uint32) for _ in range(3)]
  raw[0][0:8,10:18]=1;raw[0][30:34,30:50]=2
  labels,objects,rejected=self.combine(raw,self.plan(4,24))
  self.assertEqual(len(objects),1);self.assertTrue(np.array_equal(labels>0,raw[0]==2))
  self.assertIn('roi_border_requires_review',rejected[0]['reasons'])
 def test_legacy_envelope_unchanged(self):
  raw=np.zeros((100,100),np.uint32);raw[10,10:15]=1;raw[30:60,30:60]=2
  p=size_plan({'model':'cyto3','diameter_px':16});labels,objects,rejected=self.combine([raw],p)
  self.assertTrue(np.array_equal(labels>0,raw==1));self.assertEqual(rejected[0]['reasons'],['too_large_for_single_bubble'])
 def test_bad_scale_count_rejected(self):
  with self.assertRaises(ValueError):self.combine([np.zeros((100,100),np.uint32)],self.plan())
 def test_reference_only_preserves_small_large_and_elongated_masks(self):
  raw=[np.zeros((100,100),np.uint32) for _ in range(3)]
  raw[0][10:13,10:13]=1;raw[1][30:55,30:55]=2;raw[2][70:73,10:90]=3
  p=size_plan({'model':'cyto3','diameter_range_px':[8,20],'diameter_filter':'REFERENCE_ONLY'})
  labels,objects,rejected=self.combine(raw,p)
  self.assertEqual(len(objects),3);self.assertEqual(rejected,[])
  self.assertEqual(sum(o['outside_reference_range'] for o in objects),2)
  for o in objects:self.assertTrue(np.array_equal(labels==o['id'],raw[o['scale_index']]==o['raw_id']))
 def test_reference_still_rejects_noise_border_and_roi_dominating_mask(self):
  raw=[np.zeros((100,100),np.uint32) for _ in range(3)]
  raw[0][10:12,10:12]=1;raw[0][0:8,30:38]=2;raw[1][10:90,10:90]=3
  p=size_plan({'model':'cyto3','diameter_range_px':[8,20],'diameter_filter':'REFERENCE_ONLY'})
  labels,objects,rejected=self.combine(raw,p)
  self.assertFalse(np.any(labels));self.assertEqual(objects,[])
  self.assertEqual({v for o in rejected for v in o['reasons']},
                   {'too_small','roi_border_requires_review','dominates_roi_requires_review'})
 def test_reference_policy_for_sam_models_is_advisory(self):
  for model in ('cellpose4','microsam_lm'):
   p=size_plan({'model':model,'diameter_range_px':[8,20],'diameter_filter':'REFERENCE_ONLY'})
   self.assertEqual(p['inference_diameters_px'],[None]);self.assertEqual(p['range_applies_to'],'advisory_only')
   raw=np.zeros((100,100),np.uint32);raw[30:55,30:55]=1
   labels,objects,rejected=self.combine([raw],p)
   self.assertEqual(len(objects),1);self.assertTrue(objects[0]['outside_reference_range'])
   self.assertTrue(np.array_equal(labels,raw));self.assertEqual(rejected,[])
 def test_old_range_requests_stay_strict_and_invalid_policy_is_rejected(self):
  self.assertEqual(self.plan()['diameter_filter'],'STRICT_RANGE')
  for request in ({'model':'cyto3','diameter_range_px':[8,20],'diameter_filter':'INVALID'},
                  {'model':'cyto3','diameter_px':16,'diameter_filter':'REFERENCE_ONLY'}):
   with self.assertRaises(ValueError):size_plan(request)

if __name__=='__main__':unittest.main(verbosity=2)
