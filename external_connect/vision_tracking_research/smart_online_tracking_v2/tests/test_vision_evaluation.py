import unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from vision_evaluation import evaluate_instances,assignment,aggregate

class Evaluation(unittest.TestCase):
 def setUp(self):
  self.gt=np.zeros((60,60),int);self.gt[15:20,15:20]=1;self.gt[35:40,35:40]=2
 def objects(self,p):
  out=[]
  for n in np.unique(p):
   if n:
    y,x=np.where(p==n);out.append({'label':int(n),'area':len(x),'xy_local':[x.mean(),y.mean()]})
  return out
 def test_identity(self):
  m=evaluate_instances(self.gt,self.gt,self.objects(self.gt));self.assertEqual(m['center']['f1'],1);self.assertEqual(m['mask']['f1'],1)
 def test_missing_and_spurious(self):
  p=self.gt.copy();p[p==2]=0;p[40:44,15:19]=7
  m=evaluate_instances(self.gt,p,self.objects(p));self.assertEqual([m['center'][k] for k in ('tp','fp','fn')],[1,1,1])
 def test_duplicate_not_two_true_positives(self):
  p=self.gt.copy();p[20:22,15:20]=3
  m=evaluate_instances(self.gt,p,self.objects(p));self.assertEqual(m['center']['tp'],2);self.assertEqual(m['center']['fp'],1)
 def test_truncated_truth_excluded_and_prediction_ignored(self):
  m=evaluate_instances(self.gt,self.gt,self.objects(self.gt),full_gt_counts={1:100,2:25})
  self.assertEqual(m['excluded_gt_ids'],[1]);self.assertEqual(m['center']['tp'],1);self.assertEqual(m['center']['fp'],0)
 def test_centers_and_masks_are_distinct_metrics(self):
  p=np.zeros_like(self.gt);p[17,17]=1;p[37,37]=2
  m=evaluate_instances(self.gt,p,self.objects(p));self.assertEqual(m['center']['f1'],1);self.assertEqual(m['mask']['f1'],0)
 def test_cardinality_first_and_empty(self):
  s=np.array([[1,.01],[.01,0]]);v=np.array([[1,1],[1,0]],bool)
  self.assertEqual(len(assignment(s,v)),2);self.assertEqual(assignment(np.zeros((0,3)),np.zeros((0,3),bool)),[])
 def test_ignore_and_roi_denominator(self):
  p=self.gt.copy();p[40:44,15:19]=3;ignore=p==3
  m=evaluate_instances(self.gt,p,self.objects(p),ignore);self.assertEqual(m['center']['fp'],0)
  a=aggregate([{'image_id':'frame','metrics':m}]*2)
  self.assertEqual(a['source_frame_count'],1);self.assertEqual(a['roi_count'],2)

if __name__=='__main__':unittest.main()
