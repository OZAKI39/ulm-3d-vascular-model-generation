"""Software geometry checks, no artificial expert labels in the dataset."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sam_prompt_geometry import prompt_region,proposal_checks

class PromptGeometry(unittest.TestCase):
 def test_original_coordinate_roundtrip_and_translation(self):
  p=np.array([[200.,300.]])
  crop,box=prompt_region((600,800),p,[],[190,290,210,310])
  np.testing.assert_equal(crop,[136,236,264,364])
  np.testing.assert_equal((p-crop[:2])+crop[:2],p)
  c2,b2=prompt_region((600,800),p+[20,40],[],np.array([190,290,210,310])+[20,40,20,40])
  np.testing.assert_equal(c2-crop,[20,40,20,40]);np.testing.assert_equal(b2-box,[20,40,20,40])
 def test_box_only_and_edge_crop(self):
  c,b=prompt_region((200,300),[],[],[0,0,20,20]);np.testing.assert_equal(c[:2],[0,0])
  with self.assertRaises(ValueError):prompt_region((200,300),[],[])
  with self.assertRaises(ValueError):prompt_region((200,300),[[np.nan,1]],[])
  with self.assertRaises(ValueError):prompt_region((200,300),[[30,40]],[],[0,0,20,20])
 def test_channel_rejected_without_silhouette_changes(self):
  m=np.ones((128,128),bool);original=m.copy()
  r=proposal_checks(m,[[64,64]],[],[56,56,72,72]);self.assertFalse(r['allowed']);np.testing.assert_equal(m,original)
 def test_prompt_contradiction_rejected(self):
  m=np.zeros((128,128),bool);m[60:68,60:68]=True
  self.assertTrue(proposal_checks(m,[[64,64]],[[75,75]],[56,56,72,72])['allowed'])
  self.assertFalse(proposal_checks(m,[[64,64]],[[65,65]],[56,56,72,72])['allowed'])
  self.assertFalse(proposal_checks(m,[[80,80]],[],[56,56,72,72])['allowed'])
 def test_far_negative_point_is_not_silently_ignored(self):
  with self.assertRaises(ValueError):prompt_region((600,800),[[200,300]],[[500,600]])

if __name__=='__main__':unittest.main(verbosity=2)
