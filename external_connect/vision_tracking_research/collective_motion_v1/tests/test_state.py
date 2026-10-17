import sys,unittest,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from state import *
from affine import empty

class StateTests(unittest.TestCase):
    def state(self):return make_state('TEST','s',1,.1,10,[0,0,20,20],{'group_center_xy':[2,3],'field_validity':{'group_center_xy':True}},empty('NO_SUPPORT'),{},3)
    def test_null_preserved_tensor(self):
        x,mask=tensor_input(self.state());self.assertTrue(np.isnan(x).all());self.assertFalse(mask.any())
    def test_stale_not_current_distribution(self):
        s,r=ConsumerGate().accept(self.state(),11);self.assertFalse(r['fresh']);self.assertIsNone(s['distribution']['group_center_xy'])
    def test_duplicate_and_discontinuity(self):
        g=ConsumerGate();g.accept(self.state(),10.01);s,r=g.accept(self.state(),10.02);self.assertFalse(r['fresh'])
        a=self.state();a['frame_id']=2;a['source_PTS']=.05;s,r=g.accept(a,10.03);self.assertEqual(r['reason'],'PTS_DISCONTINUITY')
    def test_causal_gap_resets(self):
        w=CausalWindow(3);im=np.zeros((8,8));self.assertFalse(w.add(1,.1,im,[]));self.assertFalse(w.add(2,.2,im,[]));self.assertFalse(w.add(4,.3,im,[]));self.assertEqual(len(w.frames),1)
    def test_partial_validity(self):
        s=self.state();s['motion']['valid']['translation_vx']=True;s['motion']['values']['translation_vx']=3.
        x,m=tensor_input(s);self.assertEqual(x[0],3.);self.assertEqual(m.sum(),1);self.assertTrue(np.isnan(x[1:]).all())
    def test_no_new_frame_expires_without_message(self):
        g=ConsumerGate();g.accept(self.state(),10.01);s=g.current(10.2);self.assertIsNone(s['distribution']['group_center_xy']);self.assertEqual(s['motion']['reasons']['translation_vx'],'NO_NEW_FRAME_EXPIRED')

if __name__=='__main__':unittest.main()
