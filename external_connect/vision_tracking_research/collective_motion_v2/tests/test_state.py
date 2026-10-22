import sys,time,json,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import jsonschema,numpy as np
from vcommon import ROOT,NAMES
from legacy import V1,MockConsumer,empty,CausalWindow,distribution
from state_v2 import new_state,to_v1,VersionedGate

class StateTests(unittest.TestCase):
 def state(self,fid=1,status='NON_AFFINE_MOTION'):
  e=empty('NON_AFFINE','unit');e['observability_status']=status;d,_=distribution([],(256,256))
  return new_state('SYNTHETIC_MOTION_TEST','test',fid,fid*.0333,time.monotonic(),[0,0,256,256],d,e,{},dict(size=3),model_metadata=dict(name='unit',weight_sha256='0'*64,calibration_version='unit'))
 def test_schema_and_adapter(self):
  s=self.state();jsonschema.validate(s,json.loads((ROOT/'configs/collective_state_v2.schema.json').read_text()))
  with self.assertRaises(jsonschema.ValidationError):jsonschema.validate(s,json.loads((V1/'configs/collective_state.schema.json').read_text()))
  old=to_v1(s);jsonschema.validate(old,json.loads((V1/'configs/collective_state.schema.json').read_text()))
  self.assertTrue(all(x is None for x in old['motion']['values'].values()))
  s['motion']['values']['translation_vx']=0
  with self.assertRaises(jsonschema.ValidationError):jsonschema.validate(s,json.loads((ROOT/'configs/collective_state_v2.schema.json').read_text()))
 def test_gap_duplicate_ttl(self):
  gate=VersionedGate();s=self.state();_,r=gate.accept(s);self.assertTrue(r['fresh'])
  _,r=gate.accept(self.state(3));self.assertFalse(r['fresh'])
  _,r=gate.accept(self.state(3));self.assertFalse(r['fresh'])
  s=gate.current(time.monotonic()+1);self.assertTrue(s['stale_or_discontinuous']);self.assertIsNone(s['distribution']['visible_count'])
 def test_real_mock_consumer(self):
  c=MockConsumer()
  try:
   c.send(to_v1(self.state()));time.sleep(.05);self.assertTrue(c.rows);self.assertTrue(c.rows[-1]['fresh']);time.sleep(.12);self.assertTrue(c.watchdog_events)
  finally:c.close()
 def test_window_causality(self):
  w=CausalWindow(3);a=np.zeros((8,8),np.uint8)
  for i in range(3):ready=w.add(i,i*.0333,a,[])
  self.assertTrue(ready);self.assertFalse(w.add(5,.1665,a,[]));self.assertEqual(w.last_reset,'FRAME_GAP_OR_PTS_DISCONTINUITY')
if __name__=='__main__':unittest.main()
