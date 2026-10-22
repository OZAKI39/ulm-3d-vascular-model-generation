import sys,time,json,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wcommon import *
from benchmark import apply_gate
from motion_models import optimal_transport,sample_field
import jsonschema
from scipy.linalg import expm
class Contracts(unittest.TestCase):
 def test_logm_finite_nonsquare(self):
  center=np.array([128,192]);A=np.array([[.12,-.23],[.31,-.08]]);v=np.array([18.,-12.]);dt=.0666;H=expm(old.generator(v,A,center)*dt);expected=[18,-12,.04,.54,.1,.04]
  np.testing.assert_allclose(old.parameters(H[:2,:2],H[:2,2],dt,center),expected,atol=1e-10)
 def test_field_coordinates(self):
  yy,xx=np.mgrid[:384,:256];field=np.dstack([xx*2,yy*3]).astype(np.float32);xy=np.array([[21,300],[100,19]],np.float32)
  np.testing.assert_allclose(sample_field(field,xy),xy*[2,3]);self.assertTrue(np.isnan(sample_field(field,[[-10,0]])).all())
 def state(self,fid=0):
  d,_=old.distribution([],(384,256));e=old.empty('WARMUP','PIV');return new_state('PACED_RECORDING','unit',fid,fid*.0333,time.monotonic(),[650,300,256,384],d,e,{},dict(frame_ids=[fid]),model_metadata=dict(name='PIV',weight_sha256='NO_LEARNED_WEIGHTS',calibration_version='F10U frozen gate'))
 def test_schema_null_and_v1(self):
  s=self.state();v2=json.loads((ROOT.parent/'collective_motion_v2/configs/collective_state_v2.schema.json').read_text());v1=json.loads((ROOT.parent/'collective_motion_v1/configs/collective_state.schema.json').read_text());jsonschema.validate(s,v2);jsonschema.validate(to_v1(s),v1)
  s['motion']['values']['translation_vx']=0
  with self.assertRaises(jsonschema.ValidationError):jsonschema.validate(s,v2)
 def test_gap_stale_and_window(self):
  g=VersionedGate();self.assertTrue(g.accept(self.state())[1]['fresh']);s,r=g.accept(self.state(2));self.assertFalse(r['fresh']);self.assertTrue(s['stale_or_discontinuous']);self.assertIsNone(g.current(time.monotonic()+1)['distribution']['visible_count'])
  w=old.CausalWindow(3);a=np.zeros((4,4),np.uint8)
  for i in range(3):ready=w.add(i,i*.0333,a,[])
  self.assertTrue(ready);self.assertFalse(w.add(4,.1332,a,[]))
 def test_transport_translation_and_unequal_sets(self):
  rng=np.random.default_rng(42);xy=rng.uniform(20,235,(40,2));shift=np.array([1.2,-.8]);objects=lambda z:[dict(xy=p,area=20,peak=30) for p in z]
  a=objects(xy);b=objects(xy+shift);im=np.zeros((256,256),np.uint8);r=optimal_transport([im]*3,[a,a,b],[0,.0333,.0666],method='unbalanced',distance_scale=8.)
  self.assertTrue(all(r['estimate']['valid'].values()));np.testing.assert_allclose(r['estimate']['raw_estimate'][:2],shift/.0666,atol=.1)
  outliers=objects(np.r_[xy+shift,[[8,8],[240,240]]]);r=optimal_transport([im]*3,[a,a,outliers],[0,.0333,.0666],method='unbalanced',distance_scale=8.)
  self.assertFalse(any(r['estimate']['valid'].values()));self.assertEqual(r['estimate']['reasons']['translation_vx'],'OT_NUMERICAL_FAILURE')
 def test_insufficient_support_never_zero(self):
  a=[dict(xy=[10.,10.],area=20,peak=30)];im=np.zeros((32,32),np.uint8);r=optimal_transport([im]*3,[a]*3,[0,.0333,.0666]);self.assertFalse(any(r['estimate']['valid'].values()));self.assertTrue(all(x is None for x in r['estimate']['values'].values()))
 def test_gate_cannot_fill_invalid(self):
  e=old.empty('NO_SUPPORT','unit');a=apply_gate(e,1,0);self.assertFalse(any(a['valid'].values()));b=apply_gate(e,1,None);self.assertTrue(all(x is None for x in b['values'].values()))
 def test_mock_consumer_watchdog(self):
  c=old.MockConsumer()
  try:c.send(to_v1(self.state()));time.sleep(.04);self.assertTrue(c.rows);time.sleep(.13);self.assertTrue(c.watchdog_events)
  finally:c.close()
if __name__=='__main__':unittest.main()
