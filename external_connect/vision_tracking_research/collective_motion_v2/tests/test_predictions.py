import sys,unittest,numpy as np,torch
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from inference import convert_shape,shape_jacobian,estimate_from_prediction
from vcommon import NAMES
from networks import build_model
class PredictionTests(unittest.TestCase):
 def test_shape_linear_conversion(self):
  p=np.array([12,-8,.4,.7,.2,-.1]);J=shape_jacobian((384,256));np.testing.assert_allclose(convert_shape(p,.0666,(384,256)),J@p,atol=1e-10)
 def test_rejected_joint_gate_nulls(self):
  pred=dict(probability=np.array([.9,.05,.05]),algebra_valid=True,mean=np.ones(6),sigma=np.ones(6),score=2.)
  cal=dict(threshold=None,field_thresholds=[3.]*6,quantiles={})
  est=estimate_from_prediction(pred,cal,'test','abc');self.assertFalse(any(est['valid'].values()));self.assertTrue(all(v is None for v in est['values'].values()))
 def test_finite_model_gradient(self):
  torch.set_num_threads(2)
  for name in ['direct_gray','direct_fused','direct_heat','geo_fused']:
   torch.manual_seed(13);m=build_model(name);x=torch.randn(2,3,2,128,128)
   o=m(x,torch.ones(2,2)*.0333);loss=o['mean'].square().mean()+o['logits'].square().mean();loss.backward()
   self.assertTrue(torch.isfinite(loss));self.assertTrue(all(torch.isfinite(p.grad).all() for p in m.parameters() if p.grad is not None))
if __name__=='__main__':unittest.main()
