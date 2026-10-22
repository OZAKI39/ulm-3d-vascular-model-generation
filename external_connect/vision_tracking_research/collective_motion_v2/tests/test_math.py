import sys,unittest,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,torch,cv2
from scipy.linalg import logm
from geometry import transform,endpoint_six,resize_map,original_transform,apply_transform,affine_residual
from networks import weighted_affine,six_from_normalized_H
from dataset import splat,nonlinear_displacement
from heatmaps import centers_to_heatmap
class MathTests(unittest.TestCase):
 def test_exponential_log_roundtrip(self):
  rng=np.random.default_rng(81)
  for dt in [.01665,.0666,.1]:
   for _ in range(20):
    p=rng.uniform(-1,1,6)*[40,40,.6,1.2,.4,.4]
    H=transform(p,dt,[128,128]);np.testing.assert_allclose(endpoint_six(H,dt,(256,256)),p,atol=1e-10)
 def test_anisotropic_resize(self):
  p=np.array([15,-22,.4,.8,.2,-.1]);H=transform(p,.0666,[128,192]);S=resize_map((384,256),(128,128))
  G=S@H@np.linalg.inv(S);back=original_transform(G,(384,256),(128,128))
  np.testing.assert_allclose(endpoint_six(back,.0666,(384,256)),p,atol=1e-10)
  self.assertGreater(np.max(abs(endpoint_six(G,.0666,(128,128))[2:]-p[2:])),.1)
 def test_geometric_fit_and_differentiable_log(self):
  xy=torch.tensor(np.stack(np.meshgrid(np.linspace(30,226,16),np.linspace(30,226,16)),-1).reshape(-1,2),dtype=torch.float32)
  ps=np.array([[20,-15,.4,.7,.2,-.3],[-30,10,-.5,-1,.3,.2]])
  ds=np.array([apply_transform(transform(p,.1,[128,128]),xy.numpy())-xy.numpy() for p in ps])
  d=torch.tensor(ds,dtype=torch.float32,requires_grad=True)
  H,_=weighted_affine(xy,d,torch.ones(2,len(xy)))
  out=six_from_normalized_H(H,torch.tensor([.1,.1]));np.testing.assert_allclose(out.detach().numpy(),ps,atol=2e-5)
  out.square().sum().backward();self.assertTrue(torch.isfinite(d.grad).all())
 def test_nonlinear_rule_not_shear(self):
  p=np.stack(np.meshgrid(np.linspace(16,240,16),np.linspace(16,240,16)),-1).reshape(-1,2)
  self.assertGreater(affine_residual(p,nonlinear_displacement(p,.1,.1)),1.5)
  self.assertLess(affine_residual(p,np.c_[p[:,1]*.1,p[:,0]*0]),1e-10)
 def test_splat_f10u_reference(self):
  patch=np.outer(np.hanning(21),np.hanning(21)).astype(np.float32)*30
  for xy in [[101.3,120.7],[-1.2,4.1],[256.3,254.1]]:
   a=np.zeros((256,256),np.float32);splat(a,patch,xy)
   ref=-cv2.warpAffine(patch,np.array([[1,0,xy[0]-10],[0,1,xy[1]-10]],np.float32),(256,256),flags=cv2.INTER_LINEAR)
   np.testing.assert_allclose(a,ref,atol=1e-5)
 def test_heat_peak(self):
  h=centers_to_heatmap([[64,64]],(128,128));self.assertAlmostEqual(float(h[64,64]),1,places=3)
if __name__=='__main__':unittest.main()
