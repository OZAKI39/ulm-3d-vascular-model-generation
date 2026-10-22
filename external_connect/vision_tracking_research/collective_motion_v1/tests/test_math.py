import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from scipy.linalg import expm
from affine import parameters,fit_affine
from methods import piv,cpd

class MathTests(unittest.TestCase):
    def test_log_exact_finite_rotation_translation(self):
        A=np.array([[.1,-.7],[.4,-.2]]);v=np.array([13.,-8.]);c=np.array([128.,128.]);L=np.zeros((3,3));L[:2,:2]=A;L[:2,2]=v-A@c
        H=expm(L*.3);expected=[13,-8,-.1,1.1,.15,-.15]
        np.testing.assert_allclose(parameters(H[:2,:2],H[:2,2],.3,c),expected,atol=1e-9)
    def test_collinear_rejected(self):
        xy=np.c_[np.arange(20)*10,np.ones(20)*80]
        self.assertEqual(fit_affine(xy,np.ones_like(xy),.1,(256,256),'test')['observability_status'],'NOT_OBSERVABLE')
    def test_covariance_rotation_nonidentifiability(self):
        cov=np.eye(2);A=np.array([[0,-1],[1,0]])
        np.testing.assert_array_equal(A@cov+cov@A.T,np.zeros((2,2)))
    def test_official_piv_image_y_down(self):
        import cv2
        rng=np.random.default_rng(42);im=cv2.GaussianBlur(rng.uniform(0,255,(256,256)).astype(np.float32),(0,0),1)
        to=cv2.warpAffine(im,np.array([[1,0,2],[0,1,3]],np.float32),(256,256))
        r=piv(im,to,.1,'raw');vec=np.array(r['vectors']);med=np.median(vec[vec[:,4]>0,2:4],0)
        np.testing.assert_allclose(med,[2,3],atol=.15)
    def test_cpd_row_column_convention(self):
        from pycpd import AffineRegistration
        rng=np.random.default_rng(88);x=rng.normal(size=(40,2))*[2,1];F=np.array([[1.01,.02],[-.03,.99]]);b=np.array([.1,-.08]);y=x@F.T+b
        r=AffineRegistration(X=y,Y=x,max_iterations=100,tolerance=1e-8,w=0);_,(B,t)=r.register()
        np.testing.assert_allclose(B.T,F,atol=1e-5);np.testing.assert_allclose(t.ravel(),b,atol=1e-5)
    def test_multipass_not_replaced_vectors(self):
        import cv2
        rng=np.random.default_rng(43);im=cv2.GaussianBlur(rng.uniform(0,255,(256,256)).astype(np.float32),(0,0),1)
        to=cv2.warpAffine(im,np.array([[1,0,2],[0,1,3]],np.float32),(256,256))
        r=piv(im,to,.1,'local',multipass=True)
        self.assertGreater(r['valid_vector_count'],8)
        vec=np.array(r['vectors']);np.testing.assert_allclose(np.median(vec[vec[:,4]>0,2:4],0),[2,3],atol=.3)

if __name__=='__main__':unittest.main()
