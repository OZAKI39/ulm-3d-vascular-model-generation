import numpy as np
from scipy.spatial.transform import Rotation
from scipy.linalg import block_diag
from particle_3d.planar_wall_hydrodynamics import planar_wall_affine_block

def test_rotation_covariance():
    n=np.array([.3,.4,np.sqrt(.75)]);g=np.array([[3,6,2],[9,1,7],[-4,8,2.]])*100
    u=np.array([.001,.002,-.001,80,20,50.]);r,b,_=planar_wall_affine_block(1e-6,.003,1e-8,n,g,u)
    for A in [np.eye(3),Rotation.from_euler('z',90,degrees=True).as_matrix(),Rotation.from_rotvec([.4,-.7,1.1]).as_matrix()]:
        B=block_diag(A,A);rr,bb,_=planar_wall_affine_block(1e-6,.003,1e-8,A@n,A@g@A.T,B@u)
        D=np.diag([1,1,1,1e6,1e6,1e6])
        np.testing.assert_allclose(D@rr@D,D@B@r@B.T@D,rtol=1e-12,atol=1e-20)
        np.testing.assert_allclose(D@bb,D@B@b,rtol=1e-12,atol=1e-24)
