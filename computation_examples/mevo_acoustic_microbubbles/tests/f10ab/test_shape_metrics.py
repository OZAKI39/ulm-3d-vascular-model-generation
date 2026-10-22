import numpy as np
from scipy.spatial.transform import Rotation
from scripts.straight_tube_exp1.f10ab.metrics import measure, deformation, geometry
import pytest

def cloud():return np.random.default_rng(24).normal(size=(64,3))*np.array([12,7,4])*1e-6

def test_translation_rotation_not_shape():
 x=cloud();r=Rotation.from_rotvec([.4,-.2,.7]).as_matrix();y=x@r.T+[1e-3,-.002,.004]
 m=measure(x,y)
 assert m['kabsch_rmsd_3d_m']<1e-17
 assert m['similarity_rmsd_3d_m']<1e-17
 assert m['deviatoric_stretch_3d']<1e-12
 assert np.allclose(m['rigid_rotation_3d'],r,atol=1e-12)

def test_uniform_contraction_not_anisotropy():
 x=cloud();m=measure(x,.8*x)
 assert m['kabsch_rmsd_3d_m']>1e-6
 assert m['similarity_rmsd_3d_m']<1e-18
 assert abs(m['isotropic_stretch_3d']-.8)<1e-12
 assert m['deviatoric_stretch_3d']<1e-12

def test_known_stretch_and_shear():
 x=cloud();f=np.array([[1.1,.04,0],[.04,.9,0],[0,0,1]])
 m=measure(x,x@f.T)
 assert np.allclose(m['stretch_tensor_3d'],f,atol=1e-12)
 assert m['shear_norm_3d']==pytest.approx(np.sqrt(2)*.04)
 assert m['nonaffine_residual_3d_m']<1e-18
 assert m['similarity_rmsd_3d_m']>1e-7

def test_nonaffine_rearrangement():
 x=cloud();y=x.copy();y[0]+=[4e-6,-2e-6,3e-6];m=measure(x,y)
 assert m['nonaffine_residual_3d_m']>1e-7

def test_degenerate_rank_and_orientation():
 x=np.zeros((32,3));x[:,0]=np.linspace(-1,1,32)*1e-5
 assert not deformation(x,x)['affine_valid_3d']
 assert np.isnan(deformation(x,x)['affine_deformation_3d']).all()
 sphere=np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])*5e-6
 assert not geometry(sphere)['orientation_valid_3d'].any()

def test_xy_projection_does_not_define_contact():
 x=np.array([[0,0,0],[0,0,10],[10,0,0],[0,10,0]],float)*1e-6
 assert geometry(x)['hard_contact_pairs']==0

def test_fixed_identity_and_invalid_input():
 x=cloud();y=x.copy();y[[0,1]]=y[[1,0]]
 assert measure(x,y)['kabsch_rmsd_3d_m']>1e-7
 with pytest.raises(ValueError):measure(x,x[:-1])

def test_invalid_global_affine_fit_does_not_invent_axis():
 x=cloud();y=x+np.random.default_rng(9).normal(size=x.shape)*1e-5
 m=measure(x,y)
 assert not m['affine_valid_3d']
 assert np.isnan(m['stretch_tensor_3d']).all()
