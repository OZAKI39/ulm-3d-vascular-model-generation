from stage017_helpers import *
import pytest
from fem3d.adaptive_port import rim_distance,size_function

def test_exact_segment_distance_not_fitted_circle():
 polygon=np.array([[0.,0.],[2.,0.],[2.,1.],[0.,1.]])
 np.testing.assert_allclose(rim_distance([[1,.5],[.1,.6],[2.5,.5]],polygon),[.5,.1,.5])

def test_unit_slope_and_saturation():
 np.testing.assert_allclose(size_function(np.array([0,.2,1.,3.]),.1,.7,1.),[.1,.3,.7,.7])

def test_scale_equivariance_of_size_function():
 d=np.linspace(0,2,11)
 for scale in (.5,1.,4.):
  np.testing.assert_allclose(size_function(d*scale,.1*scale,.8*scale,1.)/scale,size_function(d,.1,.8,1.))
