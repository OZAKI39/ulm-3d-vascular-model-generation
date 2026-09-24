import numpy as np
from particle_3d.particle9a_diagnostics import inlet_distance

def test_oriented_plane_positive_is_inside():
    origin=np.array([1.,2.,3.]); inward=np.array([0.,0.,1.])
    assert inlet_distance(origin+2*inward,origin,inward)==2.
    assert inlet_distance(origin-2*inward,origin,inward)==-2.
    velocity=3*inward
    assert velocity@inward>0
    assert inlet_distance(origin+.01*velocity,origin,inward)>0
