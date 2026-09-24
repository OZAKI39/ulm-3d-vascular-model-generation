import numpy as np
from method_c_test_helpers import source,FixedDistribution
from particle_3d.injection_method_c import sample_flux_weighted_feasible_position

def test_symmetric_flux_and_geometry_no_side_bias():
    s=source(FixedDistribution(1e-6));sampler,mapping,bounds=s.tree.feasible_proposal(.5e-6)
    centers=[]
    for i in range(1,4001):
        e=dict(particle_id=i,species='MB',diameter_m=1e-6,radius_m=.5e-6,q=[1.,0.,0.,0.],diameter_draw_id=[123,i,920,0])
        p,_,_=sample_flux_weighted_feasible_position(e,sampler=sampler,mapping=mapping,checker=s.checker,seed=456,guard=bounds['position_guard'],bounds=bounds)
        centers.append(p)
    p=np.array(centers)
    for axis in [0,1]:
        assert abs(np.mean(p[:,axis]>0)-.5)<5*np.sqrt(.25/len(p))
        assert abs(p[:,axis].mean())<5*(.998e-6)/np.sqrt(3*len(p))
