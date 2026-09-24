import numpy as np
from method_c_test_helpers import source,FixedDistribution
from particle_3d.injection_method_c import sample_flux_weighted_feasible_position

def test_accelerated_conditional_flux_matches_analytic_rectangle():
    s=source(FixedDistribution(1e-6),flux_slope=.6)
    sampler,mapping,bounds=s.tree.feasible_proposal(.5e-6)
    x=[]
    for i in range(1,3501):
        e=dict(particle_id=i,species='MB',diameter_m=1e-6,radius_m=.5e-6,q=[1.,0.,0.,0.],diameter_draw_id=[123,i,920,0])
        p,_,_=sample_flux_weighted_feasible_position(e,sampler=sampler,mapping=mapping,checker=s.checker,seed=123,guard=bounds['position_guard'],bounds=bounds)
        x.append(p)
    x=np.array(x);L=1.5e-6-.5e-6-2e-9
    target=.6/1.5e-6*L**2/3
    se=L/np.sqrt(3*len(x))
    assert abs(x[:,0].mean()-target)<5*se
    assert abs(x[:,1].mean())<5*se
    assert np.max(abs(x[:,:2]))<=L+s.wall.roundoff_m
    assert bounds['feasible_flux_fraction_lower']<=(2*L/3e-6)**2<=bounds['feasible_flux_fraction_upper']
