from method_c_test_helpers import source,FixedDistribution

def test_global_rejection_and_certified_capacity():
    class D(FixedDistribution):
        def sample(self,rng):
            self.calls+=1
            return (3.5e-6 if self.calls==1 else 1e-6),dict(uniform=float(rng.random()))
    s=source(D(1e-6));lo,hi=s.capacity['D_geometry_max_bracket_m']
    assert lo<=3e-6-4e-9+2*s.wall.roundoff_m<=hi
    e=s.event(1)
    assert e['diameter_global_rejections']==1
    assert e['source_diameter_draws'][0]['status']=='NO_FEASIBLE_INLET_POSITION_FOR_SIZE'
    assert e['diameter_draw_id'][-1]==1 and e['diameter_m']==1e-6
