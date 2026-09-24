import numpy as np
from particle_3d.injection_method_c import MethodCSource
from particle_3d.inlet_flux import InletFluxSampler
from particle_3d.particle7_cases import channel

class FixedDistribution:
    contract_sha256='synthetic-fixed-diameter'
    def __init__(self,d):self.d=d;self.calls=0
    def sample(self,rng):self.calls+=1;return self.d,dict(uniform=float(rng.random()))

def source(distribution,width=3e-6,flux_slope=0.,seed=123):
    sampler,wall,_,_=channel(width=width,length=20e-6)
    sampler=InletFluxSampler(sampler.triangles,1+flux_slope*sampler.triangles[:,:,0]/(width/2))
    return MethodCSource(sampler=sampler,wall=wall,field=None,distribution=distribution,seed=seed,flow_sha256='synthetic-flow',geometry_sha256='synthetic-wall')

