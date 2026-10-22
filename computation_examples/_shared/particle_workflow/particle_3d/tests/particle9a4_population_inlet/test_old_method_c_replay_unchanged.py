from particle_3d.injection_method_c import MethodCSource,canonical_bytes
from particle_3d.particle7_cases import channel
from particle_3d.population_inlet_p9a4 import generate

def test_legacy_c_replay_preserved(distribution,simple_source):
 sampler,wall,_,_=channel(width=3e-6,length=20e-6)
 old=MethodCSource(sampler=sampler,wall=wall,field=None,distribution=distribution,flow_sha256='synthetic-flow',geometry_sha256='wall',seed=123)
 before=canonical_bytes(old.event(4));generate(simple_source,128)
 assert canonical_bytes(old.event(4))==before
