from types import SimpleNamespace
from particle_3d.particle82a_admission import common_event,method_b
from particle_3d.particle7_cases import channel
from particle_3d.population_inlet_p9a4 import generate

def test_legacy_b_replay_preserved(distribution,simple_source):
 sampler,wall,_,_=channel(width=3e-6,length=20e-6)
 geometry=SimpleNamespace(distances=lambda p:(wall.nearest_center_triangle(p)[1],wall.nearest_center_triangle(p)[1],0.))
 ctx=SimpleNamespace(env=SimpleNamespace(sampler=sampler,wall=wall,field=SimpleNamespace(locate=lambda p:(0,None))),geometry=geometry,distribution=distribution.original)
 event=common_event(41,ctx=ctx);before=method_b(event,ctx)
 generate(simple_source,128)
 assert method_b(event,ctx)==before
