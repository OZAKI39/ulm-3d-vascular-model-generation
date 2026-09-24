from copy import deepcopy
from particle_3d.particle_shapes import Sphere
from particle_3d.particle7_cases import channel
from particle_3d.injection_admission import FiniteSizeAdmission

def test_blocked_draws_preserve_geometry(source):
 s,w,c,v=channel(); a=FiniteSizeAdmission(s,source,wall=w,velocity=v,guard=2); e=source.next_rbc(); e.update(particle_id=1,attempt_count=0)
 original=deepcopy(e); state=deepcopy(source.rng['RBC_GEOMETRY'].bit_generator.state)
 for _ in range(2): assert a.attempt(e,{2:Sphere([0,0,0],1e-3)}) is None
 assert original['geometry']==e['geometry'] and original['q']==e['q']
 assert state==source.rng['RBC_GEOMETRY'].bit_generator.state and e['attempt_count']==2
