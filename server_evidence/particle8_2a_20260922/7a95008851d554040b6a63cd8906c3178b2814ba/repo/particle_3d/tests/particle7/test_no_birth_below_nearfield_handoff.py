from particle_3d.particle7_cases import channel
from particle_3d.injection_admission import FiniteSizeAdmission
from particle_3d.particle_shapes import Sphere

def test_wall_and_pair(source):
 s,w,c,v=channel(); a=FiniteSizeAdmission(s,source,wall=w,velocity=v); e=source.next_mb(); e.update(particle_id=2,attempt_count=0); r=e['radius_m']
 p,status,_=a.check(e,[100e-6-r-1e-9,0,0],{}); assert p is None and status=='WALL_NEARFIELD_REJECTED'
 p,status,_=a.check(e,[0,0,0],{1:Sphere([r+1e-6+1e-9,0,0],1e-6)})
 assert p is None and status=='PAIR_NEARFIELD_REJECTED'
