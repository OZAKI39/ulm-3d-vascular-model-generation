from particle_3d.particle7_cases import channel
from particle_3d.injection_admission import FiniteSizeAdmission
from particle_3d.particle_shapes import Sphere
import pytest
@pytest.mark.parametrize('species',['MB','RBC'])
def test_collision(source,species):
 s,w,c,v=channel(); a=FiniteSizeAdmission(s,source,wall=w,velocity=v); e=source.next_mb() if species=='MB' else source.next_rbc(); e.update(particle_id=2,attempt_count=0)
 p,status,_=a.check(e,[0,0,0],{1:Sphere([0,0,0],2e-6)})
 assert p is None and status=='PAIR_REJECTED'
