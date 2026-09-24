from particle_3d.particle7_cases import channel
from particle_3d.injection_admission import FiniteSizeAdmission
import numpy as np
def test_original_oblate(source):
 s,w,c,v=channel(); a=FiniteSizeAdmission(s,source,wall=w,velocity=v); e=source.next_rbc(); e.update(particle_id=1,attempt_count=0)
 p,status,_=a.check(e,[0,0,0],{}); assert p.mode_code==2 and p.position[2]==0
 assert np.array_equal(p.q,e['q']) and np.array_equal(p.axes,[e['geometry'][k] for k in ['a_m','b_m','c_m']])

def test_rbc_free_omega_preserves_p2(source):
 from types import SimpleNamespace
 from particle_3d.injection_population import geometry_from_dict
 from particle_3d.rbc_orientation import angular_velocity,short_axis
 class Field:
  def sample(self,x): return SimpleNamespace(inside_lumen=True,velocity_m_s=np.array([0.,0.,1.]),vorticity_s_inv=np.array([0.,0.,-2.]),velocity_gradient_s_inv=np.array([[0.,2.,0.],[0.,0.,0.],[0.,0.,0.]]))
 s,w,c,v=channel(); a=FiniteSizeAdmission(s,source,wall=w,field=Field()); e=source.next_rbc();e.update(particle_id=1,attempt_count=0)
 p,status,_=a.check(e,[0,0,0],{}); sample=Field().sample([0,0,0])
 expected=angular_velocity(short_axis(e['q']),geometry_from_dict(e['geometry']).jeffery_lambda,sample.velocity_gradient_s_inv,sample.vorticity_s_inv)
 assert np.array_equal(p.omega,expected)
 assert not np.array_equal(p.omega,.5*sample.vorticity_s_inv)
