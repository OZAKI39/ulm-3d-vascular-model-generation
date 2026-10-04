"""Diagnostic-only exact cylinder queries bound into unchanged production bytecode."""
from dataclasses import dataclass
from types import SimpleNamespace
import numpy as np
from particle_3d.particle_shapes import Sphere,roundoff_length
from particle_3d.wall_gap import WallGapResult
from particle_3d.kinematic_contact import ContactConstraint
from particle_3d.particle6_stepper import bind_query_dependency as bind
from particle_3d import particle65_motion as p65
from particle_3d import nearfield_handoff as handoff
from particle_3d.population_v1 import components
from particle_3d.physical_time_refinement import refine_interval
from particle_3d.nearfield_regularization import NearFieldRegularizationV1
from particle_3d.particle4_motion import ParticleWorldState

@dataclass(frozen=True)
class ExactCylinder:
 radius_m:float
 length_m:float
 open_boundary_topology:object=None

def exact_gap(s,wall,*,inside_lumen=None):
 assert type(s) is Sphere
 r=np.linalg.norm(s.center_m[:2]);assert r>0
 outward=np.r_[s.center_m[:2]/r,0.];n=-outward
 wp=s.center_m+(wall.radius_m-r)*outward;pp=s.center_m+s.radius_m*outward
 gap=float(wall.radius_m-r-s.radius_m)
 pad=roundoff_length(s.center_m,s.radius_m,wall.radius_m,wall.length_m)
 return WallGapResult(gap,'SEPARATED' if gap>pad else 'PENETRATING' if gap< -pad else 'TOUCHING',
  inside_lumen,-1,'EXACT_CYLINDER',wp,pp,n,abs(gap),s.mode,1,pad)

def exact_constraints(shapes,wall,policy):
 contacts=[];rows=[]
 for pid,s in sorted(shapes.items()):
  g=exact_gap(s,wall);lower=policy.lower_handoff_gap(policy.reference_length(s))['h_lower_m']
  if abs(g.gap_m)<=g.roundoff_m:
   contacts.append(ContactConstraint.wall(pid,g));rows.append(dict(contact_type='GEOMETRIC_HARD_CONTACT',particle_id=pid))
  if abs(g.gap_m-lower)<=g.roundoff_m:
   c=ContactConstraint(('CONTINUUM_HANDOFF_CONTACT','WALL',pid,'EXACT_CYLINDER'),pid,None,g.normal_inward,g.particle_point_m)
   contacts.append(c);rows.append(dict(contact_type='CONTINUUM_HANDOFF_CONTACT',canonical_id=c.canonical_id,h_geom_m=g.gap_m,
    h_lower_m=lower,g_nf_m=g.gap_m-lower,roundoff_m=g.roundoff_m,normal=c.normal.tolist()))
 return contacts,rows

def exact_certificate(start,end,wall,lower):
 # Norm of affine xy path is convex: maximum over a segment is at an endpoint.
 # Thus minimum gap of the SAME held-velocity path is exactly min(endpoint gaps).
 a=exact_gap(start,wall);b=exact_gap(end,wall);minimum=min(a.gap_m,b.gap_m);pad=max(a.roundoff_m,b.roundoff_m)
 return minimum-lower>=-pad,dict(proof='EXACT_CYLINDER_LINEAR_PATH_MINIMUM',minimum_h_geom_m=minimum,
  h_lower_m=lower,minimum_g_nf_m=minimum-lower,roundoff_m=pad)

def exact_swept(start,end,wall,**kwargs):
 safe,record=exact_certificate(start,end,wall,0.)
 return safe,record['proof'],None

def exact_first_handoff(shapes,velocities,dt,wall,policy):
 events=[]
 for pid,s in sorted(shapes.items()):
  lower=policy.lower_handoff_gap(policy.reference_length(s))['h_lower_m'];v=velocities[pid]
  first=exact_gap(s,wall);last=exact_gap(s.moved(s.center_m+dt*v),wall);pad=max(first.roundoff_m,last.roundoff_m)
  if last.gap_m-lower>=-pad or first.gap_m-lower<=pad:continue
  lo=0.;hi=dt
  for iteration in range(64):
   mid=lo+(hi-lo)/2
   if mid==lo or mid==hi:break
   value=exact_gap(s.moved(s.center_m+mid*v),wall).gap_m-lower
   if value>=0:lo=mid
   else:hi=mid
   if 0<=value<=pad/4:lo=mid;break
  value=exact_gap(s.moved(s.center_m+lo*v),wall).gap_m-lower
  if 0<=value<=pad:
   events.append(dict(kind='WALL',particle_id=pid,triangle_id=-1,time_from_start_s=lo,g_nf_at_root_m=value,
    roundoff_m=pad,root_iterations=iteration+1,geometry='EXACT_CYLINDER'))
 return min(events,key=lambda e:(e['time_from_start_s'],e['particle_id'])) if events else None

class DiagnosticProductionDynamics:
 def __init__(self,field,wall,mu):
  self.field=field;self.wall=wall;self.mu=mu;self.policy=NearFieldRegularizationV1()
  self.interactions=bind(handoff.exact_interactions,dict(wall_gap=exact_gap,all_pairs=lambda s:[]))
  self.require_initial=bind(handoff.require_admissible_initial,dict(exact_interactions=self.interactions))
  self.component_assembly=bind(components.assemble_components,dict(exact_interactions=self.interactions,wall_gap=exact_gap))
  # assemble_components uses the original p65 assembly, planar-wall force and component solver.
  def assembly(shapes,free,mu,wall=None,**options):
   return self.component_assembly(shapes,free,mu,wall,pairs_enabled=False,
     gradient_at=lambda x:self.field.sample(x).velocity_gradient_s_inv,**options)
  self.assembly=assembly
  self.factory=bind(p65.v1_trial,dict(assemble_v1=assembly,exact_interactions=self.interactions,handoff_constraints=exact_constraints,
   first_handoff_event=exact_first_handoff,wall_handoff_certificate=exact_certificate,
   partitioned_wall_handoff_certificate=exact_certificate,swept_clearance_certificate=exact_swept,wall_gap=exact_gap,
   solve_resistance=components.solve_components,contact_jacobian=components.sparse_jacobian,all_pairs=lambda s:[]))
  assert self.factory.__code__ is p65.v1_trial.__code__
 def provider(self,pid,s,t):
  q=self.field.sample(s.center_m)
  if not q.inside_lumen:raise ValueError('CENTER_OUTSIDE_AUTHORITATIVE_FIELD')
  return q.velocity_m_s,.5*q.vorticity_s_inv
 def rhs(self,shapes,t):
  free={i:np.concatenate(self.provider(i,s,t)) for i,s in shapes.items()}
  a=self.assembly(shapes,free,self.mu,self.wall,policy=self.policy)
  cs,rows=exact_constraints(shapes,self.wall,self.policy)
  result=components.solve_components(a,particles=shapes,constraints=cs)
  return result.velocity,rows
