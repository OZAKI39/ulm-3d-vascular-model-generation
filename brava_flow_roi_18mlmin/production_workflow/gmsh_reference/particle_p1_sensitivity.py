"""Four diagnostic exact-vs-P1 scenarios; production CPU bytecode, no CFD trajectory."""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
import json,sys,time,hashlib,traceback,argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;OUT=HERE/'error_diagnosis'
contract=json.loads((ROOT/'PARTICLE_CFD_INPUT_CONTRACT.json').read_text())
source=Path(contract['source_server_root'] if str(ROOT).startswith('/workspace') else contract['source_local_root'])
sys.path[:0]=[str(source.parent),str(ROOT/'scripts'),str(HERE)]
from particle_3d.field import FrozenFEMField,FlowSample
from particle_3d.particle_shapes import Sphere
from particle_3d.physical_time_refinement import refine_interval
from particle_3d.particle4_motion import ParticleWorldState
from poiseuille import reference
from exact_cylinder_particle import ExactCylinder,exact_gap,DiagnosticProductionDynamics
PLAN=json.loads((HERE/'GMSH_REFERENCE_MESH_V2_PLAN.json').read_text())
PROTO=json.loads((ROOT/'manufactured_pipe/PROTOCOL.json').read_text())
R=PLAN['radius_m'];L=PLAN['length_m'];MU=PROTO['mu_Pa_s'];MEAN=PROTO['mean_velocity_m_s']
DT=.000125;T_END=.2;TRANSIT=10e-6
THETA=2*np.pi*PROTO['azimuth_offset']/PROTO['azimuth_count']
SCENARIOS=[('A',.375e-6,2e-9),('B',.964365028029e-6,2e-9),('C',.964365028029e-6,2e-6),('D',.964365028029e-6,10e-6)]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def put(p,d):
 p=Path(p);q=p.with_suffix(p.suffix+'.pending');q.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n');q.replace(p)
def exact(x):return reference(np.asarray(x),R,L,MEAN,MU)
class AnalyticField:
 def sample(self,x):
  q=exact(x);g=q['G_s_inv'];curl=np.array([g[2,1]-g[1,2],g[0,2]-g[2,0],g[1,0]-g[0,1]])
  return FlowSample(q['u_m_s'],float(q['p_Pa']),g,curl,.5*(g+g.T),bool(np.linalg.norm(x[:2])<=R and 0<=x[2]<=L),-1)
def make_field(kind):
 if kind=='EXACT_ANALYTIC_FIELD':return AnalyticField()
 assert kind=='SAME_MESH_P1_ANALYTIC_FIELD'
 with np.load(ROOT/'manufactured_pipe/R1_GMSH_V2/cfd_gpu/SV_MESH/mesh_arrays.npz') as d:p=d['points_m'];t=d['tetra']
 q=exact(p)
 return FrozenFEMField(p,t,q['u_m_s'],q['p_Pa'])
def run_one(spec):
 name,a,gap,kind=spec;start=time.monotonic();field=make_field(kind);wall=ExactCylinder(R,L)
 dynamics=DiagnosticProductionDynamics(field,wall,MU)
 center=np.array([(R-a-gap)*np.cos(THETA),(R-a-gap)*np.sin(THETA),L/2])
 shape=Sphere(center,a);dynamics.require_initial({0:shape},wall,dynamics.policy,MU)
 vel,_=dynamics.rhs({0:shape},0.);world=ParticleWorldState(0.,{0:shape},{0:vel[:3]},{0:vel[3:]},[],{0:exact_gap(shape,wall)},{},{})
 ts=[0.];xs=[center.copy()];vs=[vel[:3].copy()];rhs=[vel[:3].copy()];gaps=[exact_gap(shape,wall).gap_m];cats=[];cells=[]
 records=[];ledger=[];attempts=[]
 def accepted(w):
  s=w.shapes[0];t=w.time_s;u,_=dynamics.rhs(w.shapes,t);h=exact_gap(s,wall)
  assert h.gap_m>=-h.roundoff_m
  category=dynamics.policy.interaction_state(h.gap_m,a,dynamics.policy.lower_handoff_gap(a)['h_lower_m'],h.roundoff_m)
  ts.append(t);xs.append(s.center_m.copy());vs.append(w.velocities[0].copy());rhs.append(u[:3].copy());gaps.append(h.gap_m);cats.append(category)
  cells.append(field.sample(s.center_m).tetra_id)
 trial=dynamics.factory(dynamics.provider,MU,wall,dynamics.policy,on_accept=accepted,attempts=attempts)
 status='OBSERVATION_END';failure=None
 try:
  for k in range(1,round(T_END/DT)+1):
   world,_=refine_interval(world,k*DT,trial,ledger=ledger)
 except Exception as exc:
  status='NUMERICAL_FAILURE';failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
 arr=np.asarray(xs);ts=np.asarray(ts);rv=np.asarray(rhs);vv=np.asarray(vs);gg=np.asarray(gaps)
 assert np.all(np.isfinite(arr)) and np.all(np.diff(ts)>0)
 z=arr[:,2]-center[2];crossing=np.flatnonzero(z>=TRANSIT)
 crossing_time=None
 if len(crossing):
  k=int(crossing[0]);crossing_time=float(ts[k-1]+(TRANSIT-z[k-1])/(z[k]-z[k-1])*(ts[k]-ts[k-1]))
 path=OUT/f'particle_{name}_{kind}.npz';np.savez_compressed(path,time_s=ts,position_m=arr,step_velocity_m_s=vv,RHS_velocity_m_s=rv,gap_m=gg,cell=np.array([-1]+cells))
 initial_gap=exact_gap(shape,wall)
 initial_category=dynamics.policy.interaction_state(initial_gap.gap_m,a,dynamics.policy.lower_handoff_gap(a)['h_lower_m'],initial_gap.roundoff_m)
 categories=sorted(set([initial_category]+cats));contact=[c for c in categories if 'CONTACT' in c]
 result=dict(case=name,field=kind,DIAGNOSTIC_ONLY=True,status=status,failure=failure,radius_m=a,initial_gap_requested_m=gap,initial_gap_actual_m=initial_gap.gap_m,
  initial_position_m=center.tolist(),initial_RHS_velocity_m_s=vel[:3].tolist(),dt_max_s=DT,observation_end_s=float(ts[-1]),
  accepted_steps=len(ts)-1,rejected_trials=sum(not r['accepted'] for r in attempts),maximum_refinement_depth=max((r['depth'] for r in ledger),default=0),
  minimum_gap_m=float(gg.min()),maximum_gap_m=float(gg.max()),categories=categories,contact_categories=contact,
  transit_distance_m=TRANSIT,transit_time_s=crossing_time,transit_censored=crossing_time is None,
  axial_displacement_m=float(z[-1]),radial_drift_max_m=float(np.max(np.abs(np.linalg.norm(arr[:,:2],axis=1)-np.linalg.norm(center[:2])))),
  wall_time_s=time.monotonic()-start,data_sha256=sha(path),data_path=str(path.relative_to(ROOT)))
 put(OUT/f'particle_{name}_{kind}.json',result)
 return result

def metrics(a,b):
 with np.load(ROOT/a['data_path']) as f:aa={k:f[k] for k in f.files}
 with np.load(ROOT/b['data_path']) as f:bb={k:f[k] for k in f.files}
 end=min(aa['time_s'][-1],bb['time_s'][-1]);tt=np.linspace(0,end,round(end/DT)*4+1)
 # Predeclared linear interpolation at COMMON physical times; positions are affine per accepted step.
 def align(data,key):return np.column_stack([np.interp(tt,data['time_s'],data[key][:,j]) for j in range(3)])
 xa=align(aa,'position_m');xb=align(bb,'position_m');va=align(aa,'RHS_velocity_m_s');vb=align(bb,'RHS_velocity_m_s')
 sa=np.linalg.norm(va,axis=1);sb=np.linalg.norm(vb,axis=1)
 speed=float(np.sqrt(np.trapezoid((sb-sa)**2,tt)/np.trapezoid(sa**2,tt)))
 axial=float(abs((xb[-1,2]-xb[0,2])/(xa[-1,2]-xa[0,2])-1))
 radial=np.abs(np.linalg.norm(xa[:,:2],axis=1)-np.linalg.norm(xb[:,:2],axis=1))
 ga=np.interp(tt,aa['time_s'],aa['gap_m']);gb=np.interp(tt,bb['time_s'],bb['gap_m'])
 transit=abs(b['transit_time_s']/a['transit_time_s']-1) if a['transit_time_s'] and b['transit_time_s'] else None
 contact=a['contact_categories']==b['contact_categories']
 # Persistent radial deviation: exceeds 0.05 um for every common sample in final quarter.
 persistent=bool(np.all(radial[3*len(radial)//4:]>.05e-6))
 result=dict(case=a['case'],speed_history_relative_L2=speed,axial_displacement_relative_error=axial,transit_relative_error=transit,
  maximum_radial_separation_m=float(radial.max()),maximum_gap_difference_m=float(np.max(abs(gb-ga))),contact_categories_equal=contact,
  persistent_radial_deviation=persistent,common_grid_samples=len(tt),common_grid_end_s=float(end),
  initial_RHS_velocity_relative_error=float(np.linalg.norm(np.array(b['initial_RHS_velocity_m_s'])-a['initial_RHS_velocity_m_s'])/np.linalg.norm(a['initial_RHS_velocity_m_s'])),
  exact_status=a['status'],P1_status=b['status'])
 if a['status']!='OBSERVATION_END' or b['status']!='OBSERVATION_END' or transit is None:result['classification']='INCOMPLETE_DIAGNOSTIC'
 elif max(speed,axial,transit)>.05 or not contact or persistent:result['classification']='P1_FIELD_TRAJECTORY_IMPACT_LARGE'
 elif max(speed,axial,transit)<=.02 and radial.max()<=.05e-6:result['classification']='P1_FIELD_TRAJECTORY_IMPACT_SMALL'
 elif max(speed,axial,transit)<=.05 and contact and radial.max()<=.05e-6:result['classification']='P1_FIELD_TRAJECTORY_IMPACT_MODERATE'
 else:result['classification']='INDETERMINATE_RADIAL_CRITERION'
 np.savez_compressed(OUT/f'comparison_{a["case"]}.npz',time_s=tt,exact_position_m=xa,P1_position_m=xb,exact_speed_m_s=sa,P1_speed_m_s=sb,exact_gap_m=ga,P1_gap_m=gb)
 return result

def main():
 OUT.mkdir(exist_ok=True);start=time.monotonic()
 for name,record in contract['source_identity'].items():assert sha(source/name)==record['server_sha256' if str(ROOT).startswith('/workspace') else 'local_sha256'],name
 protocol=dict(DIAGNOSTIC_ONLY=True,CFD_field_used=False,CFD_runs=0,mesh_changes=False,dt_max_s=DT,observation_horizon_s=T_END,
  transit_definition='First axial advance of 10 um from z=L/2, piecewise linear crossing time; not whole-pipe transit. If unreached, report censored, no extrapolation.',
  common_time_grid='Uniform 31.25 us through common accepted horizon; linear interpolation of instantaneous RHS speed and accepted position at physical times; no step-index matching.',
  geometry='Exact finite-radius sphere gap R-hypot(x,y)-a and inward normal; exact linear-path clearance; same in both arms',
  initial_azimuth_rad=THETA,initial_azimuth_rule='First frozen probe azimuth, chosen before running; no favorable point selection',
  backend='Existing production CPU P6.5/P9A component assembly/solve/v1_trial/refine_interval; exact analytic sampler is CPU-only; four independent CPU workers; no new CUDA kernel',
  production_physics_modified=False,maximum_worker_count=4,scenario_count=4,field_count=2,
  persistent_radial_definition='Radial separation >0.05 um at every common sample in final quarter of observation',
  limitations='Single azimuth, 0.2 s local diagnostic, one dt only; does not certify temporal convergence or whole-vessel behavior',
  source_sha256={n:sha(HERE/n) for n in ['particle_p1_sensitivity.py','exact_cylinder_particle.py']},
  mesh_sha256=sha(ROOT/'manufactured_pipe/R1_GMSH_V2/cfd_gpu/SV_MESH/mesh-complete.mesh.vtu'),source_contract_sha256=sha(ROOT/'PARTICLE_CFD_INPUT_CONTRACT.json'))
 put(OUT/'PARTICLE_DIAGNOSTIC_PROTOCOL.json',protocol)
 put(ROOT/'P1_FIELD_PARTICLE_SENSITIVITY.json',dict(status='RUNNING',protocol=protocol,cases=[]))
 cases=[]
 with ProcessPoolExecutor(max_workers=4) as pool:
  jobs=[pool.submit(run_one,(name,a,gap,kind)) for name,a,gap in SCENARIOS for kind in ['EXACT_ANALYTIC_FIELD','SAME_MESH_P1_ANALYTIC_FIELD']]
  for future in as_completed(jobs):
   result=future.result();cases.append(result);print(json.dumps({k:result[k] for k in ['case','field','status','wall_time_s','transit_time_s','initial_RHS_velocity_m_s']}),flush=True)
   put(ROOT/'P1_FIELD_PARTICLE_SENSITIVITY.json',dict(status='RUNNING',protocol=protocol,cases=cases))
 comparison=[]
 for name,_,_ in SCENARIOS:
  a=next(c for c in cases if c['case']==name and c['field']=='EXACT_ANALYTIC_FIELD');b=next(c for c in cases if c['case']==name and c['field']=='SAME_MESH_P1_ANALYTIC_FIELD')
  comparison.append(metrics(a,b))
 labels=[m['classification'] for m in comparison]
 status='P1_FIELD_TRAJECTORY_IMPACT_LARGE' if 'P1_FIELD_TRAJECTORY_IMPACT_LARGE' in labels else 'P1_FIELD_TRAJECTORY_IMPACT_SMALL' if all(v=='P1_FIELD_TRAJECTORY_IMPACT_SMALL' for v in labels) else 'P1_FIELD_TRAJECTORY_IMPACT_MODERATE' if all(v in ['P1_FIELD_TRAJECTORY_IMPACT_SMALL','P1_FIELD_TRAJECTORY_IMPACT_MODERATE'] for v in labels) else 'INCOMPLETE_DIAGNOSTIC'
 put(ROOT/'P1_FIELD_PARTICLE_SENSITIVITY.json',dict(status=status,DIAGNOSTIC_ONLY=True,protocol=protocol,cases=cases,comparison=comparison,wall_time_s=time.monotonic()-start,original_field_gate_unchanged=True,formal_trajectory_certification=False,dt_convergence_claim=False))
 print(status,flush=True)
if __name__=='__main__':main()
