"""Actual CFD and analytic-node fields through repaired production WSS; no fits."""
from pathlib import Path
import sys,json,csv,argparse,subprocess
import numpy as np
import pyvista as pv
from numpy.polynomial.legendre import leggauss
from case_common import *
from flow_solver_support.wss_case import recover,material
from flow_solver_support.flow_parser import parse_solver_log

def csvout(path,rows):
 with Path(path).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def quant(v,w):
 i=np.argsort(v);return np.interp([.05,.5,.95],(np.cumsum(w[i])-.5*w[i])/w.sum(),v[i])
def volume_error(x,t,u,ua,R,U,L):
 centers=x[t].mean(axis=1);sel=(centers[:,2]>=.3*L)&(centers[:,2]<=.7*L);t=t[sel];xx=x[t];uu=u[t];uv=ua[t]
 det=np.linalg.det(xx[:,1:]-xx[:,:1]);err=ref=inter=0.;quads=[]
 for order in [4,3,3]:
  q,w=leggauss(order);quads.append(((q+1)/2,w/2))
 for s,ws in zip(*quads[0]):
  for r,wr in zip(*quads[1]):
   for z,wz in zip(*quads[2]):
    N=np.array([(1-s)*(1-r)*(1-z),s,(1-s)*r,(1-s)*(1-r)*z]);wt=ws*wr*wz*(1-s)**2*(1-r)
    xxq=np.einsum('i,nij->nj',N,xx);uq=np.einsum('i,nij->nj',N,uu);aq=np.einsum('i,nij->nj',N,uv)
    exact=np.zeros_like(uq);exact[:,2]=2*U*(1-np.sum(xxq[:,:2]**2,axis=1)/R**2)
    err+=wt*np.dot(det,np.sum((uq-exact)**2,axis=1));ref+=wt*np.dot(det,np.sum(exact**2,axis=1));inter+=wt*np.dot(det,np.sum((uq-aq)**2,axis=1))
 return float(100*np.sqrt(err/ref)),float(100*np.sqrt(inter/ref))

def analyze(case):
 case=Path(case);policy=json.loads((case/'policy.json').read_text());execution=json.loads((case/'reports/execution.json').read_text());assert execution['status']=='PASS',case
 out=case/'wss';cmd=[sys.executable,'-B',str(C/'rotate_visualization/prepare_surface_data.py'),'--case',str(case),'--output',str(out)]
 with (case/'reports/production_wss.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
 mesh=np.load(case/'SV_MESH/mesh_arrays.npz');flow=np.load(case/'frozen_flow/flow_arrays_si.npz');x=mesh['points_m'];t=mesh['tetra'];b=mesh['boundary_triangles'];tags=mesh['facet_tags'];u=flow['velocity_m_s'];p=flow['pressure_pa'];analytic=np.load(case/'analytical_nodal.npz');ua=analytic['velocity_m_s'];pa=analytic['pressure_pa'];props=material(case);mu=props['mu_Pa_s'];R=policy['R_m'];L=policy['L_m'];U=policy['Umean_m_s'];theory=4*mu*U/R
 log=parse_solver_log((case/'run/solver.log').read_text(errors='replace'),policy['dt_s']);solves=log['linear_solves'];ratios=[]
 for s in solves:
  ms=s.get('petsc_monitor',[])
  if ms:ratios.append(ms[-1]['true_residual_norm']/max(1e-24,1e-10*ms[0]['true_residual_norm']))
 quality=dict(parsed_solves=len(solves),unparsed_rows=len(log['unparsed_rows']),failed_linear_solves=log['failed_linear_solves'],nonlinear_failure_messages=log['nonlinear_failure_messages'],recovered_attempts=log['recovered_attempts'],last_nonlinear_relative=solves[-1]['nonlinear_Ri_over_R0'],max_true_residual_over_criterion=max(ratios),final_state=execution['states'][-1])
 quality['accepted']=bool(solves and not log['unparsed_rows'] and not log['failed_linear_solves'] and not log['nonlinear_failure_messages'] and max(ratios)<1.01)
 dump(case/'reports/independent_log_checks.json',quality)
 # Snapshot methods remain useful evidence if a gate fails, but no result is labelled valid.
 err,inter=volume_error(x,t,u,ua,R,U,L)
 n=policy['radial_intervals']*policy['circumferential_intervals']+1;nz=policy['axial_intervals'];base=b[tags==4];sections=[]
 for fraction in [0,.3,.5,.7,1]:
  k=int(round(nz*fraction));tri=base+k*n;xyz=x[tri];areas=.5*np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)
  sections.append(dict(case=case.name,z_um=k*L/nz*1e6,area_um2=areas.sum()*1e12,Q_m3_s=float(np.dot(areas,u[tri,2].mean(axis=1))),Q_theory_m3_s=np.pi*R**2*U,pressure_area_mean_Pa=float(np.average(p[tri].mean(axis=1),weights=areas)),weight='section_triangle_area'))
 csvout(case/'reports/sections.csv',sections);dp=sections[1]['pressure_area_mean_Pa']-sections[3]['pressure_area_mean_Pa'];dpth=8*mu*U*(.4*L)/R**2
 rows=[]
 for kind,vel,pres in [('CFD',u,p),('analytic_nodal',ua,pa)]:
  surf,d=recover(x,t,b,tags,vel,pres,mu);cent=d['centers'];sel=(cent[:,2]>=.3*L)&(cent[:,2]<=.7*L);ww=d['magnitude'][sel];area=d['area'][sel];mean=float(np.average(ww,weights=area));qs=quant(ww,area)
  if kind=='analytic_nodal':surf.save(case/'wss/data/analytic_node_wall_wss_si.vtp')
  row=dict(case=case.name,mesh=case.name.replace('_halfdt',''),field_type=kind,region='central_z_0.3L_to_0.7L',statistical_weight='surface_area_for_WSS;volume_Duffy36_for_velocity;section_area_for_pressure',nodes=len(x),tetra=len(t),dt_s=policy['dt_s'],R_um=R*1e6,L_um=L*1e6,radial_intervals=policy['radial_intervals'],circumferential_intervals=policy['circumferential_intervals'],axial_intervals=nz,radial_step_um=R/policy['radial_intervals']*1e6,axial_step_um=L/nz*1e6,theory_WSS_Pa=theory,wss_area_mean_Pa=mean,wss_p05_Pa=float(qs[0]),wss_p50_Pa=float(qs[1]),wss_p95_Pa=float(qs[2]),wss_min_Pa=float(ww.min()),wss_max_Pa=float(ww.max()),wss_mean_signed_error_pct=100*(mean/theory-1),wss_area_L2_error_pct=float(100*np.sqrt(np.average((ww-theory)**2,weights=area))/theory),velocity_volume_L2_error_pct=err if kind=='CFD' else volume_error(x,t,ua,ua,R,U,L)[0],CFD_vs_analytic_nodal_L2_pct=inter if kind=='CFD' else 0.,central_delta_p_Pa=dp if kind=='CFD' else dpth,theory_central_delta_p_Pa=dpth,central_delta_p_signed_error_pct=100*(dp/dpth-1) if kind=='CFD' else 0.,inlet_Q_m3_s=sections[0]['Q_m3_s'] if kind=='CFD' else policy['Q_target_m3_s'],theory_Q_m3_s=np.pi*R**2*U,Q_signed_error_pct=100*(sections[0]['Q_m3_s']/(np.pi*R**2*U)-1),epsilon_mass=execution['states'][-1]['epsilon_mass'] if kind=='CFD' else 0.,solver_steady_and_log_checks=quality['accepted'] if kind=='CFD' else 'not_a_CFD_solve',reference='continuous_Poiseuille_R4um_U2mmps')
  rows.append(row)
 csvout(case/'reports/pipe_validation.csv',rows)
 grid=pv.read(case/'frozen_flow/steady_flow.vtu');xx=np.linspace(-.999*R,.999*R,161);pts=np.column_stack([xx,np.zeros_like(xx),np.full_like(xx,L/2)]);sample=pv.PolyData(pts).sample(grid)
 profile=[dict(case=case.name,x_um=float(a*1e6),y_um=0.,z_um=L*.5*1e6,velocity_z_m_s=float(v[2]),velocity_theory_m_s=float(2*U*(1-(a/R)**2)),valid=int(ok)) for a,v,ok in zip(xx,sample['Velocity'],sample['vtkValidPointMask'])]
 csvout(case/'reports/velocity_profile.csv',profile)
 print(json.dumps(rows[0],indent=2),flush=True)
 return rows

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--case',type=Path);args=parser.parse_args()
 cases=[args.case] if args.case else sorted((V/'stage2').glob('pipe_*'))
 for case in cases:
  if (case/'reports/execution.json').exists():analyze(case)
 allrows=[]
 for p in sorted((V/'stage2').glob('pipe_*/reports/pipe_validation.csv')):allrows.extend(list(csv.DictReader(p.open())))
 if allrows:csvout(V/'data/pipe_validation.csv',allrows)
