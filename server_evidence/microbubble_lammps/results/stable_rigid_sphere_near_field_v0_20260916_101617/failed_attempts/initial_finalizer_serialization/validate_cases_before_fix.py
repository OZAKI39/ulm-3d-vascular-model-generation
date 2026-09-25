from pathlib import Path
import sys,json,itertools
import numpy as np
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'src'))
from reference_rigid_sphere_lubrication import FlowGradientReference,dense_solve,assemble

def table(p):return np.atleast_1d(np.genfromtxt(p,delimiter=',',names=True,dtype=None,encoding='utf-8'))
def xyz(d):return np.column_stack([d[k] for k in ['x_m','y_m','z_m']])
def qvals(d):
 names=['vx_m_s','vy_m_s','vz_m_s']+(['omega_x','omega_y','omega_z'] if 'omega_x' in d.dtype.names else ['omega_x_rad_s','omega_y_rad_s','omega_z_rad_s'])
 return np.column_stack([d[k] for k in names])
def driving(c,ids):
 d=np.zeros((len(ids),6));s=c['drive_scale'];mode=c['drive_mode']
 if mode==1:d[:,5]=s
 if mode==2:d[:,5]=np.where(ids%2,s,-s)
 if mode==3:
  for k in range(3):d[:,k]=s*np.sin(ids*1.7+k);d[:,3+k]=1e6*s*np.cos(ids*2.3+k)
 return d

def evaluate_case(name):
 D=R/'cases'/name;c=json.loads((D/'CASE_CONTRACT.json').read_text());state=json.loads((D/'RUN_STATE.json').read_text());tr=table(D/'TRAJECTORIES.csv');stage=table(D/'INTEGRATION_STAGES.csv');hist=table(D/'SOLVER_HISTORY.csv');ids=sorted(c['ids']);n=len(ids);a=np.array([c['radii_m'][c['ids'].index(i)] for i in ids]);pairs=list(itertools.combinations(range(n),2));flow=FlowGradientReference(R/c['field']);vmax=omax=xerr=0;minswept=np.inf;contacts=0
 U,O,E,_=flow.query(xyz(stage));d=driving(c,np.array(ids))
 for step in np.unique(stage['step']):
  slots0=np.flatnonzero((stage['step']==step)&(stage['stage']==0));slots1=np.flatnonzero((stage['step']==step)&(stage['stage']==1));s0=stage[slots0];s1=stage[slots1];assert list(s0['particle_id'])==ids and list(s1['particle_id'])==ids;X0=xyz(s0);X1=xyz(s1);dt=float(s0['dt_s'][0]);q0=dense_solve(X0,a,U[slots0],O[slots0],E[slots0],pairs,base=X0,dt=dt*.5,driving=d);q1=dense_solve(X1,a,U[slots1],O[slots1],E[slots1],pairs,base=X0,dt=dt,driving=d)
  for qr,actual in [(q0,qvals(s0)),(q1,qvals(s1))]:vmax=max(vmax,float(np.max(np.abs(qr[:,:3]-actual[:,:3]))));omax=max(omax,float(np.max(np.abs(qr[:,3:]-actual[:,3:]))))
  xerr=max(xerr,float(np.max(np.abs(X1-(X0+.5*dt*q0[:,:3])))))
  new=tr[tr['step']==step+1];assert len(new)==n;Xnext=xyz(new);xerr=max(xerr,float(np.max(np.abs(Xnext-(X0+dt*q1[:,:3])))))
  for i,j in pairs:
   r=X0[i]-X0[j];v=(Xnext[i]-X0[i])-(Xnext[j]-X0[j]);t=np.clip(-r@v/(v@v),0,1) if v@v>0 else 0;gap=np.linalg.norm(r+t*v)-a[i]-a[j];minswept=min(minswept,gap);contacts+=gap< -1e-12
 # Every reported accepted velocity/omega independently re-solved.
 u,o,e,g=flow.query(xyz(tr));finalv=finalo=0
 for step in np.unique(tr['step']):
  slots=np.flatnonzero(tr['step']==step);rr=tr[slots];x=xyz(rr);dt=float(hist['dt_s'][np.flatnonzero(hist['step']==step)[0]]) if step else 0
  qr=dense_solve(x,a,u[slots],o[slots],e[slots],pairs,base=x if step else None,dt=dt,driving=d);actual=qvals(rr);finalv=max(finalv,float(np.max(np.abs(qr[:,:3]-actual[:,:3]))));finalo=max(finalo,float(np.max(np.abs(qr[:,3:]-actual[:,3:]))))
 finite=all(np.isfinite(tr[k]).all() for k in tr.dtype.names);dia=np.column_stack([tr[tr['particle_id']==i]['diameter_um'] for i in ids]);diamerr=float(np.max(np.abs(dia-2e6*a)))
 checks={'terminal':state['reason'] in ['MAX_PHYSICAL_TIME','WALL_SAFETY_STOP'],'finite':bool(finite),'stage_velocity_reference':vmax<=1e-9,'stage_omega_reference':omax<=1e-6,'accepted_velocity_reference':finalv<=1e-9,'accepted_omega_reference':finalo<=1e-6,'RK2_stage_positions':xerr<=1e-9,'continuous_swept_nonoverlap':contacts==0,'fixed_diameter':diamerr<=1e-12,'residual':state['max_residual']<=1e-10}
 r={'case':name,'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'accepted_steps':state['accepted_steps'],'time_s':state['time_s'],'reason':state['reason'],'max_stage_velocity_error_m_s':vmax,'max_stage_omega_error_rad_s':omax,'max_accepted_velocity_error_m_s':finalv,'max_accepted_omega_error_rad_s':finalo,'max_RK2_position_error_m':xerr,'independent_min_swept_gap_m':float(minswept) if np.isfinite(minswept) else None,'independent_overlap_count':contacts,'max_diameter_error_um':diamerr,'constraint_steps':state['constraint_steps'],'cross_rank_active_pair_observations':state['cross_rank_active_pair_observations']}
 (D/'INDEPENDENT_VALIDATION.json').write_text(json.dumps(r,indent=2)+'\n');return r

def gap_convergence():
 results={};selected=None
 for cg in [.4,.2,.1]:
  checks=[]
  for letter,label in [('A','normal'),('F','contact')]:
   cur=table(R/'cases'/f'{letter}_{label}_cg{int(cg*10):02d}'/'TRAJECTORIES.csv');ref=table(R/'cases'/f'{letter}_{label}_cg01'/'TRAJECTORIES.csv');maxx=maxg=0
   for tag in np.unique(cur['particle_id']):
    c=cur[cur['particle_id']==tag];r=ref[ref['particle_id']==tag];xx=xyz(c);rr=np.column_stack([np.interp(c['time_s'],r['time_s'],r[k]) for k in ['x_m','y_m','z_m']]);maxx=max(maxx,float(np.max(np.linalg.norm(xx-rr,axis=1))));maxg=max(maxg,float(np.max(np.abs(c['nearest_gap_m']-np.interp(c['time_s'],r['time_s'],r['nearest_gap_m'])))))
   checks.append({'case':letter,'max_position_diff_m':maxx,'max_gap_diff_m':maxg,'pass':maxx<=1e-8 and maxg<=1e-8})
  results[str(cg)]=checks
  if selected is None and all(x['pass'] for x in checks):selected=cg
 r={'status':'PASS' if selected else 'FAIL','selected_C_gap':selected,'reference_C_gap':.1,'gates_from':'RIGID_SPHERE_LUBRICATION_REFERENCE_CONTRACT.json','comparisons':results};(R/'validation/CGAP_CONVERGENCE.json').write_text(json.dumps(r,indent=2)+'\n');return r
if __name__=='__main__':
 if sys.argv[1]=='cgap':print(json.dumps(gap_convergence(),indent=2))
 else:
  r=evaluate_case(sys.argv[1]);print(json.dumps(r,indent=2));raise SystemExit(0 if r['status']=='PASS' else 2)
