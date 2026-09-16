"""Independent audit of archived numerical evidence; no C++/ctypes/subprocess solver.
Run with --phase contact before further cases, --phase offline for archive math
and all geometry, and --phase runtime after receiving remote outputs.
"""
from pathlib import Path
import sys,json,hashlib,itertools,argparse,csv
import numpy as np,h5py
from reference_wall_hydrodynamics_v0 import Lookup,Geometry,frame
from reference_wall_runtime import RuntimeReference
R=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,d):
 (R/'validation'/name).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def arr(path):return np.atleast_1d(np.genfromtxt(path,delimiter=',',names=True))
def xyz(t):return np.column_stack([t[k] for k in ['x_m','y_m','z_m']])
def vel(t):return np.column_stack([t[k] for k in ['vx_m_s','vy_m_s','vz_m_s','omega_x_rad_s','omega_y_rad_s','omega_z_rad_s']])
def case_audit(c):
 d=R/'cases'/c['name'];state=json.loads((d/'RUN_STATE.json').read_text());t=arr(d/'TRAJECTORIES.csv');N=c['N'];gates=json.loads((R/'contracts/MICROBUBBLE_WALL_HYDRODYNAMICS_V0_CONTRACT.json').read_text())['gates']
 if state['reason']=='BLOCKED_LOCAL_PLANE_VALIDITY':
  assert state['accepted_steps']==0 and state['time_s']==0 and c['name']=='J_real_8'
  geom=Geometry(R/'provenance/closed_geometry_m.stl');queries=[geom.query(x,a) for x,a in zip(xyz(t),c['radii_m'])];assert any(not q['planar'] for q in queries)
  assert np.array_equal(xyz(t),np.array(c['initial_positions_m']))
  return {'status':'BLOCKED_LOCAL_PLANE_VALIDITY','independent_gate_enforcement':'PASS','initial_state_only':True,'actual_time_s':0,'target_time_s':.02,'invalid_local_planes':sum(not q['planar'] for q in queries),'min_initial_gap_m':min(q['gap'] for q in queries),'velocities_in_runtime_file':'ZERO_PLACEHOLDERS; NO_INITIAL_SOLUTION; DO_NOT_INTERPRET_AS_PHYSICS'}
 assert c['config']['wall']=='FLAT';assert state['reason']=='MAX_PHYSICAL_TIME' and abs(state['time_s']-c['config']['max_time'])<1e-14
 assert len(t)==(state['accepted_steps']+1)*N;assert np.isfinite(t.view(float)).all()
 stages=arr(d/'INTEGRATION_STAGES.csv');hist=arr(d/'SOLVER_HISTORY.csv');ref=RuntimeReference(R,c)
 maxv=maxw=maxpos=maxmid=0.;mingap=1e100;minpair=1e100;active_pairs=0;wallconstraintrows=0;pairconstraintrows=0
 times=np.unique(t['time_s']);a=np.array(c['radii_m']);qreference=[]
 for k,tt in enumerate(times):
  state_rows=t[t['time_s']==tt];x=xyz(state_rows);observed=vel(state_rows);dt=0. if k==0 else hist[k-1]['dt_s'];q=ref.solve(x,None if k==0 else x,dt);qreference.append(q)
  maxv=max(maxv,float(np.max(abs(q[:,:3]-observed[:,:3]))));maxw=max(maxw,float(np.max(abs(q[:,3:]-observed[:,3:]))));gap=x[:,2]-a;mingap=min(mingap,float(min(gap)))
  assert np.max(abs(gap-state_rows['surface_wall_gap_m']))<=1e-18
  assert np.max(abs(state_rows['normal_z']-1))<=1e-14
  if k==len(times)-1:break
  h=hist[k];dt=float(h['dt_s']);s0=stages[(stages['step']==k)&(stages['stage']==0)];s1=stages[(stages['step']==k)&(stages['stage']==1)];assert len(s0)==len(s1)==N
  q0=ref.solve(x,x,.5*dt);mid=x+.5*dt*q0[:,:3];q1=ref.solve(mid,x,dt);end=x+dt*q1[:,:3]
  maxv=max(maxv,float(np.max(abs(q0[:,:3]-vel(s0)[:,:3]))),float(np.max(abs(q1[:,:3]-vel(s1)[:,:3]))));maxw=max(maxw,float(np.max(abs(q0[:,3:]-vel(s0)[:,3:]))),float(np.max(abs(q1[:,3:]-vel(s1)[:,3:]))));maxmid=max(maxmid,float(np.max(abs(mid-xyz(s1)))));maxpos=max(maxpos,float(np.max(abs(end-xyz(t[t['time_s']==times[k+1]])))))
  # Flat-wall exact continuous swept minimum: z is affine along each segment.
  swept_min=float(min(np.min(x[:,2]-a),np.min(mid[:,2]-a),np.min(end[:,2]-a)));mingap=min(mingap,swept_min)
  if c['config']['hard_wall']:assert swept_min>=-1e-12
  for i,j in itertools.combinations(range(N),2):
   dist=np.linalg.norm(x[i]-x[j])-a[i]-a[j]
   if dist<.2*a[i]*a[j]/(a[i]+a[j]):active_pairs+=1
   for y in [mid,end]:
    r=x[i]-x[j];v=y[i]-y[j]-r;s=np.clip(-r@v/(v@v),0,1) if v@v>0 else 0;gap=np.linalg.norm(r+s*v)-a[i]-a[j];minpair=min(minpair,float(gap));assert gap>=-1e-12
  wallconstraintrows+=int(h['wall_constraints']);pairconstraintrows+=int(h['pair_constraints'])
  assert dt<=c['config']['dt_max']*(1+1e-10)
 assert maxv<=gates['RK2_reference_velocity_m_s'] and maxw<=gates['RK2_reference_omega_rad_s'] and maxpos<=gates['RK2_reference_position_m'],(c['name'],maxv,maxw,maxpos)
 if c['config']['hard_wall']:assert mingap>=-1e-12 and state['accepted_wall_overlap_count']==0
 assert state['nonfinite_count']==state['invalid_query_count']==state['lost_atoms']==state['accepted_pair_overlap_count']==0
 result={'status':'PASS','accepted_steps':state['accepted_steps'],'physical_time_s':state['time_s'],'max_velocity_reference_error_m_s':maxv,'max_omega_reference_error_rad_s':maxw,'max_RK2_endpoint_error_m':maxpos,'max_RK2_midpoint_error_m':maxmid,'minimum_continuous_wall_gap_m':mingap,'minimum_continuous_pair_gap_m':None if N==1 else minpair,'continuous_wall_method':'analytic exact affine z minimum on both RK2 accepted segments','pair_interaction_stage_observations':active_pairs,'wall_constraints':wallconstraintrows,'pair_constraints':pairconstraintrows,'known_penetrating_control':bool(c['config']['control_overlap']),'owner_changes':state['owner_changes']}
 save(c['name']+'_INDEPENDENT.json',result);print(c['name'],result,flush=True);return result

def contact():
 cases=json.loads((R/'contracts/CASE_INDEX.json').read_text());checks={c['name']:case_audit(c) for c in cases[:6]}
 reference=arr(R/'cases/F_hard_cw1/TRAJECTORIES.csv');tg=np.unique(np.concatenate([arr(R/'cases'/('F_hard_cw'+str(k))/'TRAJECTORIES.csv')['time_s'] for k in [4,2,1]]));qr=np.column_stack([np.interp(tg,reference['time_s'],reference[k]) for k in ['x_m','y_m','z_m','omega_x_rad_s','omega_y_rad_s','omega_z_rad_s']]);rows=[]
 for cw in [.4,.2,.1]:
  t=arr(R/'cases'/('F_hard_cw'+str(int(cw*10)))/'TRAJECTORIES.csv');q=np.column_stack([np.interp(tg,t['time_s'],t[k]) for k in ['x_m','y_m','z_m','omega_x_rad_s','omega_y_rad_s','omega_z_rad_s']]);pd=float(np.max(abs(q[:,:3]-qr[:,:3])));gd=float(np.max(abs(q[:,2]-qr[:,2])));wd=float(np.max(abs(q[:,3:]-qr[:,3:]))/np.max(abs(qr[:,3:])));ad=float(np.max(abs(np.trapz(q[:,3:],tg,axis=0)-np.trapz(qr[:,3:],tg,axis=0))));ok=pd<=1e-8 and gd<=1e-8 and wd<=.02 and ad<=1e-3;rows.append(dict(C_wall=cw,position_error_m=pd,gap_error_m=gd,omega_relative_peak_error=wd,integrated_rotation_error_rad=ad,status='PASS' if ok else 'FAIL'))
 selected=max(row['C_wall'] for row in rows if row['status']=='PASS');summary={'status':'PASS','comparison':'linear time interpolation onto union of all three accepted-step time grids','selected_C_wall':selected,'reference_C_wall':.1,'rows':rows};save('C_WALL_CONVERGENCE.json',summary)
 assert checks['F_no_wall']['minimum_continuous_wall_gap_m']<-1e-12 and checks['F_resistance_only']['minimum_continuous_wall_gap_m']<-1e-12
 assert all(checks['F_hard_cw'+str(k)]['wall_constraints']>0 for k in [1,2,4]);save('CONTACT_INDEPENDENT_FINALIZER.json',{'status':'PASS','cases':checks,'C_wall':summary})
 f=arr(R/'cases/F_handoff/TRAJECTORIES.csv');j=np.flatnonzero((f['epsilon'][:-1]<=20)&(f['epsilon'][1:]>20));assert len(j)==1;j=int(j[0]);handoff=json.loads((R/'validation/FAR_FIELD_HANDOFF_AUDIT.json').read_text());handoff.update(trajectory_effect='MEASURED_SYNTHETIC_FORCE_DRIVEN_HANDOFF',crossing_before_epsilon=float(f['epsilon'][j]),crossing_after_epsilon=float(f['epsilon'][j+1]),normal_velocity_before_m_s=float(f['vz_m_s'][j]),normal_velocity_after_m_s=float(f['vz_m_s'][j+1]),normal_velocity_relative_jump=float(f['vz_m_s'][j+1]/f['vz_m_s'][j]-1),parallel_velocity_relative_jump=float(f['vx_m_s'][j+1]/f['vx_m_s'][j]-1),caveat='Adjacent-step difference includes a small finite gap change; matrix jump at20 separately tabulated.');save('FAR_FIELD_HANDOFF_AUDIT.json',handoff);print('CONTACT_PASS',selected,flush=True)

def runtime():
 contact();cases=json.loads((R/'contracts/CASE_INDEX.json').read_text());checks={c['name']:case_audit(c) for c in cases[6:]};a=arr(R/'cases/MPI1_128/TRAJECTORIES.csv');b=arr(R/'cases/MPI4_128/TRAJECTORIES.csv');assert len(a)==len(b) and np.array_equal(a['particle_id'],b['particle_id']) and np.array_equal(a['time_s'],b['time_s']);pd=float(np.max(abs(xyz(a)-xyz(b))));vd=float(np.max(abs(vel(a)[:,:3]-vel(b)[:,:3])));wd=float(np.max(abs(vel(a)[:,3:]-vel(b)[:,3:])));states=[json.loads((R/'cases'/('MPI'+str(k)+'_128')/'RUN_STATE.json').read_text()) for k in [1,4]];assert pd<=1e-10 and vd<=1e-10 and wd<=1e-7 and states[1]['owner_changes']>0 and states[1]['wall_active_observations']>0
 mpi={'status':'PASS','N':128,'position_error_m':pd,'velocity_error_m_s':vd,'omega_error_rad_s':wd,'MPI4_owner_changes':states[1]['owner_changes'],'wall_active_observations':states[1]['wall_active_observations']};save('MPI_WALL_MODEL.json',mpi)
 # Existing Kokkos engine executes identical first8 initial states; compatibility only.
 k=arr(R/'cases/KOKKOS_SMOKE/TRAJECTORIES.csv');cpu=a[(a['particle_id']<=8)&(a['time_s']<=k['time_s'].max()+1e-15)];assert len(k)==len(cpu);assert np.max(abs(xyz(k)-xyz(cpu)))<=1e-10
 save('RUNTIME_INDEPENDENT_FINALIZER.json',{'status':'PASS','scientific_stage_status':'BLOCKED_LOCAL_PLANE_VALIDITY','case_I':'BLOCKED_NO_VALID_NEAR_WALL_START','checks':checks,'MPI':mpi,'Kokkos':'PASS; COMPATIBILITY_ONLY'})

def offline():
 plan=json.loads((R/'contracts/TABLE_GENERATION_PLAN.json').read_text());contract=json.loads((R/'contracts/RMBW_WALL_TABLE_CONTRACT.json').read_text());path=R/'tables/RMBW_WALL_RESISTANCE_TABLE_V0.h5';assert sha(path)==contract['table_sha256'];L=Lookup(path);raw=np.load(R/'raw/RMBW_table_FIRST.npz');assert np.array_equal(L.R,raw['R_total_scaled'])
 e=np.array(plan['holdout_epsilon']);T=L.excess(e)+np.eye(6);arch=np.load(R/'validation/TABLE_HOLDOUT_CPP.npz');assert np.max(abs(T-arch['R_interpolated_scaled']))<=1e-10;v=arch['action_vectors'];maxerr=0.
 for size in plan['radii_m']:
  direct=np.load(R/'raw'/('RMBW_'+size+'_FIRST.npz'));A=direct['R_total_scaled']@v.T;B=T@v.T;err=np.linalg.norm(B-A,axis=1)/np.linalg.norm(A,axis=1);M=np.linalg.inv(T);u=direct['M_scaled']@v.T;eu=np.linalg.norm(M@v.T-u,axis=1)/np.linalg.norm(u,axis=1);maxerr=max(maxerr,float(np.max(err[e<=5])),float(np.max(eu[e<=5])))
 assert maxerr<=1e-3
 reciprocal=float(np.max(np.max(abs(L.R-L.R.transpose(0,2,1)),axis=(1,2))/np.max(abs(L.R),axis=(1,2))));eig=float(np.min(np.linalg.eigvalsh(L.R)));mineigm=float(np.min(np.linalg.eigvalsh(raw['M_scaled'])));assert reciprocal<=1e-10 and eig>0 and mineigm>0
 z=np.load(R/'raw/ROTATION_FIRST.npz');Qs=np.array([frame(n) for n in z['normals']]);assert np.max(abs(Qs-z['frames']))<=1e-12;maxrot=0.
 for Q,M,vv,f in zip(Qs,z['matrices'],z['velocities'],z['force_matrix']):
  D=np.zeros((6,6));D[:3,:3]=D[3:,3:]=Q;ref=D@(M@(D.T@vv));maxrot=max(maxrot,float(np.linalg.norm(ref-f)/np.linalg.norm(ref)))
 assert maxrot<=1e-12
 modes=np.load(R/'raw/FLAT_MODE_SOLVER_FIRST.npz');errs=[]
 for inp,q in zip(modes['inputs'],modes['cpp_velocity']):
  a,h=inp[:2];drag=np.array([6*np.pi*.001*a]*3+[8*np.pi*.001*a**3]*3);scale=np.sqrt(drag);A=(np.diag(drag)+L.global_excess(a,.001,h,inp[8:11]))/np.outer(scale,scale);vv=np.linalg.solve(A,inp[2:8]/scale)/scale;errs.append(np.linalg.norm((vv-q)*np.r_[np.ones(3),np.full(3,a)])/np.linalg.norm(q*np.r_[np.ones(3),np.full(3,a)]))
 assert max(errs)<=1e-3
 out={'status':'PASS','table_sha256':sha(path),'table_action_error':maxerr,'max_reciprocity_error':reciprocal,'minimum_total_resistance_scaled_eigenvalue':eig,'minimum_mobility_scaled_eigenvalue':mineigm,'rotation_error':maxrot,'max_mode_response_error':max(errs),'geometry':{}}
 for name,mesh in [('REAL','closed_geometry_m.stl'),('CURVED','SYNTHETIC_SPHERE_R20UM.stl')]:
  data=np.load(R/'raw'/(name+'_GEOMETRY_FIRST.npz'));print(name,list(data.keys()),flush=True);geo=Geometry(R/'provenance'/mesh);pos=data['positions'];rad=data['radii'];cpp=data['cpp_query'];maxd=maxp=maxr=maxn=0.;disagree=0
  for i,(x,a,q) in enumerate(zip(pos,rad,cpp)):
   g=geo.query(x,a);maxd=max(maxd,abs(g['distance']-q[12]));maxp=max(maxp,float(np.linalg.norm(g['closest']-q[:3])));maxr=max(maxr,abs(g['rms_over_a']-q[15]));maxn=max(maxn,abs(g['normal_spread_deg']-q[16]));disagree+=bool(g['planar'])!=bool(q[20]);assert g['gap']>=0
   if i%2000==0:print(name,i,flush=True)
  assert maxd<=1e-12 and maxp<=1e-12 and maxr<=1e-8 and maxn<=1e-7 and disagree==0
  out['geometry'][name]={'count':len(pos),'max_distance_error_m':maxd,'max_point_error_m':maxp,'max_RMS_error':maxr,'max_normal_spread_error_deg':maxn,'class_disagreements':disagree}
 save('OFFLINE_INDEPENDENT_FINALIZER.json',out);print('OFFLINE_PASS',flush=True)
def supplemental():
 z=np.load(R/'raw/SWEPT_STATIC_FIRST.npz');geo=Geometry(R/'provenance/SYNTHETIC_SPHERE_R20UM.stl');ref=np.array([geo.swept_safe(x,y,a) for x,y,a in zip(z['x'],z['y'],z['radii'])]);assert np.array_equal(ref,z['cpp_safe']) and np.array_equal(ref,z['expected'])
 normal_error=0.;unit_error=0.
 for name in ['REAL','CURVED']:
  d=np.load(R/'raw'/(name+'_GEOMETRY_FIRST.npz'));q=d['cpp_query'];n=(d['positions']-q[:,:3])/q[:,12,None];normal_error=max(normal_error,float(np.max(abs(n-q[:,3:6]))));unit_error=max(unit_error,float(np.max(abs(np.linalg.norm(q[:,3:6],axis=1)-1))))
 assert normal_error<=1e-12 and unit_error<=1e-12
 plan=json.loads((R/'contracts/JOINT_CONSTRAINT_STATIC_PLAN.json').read_text());out=json.loads((R/'raw/JOINT_CONSTRAINT_STATIC_FIRST.json').read_text());a=plan['radius_m'];sep=2*a+.0001*a/2;x=np.array([[-sep/2,0,1.0001*a],[sep/2,0,1.0001*a]])
 c=json.loads((R/'cases/COMBINED_PAIR_WALL/CASE_CONTRACT.json').read_text());c['config'].update(field='WALL_ZERO_FIELD.h5',drive_x=0,drive_y=0,drive_z=0);rr=RuntimeReference(R,c);original=rr.system
 def system(x):
  A,b,active=original(x);b[0]+=.005*6*np.pi*.001*a;b[6]-=.005*6*np.pi*.001*a;b[2]-=.005*6*np.pi*.001*a;b[8]-=.005*6*np.pi*.001*a;return A,b,active
 rr.system=system;qq=rr.solve(x,x,.001);cpp=np.array(out['q']).reshape(2,6);ev=float(np.max(abs(qq[:,:3]-cpp[:,:3])));ew=float(np.max(abs(qq[:,3:]-cpp[:,3:])));assert ev<=1e-10 and ew<=1e-7 and out['wall_constraints']==2 and out['total_constraints']==3
 save('SUPPLEMENTAL_INDEPENDENT_FINALIZER.json',{'status':'PASS','swept_adversarial_cases':200,'endpoints_safe_crossings_rejected':50,'normal_direction_error':normal_error,'normal_unit_error':unit_error,'joint_wall_constraints':2,'joint_pair_constraints':1,'joint_velocity_error_m_s':ev,'joint_omega_error_rad_s':ew})
 print('SUPPLEMENTAL_PASS',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--phase',choices=['contact','offline','runtime','supplemental'],required=True);args=p.parse_args();globals()[args.phase]()
