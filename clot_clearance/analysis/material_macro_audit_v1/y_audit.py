"""Read-only chunked audit of native full-force evidence and every transport row."""
import csv,json,hashlib,sys,time
from pathlib import Path
import numpy as np
from pd_clot.lammps_fragmentation.early_flow_coupling import DTYPE,MASS,DT
from pd_clot.lammps_fragmentation.config import ROOT

def clean(x):
 if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [clean(v) for v in x]
 if isinstance(x,np.ndarray):return x.tolist()
 if isinstance(x,np.generic):return x.item()
 return x

def run(out):
 started=time.perf_counter();base=ROOT/'verification/y_clot_point_cloud_and_early_flow_coupling_v1';case=base/'cases/C1_GAUSSIAN_A200';out.mkdir(exist_ok=True,parents=True)
 allmetrics=[];audits={};inputs=[]
 for name in ['saved_frame_force_audit','detached_force_audit']:
  path=case/(name+'.bin');raw=np.memmap(path,dtype=DTYPE,mode='r');bound=np.r_[0,np.flatnonzero(np.diff(raw['step']))+1,len(raw)]
  prevx=np.full((1927,3),np.nan);prevstep=np.full(1927,-99);maxnode={'speed_m_s':0};maxcom={'speed_m_s':0};maxdisp={'m':0};maxrel={'m_s':0}
  full_peak=0.;visited=[]
  for start,end in zip(bound[:-1],bound[1:]):
   a=raw[start:end];step=int(a['step'][0]);ids=a['atom_id']-1;v=a['v_force_half_m_s']+.5*DT*(a['force_before_flow_N']+a['drag_N']+a['numerical_N'])/MASS
   speed=np.linalg.norm(v,axis=1);q=int(speed.argmax())
   if speed[q]>maxnode['speed_m_s']:maxnode=dict(speed_m_s=float(speed[q]),step=step,time_s=step*DT,atom_id=int(a['atom_id'][q]))
   okay=prevstep[ids]==step-1
   if okay.any():
    ds=np.linalg.norm(a['x_m'][okay]-prevx[ids[okay]],axis=1);j=int(ds.argmax())
    if ds[j]>maxdisp['m']:maxdisp=dict(m=float(ds[j]),step=step,atom_id=int(a['atom_id'][okay][j]),over_dx=float(ds[j]/50e-6))
   prevx[ids]=a['x_m'];prevstep[ids]=step
   free=a['detach_step']>=0;commax=0.
   if free.any():
    _,inv=np.unique(a['component_id'][free],return_inverse=True);counts=np.bincount(inv);cv=np.column_stack([np.bincount(inv,weights=v[free,d])/counts for d in range(3)]);cs=np.linalg.norm(cv,axis=1);j=int(cs.argmax());commax=float(cs[j])
    if commax>maxcom['speed_m_s']:maxcom=dict(speed_m_s=commax,step=step,time_s=step*DT,component_id=int(np.unique(a['component_id'][free])[j]))
    relative=np.linalg.norm(v[free]-cv[inv],axis=1).max()
    if relative>maxrel['m_s']:maxrel=dict(m_s=float(relative),step=step)
   allmetrics.append([name,step,step*DT,float(speed.max()),commax,float(np.linalg.norm(a['wall_N'],axis=1).max()),float(np.linalg.norm(a['drag_N'],axis=1).max())]);visited.append(step)
  audits[name]=dict(rows=len(raw),saved_step_count=len(visited),first_step=min(visited),last_step=max(visited),max_node_full_step_speed=maxnode,max_free_COM_full_step_speed=maxcom,max_adjacent_accepted_position_change=maxdisp,max_component_internal_relative_speed=maxrel,coverage='All currently detached node records at every native step after first detachment' if name=='detached_force_audit' else 'All clot nodes only every20steps OR topology change; maximum over saved records, not certified whole-native-run maximum',velocity_reconstruction='v_full=v_force_half+0.5*dt*(force_before_flow+drag+numerical)/mass')
  inputs.append(path)
 with (out/'Y_NATIVE_SPEED_TIMELINE.csv').open('w') as f:
  w=csv.writer(f);w.writerow(['source','step','time_s','max_node_speed_m_s','max_free_COM_speed_m_s','max_wall_force_N','max_drag_force_N']);w.writerows(allmetrics)
 wall=[]
 for p in sorted(case.glob('wall_rank*.csv')):
  wall.extend(list(csv.DictReader(p.open())));inputs.append(p)
 native_wall=dict(record_rows=len(wall),unique_steps=len(set(int(r['step']) for r in wall)),minimum_center_clearance_m=min(float(r['min_clearance_m']) for r in wall),maximum_contact_nodes=max(int(r['contact_nodes']) for r in wall),outside_node_steps=sum(int(r['outside_nodes']) for r in wall),normal_contact_displacement='No active native wall-contact steps if max_contact_nodes=0; near-wall normal incremental positions not separately saved for all attached nodes')
 tr=base/'FRAGMENT_TRANSPORT_TRAJECTORIES.csv';dts=[];maxspeed={'m_s':0};maxdis={'m':0};clockerr=0.;maxdterr=0;prior={};tlist=[];vlist=[]
 with tr.open() as f:
  for r in csv.DictReader(f):
   cid=int(r['component_id']);t=float(r['relative_time_s']);dt=float(r['dt_s']);v=np.array([float(r[c]) for c in ['vx_m_s','vy_m_s','vz_m_s']]);x=np.array([float(r[c]) for c in ['x_m','y_m','z_m']]);x0=np.array([float(r[c]) for c in ['start_x_m','start_y_m','start_z_m']]);v0=np.array([float(r[c]) for c in ['start_vx_m_s','start_vy_m_s','start_vz_m_s']]);sp=float(max(np.linalg.norm(v),np.linalg.norm(v0)));ds=float(np.linalg.norm(x-x0))
   if sp>maxspeed['m_s']:maxspeed=dict(m_s=sp,component_id=cid,endpoint_relative_time_s=t,start_relative_time_s=t-dt)
   if ds>maxdis['m']:maxdis=dict(m=ds,component_id=cid,relative_time_s=t)
   clockerr=max(clockerr,abs(float(r['absolute_time_s'])-t-.00233));maxdterr=max(maxdterr,abs(t-prior.get(cid,0)-dt));prior[cid]=t;dts.append(dt);tlist.append(t);vlist.append(np.linalg.norm(v))
 te=json.loads((base/'FRAGMENT_TRANSPORT_EVENTS.json').read_text());events=te['events'];comps=te['components']
 transport=dict(rows=len(dts),components=len(prior),dt_quantiles_s=np.quantile(dts,[0,.5,.9,.99,1]),clock_offset_s=.00233,max_clock_error_s=clockerr,max_dt_reconstruction_error_s=maxdterr,max_accepted_endpoint_COM_speed=maxspeed,max_accepted_step_displacement=maxdis,maximum_recorded_wall_projection_m=max((e.get('constraint_correction_m',0) for e in events),default=0),wall_constraint_events=sum(e['event']=='WALL_CONSTRAINT' for e in events),minimum_accepted_center_clearance_m=min(c['minimum_accepted_center_clearance_m'] for c in comps),coverage='All component integrator rows, native endpoint+0 to5s; maxima accepted endpoints, no independent continuous-within-step extrema; rigid offsets, no internal modes')
 np.savez_compressed(out/'Y_TRANSPORT_SPEED_TIMELINE.npz',relative_time_s=tlist,speed_m_s=vlist,dt_s=dts)
 bound_code={}
 for name in ['early_flow_coupling.py','early_flow_transport.py','native_y_wall_contact.py','y_clot_morphology_v3.py']:
  p=ROOT/'pd_clot/lammps_fragmentation'/name;lines=p.read_text().splitlines();bound_code[name]=[{'line':i+1,'text':s} for i,s in enumerate(lines) if any(v in s for v in ['.5e-6','0.5e-6','4660','.00233'])];inputs.append(p)
 cL=math.sqrt((100000/3+4*940/3)/1000) if False else np.sqrt((100000/3+4*940/3)/1000)
 st=dict(cL_m_s=cL,cL_dt_over_dx=cL*DT/50e-6,dx_over_cL_s=50e-6/cL,wall_stiffness_N_m=50,wall_dt_omega=DT*np.sqrt(50/MASS),pair_single_contact_stiffness_N_m=15*(18*(100000/3)/(np.pi*(150.5e-6)**4))*(50e-6)**6/(150.5e-6))
 st['pair_two_equal_mass_dt_omega']=DT*np.sqrt(2*st['pair_single_contact_stiffness_N_m']/MASS);st['scope']='Scale checks only; neighbor connectivity, dilatation coupling and contact network modes not certified by single-spring bounds'
 hand=json.loads((base/'HANDOFF_CONTINUITY_AUDIT.json').read_text())
 result=dict(native=audits,native_wall=native_wall,transport=transport,handoff_existing_audit=hand,stability=st,hardcoded_clocks=bound_code,force_terms='Attached: -500*m*v numerical damping; detached: shared COM fluid drag; wall: unilateral penalty+damping. Different meanings, not combined dissipation.',no_production_changes=True,wall_seconds=time.perf_counter()-started)
 (out/'Y_AUDIT.json').write_text(json.dumps(clean(result),indent=2))
 inputs += [tr,base/'FRAGMENT_TRANSPORT_EVENTS.json',case/'accepted_trajectory.npz']
 hashes={}
 for p in inputs:
  h=hashlib.sha256()
  with p.open('rb') as f:
   for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
  hashes[str(p)]=dict(bytes=p.stat().st_size,sha256=h.hexdigest())
 (out/'Y_INPUT_SHA256.json').write_text(json.dumps(hashes,indent=2));print(json.dumps(clean(result),indent=2)[:9000])
if __name__=='__main__':run(Path(sys.argv[1]))
