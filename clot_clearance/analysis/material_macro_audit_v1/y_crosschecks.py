"""Additional read-only time-level and native normal-displacement checks."""
import sys,json,csv
from pathlib import Path
import numpy as np
from pd_clot.lammps_fragmentation.early_flow_coupling import DTYPE,DT,MASS
from pd_clot.lammps_fragmentation.config import ROOT
from pd_clot.lammps_fragmentation.native_y_wall_query import NativeYWallQuery
from pd_clot.lammps_fragmentation.y_clot_morphology_v3 import load_geometry
out=Path(sys.argv[1]);base=ROOT/'verification/y_clot_point_cloud_and_early_flow_coupling_v1';case=base/'cases/C1_GAUSSIAN_A200'
tr=np.load(case/'accepted_trajectory.npz');lookup={int(s):i for i,s in enumerate(tr['step'])};raw=np.memmap(case/'saved_frame_force_audit.bin',mode='r',dtype=DTYPE);maxve=maxxe=0.;frames=0
for chunk in raw.reshape(-1,1927):
 step=int(chunk['step'][0])
 if step not in lookup:continue
 j=lookup[step];ids=chunk['atom_id']-1;v=chunk['v_force_half_m_s']+.5*DT*(chunk['force_before_flow_N']+chunk['drag_N']+chunk['numerical_N'])/MASS
 maxve=max(maxve,float(np.max(abs(v-tr['velocity_m_s'][j,ids]))));maxxe=max(maxxe,float(np.max(abs(chunk['x_m']-tr['position_m'][j,ids]))));frames+=1
assert frames==len(lookup) and maxve<1e-10 and maxxe<1e-14
geo=load_geometry();b=geo['bonds'];e=np.load(case/'ONLINE_TOPOLOGY_EVENTS.npz')['step_and_initial_bond_index'];alive=np.ones(len(b),bool);maxrelative=0;stepmax=0
for j,step in enumerate(tr['step']):
 alive[e[e[:,0]<=step,1]]=False
 v=tr['velocity_m_s'][j];ab=b[alive];rel=np.linalg.norm(v[ab[:,0]]-v[ab[:,1]],axis=1)
 if rel.max(initial=0)>maxrelative:maxrelative=float(rel.max());stepmax=int(step)
# Match the actual contact rows; query native wall normals only at those saved positions.
te=json.loads((base/'FRAGMENT_TRANSPORT_EVENTS.json').read_text());contacts={(int(e['component_id']),float(e['time_s'])):e for e in te['events'] if e['event']=='WALL_CONSTRAINT'}
hand=json.loads((base/'FRAGMENT_HANDOFF.json').read_text());end=np.load(case/'terminal_full_precision.npz')['position_m'];members={int(c['component_id']):end[np.array(c['member_atom_ids'])-1]-np.array(c['COM_xyz_m']) for c in hand['components']}
q=NativeYWallQuery();maxnormal=0.;maxn=None;count=0
for r in csv.DictReader((base/'FRAGMENT_TRANSPORT_TRAJECTORIES.csv').open()):
 cid=int(r['component_id']);t=float(r['relative_time_s'])
 if (cid,t) not in contacts:continue
 old=np.array([float(r[k]) for k in ['start_x_m','start_y_m','start_z_m']]);new=np.array([float(r[k]) for k in ['x_m','y_m','z_m']]);sample=q.query(old+members[cid]);values=(new-old)@sample['inward_normal'].T;val=float(np.max(abs(values)));count+=1
 if val>maxnormal:maxnormal=val;maxn=dict(component_id=cid,time_s=t,dt_s=float(r['dt_s']),abs_start_wall_normal_displacement_m=val)
result=dict(native_saved_full_step_velocity_reconstruction=dict(frames=frames,max_velocity_error_m_s=maxve,max_coordinate_error_m=maxxe,status='PASS'),native_surviving_bond_max_relative_velocity=dict(m_s=maxrelative,step=stepmax,time_s=stepmax*DT,scope='234 native saved full-step trajectory frames; not every native step'),transport_contact_normal_displacement=dict(contact_rows_matched=count,expected=len(contacts),maximum=maxn,scope='Actual accepted endpoint translation projected onto native wall inward normals at starting member positions; signed normal can change duringstep; projection correction is separate'))
assert count==len(contacts)
(out/'Y_CROSSCHECKS.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
