"""Read-only analysis of the single failed V3 attempt. No continuation/rerun."""
from pathlib import Path
import json,csv,numpy as np,hashlib
from pd_clot.lammps_fragmentation.fragment_analysis import read_dump
from pd_clot.lammps_fragmentation.detached_fragment_metrics import graph_components
from pd_clot.lammps_fragmentation.fragment_semantics import component_metrics
from pd_clot.lammps_fragmentation.early_flow_coupling import DTYPE
C=json.loads(Path(__file__).with_name('config.json').read_text());R=Path.cwd();O=R/'verification'/C['task']/C['run'];P=O/'V3_replay';OLD=R/'verification/y_clot_point_cloud_and_early_flow_coupling_v1/cases/C1_GAUSSIAN_A200';G=np.load(R/'verification/y_clot_realistic_morphology_wall_contact_and_flow_transport_v1/V3_MATERIAL_POINTS_AND_BONDS.npz');b=G['bonds'];typ=G['bond_types'];n=1927;N=len(G['X']);keys=b[:,0]*N+b[:,1]
def complete(path):
 frames=[];err=None
 try:
  for f in read_dump(path):frames.append(f)
 except Exception as e:err=repr(e)
 return frames,err
new,err=complete(P/'trajectory.lammpstrj');old=np.load(OLD/'accepted_trajectory.npz');native_events,_=complete(P/'broken_bonds.local');blocks=[a for _,_,a in native_events if len(a)];events=np.concatenate(blocks) if blocks else np.empty((0,10));ends=np.sort(events[:,2:4].astype(int)-1,axis=1);idx=np.searchsorted(keys,ends[:,0]*N+ends[:,1]);assert np.array_equal(keys[idx],ends[:,0]*N+ends[:,1]);eventmap={int(ii):int(step) for ii,step in zip(idx,events[:,1])};oev=np.load(OLD/'broken_bonds.npz');old_events=oev['events'];old_idx=oev['initial_bond_indices'];alive=np.ones(len(b),bool);oldalive=alive.copy();rows=[];positions=[];velocities=[];alives=[];old_alives=[]
for step,cols,a in new:
 a=a[np.argsort(a[:,0])];assert len(a)==N and np.array_equal(a[:,0],np.arange(1,N+1));x=a[:,[cols.index(c) for c in ['x','y','z']]];v=a[:,[cols.index(c) for c in ['vx','vy','vz']]];k=int(np.flatnonzero(old['step']==step)[0]);alive[idx[events[:,1]<=step]]=False;oldalive[old_idx[old_events[:,1]<=step]]=False
 lab,count,att=graph_components(n,b,typ,alive);ol,oc,oa=graph_components(n,b,typ,oldalive)
 m=component_metrics(count,att,50e-6);om=component_metrics(oc,oa,50e-6)
 row={'step':int(step),'time_s':step*.5e-6,'max_position_difference_m':float(np.linalg.norm(x-old['position_m'][k],axis=1).max()),'max_velocity_difference_m_s':float(np.linalg.norm(v-old['velocity_m_s'][k],axis=1).max()),'broken_count':int((~alive).sum()),'old_broken_count':int((~oldalive).sum()),'topology_identical':bool(np.array_equal(alive,oldalive)),'vmax_m_s':float(np.linalg.norm(v[:n],axis=1).max()),'old_vmax_m_s':float(np.linalg.norm(old['velocity_m_s'][k,:n],axis=1).max()),'KE_J':float(.5*1.25e-10*np.sum(v[:n]**2)),'old_KE_J':float(.5*1.25e-10*np.sum(old['velocity_m_s'][k,:n]**2)),'detached_nodes':int(count[~att].sum()),'old_detached_nodes':int(oc[~oa].sum()),'detached_components':m['detached_component_count'],'old_detached_components':om['detached_component_count'],'largest_detached_nodes':m['largest_detached_nodes'],'N_detached_ge200':m['N_detached_ge200']}
 rows.append(row);positions.append(x);velocities.append(v);alives.append(alive.copy());old_alives.append(oldalive.copy())
with (O/'V3_PARTIAL_COMPARISON.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
np.savez_compressed(O/'V3_PARTIAL_COMPARISON.npz',step=[a['step'] for a in rows],new_position_m=positions,new_velocity_m_s=velocities,new_alive=alives,old_alive=old_alives,old_position_m=old['position_m'][:len(rows)],old_velocity_m_s=old['velocity_m_s'][:len(rows)],bonds=b,types=typ,X=G['X'])
# Complete binary audit rows only (truncated final write is explicitly counted).
fbytes=(P/'saved_frame_force_audit.bin').stat().st_size;a=np.fromfile(P/'saved_frame_force_audit.bin',dtype=DTYPE,count=fbytes//DTYPE.itemsize);oa=np.memmap(OLD/'saved_frame_force_audit.bin',dtype=DTYPE,mode='r');forc=[]
for step in np.unique(a['step']):
 z=a[a['step']==step];z=z[np.argsort(z['atom_id'])];oz=oa[oa['step']==step];oz=oz[np.argsort(oz['atom_id'])]
 if len(z)!=n or len(oz)!=n:continue
 forc.append({'step':int(step),'max_preflow_force_difference_N':float(np.linalg.norm(z['force_before_flow_N']-oz['force_before_flow_N'],axis=1).max()),'max_position_difference_m':float(np.linalg.norm(z['x_m']-oz['x_m'],axis=1).max())})
with (O/'V3_PARTIAL_FORCE_COMPARISON.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(forc[0]));w.writeheader();w.writerows(forc)
def first_type(ev,ii,t):
 q=ev[typ[ii]==t];return int(q[:,1].min()) if len(q) else None
hist={}
for rank in [0,1]:
 wall=np.genfromtxt(OLD/f'wall_rank{rank}.csv',delimiter=',',names=True);changes=np.flatnonzero(np.diff(wall['mobile_nodes'])!=0)+1;hist[str(rank)]={'initial_mobile_nodes':int(wall['mobile_nodes'][0]),'final_mobile_nodes':int(wall['mobile_nodes'][-1]),'ownership_count_change_steps':wall['step'][changes].astype(int).tolist(),'first_count_change_step':int(wall['step'][changes[0]]) if len(changes) else None}
report={'status':'FAILED_PARTIAL_V3_SIGSEGV','formal_attempts':1,'rerun':False,'successful_full_replay':False,'library_variant':'wvolume-only intermediate; final ownership guard NOT V3 replayed','complete_saved_frames':len(rows),'last_complete_frame_step':rows[-1]['step'],'trajectory_truncation_error':err,'last_complete_wall_row_step':300,'failure_location':'MPI rank1 BondBPM::process_broken+0x267','new_first_bulk_step':first_type(events,idx,1),'old_first_bulk_step':first_type(old_events,old_idx,1),'new_first_interface_step':first_type(events,idx,2),'old_first_interface_step':first_type(old_events,old_idx,2),'new_first_detach_status':'NOT_REACHED_IN_SAVED_PARTIAL_TRAJECTORY','old_first_detach_step':593,'first_saved_position_diff_over_1e-12_m':next((z for z in rows if z['max_position_difference_m']>1e-12),None),'first_saved_force_diff_over_1e-12_N':next((z for z in forc if z['max_preflow_force_difference_N']>1e-12),None),'last_comparable_frame':rows[-1],'historical_MPI_migration_evidence':hist,'historical_sort':'default1000 enabled; initial sort() called by Verlet setup; exact nonidentity historical permutations/times not logged, later-sort execution not directly certified','force_audit_incomplete_bytes':fbytes%DTYPE.itemsize,'claim_limits':'Partial before first detachment only. No final fragment-size improvement inferred. New m-identity observer buffered outputs were lost at native crash; no invented per-step PASS. Final guard independently tested on tiny MPI fracture only.','historical_terminal':{'detached_nodes':217,'detached_components':208,'largest_detached_nodes':4,'ge200':0},'source_input_identity':json.loads((P/'REPLAY_DESIGN_AND_IDENTITY.json').read_text())}
(O/'V3_REPLAY_COMPARISON.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k not in ['historical_MPI_migration_evidence','source_input_identity']},indent=2))
