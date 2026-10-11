"""Supplementary event interpretation and frozen-stop replay, without a new solve."""
from pathlib import Path
import sys,numpy as np
from scipy.ndimage import label
from metrics import Topology
from design import save,js

def main(base):
 q=js(base/'CRACK_PROPAGATION_AUDIT.json');detail={}
 for name in ['M5B','M5C']:
  p=base/'analysis'/name;g=dict(np.load(p/'geometry_events.npz'));T=Topology(g,g['prescribed_mask']);ev=g['events'];bt=g['break_step'];tip=ev[T.reg[ev[:,1]]==0];step=int(tip[0,0]);ids=tip[tip[:,0]==step,1];out=dict(first_notch_region_event=dict(step=step,time_s=step*.5e-6,bond_indices=ids,reference_midpoints_um=T.mid[ids]*1e6),states={})
  for key,s in q[name]['states'].items():
   alive=bt>s['step'];dead=~alive;grid=np.zeros((13,13,21),bool);kept=dead&T.band;grid[tuple(T.voxel[kept].T)]=True;lab,n=label(grid);near=kept&(abs(T.mid[:,0]+100e-6)<=50e-6+1e-12);labs=np.unique(lab[tuple(T.voxel[near].T)]);labs=labs[labs>0];root=kept&np.isin(lab[tuple(T.voxel.T)],labs);sizes=np.bincount(lab.ravel());central=(abs(T.mid[:,1])<=100e-6+1e-12);survive=T.cross&alive;out['states'][key]=dict(central_target_band_broken=int(np.sum(kept&central)),tip_connected_targetband_broken=int(root.sum()),tip_connected_targetband_central_broken=int(np.sum(root&central)),targetband_tip_connected_cluster_sizes=sizes[labs],surviving_cross_z0_bonds=int(survive.sum()),lost_cross_z0_bonds=int(np.sum(T.cross&dead)),initial_cross_z0_bonds=int(T.cross.sum()),scope='Descriptive six-neighbor broken-midpoint cloud within the previously frozen target band; connection of damage voxels is not a continuous physical free surface or proof of fracture-area advance.')
  detail[name]=out
 cps=js(base/'M5C/CHECKPOINT_TOPOLOGY.json');fail=[]
 for x in cps:
  trigger=None
  if x['window_max_step_dx_over_dx']>.2:trigger='NUMERICAL_STEP_DISPLACEMENT_LIMIT'
  elif x['broken']>=1000 and x['grip_fraction']>=.5:trigger='SEVERE_GRIP_DAMAGE_CONTAMINATION'
  elif x['broken']>=1000 and x['target_band_fraction']<.5 and x['off_plane_fraction']>.6:trigger='SINGLE_TARGET_BAND_INTERPRETATION_LOST'
  if trigger:fail.append([x['step'],trigger])
 assert fail==[[14100,'SINGLE_TARGET_BAND_INTERPRETATION_LOST']]
 detail['stop_replay']=dict(all_checkpoints=len(cps),only_trigger=fail,previous_checkpoint=cps[-2],last_checkpoint=cps[-1],future_data_used=False,all_earlier_checkpoints_below_trigger=True)
 save(base/'EVENT_AND_STOP_DETAIL.json',detail)
 E=js(base/'ENERGY_LEDGER_AUDIT.json');contact=js(base/'CONTACT_FORCE_AUDIT.json');E['contact_saved_state_audit']='CONTACT_FORCE_AUDIT.json';E['contact_force_caution']='The online non-break decomposition has no nonzero contact samples. Its zero maximum is a restricted observation, not an all-time zero force. A separately source-equivalent pair calculation agrees exactly with native contact energy on all saved states; saved-state maximum follows, while the all-step force maximum is unavailable.';E['saved_state_max_contact_node_force_N']=contact['max_force_over_saved_states_N'];E['all_step_contact_force_maximum_N']=None;save(base/'ENERGY_LEDGER_AUDIT.json',E)
 s=js(base/'FINAL_STATUS.json');s.update(overall='IMPLEMENTATION_VALID_KINEMATIC_DEFECT_REMOVED_BUT_NOTCH_AND_CRACK_OBJECTIVE_NOT_MET',first_break=q['M5C']['first_events'],next_priority='基于已有M5/M5B/M5C及Liu原文实验，开展一次只读的断裂观测量可辨识性与本构需求审查，明确哪些实验量可约束LPS及软化需求；结束同类夹具试错，不启动新参数扫描。',next_task_started=False,AI_labels=dict(allowed=['Exact numerical boundary protocol and resource identity','Native first break location/time and connected components','Native damage events and reference-space region counts','Numerical force/opening curves with specimen and discretization identifiers','Signed topology, free kinetic and numerical damping ledger with limitations'],not_certified=['Material intrinsic Gc or fracture-only energy','Mesh-independent fragment-size/strength labels','Unique physical crack area','Quasistatic full fracture trajectory','Biological/patient-specific truth','Microbubble force/AI-control ground truth']));save(base/'FINAL_STATUS.json',s)
 save(base/'POSTPROCESS_NOTES.json',dict(formal_solver_reruns=0,notes=['Initial read-only analysis corrected NumPy iinfo.max property access; no raw changes.','Removed an erroneous postprocessor equality assertion against native bond PE; V1/V2 already document its deviatoric-energy undercount. Retained original discrepancies; reused V2 independent-force-gradient and finite-difference checks unchanged.','No accepted state, simulation input, stopping rule or launch source changed after formal launch.','Signed topology energy is retained and not used as a certified fracture-only energy.']))
 print('detail and frozen-stop replay passed')

if __name__=='__main__':main(Path(sys.argv[1]).resolve())
