from common import *
from native import Native
from checks import compare,topology
import shutil,time
def decorate(row,s):
 free=s['component_sizes'][~s['component_attached']]
 for k in [14,34,114]:row['N_ge'+str(k)+'_nodes']=int(np.sum(free>=k))
 return row
class ComponentRecorder:
 def __init__(self):
  self.stream=(O/'full/component_observations.jsonl').open('w');self.max_size=0;self.old_large=set()
 def accept(self,s,prev,j):
  labels=s['component_label'];sizes=s['component_sizes'];att=s['component_attached'];large={tuple(np.flatnonzero(labels==q)+1) for q in np.flatnonzero((~att)&(sizes>=14))};maximum=int(sizes[~att].max(initial=0));changed=prev is None or not np.array_equal(prev['component_label'],labels) or not np.array_equal(prev['detached'],s['detached'])
  if large!=self.old_large or maximum>self.max_size and maximum>=3:j.key(s,'large_or_record_component_topology')
  if changed or s['step']%20==0 or abs(s['time']-END)<1e-15:
   comps=[];b=j.ref['bonds'];bt=j.ref['bond_types'];n=int(j.ref['clot_count'])
   for q in np.flatnonzero(~att):
    ids=np.flatnonzero(labels==q);mask=np.zeros(len(j.ref['X']),bool);mask[ids]=True;inside=mask[b[:,0]]&mask[b[:,1]]&(bt==1);outside=mask[b[:,0]]^mask[b[:,1]];vv=s['velocity'][ids];com=vv.mean(0)
    comps.append(dict(atom_ids=ids+1,nodes=len(ids),reference_bulk=int(inside.sum()),alive_bulk=int(np.sum(inside&(s['d']<1))),mean_internal_d=float(s['d'][inside].mean()) if inside.any() else None,external_alive=int(np.sum(outside&(s['d']<1))),COM_position_m=s['x'][ids].mean(0),COM_velocity_m_s=com,K_internal_J=float(.5*MASS*np.sum((vv-com)**2))))
   j.line(self.stream,dict(step=s['step'],time_s=s['time'],reason='TOPOLOGY_CHANGE' if changed else 'REGULAR_ACCEPTED_SAMPLE',components=comps))
  self.max_size=max(self.max_size,maximum);self.old_large=large
 def close(self):self.stream.flush();self.stream.close()
def early_checkpoint(ns,s,m,integ,j,row,first):
 assert abs(s['time']-EARLY)<1e-15 and not integ.queue
 start=time.perf_counter();j.key(s,'formal_0p8ms');cp=j.checkpoint(ns,s,'formal_0p8ms');folder=O/'early_checkpoint';folder.mkdir();target=folder/'accepted.restart';shutil.copy2(cp,target);save_state(folder/'accepted.state.npz',s)
 save_state(folder/'python_ledger_state.npz',dict(initial_x=m.initial['x'],initial_scalars=m.initial['scalars'],initial_K=m.initial['K'],initial_Uc=m.initial['Uc'],initial_Uw=m.Uw0,previous_contact_force=m.prevfc,last_material_force=m.last['force'],last_material_U=m.last['U'],defects=m.defects))
 dump(folder/'PYTHON_CONTINUATION_STATE.json',dict(work=m.work,Jc=m.Jc,scale=m.scale,maxerr=m.maxerr,flow_supply=m.flow_supply,drag_dissipation=m.drag_dissipation,integrator_queue=integ.queue,serial=integ.serial,accepts=integ.accepts,rejects=integ.rejects,base_intervals=integ.base_intervals,first_events=first,latest_row=row))
 history=np.genfromtxt(B/'V3_NATIVE_PILOT_HISTORY.csv',delimiter=',',names=True);old=history[np.flatnonzero(history['time_s']<=EARLY+1e-15)[-1]]
 comparison={k:dict(normal=float(row[k]),historical_tangent=float(old[k])) for k in ['bulk_failed','interface_failed','detached_nodes','detached_components','largest_detached_nodes','N_ge200','N_ge300','max_speed_m_s','max_contact_penetration_dx','W_load_J','W_drag_J','W_numeric_J','U_J','K_J','D_J','R_J','relative_R']}
 dump(O/'EARLY_WINDOW_COMPARISON.json',dict(status='SHORT_WINDOW_NUMERICAL_PASS',normal_time_s=s['time'],historical_actual_accepted_time_s=float(old['time_s']),matching='latest actual historical accepted state at or before 0.8 ms; no interpolation',comparison=comparison,first_events=first))
 dump(O/'CONTINUATION_DECISION.json',dict(decision='CONTINUE_SAME_CASE_TO_2P33MS',from_time_s=s['time'],reason='Normal mapping and unchanged numerical admission gates pass. A direction-controlled response is present; complete the frozen pulse and post-pulse interval to distinguish delayed fracture from no fracture. No amplitude/material/contact change.',physical_cases=1,restart_sha256=sha(target),state_sha256=sha(folder/'accepted.state.npz')))
 # Actually reload the saved 0.8 ms native checkpoint. Retain the already
 # accumulated Python work/defect ledger; do not restart it at zero.
 before=topology(ns);query=ns.ext.query;ns.close();replacement=Native('full/early_checkpoint_resume',restart=target,restart_state=s,dt=ns.dt,query=query);restored=replacement.dynamics();check=compare(s,restored);after=topology(replacement);check.update(native_topology_exact=np.array_equal(before[0],after[0]),special_exclusions_exact=before[1]==after[1],checkpoint_and_restore_wall_s=time.perf_counter()-start,same_physical_case=True,work_ledger_retained=m.work.copy())
 assert check['status']=='PASS' and check['native_topology_exact'] and check['special_exclusions_exact'],check
 dump(folder/'RESTORE_AUDIT.json',check);integ.ns=replacement;m.ns=replacement
 return replacement,restored
