from common import *
from native import Native
from mechanics import Mechanics
from journal import Journal
from integrator import Integrator
from checks import topology
import time,resource,traceback,shutil
from protocol_support import decorate,ComponentRecorder,early_checkpoint
freeze=json.loads((O/'LOAD_DIRECTION_FREEZE.json').read_text());assert sha(O/'LOAD_DIRECTION_FREEZE.json')==(O/'LOAD_DIRECTION_FREEZE.sha256').read_text().split()[0]
for name,h in freeze['driver_sha256'].items():assert sha(A/name)==h,name
for name,h in freeze['historical_driver_sha256'].items():assert sha(HIST_A/name)==h,name
assert sha(O/'REFERENCE_NORMAL_FIELD.npz')==freeze['reference_normal_sha256']
for name,h in freeze['input_sha256'].items():assert sha(O/'provenance'/name)==h,name
assert sha(LIB)==freeze['library_sha256']
with (O/'FORMAL_PILOT_LAUNCHED.json').open('x') as f:json.dump(dict(task=TASK,run_id=RUN,allowed_cases=1,from_real_initial_state=True,epoch=time.time(),freeze_sha256=sha(O/'LOAD_DIRECTION_FREEZE.json')),f)
start=time.perf_counter();cpu=time.process_time();ns=j=integ=m=None;status='NOT_STARTED';failure=None;first={k:None for k in ['bulk_softening','bulk_failure','interface_softening','interface_failure','detachment','particle_contact','wall_contact']};detach_times=np.full(1927,np.nan);drag_times=np.full(1927,np.nan);accepted=0;maxes={};newstates=[];next_fixed=.00025;final_match=False;early_done=False;recorder=None
try:
 ns=Native('full/native');s=ns.dynamics();assert np.array_equal(s['x'],ns.ref['X']) and not s['d'].any() and not s['velocity'].any()
 dump(O/'LOADED_NATIVE_IDENTITY.json',dict(paths=ns.maps,sha256=sha(LIB),MPI_ranks=ns.comm.size,threads=1,original_atom_ids_exact=True))
 m=Mechanics(ns,s,freeze['gates']);row,pairs,newsoft,newfail=m.check(s,s,0);decorate(row,s);j=Journal(list(row),ns.ref);j.accept(s,row,pairs=pairs);j.key(s,'initial');j.checkpoint(ns,s,'initial');integ=Integrator(ns,j);status='RUNNING';recorder=ComponentRecorder();recorder.accept(s,None,j)
 while s['time']<END-1e-15:
  prev=s;s,h=integ.next(prev);ns=integ.ns;m.ns=ns
  try:row,pairs,newsoft,newfail=m.check(prev,s,h)
  except BaseException:
   save_state(O/'full/unaccepted_mechanical_candidate.npz',s);raise
  decorate(row,s);ns.ext.accept(s);j.accept(s,row,prev,pairs,newsoft,newfail);recorder.accept(s,prev,j);accepted+=1
  b=ns.ref['bonds'];bt=ns.ref['bond_types'];trigger=[]
  for kind,ids in [('softening',newsoft),('failure',newfail)]:
   for typ,label in [(1,'bulk'),(2,'interface')]:
    key=label+'_'+kind;which=ids[bt[ids]==typ]
    if len(which) and first[key] is None:
     first[key]=dict(step=s['step'],time_s=s['time'],bond_ids=which+1,atom_ids=b[which]+1,reference_midpoints_m=ns.ref['X'][b[which]].mean(1),actual_midpoints_m=s['x'][b[which]].mean(1),kappa=s['kappa'][which],d=s['d'][which],stretch=np.linalg.norm(s['x'][b[which,1]]-s['x'][b[which,0]],axis=1)/ns.ref['r0'][which]-1);trigger.append('first_'+key)
  newdet=np.flatnonzero(s['detached']&~prev['detached']);detach_times[newdet]=s['time'];on=np.flatnonzero(s['valid_flow']&np.isnan(drag_times));drag_times[on]=s['time']
  if len(newdet):
   j.event('DETACHMENT',step=s['step'],time_s=s['time'],atom_ids=newdet+1,valid_drag=s['valid_flow'][newdet],components=s['component_sizes'][~s['component_attached']])
   if first['detachment'] is None:first['detachment']=dict(step=s['step'],time_s=s['time'],atom_ids=newdet+1,component_sizes=s['component_sizes'][~s['component_attached']]);trigger.append('first_detachment')
  for key,condition in [('particle_contact',bool(len(pairs) and row['max_contact_penetration_dx']>1e-12)),('wall_contact',row['wall_contact_nodes']>0)]:
   if condition and first[key] is None:first[key]=dict(step=s['step'],time_s=s['time']);trigger.append('first_'+key)
  if s['time']>=next_fixed-1e-15:trigger.append('fixed_observation');next_fixed+=.00025
  for reason in trigger:j.key(s,reason)
  if trigger:dump(O/'FIRST_EVENTS.json',first)
  for k in ['max_speed_m_s','max_dd','relative_R','max_contact_penetration_dx','max_wall_overlap_m','wall_contact_nodes','detached_nodes','detached_components','largest_detached_nodes','largest_detached_diameter_m','N_ge200','N_ge300','momentum_residual_Ns','x_equation_m','v_equation_m_s']:
   if k not in maxes or row[k]>maxes[k]['value']:maxes[k]=dict(value=row[k],step=s['step'],time_s=s['time'])
  if not early_done and s['time']>=EARLY-1e-15:
   ns,s=early_checkpoint(ns,s,m,integ,j,row,first);early_done=True;print('EARLY_CHECKPOINT_RESTORED_CONTINUING',s['time'],flush=True)
  if accepted%400==0:j.checkpoint(ns,s,'periodic_400_accepted');assert sum(p.stat().st_size for p in O.rglob('*') if p.is_file())<2*1024**3,'OUTPUT_LIMIT'
  if accepted%200==0 or trigger:
   dump(O/'PROGRESS.json',dict(status=status,step=s['step'],time_s=s['time'],wall_s=time.perf_counter()-start,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,bulk_failed=row['bulk_failed'],interface_failed=row['interface_failed'],detached_nodes=row['detached_nodes'],detached_components=row['detached_components'],max_relative_R_percent=maxes['relative_R']['value']*100,rejects=integ.rejects))
   print('ACCEPTED',s['step'],format(s['time']*1e3,'.6f'),'ms',row['bulk_failed'],row['interface_failed'],'free',row['detached_nodes'],'components',row['detached_components'],'R%',format(row['relative_R']*100,'.6g'),flush=True)
 codes,sp=topology(ns);alive=b[s['d']<1];expected=np.sort((alive[:,0]+1)*ns.N+(alive[:,1]+1));assert np.array_equal(codes,expected);final_match=True
 dump(O/'NATIVE_FINAL_TOPOLOGY_CHECK.json',dict(status='PASS',all_original_surviving_endpoint_pairs_exact=True,special_exclusions_sha256=sp,native_bond_count=len(codes)))
 j.key(s,'final');j.checkpoint(ns,s,'completed_protocol');status='FULL_PROTOCOL_NUMERICAL_PASS'
except BaseException as exc:
 status='STOPPED_NUMERICAL_OR_IMPLEMENTATION_GATE';failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc());print(failure['traceback'],flush=True)
 if j:j.failure(exc,failure['traceback'])
finally:
 if recorder:recorder.close()
 if integ:integ.save()
 if j:j.finish();shutil.copy2(O/'V3_NATIVE_PILOT_HISTORY.csv',O/'V3_DIRECTION_CONTROL_DYNAMIC_HISTORY.csv')
 final=None if j is None else j.last_row
 dump(O/'FIRST_EVENTS.json',first)
 np.savez_compressed(O/'DETACHMENT_DRAG_TIMES.npz',detach_time_s=detach_times,drag_time_s=drag_times)
 dump(O/'NUMERICAL_AUDIT.json',dict(status=status,accepted_steps=accepted,damage_gate_accepted=0 if integ is None else integ.accepts,rejected_damage_trials=0 if integ is None else integ.rejects,maxima=maxes,max_independent_errors={} if m is None else m.maxerr,final_native_topology_match=final_match,exception=failure,last_accepted=final))
 q=resource.getrusage(resource.RUSAGE_SELF);dump(O/'RUN_RESOURCES.json',dict(main_wall_s=time.perf_counter()-start,main_CPU_s=time.process_time()-cpu,peak_RSS_bytes=q.ru_maxrss*1024,MPI_ranks=1,threads=1,GPU_used=False,wall_limit_s=None,RSS_limit_GiB=None,native_candidate_wall_s=None if integ is None else integ.native_wall,restart_wall_s=None if integ is None else integ.restore_wall,checkpoints=None if j is None else j.checkpoints,checkpoint_wall_s=None if j is None else j.checkpoint_wall_s,output_bytes=sum(p.stat().st_size for p in O.rglob('*') if p.is_file())))
 dump(O/'FINAL_STATUS.json',dict(task=TASK,run_id=RUN,status=status,scope='EXPLORATORY_APPLICATION',scientific_label='EXPLORATORY_UNCALIBRATED',full_protocol_completed=status=='FULL_PROTOCOL_NUMERICAL_PASS',planned_end_s=END,last_accepted_time_s=None if final is None else final['time_s'],accepted_steps=accepted,formal_cases=1,INTERFACE_SOFTENING_ASSUMPTION='EXPLORATORY_UNCALIBRATED',GC_EFFECTIVE='UNAVAILABLE',CRITICAL_LOAD='UNAVAILABLE',macro_strength='UNVALIDATED',production_unchanged=True,delivery_status='POSTPROCESSING'))
 if ns:
  try:ns.close()
  except BaseException:pass
print('FORMAL_PILOT_FINISHED',status,flush=True)

raise SystemExit(0 if status=='FULL_PROTOCOL_NUMERICAL_PASS' else 1)
