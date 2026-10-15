from common import *
import gzip,csv,time,resource
class Journal:
 def __init__(self,columns,ref):
  self.folder=O/'full';self.folder.mkdir(exist_ok=True);(self.folder/'checkpoints').mkdir(exist_ok=True);(self.folder/'keyframes').mkdir(exist_ok=True);(self.folder/'rejected').mkdir(exist_ok=True)
  self.n=len(ref['X']);self.ref=ref;self.columns=columns;self.last=None;self.last_row=None;self.latest=None;self.attempt=None;self.records=0;self.checkpoints=0;self.checkpoint_wall_s=0.;self.keys=[];self.start=time.perf_counter()
  self.events=(self.folder/'events.jsonl').open('w');self.trials=(self.folder/'trials.jsonl').open('w');self.table=(O/'NATIVE_DYNAMIC_HISTORY.csv').open('w',newline='');self.writer=csv.DictWriter(self.table,fieldnames=columns);self.writer.writeheader()
  self.nodes=(self.folder/'native_frames.bin').open('wb');self.frame_index=(self.folder/'native_frames.csv').open('w',newline='');self.fw=csv.writer(self.frame_index);self.fw.writerow(['record','step','time_s']);self.frames=0
  self.damage=gzip.open(self.folder/'damage_updates.bin.gz','wb',compresslevel=1);self.contacts=gzip.open(self.folder/'contact_pairs.bin.gz','wb',compresslevel=1)
  self.damage_dtype=np.dtype([('step','<i4'),('bond_id','<i4'),('kappa','<f8'),('d','<f8')]);self.contact_dtype=np.dtype([('step','<i4'),('i','<i4'),('j','<i4'),('distance_m','<f8')])
  dump(self.folder/'STREAM_LAYOUT.json',dict(nodes_dtype='<f8',nodes_per_record_shape=[4,self.n,3],node_fields=['x','velocity','force','material_force'],record_index='native_frames.csv',sampling='Every 20 accepted steps plus final; NOT every-step nodal trajectory. Every accepted scalar/event is saved; full native kappa/d at selected keys and rolling checkpoints.',damage_dtype=self.damage_dtype.descr,damage_scope='Every accepted change in d with actual kappa; gzip binary; initial precut excluded. Subthreshold kappa only at full saved states.',contact_dtype=self.contact_dtype.descr))
 def line(self,f,z):f.write(json.dumps(z,default=cv,allow_nan=False,separators=(',',':'))+'\n');f.flush()
 def event(self,kind,**data):self.line(self.events,dict(kind=kind,**data))
 def trial_start(self,serial,s,h,depth,clearance):
  self.attempt=dict(trial=serial,last_accepted_step=s['step'],last_accepted_time_s=s['time'],candidate_step=s['step']+1,candidate_time_s=s['time']+h,h_s=h,depth=depth,predicted_domain_clearance_m=clearance);self.line(self.trials,dict(kind='START',**self.attempt))
 def frame(self,s):
  self.nodes.write(np.asarray([s[k] for k in ['x','velocity','force','material_force']],dtype='<f8').tobytes());self.fw.writerow([self.frames,s['step'],s['time']]);self.frames+=1;self.nodes.flush();self.frame_index.flush()
 def key(self,s,reason):
  path=self.folder/'keyframes'/f'step_{s["step"]:06d}.npz'
  if not path.exists():save_state(path,s)
  record=dict(step=s['step'],time_s=s['time'],reason=reason,path=str(path.relative_to(self.folder)));self.keys.append(record);dump(self.folder/'KEYFRAME_INDEX.json',self.keys)
 def accept(self,s,row,prev=None,pairs=None,newsoft=None,newfail=None):
  assert np.array_equal(s['atom_ids'],np.arange(1,self.n+1))
  if prev is not None:
   changes=np.flatnonzero(s['d']!=prev['d']);z=np.empty(len(changes),self.damage_dtype);z['step']=s['step'];z['bond_id']=changes+1;z['kappa']=s['kappa'][changes];z['d']=s['d'][changes];self.damage.write(z.tobytes())
  if pairs is not None and len(pairs):
   z=np.empty(len(pairs),self.contact_dtype);z['step']=s['step'];z['i']=pairs[:,0]+1;z['j']=pairs[:,1]+1;z['distance_m']=np.linalg.norm(s['x'][pairs[:,1]]-s['x'][pairs[:,0]],axis=1);self.contacts.write(z.tobytes())
  self.writer.writerow(row);self.table.flush();self.last=s;self.last_row=row;self.records+=1
  self.event('ACCEPTED',step=s['step'],time_s=s['time'],atom_count=self.n,new_soft_ids=[] if newsoft is None else newsoft+1,new_failed_ids=[] if newfail is None else newfail+1,native_bonds=s['native_bonds'])
  if s['step']%20==0:self.frame(s)
  dump(self.folder/'LAST_ACCEPTED.json',dict(step=s['step'],time_s=s['time'],ledger=row,nearest_checkpoint=self.latest))
 def checkpoint(self,ns,s,reason):
  start=time.perf_counter();slot=self.checkpoints%2;path=self.folder/'checkpoints'/f'accepted_slot_{slot}.restart';tmp=Path(str(path)+'.pending');ns.l.command(f'write_restart {tmp}');os.replace(tmp,path);state=Path(str(path)+'.state.npz');save_state(state,s)
  self.checkpoints+=1;self.latest=dict(restart=str(path),restart_sha256=sha(path),state=str(state),state_sha256=sha(state),step=s['step'],time_s=s['time'],reason=reason);dump(self.folder/'LATEST_NATIVE_CHECKPOINT.json',self.latest);self.event('CHECKPOINT',**self.latest);self.checkpoint_wall_s+=time.perf_counter()-start
  self.damage.flush();self.contacts.flush();self.events.flush();self.trials.flush();return path
 def failure(self,exc,trace=''):
  if self.last is not None:save_state(self.folder/'last_accepted.npz',self.last)
  dump(self.folder/'FAILURE_RECORD.json',dict(status='FAILED',exception_type=type(exc).__name__,exception=str(exc),attempt=self.attempt,last_accepted_step=None if self.last is None else self.last['step'],last_accepted_time_s=None if self.last is None else self.last['time'],nearest_native_checkpoint=self.latest,native_calls_after_exception=0,candidate_not_admitted=True,unavailable='Unaccepted native candidate deliberately not queried',peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024));(self.folder/'EXCEPTION.txt').write_text(trace);self.event('FAILURE',message=str(exc))
 def finish(self):
  if self.last is not None:
   save_state(self.folder/'last_accepted.npz',self.last)
   if self.last['step']%20:self.frame(self.last)
  for f in [self.damage,self.contacts,self.nodes,self.frame_index,self.table,self.events,self.trials]:f.flush();f.close()
