"""Read-only diagnostics. The only solver mutations are the unmodified original callbacks.
Native END_OF_STEP supplies actual accepted x/v/f (not reconstructed half velocities).
"""
from common import *
import numpy as np,ctypes as ct,json,csv,traceback,os,time
from pd_clot.lammps_fragmentation import early_flow_coupling as E
from pd_clot.lammps_fragmentation.fragment_semantics import component_metrics
DT=E.DT;V=(50e-6)**3;M=E.MASS
class Flow(E.EarlyFlow):
 def __init__(self,*args):
  super().__init__(*args)
  if self.comm.rank==0:
   self.file.close();self.fullfile.close();self.file=open(os.devnull,'wb');self.fullfile=open(os.devnull,'wb')
 def compute(self,step,nlocal,tags,x,f):
  self.half_v=self.lmp.numpy.extract_atom('v')[:nlocal].copy()
  self.before_flow=self.lmp.numpy.extract_atom('f')[:nlocal].copy()
  super().compute(step,nlocal,tags,x,f)
  self.wall.file.flush()
class Audit:
 def __init__(self,L,comm,flow,wall):
  self.L=L;self.comm=comm;self.rank=comm.rank;self.flow=flow;self.wall=wall;self.g=flow.g
  self.X=self.g['X'];self.b=self.g['bonds'];self.types=self.g['bond_types'];self.N=len(self.X);self.n=flow.n;self.keys=flow.keys
  r=np.linalg.norm(self.X[self.b[:,1]]-self.X[self.b[:,0]],axis=1);h=np.where(abs(r-150.5e-6)<=25e-6,(175.5e-6-r)/50e-6,1.)
  self.m=np.bincount(self.b.ravel(),weights=np.repeat(r*h*V,2),minlength=self.N);assert (self.m>0).all()
  self.probe=ct.CDLL(str(O/'observer.so'));self.probe.observe.argtypes=[ct.c_void_p,ct.c_int];self.probe.observe.restype=ct.POINTER(ct.c_double)
  self.probe.counter.argtypes=[ct.c_void_p,ct.c_int];self.probe.counter.restype=ct.c_longlong
  self.probe.permutation.argtypes=[ct.c_void_p];self.probe.permutation.restype=ct.POINTER(ct.c_int)
  self.probe.endpoint_records.argtypes=[ct.c_void_p,ct.c_void_p,ct.c_int];self.probe.endpoint_records.restype=ct.c_int
  self.probe.contact_readout.argtypes=[ct.c_void_p,ct.c_void_p,ct.c_void_p,ct.c_int];self.probe.contact_readout.restype=ct.c_int
  self.buf=np.empty((2*len(self.b),3),np.int32)
  self.cbtype=ct.CFUNCTYPE(None,ct.c_longlong,ct.c_int);self.callback=self.cbtype(self.safe);self.probe.register_observer.argtypes=[self.cbtype];self.probe.register_observer(self.callback)
  self.last_ids=None;self.lastghost=None;self.lastsort=None;self.lastmax=None;self.prevalive=np.ones(len(self.b),bool);self.prevowner=None;self.prevx=None;self.first_flags=set();self.lastcomplete=-1
  self.rankfile=(P/f'identity_rank{self.rank}.csv').open('w',buffering=1);self.rankcsv=csv.writer(self.rankfile);self.rankcsv.writerow(['step','phase','nlocal','nghost','nmax','nextsort','neighbor_ago','sort_called','nonidentity_sort_count','local_order_changed','incoming_owned','ghost_list_changed','m_abs_owned','m_rel_owned','m_abs_ghost','m_rel_ghost','worst_id_owned','worst_id_ghost','owned_endpoint_records'])
  self.eventfile=(P/f'identity_events_rank{self.rank}.jsonl').open('w',buffering=1)
  if not self.rank:
   self.scalarf=(P/'accepted_step_scalars.csv').open('w',buffering=1);self.scalarw=None
   self.breakf=(P/'break_events.csv').open('w',buffering=1);self.breakw=csv.writer(self.breakf);self.breakw.writerow(['step','initial_bond_index','id1','id2','type','owner1','owner2','cross_rank','index1_rank0','index2_rank0','index1_rank1','index2_rank1','nlocal_rank0','nlocal_rank1'])
   self.comparef=(P/'old_new_force_comparison.csv').open('w',buffering=1);self.comparew=None
   self.old=np.memmap(OLD/'saved_frame_force_audit.bin',dtype=E.DTYPE,mode='r');assert len(self.old)%self.n==0;self.old=self.old.reshape(-1,self.n);self.oldsteps=self.old['step'][:,0].copy();assert np.all(np.diff(self.oldsteps)>0);self.oldlookup={int(s):i for i,s in enumerate(self.oldsteps)}
   (P/'states').mkdir(exist_ok=False);np.savez_compressed(P/'reference.npz',X=self.X,bonds=self.b,bond_types=self.types,reference_m=self.m,clot_count=self.n,volume_m3=V,mass_kg=M)
 def safe(self,step,phase):
  try:self.observe(int(step),int(phase))
  except BaseException:
   dump(P/f'OBSERVATION_FAILURE_rank{self.rank}.json',{'step':int(step),'phase':int(phase),'last_complete':self.lastcomplete,'traceback':traceback.format_exc()})
   self.rankfile.flush();self.eventfile.flush();self.comm.Abort(51)
 def observe(self,step,phase):
  L=self.L;nl=L.extract_global('nlocal');ng=L.extract_global('nghost');nt=nl+ng
  allid=L.numpy.extract_atom('id')[:nt].astype(int).copy();ids=allid[:nl];ghost=allid[nl:];assert np.all((allid>=1)&(allid<=self.N))
  mm=np.ctypeslib.as_array(self.probe.observe(L.lmp,0),shape=(nt,)).copy();theta=np.ctypeslib.as_array(self.probe.observe(L.lmp,1),shape=(nt,)).copy();err=abs(mm-self.m[allid-1]);rel=err/self.m[allid-1]
  cs=[int(self.probe.counter(L.lmp,k)) for k in range(5)];nextsort,ago,nmax,nbins,sortfreq=cs;assert sortfreq==1000
  sort_called=self.lastsort is not None and nextsort!=self.lastsort
  initialsort=self.lastsort is None
  perm=np.ctypeslib.as_array(self.probe.permutation(L.lmp),shape=(nl,)).copy() if (sort_called or initialsort) and nbins>1 else np.arange(nl)
  assert np.array_equal(np.sort(perm),np.arange(nl))
  sortedmask=perm!=np.arange(nl);before_sort=ids.copy();before_sort[perm]=ids
  changed=self.last_ids is not None and not np.array_equal(ids,self.last_ids);incoming=[] if self.last_ids is None else np.setdiff1d(ids,self.last_ids).tolist();outgoing=[] if self.last_ids is None else np.setdiff1d(self.last_ids,ids).tolist()
  ghostchanged=self.lastghost is not None and not np.array_equal(ghost,self.lastghost);growth=self.lastmax is not None and nmax!=self.lastmax
  event=(changed or ghostchanged or sort_called or initialsort or growth or ago==0)
  if event:
   self.eventfile.write(json.dumps({'step':step,'rank':self.rank,'nextsort':nextsort,'sort_called':bool(sort_called),'initial_sort_observed':initialsort,'nonidentity_sort_count':int(sortedmask.sum()),'sorted_ids':ids[sortedmask].tolist(),'pre_sort_ids':before_sort.tolist() if sort_called or initialsort else None,'post_ids':ids.tolist(),'previous_ids':None if self.last_ids is None else self.last_ids.tolist(),'incoming_owned_ids':incoming,'outgoing_owned_ids':outgoing,'ghost_ids':ghost.tolist(),'ghost_list_changed':ghostchanged,'neighbor_ago':ago,'nmax':nmax,'storage_growth':growth,'m_abs_max':float(err.max()),'m_relative_max':float(rel.max()),'worst_id':int(allid[np.argmax(rel)])},allow_nan=False)+'\n')
  assert np.isfinite(mm).all() and np.max(rel)<=2e-12,('m_identity',step,self.rank,int(allid[np.argmax(rel)]),float(rel.max()))
  alive=self.flow.alive;broken=np.flatnonzero(self.prevalive&~alive);changed_top=bool(len(broken)) or step==0
  endpoint_n=-1
  if changed_top:
   endpoint_n=int(self.probe.endpoint_records(L.lmp,self.buf.ctypes.data,len(self.buf)));assert endpoint_n>=0,('illegal_owned_index',endpoint_n)
   ar=self.buf[:endpoint_n];assert (ar[:,:2]>0).all() and (ar[:,:2]<=self.N).all() and np.isin(ar[:,0],ids).all()
   ends=np.sort(ar[:,:2].astype(np.int64)-1,axis=1);key=ends[:,0]*self.N+ends[:,1];idx=np.searchsorted(self.keys,key)
   assert np.all(idx<len(self.b)) and np.array_equal(self.keys[idx],key) and np.array_equal(self.types[idx],ar[:,2]) and alive[idx].all()
   assert len(np.unique((ar[:,0].astype(np.int64)-1)*self.N+ar[:,1]-1))==endpoint_n,'duplicate owned endpoint'
   expected_degree=np.bincount(self.b[alive].ravel(),minlength=self.N);actualdegree=np.bincount(ar[:,0]-1,minlength=self.N)
   assert np.array_equal(actualdegree[ids-1],expected_degree[ids-1]),'missing owned endpoint'
  # All numerical arrays below are copied, never written through LAMMPS pointers.
  x=L.numpy.extract_atom('x')[:nl].copy();v=L.numpy.extract_atom('v')[:nl].copy();force=L.numpy.extract_atom('f')[:nl].copy()
  drag=np.column_stack([L.numpy.extract_atom('d_drag'+a)[:nl].copy() for a in 'xyz']);numeric=np.column_stack([L.numpy.extract_atom('d_num'+a)[:nl].copy() for a in 'xyz'])
  block=np.column_stack([x,v,force,self.flow.half_v,self.flow.before_flow,drag,numeric,self.wall.last_force,mm[:nl],theta[:nl]])
  assert np.isfinite(block).all(),'nonfinite accepted state'
  indexbyid=np.full(self.N,-1,int);indexbyid[allid-1]=np.arange(nt)
  localmeta={'rank':self.rank,'nlocal':nl,'nghost':ng,'nmax':nmax,'ago':ago,'sort_called':sort_called,'sorted_count':int(sortedmask.sum()),'incoming':incoming,'endpoint_count':endpoint_n,'m_rel':float(rel.max()),'m_abs':float(err.max()),'indexbyid':indexbyid}
  contactbuf=np.zeros((nt,3));contactends=np.empty((len(self.b),2),np.int32)
  ncontact=int(self.probe.contact_readout(L.lmp,contactbuf.ctypes.data,contactends.ctypes.data,len(contactends)));assert ncontact>=0
  contact_global=np.zeros((self.N,3));np.add.at(contact_global,allid-1,contactbuf*V)
  packets=self.comm.gather((ids,block,localmeta,contact_global,contactends[:ncontact].copy()),root=0)
  self.rankcsv.writerow([step,phase,nl,ng,nmax,nextsort,ago,int(sort_called),int(sortedmask.sum()),int(changed),len(incoming),int(ghostchanged),float(err[:nl].max()),float(rel[:nl].max()),float(err[nl:].max(initial=0)),float(rel[nl:].max(initial=0)),int(ids[np.argmax(rel[:nl])]),int(ghost[np.argmax(rel[nl:])]) if ng else -1,endpoint_n])
  if not self.rank:self.global_record(step,phase,packets,broken,changed_top)
  self.comm.Barrier() # Global audit for this accepted step completed and flushed before advancing.
  self.last_ids=ids;self.lastghost=ghost;self.lastsort=nextsort;self.lastmax=nmax;self.prevalive=alive.copy();self.lastcomplete=step
 def global_record(self,step,phase,packets,broken,changed_top):
  allids=np.concatenate([p[0] for p in packets]);assert np.array_equal(np.sort(allids),np.arange(1,self.N+1)),'ownership coverage'
  block=np.zeros((self.N,26));owner=np.zeros(self.N,int)
  for ids,bb,meta,_,_ in packets:block[ids-1]=bb;owner[ids-1]=meta['rank']
  contact_force=sum(p[3] for p in packets);contact_edges=np.concatenate([p[4] for p in packets])
  x=block[:,:3];v=block[:,3:6];force=block[:,6:9];half=block[:,9:12];before=block[:,12:15];drag=block[:,15:18];numeric=block[:,18:21];wall=block[:,21:24];m=block[:,24];theta=block[:,25]
  if changed_top:assert sum(p[2]['endpoint_count'] for p in packets)==2*self.flow.nbonds
  moving=np.flatnonzero(owner!=self.prevowner) if self.prevowner is not None else np.array([],int)
  for ib in broken:
   a,z=self.b[ib];self.breakw.writerow([step,int(ib),int(a+1),int(z+1),int(self.types[ib]),int(owner[a]),int(owner[z]),int(owner[a]!=owner[z]),int(packets[0][2]['indexbyid'][a]),int(packets[0][2]['indexbyid'][z]),int(packets[1][2]['indexbyid'][a]),int(packets[1][2]['indexbyid'][z]),packets[0][2]['nlocal'],packets[1][2]['nlocal']])
  det=self.flow.detached;sz=self.flow.sizes;att=self.flow.attached
  assert int(sz.sum())==self.n and int(sz[~att].sum())==int(det.sum())
  assert not np.any(drag[:self.n][~det]) and not np.any(numeric[:self.n][det])
  speed=np.linalg.norm(v[:self.n],axis=1);peak=int(np.argmax(speed));dstep=np.linalg.norm(x[:self.n]-self.prevx[:self.n],axis=1) if self.prevx is not None else np.zeros(self.n)
  free=sz[~att];largest=int(free.max(initial=0));tiny=int(free[free<=2].sum())
  row=dict(step=step,time_s=step*DT,accepted_after_final_integrate=phase,native_bonds=int(self.flow.nbonds),broken_bulk=int(np.sum(~self.flow.alive&(self.types==1))),broken_interface=int(np.sum(~self.flow.alive&(self.types==2))),broken_anchor=int(np.sum(~self.flow.alive&(self.types==3))),new_breaks=len(broken),new_cross_rank_breaks=int(np.sum(owner[self.b[broken,0]]!=owner[self.b[broken,1]])),detached_nodes=int(det.sum()),detached_components=len(free),attached_nodes=int((~det).sum()),tiny_detached_nodes=tiny,resolved_detached_components=int(np.sum(free>=3)),largest_detached_nodes=largest,largest_detached_diameter_m=float((6*largest*V/np.pi)**(1/3)),N_ge200=int(np.sum(free*V>=np.pi*(200e-6)**3/6)),N_ge300=int(np.sum(free*V>=np.pi*(300e-6)**3/6)),vmax_m_s=float(speed[peak]),vmax_atom_id=peak+1,max_step_displacement_m=float(dstep.max()),max_displacement_atom_id=int(np.argmax(dstep))+1,KE_material_J=float(.5*M*np.sum(v[:self.n]**2)),KE_detached_J=float(.5*M*np.sum(v[:self.n][det]**2)),drag_norm_sum_N=float(np.linalg.norm(drag,axis=1).sum()),numerical_norm_sum_N=float(np.linalg.norm(numeric,axis=1).sum()),wall_norm_sum_N=float(np.linalg.norm(wall,axis=1).sum()),drag_power_halfvelocity_W=float(np.sum(drag*half)),numerical_power_halfvelocity_W=float(np.sum(numeric*half)),wall_power_halfvelocity_W=float(np.sum(wall*half)),m_max_relative=max(p[2]['m_rel'] for p in packets),m_max_absolute=max(p[2]['m_abs'] for p in packets),ownership_changes=len(moving),sort_ranks=sum(p[2]['sort_called'] for p in packets),sort_nonidentity_nodes=sum(p[2]['sorted_count'] for p in packets),neighbor_rebuild_ranks=sum(p[2]['ago']==0 for p in packets))
  row.update(contact_pairs=len(contact_edges),contact_force_norm_sum_N=float(np.linalg.norm(contact_force[:self.n],axis=1).sum()),contact_power_fullvelocity_W=float(np.sum(contact_force*v)))
  row.update({k:self.flow.history[-1][k] for k in ['valid_drag_nodes','invalid_detached_nodes','drag_on_attached_nodes','numerical_damping_on_detached_nodes']})
  if self.scalarw is None:self.scalarw=csv.DictWriter(self.scalarf,fieldnames=list(row));self.scalarw.writeheader()
  self.scalarw.writerow(row)
  reasons=[]
  if step in [0,240,241,320,455,593,621,1000,2000,4660] or step%20==0 and step<=700:reasons.append('requested_or_early_alignment')
  flags={'first_bulk_break':bool(np.any(self.types[broken]==1)),'first_interface_break':bool(np.any(self.types[broken]==2)),'first_detachment':bool(det.any()),'first_migration':bool(len(moving)),'first_cross_rank_break':bool(row['new_cross_rank_breaks']),'first_ge200':bool(row['N_ge200']),'first_ge300':bool(row['N_ge300'])}
  for k,yes in flags.items():
   if yes and k not in self.first_flags:self.first_flags.add(k);reasons.append(k)
  if row['sort_ranks']:reasons.append('actual_sort')
  if len(moving):reasons.append('actual_migration')
  if len(contact_edges) and 'first_particle_contact' not in self.first_flags:self.first_flags.add('first_particle_contact');reasons.append('first_particle_contact')
  if row['wall_norm_sum_N'] and 'first_wall_contact' not in self.first_flags:self.first_flags.add('first_wall_contact');reasons.append('first_wall_contact')
  if step in self.oldlookup:
   old=self.old[self.oldlookup[step]];assert np.all(old['step']==step) and np.array_equal(old['atom_id'],np.arange(1,self.n+1))
   cr={'step':step,'time_s':step*DT}
   for k,arr in [('x_m',x[:self.n]),('v_force_half_m_s',half[:self.n]),('force_before_flow_N',before[:self.n]),('drag_N',drag[:self.n]),('numerical_N',numeric[:self.n]),('wall_N',wall[:self.n])]:
    dif=np.linalg.norm(arr-old[k],axis=1);cr[k+'_max_delta']=float(dif.max());cr[k+'_worst_id']=int(np.argmax(dif))+1;cr[k+'_old_scale']=float(np.linalg.norm(old[k],axis=1).max());cr[k+'_bitwise_equal']=bool(np.array_equal(arr,old[k]))
   cr['different_component_nodes']=int(np.count_nonzero(old['component_id']!=self.flow.cid[self.flow.lab]));cr['different_detachment_nodes']=int(np.count_nonzero((old['detach_step']>=0)!=det))
   if self.comparew is None:self.comparew=csv.DictWriter(self.comparef,fieldnames=list(cr));self.comparew.writeheader()
   self.comparew.writerow(cr)
   meaningful=cr['x_m_max_delta']>5e-12 or cr['v_force_half_m_s_max_delta']>1e-8 or cr['force_before_flow_N_max_delta']>1e-12
   if meaningful and 'first_meaningful_saved_difference' not in self.first_flags:
    self.first_flags.add('first_meaningful_saved_difference');reasons.append('first_meaningful_saved_difference');dump(P/'FIRST_SAVED_DIFFERENCE.json',cr)
  if reasons:
   np.savez_compressed(P/'states'/f'step_{step:05d}.npz',step=step,position_m=x,velocity_m_s=v,force_N=force,velocity_half_m_s=half,force_before_flow_N=before,drag_N=drag,numerical_N=numeric,wall_N=wall,m=m,theta=theta,owner=owner,alive=self.flow.alive,component_label=self.flow.lab,component_sizes=sz,component_attached=att,detached=det,migrating_ids=moving+1,contact_force_N=contact_force,contact_edges=contact_edges,reasons=np.array(reasons))
  if step%100==0 or reasons:dump(P/'LAST_COMPLETED_STEP.json',{'step':step,'time_s':step*DT,'phase':phase,'m_max_relative':row['m_max_relative'],'native_bonds':row['native_bonds']})
  self.prevowner=owner.copy();self.prevx=x.copy()
 def close(self):
  self.rankfile.close();self.eventfile.close()
  if not self.rank:
   self.scalarf.close();self.breakf.close();self.comparef.close()
