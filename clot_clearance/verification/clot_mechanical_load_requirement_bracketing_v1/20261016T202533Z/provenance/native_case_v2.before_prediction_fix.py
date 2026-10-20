"""Native small-specimen subprocess. Reuses immutable constitutive/native code.

The existing precracked small geometry is copied byte-for-byte. The old top
displacement grip is replaced by physical shear traction on the flat top.
All material nodes are integrated; only the original bottom anchors are fixed.
"""
from common import *
from prescribed_traction import Load, END, RISE
import gzip
import ast,ctypes,types,shutil,time,resource,traceback,csv,argparse
from mpi4py import MPI
from lammps import lammps
graph=module('accepted_graph_ops',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')
def native_factory(case,flow):
 source=R/'analysis/v3_softening_pilot_v1/native.py'
 tree=ast.parse(source.read_text());classes=[x for x in tree.body if isinstance(x,ast.ClassDef) and x.name in ['Native','Contact']]
 namespace=dict(globals(),O=case);exec(compile(ast.Module(body=classes,type_ignores=[]),str(source),'exec'),namespace)
 Base=namespace['Native'];Contact=namespace['Contact']
 class Native(Base):
  def __init__(self):
   super().__init__('native',dt=BASE)
   d=self.ref['initial_cut'].astype(float);k=SC[self.ref['bond_types']]*d
   assert self.lib.lpsvar_action(self.l.lmp,1,0,self.ptr(k),self.ptr(d))==0
  def dynamics(self):
   self.l.commands_string(f'group clot type 1\nfix integrate clot nve\nfix coordinator all LPS/VARIATIONAL/DYNAMIC {RHO:.17g}\nfix v3forces clot external pf/callback 1 1')
   self.ext=types.SimpleNamespace(error=None,last={})
   def callback(caller,step,nlocal,tags,x,f):
    f[:]=0
    try:
     ids=tags.astype(int)-1;assert nlocal==self.N
     t=float(self.l.get_thermo('time'));force=flow.force(t);f[:]=force[ids]
     self.ext.last={'external_force':force,'external_time':np.array(t)}
    except BaseException:self.ext.error=traceback.format_exc();f[:]=np.nan
   self.callback=callback;self.l.set_fix_external_callback('v3forces',self.callback);self.dynamic=True;self.l.command('run 0');return self.state()
 return Native,Contact
def predict(model,ref,x,k):
 e,_=model.geometry(x);k=np.maximum(k,np.maximum(e/ref['r0'],0));on=S0[ref['bond_types']];off=SC[ref['bond_types']];d=np.zeros(len(k));d[k>=off]=1;part=(k>on)&(k<off);d[part]=off[part]*(k[part]-on[part])/(k[part]*(off[part]-on[part]));return k,d
def components(ref,s,cached):
 labels,sizes,attached=cached;out=[];bt=ref['bond_types'];cut=ref['initial_cut']
 for q in np.flatnonzero(~attached):
  ids=np.flatnonzero(labels==q);_,inside=graph.boundary(ref,ids);bulk=inside[(bt[inside]==1)&~cut[inside]]
  out.append(dict(atom_ids=ids+1,nodes=len(ids),equivalent_diameter_m=(6*len(ids)*V/np.pi)**(1/3),mass_kg=len(ids)*RHO*V,internal_alive_bonds=int(np.sum(s['d'][bulk]<1)),internal_initial_alive_bonds=len(bulk),internal_alive_fraction=float(np.mean(s['d'][bulk]<1)) if len(bulk) else None))
 return out

def run(case,params):
 start=time.perf_counter();cpu=time.process_time();case=Path(case)
 if (case/'SUMMARY.json').exists():return json.loads((case/'SUMMARY.json').read_text())
 assert not (case/'history.csv').exists(),'Interrupted case retained; use a new attempt directory'
 (case/'provenance').mkdir(parents=True,exist_ok=True);(case/'keys').mkdir(exist_ok=True)
 for name in ['reference.npz','reference.native','body.data']:
  target=case/'provenance'/name
  if not target.exists():target.symlink_to(Path('../../../provenance')/name)
 ref=load_npz(O/'provenance/reference.npz');flow=Load(ref,params)
 dump(case/'REQUEST.json',dict(parameters=params,load_audit=flow.audit(ref),startup_recovery_sha256=sha(O/'STARTUP_RECOVERY_FREEZE.json'),executed_driver_sha256=sha(Path(__file__)),freeze_sha256=sha(O/'BASELINE_AND_SPECIMEN_FREEZE.json'),native_sha256=sha(LIB),label='PRESCRIBED_TRACTION_MECHANICAL_REQUIREMENT',boundary='IDEAL_MECHANICAL_BOUNDARY' if params['mode']=='IDEAL_LOCAL_PULL' else 'NOMINAL_DEAD_TRACTION'))
 Native,Contact=native_factory(case,flow);ns=Native();n=ns.n;b=ref['bonds'];bt=ref['bond_types'];cut=ref['initial_cut'];model=energy_model(ref);contact=Contact(ref)
 s=ns.dynamics();initial=s;prev=s;W=0.;Jc=0.;impulse=np.zeros(3);rows=[];err={};maxfree=0;events={};substeps=0;rejections=0;scale=1e-24;stop=[];secondary=[];regimes=['NO_DAMAGE'];hist=[]
 def error(key,val):err[key]=max(err.get(key,0.),float(val))
 def save(state,name):np.savez_compressed(case/'keys'/name,**state)
 def first(name,ids,state):
  if name in events or not len(ids):return
  ends=b[ids];delta=state['x'][ends[:,1]]-state['x'][ends[:,0]];stretch=np.linalg.norm(delta,axis=1)/ref['r0'][ids]-1
  X=ref['X'];mid=X[ends].mean(1);limits=np.array([X[:n].min(0),X[:n].max(0)])
  events[name]=dict(time_s=state['time'],step=state['step'],original_bond_ids=ids+1,atom_ids=ends+1,bond_types=bt[ids],reference_endpoints_m=X[ends],actual_endpoints_m=state['x'][ends],stretch=stretch,kappa=state['kappa'][ids],d=state['d'][ids],near_anchor=(mid[:,2]-X[:n,2].min())<=HORIZON,near_material_bbox_surface=np.min(np.minimum(abs(mid-limits[0]),abs(mid-limits[1])),axis=1)<=DX,evidence=f'keys/{name}.npz')
  save(state,name+'.npz')
 save(s,'initial.npz');stream=(case/'history.csv').open('w');writer=None;changedfile=gzip.open(case/'damage_events.jsonl.gz','wt')
 cached=graph.topology(ref,s['d']);oldfree=[]
 try:
  while prev['time']<END-1e-15:
   # Resolve the near-zero-energy smooth startup without relaxing the 1% gate.
   startup_limit=BASE/8 if prev['time']<2e-6-1e-15 else BASE/4 if prev['time']<4e-6-1e-15 else BASE/2 if prev['time']<8e-6-1e-15 else BASE
   h=min(startup_limit,END-prev['time']);depth=0
   while True:
    xp=prev['x'].copy();xp[:n]+=h*prev['velocity'][:n]+h*h/(2*RHO*V)*prev['force'][:n]
    kk,dp=predict(model,ref,xp,prev['kappa'])
    if np.max(dp-prev['d'])<=.1:break
    assert h>BASE/8*(1+1e-12),'MINIMUM_EVENT_SUBSTEP_EXHAUSTED'
    h/=2;depth+=1;rejections+=1
   substeps+=bool(depth)
   if h!=ns.dt:ns.l.command(f'timestep {h:.17g}');ns.dt=h
   ns.l.command('run 1 pre no post no');s=ns.state();p=prev;dx=s['x']-p['x'];dd=s['d']-p['d'];sc=s['scalars'];nf=np.flatnonzero((s['d']==1)&(p['d']<1)&~cut)
   error('predicted_x_m',np.max(abs(xp-s['x'])));error('predicted_damage',np.max(abs(dp-s['d'])));error('predicted_kappa',np.max(abs(kk-s['kappa'])));error('irreversibility',max(np.max(-dd),np.max(p['kappa']-s['kappa'])))
   error('v_VV_m_s',np.max(abs(s['velocity'][:n]-p['velocity'][:n]-h/(2*RHO*V)*(p['force'][:n]+s['force'][:n]))))
   assert err['predicted_x_m']<2e-18 and err['v_VV_m_s']<1e-10 and max(err['predicted_damage'],err['predicted_kappa'],err['irreversibility'])<1e-10
   assert np.max(dd)<=.1+1e-10 and np.array_equal(s['x'][n:],ref['X'][n:]) and not s['velocity'][n:].any()
   assert s['native_bonds']==int(np.sum(s['d']<1)) and sc[3]==int(np.sum((s['d']==1)&~cut)) and not s['d'][bt==3].any()
   assert np.isfinite(s['force']).all() and np.array_equal(s['atom_ids'],np.arange(1,ns.N+1))
   wnew=W+float(.5*np.sum((p['external_force']+s['external_force'])*dx));jnew=Jc
   if len(nf):
    rr=np.linalg.norm(s['x'][b[nf,1]]-s['x'][b[nf,0]],axis=1);jnew+=float(.5*contact.c*np.sum(np.minimum(rr-DX,0)**2))
   res=wnew-(sc[0]-initial['scalars'][0]+s['K']-initial['K']+sc[1]+s['Uc']-initial['Uc']-jnew);newscale=max(scale,abs(wnew),abs(sc[0]-initial['scalars'][0])+s['K']+sc[1]+s['Uc']);rel=abs(res)/newscale
   assert rel<=.01,('ENERGY_RESIDUAL',rel,res,newscale)
   if s['step']%100==0 or np.any(dd):
    new=model.evaluate(s['x'],s['d']);old=model.evaluate(s['x'],p['d']) if np.any(dd) else new;uc,fc,pairs=contact.evaluate(s['x'],s['d'])
    error('native_force_N',np.max(abs(new['force']-s['material_force'])));error('external_port_N',np.max(abs(s['force']-s['material_force']-fc-s['external_force'])));error('native_U_J',abs(new['U']-sc[0]));error('dissipation_increment_J',abs(old['U']-new['U']-(sc[1]-p['scalars'][1])));error('contact_U_J',abs(uc-s['Uc']))
    assert max(err['native_force_N'],err['external_port_N'])<1e-12 and max(err['native_U_J'],err['dissipation_increment_J'],err['contact_U_J'])<1e-17
   # All gates passed: commit the Python ledger and evidence for this native state.
   W=wnew;Jc=jnew;scale=newscale;impulse+=.5*h*(p['external_force'].sum(0)+s['external_force'].sum(0))
   first('first_kappa_exceeds_s0',np.flatnonzero((s['kappa']>S0[bt])&~cut),s)
   first('first_softening',np.flatnonzero((s['d']>0)&~cut),s);first('first_failure',nf,s)
   if len(nf):cached=graph.topology(ref,s['d'])
   labels,sizes,attached=cached;free=sizes[~attached];largest=int(free.max(initial=0));natural=(s['d']==1)&~cut
   if len(nf):
    cs=components(ref,s,cached);sets=[set(c['atom_ids'].tolist()) for c in cs]
    for before in oldfree:
     children=[sorted(before & q) for q in sets if before & q]
     if len(children)>1:
      secondary.append(dict(time_s=s['time'],step=s['step'],parent_atom_ids=sorted(before),children_atom_ids=children));
      if len(secondary)==1:save(s,'first_secondary.npz')
    oldfree=sets
    hist.append(dict(time_s=s['time'],step=s['step'],components=cs))
   if largest>maxfree:maxfree=largest;save(s,'max_free.npz')
   if largest>=2 and 'first_multinode' not in events:
    events['first_multinode']=dict(time_s=s['time'],step=s['step'],components=components(ref,s,cached),evidence='keys/first_multinode.npz');save(s,'first_multinode.npz')
   classification='SECONDARY_FRAGMENTATION' if secondary else 'COHESIVE_MULTINODE_DETACHMENT' if largest>=2 else 'DETACHED_SINGLETON_DEBRIS' if largest else 'BOND_FAILURE_NO_DETACHMENT' if natural.any() else 'SOFTENING_ONLY' if np.any((s['d']>0)&~cut) else 'NO_DAMAGE'
   if classification!=regimes[-1]:regimes.append(classification)
   diameter=(6*free*V/np.pi)**(1/3)
   row=dict(step=s['step'],time_s=s['time'],dt_s=h,traction_pa=params['peak_pa']*float(__import__('prescribed_traction').envelope(s['time'])),force_norm_n=float(np.linalg.norm(s['external_force'].sum(0))),impulse_ns=float(np.linalg.norm(impulse)),bulk_failed=int(np.sum(natural&(bt==1))),interface_failed=int(np.sum(natural&(bt==2))),new_softened=int(np.sum((s['d']>0)&~cut)),largest_detached_nodes=largest,detached_nodes=int(free.sum()),detached_components=len(free),N_ge14=int(np.sum(free>=14)),N_ge34=int(np.sum(free>=34)),N_ge114=int(np.sum(free>=114)),N_ge200um=int(np.sum(diameter>=200e-6)),N_ge300um=int(np.sum(diameter>=300e-6)),debris_1_2_mass_kg=float(free[free<=2].sum()*RHO*V),W_J=W,U_J=float(sc[0]),K_J=s['K'],D_J=float(sc[1]),Uc_J=s['Uc'],Jc_J=Jc,R_J=res,relative_R=rel,max_dd=float(np.max(dd)),max_speed_m_s=float(np.linalg.norm(s['velocity'],axis=1).max()),max_tensile_kappa=float(s['kappa'][~cut].max()),max_damage=float(s['d'][~cut].max()),classification=classification)
   if writer is None:writer=csv.DictWriter(stream,fieldnames=list(row));writer.writeheader()
   writer.writerow(row);rows.append(row)
   changed=np.flatnonzero(dd>0)
   if len(changed):changedfile.write(json.dumps(dict(step=s['step'],time_s=s['time'],original_bond_ids=(changed+1).tolist(),d=s['d'][changed].tolist()))+'\n')
   if s['step']%1000==0:stream.flush();print('STEP',case.name,s['step'],s['time'],'R',rel,'maxfree',maxfree,'wall',time.perf_counter()-start,flush=True)
   prev=s
  save(s,'final.npz')
 except BaseException as exc:
  stop.append(repr(exc));(case/'FAILURE_TRACEBACK.txt').write_text(traceback.format_exc());save(s,'failed_unaccepted.npz');s=prev;save(s,'last_accepted.npz')
 finally:stream.close();changedfile.close()
 cached=graph.topology(ref,s['d']);cs=components(ref,s,cached)
 # Native may be at an unaccepted trial. Never compare its graph with a different time layer.
 count,raw=ns.l.gather_bonds();actual=np.ctypeslib.as_array(raw,shape=(3*count,)).reshape(-1,3)[:,1:].copy();actual.sort(1);native_d=ns.field(3);wanted=b[native_d<1]+1
 ok=np.array_equal(np.sort(actual[:,0]*ns.N+actual[:,1]),np.sort(wanted[:,0]*ns.N+wanted[:,1]));assert ok
 ns.close();summary=dict(status='COMPLETED' if not stop else 'NUMERICALLY_UNRESOLVED',physical_status='MECHANICAL_REQUIREMENT_ONLY',params=params,time_s=float(s['time']),accepted_steps=len(rows),stop_reasons=stop,max_independent_errors=err,last_row=rows[-1] if rows else None,max_relative_energy_residual=max((r['relative_R'] for r in rows),default=0),history_max_free_nodes=maxfree,final_free_components=cs,topology_matches_current_native=ok,topology_check_state='last_native_candidate' if stop else 'final_accepted',events=events,secondary_fracture_events=secondary,classification_history=regimes,classification='NUMERICALLY_UNRESOLVED' if stop else regimes[-1],predictive_substeps=substeps,predictive_rejections=rejections,wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,library_sha256=sha(LIB),loaded_area_m2=flow.area,impulse_vector_ns=impulse)
 dump(case/'FRACTURE_EVENTS.json',events);dump(case/'DETACHED_COMPONENT_HISTORY.json',hist);dump(case/'SUMMARY.json',summary);print('DONE',case.name,summary['status'],summary['wall_s'],flush=True);return summary
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--request',required=True);args=p.parse_args();q=json.loads(Path(args.request).read_text());z=run(q['case'],q['params']);raise SystemExit(0 if z['status']=='COMPLETED' else 2)
