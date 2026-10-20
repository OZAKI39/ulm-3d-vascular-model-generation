"""Native small-specimen subprocess. Reuses immutable constitutive/native code.

The existing precracked small geometry is copied byte-for-byte. The old top
displacement grip is replaced by physical shear traction on the flat top.
All material nodes are integrated; only the original bottom anchors are fixed.
"""
from common import *
from traction_library import Couette
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
   self.area=np.zeros(self.N);top=self.ref['top_ids']-1;X=self.ref['X'];weights=np.ones(len(top))
   for j in [0,1]:weights[(X[top,j]==X[top,j].min())|(X[top,j]==X[top,j].max())]*=.5
   self.area[top]=weights*DX**2
  def dynamics(self):
   self.l.commands_string(f'group clot type 1\nfix integrate clot nve\nfix coordinator all LPS/VARIATIONAL/DYNAMIC {RHO:.17g}\nfix v3forces clot external pf/callback 1 1')
   self.ext=types.SimpleNamespace(error=None,last={})
   def callback(caller,step,nlocal,tags,x,f):
    f[:]=0
    try:
     ids=tags.astype(int)-1;assert nlocal==self.N
     t=float(self.l.get_thermo('time'));force=np.zeros((self.N,3));force[:,0]=self.area*float(flow.shear(t));f[:]=force[ids]
     self.ext.last={'external_force':force,'external_time':np.array(t)}
    except BaseException:self.ext.error=traceback.format_exc();f[:]=np.nan
   self.callback=callback;self.l.set_fix_external_callback('v3forces',self.callback);self.dynamic=True;self.l.command('run 0');return self.state()
 return Native,Contact
def predict(model,ref,x,k):
 e,_=model.geometry(x);k=np.maximum(k,np.maximum(e/ref['r0'],0));on=S0[ref['bond_types']];off=SC[ref['bond_types']];d=np.zeros(len(k));d[k>=off]=1;part=(k>on)&(k<off);d[part]=off[part]*(k[part]-on[part])/(k[part]*(off[part]-on[part]));return k,d
def run(case,params):
 start=time.perf_counter();cpu=time.process_time();case=Path(case);case.mkdir(parents=True,exist_ok=True)
 if (case/'SUMMARY.json').exists():return json.loads((case/'SUMMARY.json').read_text())
 (case/'provenance').mkdir(exist_ok=True);(case/'keys').mkdir(exist_ok=True)
 for name in ['reference.npz','reference.native','body.data']:shutil.copy2(SPEC/name,case/'provenance'/name)
 flow=Couette(params['drive_speed_m_s'],params['hold_until_s']);dump(case/'REQUEST.json',dict(parameters=params,fluid=flow.provenance(),design_freeze_sha256=sha(O/'DESIGN_FREEZE.json'),native_sha256=sha(LIB),reference_hashes={p.name:sha(p) for p in (case/'provenance').iterdir()}))
 Native,Contact=native_factory(case,flow);ns=Native();ref=ns.ref;n=ns.n;b=ref['bonds'];bt=ref['bond_types'];cut=ref['initial_cut'];model=energy_model(ref);contact=Contact(ref)
 s=ns.dynamics();initial=s;prev=s;W=0.;Jc=0.;rows=[];err={};maxfree=0;events=[];substeps=0;scale=1e-24;maxproxy=0.;largest_step=0;stop=[]
 # Frozen candidate = upper half. Boundary = remaining z=0 ligament, no initial cut.
 candidate=np.flatnonzero(ref['X'][:n,2]>0);out,inside=graph.boundary(ref,candidate);bridge=out[(bt[out]==1)&~cut[out]];internal=inside[(bt[inside]==1)&~cut[inside]]
 def error(key,val):err[key]=max(err.get(key,0.),float(val))
 def save(state,name):np.savez_compressed(case/'keys'/name,**state)
 save(s,'initial.npz');stream=(case/'history.csv').open('w');writer=None;changedfile=(case/'damage_events.jsonl').open('w')
 cached=graph.topology(ref,s['d']);history_components=[]
 try:
  while prev['time']<flow.end-1e-15:
   h=min(BASE,flow.end-prev['time']);depth=0
   while True:
    xp=prev['x'].copy();xp[:n]+=h*prev['velocity'][:n]+h*h/(2*RHO*V)*prev['force'][:n]
    kk,dp=predict(model,ref,xp,prev['kappa'])
    if np.max(dp-prev['d'])<=.1:break
    assert h>BASE/8*(1+1e-12),'MINIMUM_EVENT_SUBSTEP_EXHAUSTED'
    h/=2;depth+=1
   substeps+=bool(depth)
   if h!=ns.dt:ns.l.command(f'timestep {h:.17g}');ns.dt=h
   ns.l.command('run 1 pre no post no');s=ns.state();p=prev;dx=s['x']-p['x'];dd=s['d']-p['d'];sc=s['scalars'];nf=np.flatnonzero((s['d']==1)&(p['d']<1))
   error('predicted_x_m',np.max(abs(xp-s['x'])));error('predicted_damage',np.max(abs(dp-s['d'])));error('predicted_kappa',np.max(abs(kk-s['kappa'])));error('irreversibility',max(np.max(-dd),np.max(p['kappa']-s['kappa'])))
   error('v_VV_m_s',np.max(abs(s['velocity'][:n]-p['velocity'][:n]-h/(2*RHO*V)*(p['force'][:n]+s['force'][:n]))))
   assert err['predicted_x_m']<2e-18 and err['v_VV_m_s']<1e-10 and max(err['predicted_damage'],err['predicted_kappa'],err['irreversibility'])<1e-10
   assert np.max(dd)<=.1+1e-10 and np.array_equal(s['x'][n:],ref['X'][n:]) and not s['velocity'][n:].any()
   assert s['native_bonds']==int(np.sum(s['d']<1)) and sc[3]==int(np.sum((s['d']==1)&~cut)) and not s['d'][bt==3].any()
   assert np.isfinite(s['force']).all() and np.array_equal(s['atom_ids'],np.arange(1,ns.N+1))
   W+=float(.5*np.sum((p['external_force']+s['external_force'])*dx))
   if len(nf):
    rr=np.linalg.norm(s['x'][b[nf,1]]-s['x'][b[nf,0]],axis=1);Jc+=float(.5*contact.c*np.sum(np.minimum(rr-DX,0)**2));cached=graph.topology(ref,s['d'])
   res=W-(sc[0]-initial['scalars'][0]+s['K']-initial['K']+sc[1]+s['Uc']-initial['Uc']-Jc);scale=max(scale,abs(W),abs(sc[0]-initial['scalars'][0])+s['K']+sc[1]+s['Uc']);rel=abs(res)/scale
   assert rel<=.01,('ENERGY_RESIDUAL',rel,res,scale)
   if s['step']%100==0 or np.any(dd):
    new=model.evaluate(s['x'],s['d']);old=model.evaluate(s['x'],p['d']) if np.any(dd) else new;uc,fc,pairs=contact.evaluate(s['x'],s['d'])
    error('native_force_N',np.max(abs(new['force']-s['material_force'])));error('external_port_N',np.max(abs(s['force']-s['material_force']-fc-s['external_force'])));error('native_U_J',abs(new['U']-sc[0]));error('dissipation_increment_J',abs(old['U']-new['U']-(sc[1]-p['scalars'][1])));error('contact_U_J',abs(uc-s['Uc']))
    assert max(err['native_force_N'],err['external_port_N'])<1e-12 and max(err['native_U_J'],err['dissipation_increment_J'],err['contact_U_J'])<1e-17
   labels,sizes,attached=cached;free=sizes[~attached];largest=int(free.max(initial=0));natural=(s['d']==1)&~cut
   boundary_d=float(np.mean(s['d'][bridge]));internal_d=float(np.mean(s['d'][internal]));precursor=float(np.mean(np.clip(s['kappa'][bridge]/S0[bt[bridge]],0,1)))
   proxy=boundary_d-internal_d+.01*precursor;maxproxy=max(maxproxy,proxy)
   if largest>maxfree:maxfree=largest;largest_step=s['step'];save(s,'max_free.npz')
   row=dict(step=s['step'],time_s=s['time'],dt_s=h,shear_Pa=float(flow.shear(s['time'])),bulk_failed=int(np.sum(natural&(bt==1))),interface_failed=int(np.sum(natural&(bt==2))),new_softened=int(np.sum((s['d']>0)&~cut)),largest_detached_nodes=largest,detached_nodes=int(free.sum()),detached_components=len(free),N_ge14=int(np.sum(free>=14)),N_ge34=int(np.sum(free>=34)),N_ge114=int(np.sum(free>=114)),debris_1_2_mass_kg=float(free[free<=2].sum()*RHO*V),W_J=W,U_J=float(sc[0]),K_J=s['K'],D_J=float(sc[1]),Uc_J=s['Uc'],Jc_J=Jc,R_J=res,relative_R=rel,max_dd=float(np.max(dd)),max_speed_m_s=float(np.linalg.norm(s['velocity'],axis=1).max()),boundary_damage=boundary_d,internal_damage=internal_d,boundary_bridge_count=int(np.sum(s['d'][bridge]<1)),boundary_integrity=float(np.mean(1-s['d'][bridge])),precursor=precursor,proxy=proxy,max_tensile_kappa=float(s['kappa'][~cut].max()))
   if writer is None:writer=csv.DictWriter(stream,fieldnames=list(row));writer.writeheader()
   writer.writerow(row);rows.append(row)
   changed=np.flatnonzero(dd>0)
   if len(changed):changedfile.write(json.dumps(dict(step=s['step'],time_s=s['time'],original_bond_ids=(changed+1).tolist(),d=s['d'][changed].tolist()))+'\n')
   if len(nf):
    history_components.append(dict(step=s['step'],free_components=[(np.flatnonzero(labels==j)+1).tolist() for j in np.flatnonzero(~attached)]));save(s,f"event_{s['step']:06d}.npz")
   if s['step']%1000==0:stream.flush();print('STEP',case.name,s['step'],s['time'],'R',rel,'proxy',proxy,'wall',time.perf_counter()-start,flush=True)
   prev=s
  save(s,'final.npz')
 except BaseException as exc:
  stop.append(repr(exc));(case/'FAILURE_TRACEBACK.txt').write_text(traceback.format_exc());save(s,'failed_unaccepted.npz');s=prev;save(s,'last_accepted.npz')
 finally:stream.close();changedfile.close()
 labels,sizes,attached=graph.topology(ref,s['d']);comps=[]
 for q in np.flatnonzero(~attached):
  ids=np.flatnonzero(labels==q);out,inside=graph.boundary(ref,ids);bulk=inside[bt[inside]==1];comps.append(dict(atom_ids=ids+1,nodes=len(ids),equivalent_diameter_m=(6*len(ids)*V/np.pi)**(1/3),mass_kg=len(ids)*RHO*V,internal_live_bonds=int(np.sum(s['d'][bulk]<1)),internal_reference_bonds=len(bulk),internal_alive_fraction=float(np.mean(s['d'][bulk]<1)) if len(bulk) else None))
 count,raw=ns.l.gather_bonds();actual=np.ctypeslib.as_array(raw,shape=(3*count,)).reshape(-1,3)[:,1:].copy();actual.sort(1);wanted=b[s['d']<1]+1;ok=np.array_equal(np.sort(actual[:,0]*ns.N+actual[:,1]),np.sort(wanted[:,0]*ns.N+wanted[:,1]));assert ok
 ns.close();summary=dict(status='COMPLETED' if not stop else 'NUMERICAL_FAILURE',physical_status='VALID_LOCAL_STOKES_APPROXIMATION',fidelity='LOW_FIDELITY_PROXY',params=params,time_s=float(s['time']),accepted_steps=len(rows),stop_reasons=stop,objective=maxproxy,objective_kind='Frozen history maximum of boundary mean damage - internal mean damage + 0.01 mean clipped ligament kappa/s0',max_independent_errors=err,last_row=rows[-1] if rows else None,max_relative_energy_residual=max((r['relative_R'] for r in rows),default=0),history_max_free_nodes=maxfree,history_max_free_step=largest_step,final_free_components=comps,topology_matches_native=ok,secondary_fracture_events=history_components,secondary_fracture_scope='Original-ID components stored at every new failure; none if no failure',predictive_substeps=substeps,wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,library_sha256=sha(LIB),top_loaded_area_m2=float(ns.area.sum()))
 dump(case/'SUMMARY.json',summary);print('DONE',case.name,summary['status'],summary['wall_s'],flush=True);return summary
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--request',required=True);args=p.parse_args();q=json.loads(Path(args.request).read_text());z=run(q['case'],q['params']);raise SystemExit(0 if z['status']=='COMPLETED' else 2)
