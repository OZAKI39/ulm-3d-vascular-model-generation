from common import *
from native import native,predict
from indenter import *
import time,resource,traceback,csv,gzip,argparse,ctypes
graph=module('crushing_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')
def save(case,name,s):
 p=case/'keys'/name;tmp=p.with_suffix('.tmp.npz');np.savez_compressed(tmp,**s);os.replace(tmp,p)
def run(name,heterogeneous):
 start=time.perf_counter();cpu=time.process_time();case=O/'cases'/name
 assert not (case/'history.csv').exists(),'Preserve interrupted attempts; use a new output name for a justified recovery'
 ns,C=native(case,heterogeneous);ref=ns.ref;n=ns.n;N=ns.N;b=ref['bonds'];bt=ref['bond_types'];cut=ref['initial_cut'];model=energy_model(ref);contact=C(ref);on=ns.field(13)
 assert json.loads((O/'KEY_CONSISTENCY_CHECKS.json').read_text())['status']=='PASS'
 dump(case/'REQUEST.json',dict(name=name,heterogeneous=heterogeneous,label='KINEMATIC_IMPACT_SURROGATE',library_sha256=sha(LIB),source_sha256={p.name:sha(p) for p in A.glob('*.py')},input_sha256={name:sha(O/name) for name in ['SPECIMEN_AND_MATERIAL_FREEZE.json','INDENTER_GEOMETRY.json','INDENTER_LOADING_PROTOCOL.json','HETEROGENEITY_FIELD.npz','HETEROGENEITY_FREEZE.json','ANALYSIS_CRITERIA_FREEZE.json']}))
 s=ns.dynamics();initial=s;prev=s;save(case,'initial.npz',s);W=0.;Jc=0.;impulse=np.zeros(3);rows=[];events={};err={};stop=[];substeps=0;rejects=0;scale=1e-24;secondary=[];maxscore=(-1,);contact_duration=0.;snap_next=50e-6
 cached=graph.topology(ref,s['d']);registry={};lookup={};current=[];previous_free={}
 history=gzip.open(case/'CONNECTED_COMPONENT_HISTORY.jsonl.gz','wt');changes=gzip.open(case/'damage_events.jsonl.gz','wt');stream=(case/'history.csv').open('w');writer=None
 def error(key,value):err[key]=max(err.get(key,0.),float(value))
 def event(name,ids,state):
  if name in events or not len(ids):return
  ends=b[ids];mid=ref['X'][ends].mean(1)
  events[name]=dict(time_s=state['time'],step=state['step'],original_bond_ids=ids+1,atom_ids=ends+1,bond_types=bt[ids],reference_midpoints_m=mid,actual_midpoints_m=state['x'][ends].mean(1),stretch=np.linalg.norm(state['x'][ends[:,1]]-state['x'][ends[:,0]],axis=1)/ref['r0'][ids]-1,kappa=state['kappa'][ids],damage=state['d'][ids],s0=on[ids],near_anchor=(mid[:,2]-ref['X'][:n,2].min())<=HORIZON,evidence='keys/'+name+'.npz');save(case,name+'.npz',state)
 def component_row(state,topology_changed):
  nonlocal current,previous_free
  lab,sizes,attached=cached;cs=[]
  if topology_changed:
   current=[]
   for label in np.flatnonzero(~attached):
    ids=np.flatnonzero(lab==label);key=tuple((ids+1).tolist())
    if key not in lookup:
     cid=len(registry)+1;lookup[key]=cid;registry[cid]=dict(atom_ids=ids+1,nodes=len(ids),reference_volume_m3=len(ids)*V,equivalent_diameter_m=(6*len(ids)*V/np.pi)**(1/3),birth_time_s=state['time'],birth_step=state['step'])
    cid=lookup[key];current.append((label,cid,ids))
   for old,oldset in previous_free.items():
    children=[cid for _,cid,ids in current if oldset.intersection((ids+1).tolist())]
    if len(children)>1:
     secondary.append(dict(time_s=state['time'],step=state['step'],parent_component_id=old,children_component_ids=children))
     if len(secondary)==1:save(case,'first_secondary.npz',state)
   previous_free={cid:set((ids+1).tolist()) for _,cid,ids in current}
  bulk=(bt==1)&~cut;bb=b[bulk];same=lab[bb[:,0]]==lab[bb[:,1]];ilab=lab[bb[same,0]];dam=state['d'][bulk][same];counts=np.bincount(ilab,minlength=len(sizes));alive=np.bincount(ilab,weights=dam<1,minlength=len(sizes));dsum=np.bincount(ilab,weights=dam,minlength=len(sizes))
  for label,cid,ids in current:
   v=state['velocity'][ids];vc=v.mean(0);relativeKE=.5*RHO*V*np.sum((v-vc)**2)
   cs.append(dict(component_id=cid,internal_initial_bonds=int(counts[label]),internal_alive_bonds=int(alive[label]),internal_survival=float(alive[label]/counts[label]) if counts[label] else None,mean_internal_damage=float(dsum[label]/counts[label]) if counts[label] else None,centroid_m=state['x'][ids].mean(0),centroid_velocity_m_s=vc,internal_relative_kinetic_J=float(relativeKE),remaining_live_external_connections=0))
  history.write(json.dumps(dict(time_s=state['time'],step=state['step'],attached_nodes=int(sizes[attached].sum()),detached_nodes=int(sizes[~attached].sum()),components=cs),default=cv)+'\n')
  return cs
 component_row(s,True)
 try:
  while prev['time']<END-1e-15:
   boundary=min([t for t in [APPROACH,PEAK,RETRACT,CLEAR,END] if t>prev['time']+1e-15]);h=min(BASE,boundary-prev['time']);depth=0
   while True:
    xp=prev['x'].copy();xp[:n]+=h*(prev['velocity'][:n]+(.5*h/(RHO*V))*prev['force'][:n]);kk,dp=predict(model,ref,xp,prev['kappa'],on)
    if np.max(dp-prev['d'])<=.1:break
    assert h>MIN_H*(1+1e-12),'MINIMUM_EVENT_SUBSTEP_EXHAUSTED'
    h/=2;depth+=1;rejects+=1
   substeps+=bool(depth)
   if h!=ns.dt:ns.l.command(f'timestep {h:.17g}');ns.dt=h
   ns.l.command('run 1 pre no post no');s=ns.state();p=prev;dd=s['d']-p['d'];sc=s['scalars'];nf=np.flatnonzero((s['d']==1)&(p['d']<1)&~cut)
   error('predicted_x_m',np.max(abs(xp-s['x'])));error('predicted_damage',np.max(abs(dp-s['d'])));error('predicted_kappa',np.max(abs(kk-s['kappa'])));error('irreversibility',max(np.max(-dd),np.max(p['kappa']-s['kappa'])));error('velocity_VV_m_s',np.max(abs(s['velocity'][:n]-p['velocity'][:n]-h/(2*RHO*V)*(p['force'][:n]+s['force'][:n]))))
   assert err['predicted_x_m']<2e-18 and err['velocity_VV_m_s']<1e-10 and max(err['predicted_damage'],err['predicted_kappa'],err['irreversibility'])<1e-10
   assert np.max(dd)<=.1+1e-10 and np.array_equal(s['x'][n:],ref['X'][n:]) and not s['velocity'][n:].any()
   assert s['native_bonds']==int(np.sum(s['d']<1)) and not s['d'][bt==3].any() and np.isfinite(s['force']).all()
   error('indenter_energy_J',abs(float(s['indenter_U']-s['indenter_U_oracle'])));error('action_reaction_N',np.max(abs(s['indenter_reaction']+s['indenter_force'].sum(0))))
   assert err['indenter_energy_J']<1e-17 and err['action_reaction_N']<1e-12
   winc=ns.tool.work(p['x'],s['x'],p['time'],s['time']);wnew=W+winc;jnew=Jc
   if len(nf):
    rr=np.linalg.norm(s['x'][b[nf,1]]-s['x'][b[nf,0]],axis=1);jnew+=float(.5*contact.c*np.sum(np.minimum(rr-DX,0)**2))
   res=wnew-(sc[0]-initial['scalars'][0]+s['K']-initial['K']+sc[1]+s['Uc']-initial['Uc']-jnew+float(s['indenter_U']-initial['indenter_U']));newscale=max(scale,abs(wnew),abs(sc[0]-initial['scalars'][0])+s['K']+sc[1]+s['Uc']+float(s['indenter_U']));rel=abs(res)/newscale
   assert rel<=.01,('ENERGY_RESIDUAL',rel,res,newscale)
   if s['step']%100==0 or np.any(dd):
    new=model.evaluate(s['x'],s['d']);old=model.evaluate(s['x'],p['d']) if np.any(dd) else new;uc,fc,pairs=contact.evaluate(s['x'],s['d']);error('native_force_N',np.max(abs(new['force']-s['material_force'])));error('force_port_N',np.max(abs(s['force']-s['material_force']-fc-s['indenter_force'])));error('native_U_J',abs(new['U']-sc[0]));error('dissipation_increment_J',abs(old['U']-new['U']-(sc[1]-p['scalars'][1])));error('pair_contact_U_J',abs(uc-s['Uc']));error('work_quadrature_4_vs_8_J',abs(winc-ns.tool.work(p['x'],s['x'],p['time'],s['time'],8)))
    assert max(err['native_force_N'],err['force_port_N'])<1e-12 and max(err['native_U_J'],err['dissipation_increment_J'],err['pair_contact_U_J'])<1e-17
   W=wnew;Jc=jnew;scale=newscale;impulse+=.5*h*(p['indenter_reaction']+s['indenter_reaction']);contact_duration+=h*bool(len(p['contact_atom_ids']) or len(s['contact_atom_ids']))
   event('first_kappa_exceeds_s0',np.flatnonzero((s['kappa']>on)&~cut),s);event('first_softening',np.flatnonzero((s['d']>0)&~cut),s);event('first_failure',nf,s)
   if len(nf):cached=graph.topology(ref,s['d'])
   lab,sizes,attached=cached;free=sizes[~attached];large=int(free.max(initial=0));medium=(free>=20)&(free<=200);debris=free[free<=2].sum();cs=component_row(s,bool(len(nf)))
   if large>=2 and 'first_multinode' not in events:
    events['first_multinode']=dict(time_s=s['time'],step=s['step'],component_ids=[q['component_id'] for q in cs],evidence='keys/first_multinode.npz');save(case,'first_multinode.npz',s)
   score=(int(medium.sum()),int(free[medium].sum()),int(np.sum(free>=3)),int(free.sum()))
   if score>maxscore:maxscore=score;save(case,'max_multifragment.npz',s)
   if len(s['contact_atom_ids']) and 'first_contact' not in events:
    events['first_contact']=dict(time_s=s['time'],step=s['step'],atom_ids=s['contact_atom_ids'],evidence='keys/first_contact.npz');save(case,'first_contact.npz',s)
   for t,key in [(PEAK,'peak_indent'),(RETRACT,'retract_start'),(CLEAR,'withdrawal_end')]:
    if abs(s['time']-t)<1e-14:save(case,key+'.npz',s)
   if s['time']>=snap_next-1e-14:save(case,f'frame_{s["step"]:06d}.npz',s);snap_next+=50e-6
   diameter=(6*free*V/np.pi)**(1/3);dfrac=float(debris/free.sum()) if free.sum() else None;instant=bool(medium.sum()>=5 and large/n<=.3 and dfrac is not None and dfrac<=.2)
   if instant and 'first_instant_target' not in events:
    events['first_instant_target']=dict(time_s=s['time'],step=s['step'],evidence='keys/first_instant_target.npz');save(case,'first_instant_target.npz',s)
   row=dict(step=s['step'],time_s=s['time'],dt_s=h,phase=phase(s['time']),indent_displacement_m=displacement(s['time']),indent_speed_m_s=speed(s['time']),reaction_x_N=float(s['indenter_reaction'][0]),reaction_y_N=float(s['indenter_reaction'][1]),reaction_z_N=float(s['indenter_reaction'][2]),reaction_norm_N=float(np.linalg.norm(s['indenter_reaction'])),contact_nodes=len(s['contact_atom_ids']),contact_duration_s=contact_duration,max_overlap_m=float(s['indenter_overlap'].max(initial=0)),impulse_x_Ns=impulse[0],impulse_y_Ns=impulse[1],impulse_z_Ns=impulse[2],bulk_failed=int(np.sum((s['d']==1)&~cut&(bt==1))),interface_failed=int(np.sum((s['d']==1)&~cut&(bt==2))),softened=int(np.sum((s['d']>0)&~cut)),attached_nodes=int(sizes[attached].sum()),detached_nodes=int(free.sum()),components=len(free),largest_nodes=large,largest_fraction=large/n,N20_200=int(medium.sum()),mass_fraction20_200=float(free[medium].sum()/n),N_ge34=int(np.sum(free>=34)),N_ge114=int(np.sum(free>=114)),N_ge200um=int(np.sum(diameter>=200e-6)),N_ge300um=int(np.sum(diameter>=300e-6)),debris12_nodes=int(debris),debris12_detached_fraction=dfrac,instant_target=instant,W_J=W,U_J=float(sc[0]),K_J=s['K'],D_J=float(sc[1]),Uc_J=s['Uc'],U_indenter_J=float(s['indenter_U']),Jc_J=Jc,numerical_damping_work_J=0.,R_J=res,relative_R=rel,max_dd=float(dd.max()),max_speed_m_s=float(np.linalg.norm(s['velocity'],axis=1).max()),max_kappa=float(s['kappa'][~cut].max()))
   if writer is None:writer=csv.DictWriter(stream,fieldnames=list(row));writer.writeheader()
   writer.writerow(row);rows.append(row);changed=np.flatnonzero(dd>0)
   if len(changed):changes.write(json.dumps(dict(step=s['step'],time_s=s['time'],original_bond_ids=changed+1,d=s['d'][changed]),default=cv)+'\n')
   if s['step']%1000==0:stream.flush();history.flush();changes.flush();print('STEP',name,s['step'],s['time'],'contact',len(s['contact_atom_ids']),'soft',row['softened'],'N20',row['N20_200'],'largest',large,'R',rel,'wall',time.perf_counter()-start,flush=True)
   prev=s
  save(case,'final.npz',s)
 except BaseException as exc:
  stop.append(repr(exc));(case/'FAILURE_TRACEBACK.txt').write_text(traceback.format_exc());save(case,'failed_unaccepted.npz',s);s=prev;save(case,'last_accepted.npz',s)
 finally:stream.close();history.close();changes.close()
 count,raw=ns.l.gather_bonds();actual=np.ctypeslib.as_array(raw,shape=(3*count,)).reshape(-1,3)[:,1:].copy();actual.sort(1);wanted=b[ns.field(3)<1]+1;native_topology=bool(np.array_equal(np.sort(actual[:,0]*N+actual[:,1]),np.sort(wanted[:,0]*N+wanted[:,1])));assert native_topology
 assert np.array_equal(ns.field(13),on);ns.close()
 summary=dict(status='COMPLETED' if not stop else 'NUMERICALLY_UNRESOLVED',name=name,heterogeneous=heterogeneous,accepted_steps=len(rows),time_s=float(s['time']),stop_reasons=stop,events=events,last_row=rows[-1] if rows else None,max_independent_errors=err,max_relative_energy_residual=max((q['relative_R'] for q in rows),default=0),max_medium_components=max((q['N20_200'] for q in rows),default=0),history_max_free_nodes=max((q['largest_nodes'] for q in rows),default=0),secondary_events=secondary,predictive_substeps=substeps,predictive_rejections=rejects,wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,library_sha256=sha(LIB),native_topology_matches_current_native=True,original_onsets_unchanged=True)
 dump(case/'COMPONENT_REGISTRY.json',registry);dump(case/'FRACTURE_EVENTS.json',events);dump(case/'SUMMARY.json',summary);print('DONE',name,summary['status'],summary['wall_s'],flush=True);return summary
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--name',required=True);p.add_argument('--heterogeneous',action='store_true');args=p.parse_args();result=run(args.name,args.heterogeneous);raise SystemExit(0 if result['status']=='COMPLETED' else 2)
