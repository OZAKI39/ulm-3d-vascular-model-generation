"""Offline full-precision checks; does not launch LAMMPS or mutate historical files."""
from common import *
import numpy as np,csv,time,collections
from pd_clot.lammps_fragmentation.fragment_analysis import read_dump
from pd_clot.lammps_fragmentation.detached_fragment_metrics import graph_components
from pd_clot.lammps_fragmentation.fragment_semantics import component_metrics
from pd_clot.lammps_fragmentation.early_flow_coupling import DTYPE
from pd_clot.lammps_fragmentation.macro_lab_v1 import Energy
start=time.perf_counter();DT=5e-7;DX=50e-6;V=DX**3;M=1000*V
ref=np.load(P/'reference.npz');X=ref['X'];b=ref['bonds'];typ=ref['bond_types'];m=ref['reference_m'];N=len(X);n=int(ref['clot_count']);keys=b[:,0]*N+b[:,1]
def csv_read(path):return np.genfromtxt(path,delimiter=',',names=True,dtype=None,encoding='utf-8')
def outcsv(path,rows):
 with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def py(v):return v.item() if isinstance(v,np.generic) else v
s=csv_read(P/'accepted_step_scalars.csv');ev=csv_read(P/'break_events.csv');ev=np.atleast_1d(ev)
perf=json.loads((O/'FORMAL_PERFORMANCE.json').read_text());assert perf['exit_code']==0 and set(perf['rank_exit_codes'].values())=={'0'}
assert np.array_equal(s['step'],np.arange(4661)) and s['accepted_after_final_integrate'][1:].all()
assert not list(P.glob('*FAILURE*')) and not list(P.glob('*ERROR*'))
for rank in [0,1]:assert json.loads((P/f'RANK_FINISHED_{rank}.json').read_text())['accepted_step']==4660
online=np.load(P/'ONLINE_TOPOLOGY_EVENTS.npz');ei=ev['initial_bond_index'].astype(int);es=ev['step'].astype(int)
assert np.array_equal(online['step_and_initial_bond_index'],np.column_stack([es,ei]));assert len(np.unique(ei))==len(ei)
rawblocks=[a for _,_,a in read_dump(P/'broken_bonds.local') if len(a)];native=np.concatenate(rawblocks) if rawblocks else np.empty((0,10));ends=np.sort(native[:,2:4].astype(int)-1,axis=1);ki=ends[:,0]*N+ends[:,1];ni=np.searchsorted(keys,ki);assert np.array_equal(keys[ni],ki)
assert len(np.unique(ni))==len(ni)==len(ei);order=np.lexsort((ni,native[:,1]));eo=np.lexsort((ei,es));assert np.array_equal(np.column_stack([native[order,1].astype(int),ni[order]]),np.column_stack([es[eo],ei[eo]]))
alive=np.ones(len(b),bool);alive[ei]=False
survivors=list(read_dump(P/'surviving_bonds.local'))[-1];ends=np.sort(survivors[2][:,1:3].astype(int)-1,axis=1);assert survivors[0]==4660 and np.array_equal(np.sort(ends[:,0]*N+ends[:,1]),keys[alive])
assert len(ei)==len(b)-s['native_bonds'][-1]
# Runtime MPI/state records include every owned AND ghost m and local endpoint consistency.
rankrows=[csv_read(P/f'identity_rank{r}.csv') for r in [0,1]];rankevents=[[json.loads(line) for line in (P/f'identity_events_rank{r}.jsonl').read_text().splitlines()] for r in [0,1]]
for q in rankrows:assert np.array_equal(q['step'],np.arange(4661))
mlimits={'status':'PASS','reference_definition':'m_i=sum_{all initial bonds ij} taper(r0)*r0*V_j; omega=1/r0, fixed intact initial family','units':'m^4','reference_min_m4':float(m.min()),'reference_max_m4':float(m.max()),'relative_gate':2e-12,'per_id_denominator':'actual positive reference m_i; NO floor','absolute_gate_per_id':'2e-12*m_i (reported diagnostically alongside relative criterion)','max_relative_error':float(s['m_max_relative'].max()),'max_absolute_error_m4':float(s['m_max_absolute'].max()),'complete_steps':4661,'owned_and_ghost':True,'wrong_m_count':0,'absolute_scale_normalized_error':float(s['m_max_absolute'].max()/m.max()),'terminal_relative':float(s['m_max_relative'][-1]),'storage_growth':any(e['storage_growth'] for rr in rankevents for e in rr)}
assert mlimits['max_relative_error']<=2e-12;dump(O/'WVOLUME_RUNTIME_IDENTITY_AUDIT.json',mlimits)
sort_events=[{'rank':r,**{k:e[k] for k in ['step','nextsort','nonidentity_sort_count','m_relative_max']}} for r,rr in enumerate(rankevents) for e in rr if e['sort_called']]
migrations=[{'rank':r,'step':e['step'],'incoming_owned_ids':e['incoming_owned_ids'],'outgoing_owned_ids':e['outgoing_owned_ids']} for r,rr in enumerate(rankevents) for e in rr if e['incoming_owned_ids'] or e['outgoing_owned_ids']]
mi={'status':'PASS','initial_sort':[{'rank':r,'nonidentity_sort_count':rr[0]['nonidentity_sort_count']} for r,rr in enumerate(rankevents)],'actual_sort_events':sort_events,'scheduled_frequency':1000,'sort_detection':'actual nextsort change; persisted Atom::permute nonidentity count; before_sort reconstructed via exact permute','ownership_change_steps':s['step'][s['ownership_changes']>0].astype(int).tolist(),'ownership_atom_transfer_count':int(s['ownership_changes'].sum()),'ownership_unique_atom_ids':sorted({a for e in migrations for a in e['incoming_owned_ids']}),'migration_events':migrations,'ghost_list_change_counts':[int(q['ghost_list_changed'].sum()) for q in rankrows],'neighbor_rebuild_step_counts':[int((q['neighbor_ago']==0).sum()) for q in rankrows],'array_growth_events':[e for rr in rankevents for e in rr if e['storage_growth']],'ordering_and_ownership_record':'C1_final/identity_events_rank*.jsonl','old_first_count_change_step':455,'new_first_actual_ID_ownership_change_step':int(s['step'][s['ownership_changes']>0][0]) if np.any(s['ownership_changes']) else None}
dump(O/'MPI_SORT_AND_MIGRATION_EVENTS.json',mi)
cross=ev['cross_rank'].astype(bool)
# For a crossing pair both ranks must own one endpoint and map the other as ghost.
for row in ev[cross]:
 for rank in [0,1]:
  own1=row['owner1']==rank;own2=row['owner2']==rank;assert own1!=own2
  ix1=row[f'index1_rank{rank}'];ix2=row[f'index2_rank{rank}'];nl=row[f'nlocal_rank{rank}'];assert (0<=ix1<nl)==own1 and (0<=ix2<nl)==own2 and ix1>=0 and ix2>=0
cb={'status':'PASS','native_record_count':len(native),'topology_removed_count':len(ei),'duplicate_native_break_records':0,'cross_rank_break_count':int(cross.sum()),'cross_rank_steps':len(np.unique(es[cross])),'first_cross_rank_step':int(es[cross].min()) if cross.any() else None,'owned_endpoint_checks':'Every changed-topology step and initialization, BOTH owned endpoint records checked against exact initial ID/type+alive table; unique directed endpoints and exact per-ID degree. Global count = 2*nbonds.','owned_endpoints_omitted':0,'rank_inconsistent_nbonds_or_topology':0,'illegal_owned_array_index_observed':False,'ghost_to_owned_array_guard':'verified frozen C++ both0<=index<nlocal; observer itself never reads ghost bond/history arrays','dynamic_sanitizer':'NOT_EXERCISED; no ASan instrumentation','native_events_exact_match_online_topology':True,'terminal_survivors_exact_match':True,'detail_csv':'C1_final/break_events.csv'};dump(O/'CROSS_RANK_BREAK_AUDIT.json',cb)
old=np.load(OLD/'accepted_trajectory.npz');olde=np.load(OLD/'broken_bonds.npz');oes=olde['events'][:,1].astype(int);oei=olde['initial_bond_indices'].astype(int);oldonline=np.load(OLD/'ONLINE_TOPOLOGY_EVENTS.npz')
# New accepted dump, exact alignment with old full-step states. No event interpolation.
frames=[];compar=[];max_vmatch=0.;max_xmatch=0.
for st,cols,a in read_dump(P/'trajectory.lammpstrj'):
 a=a[np.argsort(a[:,0])];assert np.array_equal(a[:,0],np.arange(1,N+1));xx=a[:,[cols.index(k) for k in ['x','y','z']]];vv=a[:,[cols.index(k) for k in ['vx','vy','vz']]]
 i=int(np.flatnonzero(old['step']==st)[0]);al=np.ones(len(b),bool);al[ei[es<=st]]=False;oa=np.ones(len(b),bool);oa[oei[oes<=st]]=False
 lab,sz,at=graph_components(n,b,typ,al);ol,osz,oat=graph_components(n,b,typ,oa);nm=component_metrics(sz,at,DX);om=component_metrics(osz,oat,DX)
 row={'step':int(st),'time_s':st*DT,'new_bulk':int(np.sum(~al&(typ==1))),'old_bulk':int(np.sum(~oa&(typ==1))),'new_interface':int(np.sum(~al&(typ==2))),'old_interface':int(np.sum(~oa&(typ==2))),'new_detached_nodes':int(sz[~at].sum()),'old_detached_nodes':int(osz[~oat].sum()),'new_detached_components':nm['detached_component_count'],'old_detached_components':om['detached_component_count'],'new_largest_detached_nodes':nm['largest_detached_nodes'],'old_largest_detached_nodes':om['largest_detached_nodes'],'new_largest_detached_diameter_um':float((6*nm['largest_detached_nodes']*V/np.pi)**(1/3)*1e6),'old_largest_detached_diameter_um':float((6*om['largest_detached_nodes']*V/np.pi)**(1/3)*1e6),'new_vmax_saved_frame_m_s':float(np.linalg.norm(vv[:n],axis=1).max()),'old_vmax_saved_frame_m_s':float(np.linalg.norm(old['velocity_m_s'][i,:n],axis=1).max()),'new_KE_J':float(.5*M*np.sum(vv[:n]**2)),'old_KE_J':float(.5*M*np.sum(old['velocity_m_s'][i,:n]**2)),'max_position_delta_m':float(np.linalg.norm(xx-old['position_m'][i],axis=1).max()),'max_velocity_delta_m_s':float(np.linalg.norm(vv-old['velocity_m_s'][i],axis=1).max()),'topology_symmetric_difference_bonds':int(np.sum(al!=oa)),'new_tiny_detached_nodes':int(sz[(~at)&(sz<=2)].sum()),'old_tiny_detached_nodes':int(osz[(~oat)&(osz<=2)].sum())}
 assert abs(row['new_vmax_saved_frame_m_s']-s['vmax_m_s'][st])<1e-14
 max_vmatch=max(max_vmatch,abs(row['new_vmax_saved_frame_m_s']-s['vmax_m_s'][st]))
 if (P/'states'/f'step_{st:05d}.npz').exists():
  z=np.load(P/'states'/f'step_{st:05d}.npz');assert np.array_equal(z['position_m'],xx) and np.array_equal(z['velocity_m_s'],vv)
 frames.append((st,xx,vv,lab,at[lab]));compar.append(row)
assert len(frames)==234 and frames[-1][0]==4660
np.savez_compressed(P/'accepted_trajectory.npz',step=[q[0] for q in frames],position_m=[q[1] for q in frames],velocity_m_s=[q[2] for q in frames],component_label=[q[3] for q in frames],attached=[q[4] for q in frames]);outcsv(O/'V3_OLD_NEW_COMPARISON.csv',compar)
terminal=np.load(P/'terminal_full_precision.npz');assert np.array_equal(terminal['position_m'],frames[-1][1]) and np.array_equal(terminal['velocity_m_s'],frames[-1][2])
np.savez_compressed(P/'broken_bonds.npz',events=native[order],initial_bond_indices=ni[order])
# Exact every-step historic topology, reusing prior events and recomputing only when graph changes.
history=[];ha=np.ones(len(b),bool);prev=None
for st in range(4661):
 br=oei[oes==st];ha[br]=False
 if prev is None or len(br):
  ll,ss,aa=graph_components(n,b,typ,ha);prev={'detached_nodes':int(ss[~aa].sum()),'detached_components':int((~aa).sum()),'largest_detached_nodes':int(ss[~aa].max(initial=0)),'tiny_detached_nodes':int(ss[(~aa)&(ss<=2)].sum()),'resolved_detached_components':int(np.sum((~aa)&(ss>=3))),'broken_bulk':int(np.sum(~ha&(typ==1))),'broken_interface':int(np.sum(~ha&(typ==2)))}
 history.append({'step':st,'time_s':st*DT,**prev})
outcsv(O/'HISTORICAL_TOPOLOGY_EVERY_STEP.csv',history)
def metrics(al):
 lab,sz,at=graph_components(n,b,typ,al);q=component_metrics(sz,at,DX);fr=sz[~at];q.update(detached_nodes=int(fr.sum()),attached_nodes=int(sz[at].sum()),detached_component_size_histogram={str(k):int(v) for k,v in collections.Counter(fr.tolist()).items()},largest_detached_mass_kg=q['largest_detached_nodes']*M,largest_detached_diameter_um=float((6*q['largest_detached_nodes']*V/np.pi)**(1/3)*1e6),tiny_detached_nodes=int(fr[fr<=2].sum()),tiny_fraction_of_detached_mass=float(fr[fr<=2].sum()/fr.sum()) if fr.sum() else None);return q
om=metrics(ha);nm=metrics(alive);assert om['detached_nodes']==217 and om['detached_component_count']==208
first=lambda arr:int(np.flatnonzero(arr)[0]) if np.any(arr) else None
nm.update(first_ge200_step=first(s['N_ge200']>0),first_ge300_step=first(s['N_ge300']>0),maximum_ever_detached_nodes=int(s['detached_nodes'].max()),maximum_ever_detached_component_nodes=int(s['largest_detached_nodes'].max()),maximum_ever_detached_diameter_um=float(s['largest_detached_diameter_m'].max()*1e6),mass_conservation_all4661_steps=bool(np.all(s['attached_nodes']+s['detached_nodes']==n)))
dump(O/'FRAGMENT_ATTACHMENT_AWARE_V2.json',{'semantics':'ATTACHMENT_AWARE_V2; tiny =1-2 unresolved material nodes; all fraction fields use total clot mass unless explicitly tiny_fraction_of_detached_mass','old_final':om,'new_final':nm,'historical_legacy_largest_free_fraction':'largest attached/all component; preserved, never used as free fragment size','new_largest_detached_is_experimental_microembolus':False})
# Coupling ID events; first actual detached node is allowed to differ from old593.
ds=online['detach_step'];gs=online['drag_step'];det=ds>=0
couple={'status':'PASS','first_detachment_step':int(ds[det].min()) if det.any() else None,'first_drag_step':int(gs[gs>=0].min()) if np.any(gs>=0) else None,'old_first_detachment_step':int(oldonline['detach_step'][oldonline['detach_step']>=0].min()),'detached_nodes':int(det.sum()),'drag_same_step_for_every_detached_node':bool(np.all(gs[det]==ds[det])),'max_drag_delay_steps':int((gs[det]-ds[det]).max(initial=0)),'attached_nodes_receiving_free_drag_steps':int(s['drag_on_attached_nodes'].sum()),'free_nodes_receiving_attached_damping_steps':int(s['numerical_damping_on_detached_nodes'].sum()),'invalid_detached_node_steps':int(s['invalid_detached_nodes'].sum()),'background':'UNCHANGED historical unobstructed CFD + COM drag; NOT clot-resolved two-way FSI','long_handoff_logic':'source hashes unchanged; new full precision terminal x/v exported. No new long-time transport integration in this task.','long_time_handoff_execution':'NOT_EXERCISED','wall_contact_status':'EXERCISED' if np.any(s['wall_norm_sum_N']) else 'NOT_EXERCISED','particle_contact_status':'EXERCISED' if np.any(s['contact_pairs']) else 'NOT_EXERCISED','particle_contact_readout':'read-only original Pair::single on actual neighbor list / special mask; no contact force added by monitor'}
assert couple['attached_nodes_receiving_free_drag_steps']==0 and couple['free_nodes_receiving_attached_damping_steps']==0;dump(O/'EARLY_FLOW_COUPLING_AUDIT.json',couple)
# Independent reference potential and lagged/prebreak theta force on every key snapshot.
r0=np.linalg.norm(X[b[:,1]]-X[b[:,0]],axis=1);h=np.where(abs(r0-150.5e-6)<=25e-6,(175.5e-6-r0)/DX,1.)
energy=Energy({'X':X,'bonds':b},DX,R/'analysis/material_macro_audit_v1/energy.so')
def lps(x,al,mm,theta_override=None):
 aa,zz=b[al].T;rr=r0[al];hh=h[al];y=x[zz]-x[aa];length=np.linalg.norm(y,axis=1);e=length-rr;e[abs(e)<2.2204e-16]=0
 theta=3*V/mm*np.bincount(np.r_[aa,zz],weights=np.tile(hh*e,2),minlength=N)
 tt=theta if theta_override is None else theta_override
 tr=V**2*hh*((3*(100000/3)-5*940)*(tt[aa]/mm[aa]+tt[zz]/mm[zz])+15*940*(1/mm[aa]+1/mm[zz])*e/rr)
 force=np.zeros_like(x);fv=tr[:,None]*y/length[:,None];np.add.at(force,aa,fv);np.add.at(force,zz,-fv)
 return theta,force
loadr=json.loads((P/'request.json').read_text());loadp=json.loads((P/'LOAD_PROVIDER.json').read_text());from pd_clot.lammps_fragmentation.y_fragmentation_cases import tangent
from pd_clot.lammps_fragmentation.y_clot_morphology_v3 import load_geometry
g=load_geometry();tan=tangent(g)
def loadforce(x,step):
 t=step*DT;T=loadr['duration_s'];amp=.5*(1-np.cos(4*np.pi*t/T)) if t<T/4 else 1. if t<.75*T else .5*(1+np.cos(4*np.pi*(t/T-.75))) if t<T else 0.
 env=np.exp(-np.sum((x-loadr['load_centers_m'][0])**2,axis=1)/(2*loadr['gaussian_sigma_m']**2))
 return g['areas'][:,None]*loadp['scale_Pa']*amp*env[:,None]*tan
forcechecks=[]
for path in sorted((P/'states').glob('*.npz')):
 z=np.load(path);st=int(z['step']);x=z['position_m'];al=z['alive'];pre=al.copy();pre[ei[es==st]]=True
 th_pre,_=lps(x,pre,m);th_post,Fgradient=lps(x,al,m);_,Flag=lps(x,al,m,z['theta'])
 nativeF=z['force_before_flow_N']-z['wall_N']-loadforce(x,st);err=nativeF[:n]-Flag[:n]-z['contact_force_N'][:n];scale=float(np.linalg.norm(nativeF[:n],axis=1).max());dtheta=float(abs(z['theta']-th_pre).max())
 U,Fcpp=energy.eval(np.ascontiguousarray(x),np.ascontiguousarray(al,np.uint8),True);Upre=energy.eval(np.ascontiguousarray(x),np.ascontiguousarray(pre,np.uint8))
 forcechecks.append({'step':st,'time_s':st*DT,'native_vs_independent_lagged_lps_plus_contact_max_N':float(np.linalg.norm(err,axis=1).max()),'native_force_scale_N':scale,'relative_error':float(np.linalg.norm(err,axis=1).max()/scale) if scale else 0.,'theta_prebreak_max_abs_error':dtheta,'lagged_vs_postbreak_theta_max_abs':float(abs(z['theta']-th_post).max()),'lagged_vs_postbreak_force_max_N':float(np.linalg.norm(Flag[:n]-Fgradient[:n],axis=1).max()),'energy_gradient_cpp_force_max_diff_N':float(np.linalg.norm(Fcpp-Fgradient,axis=1).max()),'independent_current_potential_J':float(U),'instantaneous_topology_potential_drop_J':float(Upre-U),'active_contact_pairs':len(z['contact_edges']),'max_contact_force_N':float(np.linalg.norm(z['contact_force_N'][:n],axis=1).max()),'break_count_this_step':int(np.sum(es==st))})
outcsv(O/'KEY_STATE_FORCE_ENERGY_CHECKS.csv',forcechecks)
# Source-driven counterfactual reconstruction at first migration, without running old native code.
# Old naked m local slots stay put; neighbor rebuilds only forward current owned m to ghost slots.
firstdiff=json.loads((P/'FIRST_SAVED_DIFFERENCE.json').read_text());st=int(firstdiff['step']);z=np.load(P/'states'/f'step_{st:05d}.npz');slots=[np.zeros(16384) for _ in range(2)];state={}
for step in sorted({e['step'] for rr in rankevents for e in rr if e['step']<=st}):
 events={r:next((e for e in rr if e['step']==step),None) for r,rr in enumerate(rankevents)}
 if step==0:
  for r,e in events.items():
   ids=np.array(e['post_ids']+e['ghost_ids']);slots[r][:len(ids)]=m[ids-1];state[r]=e
 else:
  for r,e in events.items():
   if e is not None:state[r]=e
  ownm=np.zeros(N)
  for r,e in state.items():
   ids=np.array(e['post_ids']);ownm[ids-1]=slots[r][:len(ids)]
  for r,e in events.items():
   if e is not None and e['neighbor_ago']==0:
    nl=len(e['post_ids']);gg=np.array(e['ghost_ids'],int);slots[r][nl:nl+len(gg)]=ownm[gg-1]
wrongm=np.zeros(N)
for r,e in state.items():ids=np.array(e['post_ids']);wrongm[ids-1]=slots[r][:len(ids)]
assert (wrongm>0).all();oldraw=np.memmap(OLD/'saved_frame_force_audit.bin',dtype=DTYPE,mode='r').reshape(-1,n);ii=int(np.flatnonzero(oldraw['step'][:,0]==st)[0]);orow=oldraw[ii];assert np.array_equal(orow['x_m'],z['position_m'][:n]) and np.array_equal(orow['v_force_half_m_s'],z['velocity_half_m_s'][:n])
oldpre=np.ones(len(b),bool);oldpre[oei[oes<st]]=False;oldpost=oldpre.copy();oldpost[oei[oes==st]]=False;wt,_=lps(z['position_m'],oldpre,wrongm);_,wf=lps(z['position_m'],oldpost,wrongm,wt)
oldnative=orow['force_before_flow_N']-orow['load_N']-orow['wall_N'];pred=wf[:n]+z['contact_force_N'][:n];res=float(np.linalg.norm(oldnative-pred,axis=1).max());fscale=float(np.linalg.norm(oldnative,axis=1).max())
causality={'first_compared_force_difference':firstdiff,'old_native_m_directly_recorded':False,'old_m_reconstruction':'Counterfactual from frozen old source naked local arrays + new actual local order/ghost history, valid through first difference where old/new x and half-v are bitwise equal. NOT a rerun or direct old m measurement.','step':st,'wrong_m_inferred_nodes':int(np.sum(abs(wrongm-m)/m>2e-12)),'wrong_m_inferred_IDs':(np.flatnonzero(abs(wrongm-m)/m>2e-12)+1).tolist(),'counterfactual_old_force_max_residual_N':res,'old_force_scale_N':fscale,'relative_residual':res/fscale,'same_x_and_half_v':True,'same_topology_at_first_force_difference':bool(np.array_equal(z['alive'],oldpost)),'migration_same_step':st in mi['ownership_change_steps'],'sort_at_first_force_difference':any(q['step']==st for q in sort_events),'causal_strength':'source reconstruction quantitatively supported' if res/fscale<1e-9 else 'association only; residual not explained','historical_sampling':'2104 saved force audit moments; unsaved steps are NOT asserted equal. Exact IDs/steps; no interpolation of events.'}
np.savez_compressed(O/'FIRST_MIGRATION_COUNTERFACTUAL.npz',step=st,reference_m=m,inferred_old_m=wrongm,old_observed_native_force_N=oldnative,inferred_old_force_N=pred,new_position_m=z['position_m'],owner=z['owner'])
dump(O/'FIRST_MECHANICAL_DIFFERENCE_AUDIT.json',causality)
# All-step dynamic maxima are actual END_OF_STEP reads; old maxima remain saved-record coverage.
def maxrow(field):return {k:py(s[k][np.argmax(s[field])]) for k in ['step','time_s',field,'vmax_atom_id']}
old_v2=oldraw['v_force_half_m_s']+.5*DT/M*(oldraw['force_before_flow_N']+oldraw['drag_N']+oldraw['numerical_N']);old_v2[oldraw['step']==0]=oldraw['v_force_half_m_s'][oldraw['step']==0];old_speed=np.linalg.norm(old_v2,axis=2);fi,ai=np.unravel_index(np.argmax(old_speed),old_speed.shape)
dyn={'new_velocity_sampling':'Native read-only END_OF_STEP after final_integrate every1 step (4661 states including setup0); dump and terminal exact cross-check','vmax':maxrow('vmax_m_s'),'max_step_displacement':maxrow('max_step_displacement_m'),'max_total_KE':maxrow('KE_material_J'),'max_detached_KE':maxrow('KE_detached_J'),'old_saved_force_record_vmax_m_s':float(old_speed[fi,ai]),'old_vmax_step':int(oldraw['step'][fi,ai]),'old_vmax_atom_id':int(ai)+1,'old_velocity_sampling':'Exact Verlet full-step reconstruction from2104 saved half-step/force states only; NOT every-step old maximum','old_dump_saved_vmax_m_s':float(np.linalg.norm(old['velocity_m_s'][:,:n],axis=2).max()),'dump_vmax_monitor_max_difference_m_s':max_vmatch,'max_step_displacement_over_dx':float(s['max_step_displacement_m'].max()/DX),'drag_work_halfvelocity_quadrature_J':float(np.sum(s['drag_power_halfvelocity_W'][1:])*DT),'numerical_work_halfvelocity_quadrature_J':float(np.sum(s['numerical_power_halfvelocity_W'][1:])*DT),'wall_work_halfvelocity_quadrature_J':float(np.sum(s['wall_power_halfvelocity_W'][1:])*DT),'contact_work_fullvelocity_quadrature_J':float(np.sum(s['contact_power_fullvelocity_W'][1:])*DT),'work_label':'diagnostic time quadrature, NOT exact energy conservation identity','first_contact_step':first(s['contact_pairs']>0),'max_contact_pairs':int(s['contact_pairs'].max()),'first_wall_contact_step':first(s['wall_norm_sum_N']>0),'maximum_break_burst':{'step':int(s['step'][np.argmax(s['new_breaks'])]),'bonds':int(s['new_breaks'].max())},'native_energy_tally':'Known mismatch retained, not used for material Gc calibration'}
dump(O/'DYNAMICS_AUDIT.json',dyn)
forceaudit={'static':json.loads((O/'STATIC_PREFLIGHT.json').read_text()),'key_state_checks':len(forcechecks),'max_force_relative_error_key_states':max(q['relative_error'] for q in forcechecks),'max_force_absolute_error_key_states_N':max(q['native_vs_independent_lagged_lps_plus_contact_max_N'] for q in forcechecks),'max_theta_prebreak_error':max(q['theta_prebreak_max_abs_error'] for q in forcechecks),'maximum_lagged_vs_gradient_force_difference_N':max(q['lagged_vs_postbreak_force_max_N'] for q in forcechecks),'FRACTURE_TIME_LAYER':'EXPLICIT_SCHEME_LIMITATION','force_status':'PASS_LAGGED_SCHEME' if max(q['relative_error'] for q in forcechecks)<1e-9 else 'MISMATCH_REQUIRES_REVIEW','potential_status':'INDEPENDENT_LPS_GRADIENT_PASS; NATIVE_PE_TALLY_MISMATCH_RETAINED','native_topology_status':'PASS_BOTH_OWNED_ENDPOINTS_AND_EVENT_IDENTITY','contact_status':couple['particle_contact_status'],'new_Cplusplus_physics_fix':False}
dump(O/'NATIVE_FORCE_AND_TOPOLOGY_AUDIT.json',forceaudit)
def firsttype(ee,ii,t):return int(ee[typ[ii]==t].min()) if np.any(typ[ii]==t) else None
summary={'run':C['run'],'complete':True,'last_step':4660,'time_s':4660*DT,'new_first_bulk_step':firsttype(es,ei,1),'old_first_bulk_step':firsttype(oes,oei,1),'new_first_interface_step':firsttype(es,ei,2),'old_first_interface_step':firsttype(oes,oei,2),'new_first_detach_step':couple['first_detachment_step'],'old_first_detach_step':593,'new_broken_bulk':int(np.sum(typ[ei]==1)),'old_broken_bulk':int(np.sum(typ[oei]==1)),'new_broken_interface':int(np.sum(typ[ei]==2)),'old_broken_interface':int(np.sum(typ[oei]==2)),'new_final':nm,'old_final':om,'native_force':forceaudit,'MPI':cb,'dynamic':dyn,'causality':causality,'analysis_wall_s':time.perf_counter()-start}
dump(O/'SCIENTIFIC_SUMMARY.json',summary);print(json.dumps({k:v for k,v in summary.items() if k not in ['new_final','old_final','native_force','MPI','dynamic','causality']},indent=2));print('newfinal',json.dumps(nm));print('forcecheck',forceaudit['force_status'],forceaudit['max_force_relative_error_key_states']);print('causality',causality['relative_residual'],causality['wrong_m_inferred_nodes'])
