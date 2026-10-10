"""Minimal event/force/continuity checks and one-case transport continuation."""
import json, csv, shutil, sys, time
import numpy as np
from .config import dump
from .early_flow_coupling import OUT, RESULTS, CASE, OLD, PRIOR, DT, MASS, DTYPE, drag_law
from .y_clot_morphology_v3 import load_geometry
from .detached_fragment_metrics import graph_components
from .failure_map import csv_out
from .background_flow_sampler import BackgroundFlowSampler

def audit_native():
    g=load_geometry();n=g['clot_count'];tr=np.load(CASE/'accepted_trajectory.npz');ev=np.load(CASE/'broken_bonds.npz');online=np.load(CASE/'ONLINE_TOPOLOGY_EVENTS.npz')
    events=np.c_[ev['events'][:,1].astype(int),ev['initial_bond_indices']];logged=online['step_and_initial_bond_index']
    sort=lambda a:a[np.lexsort((a[:,1],a[:,0]))]
    assert np.array_equal(sort(events),sort(logged)),'Online current-step topology differs from independent native break events'
    alive=np.ones(len(g['bonds']),bool);first=np.full(n,-1,int);max_comp=0
    for step in np.unique(events[:,0]):
        alive[events[events[:,0]==step,1]]=False;l,c,a=graph_components(n,g['bonds'],g['bond_types'],alive)
        new=(~a[l])&(first<0);first[new]=step;max_comp=max(max_comp,int(c[~a].max(initial=0)))
    assert np.array_equal(first,online['detach_step'])
    ids=np.flatnonzero(first>=0);delay=online['drag_step'][ids]-first[ids];assert (delay==0).all()
    raw=np.memmap(CASE/'detached_force_audit.bin',dtype=DTYPE,mode='r');full=np.memmap(CASE/'saved_frame_force_audit.bin',dtype=DTYPE,mode='r')
    max_balance=0.;invalid=0;near=0;uncertain=0
    for start in range(0,len(raw),100000):
        a=raw[start:start+100000];assert np.all(a['step']>=a['detach_step']) and np.all(a['detach_step']>=0)
        assert not a['numerical_N'].any();valid=a['valid_CFD'].astype(bool);assert not a['drag_N'][~valid].any()
        assert np.isfinite(a['drag_N']).all();invalid+=int((~valid).sum());near+=int(a['near_wall_uncertain'].sum());uncertain+=int(a['clot_present_flow_uncertain'].sum())
    for start in range(0,len(full),n):
        a=full[start:start+n];assert len(a)==n and np.array_equal(a['atom_id'],np.arange(1,n+1))
        attached=a['detach_step']<0;assert not a['drag_N'][attached].any() and not a['numerical_N'][~attached].any()
        assert np.max(abs(a['numerical_N'][attached]+500*MASS*a['v_force_half_m_s'][attached]),initial=0)<1e-20
        for cid in np.unique(a['component_id'][~attached]):
            z=a[a['component_id']==cid]
            if z['valid_CFD'][0]:
                expected,tau,_=drag_law(np.array([len(z)]),z['v_force_half_m_s'].mean(0)[None,:],z['u_COM_m_s'][0,None,:])
                max_balance=max(max_balance,float(np.max(abs(z['drag_N'].sum(0)-expected[0]))))
    assert max_balance<1e-18
    native=np.load(CASE/'terminal_full_precision.npz');assert np.array_equal(native['position_m'],tr['position_m'][-1]) and np.array_equal(native['velocity_m_s'],tr['velocity_m_s'][-1])
    wall=json.loads((CASE/'WALL_CONTACT_DIAGNOSTICS.json').read_text());assert wall['maximum_center_penetration_m']<=2e-6 and wall['swept_crossings']==0
    old=json.loads((OLD/'result.json').read_text());new=json.loads((CASE/'result.json').read_text());comparison=[]
    for name,case,z in [('old',OLD,old),('early_flow',CASE,new)]:
        states=np.load(case/'accepted_trajectory.npz');labels=states['component_id'][-1];attached=states['attached'][-1];counts=np.unique(labels[~attached],return_counts=True)[1]
        core=max(np.unique(labels[attached],return_counts=True)[1]);events=np.load(case/'broken_bonds.npz')
        comparison.append(dict(case=name,first_bulk_break_s=z['first_bulk_break_time'],first_interface_break_s=z['first_interface_break_time'],detached_nodes=int((~attached).sum()),detached_components=len(counts),native_broken_bonds=len(events['events']),attached_nodes=int(attached.sum()),attached_mass_kg=int(attached.sum())*MASS,largest_attached_core_mass_kg=int(core)*MASS,unresolved_detached_nodes=int(counts[counts<=2].sum()),unresolved_mass_kg=int(counts[counts<=2].sum())*MASS,max_final_detached_component_nodes=int(counts.max(initial=0)),wall_contact_nodes=z['wall_contact']['unique_wall_contact_nodes'],maximum_center_penetration_m=z['wall_contact']['maximum_center_penetration_m'],wall_force_max_N=z['wall_contact']['maximum_nodal_wall_force_N'],solver_wall_seconds=z['runtime']['wall_seconds']))
    csv_out(OUT/'OLD_NEW_CASE_COMPARISON.csv',comparison)
    dump(OUT/'NATIVE_EARLY_FLOW_VALIDATION.json',dict(status='PASS',native_events_equal_online_events=True,independent_detachment_reconstruction=True,event_count=len(events['events']),detached_nodes=len(ids),all_first_drag_delays_steps=0,granularity_s=DT,force_audit_rows=len(raw),full_snapshot_rows=len(full),invalid_CFD_detached_node_steps=invalid,near_wall_uncertain_node_steps=near,clot_present_flow_uncertain_node_steps=uncertain,maximum_COM_force_distribution_error_N=max_balance,velocity_write=False,double_drag=False,max_detached_component_nodes_all_event_steps=max_comp,macro_fragment_ge200um=max_comp>=34,wall=wall,comparison=comparison))
    print('NATIVE_AUDIT_PASS',comparison,flush=True)

def handoff_and_transport():
    assert json.loads((OUT/'NATIVE_EARLY_FLOW_VALIDATION.json').read_text())['status']=='PASS'
    from . import fragment_handoff, early_flow_transport
    # All exports explicitly redirected in this short-lived process. Old assets untouched.
    fragment_handoff.OUT=OUT
    shutil.copy2(PRIOR/'CASE_CONFIG_FROZEN.json',OUT/'CASE_CONFIG_FROZEN.json')
    fragment_handoff.export()
    config={**early_flow_transport.CONFIG,'initial_flow_already_active_in_LAMMPS':True,'instantaneous_force_function':'early_flow_coupling.drag_law','output_early_grid':'10us through500us,50us through5ms; then20ms through5s','freeze_epoch_s':time.time()}
    dump(OUT/'TRANSPORT_CONFIG_FROZEN.json',config);early_flow_transport.run()

def audit_handoff():
    tr=np.load(CASE/'accepted_trajectory.npz');hand=json.loads((OUT/'FRAGMENT_HANDOFF.json').read_text())['components'];trans=np.load(OUT/'FRAGMENT_TRANSPORT_TRAJECTORIES.npz');flow=BackgroundFlowSampler();raw=np.memmap(CASE/'detached_force_audit.bin',dtype=DTYPE,mode='r');last=raw[raw['step']==4660]
    rows=[];comparison=[];continuation=list(csv.DictReader((OUT/'FRAGMENT_TRANSPORT_TRAJECTORIES.csv').open()));first={}
    for r in continuation:first.setdefault(int(r['component_id']),r)
    for j,r in enumerate(hand):
        ids=np.array(r['member_atom_ids'])-1;X=tr['position_m'][-1,ids];V=tr['velocity_m_s'][-1,ids];x=X.mean(0);v=V.mean(0);u=flow.sample([x]);assert u['valid'][0]
        F,tau,re=drag_law(np.array([len(ids)]),v[None,:],u['velocity_m_s']);rec=first[r['component_id']]
        ft=np.array([float(rec['drag_start_'+k+'_N']) for k in 'xyz']);xt=np.array([float(rec['start_'+k+'_m']) for k in 'xyz']);vt=np.array([float(rec['start_v'+k+'_m_s']) for k in 'xyz'])
        posjump=float(np.max(abs(trans['position_m'][0,j]-x)));veljump=float(np.max(abs(trans['velocity_m_s'][0,j]-v)));assert posjump==0 and veljump==0 and np.array_equal(xt,x) and np.array_equal(vt,v)
        lawjump=float(np.max(abs(ft-F[0])));assert lawjump<1e-20
        native_last=last[np.isin(last['atom_id'],ids+1)];half_force=native_last['drag_N'].sum(0);dv=v-native_last['v_force_half_m_s'].mean(0)
        # Verify native final half-kick, keeping integration time levels distinct.
        expected=native_last['v_force_half_m_s']+(native_last['force_before_flow_N']+native_last['drag_N']+native_last['numerical_N'])*(DT/(2*MASS))
        halfkick_error=float(np.max(abs(expected-V)));assert halfkick_error<1e-12
        rows.append(dict(component_id=r['component_id'],node_count=len(ids),member_atom_ids=json.dumps(r['member_atom_ids']),position_jump_m=posjump,velocity_jump_m_s=veljump,fluid_force_law_jump_N=lawjump,last_half_to_full_force_change_N=float(np.linalg.norm(F[0]-half_force)),native_halfkick_velocity_error_m_s=halfkick_error,internal_relative_speed_max_m_s=float(np.linalg.norm(V-v,axis=1).max()),internal_modes_after_handoff='TRUNCATED: rigid component COM model; not pointwise velocity continuity',native_COM_speed_m_s=float(np.linalg.norm(v)),CFD_speed_m_s=float(np.linalg.norm(u['velocity_m_s'][0])),slip_m_s=float(np.linalg.norm(v-u['velocity_m_s'][0])),tau_s=float(tau[0])))
        xe=np.array([float(rec[k+'_m']) for k in 'xyz']);ve=np.array([float(rec['v'+k+'_m_s']) for k in 'xyz']);ue=flow.sample([xe]);assert ue['valid'][0]
        fe=drag_law(np.array([len(ids)]),ve[None,:],ue['velocity_m_s'])[0][0]
        for side,xx,vv,force in [('LAMMPS_FULL_STEP_FORCE_REEVALUATED',x,v,F[0]),('TRANSPORT_INITIAL',xt,vt,ft),('TRANSPORT_FIRST_ACCEPTED',xe,ve,fe)]:
            row=dict(component_id=r['component_id'],phase=side,time_s=.00233 if side!='TRANSPORT_FIRST_ACCEPTED' else float(rec['absolute_time_s']))
            for prefix,array in [('x_m',xx),('v_m_s',vv),('Fdrag_N',force)]:
                for ax,val in zip('xyz',array):row[prefix+'_'+ax]=float(val) if np.isfinite(val) else ''
            comparison.append(row)
    csv_out(OUT/'HANDOFF_POSITION_VELOCITY_FORCE_AUDIT.csv',rows);csv_out(OUT/'PHASE_BOUNDARY_STATES.csv',comparison)
    transport=json.loads((OUT/'FRAGMENT_TRANSPORT_EVENTS.json').read_text());assert transport['no_accepted_center_wall_penetration']
    dump(OUT/'HANDOFF_CONTINUITY_AUDIT.json',dict(status='PASS',scope='Component COM position/velocity and identical instantaneous drag law at accepted full-step state',max_position_jump_m=max(r['position_jump_m'] for r in rows),max_velocity_jump_m_s=max(r['velocity_jump_m_s'] for r in rows),max_force_law_jump_N=max(r['fluid_force_law_jump_N'] for r in rows),max_halfstep_to_fullstep_force_change_N=max(r['last_half_to_full_force_change_N'] for r in rows),max_internal_relative_speed_truncated_m_s=max(r['internal_relative_speed_max_m_s'] for r in rows),native_and_transport_integrators_differ=True,wall_model_switch='Native20um penalty -> transport5um unilateral rigid-member constraint; retained existing model, not contact-law continuity',transport_wall_seconds=transport['wall_seconds'],transport_no_illegal_center_penetration=True,component_count=len(rows)))
    print('HANDOFF_AUDIT_PASS',len(rows),flush=True)

def audit_integrator():
    """Verify no hidden x/v replacement at onset or any recorded detached step."""
    raw=np.memmap(CASE/'detached_force_audit.bin',dtype=DTYPE,mode='r');full=np.memmap(CASE/'saved_frame_force_audit.bin',dtype=DTYPE,mode='r');n=1927
    order=np.lexsort((raw['step'],raw['atom_id']));a=raw[order]
    same=(a['atom_id'][1:]==a['atom_id'][:-1])&(a['step'][1:]==a['step'][:-1]+1)
    prev=a[:-1][same];next_=a[1:][same]
    force=prev['force_before_flow_N']+prev['drag_N']+prev['numerical_N']
    expected_v=prev['v_force_half_m_s']+DT/MASS*force
    expected_x=prev['x_m']+DT*expected_v
    v_err=float(np.max(abs(expected_v-next_['v_force_half_m_s'])));x_err=float(np.max(abs(expected_x-next_['x_m'])))
    assert v_err<1e-12 and x_err<1e-15
    snapshots={int(full[i]['step']):i for i in range(0,len(full),n)}
    online=np.load(CASE/'ONLINE_TOPOLOGY_EVENTS.npz');rows=[]
    for atom in np.flatnonzero(online['detach_step']>=0):
        step=int(online['detach_step'][atom]);current=full[snapshots[step]+atom]
        if step-1 not in snapshots:
            rows.append(dict(atom_id=int(atom+1),detach_step=step,previous_step_saved=False,velocity_error_m_s=None,position_error_m=None));continue
        previous=full[snapshots[step-1]+atom];F=previous['force_before_flow_N']+previous['drag_N']+previous['numerical_N'];v=previous['v_force_half_m_s']+DT/MASS*F;x=previous['x_m']+DT*v
        ve=float(np.max(abs(v-current['v_force_half_m_s'])));xe=float(np.max(abs(x-current['x_m'])));assert ve<1e-12 and xe<1e-15
        rows.append(dict(atom_id=int(atom+1),detach_step=step,previous_step_saved=True,velocity_error_m_s=ve,position_error_m=xe))
    csv_out(OUT/'DETACHMENT_ONSET_INTEGRATOR_CHECK.csv',rows)
    dump(OUT/'NO_VELOCITY_RESET_INTEGRATOR_AUDIT.json',dict(status='PASS',consecutive_detached_node_steps=int(same.sum()),maximum_velocity_recurrence_error_m_s=v_err,maximum_position_recurrence_error_m=x_err,onsets_with_exact_previous_step_saved=sum(r['previous_step_saved'] for r in rows),total_detachments=len(rows),missing_previous_step_scope='If any, only sparse pre-detachment full snapshots were unavailable. Every detached transition thereafter and the read-only callback contract are checked.',formula='v_half_next = v_half_previous + dt*F_previous/m; x_next = x_previous + dt*v_half_next. Actual native forces include mechanics, load, wall, numeric and fluid.',no_new_solve=True))
    print('INTEGRATOR_AUDIT_PASS',v_err,x_err,flush=True)

if __name__=='__main__':{'native':audit_native,'transport':handoff_and_transport,'handoff':audit_handoff,'integrator':audit_integrator}[sys.argv[1]]()
