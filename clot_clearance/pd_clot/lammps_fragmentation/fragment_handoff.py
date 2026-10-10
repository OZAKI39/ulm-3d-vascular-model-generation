"""Single frozen native C1 checkpoint; exact graph identities and inherited momenta."""
import csv,json,time
import numpy as np
from .config import dump
from .y_clot_morphology_v3 import OUT,load_geometry
from .y_bifurcation_sources import YSources,digest
from .detached_fragment_metrics import graph_components,equivalent_diameter_um
from .failure_map import csv_out

def export():
    assert not (OUT/'FRAGMENT_HANDOFF.json').exists()
    g=load_geometry();n=g['clot_count'];p=OUT/'cases/C1_GAUSSIAN_A200';tr=np.load(p/'accepted_trajectory.npz');z=np.load(p/'broken_bonds.npz');ev=z['events'];bi=z['initial_bond_indices']
    cfg=json.loads((OUT/'CASE_CONFIG_FROZEN.json').read_text());assert tr['step'][-1]==cfg['handoff_step']==4660
    alive=np.ones(len(g['bonds']),bool);alive[bi]=False;lab,count,attached=graph_components(n,g['bonds'],g['bond_types'],alive)
    graph=np.load(p/'lineage/graph_reconstruction.npz');membership=json.loads((p/'lineage/membership_and_shapes.json').read_text())['records'];life={int(a['track_id']):a for a in csv.DictReader((p/'lineage/lifetimes.csv').open())}
    s=YSources();rows=[];X=tr['position_m'][-1];v=tr['velocity_m_s'][-1]
    for j in np.flatnonzero(~attached):
        ids=np.flatnonzero(lab==j);tid=int(graph['track_id'][-1,ids[0]]);rec=membership[str(tid)];assert rec['node_ids']==(ids+1).tolist()
        rows.append(dict(component_id=int(ids.min()+1),track_id=tid,topology_parent_id=rec['parent_track_ids'],member_atom_ids=(ids+1).tolist(),detached_from_wall=True,node_count=len(ids),V_m3=len(ids)*(50e-6)**3,mass_kg=len(ids)*1000*(50e-6)**3,COM_xyz_m=X[ids].mean(0),COM_velocity_mps=v[ids].mean(0),equivalent_diameter_m=float(equivalent_diameter_um(len(ids),50e-6)*1e-6),diameter_scope='REFERENCE_VOLUME_GRAPH_PROXY_NOT_CLINICAL_SIZE',category='RESOLVED_GRAPH_COMPONENT' if len(ids)>=3 else 'NUMERICALLY_UNRESOLVED_DEBRIS',actual_handoff_time_s=cfg['handoff_time_s'],handoff_step=4660,first_detached_time_s=float(life[tid]['first_detached_time_s']),exact_membership_observed_span_s=float(life[tid]['observed_span_s']),geometry_id='CURVED_MURAL_MULTILOBE_CLOT_V3',wall_geometry_sha256=s.geometry_id,flow_sha256=s.flow_id,contact_model_id='MODIFIED_1P0DX + NATIVE_Y_WALL_REPULSION_V1'))
    assert sum(r['node_count'] for r in rows)==int((~tr['attached'][-1]).sum())
    hand=dict(status='PASS' if rows else 'NO_DETACHED_FRAGMENTS',source_case=str(p),source_trajectory_sha256=digest(p/'accepted_trajectory.npz'),source_broken_bonds_sha256=digest(p/'broken_bonds.npz'),handoff_time_s=.00233,handoff_step=4660,components=rows,initial_clot_mass_kg=n*1000*(50e-6)**3,detached_mass_kg=sum(r['mass_kg'] for r in rows),attached_mass_kg=int(tr['attached'][-1].sum())*1000*(50e-6)**3,resolved_components=sum(r['node_count']>=3 for r in rows),unresolved_packets=sum(r['node_count']<=2 for r in rows),no_velocity_replacement=True,no_attached_component_transport=True)
    dump(OUT/'FRAGMENT_HANDOFF.json',hand);csv_out(OUT/'FRAGMENT_HANDOFF.csv',[{**r,'member_atom_ids':json.dumps(r['member_atom_ids'])} for r in rows])
    # Event-step size audit and exact first-detached observation for final members.
    alive[:]=True;cursor=0;audit=[];first={r['component_id']:None for r in rows};started=time.perf_counter()
    for step in np.unique(ev[:,1]).astype(int):
        end=np.searchsorted(ev[:,1],step,side='right');alive[bi[cursor:end]]=False;cursor=end;l,c,a=graph_components(n,g['bonds'],g['bond_types'],alive);free=c[~a]
        audit.append(dict(step=int(step),time_s=float(step*.5e-6),maximum_detached_nodes=int(free.max(initial=0)),detached_ge200_count=int((free>=34).sum()),detached_ge300_count=int((free>=114).sum())))
        for r in rows:
            if first[r['component_id']] is None:
                j=l[r['member_atom_ids'][0]-1]
                if not a[j] and c[j]==r['node_count'] and np.array_equal(np.flatnonzero(l==j)+1,r['member_atom_ids']):first[r['component_id']]=int(step)
    csv_out(OUT/'NATIVE_EVENT_SIZE_THRESHOLDS.csv',audit)
    dump(OUT/'NATIVE_EVENT_TRANSIENT_SIZE_AUDIT.json',dict(status='PASS',distinct_break_steps=len(audit),maximum_detached_nodes=max((r['maximum_detached_nodes'] for r in audit),default=0),ge200_ever=any(r['detached_ge200_count'] for r in audit),ge300_ever=any(r['detached_ge300_count'] for r in audit),first_exact_detached_step_for_handoff_components=first,wall_seconds=time.perf_counter()-started))
    for r in rows:
        r['first_detached_saved_frame_time_s']=r['first_detached_time_s'];r['first_exact_detached_native_step']=first[r['component_id']];r['first_exact_detached_native_time_s']=first[r['component_id']]*.5e-6
    dump(OUT/'FRAGMENT_HANDOFF.json',hand);csv_out(OUT/'FRAGMENT_HANDOFF.csv',[{**r,'member_atom_ids':json.dumps(r['member_atom_ids'])} for r in rows])
    print('HANDOFF',len(rows),hand['resolved_components'],hand['unresolved_packets'],hand['detached_mass_kg'],flush=True)

if __name__=='__main__':export()
