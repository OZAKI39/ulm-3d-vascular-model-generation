"""Same-time old/new H comparison using fixed original material cohorts."""
from common import *
graph=module('motion_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')

def cohort(s,ids,ref,top):
    labels,sizes,attached=top;v=s['velocity'][ids];com=s['x'][ids].mean(0);velocity=v.mean(0)
    unique,counts=np.unique(labels[ids],return_counts=True)
    return dict(COM_m=com,COM_velocity_m_s=velocity,COM_speed_m_s=float(np.linalg.norm(velocity)),
        total_kinetic_J=float(.5*RHO*V*np.sum(v*v)),relative_kinetic_J=float(.5*RHO*V*np.sum((v-velocity)**2)),
        radius_gyration_m=float(np.sqrt(np.sum((s['x'][ids]-com)**2)/len(ids))),
        actual_node_fluid_force_N=s['external_force'][ids].sum(0),actual_node_fluid_power_W=float(np.sum(s['external_force'][ids]*v)),
        component_memberships=[dict(current_label=int(label),total_nodes_in_component=int(sizes[label]),original_cohort_members=int(count),anchor_reachable=bool(attached[label])) for label,count in zip(unique,counts)])

def main():
    result=json.loads((O/'NATIVE_V3_RESULT.json').read_text());ref=load_npz(O/'provenance/reference.npz');n=int(ref['clot_count'])
    oldreg=json.loads((V2_O/'cases/biphasic/COMPONENT_REGISTRY.json').read_text())
    cohorts=[dict(definition=f'V2 final target {cid}',atom_ids=oldreg[str(cid)]['atom_ids']) for cid in [3769,3809,3846]]
    cohorts += [dict(definition=f'V3 final target {c["component_id"]}',atom_ids=c['original_atom_ids']) for c in result['final_target_components']]
    stages=[]
    for key in ['initial','rise_end','peak_load','upper_drive_off','final']:
        paths=[base/'cases/biphasic/keys'/(key+'.npz') for base in [V2_O,O]]
        if not all(p.exists() for p in paths):continue
        old,new=[load_npz(p) for p in paths];assert abs(float(old['time'])-float(new['time']))<1e-14
        tops=[graph.topology(ref,s['d']) for s in [old,new]];records=[]
        for c in cohorts:
            ids=np.asarray(c['atom_ids'])-1
            records.append(dict(**c,nodes=len(ids),old_H=cohort(old,ids,ref,tops[0]),new_H=cohort(new,ids,ref,tops[1]),
                member_RMS_position_difference_m=float(np.sqrt(np.mean(np.sum((old['x'][ids]-new['x'][ids])**2,axis=1)))),
                member_RMS_velocity_difference_m_s=float(np.sqrt(np.mean(np.sum((old['velocity'][ids]-new['velocity'][ids])**2,axis=1))))))
        stages.append(dict(key=key,time_s=float(new['time']),old_source=str(paths[0].relative_to(R)),new_source=str(paths[1].relative_to(R)),
            all_material_RMS_position_difference_m=float(np.sqrt(np.mean(np.sum((old['x'][:n]-new['x'][:n])**2,axis=1)))),
            all_material_RMS_velocity_difference_m_s=float(np.sqrt(np.mean(np.sum((old['velocity'][:n]-new['velocity'][:n])**2,axis=1)))),cohorts=records))
    passive=[]
    off=next((s for s in stages if s['key']=='upper_drive_off'),None);final=next((s for s in stages if s['key']=='final'),None)
    if off and final:
        for first,last in zip(off['cohorts'],final['cohorts']):
            assert first['atom_ids']==last['atom_ids']
            values={}
            for model in ['old_H','new_H']:
                a,b=first[model],last[model]
                values[model]=dict(COM_displacement_m=np.asarray(b['COM_m'])-a['COM_m'],
                    kinetic_change_J=b['total_kinetic_J']-a['total_kinetic_J'],relative_kinetic_change_J=b['relative_kinetic_J']-a['relative_kinetic_J'],
                    drive_off_speed_m_s=a['COM_speed_m_s'],final_speed_m_s=b['COM_speed_m_s'])
            passive.append(dict(definition=first['definition'],nodes=first['nodes'],**values))
    dump(O/'SURFACE_MOTION_COMPARISON.json',dict(status='PASS',stages=stages,passive_window_fixed_cohorts=passive,
        definitions='Fixed original atom cohorts, not assumed to be the same connected fragments across versions or times. Component memberships are measured separately. Exact physical boundary times; no trajectory interpolation.',
        causal_scope='Same frozen geometry/material/load, altered surface discretization; association with global trajectory changes does not isolate drag as the cause of any particular secondary split.'))
    print('SAME_TIME_MOTION_COMPARISON',len(stages),len(cohorts),flush=True)

if __name__=='__main__':main()
