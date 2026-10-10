"""Fragments are connected components of the actual surviving bond graph."""
from pathlib import Path
import csv,json
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from .config import dump

def read_dump(path):
    with Path(path).open() as f:
        while True:
            line=f.readline()
            if not line:break
            if line.strip()!='ITEM: TIMESTEP':raise ValueError(f'Malformed dump header: {line}')
            step=int(f.readline());f.readline();count=int(f.readline());line=f.readline()
            if line.startswith('ITEM: BOX BOUNDS'):
                for _ in range(3):f.readline()
                line=f.readline()
            columns=line.split()[2:]
            values=np.fromstring(''.join(f.readline() for _ in range(count)),sep=' ').reshape(count,len(columns))
            yield step,columns,values

def _events(path):
    blocks=[a for _,_,a in read_dump(path) if len(a)]
    return np.concatenate(blocks) if blocks else np.empty((0,10))

def analyze_fragments(case):
    p=Path(case['path']);r=case['request'];g=case['geometry'];N=len(g['X']);bonds=g['bonds'];M=len(bonds)
    initial_keys=bonds[:,0]*N+bonds[:,1];events=_events(p/'broken_bonds.local')
    if len(events):
        endpoints=np.sort(events[:,2:4].astype(np.int64)-1,axis=1);keys=endpoints[:,0]*N+endpoints[:,1]
        order=np.argsort(events[:,1],kind='stable');_,first=np.unique(keys[order],return_index=True);events=events[order[first]];events=events[np.argsort(events[:,1],kind='stable')]
        endpoints=np.sort(events[:,2:4].astype(np.int64)-1,axis=1);keys=endpoints[:,0]*N+endpoints[:,1];event_indices=np.searchsorted(initial_keys,keys)
        if np.any(event_indices>=M) or not np.array_equal(initial_keys[event_indices],keys):raise ValueError('Break event is not an initial bond')
    else:event_indices=np.array([],int)
    steps=[];positions=[];velocities=[];damage=[];labels_all=[];rows=[];fragments=[];alive=np.ones(M,bool);cursor=0
    initial_length=np.linalg.norm(g['X'][bonds[:,1]]-g['X'][bonds[:,0]],axis=1);maximum_sampled_stretch=0.;min_contact=np.inf;contact_samples=[]
    with (p/'fragment_components.csv').open('w',newline='') as ff:
        fields=['step','time_s','component_id','node_count','volume_m3','mass_kg','com_x_m','com_y_m','com_z_m','vx_m_s','vy_m_s','vz_m_s','speed_m_s','contains_base'];fw=csv.DictWriter(ff,fieldnames=fields);fw.writeheader()
        for step,cols,a in read_dump(p/'trajectory.lammpstrj'):
            a=a[np.argsort(a[:,cols.index('id')])]
            if not np.array_equal(a[:,cols.index('id')],np.arange(1,N+1)):raise ValueError('Lost/duplicated atoms; mass conservation violated')
            X=a[:,[cols.index(k) for k in ['x','y','z']]];vel=a[:,[cols.index(k) for k in ['vx','vy','vz']]]
            if not np.isfinite(a).all():raise ValueError('Nonfinite production state')
            if np.max(abs(X[g['base']]-g['X'][g['base']]))>1e-12:raise ValueError('Anchored base moved')
            while cursor<len(events) and events[cursor,1]<=step:alive[event_indices[cursor]]=False;cursor+=1
            e=bonds[alive];graph=coo_matrix((np.ones(2*len(e)),(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),shape=(N,N)).tocsr()
            nc,lab=connected_components(graph);counts=np.bincount(lab,minlength=nc);volume=np.bincount(lab,weights=g['volume'],minlength=nc);mass=np.bincount(lab,weights=g['mass'],minlength=nc)
            com=np.stack([np.bincount(lab,weights=g['mass']*X[:,j],minlength=nc)/mass for j in range(3)],axis=1)
            vcom=np.stack([np.bincount(lab,weights=g['mass']*vel[:,j],minlength=nc)/mass for j in range(3)],axis=1)
            ids=np.full(nc,N+1,int);np.minimum.at(ids,lab,np.arange(1,N+1));anchored=np.bincount(lab,weights=g['base'],minlength=nc)>0
            stable=ids[lab];largest=float(volume.max()/g['volume'].sum());frag=[]
            for j in range(nc):
                z=dict(step=int(step),time_s=step*r['dt_s'],component_id=int(ids[j]),node_count=int(counts[j]),volume_m3=float(volume[j]),mass_kg=float(mass[j]),com_x_m=float(com[j,0]),com_y_m=float(com[j,1]),com_z_m=float(com[j,2]),vx_m_s=float(vcom[j,0]),vy_m_s=float(vcom[j,1]),vz_m_s=float(vcom[j,2]),speed_m_s=float(np.linalg.norm(vcom[j])),contains_base=bool(anchored[j]));fw.writerow(z);frag.append(z)
            pairs=cKDTree(X).query_pairs(1.35*r['dx_m'],output_type='ndarray')
            if len(pairs):pairs=pairs[stable[pairs[:,0]]!=stable[pairs[:,1]]]
            if len(pairs):
                rel=X[pairs[:,1]]-X[pairs[:,0]];dist=np.linalg.norm(rel,axis=1);closing=np.sum(rel*(vel[pairs[:,1]]-vel[pairs[:,0]]),axis=1)<0
                minimum=float(dist.min());min_contact=min(min_contact,minimum);contact_samples.append(dict(step=int(step),pairs=len(pairs),minimum_distance_m=minimum,approaching_pairs=int(closing.sum())))
            sampled=float(np.max(np.linalg.norm(X[e[:,1]]-X[e[:,0]],axis=1)/initial_length[alive]-1)) if len(e) else 0.
            maximum_sampled_stretch=max(maximum_sampled_stretch,sampled)
            amp=.5*(1-np.cos(np.pi*min(step*r['dt_s']/r['ramp_time_s'],1.)))
            rows.append(dict(step=int(step),time_s=step*r['dt_s'],load_fraction=float(amp),broken_bond_count=int((~alive).sum()),broken_bond_fraction=float((~alive).mean()),N_fragments=int(nc),resolved_fragments=int((counts>=r['resolved_fragment_min_nodes']).sum()),largest_fragment_fraction=largest,FI=1-largest,median_fragment_volume_m3=float(np.median(volume)),sampled_max_surviving_bond_stretch=sampled,max_speed_m_s=float(np.linalg.norm(vel,axis=1).max())))
            steps.append(step);positions.append(X);velocities.append(vel);damage.append(a[:,cols.index('d_damage')]);labels_all.append(stable);fragments=frag
    trajectory=dict(step=np.array(steps),position_m=np.array(positions),velocity_m_s=np.array(velocities),damage=np.array(damage),component_id=np.array(labels_all),base=g['base'],surface=g['surface'],volume_m3=g['volume'],mass_kg=g['mass'])
    np.savez_compressed(p/'accepted_trajectory.npz',**trajectory)
    np.savez_compressed(p/'broken_bonds.npz',events=events,initial_bond_indices=event_indices,initial_bonds=bonds)
    with (p/'fragment_metrics.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    # Cross-check with the topology actually reported by LAMMPS, not damage colors.
    topo=list(read_dump(p/'surviving_bonds.local'));last=topo[-1][2];ends=np.sort(last[:,1:3].astype(int)-1,axis=1);actual_keys=np.unique(ends[:,0]*N+ends[:,1]);topology_pass=np.array_equal(actual_keys,initial_keys[alive]) and topo[-1][0]==steps[-1]
    motion=[];labels=np.array(labels_all);XX=np.array(positions)
    for frag in fragments:
        if frag['contains_base']:continue
        nodes=labels[-1]==frag['component_id'];birth=None
        for i,l in enumerate(labels):
            groups=np.unique(l[nodes]);base_labels=np.unique(l[g['base']])
            if not np.intersect1d(groups,base_labels).size:birth=i;break
        if birth is None:continue
        com0=np.average(XX[birth,nodes],axis=0,weights=g['mass'][nodes]);com1=np.average(XX[-1,nodes],axis=0,weights=g['mass'][nodes]);travel=float(np.linalg.norm(com1-com0))
        motion.append(dict(component_id=frag['component_id'],nodes=frag['node_count'],first_detached_frame=birth,first_detached_time_s=steps[birth]*r['dt_s'],post_detachment_COM_travel_m=travel,final_speed_m_s=frag['speed_m_s'],resolved=frag['node_count']>=r['resolved_fragment_min_nodes']))
    motion_pass=any(m['post_detachment_COM_travel_m']>.2*r['dx_m'] for m in motion)
    contact='NOT_EXERCISED' if not contact_samples else ('OBSERVED_AND_PASS' if min_contact>=.8*r['dx_m'] else 'FAIL')
    first=events[0] if len(events) else None;first_time=float(first[1]*r['dt_s']) if first is not None else None;final=rows[-1];perf=json.loads((p/'PERFORMANCE.json').read_text())
    result={k:r[k] for k in ['model_level','method_id','geometry_id','load_provider','pattern_id','scope']}
    result.update(material_id=case['material']['material_id'],dx=r['dx_m'],delta=r['delta_m'],node_count=N,bond_count=M,K=r['K_Pa'],G=r['G_Pa'],rho=r['rho_kg_m3'],s00=r['s00'],alpha=r['alpha'],
        load_budget=r['load_budget_N'],peak_traction=case['load']['info']['peak_local_traction_Pa'],time_to_first_break=first_time,first_break_step=int(first[1]) if first is not None else None,
        first_break_location_m=first[4:7].tolist() if first is not None else None,first_break_reference_location_m=first[7:10].tolist() if first is not None else None,
        first_break_load_fraction=float(.5*(1-np.cos(np.pi*min(first_time/r['ramp_time_s'],1)))) if first_time is not None else None,
        broken_bond_count=final['broken_bond_count'],broken_bond_fraction=final['broken_bond_fraction'],fragment_count=final['N_fragments'],resolved_fragment_count=final['resolved_fragments'],
        largest_fragment_fraction=final['largest_fragment_fraction'],FI=final['FI'],fragment_size_statistics=dict(min_volume_m3=min(f['volume_m3'] for f in fragments),median_volume_m3=float(np.median([f['volume_m3'] for f in fragments])),max_volume_m3=max(f['volume_m3'] for f in fragments),raw_components_retained=True,resolved_fragment_min_nodes=r['resolved_fragment_min_nodes']),
        final_fragments=fragments,fragment_motion=motion,first_break_ramp_fraction=first_time/r['ramp_time_s'] if first_time is not None else None,
        maximum_sampled_surviving_bond_stretch=maximum_sampled_stretch,
        runtime_seconds=perf['wall_seconds'],steps_per_second=perf['steps_per_second'],
        validation_status=dict(geometry_mass=g['info']['validation'],all_atoms_retained=True,anchored_base_fixed=True,finite_states=True,actual_surviving_topology_matches=bool(topology_pass),true_bond_break=len(events)>0,multiple_components=final['N_fragments']>=2,fragment_motion=motion_pass,contact=contact),
        contact_detail=dict(status=contact,minimum_sampled_intercomponent_distance_m=min_contact if np.isfinite(min_contact) else None,samples=contact_samples,scope='Post-fracture contacts in this production trajectory; no artificial collision case; frame-sampled sanity, not a penetration convergence study'),
        ranking_status='PENDING_THREE_PATTERN_COMPARISON',frames=len(rows),explicit_steps=steps[-1],biological_quantitative_validation='NOT_ESTABLISHED',patient_specific='NOT_ESTABLISHED')
    if not topology_pass:raise ValueError('Actual surviving topology disagrees with break-event reconstruction')
    export_result(case,result);return result

def export_result(case,result):
    dump(Path(case['path'])/'result.json',result);return result
