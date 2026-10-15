"""Native-event clot graphs, explicit attachment and exact-membership lifetimes.

Never edits a source case. Sizes use reference volume, not displayed glyphs.
"""
import csv,json,math,time,hashlib
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from .config import dump
from .fragment_analysis import read_dump,_events
from .failure_map import csv_out
from .failure_classifier import classify

def equivalent_diameter_um(nodes,dx):
    return (6*np.asarray(nodes)*dx**3/np.pi)**(1/3)*1e6

def minimum_nodes(diameter_um,dx):
    return math.ceil(np.pi*(diameter_um*1e-6)**3/(6*dx**3))

def graph_components(n,bonds,types,alive):
    bulk=bonds[alive&(types==1)];assert not len(bulk) or bulk.max()<n
    graph=coo_matrix((np.ones(2*len(bulk)),(np.r_[bulk[:,0],bulk[:,1]],np.r_[bulk[:,1],bulk[:,0]])),shape=(n,n)).tocsr()
    nc,labels=connected_components(graph);counts=np.bincount(labels,minlength=nc)
    interface=bonds[alive&(types==2)];touched=interface[interface<n]
    attached=np.zeros(nc,bool);attached[labels[touched]]=True
    return labels,counts,attached

def native_events(path,g):
    N=len(g['X']);b=g['bonds'];keys=b[:,0]*N+b[:,1];assert np.all(np.diff(keys)>0)
    raw=_events(path/'broken_bonds.local')
    if not len(raw):return raw,np.empty(0,int)
    ends=np.sort(raw[:,2:4].astype(np.int64)-1,axis=1);key=ends[:,0]*N+ends[:,1]
    order=np.argsort(raw[:,1],kind='stable');_,first=np.unique(key[order],return_index=True);ev=raw[order[first]];ev=ev[np.argsort(ev[:,1],kind='stable')]
    ends=np.sort(ev[:,2:4].astype(np.int64)-1,axis=1);ek=ends[:,0]*N+ends[:,1];idx=np.searchsorted(keys,ek)
    assert np.all(idx<len(keys)) and np.array_equal(keys[idx],ek)
    cached=np.load(path/'broken_bonds.npz');assert np.array_equal(idx,cached['initial_bond_indices']);assert np.array_equal(ev,cached['events'])
    return ev,idx

def shape_metrics(points,dx):
    centered=points-points.mean(0);_,axes=np.linalg.eigh(centered.T@centered);axes=axes[:,::-1]
    projected=centered@axes;span=np.ptp(projected,axis=0);extent=span+dx*np.abs(axes).sum(axis=0)
    return dict(reference_PCA_center_spans_um=(span*1e6).tolist(),reference_PCA_voxel_bbox_um=(extent*1e6).tolist(),reference_PCA_axes=axes.tolist(),shortest_axis_over_dx=float(extent.min()/dx),geometric_label='GEOMETRICALLY_THIN_CANDIDATE' if extent.min()<2*dx*(1-1e-10) else 'THREE_AXIS_EXTENT_GE_2DX')

def analyze(path,out,name,g,factor,*,shape_min_nodes=114,classify_fn=classify):
    start=time.perf_counter();path=Path(path);out=Path(out);out.mkdir(parents=True,exist_ok=False)
    r=json.loads((path/'request.json').read_text());old=json.loads((path/'result.json').read_text());n=g['clot_count'];N=len(g['X']);dx=r['dx_m'];vol=dx**3;mass=r['rho_kg_m3']*vol
    b=g['bonds'];typ=g['bond_types'];keys=b[:,0]*N+b[:,1];ev,idx=native_events(path,g);alive=np.ones(len(b),bool);cursor=0
    native=np.load(path/'accepted_trajectory.npz');rows=[];records={};previous={};prev_tracks=None;track_frames=[];label_frames=[];attachment_frames=[];nextid=1;split_edges=[];shapes={}
    thresholds={str(d):minimum_nodes(d,dx) for d in [200,300]};assert thresholds=={'200':34,'300':114}
    fields=['case_id','frame','step','time_s','track_id','membership_hash','minimum_node_id','node_count','volume_m3','mass_kg','d_eq_um','attachment','ge200','ge300','COM_x_m','COM_y_m','COM_z_m','COM_vx_m_s','COM_vy_m_s','COM_vz_m_s','COM_displacement_from_reference_m','COM_displacement_since_detachment_m','COM_KE_J','tracking_ambiguity','geometric_label']
    with (out/'components.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for k,(step,cols,a) in enumerate(read_dump(path/'trajectory.lammpstrj')):
            a=a[np.argsort(a[:,cols.index('id')])];assert np.array_equal(a[:,cols.index('id')],np.arange(1,N+1)) and np.isfinite(a).all()
            X=a[:,[cols.index(d) for d in ['x','y','z']]];V=a[:,[cols.index(d) for d in ['vx','vy','vz']]];t=step*r['dt_s']
            assert step==native['step'][k] and np.array_equal(X,native['position_m'][k]) and np.array_equal(V,native['velocity_m_s'][k])
            assert np.max(abs(X[n:]-g['X'][n:]))<1e-12
            while cursor<len(ev) and ev[cursor,1]<=step:alive[idx[cursor]]=False;cursor+=1
            labels,counts,attached=graph_components(n,b,typ,alive);nc=len(counts);assert counts.sum()==n and alive[typ==3].all()
            stable=np.full(nc,n+1,int);np.minimum.at(stable,labels,np.arange(1,n+1))
            assert np.array_equal(stable[labels],native['component_id'][k]) and np.array_equal(attached[labels],native['attached'][k])
            members=np.split(np.argsort(labels,kind='stable'),np.cumsum(counts)[:-1]);current={};track_per_node=np.empty(n,int)
            com=np.stack([np.bincount(labels,weights=X[:n,d])/counts for d in range(3)],axis=1);vc=np.stack([np.bincount(labels,weights=V[:n,d])/counts for d in range(3)],axis=1)
            for j,nodes in enumerate(members):
                key=nodes.astype('<i4').tobytes();uid=previous.get(key)
                if uid is None:
                    uid=nextid;nextid+=1;parents=[] if prev_tracks is None else np.unique(prev_tracks[nodes]).tolist();ambiguous=len(parents)>1
                    rec=dict(track_id=uid,case_id=name,membership_hash=hashlib.sha256(key).hexdigest(),node_ids=(nodes+1).tolist(),node_count=len(nodes),d_eq_um=float(equivalent_diameter_um(len(nodes),dx)),first_frame=k,last_frame=k,detached_frames=[],parent_track_ids=parents,tracking_ambiguity=ambiguous,reference_COM=g['X'][nodes].mean(0).tolist(),detachment_COM=None)
                    records[uid]=rec
                    for parent in parents:split_edges.append(dict(case_id=name,parent_track_id=parent,child_track_id=uid,first_child_time_s=t,overlap_nodes=int((prev_tracks[nodes]==parent).sum()),event='SPLIT',tracking_ambiguity=ambiguous))
                    if len(nodes)>=shape_min_nodes:shapes[uid]=shape_metrics(g['X'][nodes],dx)
                rec=records[uid];rec['last_frame']=k
                if not attached[j]:
                    rec['detached_frames'].append(k)
                    if rec['detachment_COM'] is None:rec['detachment_COM']=com[j].tolist()
                current[key]=uid;track_per_node[nodes]=uid
                row=dict(case_id=name,frame=k,step=step,time_s=t,track_id=uid,membership_hash=rec['membership_hash'],minimum_node_id=int(nodes[0]+1),node_count=len(nodes),volume_m3=len(nodes)*vol,mass_kg=len(nodes)*mass,d_eq_um=rec['d_eq_um'],attachment='ATTACHED' if attached[j] else 'DETACHED',ge200=len(nodes)>=34,ge300=len(nodes)>=114,COM_displacement_from_reference_m=float(np.linalg.norm(com[j]-rec['reference_COM'])),COM_displacement_since_detachment_m=None if rec['detachment_COM'] is None else float(np.linalg.norm(com[j]-rec['detachment_COM'])),COM_KE_J=float(.5*len(nodes)*mass*np.sum(vc[j]**2)),tracking_ambiguity=rec['tracking_ambiguity'],geometric_label=shapes[uid]['geometric_label'] if uid in shapes else f'BELOW_{shape_min_nodes}_NODE_SHAPE_SCREEN')
                for d,value in zip('xyz',com[j]):row['COM_'+d+'_m']=float(value)
                for d,value in zip('xyz',vc[j]):row['COM_v'+d+'_m_s']=float(value)
                writer.writerow(row)
            previous=current;prev_tracks=track_per_node;track_frames.append(track_per_node);label_frames.append(stable[labels]);attachment_frames.append(attached[labels])
            free=~attached;largest=int(counts[free].max()) if free.any() else 0
            row=dict(case_id=name,step=step,time_s=t,bulk_broken_fraction=float((~alive[typ==1]).mean()),interface_broken_fraction=float((~alive[typ==2]).mean()),debris_fraction=float(counts[counts<=2].sum()/n),debris_1_2node_fraction=float(counts[counts<=2].sum()/n),median_fragment_nodes=float(np.median(counts)),largest_free_fraction=float(counts.max()/n),resolved_fragment_count=int((counts>=3).sum()),largest_attached_component_fraction=float(counts[attached].max()/n) if attached.any() else 0.,total_attached_fraction=float(counts[attached].sum()/n),largest_detached_nodes=largest,largest_detached_d_eq_um=float(equivalent_diameter_um(largest,dx)),largest_detached_mass_kg=largest*mass,clot_KE_J=float(.5*mass*np.sum(V[:n]**2)))
            row['failure_mode']=classify_fn(row)
            for diameter,limit in thresholds.items():
                selected=free&(counts>=limit);row['N_detached_ge'+diameter]=int(selected.sum());row['mass_detached_ge'+diameter+'_kg']=float(counts[selected].sum()*mass);row['mass_fraction_detached_ge'+diameter]=float(counts[selected].sum()/n)
            pairs=cKDTree(X).query_pairs(1.35*dx,output_type='ndarray');pk=pairs[:,0]*N+pairs[:,1];loc=np.searchsorted(keys,pk);valid=loc<len(keys);bonded=np.zeros(len(pk),bool);ii=np.flatnonzero(valid);bonded[ii]=(keys[loc[ii]]==pk[ii])&alive[loc[ii]];pairs=pairs[(~bonded)&(pairs[:,0]<n)]
            dist=np.linalg.norm(X[pairs[:,1]]-X[pairs[:,0]],axis=1);deep=pairs[dist<.5*dx];affected=np.unique(deep.ravel());affected=affected[affected<n]
            row.update(active_contact_pairs=int((dist<factor*dx).sum()),min_nonbonded_distance_over_dx=float(dist.min()/dx) if len(dist) else None,deep_overlap_pairs=len(deep),deep_overlap_clot_fraction=len(affected)/n)
            rows.append(row)
    assert len(rows)==len(native['step']) and rows[-1]['step']==r['explicit_steps']
    top=list(read_dump(path/'surviving_bonds.local'))[-1];ends=np.sort(top[2][:,1:3].astype(int)-1,axis=1);assert np.array_equal(np.unique(ends[:,0]*N+ends[:,1]),keys[alive])
    times=np.array([a['time_s'] for a in rows]);life=[]
    for uid,rec in records.items():
        a=rec['first_frame'];z=rec['last_frame'];df=rec['detached_frames'];right=z==len(rows)-1;end=times[z] if right else times[z+1]
        item={k:rec[k] for k in ['case_id','track_id','membership_hash','node_count','d_eq_um','tracking_ambiguity']};item.update(parent_track_ids=json.dumps(rec['parent_track_ids']),first_time_s=times[a],last_seen_time_s=times[z],end_exclusive_or_censor_time_s=end,frame_persistence_count=z-a+1,observed_span_s=times[z]-times[a],left_hold_duration_s=end-times[a],right_censored=right,first_detached_time_s=times[df[0]] if df else None,last_detached_seen_s=times[df[-1]] if df else None,detached_frame_count=len(df),detached_observed_span_s=times[df[-1]]-times[df[0]] if df else 0.,detached_left_hold_duration_s=end-times[df[0]] if df else 0.,persistent_detached_ge200=len(df)>=3 and rec['node_count']>=34,persistent_detached_ge300=len(df)>=3 and rec['node_count']>=114,geometric_label=shapes[uid]['geometric_label'] if uid in shapes else f'BELOW_{shape_min_nodes}_NODE_SHAPE_SCREEN')
        life.append(item)
    csv_out(out/'timeseries.csv',rows);csv_out(out/'lifetimes.csv',life);csv_out(out/'split_relations.csv',split_edges)
    dump(out/'membership_and_shapes.json',dict(tracking='Exact particle-set membership episodes. Parent-child maximum-overlap relation retained; every split starts new child episodes, including the largest daughter. No persistent identity inferred from minimum node ID.',records=records,shapes=shapes))
    np.savez_compressed(out/'graph_reconstruction.npz',step=native['step'],track_id=track_frames,component_id=label_frames,attached=attachment_frames)
    summaries={}
    for d in thresholds:
        count=np.array([a['N_detached_ge'+d] for a in rows]);fraction=np.array([a['mass_fraction_detached_ge'+d] for a in rows]);where=np.flatnonzero(count>0);episodes=[a for a in life if a['node_count']>=thresholds[d] and a['detached_frame_count']]
        summaries[d]=dict(EVER_OBSERVED=bool(len(where)),MAX_SIMULTANEOUS_COUNT=int(count.max()),PEAK_DETACHED_MASS_FRACTION=float(fraction.max()),AT_FINAL_FRAME=int(count[-1]),FINAL_MASS_FRACTION=float(fraction[-1]),FIRST_APPEARANCE_TIME=float(times[where[0]]) if len(where) else None,LAST_SEEN_TIME=float(times[where[-1]]) if len(where) else None,FRAME_PERSISTENCE_COUNT=len(where),SAVED_GRID_PRESENCE_DURATION_S=float(np.diff(times,append=times[-1])[count>0].sum()),membership_episode_count=len(episodes),persistent_membership_episodes=sum(a['detached_frame_count']>=3 for a in episodes),longest_exact_membership_detached_span_s=max([a['detached_observed_span_s'] for a in episodes],default=0),right_censored=bool(count[-1]>0))
    for row in rows:row['failure_mode']=str(row['failure_mode'])
    cat={}
    for mode in sorted(set(a['failure_mode'] for a in rows)):
        present=np.array([a['failure_mode']==mode for a in rows]);cat[mode]=dict(first_time_s=float(times[present][0]),saved_grid_duration_s=float(np.diff(times,append=times[-1])[present].sum()),at_final=bool(present[-1]))
    deep=np.array([a['deep_overlap_clot_fraction']>=.01 for a in rows]);conv=np.convolve(deep.astype(int),np.ones(3,int),mode='valid');valid_dist=[a['min_nonbonded_distance_over_dx'] for a in rows if a['min_nonbonded_distance_over_dx'] is not None]
    summary=dict(case_id=name,source_path=str(path),contact_model_id='MODIFIED_1P0DX' if factor==1 else 'NATIVE_1P35DX',material_id=r['material_id'],interface_id='PORCINE_SKIN_PEELING_MODEL_PROXY',pattern_id=r['pattern_id'],tau_budget_Pa=r['tau_budget_Pa'],pulse_duration_s=r['duration_s'],contact_distance_factor=factor,dx_m=dx,horizon_m=r['delta_m'],time_end_s=float(times[-1]),n_steps=int(rows[-1]['step']),initial_clot_mass_kg=n*mass,thresholds=summaries,final=rows[-1],largest_ever_detached_d_eq_um=max(a['largest_detached_d_eq_um'] for a in rows),largest_final_detached_d_eq_um=rows[-1]['largest_detached_d_eq_um'],peak_clot_KE_J=max(a['clot_KE_J'] for a in rows),max_active_contact_pairs=max(a['active_contact_pairs'] for a in rows),min_nonbonded_distance_over_dx=min(valid_dist) if valid_dist else None,max_deep_overlap_clot_fraction=max(a['deep_overlap_clot_fraction'] for a in rows),sustained_severe_overlap=bool((conv>=3).any()),severe_overlap_saved_duration_s=float(np.diff(times,append=times[-1])[deep].sum()),dynamic_failure_classification=cat,tracking_ambiguities=sum(a['tracking_ambiguity'] for a in life),thin_detached_ge300_episodes=sum(a['geometric_label']=='GEOMETRICALLY_THIN_CANDIDATE' and a['detached_frame_count']>0 for a in life),reference_peak_local_traction_Pa=old['peak_local_traction'],actual_peak_local_traction_Pa=old['actual_peak_local_traction_Pa'],peak_actual_mean_traction_budget_Pa=old['peak_actual_mean_traction_budget_Pa'],wall_seconds=old['runtime']['wall_seconds'],original_postprocess_wall_s=old['postprocess_wall_s'],size_postprocess_wall_s=time.perf_counter()-start,provenance=dict(native_events_reconstructed=True,cached_events_equal=True,all_native_frames_reread=True,cached_components_and_attachment_equal=True,native_final_topology_equal=True,wall_excluded=True,source_case_unchanged=True),EXPLORATORY_ONLY=True,component_lifetimes=str(out/'lifetimes.csv'),debris_identity='NUMERICALLY_UNRESOLVED_DEBRIS')
    for typid,label in [(1,'bulk'),(2,'interface')]:
        sub=ev[typ[idx]==typid];summary['first_'+label+'_break_time_s']=float(sub[:,1].min()*r['dt_s']) if len(sub) else None
    dump(out/'summary.json',summary);print('SIZE_ANALYZED',name,summary['final']['failure_mode'],summaries,flush=True);return summary

def focused_checks():
    dx=50e-6;assert minimum_nodes(200,dx)==34 and minimum_nodes(300,dx)==114
    sizes={n:float(equivalent_diameter_um(n,dx)) for n in [1,2,3,5,34,49,113,114]}
    assert sizes[34]>=200 and sizes[49]>=200 and sizes[113]<300<=sizes[114]
    # Four clot nodes, two wall nodes: two bulk pairs and one interface only.
    bonds=np.array([[0,1],[1,4],[2,3],[4,5]]);types=np.array([1,2,1,3]);alive=np.ones(4,bool)
    lab,count,att=graph_components(4,bonds,types,alive);assert sorted(count)==[2,2] and att[lab[0]] and not att[lab[2]] and count.sum()==4
    alive[1]=False;lab,count,att=graph_components(4,bonds,types,alive);assert not att.any()
    thin=shape_metrics(np.c_[np.arange(114)*dx,np.zeros(114),np.zeros(114)],dx);assert thin['geometric_label']=='GEOMETRICALLY_THIN_CANDIDATE'
    return dict(status='PASS',diameters_um=sizes,threshold_nodes={'200':34,'300':114},interface_graph_and_wall_exclusion='PASS',thin_geometry_flag='PASS',scope='Small formula/topology unit checks only; no mechanics or stress benchmark')
