"""Accepted-state V3 analysis; historical V2 evidence stays read-only."""
from common import *
from current_surface import CurrentSurface, MIN_WEIGHT, MAX_WEIGHT_L1
from multifragment_shear_objective import metrics,shape,BINS
import csv,gzip,time,resource
graph=module('graph_v3_analysis',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')

def read_history(path):
    with open(path) as f:rows=list(csv.DictReader(f))
    for r in rows:
        for k,v in r.items():
            if v=='':r[k]=None
            elif v in ['True','False']:r[k]=(v=='True')
            else:
                try:r[k]=float(v)
                except ValueError:pass
    return rows

def write_csv(path,rows):
    path=Path(path)
    if not rows:path.write_text('status\nNOT_OBSERVED\n');return
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader()
        for row in rows:w.writerow({k:json.dumps(v,default=cv) if isinstance(v,(list,dict,np.ndarray)) else v for k,v in row.items()})

def l2_resolution(surf,z):
    """Partition actual L2 candidates by their first unmet discrete condition."""
    top=z['topology'];owner=surf.owner;part=np.full(surf.M,-1,np.int8)
    conditions=[
        ('LIVE_ORIGINAL_BRIDGE',top['bridge_count']>0),
        ('REFERENCE_FIT_UNRESOLVED',~top['reference_fit_valid'][owner]),
        ('MAPPING_EXTRAPOLATION_UNRESOLVED',~top['weight_ok']),
        ('CURRENT_DEFORMATION_UNRESOLVED',~z['l2_cell_valid'][owner]),
        ('AREA_UNRESOLVED',~z['geometric'][surf.L:]),
        ('L1_OVERLAP_EXCLUDED',z['l2_overlap_guard']),
        ('CLOSED_OR_OCCLUDED_PROXY',z['face_class'][surf.L:]==4),
        ('OPEN_BUT_INACCESSIBLE_PROXY',z['face_class'][surf.L:]==3),
        ('ACTIVE_ORIGINAL',z['active'][surf.L:]&(z['face_class'][surf.L:]==1)),
        ('ACTIVE_NEW',z['active'][surf.L:]&(z['face_class'][surf.L:]==5)),
        ('OTHER_UNRESOLVED',np.ones(surf.M,bool))]
    for i,(_,condition) in enumerate(conditions):part[(part<0)&top['boundary']&condition]=i
    assert np.array_equal(part>=0,top['boundary'])
    assert np.array_equal((part==8)|(part==9),z['active'][surf.L:])
    margins=[];near=np.zeros(surf.M,bool)
    for f in np.flatnonzero(part==2):
        i=owner[f];g=top['gradients'][surf.axis[f]];lo,hi=g.indptr[i:i+2]
        js=np.r_[g.indices[lo:hi],i];vs=np.r_[.5*surf.sign[f]*g.data[lo:hi],1.]
        merged={}
        for j,v in zip(js,vs):merged[int(j)]=merged.get(int(j),0.)+float(v)
        vs=np.array(list(merged.values()));minimum=float(vs.min(initial=0));norm=float(abs(vs).sum())
        assert minimum<MIN_WEIGHT or norm>MAX_WEIGHT_L1
        near[f]=minimum>=MIN_WEIGHT-1e-12 and norm<=MAX_WEIGHT_L1+1e-12
        margins.append(dict(L2_patch_index=int(f),original_owner_atom_id=int(i+1),reference_axis=int(surf.axis[f]),
            reference_sign=float(surf.sign[f]),minimum_proposed_weight=minimum,sum_absolute_proposed_weights=norm,
            inside_bounds_expanded_by_1e_12=bool(near[f]),actual_activation=False))
    return part,[name for name,_ in conditions],margins,near

def main():
    begin=time.perf_counter();cpu=time.process_time();case=O/'cases/biphasic'
    summary=json.loads((case/'SUMMARY.json').read_text());reg=json.loads((case/'COMPONENT_REGISTRY.json').read_text());rows=read_history(case/'history.csv')
    freeze=json.loads((O/'SOURCE_FREEZE.json').read_text())
    for p,h in freeze['implementation_sha256'].items():assert sha(R/p)==h,p
    for p,h in freeze['input_sha256'].items():assert sha(O/p)==h,p
    ref=load_npz(O/'provenance/reference.npz');n=int(ref['clot_count']);N=len(ref['X']);b=ref['bonds'];bt=ref['bond_types'];phase=load_npz(O/'BIPHASIC_FIELD.npz')['material_phase']
    initial=load_npz(case/'keys/initial.npz');endpath=case/'keys'/('final.npz' if summary['status']=='COMPLETED' else 'last_accepted.npz');end=load_npz(endpath)
    old_initial=load_npz(V2_O/'cases/biphasic/keys/initial.npz')
    for k in ['x','velocity','kappa','d','atom_ids','material_force']:assert np.array_equal(initial[k],old_initial[k]),k
    assert np.max(abs(initial['external_force']-old_initial['external_force']))<1e-18
    lab,sizes,attached=graph.topology(ref,end['d']);free=sizes[~attached];objective=metrics(free,n)
    step_times={0:0.,**{int(r['step']):r['time_s'] for r in rows}}
    death={int(e['parent_component_id']):e['time_s'] for e in summary['secondary_events']}
    parent={int(c):int(e['parent_component_id']) for e in summary['secondary_events'] for c in e['children_component_ids']}
    first={};last={};off={};final={};targets=[];fluid_rows=[];records=0;max_component_power_sum_error=0.
    with gzip.open(case/'CONNECTED_COMPONENT_HISTORY.jsonl.gz','rt') as f:
        for line in f:
            h=json.loads(line);records+=1;current={};members=[]
            assert int(h['step']) in step_times and abs(h['time_s']-step_times[int(h['step'])])<1e-14
            for c in h['components']:
                cid=int(c['component_id']);nn=reg[str(cid)]['nodes'];current[cid]=c;first.setdefault(cid,c);last[cid]=c;members.extend(reg[str(cid)]['atom_ids'])
                assert c['internal_alive_bonds']>=nn-1 and c['remaining_live_external_connections']==0
                assert np.isfinite(c['fluid_power_W'])
                if h['time_s']>=OFF-1e-14:assert c['fluid_power_W']<=1e-17,'Passive drag injected net power into a disconnected component'
                if c['fluid_effective_area_m2']==0:assert c['fluid_loading_status']=='FLUID_LOADING_UNRESOLVED'
                if 20<=nn<=200:
                    fluid_rows.append(dict(step=h['step'],time_s=h['time_s'],component_id=cid,nodes=nn,**{k:v for k,v in c.items() if k.startswith(('fluid_','member_integrated_','work_since_','surface_method'))}))
            assert len(members)==len(set(members))==h['detached_nodes'] and h['attached_nodes']+h['detached_nodes']==n
            if abs(h['time_s']-OFF)<1e-14:off=current
            final=current
    assert records==len(rows)+1
    for cid,c in last.items():
        q=reg[str(cid)];nn=q['nodes']
        if not 20<=nn<=200:continue
        targets.append(dict(component_id=cid,nodes=nn,atom_ids=q['atom_ids'],birth_s=q['birth_time_s'],end_s=death.get(cid,float(end['time'])),
            lifetime_s=death.get(cid,float(end['time']))-q['birth_time_s'],right_censored=cid not in death,birth_parent_component_id=parent.get(cid),
            survives_exact_set_to_end=cid in final,at_drive_off=cid in off,first_fluid=first[cid],last_fluid=c,reference_shape=shape(ref['X'][np.array(q['atom_ids'])-1])))
    descendants=[]
    final_sets={cid:set(reg[str(cid)]['atom_ids']) for cid in final}
    for cid,c in off.items():
        q=reg[str(cid)];nn=q['nodes']
        if not 20<=nn<=200:continue
        ids=set(q['atom_ids']);children=[];seen=set()
        for child,member in final_sets.items():
            if not ids.intersection(member):continue
            assert member<=ids and not seen.intersection(member);seen.update(member)
            children.append(dict(component_id=child,nodes=len(member),member_retention_fraction=len(member)/nn,target_20_200=20<=len(member)<=200,last_fluid=final[child]))
        assert seen==ids
        descendants.append(dict(component_id_at_drive_off=cid,nodes_at_drive_off=nn,exact_set_survives=cid in final,original_atom_ids=q['atom_ids'],
            max_member_retention=max(z['member_retention_fraction'] for z in children),target_descendants=sum(z['target_20_200'] for z in children),
            drive_off_state=c,descendants=sorted(children,key=lambda z:-z['nodes'])))
    splits=[]
    for event in summary['secondary_events']:
        children=list(map(int,event['children_component_ids']));pid=int(event['parent_component_id'])
        if sum(reg[str(c)]['nodes']>=20 for c in children)<2:continue
        original=set(reg[str(pid)]['atom_ids']);seen=set()
        for child in children:
            mem=set(reg[str(child)]['atom_ids']);assert mem<=original and not seen.intersection(mem);seen|=mem
        assert seen==original
        splits.append(dict(**event,parent_nodes=reg[str(pid)]['nodes'],parent_last_accepted_observation=last[pid],
            children=[dict(component_id=c,nodes=reg[str(c)]['nodes'],first_accepted_observation=first[c]) for c in children],
            inference='Temporal load/damage evidence; no isolated causal claim that drag produced this split'))
    morphology=[];target_surface=[];resolution=[];surf=CurrentSurface(ref,'uniform');z=surf.evaluate(end['x'],end['velocity'],end['d'],float(end['time']))
    resolution_code,resolution_names,weight_margins,near_limit=l2_resolution(surf,z)
    lookup={tuple(q['atom_ids']):int(cid) for cid,q in reg.items()}
    component_force_error=component_moment_error=component_power_error=0.
    for label in np.flatnonzero(~attached):
        ids=np.flatnonzero(lab==label);cid=lookup[tuple((ids+1).tolist())];member=np.zeros(N,bool);member[ids]=True
        internal=member[b].all(1)&(end['d']<1);external=(member[b].sum(1)==1)&(end['d']<1)
        assert not external.any()
        mask=member[surf.tri[:,0]]&z['active'];sites=np.unique(np.r_[surf.base.tri[mask[:surf.L]].ravel(),surf.owner[mask[surf.L:]]])
        com=end['x'][ids].mean(0);ff=end['surface_force'][mask];nf=end['external_force'][ids]
        mapping_error=dict(force_N=float(np.max(abs(ff.sum(0)-nf.sum(0)))),
            moment_about_COM_Nm=float(np.max(abs(np.cross(end['surface_centroid'][mask]-com,ff).sum(0)-np.cross(end['x'][ids]-com,nf).sum(0)))),
            power_W=float(abs(np.sum(end['surface_velocity'][mask]*ff)-np.sum(end['velocity'][ids]*nf))))
        component_force_error=max(component_force_error,mapping_error['force_N']);component_moment_error=max(component_moment_error,mapping_error['moment_about_COM_Nm']);component_power_error=max(component_power_error,mapping_error['power_W'])
        assert mapping_error['force_N']<1e-13 and mapping_error['moment_about_COM_Nm']<1e-16 and mapping_error['power_W']<1e-13
        morphology.append(dict(component_id=cid,nodes=len(ids),target_20_200=20<=len(ids)<=200,original_atom_ids=ids+1,
            reference_bbox_min_m=ref['X'][ids].min(0),reference_bbox_max_m=ref['X'][ids].max(0),current_bbox_min_m=end['x'][ids].min(0),current_bbox_max_m=end['x'][ids].max(0),
            COM_m=end['x'][ids].mean(0),COM_velocity_m_s=end['velocity'][ids].mean(0),internal_live_bulk_bonds=int((internal&(bt==1)).sum()),internal_live_interface_bonds=int((internal&(bt==2)).sum()),
            strong_nodes=int((phase[ids]==1).sum()),quadrature_owner_or_vertex_fraction=len(sites)/len(ids),internal_edge_density_per_m3=float(internal.sum()/(len(ids)*V)),**shape(end['x'][ids]),**{k:v for k,v in final[cid].items() if k.startswith('fluid_') or k=='surface_method'}))
        if 20<=len(ids)<=200:target_surface.append(morphology[-1])
        dual=member[surf.owner]&(resolution_code>=0)
        resolution.append(dict(component_id=cid,nodes=len(ids),target_20_200=20<=len(ids)<=200,
            L1_active=int(mask[:surf.L].sum()),L2_active=int(mask[surf.L:].sum()),L2_candidates=int(dual.sum()),
            L2_first_condition_counts={name:int(np.sum(dual&(resolution_code==i))) for i,name in enumerate(resolution_names)},
            rejected_mapping_near_strict_limit=int(np.sum(dual&near_limit)),
            effective_loaded_area_m2=final[cid]['fluid_effective_area_m2'],loading_status=final[cid]['fluid_loading_status'],mapping_error=mapping_error))
    force_error=moment_error=power_error=0.;snapshots=[]
    for p in sorted((case/'keys').glob('*.npz')):
        if p.stem=='failed_unaccepted':continue
        s=load_npz(p);assert s['x'].dtype==s['d'].dtype==np.float64
        assert int(s['step']) in step_times and abs(float(s['time'])-step_times[int(s['step'])])<1e-14
        assert np.array_equal(s['atom_ids'],np.arange(1,N+1)) and np.array_equal(s['x'][n:],ref['X'][n:]) and not s['velocity'][n:].any()
        ff=s['surface_force'];nf=s['external_force'];force_error=max(force_error,float(np.max(abs(ff.sum(0)-nf.sum(0)))))
        moment_error=max(moment_error,float(np.max(abs(np.cross(s['surface_centroid'],ff).sum(0)-np.cross(s['x'],nf).sum(0)))))
        power_error=max(power_error,float(abs(np.sum(s['surface_velocity']*ff)-np.sum(s['velocity']*nf))))
        snapshots.append(dict(key=p.name,step=int(s['step']),time_s=float(s['time']),sha256=sha(p)))
    assert force_error<1e-13 and moment_error<1e-16 and power_error<1e-13
    damage=ref['initial_cut'].astype(float)
    with gzip.open(case/'damage_events.jsonl.gz','rt') as f:
        for line in f:
            e=json.loads(line);ids=np.asarray(e['original_bond_ids'])-1;d=np.asarray(e['d']);assert (d>=damage[ids]).all();damage[ids]=d
    assert np.array_equal(damage,end['d'])
    initialU=float(initial['scalars'][0]);initialK=float(initial['K']);initialUc=float(initial['Uc']);ledger=[]
    for r in rows:
        deltaU=r['U_J']-initialU;deltaK=r['K_J']-initialK;deltaUc=r['Uc_J']-initialUc
        error=r['W_drive_J']+r['W_drag_J']-(deltaU+deltaK+r['D_J']+deltaUc-r['Jc_J']+r['R_J'])
        assert abs(error)<1e-18 and r['relative_R']<=.01 and r['max_dd']<=.1+1e-10
        ledger.append({k:r[k] for k in ['step','time_s','W_drive_J','W_drag_J','W_J','U_J','K_J','D_J','Uc_J','Jc_J','R_J','relative_R','W_accessible_new_faces_J','W_L2_surface_J']}|dict(DeltaU_J=deltaU,DeltaK_J=deltaK,DeltaUc_J=deltaUc,identity_error_J=error))
    if off:assert all(abs(r['drive_power_W'])<1e-20 for r in rows if r['time_s']>=OFF-1e-14)
    atomwork=load_npz(case/'accepted_fluid_work.npz')['node_work_drive_drag_J']
    assert np.max(abs(atomwork.sum(0)-np.array([rows[-1]['W_drive_J'],rows[-1]['W_drag_J']])))<1e-17
    objective.update(full_protocol_completed=summary['status']=='COMPLETED',target_count_ge5_at_end=int(objective['n_fragments_20_200'])>=5,
        max_target_count=int(max(r['N20_200'] for r in rows)),exact_target_sets_surviving_drive_off=sum(q['exact_set_survives'] for q in descendants),
        target_descendants_of_drive_off_targets=sum(q['target_descendants'] for q in descendants),
        minimum_target_count_in_drive_off=min([int(r['N20_200']) for r in rows if r['time_s']>=OFF-1e-14],default=None))
    historical={m:json.loads((V2_O/f'NATIVE_{m.upper()}_RESULT.json').read_text()) for m in ['uniform','biphasic']}
    offrows=[r for r in rows if r['time_s']>=OFF-1e-14]
    energy_off=dict(start=offrows[0],end=rows[-1]) if offrows else None
    result=dict(task=TASK,run=RUN,status=summary['status'],native=summary,objective=objective,last_row=rows[-1],final_free_sizes=sorted(free.tolist(),reverse=True),
        final_target_components=target_surface,unload=energy_off,selected_final_state=str(endpath.relative_to(R)),
        scientific_labels=[HET_LABEL,LABEL,'CURRENT_SURFACE_APPROXIMATE_COUPLING','MICROSTRUCTURE_SCALE_UNCALIBRATED'],
        comparison_scope='NewH vs oldH holds geometry/material/load fixed and changes surface representation. NewH vs oldU changes both material and surface model; it is not a material-only causal comparison.')
    dump(O/'NATIVE_V3_RESULT.json',result)
    dump(O/'FINAL_SURFACE_RESOLUTION_DIAGNOSIS.json',dict(source=str(endpath.relative_to(R)),source_sha256=sha(endpath),
        time_s=float(end['time']),components=resolution,reason_order=resolution_names,
        mapping_limits=dict(minimum_weight=MIN_WEIGHT,maximum_absolute_sum=MAX_WEIGHT_L1,rejected_candidates=len(weight_margins),
            rejected_but_within_1e_12_of_limits=int(near_limit.sum()),
            diagnostic_scope='Proximity to strict floating-point bounds only; not a new admission rule or proof that any omitted patch is physically valid. No face is enabled by this postprocessing.'),
        interpretation='Disjoint first-unmet-condition counts for L2 candidates; later conditions may also fail. Counts are not physical surface areas. L1 overlap suppression is a conservative duplicate-area guard, not proof of complete L1 coverage. Proxy closed/inaccessible classes do not establish zero physical fluid drag. No load rule or native dynamics changed.'))
    write_csv(O/'MAPPING_ADMISSIBILITY_MARGINS.csv',weight_margins)
    dump(O/'FRAGMENT_LINEAGE.json',dict(target_tracks=targets,unload_targets=descendants,cohesive_secondary_splits=splits,all_secondary_event_count=len(summary['secondary_events']),
        definitions='Exact original member-set identity and surviving original material connectivity. Contact proximity and color do not define fragments. 20 and 200 are included.'))
    write_csv(O/'FRAGMENT_MORPHOLOGY.csv',morphology);write_csv(O/'FRAGMENT_FORCE_AND_POWER.csv',fluid_rows)
    write_csv(O/'FRAGMENT_SIZE_HISTORY.csv',[{k:r[k] for k in r if k in ['step','time_s','attached_nodes','detached_nodes','components','largest_nodes','N20_200','mass_fraction20_200','mass_target_kg','debris12_detached_fraction'] or k.startswith(('count_bin_','nodes_bin_','mass_fraction_bin_'))} for r in rows])
    exposurekeys=['step','time_s','active_surface_faces','active_current_area_m2','effective_loaded_area_m2','original_exposed_faces','fluid_accessible_new_faces','closed_or_unresolved_faces','active_free_faces','free_nodes_with_retained_surface','free_nodes_without_resolved_surface','L2_active_patches','L2_effective_area_m2','L2_unresolved_patches','L2_candidate_patches','shear_peak_Pa','drive_power_W','drag_power_W']
    write_csv(O/'FRAGMENT_EXPOSURE_HISTORY.csv',[{k:r[k] for k in exposurekeys} for r in rows]);write_csv(O/'ENERGY_LEDGER.csv',ledger)
    dump(O/'DATA_INVENTORY.json',dict(force_csv_scope='20–200-node target components at every accepted state; full all-size per-fragment fluid quantities are losslessly retained in native component history',
        all_fragment_force_power_area_work_source='cases/biphasic/CONNECTED_COMPONENT_HISTORY.jsonl.gz',all_original_member_ids='cases/biphasic/COMPONENT_REGISTRY.json',
        surface_index='SURFACE_QUADRATURE_INDEX.npz: L1 has actual triangular vertex IDs; L2 has original control-volume owner, reference axis/sign/offset. Method 0=L1 and 1=L2.',
        legacy_tri_atom_ids_warning='For method=1, accepted_surface_history.npz tri_atom_ids repeats the owner three times for indexing compatibility. Those entries are NOT degenerate triangles or material edges. Render L2 from its current affine dual-face corners, or use the recorded quadrature centroid/normal/area.',
        cumulative_work_definition='member_integrated_* follows original atom ancestry from t=0; work_since_exact_member_set_birth starts at that exact set birth; no redistribution of work at splitting',snapshots=snapshots))
    np.savez_compressed(O/'SURFACE_QUADRATURE_INDEX.npz',method=surf.method,L1_tri_atom_ids=surf.base.tri+1,
        L2_original_owner_atom_ids=surf.owner+1,L2_reference_axis=surf.axis,L2_reference_sign=surf.sign,L2_reference_offset_m=surf.offset,
        L1_count=np.array(surf.L),L2_count=np.array(surf.M))
    costs={k:summary[k] for k in ['wall_s','cpu_s','peak_RSS_bytes','accepted_steps','predictive_rejections','predictive_substeps','saved_state_count']}
    costs.update(formal_runs=1,new_uniform_runs=0,old_H_reference_cpu_s=historical['biphasic']['native']['cpu_s'],
        analysis_wall_s=time.perf_counter()-begin,analysis_cpu_s=time.process_time()-cpu,analysis_peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    dump(O/'COMPUTE_COST.json',costs)
    dump(O/'DELIVERY_VALIDATION.json',dict(status='PASS',history_records=records,accepted_rows=len(rows),all_saved_FP64_states=len(snapshots),
        max_resultant_force_error_N=force_error,max_moment_error_Nm=moment_error,max_power_pair_error_W=power_error,damage_event_reconstruction_exact=True,
        final_free_components_checked=len(resolution),max_final_component_force_error_N=component_force_error,
        max_final_component_moment_COM_error_Nm=component_moment_error,max_final_component_power_pair_error_W=component_power_error,
        initial_V2_material_and_geometry_exact=True,complete_original_member_closure=True,all_target_components_have_real_internal_bonds=True,member_work_matches_global_ledger=True))
    dump(O/'FINAL_STATUS.json',dict(status='COUNT_TARGET_OBSERVED_WITH_UNCALIBRATED_SURFACE_AND_MATERIAL' if summary['status']=='COMPLETED' and objective['target_count_ge5_at_end'] else 'TARGET_NOT_MET_SINGLE_FORMAL_RUN' if summary['status']=='COMPLETED' else 'NUMERICALLY_UNRESOLVED_SINGLE_FORMAL_RUN',
        formal_status=summary['status'],time_s=float(end['time']),target_count=objective['n_fragments_20_200'],max_target_count=objective['max_target_count'],
        passive_window_target_exact_survivors=objective['exact_target_sets_surviving_drive_off'],physical_surface_validation='PARTIAL_APPROXIMATE_NOT_CALIBRATED',
        formal_runs_used=1,automatic_additional_runs=False,final_delivery_complete=False))
    print('V3_ANALYSIS_COMPLETE',summary['status'],objective,flush=True)

if __name__=='__main__':main()
