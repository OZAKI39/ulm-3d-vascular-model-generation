"""Read-only reductions of accepted states; no topology or loading mutations."""
from common import *
from multifragment_shear_objective import metrics,BINS,SPEC,shape
from current_surface import CurrentSurface
from shear_field_library import plateau
from bridge_audit import inspect,candidate_regions
import ast,csv,gzip,time

graph=module('result_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')
# Reuse the already reviewed failed-edge lineage and local deformation algorithms.
_tree=ast.parse((HA/'analyze.py').read_text())
_names={'histories','crack_audit','deformation'}
exec(compile(ast.Module(body=[q for q in _tree.body if isinstance(q,ast.FunctionDef) and q.name in _names],type_ignores=[]),str(HA/'analyze.py'),'exec'))
def read(p):return json.loads(Path(p).read_text())
def csvout(path,rows,columns=None):
    with Path(path).open('w') as f:
        w=csv.DictWriter(f,fieldnames=columns or list(rows[0]));w.writeheader();w.writerows(rows)

def all_material_failed_networks(case,ref):
    """Inclusive failed-edge graph, with original bulk/interface types intact.

    Fixed-anchor endpoints remain vertices, never remapped to material nodes.
    This separate view complements the historical bulk-only lineage definition.
    """
    b=ref['bonds'];bt=ref['bond_types'];parent=np.arange(len(ref['X']));clusters={};births=[];merges=[];history=[];failed=set();peak=0
    def root(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return int(i)
    def join(i,j):
        a,c=root(i),root(j)
        if a!=c:parent[max(a,c)]=min(a,c)
    with gzip.open(case/'damage_events.jsonl.gz','rt') as f:
        for line in f:
            e=json.loads(line);ids=np.asarray(e['original_bond_ids'])-1;dam=np.asarray(e['d']);new=ids[(dam==1)&np.isin(bt[ids],[1,2])&~ref['initial_cut'][ids]]
            if not len(new):continue
            old={root(key):set(value) for key,value in clusters.items()}
            for a,c in b[new]:join(int(a),int(c))
            merged={}
            for key,value in old.items():merged.setdefault(root(key),set()).update(value)
            born={root(int(i)) for i in b[new].ravel()}-set(merged)
            for key in sorted(born):
                cid=len(births)+1;merged[key]={cid};group=np.asarray([i for i in new if root(int(b[i,0]))==key],int)
                births.append(dict(cluster_id=cid,time_s=e['time_s'],step=e['step'],first_original_bond_ids=group+1,original_bond_types=bt[group],reference_centroid_m=ref['X'][b[group]].mean((0,1))))
            for key,value in merged.items():
                ancestors=[p for p in old if root(p)==key]
                if len(ancestors)>1:merges.append(dict(time_s=e['time_s'],step=e['step'],n_merged_components=len(ancestors),original_nucleation_ids=sorted(value),new_bridge_original_bond_ids=[int(i+1) for i in new if root(int(b[i,0]))==key]))
            clusters=merged;failed.update(new.tolist());peak=max(peak,len(clusters))
            history.append(dict(time_s=e['time_s'],step=e['step'],separate_failed_material_networks=len(clusters),new_bulk_failures=int(np.sum(bt[new]==1)),new_interface_failures=int(np.sum(bt[new]==2)),new_original_bond_ids=new+1,original_bond_types=bt[new]))
    ids=np.array(sorted(failed),int);edge=b[ids]
    matrix=graph.coo_matrix((np.ones(2*len(edge)),(np.r_[edge[:,0],edge[:,1]],np.r_[edge[:,1],edge[:,0]])),shape=(len(ref['X']),len(ref['X']))).tocsr()
    _,labels=graph.csgraph.connected_components(matrix,directed=False)
    assert len(set(labels[np.unique(edge)].tolist()))==len(clusters)
    return dict(definition='All naturally failed original bulk and interface bonds, shared original endpoints; fixed-anchor endpoints retained and labelled, initial precuts excluded. Contact never enters either graph. This inclusive graph supplements, not replaces, historical bulk-only birth counts.',n_births=len(births),n_merge_events=len(merges),max_simultaneous_networks=peak,final_networks=len(clusters),independent_final_graph_check='PASS',new_failed_bulk_bonds=int(np.sum(bt[ids]==1)),new_failed_interface_bonds=int(np.sum(bt[ids]==2)),births=births,merges=merges,history=history,original_failed_bond_ids=ids+1,original_failed_bond_types=bt[ids])

def analyze_case(name,ref):
    case=O/'cases'/name;summary=read(case/'SUMMARY.json');rows=histories(case)
    registry={int(k):v for k,v in read(case/'COMPONENT_REGISTRY.json').items()}
    final=load_npz(case/'keys'/('final.npz' if summary['status']=='COMPLETED' else 'last_accepted.npz'))
    n=int(ref['clot_count']);lab,sizes,att=graph.topology(ref,final['d']);free=sizes[~att]
    first={};last={};trace={};cohort=None;endhistory=None;max_counts={key:0 for key in ['N_ge2','N_ge14','N_ge20','N_ge34','N_ge114','N_ge200um','N_ge300um']}
    with gzip.open(case/'CONNECTED_COMPONENT_HISTORY.jsonl.gz','rt') as stream:
        for line in stream:
            z=json.loads(line);t=z['time_s'];endhistory=z
            if abs(t-SPEC['common_cohort_time_s'])<1e-13:cohort=z
            for q in z['components']:
                cid=q['component_id'];first.setdefault(cid,dict(time_s=t,**q));last[cid]=dict(time_s=t,**q)
                if registry[cid]['nodes']>=2:
                    # Full shape/coordinates remain in the compressed accepted
                    # history. Avoid a second uncompressed copy for plotting.
                    trace.setdefault(cid,[]).append(dict(time_s=t,**{key:q[key] for key in ['internal_survival','mean_internal_damage','internal_alive_bonds','internal_initial_bonds','internal_relative_kinetic_J','centroid_velocity_m_s','remaining_live_external_connections']}))
    for key in max_counts:max_counts[key]=int(max(q[key] for q in rows))
    stable=[];splits=summary['secondary_events'];split_by_parent={q['parent_component_id']:q for q in splits}
    for cid,c in registry.items():
        ids=np.asarray(c['atom_ids'])-1
        same=bool(np.all(lab[ids]==lab[ids[0]]) and not att[lab[ids[0]]] and sizes[lab[ids[0]]]==len(ids))
        lasttime=last[cid]['time_s'];split=split_by_parent.get(cid)
        stable.append(dict(fragment_id=cid,**c,mass_kg=c['nodes']*RHO*V,
            first_state=first[cid],last_state=last[cid],last_observed_s=lasttime,
            observed_lifetime_s=lasttime-c['birth_time_s'],same_original_members_connected_at_final=same,
            status='CONNECTED_THROUGH_OBSERVED_WINDOW' if same else 'SPLIT_WITHIN_OBSERVED_WINDOW',
            secondary_split_history=split,split_time_bracket_s=[lasttime,split['time_s']] if split else None,
            common_post_unload_observed_s=max(0,lasttime-max(c['birth_time_s'],SPEC['common_cohort_time_s'])),
            lifetime_scope='Accepted-state observation, right-censored at final time; not permanent stability'))
    stableids={q['fragment_id'] for q in stable if q['same_original_members_connected_at_final']}
    cohort_ids=[q['component_id'] for q in cohort['components'] if 20<=registry[q['component_id']]['nodes']<=200] if cohort else []
    cohort_metrics=metrics([registry[q['component_id']]['nodes'] for q in cohort['components']],n) if cohort else None
    final_metrics=metrics(free,n)
    criteria=dict(common_cohort=cohort_metrics['independent_criteria'] if cohort_metrics else None,
        final=final_metrics['independent_criteria'],five_cohort_memberships_survive=sum(cid in stableids for cid in cohort_ids)>=5,
        numerical_completion=summary['status']=='COMPLETED')
    target=bool(cohort_metrics and cohort_metrics['instant_target'] and final_metrics['instant_target'] and criteria['five_cohort_memberships_survive'] and criteria['numerical_completion'])
    finalcs=[dict(**registry[q['component_id']],**q,mass_kg=registry[q['component_id']]['nodes']*RHO*V) for q in endhistory['components']]
    networks=crack_audit(case,ref);networks['all_material_including_interface']=all_material_failed_networks(case,ref);diagnostics={}
    inclusive=networks['all_material_including_interface'];wanted=np.flatnonzero((final['d']==1)&~ref['initial_cut']&np.isin(ref['bond_types'],[1,2]))+1
    assert np.array_equal(inclusive['original_failed_bond_ids'],wanted)
    assert inclusive['new_failed_bulk_bonds']==rows[-1]['bulk_failed'] and inclusive['new_failed_interface_bonds']==rows[-1]['interface_failed']
    centers=np.asarray([q['reference_centroid_m'] for q in networks['births']])
    separation=float(np.linalg.norm(centers[:,None]-centers[None,:],axis=2).max()) if len(centers)>1 else 0.
    for ev in ['first_softening','first_failure']:
        if ev in summary['events']:
            e=summary['events'][ev];state=load_npz(case/e['evidence']);diagnostics[ev]=dict(event=e,local_affine=deformation(ref,state,np.asarray(e['atom_ids']).ravel()-1))
    # These labels describe the data, never replace the independent frozen gates.
    classification=[]
    if not rows[-1]['softened']:classification.append('NO_DAMAGE')
    elif not free.sum():classification.append('DAMAGE_WITHOUT_DETACHMENT')
    if free.sum() and final_metrics['fine_debris_mass_fraction_detached']>.2:
        classification.append('SURFACE_EROSION' if free.sum()/n<.1 else 'DEBRIS_DOMINATED_DETACHMENT')
    if final_metrics['fragment_count_target_band']>=2:classification.append('MULTIPLE_TARGET_BAND_COMPONENTS')
    if free.sum() and final_metrics['largest_free_fragment_fraction_detached']>.5:classification.append('DOMINANT_FREE_COMPONENT')
    if splits:classification.append('SECONDARY_FRAGMENTATION')
    if summary['status']!='COMPLETED':classification.append('NUMERICALLY_UNRESOLVED')
    if not classification:classification.append('DETACHMENT_WITHOUT_EXPLORATORY_TARGET')
    result=dict(name=name,status=summary['status'],classification=classification,exploratory_target_met=target,
        target_criteria=criteria,final=rows[-1],final_metrics=final_metrics,final_components=finalcs,
        first_softening_s=summary['events'].get('first_softening',{}).get('time_s'),
        first_failure_s=summary['events'].get('first_failure',{}).get('time_s'),
        first_multinode_s=summary['events'].get('first_multinode',{}).get('time_s'),
        history_max_N20_200=summary['max_medium_components'],history_max_free_nodes=summary['history_max_free_nodes'],
        history_max_size_threshold_counts=max_counts,common_cohort_medium_ids=cohort_ids,
        common_cohort_medium_surviving=sum(cid in stableids for cid in cohort_ids),common_cohort_metrics=cohort_metrics,
        cohort_time_s=SPEC['common_cohort_time_s'],upper_drive_stop_s=OFF,remaining_drive_fraction_at_cohort=0.,passive_drag_remains=True,
        secondary_events=splits,secondary_after_common_cohort=sum(q['time_s']>=SPEC['common_cohort_time_s'] for q in splits),
        crack_networks=dict(n_births=networks['n_births'],n_merges=networks['n_merge_events'],max_simultaneous=networks['max_simultaneous_networks'],max_birth_centroid_separation_m=separation,separated_birth_centroids_ge_horizon=separation>=HORIZON),
        inclusive_failed_material_networks={key:networks['all_material_including_interface'][key] for key in ['n_births','n_merge_events','max_simultaneous_networks','new_failed_bulk_bonds','new_failed_interface_bonds']},
        max_relative_energy_residual=summary['max_relative_energy_residual'],
        peak_shear_Pa=max(q['shear_peak_Pa'] for q in rows),peak_sum_nodal_force_N=max(q['sum_force_norm_N'] for q in rows),geometry_first_diagnostic=summary['layer_geometry_first_diagnostic'],
        source_sha256={p.name:sha(p) for p in [case/'SUMMARY.json',case/'history.csv',case/'CONNECTED_COMPONENT_HISTORY.jsonl.gz',case/'COMPONENT_REGISTRY.json']},
        scope='Original surviving bulk/interface bonds and anchor reachability; contact excluded. Same frozen material, nominal shear cap, time envelope and follower coupling; different spatial mode. Driver is zero after 2.8ms but passive slip drag remains on retained outer facets.')
    dump(case/'ANALYSIS.json',result);dump(case/'CRACK_NETWORK_AUDIT.json',networks)
    dump(case/'LOCAL_DEFORMATION_AUDIT.json',diagnostics)
    dump(case/'FRAGMENT_HISTORY.json',dict(fragments=stable,candidates_ge14=[q for q in stable if q['nodes']>=14],secondary_events=splits,common_cohort_ids=cohort_ids))
    dump(case/'FRAGMENT_TRACE.json',trace)
    return result,summary,rows,networks

def candidate_audit(ref,s):
    proposals,screen=graph.proposals(ref,s);out=[];n=int(ref['clot_count']);lab,sizes,att=graph.topology(ref,s['d'])
    for members,origins in proposals.items():
        ids=np.asarray(members);boundary,inside=graph.boundary(ref,ids)
        boundary=boundary[~ref['initial_cut'][boundary]]
        inside=inside[(ref['bond_types'][inside]==1)&~ref['initial_cut'][inside]]
        shape0=graph.compactness(ref['X'][ids],DX);survival=float(np.mean(s['d'][inside]<1)) if len(inside) else 0.
        loss=float(np.mean(s['d'][boundary])) if len(boundary) else 1.
        internal=float(np.mean(s['d'][inside])) if len(inside) else 0.
        cap=float(np.sum(1-s['d'][boundary]));degree=cap+2*float(np.sum(1-s['d'][inside]));ratio=cap/max(degree,1e-30)
        gates=dict(size=3<=len(ids)<=n//2,reference_compact=shape0['voxel_bbox_fill']>=.2,internal_survival=survival>=.8,boundary_weaker=loss>=.5 and loss>internal,low_conductance=ratio<=.15)
        out.append(dict(original_atom_ids=ids+1,nodes=len(ids),origins=origins,reference_compactness=shape0,current_shape=shape(s['x'][ids]),internal_survival=survival,internal_mean_damage=internal,boundary_integrity_loss=loss,weighted_boundary_to_degree=ratio,remaining_live_boundary_bonds=int(np.sum(s['d'][boundary]<1)),anchor_reachable_nodes=int(att[lab[ids]].sum()),gates=gates,near_separation_proxy_pass=all(gates.values())))
    for q in out:
        ids=np.asarray(q['original_atom_ids'])-1
        q['internal_live_components']=int(graph.internal_components(ref,ids,s['d']))
        q['attached_target_band_3D_candidate']=bool(q['near_separation_proxy_pass'] and 20<=q['nodes']<=200 and q['anchor_reachable_nodes']>0 and q['internal_live_components']==1 and q['current_shape']['shape_class']=='THREE_DIMENSIONAL_POINT_CLOUD')
    return dict(candidates=out,integrity_threshold_screen=screen,passing=sum(q['near_separation_proxy_pass'] for q in out),attached_target_band_3D_candidate_count=sum(q['attached_target_band_3D_candidate'] for q in out),scope='Read-only threshold and spectral partitions of actual damage. Frozen historical heuristic, not true fragments, not material bond deletion and not a load mask. Optimizer review further requires attached, internally connected, 20–200 nodes and a 3D point cloud; duplicate/overlapping candidates are not multiple fragments.')

def main():
    start=time.perf_counter();ref=load_npz(O/'provenance/reference.npz');n=int(ref['clot_count'])
    results={};summaries={};hist={};networks={};bridges={};candidates={};distribution=[];bridge_rows=[];surface={}
    for name in CASES:
        q,summary,rows,network=analyze_case(name,ref);results[name]=q;summaries[name]=summary;hist[name]=rows;networks[name]=network
        dump(O/('UNIFORM_SHEAR_RESULTS.json' if name=='uniform' else 'COUNTER_ROTATING_SHEAR_RESULTS.json'),q)
        tool=CurrentSurface(ref,name);bridges[name]={};candidates[name]={};surface[name]=[]
        for stage in ['max_multifragment','common_cohort','final']:
            path=O/'cases'/name/'keys'/(stage+'.npz')
            if not path.exists():continue
            s=load_npz(path);_,sizes,att=graph.topology(ref,s['d']);free=sizes[~att]
            for lo,hi in BINS:
                vals=free[(free>=lo)&((free<=hi) if hi else True)];total=int(vals.sum())
                distribution.append(dict(case=name,state=stage,time_s=float(s['time']),size_bin=f'{lo}-{hi}' if hi else '201+',component_count=len(vals),material_nodes=total,mass_kg=total*RHO*V,mass_fraction_initial=total/n,mass_fraction_detached=total/int(free.sum()) if free.sum() else None))
        for stage in ['peak_load','final']:
            path=O/'cases'/name/'keys'/(stage+'.npz')
            if not path.exists():continue
            s=load_npz(path);items=[]
            for region in candidate_regions(ref):
                info,rr=inspect(ref,s,region);items.append(info);bridge_rows.extend(dict(case=name,state=stage,**row) for row in rr)
            bridges[name][stage]=dict(time_s=float(s['time']),source_sha256=sha(path),regions=items)
            candidates[name][stage]=candidate_audit(ref,s)
        for path in sorted((O/'cases'/name/'keys').glob('*.npz')):
            if 'failed_unaccepted' in path.name:continue
            s=load_npz(path);z=tool.evaluate(s['x'],s['velocity'],s['d'],float(s['time']))
            F=z['force'];ff=z['face_force'];lab,sizes,att=graph.topology(ref,s['d']);supported=np.unique(tool.tri[z['active']]);free=np.flatnonzero(~att[lab])
            budget=np.asarray(tool.W@np.linalg.norm(ff,axis=1));mag=np.linalg.norm(F,axis=1)
            surface[name].append(dict(key=path.stem,time_s=float(s['time']),**tool.metrics(z,s['x']),force_error_N=float(np.max(abs(F.sum(0)-ff.sum(0)))),moment_error_Nm=float(np.max(abs(np.cross(s['x'],F).sum(0)-np.cross(z['centroid'],ff).sum(0)))),power_error_W=float(abs(np.sum(F*s['velocity'])-np.sum(ff*z['face_velocity']))),local_force_vs_positive_budget_max=float(np.max(mag/np.maximum(budget,1e-300))),free_nodes_with_retained_surface=int(np.isin(free,supported).sum()),free_nodes_without_resolved_surface=int((~np.isin(free,supported)).sum()),native_force_replay_error_N=float(np.max(abs(F-s['external_force'])))))
            offset=z['centroid']-tool.field.center;chi=plateau(offset@tool.field.e1,400e-6,700e-6)*plateau(offset@tool.field.e2,250e-6,450e-6)*plateau(offset[:,2],100e-6,300e-6)
            active=z['active'];vt=z['face_velocity']-z['normal']*np.sum(z['face_velocity']*z['normal'],axis=1)[:,None]
            surface[name][-1].update(retained_faces_outside_drive_region=int(np.sum(active&(chi==0))),retained_faces_inside_full_plateau=int(np.sum(active&(chi==1))),retained_face_window_max=float(chi[active].max(initial=0)),outer_speed_max_m_s=float(np.linalg.norm(z['outer_velocity'][active],axis=1).max(initial=0)),tangential_surface_speed_max_m_s=float(np.linalg.norm(vt[active],axis=1).max(initial=0)),sum_drive_node_force_N=float(np.linalg.norm(z['drive_force'],axis=1).sum()),sum_drag_node_force_N=float(np.linalg.norm(z['drag_force'],axis=1).sum()),sum_actual_node_force_N=float(np.linalg.norm(z['force'],axis=1).sum()))
    csvout(O/'FRAGMENT_SIZE_DISTRIBUTION.csv',distribution)
    csvout(O/'SURVIVING_BRIDGE_BONDS.csv',bridge_rows)
    dump(O/'CRACK_BRIDGE_AUDIT.json',dict(cases=bridges,diagnostic_only=True,pinned_min_cut_calls=12))
    dump(O/'NEAR_SEPARATION_AUDIT.json',candidates);dump(O/'CURRENT_SURFACE_TRAJECTORY_AUDIT.json',surface)
    dump(O/'ANALYSIS.json',results);dump(O/'NATIVE_SHEAR_CASE_RESULTS.json',summaries);dump(O/'DAMAGE_NETWORK_HISTORY.json',networks)
    dump(O/'FRACTURE_EVENTS.json',{k:v['events'] for k,v in summaries.items()})
    dump(O/'CONNECTED_COMPONENT_HISTORY.json',dict(cases={k:dict(history=str((O/'cases'/k/'CONNECTED_COMPONENT_HISTORY.jsonl.gz').relative_to(O)),registry=read(O/'cases'/k/'COMPONENT_REGISTRY.json'),fragment_lifetimes=read(O/'cases'/k/'FRAGMENT_HISTORY.json'),records=len(hist[k])+1) for k in results},definition='Index and exact original memberships; every accepted state in lossless JSONL gzip. Contact does not connect fragments.'))
    ledger={};requirements=[];windows={}
    for name,q in results.items():
        rows=hist[name];end=q['final'];initial=load_npz(O/'cases'/name/'keys/initial.npz')
        ledger[name]={key:end[key] for key in ['W_J','W_drive_J','W_drag_J','W_endpoint_power_J','work_path_vs_endpoint_J','U_J','K_J','D_J','Uc_J','Jc_J','R_J']}
        ledger[name].update(initial_U_J=float(initial['scalars'][0]),initial_K_J=float(initial['K']),initial_Uc_J=float(initial['Uc']),max_relative_residual=q['max_relative_energy_residual'])
        tests={'Re_h_outer_gt_1':lambda r:r['Re_h_outer_max']>1,'Re_h_slip_gt_1':lambda r:r['Re_h_slip_max']>1,'displacement_gt_0p1h':lambda r:r['max_loaded_displacement_m']>.1*H_EFF,'normal_rotation_gt_15deg':lambda r:r['max_normal_rotation_deg']>15,'surface_support_loss':lambda r:r['inactive_broken_edge_faces']>0,'area_ratio_outside_half_to_two':lambda r:r['max_active_area_ratio']>2 or (r['min_active_area_ratio'] is not None and r['min_active_area_ratio']<.5)}
        windows[name]={key:next((dict(time_s=r['time_s'],step=r['step']) for r in rows if pred(r)),None) for key,pred in tests.items()}
        windows[name].update(scope='Independent approximation diagnostics, not a calibrated clinical validity boundary. Numerical completion admits only this idealized mechanical model. New wetting and hydrodynamics of unrepresented fragment surfaces remain unresolved.',retained_reference_area_fraction_final=end['active_reference_area_fraction'],free_nodes_not_vertices_of_retained_facets_final=end['free_nodes_without_resolved_surface'])
        requirements.append(dict(case=name,nominal_stationary_shear_cap_Pa=1000,outer_velocity_cap_m_s=25,mu_Pa_s=MU,h_eff_m=H_EFF,actual_peak_slip_shear_Pa=q['peak_shear_Pa'],peak_drive_shear_Pa=max(r['actual_peak_drive_shear_Pa'] for r in rows),peak_sum_nodal_force_N=q['peak_sum_nodal_force_N'],max_Re_h_outer=max(r['Re_h_outer_max'] for r in rows),max_Re_h_slip=max(r['Re_h_slip_max'] for r in rows),work_J=end['W_J'],drive_work_J=end['W_drive_J'],drag_work_J=end['W_drag_J'],history_max_N20_200=q['history_max_N20_200'],final_N20_200=end['N20_200'],final_detached_nodes=end['detached_nodes'],final_fines_fraction=end['fine_debris_mass_fraction_detached'],label=LABEL))
    dump(O/'MODEL_VALIDITY_WINDOWS.json',windows);csvout(O/'SHEAR_REQUIREMENTS.csv',requirements)
    dump(O/'ENERGY_LEDGER.json',dict(identity='W_drive + W_drag = Delta U_LPS + Delta K + D + Delta U_contact - J_contact_activation + R',work='Current surface force dotted into actual accepted native drift; Gauss4 quadrature with edge-failure root splits. Gauss8 checks and independent endpoint-power sum retained.',cases=ledger,numerical_residual_is_not_dissipation=True,no_extra_numerical_damping=True,passive_drag_after_driver_off=True))
    target=any(q['exploratory_target_met'] for q in results.values());band=any(q['history_max_N20_200']>0 for q in results.values());near=sum(stage['passing'] for case in candidates.values() for stage in case.values())
    credible=sum(stage['attached_target_band_3D_candidate_count'] for case in candidates.values() for stage in case.values())
    readiness='BOUNDED_FOLLOWUP_CANDIDATE_REVIEW' if band or credible else 'NO_USEFUL_TARGET_OR_NEAR_SEPARATION_SIGNAL'
    status=dict(SHEAR_LIBRARY_STATUS='TWO_EXECUTABLE_MODES',CURRENT_SURFACE_COUPLING_STATUS='ORIGINAL_OUTER_FACETS_FOLLOW_CURRENT_STATE',LOCAL_TRACTION_MAPPING_STATUS='POSITIVE_VERTEX_MAP_CHECKED_SMALL_SPECIMEN_ONLY',NATIVE_SHEAR_DYNAMICS_STATUS={k:v['status'] for k,v in results.items()},MULTIPLE_CRACK_NUCLEATION_STATUS={k:v['crack_networks'] for k,v in results.items()},CRACK_PERCOLATION_STATUS='FAILED_EDGE_MERGING_IS_NOT_CLOSED_CRACK_SURFACE_PROOF',MULTIFRAGMENT_FORMATION_STATUS='EXPLORATORY_TARGET_MET' if target else 'EXPLORATORY_TARGET_NOT_MET',TARGET_FRAGMENT_COUNT={k:dict(maximum=v['history_max_N20_200'],final=v['final']['N20_200']) for k,v in results.items()},FINE_DEBRIS_FRACTION={k:v['final_metrics']['fine_debris_mass_fraction_detached'] for k,v in results.items()},FRAGMENT_COHESION_STATUS='EXACT_MEMBERSHIP_INTERNAL_BONDS_AND_OBSERVED_LIFETIMES_REPORTED',ENERGY_ADMISSION_STATUS={k:dict(max_relative_R=v['max_relative_energy_residual'],within_1_percent=v['max_relative_energy_residual']<=.01) for k,v in results.items()},OPTIMIZATION_READINESS=readiness,ACOUSTIC_FEASIBILITY_STATUS='NOT_ESTABLISHED',SCIENTIFIC_CONCLUSION='Two frozen follower-shear modes compared using true native fracture. No acoustic realizability or optimum is inferred; target success and near-separation evidence must be assessed separately.',near_separation_proxy_count=near,new_optimizer_trials=0,new_main_native_cases=2,next_V3_admitted=False,scientific_label=LABEL,new_crack_wetting='NEW_CRACK_FACE_WETTING_UNRESOLVED')
    status['attached_target_band_3D_candidate_observations']=credible
    a=results['uniform'];b=results['counter_rotating']
    if not band:
        status['SCIENTIFIC_CONCLUSION']=f'两个冻结工况均未形成20–200节点凝聚性游离块。单向剪切历史最大块{a["history_max_free_nodes"]}节点，双涡旋{b["history_max_free_nodes"]}节点；本次双涡旋没有展示多目标碎块优势，尚未识别成功剪切需求。表面支撑、空间覆盖、相对速度和存活三维桥接限制推断；不支持盲目扩大BoTorch搜索。'
    dump(O/'FINAL_STATUS.json',status)
    dump(O/'RUN_RESOURCES.json',dict(cases={k:{key:v[key] for key in ['wall_s','cpu_s','peak_RSS_bytes','accepted_steps','predictive_substeps','predictive_rejections']} for k,v in summaries.items()},native_wall_s_sum=sum(v['wall_s'] for v in summaries.values()),native_cpu_s=sum(v['cpu_s'] for v in summaries.values()),max_single_case_peak_RSS_bytes=max(v['peak_RSS_bytes'] for v in summaries.values()),sum_individual_peak_RSS_upper_bound_bytes=sum(v['peak_RSS_bytes'] for v in summaries.values()),parallel_cases=True,ranks_per_case=1,threads_per_case=1,new_main_native_cases=2,new_full_V3_cases=0,optimizer_trials=0,analysis_wall_s=time.perf_counter()-start))
    print('ANALYSIS_DONE',status)

if __name__=='__main__':main()
