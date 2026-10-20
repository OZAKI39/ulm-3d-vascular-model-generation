"""Read-only reductions of accepted states; no topology or loading mutations."""
from common import *
from multifragment_objective import metrics,BINS,SPEC,shape
from distributed_multisite_loading import DistributedLoading
from offline_crack_bridge_audit import inspect,candidate_regions
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
    clear=DistributedLoading(ref,name).clear_time
    contacts=[q for q in rows if q['contact_nodes']>0];lastcontact=contacts[-1]['time_s'] if contacts else None
    contact_after=any(q['time_s']>=clear-1e-14 for q in contacts)
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
        cohort_time_s=SPEC['common_cohort_time_s'],scheduled_return_s=clear,
        contact_after_scheduled_return=contact_after,last_actual_contact_s=lastcontact,
        actual_contact_free_tail_s=max(0,float(final['time'])-max(clear,lastcontact or 0)),
        secondary_events=splits,secondary_after_common_cohort=sum(q['time_s']>=SPEC['common_cohort_time_s'] for q in splits),
        crack_networks=dict(n_births=networks['n_births'],n_merges=networks['n_merge_events'],max_simultaneous=networks['max_simultaneous_networks'],max_birth_centroid_separation_m=separation,separated_birth_centroids_ge_horizon=separation>=HORIZON),
        inclusive_failed_material_networks={key:networks['all_material_including_interface'][key] for key in ['n_births','n_merge_events','max_simultaneous_networks','new_failed_bulk_bonds','new_failed_interface_bonds']},
        max_relative_energy_residual=summary['max_relative_energy_residual'],
        peak_site_reactions_N=[max(q[f'site{j}_reaction_norm_N'] for q in rows) for j in range(2)],
        max_site_overlap_dx=[max(q[f'site{j}_max_overlap_m'] for q in rows)/DX for j in range(2)],
        source_sha256={p.name:sha(p) for p in [case/'SUMMARY.json',case/'history.csv',case/'CONNECTED_COMPONENT_HISTORY.jsonl.gz',case/'COMPONENT_REGISTRY.json']},
        scope='Original surviving bulk/interface bonds and anchor reachability; contact excluded. Same prescribed kinematics, not equal work. Post-return does not imply contact-free until actual last contact.')
    dump(case/'ANALYSIS.json',result);dump(case/'CRACK_NETWORK_AUDIT.json',networks)
    dump(case/'LOCAL_DEFORMATION_AUDIT.json',diagnostics)
    dump(case/'FRAGMENT_HISTORY.json',dict(fragments=stable,candidates_ge14=[q for q in stable if q['nodes']>=14],secondary_events=splits,common_cohort_ids=cohort_ids))
    dump(case/'FRAGMENT_TRACE.json',trace)
    return result,summary,rows,networks

def main():
    start=time.perf_counter();ref=load_npz(O/'provenance/reference.npz');results={};summaries={};hist={};networks={}
    for name in ['simultaneous','staggered']:results[name],summaries[name],hist[name],networks[name]=analyze_case(name,ref)
    records=[];distribution=[]
    for name,result in results.items():
        fm=result['final_metrics'];end=result['final']
        records.append(dict(case=name,source='new_uniform_two_site',status=result['status'],final_time_s=end['time_s'],detached_nodes=int(end['detached_nodes']),largest_nodes=int(end['largest_nodes']),**{k:v for k,v in fm.items() if k!='independent_criteria'},history_max_N20_200=result['history_max_N20_200'],history_max_free_nodes=result['history_max_free_nodes'],target_met=result['exploratory_target_met'],W_J=end['W_J'],D_J=end['D_J'],K_J=end['K_J'],max_relative_R=result['max_relative_energy_residual']))
        for phase,key in [('common_cohort','common_cohort'),('final','final')]:
            state=load_npz(O/'cases'/name/'keys'/(key+'.npz'));_,sizes,att=graph.topology(ref,state['d']);free=sizes[~att]
            for lo,hi in BINS:
                values=free[(free>=lo)&((free<=hi) if hi else True)];count=int(values.sum())
                distribution.append(dict(case=name,state=phase,time_s=float(state['time']),size_bin=f'{lo}-{hi}' if hi else f'{lo}+',component_count=len(values),material_nodes=count,mass_kg=count*RHO*V,mass_fraction_initial=count/int(ref['clot_count']),mass_fraction_detached=count/int(free.sum()) if free.sum() else None))
    old=read(HISTORY/'cases/uniform/ANALYSIS.json');oldsum=read(HISTORY/'cases/uniform/SUMMARY.json');oldfinal=old['final'];oldstate=load_npz(HISTORY/'cases/uniform/keys/final.npz');_,ss,aa=graph.topology(ref,oldstate['d']);om=metrics(ss[~aa],int(ref['clot_count']))
    records.append(dict(case='historical_single_uniform',source='read_only_historical_different_contact_geometry',status=old['status'],final_time_s=oldfinal['time_s'],detached_nodes=int(oldfinal['detached_nodes']),largest_nodes=int(oldfinal['largest_nodes']),**{k:v for k,v in om.items() if k!='independent_criteria'},history_max_N20_200=old['history_max_N20_200'],history_max_free_nodes=old['history_max_free_nodes'],target_met=False,W_J=oldfinal['W_J'],D_J=oldfinal['D_J'],K_J=oldfinal['K_J'],max_relative_R=old['max_relative_energy_residual']))
    csvout(O/'MULTIFRAGMENT_RESULTS.csv',records);csvout(O/'FRAGMENT_SIZE_DISTRIBUTIONS.csv',distribution)
    fragment_history={k:read(O/'cases'/k/'FRAGMENT_HISTORY.json') for k in results}
    dump(O/'ANALYSIS.json',results);dump(O/'FRAGMENT_HISTORY.json',fragment_history)
    measurements=[]
    for name,history in fragment_history.items():
        for q in history['fragments']:
            last=q['last_state']
            measurements.append(dict(case=name,fragment_id=q['fragment_id'],nodes=q['nodes'],birth_time_s=q['birth_time_s'],last_observed_s=q['last_observed_s'],observed_lifetime_s=q['observed_lifetime_s'],mass_kg=q['mass_kg'],reference_volume_m3=q['reference_volume_m3'],equivalent_diameter_m=q['equivalent_diameter_m'],internal_alive_bonds=last['internal_alive_bonds'],internal_bond_survival=last['internal_survival'],external_surviving_bonds=last['remaining_live_external_connections'],anchor_reachable=last['anchor_reachable'],shape_class=last['shape_class'],radius_gyration_m=last['radius_gyration_m'],span_x_m=last['actual_span_m'][0],span_y_m=last['actual_span_m'][1],span_z_m=last['actual_span_m'][2],principal_axes_json=json.dumps(last['principal_axes_columns']),principal_rms_lengths_json=json.dumps(last['principal_rms_lengths_m']),internal_edge_density_m3=last['internal_alive_edge_density_m3'],original_atom_ids=';'.join(map(str,q['atom_ids'])),connected_at_final=q['same_original_members_connected_at_final']))
    if measurements:csvout(O/'FRAGMENT_MEASUREMENTS.csv',measurements)
    dump(O/'DAMAGE_NETWORK_HISTORY.json',networks)
    csvout(O/'NEW_DAMAGE_NETWORK_HISTORY.csv',[dict(case=k,**q) for k,v in networks.items() for q in v['history']])
    csvout(O/'ALL_MATERIAL_DAMAGE_NETWORK_HISTORY.csv',[dict(case=k,**{field:q[field] for field in ['time_s','step','separate_failed_material_networks','new_bulk_failures','new_interface_failures']}) for k,v in networks.items() for q in v['all_material_including_interface']['history']])
    historical_identity={}
    for name in ['uniform','heterogeneous']:
        p=HISTORY/'cases'/name/'keys/final.npz';state=load_npz(p);ids=np.flatnonzero((state['d']==1)&~ref['initial_cut']&np.isin(ref['bond_types'],[1,2]));interface=int(np.sum(ref['bond_types'][ids]==2));oldsummary=read(HISTORY/'cases'/name/'SUMMARY.json')
        assert interface==oldsummary['last_row']['interface_failed']==0
        historical_identity[name]=dict(source=str(p),source_sha256=sha(p),original_failed_bond_ids=ids+1,original_bond_types=ref['bond_types'][ids],interface_failed=interface,bulk_only_equals_all_naturally_failed_material_edges=True)
    dump(O/'FAILURE_GRAPH_IDENTITY_AUDIT.json',dict(historical=historical_identity,new_cases={k:v['inclusive_failed_material_networks'] for k,v in results.items()},scope='Both historical cases have zero natural interface failures, so their bulk failed-edge graphs are already complete. New interface failures retain original identity in the inclusive graphs.'))
    dump(O/'ENERGY_LEDGER.json',dict(identity='W0+W1 = Delta U_LPS + Delta K + D + Delta U_pair - J_pair_activation + Delta U_tools + R',no_damping=True,history_includes_independent_per_source_forces_and_work=True,normalization='abs(R)/running_max(abs(W),abs(Delta U_LPS)+K+D+U_pair+U_tools), gate 0.01; not normalized by final work',cases={k:dict(final=v['final'],peak_site_reactions_N=v['peak_site_reactions_N'],max_relative_R=v['max_relative_energy_residual']) for k,v in results.items()}))
    after={};bridge_rows=[]
    for name in results:
        after[name]={}
        for key in ['max_multifragment','final']:
            p=O/'cases'/name/'keys'/(key+'.npz');s=load_npz(p);items=[]
            for region in candidate_regions(ref):
                info,rr=inspect(ref,s,region)
                ids=np.asarray(region['original_atom_ids'])-1;_,inside=graph.boundary(ref,ids);inside=inside[(ref['bond_types'][inside]==1)&~ref['initial_cut'][inside]]
                info.update(candidate_label='TETHERED_DIAGNOSTIC_REGION' if info['anchor_reachable_nodes'] else 'UNANCHORED_REGION_NOT_NECESSARILY_ONE_COMPONENT',internal_initial_bonds=len(inside),internal_alive_bonds=int(np.sum(s['d'][inside]<1)),mean_internal_damage=float(s['d'][inside].mean()) if len(inside) else None)
                items.append(info);bridge_rows.extend(dict(case=name,state=key,time_s=float(s['time']),**q) for q in rr)
            after[name][key]=dict(time_s=float(s['time']),source_sha256=sha(p),regions=items)
    dump(O/'CRACK_BRIDGE_AFTER_LOADING.json',dict(cases=after,pinned_min_cut_calls=12,diagnostic_only=True));csvout(O/'SURVIVING_BRIDGE_BONDS_AFTER_LOADING.csv',bridge_rows)
    resources=dict(cases={k:{q:v[q] for q in ['wall_s','cpu_s','peak_RSS_bytes','accepted_steps','predictive_substeps','predictive_rejections']} for k,v in summaries.items()},native_wall_s=sum(v['wall_s'] for v in summaries.values()),native_cpu_s=sum(v['cpu_s'] for v in summaries.values()),peak_RSS_bytes=max(v['peak_RSS_bytes'] for v in summaries.values()),rank_per_case=1,threads_per_case=1,new_formal_cases=2,additional_formal_cases=0,full_V3_runs=0,CFD_runs=0,offline_initial_wall_s=read(O/'CRACK_BRIDGE_AUDIT.json')['wall_s'],build_wall_s=read(O/'BUILD_AUDIT.json')['wall_s'],analysis_wall_s=time.perf_counter()-start,scope='Formal driver wall/CPU includes native stepping, per-step validation and output. Compilation/offline reduction separate; no imposed wall/RSS budget.')
    resources['key_interface_check_wall_s']=read(O/'KEY_CHECKS.json')['wall_s']
    resources['historical_physics_reruns']=0
    dump(O/'RUN_RESOURCES.json',resources)
    print('ANALYZED',[(k,v['final']['detached_nodes'],v['final_metrics']['fragment_count_target_band'],v['exploratory_target_met']) for k,v in results.items()])
if __name__=='__main__':main()
