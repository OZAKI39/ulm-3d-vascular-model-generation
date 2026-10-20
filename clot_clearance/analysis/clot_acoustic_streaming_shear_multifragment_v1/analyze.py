"""Read-only reductions of accepted states; no topology or loading mutations."""
from common import *
from multifragment_objective import metrics,BINS,SPEC,shape
from shear_loading import ShearLoad,OFF,COHORT,WAVENUMBER
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
        cohort_time_s=SPEC['common_cohort_time_s'],upper_drive_stop_s=OFF,remaining_shear_fraction_at_cohort=ShearLoad(ref,name).peak(COHORT)/(summary['physical_template_peak_Pa']*summary['traction_multiplier']),
        secondary_events=splits,secondary_after_common_cohort=sum(q['time_s']>=SPEC['common_cohort_time_s'] for q in splits),
        crack_networks=dict(n_births=networks['n_births'],n_merges=networks['n_merge_events'],max_simultaneous=networks['max_simultaneous_networks'],max_birth_centroid_separation_m=separation,separated_birth_centroids_ge_horizon=separation>=HORIZON),
        inclusive_failed_material_networks={key:networks['all_material_including_interface'][key] for key in ['n_births','n_merge_events','max_simultaneous_networks','new_failed_bulk_bonds','new_failed_interface_bonds']},
        max_relative_energy_residual=summary['max_relative_energy_residual'],
        peak_shear_Pa=max(q['shear_peak_Pa'] for q in rows),peak_sum_nodal_force_N=max(q['sum_force_norm_N'] for q in rows),geometry_first_diagnostic=summary['frozen_geometry_first_diagnostic'],
        source_sha256={p.name:sha(p) for p in [case/'SUMMARY.json',case/'history.csv',case/'CONNECTED_COMPONENT_HISTORY.jsonl.gz',case/'COMPONENT_REGISTRY.json']},
        scope='Original surviving bulk/interface bonds and anchor reachability; contact excluded. Same frozen shear shape and duration, scaled stress. Decaying shear is not zero immediately after upper-wall stop.')
    dump(case/'ANALYSIS.json',result);dump(case/'CRACK_NETWORK_AUDIT.json',networks)
    dump(case/'LOCAL_DEFORMATION_AUDIT.json',diagnostics)
    dump(case/'FRAGMENT_HISTORY.json',dict(fragments=stable,candidates_ge14=[q for q in stable if q['nodes']>=14],secondary_events=splits,common_cohort_ids=cohort_ids))
    dump(case/'FRAGMENT_TRACE.json',trace)
    return result,summary,rows,networks

def main():
 start=time.perf_counter();ref=load_npz(O/'provenance/reference.npz');n=int(ref['clot_count']);results={};summaries={};hist={};networks={};requirements=[];distribution=[];after={};bridge_rows=[]
 for name in ['stokes','inverse']:
  results[name],summaries[name],hist[name],networks[name]=analyze_case(name,ref);q=results[name];s=summaries[name];end=q['final'];load=ShearLoad(ref,name)
  requirements.append(dict(case=name,specimen_id='original_1792_material_336_anchor_notched_3D_dx50um',label=load.label,evidence=s['mechanical_status'],peak_shear_Pa=q['peak_shear_Pa'],spatial_mode='sin(2*pi*y/700um)*ex on top reference triangles',area_m2=float(load.area.sum()),time_s=END,upper_drive_stop_s=OFF,peak_sum_nodal_force_N=q['peak_sum_nodal_force_N'],peak_abs_net_force_N=max(abs(r['net_force_x_N']) for r in hist[name]),peak_abs_reference_torque_Nm=max(abs(r['net_moment_z_Nm']) for r in hist[name]),absolute_impulse_Ns=end['absolute_impulse_Ns'],net_impulse_x_Ns=end['impulse_x_Ns'],reference_angular_impulse_z_Nms=end['moment_impulse_z_Nms'],W_J=end['W_J'],first_softening_s=q['first_softening_s'],first_failure_s=q['first_failure_s'],first_softening_shear_Pa=next((r['shear_peak_Pa'] for r in hist[name] if r['softened']>0),None),first_failure_shear_Pa=next((r['shear_peak_Pa'] for r in hist[name] if r['bulk_failed']+r['interface_failed']>0),None),first_target_fragment_s=next((r['time_s'] for r in hist[name] if r['N20_200']>0),None),n_fragments_20_200_final=end['N20_200'],n_fragments_20_200_history_max=q['history_max_N20_200'],n_fragments_34_200=end['n_fragments_34_200'],n_fragments_ge_114=end['n_fragments_ge_114'],mass_fraction_target_initial=end['mass_fraction_target_initial'],mass_fraction_target_detached=end['mass_fraction_target_detached'],fines_fraction_detached=end['debris12_detached_fraction'],largest_fraction_detached=end['largest_free_fragment_fraction_detached'],total_detached_fraction_initial=end['total_detached_mass_fraction'],target_met=q['exploratory_target_met'],limitation='Frozen one-way reference load; no bubble realization; material uncalibrated',max_relative_R=q['max_relative_energy_residual']))
  for stage in ['common_cohort','final']:
   state=load_npz(O/'cases'/name/'keys'/(stage+'.npz'));_,sizes,attached=graph.topology(ref,state['d']);free=sizes[~attached]
   for lo,hi in BINS:
    vals=free[(free>=lo)&((free<=hi) if hi else True)];total=int(vals.sum())
    distribution.append(dict(case=name,state=stage,time_s=float(state['time']),size_bin=f'{lo}-{hi}' if hi else '201+',component_count=len(vals),material_nodes=total,mass_kg=total*RHO*V,mass_fraction_initial=total/n,mass_fraction_detached=total/int(free.sum()) if free.sum() else None))
  after[name]={}
  for key in ['peak_load','final']:
   path=O/'cases'/name/'keys'/(key+'.npz');state=load_npz(path);items=[]
   for region in candidate_regions(ref):
    info,rows=inspect(ref,state,region);items.append(info);bridge_rows.extend(dict(case=name,state=key,**row) for row in rows)
   after[name][key]=dict(time_s=float(state['time']),source_sha256=sha(path),regions=items)
 csvout(O/'MULTIFRAGMENT_SHEAR_REQUIREMENTS.csv',requirements);csvout(O/'FRAGMENT_SIZE_DISTRIBUTION.csv',distribution);csvout(O/'SURVIVING_BRIDGE_BONDS.csv',bridge_rows)
 dump(O/'CRACK_BRIDGE_AUDIT.json',dict(cases=after,diagnostic_only=True,pinned_min_cut_calls=12))
 dump(O/'ANALYSIS.json',results);dump(O/'NATIVE_SHEAR_CASE_RESULTS.json',summaries);dump(O/'DAMAGE_NETWORK_HISTORY.json',networks)
 dump(O/'FRACTURE_EVENTS.json',{k:v['events'] for k,v in summaries.items()})
 dump(O/'CONNECTED_COMPONENT_HISTORY.json',dict(cases={k:dict(history=str((O/'cases'/k/'CONNECTED_COMPONENT_HISTORY.jsonl.gz').relative_to(O)),registry=read(O/'cases'/k/'COMPONENT_REGISTRY.json'),fragment_lifetimes=read(O/'cases'/k/'FRAGMENT_HISTORY.json'),records=len(hist[k])+1) for k in results},definition='Index and exact original membership registry; full every-accepted-state histories are in lossless JSONL gzip, not contact-based clustering.'))
 ledger={}
 for k,q in results.items():
  end=q['final'];ledger[k]={key:end[key] for key in ['W_J','U_J','K_J','D_J','Uc_J','Jc_J','R_J']};ledger[k]['max_relative_residual']=q['max_relative_energy_residual']
 dump(O/'ENERGY_LEDGER.json',dict(identity='W = Delta U_LPS + Delta K + D + Delta U_pair - J_pair_activation + R',work='Gauss integration of prescribed shear over each accepted displacement; actual native material motion, no damping',numerical_residual_not_dissipation=True,cases=ledger,normalization='running max(abs(W),abs(Delta U)+K+D+Uc), tolerance 1%',no_rigid_indenter=True))
 audit=read(O/'ACOUSTIC_AND_BUBBLE_INPUT_AUDIT.json');a=requirements[0];flow=ShearLoad(ref,'stokes');ts=np.linspace(0,END,14401);gradient=float(np.max(abs(flow.amplitude(ts)))*WAVENUMBER)
 physical=[dict(source='MEVO_F10S_APECSS_actual_files',evidence='ACTUAL_ACOUSTIC_AND_RADIAL_DATA_ONLY',frequency_Hz=audit['current_acoustics']['frequency_Hz'],acoustic_conditions=audit['current_acoustics']['model'],bubble_positions='Own F10S coordinates, unregistered to clot',radial_response='APECSS volume harmonics in valid adaptive leaves',geometry='F10S own channel; clot gap UNAVAILABLE',mu_Pa_s='UNAVAILABLE_FOR_REGISTERED_CLOT_MODEL',local_fluid_velocity_m_s='UNAVAILABLE',peak_shear_Pa='UNAVAILABLE',max_spatial_shear_gradient_Pa_m='UNAVAILABLE',duration_s='UNAVAILABLE_FOR_STREAMING',peak_sum_force_N='UNAVAILABLE',impulse_Ns='UNAVAILABLE',scope='STREAMING_LOAD_NOT_YET_DETERMINED'),dict(source='Analytic-mode unsteady viscous diffusion',evidence='PHYSICS_APPROXIMATION / PRESCRIBED_STOKES_SHEAR_TEMPLATE',frequency_Hz='UNAVAILABLE_NOT_ACOUSTIC_DRIVE',acoustic_conditions='No acoustic input converted',bubble_positions='UNAVAILABLE',radial_response='UNAVAILABLE_NOT_USED',geometry='25um gap, 700um transverse periodic pattern, stationary bottom, prescribed upper no-slip velocity',mu_Pa_s=.001,local_fluid_velocity_m_s=.004,peak_shear_Pa=a['peak_shear_Pa'],max_spatial_shear_gradient_Pa_m=gradient,duration_s=END,peak_sum_force_N=a['peak_sum_nodal_force_N'],impulse_Ns=a['absolute_impulse_Ns'],scope='Infinite-plane fluid solution; small specimen frozen one-way traction approximation')]
 csvout(O/'STREAMING_PHYSICAL_CAPABILITY.csv',physical)
 target=any(q['exploratory_target_met'] for q in results.values());observed=any(q['history_max_N20_200']>=2 for q in results.values());target_band_observed=any(q['history_max_N20_200']>0 for q in results.values())
 status=dict(software_status='COMPLETED',STREAMING_LOAD_COMPUTED=False,STREAMING_LOAD_NOT_YET_DETERMINED=True,SHEAR_REQUIREMENT_IDENTIFIED=target,MULTIFRAGMENT_NOT_OBSERVED=not observed,MULTIFRAGMENT_OBSERVED=observed,TARGET_BAND_FRAGMENT_OBSERVED=target_band_observed,five_fragment_target_met=target,ACOUSTIC_FEASIBILITY_NOT_ESTABLISHED=True,template_traction_ratio=requirements[1]['peak_shear_Pa']/requirements[0]['peak_shear_Pa'],ratio_scope='Same geometry/mode/time between prescribed templates only; not real bubble shortage nor a target-fragment threshold',cases={k:dict(status=v['status'],history_max_target_fragments=v['history_max_N20_200'],final_target_fragments=v['final']['N20_200'],detached_nodes=v['final']['detached_nodes'],history_max_free_nodes=v['history_max_free_nodes']) for k,v in results.items()},next_V3_admitted=False,new_optimizer_trials=0,native_library_unchanged=True)
 status['SCIENTIFIC_CONCLUSION']='Shear-to-native-fracture software chain restored. One continuous counter-shear mode was measured, not optimized. Bubble streaming is undetermined; any discrete target observations remain frozen-surface mechanical evidence, not acoustic feasibility or a best-mode ranking.'
 dump(O/'FINAL_STATUS.json',status)
 dump(O/'RUN_RESOURCES.json',dict(cases={k:{key:v[key] for key in ['wall_s','cpu_s','peak_RSS_bytes','accepted_steps','predictive_substeps','predictive_rejections']} for k,v in summaries.items()},native_wall_s=sum(v['wall_s'] for v in summaries.values()),native_cpu_s=sum(v['cpu_s'] for v in summaries.values()),peak_RSS_bytes=max(v['peak_RSS_bytes'] for v in summaries.values()),ranks=1,threads_per_case=1,new_main_native_cases=2,new_full_V3_cases=0,optimizer_trials=0,analysis_wall_s=time.perf_counter()-start))
 print('ANALYSIS_DONE',status)
if __name__=='__main__':main()

