"""Read-only analysis of the two frozen protocols, including failed completion."""
from common import *
from current_surface import CurrentSurface
from multifragment_shear_objective import metrics,BINS
import csv,gzip,time,resource
graph=module('graph_v2_analysis',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')

def read_history(path):
    with open(path) as f:
        rows=list(csv.DictReader(f))
    for r in rows:
        for k,v in r.items():
            if v=='':r[k]=None
            elif v in ['True','False']:r[k]=(v=='True')
            else:
                try:r[k]=float(v)
                except ValueError:pass
    return rows

def write_csv(path,rows):
    with open(path,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def damage_regions(case,ref):
    N=len(ref['X']);parent=np.arange(N);active=np.zeros(N,bool);births=[];merges=0;count=0;maximum=0
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    with gzip.open(case/'damage_events.jsonl.gz','rt') as f:
        for line in f:
            q=json.loads(line);ids=np.asarray(q['original_bond_ids'])-1
            ids=ids[np.asarray(q['d'])==1]
            for bid in ids:
                i,j=ref['bonds'][bid];ai,aj=active[i],active[j]
                if not ai and not aj:
                    count+=1;births.append(dict(time_s=q['time_s'],bond_id=int(bid+1),reference_center_m=ref['X'][[i,j]].mean(0)))
                elif ai and aj and find(i)!=find(j):count-=1;merges+=1
                active[i]=active[j]=True;parent[find(j)]=find(i)
            maximum=max(maximum,count)
    return dict(birth_count=len(births),merge_count=merges,max_concurrent_regions=maximum,final_regions=count,births=births,
        definition='Naturally fully failed original bonds connected by shared endpoints; not a reconstructed complete crack surface or free material block')

def main():
    started=time.perf_counter();cpu=time.process_time();ref=load_npz(O/'provenance/reference.npz');n=int(ref['clot_count']);N=len(ref['X']);b=ref['bonds'];cut=ref['initial_cut'];all_results={};surfaces={};energies={};secondary={};topologies={};workrows=[];bins=[];loading=[];validation={};resources={};work_diagnostics={}
    freeze=json.loads((O/'SOURCE_FREEZE.json').read_text())
    for p,h in freeze['implementation_sha256'].items():assert sha(R/p)==h,p
    for p,h in freeze['input_sha256'].items():assert sha(O/p)==h,p
    for mode in CASES:
        case=O/'cases'/mode;summary=json.loads((case/'SUMMARY.json').read_text());rows=read_history(case/'history.csv');last=rows[-1]
        initial=load_npz(case/'keys/initial.npz');endpath=case/'keys'/('final.npz' if summary['status']=='COMPLETED' else 'last_accepted.npz');end=load_npz(endpath)
        assert np.array_equal(initial['x'],ref['X']) and not initial['velocity'].any() and np.array_equal(initial['d'],cut.astype(float))
        labels,sizes,attached=graph.topology(ref,end['d']);free=sizes[~attached];objective=metrics(free,n)
        registry=json.loads((case/'COMPONENT_REGISTRY.json').read_text());history_records=0;live_ids=set();cohort_ids=set();last_components=[]
        with gzip.open(case/'CONNECTED_COMPONENT_HISTORY.jsonl.gz','rt') as f:
            for line in f:
                h=json.loads(line);history_records+=1
                ids=[c['component_id'] for c in h['components']]
                members=[int(a) for cid in ids for a in registry[str(cid)]['atom_ids']]
                assert len(members)==len(set(members))==h['detached_nodes'] and h['detached_nodes']+h['attached_nodes']==n
                if abs(h['time_s']-COHORT)<1e-14:cohort_ids={cid for cid in ids if 20<=registry[str(cid)]['nodes']<=200}
                live_ids=set(ids);last_components=h['components']
        assert history_records==len(rows)+1
        survivors=cohort_ids&live_ids
        objective['cohort_target_survivors_to_end']=len(survivors)
        objective['full_target']=bool(objective['instant_target'] and len(survivors)>=5 and summary['status']=='COMPLETED')
        bounds_error=0.;replay_error=0.;map_force=0.;map_moment=0.;map_power=0.;snapshots=[]
        surf=CurrentSurface(ref,'checkerboard')
        initial_mask=np.sum(ref['initial_cut'][surf.te]*(1<<np.arange(6)),axis=1)
        initially_disconnected=surf.disconnected[initial_mask]
        precut_source=~surf.original&(initially_disconnected[surf.owners[:,0]]|initially_disconnected[np.maximum(surf.owners[:,1],0)])
        for p in sorted((case/'keys').glob('*.npz')):
            if p.stem=='failed_unaccepted':continue
            s=load_npz(p);assert s['x'].dtype==np.float64 and s['d'].dtype==np.float64
            assert np.array_equal(s['atom_ids'],np.arange(1,N+1))
            assert np.array_equal(s['x'][n:],ref['X'][n:]) and not s['velocity'][n:].any()
            assert np.isfinite(s['x']).all() and np.isfinite(s['force']).all() and (s['d']>=0).all() and (s['d']<=1).all()
            assert abs(float(s['field_amplitude_Pa'])-surf.field.amplitude(float(s['time'])))<1e-8
            z=surf.evaluate(s['x'],s['velocity'],s['d'],float(s['time']))
            replay_error=max(replay_error,float(abs(z['force']-s['external_force']).max()))
            assert np.array_equal(z['active'],s['surface_active'])
            assert not np.any(z['original']&z['accessible'])
            assert not z['face_force'][~z['active']].any()
            map_force=max(map_force,float(abs(z['force'].sum(0)-z['face_force'].sum(0)).max()))
            map_moment=max(map_moment,float(abs(np.cross(s['x'],z['force']).sum(0)-np.cross(z['centroid'],z['face_force']).sum(0)).max()))
            map_power=max(map_power,abs(float(np.sum(z['force']*s['velocity']))-float(np.sum(z['face_force']*z['face_velocity']))))
            lab,sz,att=graph.topology(ref,s['d']);nodefree=np.r_[~att[lab],np.zeros(N-n,bool)]
            freefaces=z['active']&nodefree[surf.tri[:,0]]
            supported_labels=np.unique(lab[surf.tri[freefaces].ravel()])
            unknown_labels=np.setdiff1d(np.flatnonzero(~att),supported_labels)
            candidate_area=float(z['area'][z['candidate']].sum());unknown_area=float(z['area'][z['unresolved']].sum())
            snapshots.append(dict(key=p.stem,time_s=float(s['time']),original_faces=int(z['original'].sum()),new_open_faces=int(z['new_open'].sum()),accessible_new_faces=int(z['accessible'].sum()),accessible_reopened_precut_faces=int((z['accessible']&precut_source).sum()),accessible_natural_new_fracture_faces=int((z['accessible']&~precut_source).sum()),new_area_m2=float(z['area'][z['accessible']].sum()),closed_unresolved_faces=int(z['unresolved'].sum()),unresolved_candidate_area_m2=unknown_area,unresolved_candidate_area_fraction=unknown_area/candidate_area if candidate_area else None,
                occluded_faces=int(z['occluded'].sum()),true_free_components=int((~att).sum()),free_components_with_surface=len(supported_labels),free_components_without_surface=len(unknown_labels),free_nodes_in_components_without_surface=int(sz[unknown_labels].sum()),free_surface_faces=int(freefaces.sum()),
                accessible_face_ids=(np.flatnonzero(z['accessible'])+1),accessible_aperture_m=z['aperture_m'][z['accessible']],fluid_region_ids=z['fluid_region'][z['accessible']],unresolved_face_ids=np.flatnonzero(z['unresolved'])+1))
        assert replay_error<1e-15 and map_force<1e-16 and map_moment<1e-19 and map_power<1e-16
        scal=np.array([[r[k] for k in ['W_J','W_drive_J','W_drag_J','U_J','K_J','D_J','Uc_J','Jc_J','R_J','relative_R']] for r in rows])
        minimum_D_increment=float(np.min(np.diff(np.r_[0.,scal[:,5]])));assert minimum_D_increment>=-1e-20
        decomp=float(abs(scal[:,0]-scal[:,1]-scal[:,2]).max())
        identity=float(abs(scal[:,0]-(scal[:,3]-initial['scalars'][0]+scal[:,4]-initial['K']+scal[:,5]+scal[:,6]-initial['Uc']-scal[:,7])-scal[:,8]).max())
        assert decomp<1e-20 and identity<1e-20 and scal[:,9].max()<=.01
        assert all(abs(r['drive_power_W'])<1e-20 for r in rows if r['time_s']>=OFF+1e-14)
        assert max(r['max_dd'] for r in rows)<=.1+1e-10
        # Damage stream plus original IDs reconstructs every accepted damage state.
        damage=cut.astype(float);changed_events=0
        with gzip.open(case/'damage_events.jsonl.gz','rt') as f:
            for line in f:
                q=json.loads(line);ids=np.asarray(q['original_bond_ids'])-1;value=np.array(q['d']);assert (value>=damage[ids]).all();damage[ids]=value;changed_events+=1
        assert np.array_equal(damage,end['d'])
        final_surface=next(q for q in snapshots if q['key']==endpath.stem)
        exposure=load_npz(case/'accepted_surface_history.npz')['ever_loaded_faces']
        surface_result=dict(status='PROXY_IMPLEMENTED_WITH_EXPLICIT_UNRESOLVED_GEOMETRY',snapshots=snapshots,final=final_surface,
            unique_ever_loaded_new_faces=int((exposure&~surf.original).sum()),unique_ever_loaded_natural_new_fracture_faces=int((exposure&~surf.original&~precut_source).sum()),unique_ever_loaded_reopened_precut_faces=int((exposure&precut_source).sum()),exposure_channel_note='New interior exposure includes reopening of existing precut and naturally formed crack faces; distinguish counts. Exact recorded new-surface work is their combined channel, not solely newly fractured material.',initial_area_m2=last['initial_loaded_area_m2'],cumulative_unique_area_at_first_loading_m2=last['cumulative_unique_loaded_area_m2'],max_accessible_new_faces=int(max(r['fluid_accessible_new_faces'] for r in rows)),max_new_area_m2=max(r['new_accessible_area_m2'] for r in rows),
            true_unrepresented_continuum_area_fraction=None,true_area_status='UNDETERMINED_FROM_POINT_DATA; candidate-area and component-coverage metrics are not true continuum surface completeness',mapping_duplicate_triangles=0,
            contribution=dict(actual_net_work_on_new_accessible_faces_J=last['W_accessible_new_faces_J'],causal_claim='No same-model wetting-off control was run, so this work channel is direct evidence of load contribution, not isolated causal proof of improved breakup'))
        workrow=dict(case=mode,status='SAME_FROZEN_OPEN_LOOP_DRIVE',actual_net_J=last['W_J'],drive_J=last['W_drive_J'],drag_J=last['W_drag_J'],actual_peak_shear_Pa=max(r['shear_peak_Pa'] for r in rows),max_surface_speed_m_s=max(r['max_surface_speed_m_s'] for r in rows),max_outer_speed_m_s=max(r['max_outer_speed_m_s'] for r in rows),max_Re_h=max(r['Re_h_outer_max'] for r in rows),initial_area_m2=last['initial_loaded_area_m2'],cumulative_unique_area_m2=last['cumulative_unique_loaded_area_m2'],drive_duration_s=OFF,protocol_time_s=summary['time_s'])
        death={int(e['parent_component_id']):e['time_s'] for e in summary['secondary_events']}
        lifetime=[dict(component_id=int(cid),nodes=q['nodes'],birth_s=q['birth_time_s'],last_s=death.get(int(cid),summary['time_s']),lifetime_s=death.get(int(cid),summary['time_s'])-q['birth_time_s'],right_censored=int(cid) not in death,atom_ids=q['atom_ids']) for cid,q in registry.items()]
        topologies[mode]=dict(final_components=last_components,final_metrics=objective,fragment_lifetimes=lifetime,registry_file=str((case/'COMPONENT_REGISTRY.json').relative_to(O)),accepted_history_file=str((case/'CONNECTED_COMPONENT_HISTORY.jsonl.gz').relative_to(O)),historical_max_target_count=summary['max_medium_components'],historical_largest_component=summary['history_max_free_nodes'],definition='Original material IDs connected by currently live bulk/interface bonds; anchors determine attachment. Contact and proximity never form material edges.',nonlocal_connections='Any region still in an attached or larger component retains live material paths; no min-cut/bridge optimization performed')
        secondary[mode]=dict(count=len(summary['secondary_events']),events=summary['secondary_events'])
        energies[mode]=dict(minimum_damage_energy_increment_J=minimum_D_increment,identity='W_drive+W_drag = DeltaU+DeltaK+D+DeltaUc-Jc+R',last=last,max_relative_R=float(scal[:,9].max()),max_decomposition_error_J=decomp,max_identity_error_J=identity,work_scheme='Accepted endpoint force trapezoid dot actual native drift; no nominal area estimate, curve rescaling or absorption of R into D',full_history=str((case/'history.csv').relative_to(O)))
        resources[mode]={k:summary[k] for k in ['wall_s','cpu_s','peak_RSS_bytes','accepted_steps','time_s']}
        validation[mode]=dict(status='PASS',accepted_rows=len(rows),component_records=history_records,damage_events=changed_events,FP64_keys=len(snapshots),force_replay_error_N=replay_error,resultant_error_N=map_force,moment_error_Nm=map_moment,virtual_power_error_W=map_power,work_channel_error_J=decomp,energy_identity_error_J=identity,initial_state_matches_frozen=True,original_ID_and_anchor_checks=True,native_irreversibility_and_force_checks=summary['max_independent_errors'])
        for lo,hi in BINS:
            selected=free[(free>=lo)&((free<=hi) if hi else True)]
            bins.append(dict(case=mode,bin=f'{lo}-{hi}' if hi else '>200',count=len(selected),nodes=int(selected.sum()),reference_volume_m3=float(selected.sum()*V),mass_kg=float(selected.sum()*V*RHO),fraction_initial=float(selected.sum()/n),fraction_detached=float(selected.sum()/free.sum()) if free.sum() else None))
        all_results[mode]=dict(native=summary,work=workrow,objective=objective,final_free_sizes=sorted(free.tolist(),reverse=True),surface=surface_result,damage_regions=damage_regions(case,ref),secondary_count=len(summary['secondary_events']))
        surfaces[mode]=surface_result;workrows.append(workrow)
        loading.extend([dict(case=mode,**{k:r[k] for k in ['step','time_s','nominal_amplitude_Pa','shear_peak_Pa','max_node_force_N','max_surface_speed_m_s','max_outer_speed_m_s','Re_h_outer_max','active_current_area_m2','effective_loaded_area_m2','fluid_accessible_new_faces','new_accessible_area_m2','closed_or_unresolved_faces','occluded_faces','facet_shear_impulse_Ns','W_J','W_drive_J','W_drag_J','W_accessible_new_faces_J','drive_power_W','drag_power_W','phase','dt_s']}) for r in rows])
    write_csv(O/'ACTUAL_WORK_COMPARISON.csv',workrows);write_csv(O/'FRAGMENT_SIZE_DISTRIBUTION.csv',bins);write_csv(O/'ACTUAL_LOADING_HISTORY.csv',loading)
    dump(O/'UNIFORM_NATIVE_RESULT.json',all_results[CASES[0]]);dump(O/'HETEROGENEOUS_NATIVE_RESULT.json',all_results[CASES[1]])
    dump(O/'FLUID_ACCESSIBLE_CRACK_FACES.json',surfaces);dump(O/'TRUE_FRAGMENT_HISTORY.json',topologies);dump(O/'SECONDARY_FRAGMENTATION.json',secondary);dump(O/'ENERGY_LEDGER.json',energies);dump(O/'DELIVERY_VALIDATION.json',dict(status='PASS',cases=validation))
    resources['new_checks']=json.loads((O/'VARIATIONAL_ENERGY_AUDIT.json').read_text())
    resources['notes']='Two new full native protocols, same frozen drive, no outcome-based retry or search. Per-process CPU/wall/RSS; parallel elapsed time is not sum of case walls.'
    resources['analysis']=dict(wall_s=time.perf_counter()-started,cpu_s=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    dump(O/'RUN_RESOURCES.json',resources)
    from spatial_analysis import analyze_spatial
    spatial=analyze_spatial(ref,all_results,topologies)
    write_csv(O/'FRAGMENT_MASS_DISTRIBUTION.csv',bins)
    dump(O/'COMPARISON.json',dict(cases=all_results,spatial=spatial))
    print('ANALYSIS_COMPLETE',[(m,q['work']['actual_net_J'],q['final_free_sizes']) for m,q in all_results.items()],flush=True)

if __name__=='__main__':main()
