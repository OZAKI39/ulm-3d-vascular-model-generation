"""Joint orientation/contact/topology planning in the existing s1-8 workflow."""
from __future__ import annotations
import copy,hashlib,json,shutil,time,traceback
from pathlib import Path
import numpy as np
import trimesh
import yaml
from .support_inputs import ROOT,save,table,verify,transform,digest
from .support_v2_inputs import load,oriented,public
from .support_v2_orientation import coarse,diverse_top,fits,rigid
from .support_custom_detection import fused_detect,intervals
from .support_v2_layers import analyse,coverage
from .support_v2_planner import candidate_roofs,select_roofs,scaffold,self_buildability
from .support_custom_qc import canonical_stl,intersection_volume

def plan_candidate(out,row,inputs,profile,cfg,mold,vascular,box,branches,topology,reuse_scaffold=False):
    folder=out/'candidates'/row['orientation_id'];folder.mkdir(parents=True,exist_ok=True);start=time.monotonic();m=np.array(row['print_transform']);restore=np.linalg.inv(m)
    model=mold.copy();model.apply_transform(m);v=vascular.copy();v.apply_transform(m);b=box.copy();b.apply_transform(m)
    bs,top=oriented(branches,topology,m);corners=np.array(np.meshgrid(*np.array(inputs['inner_bounds']).T)).reshape(3,-1).T;corners=transform(corners,m);inner=np.array([corners.min(0),corners.max(0)])
    det=fused_detect(v,b,bs,top,inner,profile,cfg,folder)
    field=analyse(v,profile,'vascular',b,folder);boxfield=analyse(b,profile,'box',None,folder)
    print('V2 DETECT',row['orientation_id'],len(field['rows']),'vascular layer regions',len(boxfield['rows']),'box regions',flush=True)
    roofs=candidate_roofs(det,field,boxfield,v,b,bs,profile,profile['roof_gap_mm'],top);selected,planning=select_roofs(roofs,[field,boxfield],profile)
    print('V2 CONTACTS',row['orientation_id'],len(roofs),'raw',len(selected),'selected',flush=True)
    save(folder/'roof_candidates_pre_graph.json',public(roofs));save(folder/'selected_roofs_pre_graph.json',public(selected))
    cache=folder/'scaffold_local_cache.pkl'
    if reuse_scaffold and cache.exists():
        import pickle
        from .support_v2_mesh import support_stl
        with cache.open('rb') as f:graph,connected,pieces,approved=pickle.load(f)
        support=support_stl(trimesh.boolean.union(pieces,engine='manifold'))
        graph['numerical_support_cleanup']=support.metadata.get('support_only_numerical_cleanup_tolerance_mm',0.)
        graph['export_rebuilt_from_saved_primitives']=True
    else:support,graph,connected,pieces,approved=scaffold(selected,v,b,profile,restore,np.array(inputs['inner_bounds']),folder)
    for roof in connected:
        if roof['geometry_class']!='box' and len(top.get('bifurcations',[])):
            roof['distance_to_bifurcation']=float(np.linalg.norm(np.array(top['bifurcations'])-roof['xyz'],axis=1).min())
            roof['bifurcation_distance_method']='Euclidean distance to transformed SWC bifurcation node'
    table(folder/'roof_candidates.csv',public(roofs));table(folder/'selected_roofs.csv',public(connected));save(folder/'roofs.json',public(connected));save(folder/'support_graph_v2.json',graph)
    table(folder/'support_primitives.csv',graph['primitives']);table(folder/'support_roots.csv',graph['roots']);table(folder/'merge_attempts.csv',graph['merge_attempts']);table(folder/'rejected_roofs.csv',graph['rejected'])
    # New critical intervals use the actual FDM requirement classification.
    from shapely.geometry import Point
    from shapely.strtree import STRtree
    required_tree=STRtree(field['polygons'])
    for p in det['raw']:
        point=Point(p['xyz'][:2]);near=required_tree.query(point.buffer(profile['line_width']),predicate='intersects')
        p['risk_level']='CRITICAL' if any(field['rows'][i]['classification']=='SUPPORT_REQUIRED' and abs(field['rows'][i]['z_mm']-p['xyz'][2])<max(.6,p['local_radius']) for i in near) else 'ATTENTION'
    critical=intervals(det['raw'],cfg);table(folder/'critical_intervals.csv',critical)
    cov=coverage(field,connected,profile);bcov=coverage(boxfield,connected,profile);table(folder/'fdm_coverage.csv',cov['rows']);table(folder/'box_fdm_coverage.csv',bcov['rows'])
    model.export(folder/'mold_without_support.stl');v.export(folder/'vascular_with_ports.stl');b.export(folder/'box.stl')
    result=dict(candidate_id=row['orientation_id'],folder=str(folder),orientation=row,architecture='FDM_NATIVE_VASCULAR_SUPPORT_V2',scale_factor=1.,
        critical_intervals=len(critical),raw_roof_count=len(roofs),selected_roof_count=len(connected),planning=planning,vascular_coverage={k:val for k,val in cov.items() if k!='rows'},box_coverage={k:val for k,val in bcov.items() if k!='rows'},
        merge_count=graph['merge_count'],root_count=graph['root_count'],rail_count=graph['rail_count'],bridge_count=graph['bridge_count'],topology_warning=graph['topology_warning'],
        primitives_by_kind={k:sum(p['kind']==k for p in graph['primitives']) for k in {p['kind'] for p in graph['primitives']}},slice={},removal={},physical_gap_status=profile['gap_status'],bridge_status=profile['bridge_calibration_status'])
    if len(support.faces):
        support.export(folder/'custom_support.stl');assembly=canonical_stl(trimesh.boolean.union([model,support],engine='manifold'));assembly.export(folder/'mold_with_custom_support.stl')
        result['semantic_part_export_audit']=[]
        for name,kinds in [('support_roofs.stl',{'ROOF_PAD'}),('support_trunks.stl',{'VERTICAL_PILLAR','SELF_SUPPORTING_BRANCH','LONGITUDINAL_RAIL','SHORT_BRIDGE'}),('support_breakaway_parts.stl',{'BREAKAWAY_FOOT','BREAKAWAY_SEGMENT_NECK','BUTTRESS'})]:
            shapes=[p for p,r in zip(pieces,graph['primitives']) if r['kind'] in kinds]
            if shapes:
                from .support_v2_mesh import support_stl
                raw=trimesh.boolean.union(shapes,engine='manifold')
                try:
                    support_stl(raw).export(folder/name);result['semantic_part_export_audit'].append(dict(file=name,closed_positive_volume=True))
                except ValueError as error:
                    # Preserve a reviewable failure artifact rather than lose
                    # the complete, closed support/assembly and its checks.
                    raw.export(folder/name);raw.export(folder/(name+'.source.ply'))
                    result['semantic_part_export_audit'].append(dict(file=name,closed_positive_volume=False,error=str(error),purpose='DIAGNOSTIC_ONLY_NOT_SEPARATE_PRINT_PART'))
        allowed=trimesh.boolean.union(approved,engine='manifold') if approved else trimesh.Trimesh()
        overlap=trimesh.boolean.intersection([support,b],engine='manifold');excess=trimesh.boolean.difference([overlap,allowed],engine='manifold') if len(overlap.faces) and len(allowed.faces) else overlap
        audit=dict(vascular_intersection_mm3=intersection_volume(support,v),illegal_box_intersection_mm3=abs(float(excess.volume)) if len(excess.faces) else 0.,frozen_model_removed_mm3=intersection_volume(model,model)-intersection_volume(model,assembly))
        audit['passed']=audit['vascular_intersection_mm3']<1e-6 and audit['illegal_box_intersection_mm3']<1e-5 and abs(audit['frozen_model_removed_mm3'])<1e-4
        result.update(support_volume_mm3=float(support.volume),union_audit=audit,self_buildability=self_buildability(support,b,profile),build_volume_fit=fits(assembly.vertices,np.eye(4),profile['build_bounds_mm']))
        heights=[p['bbox_mm'][1][2]-p['bbox_mm'][0][2] for p in graph['primitives'] if p['kind']=='VERTICAL_PILLAR'];widths=[p['width_mm'] for p in graph['primitives'] if 'width_mm' in p and p['kind']!='ROOF_PAD']
        result['mechanical_proxy']=dict(max_pillar_height_mm=max(heights,default=0.),min_width_mm=min(widths,default=0.),max_slenderness=max([height/max(width,1e-9) for height,width in zip(heights,[p['width_mm'] for p in graph['primitives'] if p['kind']=='VERTICAL_PILLAR'])],default=0.),bridge_count=graph['bridge_count'],bridge_max_span_mm=max([p.get('span_mm',0) for p in graph['primitives']],default=0.),supported_load_center=np.mean([r['xyz'] for r in connected],axis=0).tolist(),physical_strength_validated=False)
    else:result.update(support_volume_mm3=0.,union_audit={'passed':False},self_buildability={'passed':False},build_volume_fit=row['build_volume_fit'],error='NO_LEGAL_FDM_SCAFFOLD')
    result['planning_runtime_seconds']=time.monotonic()-start;save(folder/'candidate.json',result)
    print('V2 PLAN',row['orientation_id'],'roofs',len(connected),'merges',graph['merge_count'],'coverage',round(cov['percent'],2),'self',result['self_buildability']['passed'],flush=True)
    return result

def gates(r):
    f=[]
    if any(not row['closed_positive_volume'] for row in r.get('semantic_part_export_audit',[])):f.append('SEMANTIC_PART_STL_PRECISION_FAILURE')
    for key,label in [('vascular_coverage','UNSUPPORTED_VASCULAR_OR_PORT_REGIONS'),('box_coverage','UNSUPPORTED_BOX_REGIONS')]:
        if not r.get(key,{}).get('all_required_covered'):f.append(label)
    if not r.get('self_buildability',{}).get('passed'):f.append('SUPPORT_REQUIRES_SUPPORT')
    if not r.get('union_audit',{}).get('passed'):f.append('ILLEGAL_GEOMETRY_INTERSECTION')
    if not r.get('build_volume_fit'):f.append('BUILD_VOLUME_EXCEEDED')
    if r.get('removal',{}).get('trapped_pieces')!=0:f.append('TRAPPED_OR_UNASSESSED_REMOVAL_SEGMENTS')
    if r.get('removal',{}).get('break_design_audit',{}).get('passed') is False:f.append('BREAK_SECTION_CROSSES_UNWEAKENED_SUPPORT')
    if r.get('removal',{}).get('precision_failed_fragments',0):f.append('FRAGMENT_STL_PRECISION_FAILURE')
    if r.get('removal',{}).get('opening_exit_audit',{}).get('passed') is False:f.append('REMOVAL_PATH_DOES_NOT_FULLY_EXIT_CASTING_OPENING')
    s=r.get('slice',{})
    if not s.get('slicer_success'):f.append('SLICER_VALIDATION_PENDING')
    elif not (s.get('auto_support_off') and s.get('all_features_preserved') and s.get('all_roots_have_start_paths') and s.get('vascular_toolpaths_preserved')):f.append('SLICER_FEATURE_OR_ROOT_NOT_PRESERVED')
    if s.get('slicer_success'):
        if not s.get('actual_vascular_coverage',{}).get('all_required_covered') or not s.get('actual_box_coverage',{}).get('all_required_covered'):f.append('ACTUAL_ROOF_COVERAGE_INCOMPLETE')
        if s.get('native_support_toolpaths',-1)!=0:f.append('UNEXPECTED_NATIVE_SUPPORT')
        if not s.get('support_only_slicer_success'):f.append('SUPPORT_ONLY_SLICE_FAILED')
        if not s.get('bridge_toolpaths',{}).get('passed'):f.append('BRIDGE_TOOLPATH_NOT_CONFIRMED')
    if s.get('unresolved_warning_count',0):f.append('UNRESOLVED_FLOATING_WARNING')
    return f

def rank(r):
    failed=gates(r);geometric=[f for f in failed if f!='SLICER_VALIDATION_PENDING' and f!='TRAPPED_OR_UNASSESSED_REMOVAL_SEGMENTS']
    vc=r.get('slice',{}).get('actual_vascular_coverage',r.get('vascular_coverage',{}));bc=r.get('slice',{}).get('actual_box_coverage',r.get('box_coverage',{}))
    # Failed candidates are diagnostic only. Counting failed gates rewards a
    # nearly empty scaffold for being self-printable; use coverage to order
    # failed but geometrically legal candidates instead.
    return (bool(geometric),not r.get('union_audit',{}).get('passed') or not r.get('build_volume_fit'),-vc.get('percent',0),-bc.get('percent',0),not r.get('self_buildability',{}).get('passed'),bool(r.get('topology_warning')),
        r.get('slice',{}).get('contact_area_mm2',r.get('near_vessel_roof_area_mm2',1e9)),r.get('mechanical_proxy',{}).get('max_slenderness',1e9),r.get('removal',{}).get('trapped_pieces') if r.get('removal',{}).get('trapped_pieces') is not None else 1e9,r.get('support_volume_mm3',1e9),r.get('slice',{}).get('estimated_print_time_seconds') or 1e12,r['candidate_id'])

def removal_worker(record,inputs):
    from .support_v2_removal import run_removal
    try:return run_removal(record,inputs)
    except Exception as exc:
        save(Path(record['folder'])/'removal_failure.json',dict(error=str(exc),traceback=traceback.format_exc()))
        return dict(trapped_pieces=None,error=str(exc),status='REMOVAL_VALIDATION_FAILED')

def complete_contact_semantics(record,topology):
    """Fill old cached contact metadata without pretending to re-plan it."""
    points=topology.get('bifurcations',[])
    if not len(points):return
    from scipy.spatial import cKDTree
    tree=cKDTree(transform(points,np.asarray(record['orientation']['print_transform'])))
    folder=Path(record['folder']);changes=[]
    for name,csv_name in [('roofs.json','selected_roofs.csv'),('roof_candidates_pre_graph.json','roof_candidates.csv'),('selected_roofs_pre_graph.json',None)]:
        path=folder/name
        if not path.exists():continue
        rows=json.loads(path.read_text());count=0
        for row in rows:
            if row['geometry_class']!='box' and row.get('distance_to_bifurcation') is None:
                row['distance_to_bifurcation']=float(tree.query(row['xyz'])[0])
                row['bifurcation_distance_method']='Euclidean distance to transformed SWC bifurcation node'
                row['distance_metadata_backfilled_after_planning']=True;count+=1
        if count:
            save(path,rows)
            if csv_name:table(folder/csv_name,rows)
            changes.append(dict(file=name,filled_rows=count))
    audit=folder/'contact_semantic_backfill.json'
    if changes:save(audit,dict(changes=changes,geometry_changed=False,contact_selection_reoptimized=False,original_missing_values_were_not_used_as_bifurcation_penalties=True))
    if audit.exists():record['contact_semantic_backfill']=json.loads(audit.read_text())

def run(source=None,output=None,stage='all',config=None):
    from .support_inputs import discover
    inputs=discover(source=source);out=Path(output).resolve() if output else inputs['root'].parent/'vascular_support_planner/v2';out.mkdir(parents=True,exist_ok=True)
    prior=None
    if stage=='resume':
        names=['source_discovery.json','profile_runtime.json','resolved_config.json']
        if not all((out/n).exists() for n in names):raise ValueError('V2_RESUME_INPUT_MANIFEST_MISSING')
        previous=[json.loads((out/n).read_text()) for n in names]
        prior=hashlib.sha256(json.dumps(dict(inputs=previous[0]['hashes'],profile=previous[1],config=previous[2]),sort_keys=True,default=str).encode()).hexdigest()
    start=time.monotonic();inputs,profile,discovery,mold,v,box,branches,topology=load(out,source)
    legacy=out.parent/'legacy_explicit_support_v1';before=json.loads((legacy/'protected_before.json').read_text()) if (legacy/'protected_before.json').exists() else inputs['hashes'];verify(before);save(out/'protected_before.json',before)
    if stage=='all' and (legacy/'vascular_support_report.md').exists():
        prior_replay=out/'v1_merge_replay.json'
        current=legacy/'vascular_support_report.md'
        if not prior_replay.exists() or json.loads(prior_replay.read_text()).get('historical_report_sha256')!=digest(current):
            from .support_v2_postmortem import run as postmortem
            postmortem(out,inputs)
    cfg=yaml.safe_load(Path(config or ROOT/'config/vascular_support_planner.yaml').read_text());save(out/'resolved_config.json',cfg)
    key=hashlib.sha256(json.dumps(dict(inputs=inputs['hashes'],profile=profile,config=cfg),sort_keys=True,default=str).encode()).hexdigest()
    if prior is not None and key!=prior:raise ValueError('V2_RESUME_INPUTS_CHANGED_RUN_PLAN')
    if stage in ['all','orientation'] or not (out/'top_orientations.json').exists():
        rows=coarse(mold,v,box,branches,profile);table(out/'orientation_candidates.csv',rows);save(out/'orientation_candidates.json',rows);save(out/'top_orientations.json',diverse_top(rows))
    orientations=json.loads((out/'top_orientations.json').read_text())
    orientations.sort(key=lambda row:(not np.allclose(row['gravity_in_casting'],[0,0,1]),row['coarse_rank_value'],row['orientation_id']))
    if stage=='orientation':return dict(output=str(out),status='NEEDS_ADJUSTMENT')
    if stage in ['all','plan','geometry','research','resume']:
        results=[]
        from concurrent.futures import ProcessPoolExecutor,as_completed
        import multiprocessing
        # Avoid inheriting native VTK/FCL/TBB worker state through fork.
        pending=[]
        for row in orientations:
            path=out/'candidates'/row['orientation_id']/'candidate.json'
            if stage=='resume' and path.exists():
                r=json.loads(path.read_text())
                if r.get('self_buildability') and r.get('orientation')==row:
                    results.append(r);continue
            pending.append(row)
        with ProcessPoolExecutor(max_workers=3,mp_context=multiprocessing.get_context('spawn')) as pool:
            futures={pool.submit(plan_candidate,out,row,inputs,profile,cfg,mold,v,box,branches,topology,stage=='resume'):row for row in pending}
            for future in as_completed(futures):
                row=futures[future]
                try:r=future.result()
                except Exception as exc:
                    folder=out/'candidates'/row['orientation_id'];folder.mkdir(parents=True,exist_ok=True);r=dict(candidate_id=row['orientation_id'],folder=str(folder),orientation=row,error=str(exc),rejection_reason='PLANNING_EXCEPTION');save(folder/'failure.json',dict(error=str(exc),traceback=traceback.format_exc()));print('V2 FAILURE',r['candidate_id'],str(exc),flush=True)
                results.append(r);results.sort(key=lambda r:r['candidate_id']);save(out/'planning_results.json',results)
        save(out/'planning_cache.json',dict(key=key))
    else:
        if json.loads((out/'planning_cache.json').read_text())['key']!=key:raise ValueError('V2_CACHE_MISMATCH_RUN_PLAN')
        results=json.loads((out/'planning_results.json').read_text())
    # Planning failures remain in the results. Select four actual geometry
    # candidates for final machine slicing, never a native/V1 benchmark.
    from .support_v2_sensitivity import add_metrics
    for r in results:
        if r.get('self_buildability') and (Path(r['folder'])/'custom_support.stl').exists():
            complete_contact_semantics(r,topology);add_metrics(r)
    finalists=sorted([r for r in results if r.get('self_buildability') and (Path(r['folder'])/'mold_with_custom_support.stl').exists()],key=lambda r:rank({**r,'slice':{},'removal':{}}))[:4]
    frozen_top4=out/'top4_candidate_ids.json'
    if stage not in ['all','plan','geometry','research','resume'] and frozen_top4.exists():
        identifiers=json.loads(frozen_top4.read_text());finalists=[next(r for r in results if r['candidate_id']==identity) for identity in identifiers]
    else:save(frozen_top4,[r['candidate_id'] for r in finalists])
    if stage in ['all','refine']:
        from .support_v2_refinement import refine
        for index,r in enumerate(finalists):
            try:new=refine(r,inputs,profile,cfg,mold,v,box,branches,topology,out)
            except Exception as exc:
                r['refinement']=dict(error=str(exc),status='REFINEMENT_FAILED');save(Path(r['folder'])/'refinement_failure.json',dict(error=str(exc),traceback=traceback.format_exc()));continue
            results[results.index(r)]=new;finalists[index]=new;add_metrics(new)
        save(out/'planning_results.json',results)
    if stage in ['all','removal']:
        from concurrent.futures import ProcessPoolExecutor,as_completed
        import multiprocessing
        with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn')) as pool:
            pending={pool.submit(removal_worker,r,inputs):r for r in finalists}
            for future in as_completed(pending):
                r=pending[future];r['removal']=future.result();save(Path(r['folder'])/'candidate.json',r);save(out/'planning_results.json',results)
    if stage in ['all','slice']:
        from .support_v2_validation import slice_candidate
        from .support_v2_sensitivity import run as sensitivity
        for index,r in enumerate(finalists):
            try:r['sensitivity']=sensitivity(r,profile)
            except Exception as exc:
                r['sensitivity']=dict(error=str(exc),status='SENSITIVITY_FAILED');save(Path(r['folder'])/'sensitivity_failure.json',dict(error=str(exc),traceback=traceback.format_exc()))
            try:r['slice']=slice_candidate(r,profile,discovery,secondary=index==0)
            except Exception as exc:
                native=Path(r['folder'])/'bambu/slice_result.json';observed=json.loads(native.read_text()) if native.exists() else {}
                r['slice']=dict(slicer_success=observed.get('slicer_success',False),auto_support_off=None,error=str(exc),status='SLICER_AUDIT_FAILED');save(Path(r['folder'])/'slice_failure.json',dict(error=str(exc),traceback=traceback.format_exc()))
            save(Path(r['folder'])/'candidate.json',r);save(out/'planning_results.json',results)
    if stage in ['all','coupon']:
        from .support_v2_coupon import run_coupon
        run_coupon(out/'coupon',profile,discovery)
    # Removal retries may have produced a newer independently saved diagnosis.
    # Merge those completed records instead of retaining a stale failure stub.
    for r in results:
        path=Path(r['folder'])/'support_removal_plan.json'
        if path.exists():
            r['removal']=json.loads(path.read_text())
            if 'opening_exit_audit' not in r['removal']:
                from .support_v2_removal import opening_exit_audit
                r['removal']['opening_exit_audit']=opening_exit_audit(r['folder'],r['removal'],inputs['inner_bounds'],box.bounds)
                save(path,r['removal'])
    for r in results:
        r['hard_gate_failures']=gates(r);save(Path(r['folder'])/'candidate.json',r)
    viable=[r for r in finalists if not r['hard_gate_failures']];best=min(viable or finalists,key=rank) if finalists else None
    summary=dict(architecture='FDM_NATIVE_VASCULAR_SUPPORT_V2',status='READY_FOR_HUMAN_REVIEW_AND_FDM_COUPON' if viable else 'NEEDS_ADJUSTMENT',production_selected=best['candidate_id'] if viable else None,
        diagnostic_candidate=best['candidate_id'] if best else None,scale_factor=1.,source_protection=verify(before),profile=profile,coarse_count=len(json.loads((out/'orientation_candidates.json').read_text())),full_planning_count=len(results),top4=[r['candidate_id'] for r in finalists],selected=best,physical_validation=False,runtime_seconds=time.monotonic()-start,runtime_scope=f'Current invocation stage={stage}; not cumulative project elapsed time',printer_upload=False,printer_start=False)
    save(out/'planning_results.json',results);save(out/'vascular_support_v2_summary.json',summary)
    for r in finalists:
        meta=dict(candidate_id=r['candidate_id'],selection_role='REVIEW_CANDIDATE' if not r['hard_gate_failures'] else 'DIAGNOSTIC_ONLY',scale_factor=1.,physical_validation=False)
        folder=Path(r['folder'])
        save(folder/'print_transform_v2.json',dict(**meta,matrix=r['orientation']['print_transform'],from_frame='FROZEN_CASTING',to_frame='X1C_PRINT'))
        save(folder/'casting_restore_transform_v2.json',dict(**meta,matrix=r['orientation']['casting_restore_transform'],from_frame='X1C_PRINT',to_frame='FROZEN_CASTING'))
    if best:
        shutil.copy2(Path(best['folder'])/'print_transform_v2.json',out/'print_transform_v2.json')
        shutil.copy2(Path(best['folder'])/'casting_restore_transform_v2.json',out/'casting_restore_transform_v2.json')
    save(out/'source_protection.json',verify(before))
    if stage in ['all','report','slice','coupon','removal']:
        from .support_v2_review import generate
        generate(out,summary,results)
    print('V2 RESULT',summary['status'],summary['diagnostic_candidate'],flush=True);return summary
