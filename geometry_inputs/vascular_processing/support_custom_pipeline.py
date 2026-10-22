"""s1-8 explicit production route; prior slicer-generated trees remain benchmarks."""
from __future__ import annotations
import copy,csv,hashlib,json,shutil,time,traceback,zipfile
from pathlib import Path
import numpy as np
import trimesh
import yaml
from .support_inputs import ROOT,discover,resolve_profile,swc_branches,port_branches,transform,save,table,verify,digest
from .support_custom_detection import fused_detect,span_sensitive_raw
from .support_contact_planner import prune,coverage,write_contacts
from . import support_custom_graph as graphing
from . import support_custom_qc as qc
from . import support_slicer_adapter as slicer
from .support_breakaway import coupon
from .support_removal import evaluate

def cache_identity(inputs,profile,cfg):
    return hashlib.sha256(json.dumps(dict(source=inputs['hashes'],orientation=inputs['print_transform'],profile=profile,config=cfg),sort_keys=True,default=str).encode()).hexdigest()

def bottom_mesh(box,inner):
    lo,hi=box.bounds;size=hi-lo+2;size[2]=inner[0,2]-lo[2]+1
    slab=trimesh.creation.box(size);slab.apply_translation([(lo[0]+hi[0])/2,(lo[1]+hi[1])/2,inner[0,2]-size[2]/2])
    return trimesh.boolean.intersection([box,slab],engine='manifold')

def write_components(folder,components):
    expected={f'broken_component_{i:03}.stl' for i in range(len(components))}
    for previous in folder.glob('broken_component_*.stl'):
        if previous.name not in expected:
            history=folder/'history/obsolete_components';history.mkdir(parents=True,exist_ok=True)
            previous.rename(history/(previous.stem+'_'+digest(previous)[:12]+previous.suffix))
    for i,part in enumerate(components):qc.canonical_stl(part).export(folder/f'broken_component_{i:03}.stl')

def gap_audit(support,graph,vascular,profile):
    tips=[t for t in graph['tips'] if t['root_connected']]
    rows=[dict(tip_id=t['tip_id'],target_mm=profile['support_top_gap'],measured_triangle_gap_mm=t['measured_top_gap_mm'],tolerance_mm=.08,passed=abs(t['measured_top_gap_mm']-profile['support_top_gap'])<=.08) for t in tips]
    body=[]
    for e in graph['edges']:
        if e['kind']=='taper':continue
        part=graphing.edge_geometry(graph,e)
        distance=graphing.minimum_distance(part,vascular);body.append(dict(edge_id=e['edge_id'],measured_triangle_distance_mm=distance,target_mm=profile['support_xy_gap'],tolerance_mm=.015,passed=distance>=profile['support_xy_gap']-.015))
    return dict(passed=bool(rows) and all(r['passed'] for r in rows+body),tips=rows,noncontact_edges=body,method='FCL triangle mesh minimum distance; finite tip cap included')

def geometry_candidate(out,index,contacts,raw,mesh,box,mold,inner,profile,cfg,tip,neck,pillar):
    folder=out/'custom_candidates'/f'custom_{index:02}';folder.mkdir(parents=True,exist_ok=True);began=time.monotonic()
    support,graph,feet=graphing.build(contacts,mesh,box,inner,profile,cfg,tip,neck,pillar)
    if not len(support.faces):raise ValueError('NO_ROOTED_SUPPORT_GEOMETRY')
    support=qc.canonical_stl(support)
    support.export(folder/'custom_support.stl');save(folder/'support_graph.json',graph)
    rooted={c for t in graph['tips'] if t['root_connected'] for c in t['contact_ids']};actual_contacts=[c for c in contacts if c['contact_id'] in rooted]
    cov=coverage(raw,actual_contacts,cfg);union=qc.root_union_audit(support,mesh,box,feet,graph['roots'],inner)
    bottom=bottom_mesh(box,inner);selftest=qc.self_printability(support,bottom,profile,folder)
    # Repair an unstable merged topology by falling back to rooted pillars.
    # The vascular contact set stays fixed throughout this repair.
    if not selftest['passed']:
        repaired,repaired_graph,repaired_feet=graphing.build(contacts,mesh,box,inner,profile,{**cfg,'merge_max_distance_mm':0},tip,neck,pillar)
        repaired=qc.canonical_stl(repaired)
        repair_check=qc.self_printability(repaired,bottom,profile)
        if repair_check['unsupported_regions']<selftest['unsupported_regions']:
            support,graph,feet=repaired,repaired_graph,repaired_feet;graph['stability_repairs'].append(dict(action='DISABLE_UNPRINTABLE_MERGES_USE_VERTICAL_PILLARS',added_vascular_tips=0))
            support.export(folder/'custom_support.stl');save(folder/'support_graph.json',graph);selftest=qc.self_printability(support,bottom,profile,folder)
            rooted={c for t in graph['tips'] if t['root_connected'] for c in t['contact_ids']};actual_contacts=[c for c in contacts if c['contact_id'] in rooted]
            cov=coverage(raw,actual_contacts,cfg);union=qc.root_union_audit(support,mesh,box,feet,graph['roots'],inner)
    components,removal=evaluate(support,graph['roots'],graph['tips'],mesh,box)
    write_components(folder,components)
    save(folder/'removal.json',removal);table(folder/'support_component_removal.csv',removal['rows']);gaps=gap_audit(support,graph,mesh,profile);save(folder/'gap_audit.json',gaps)
    assembly=qc.canonical_stl(trimesh.boolean.union([mold,support],engine='manifold'));assembly.export(folder/'final_mold_with_breakaway_support.stl')
    lost=trimesh.boolean.difference([mold,assembly],engine='manifold');union['frozen_mold_removed_volume_mm3']=abs(float(lost.volume)) if len(lost.faces) else 0.
    union['passed'] &= union['frozen_mold_removed_volume_mm3']<1e-5
    foundation=qc.canonical_stl(trimesh.boolean.union([bottom,support],engine='manifold'));foundation.export(folder/'support_only_with_approved_bottom.stl')
    result=dict(candidate_id=f'custom_{index:02}',folder=str(folder),source_unchanged=True,tip_width_mm=tip,neck_width_mm=neck,pillar_diameter_mm=pillar,
        planned_tip_count=sum(t['root_connected'] for t in graph['tips']),planned_coverage=cov,self_printability=selftest,union_audit=union,gap_audit=gaps,
        removal=removal,edges=graph['edges'],bifurcation_contacts=sum(t['root_connected'] and t['distance_to_bifurcation']<cfg['bifurcation_keepout_mm'] for t in graph['tips']),
        support_volume_mm3=float(support.volume),support_mass_estimate_g=float(support.volume)*.00104,graph_runtime_seconds=time.monotonic()-began,
        rejected_routes=graph['rejected_routes'],total_path_length_mm=sum(e['length_mm'] for e in graph['edges']),slice={})
    save(folder/'candidate.json',result)
    print('CUSTOM GEOMETRY',index,'tips',result['planned_tip_count'],'coverage',cov['percent'],'self',selftest['passed'],'trapped',removal['trapped_support_count'],flush=True)
    return result

def profiles_support_off(discovery):
    profiles=slicer.engine_profiles(discovery,'Bambu');profiles['process']['enable_support']='0'
    return profiles

def slice_candidate(r,inputs,discovery,profile,cfg,mesh,raw,contacts):
    from .vascular_support_pipeline import sliced
    folder=Path(r['folder']);engine=slicer.discover_engine('Bambu',discovery,folder);profiles=profiles_support_off(discovery);graph=json.loads((folder/'support_graph.json').read_text())
    for name,filename in [('assembly','final_mold_with_breakaway_support.stl'),('support_only','support_only_with_approved_bottom.stl')]:
        sliced_record=sliced(engine,profiles,folder/filename,folder/filename,folder/('bambu_'+name),cfg['slice_timeout_seconds'])
        support=trimesh.load_mesh(folder/'custom_support.stl');metrics=qc.audit_slice(sliced_record,support,mesh,raw,contacts,graph,profile,cfg,folder/('bambu_'+name))
        r['slice' if name=='assembly' else 'support_only_slice']=metrics
        if name=='assembly':r['actual_coverage']=metrics.get('coverage');r['estimated_print_time_seconds']=metrics.get('estimated_print_time_seconds')
    return r

def write_vtp(folder,graph):
    import pyvista as pv
    tips=[p for p in graph['tips'] if p['root_connected']];p=pv.PolyData(np.array([t['xyz'] for t in tips]));p['tip_width_mm']=[t['width_mm'] for t in tips];p.save(folder/'support_tips.vtp')
    lines=[];width=[]
    for e in graph['edges']:lines.extend([graph['nodes'][e['start']]['xyz'],graph['nodes'][e['end']]['xyz']]);width.append(e['diameter_mm'])
    v=pv.PolyData(np.array(lines));v.lines=np.c_[np.full(len(width),2),np.arange(len(lines)).reshape(-1,2)].ravel();v.cell_data['diameter_mm']=np.r_[np.zeros(len(lines)),width] if v.n_cells>len(width) else width;v.save(folder/'support_graph.vtp')

def run(source=None,output=None,stage='all',config=None):
    raise RuntimeError('LEGACY_EXPLICIT_SUPPORT_V1_RETIRED: use s1-8 default V2; frozen V1 outputs and source snapshot remain benchmark evidence')
    began=time.monotonic();inputs=discover(source=source);out=Path(output).resolve() if output else inputs['root'].parent/'vascular_support_planner';out.mkdir(parents=True,exist_ok=True)
    before=json.loads((out/'protected_before.json').read_text()) if (out/'protected_before.json').exists() else inputs['hashes'];verify(before);before.update(inputs['hashes']);save(out/'protected_before.json',before)
    cfg=yaml.safe_load(Path(config or ROOT/'config/vascular_support_planner.yaml').read_text());cfg['selection_priority']='vascular';cfg['merge_max_distance_mm']=cfg['custom_support']['merge_max_distance_mm']
    profile,discovery=resolve_profile(inputs,out/'profiles');save(out/'profile_runtime.json',profile);save(out/'source_discovery.json',inputs);save(out/'resolved_config.json',cfg)
    key=cache_identity(inputs,profile,cfg);mesh=trimesh.load_mesh(inputs['core']);mesh.apply_transform(inputs['print_transform']);box=trimesh.load_mesh(inputs['box']);box.apply_transform(inputs['print_transform']);mold=trimesh.load_mesh(inputs['oriented_mold']);inner=transform(inputs['inner_bounds'],inputs['print_transform'])
    try:
        if stage in ['all','plan']:
            start=time.monotonic();casting=np.array(json.loads(inputs['casting_transform'].read_text())['transform_4x4']);branches,topology=swc_branches(inputs['smoothed_swc'],np.array(inputs['print_transform'])@casting,inputs['branch_mapping'],cfg['sample_step_mm']);branches+=port_branches(inputs,cfg['sample_step_mm'])
            det=fused_detect(mesh,box,branches,topology,inner,profile,cfg,out);contacts,cov=prune(det['raw'],cfg)
            for p in contacts:
                if p['keepout_penalty']:p['status']='SUPPORT_NEAR_BIFURCATION_EXCEPTION'
            write_contacts(out,contacts,cov,cfg);save(out/'branches.json',branches);save(out/'topology.json',topology);table(out/'tables/selected_support_tips.csv',contacts)
            save(out/'custom_detection_cache.json',dict(key=key,seconds=time.monotonic()-start));sensitivity=[]
            for span in cfg['span_sensitivity_mm']:
                altered={**cfg,'critical_max_span_mm':span};raw=span_sensitive_raw(det['raw'],span)
                c,v=prune(raw,altered);sensitivity.append(dict(critical_max_span_mm=span,selected_tips=len(c),planned_coverage=v['percent']))
            table(out/'tables/span_sensitivity.csv',sensitivity);print('CUSTOM PLAN',len(det['critical_intervals']),'intervals,',len(contacts),'contacts',flush=True)
        else:
            if json.loads((out/'custom_detection_cache.json').read_text())['key']!=key:raise ValueError('CUSTOM_PLANNING_CACHE_MISMATCH_RUN_PLAN')
            det=json.loads((out/'detection.json').read_text());contacts=json.loads((out/'contact_plan.json').read_text())['contacts']
        if stage=='plan':return dict(status='NEEDS_ADJUSTMENT',stage=stage,output=str(out))
        if stage in ['all','geometry','research']:
            candidates=[];w=max(profile['line_width'],profile['nozzle_diameter'])
            for i,(t,n,d) in enumerate(zip(cfg['custom_support']['tip_width_lines'],cfg['custom_support']['neck_width_lines'],cfg['custom_support']['pillar_width_lines'])):
                candidates.append(geometry_candidate(out,i,contacts,det['raw'],mesh,box,mold,inner,profile,cfg,t*w,n*w,d*w));save(out/'custom_search.json',candidates)
        else:candidates=json.loads((out/'custom_search.json').read_text())
        if stage=='removal':
            for r in candidates:
                folder=Path(r['folder']);g=json.loads((folder/'support_graph.json').read_text());s=trimesh.load_mesh(folder/'custom_support.stl');parts,removal=evaluate(s,g['roots'],g['tips'],mesh,box)
                r['removal']=removal
                write_components(folder,parts)
                save(folder/'removal.json',removal);table(folder/'support_component_removal.csv',removal['rows']);save(folder/'candidate.json',r);save(out/'custom_search.json',candidates)
                print('CUSTOM REMOVAL',r['candidate_id'],removal['component_count'],removal['trapped_support_count'],flush=True)
        if stage in ['all','coupon']:
            cm,rows,features=coupon(out/'coupon',profile)
            from .vascular_support_pipeline import sliced
            engine=slicer.discover_engine('Bambu',discovery,out/'coupon');profiles=profiles_support_off(discovery)
            result=sliced(engine,profiles,out/'coupon/support_breakaway_coupon.stl',out/'coupon/support_breakaway_coupon.stl',out/'coupon/bambu',cfg['slice_timeout_seconds'])
            if result.get('slicer_success'):
                paths=qc.parse_actual_paths(result['archive']);points=qc.sample_segments(paths['segments'],.10);preserve=qc.feature_preservation(points,features,profile);save(out/'coupon/feature_preservation.json',preserve)
                with zipfile.ZipFile(result['archive']) as z:settings=json.loads(z.read('Metadata/project_settings.config'))
                result['auto_support_off']=str(settings.get('enable_support')) in ['0','false','False']
                result['native_support_toolpath_count']=int((paths['kind']==1).sum());result['all_coupon_features_preserved']=preserve['passed']
            save(out/'coupon/slice_validation.json',result)
        if stage in ['all','slice']:
            for r in candidates:
                slice_candidate(r,inputs,discovery,profile,cfg,mesh,det['raw'],contacts);save(Path(r['folder'])/'candidate.json',r);save(out/'custom_search.json',candidates)
                print('CUSTOM SLICE',r['candidate_id'],r['slice'].get('slicer_success'),'OFF',r['slice'].get('auto_support_off'),flush=True)
        benchmark=out/'benchmark/guided_v1/baseline_metrics.json';baseline=json.loads(benchmark.read_text()) if benchmark.exists() else {}
        base_area=baseline.get('BAMBU_NATIVE_BASELINE',{}).get('contact_area_mm2')
        for r in candidates:
            area=r['slice'].get('contact_area_mm2');r['contact_area_reduction_percent']=100*(1-area/base_area) if area is not None and base_area else -1
            r['hard_gate_failures']=qc.hard_gates(r,cfg);save(Path(r['folder'])/'candidate.json',r)
        # Failed variants remain diagnostics. A successful native/guided tree
        # can NEVER be promoted into production through this selection.
        feasible=[r for r in candidates if not r['hard_gate_failures']]
        rank=lambda r:(r['slice'].get('contact_area_mm2') or 1e9,r['bifurcation_contacts'],r['planned_tip_count'],r['removal']['trapped_support_count'],r['support_volume_mm3'],r.get('estimated_print_time_seconds') or 1e9,r['candidate_id'])
        best=min(feasible,key=rank) if feasible else min(candidates,key=lambda r:(-(r.get('actual_coverage') or r['planned_coverage'])['percent'],len(r['hard_gate_failures']),rank(r)))
        geometry=out/'geometry';geometry.mkdir(exist_ok=True);folder=Path(best['folder']);graph=json.loads((folder/'support_graph.json').read_text())
        for name in ['custom_support.stl','final_mold_with_breakaway_support.stl']:shutil.copy2(folder/name,geometry/name)
        shutil.copy2(folder/'support_graph.json',out/'support_graph.json');shutil.copy2(inputs['oriented_mold'],geometry/'frozen_mold.stl');mesh.export(geometry/'frozen_vascular_with_ports.stl');write_vtp(geometry,graph)
        for path in ['bambu_assembly/sliced.3mf','bambu_support_only/sliced.3mf']:
            keyname='slice' if path.startswith('bambu_assembly') else 'support_only_slice'
            destination=out/('custom_support_auto_support_off.3mf' if keyname=='slice' else 'custom_support_only_auto_support_off.3mf')
            if best.get(keyname,{}).get('slicer_success') and (folder/path).exists():shutil.copy2(folder/path,destination)
            elif destination.exists():
                archived=out/'history/stale_geometry_slices';archived.mkdir(parents=True,exist_ok=True)
                destination.rename(archived/(destination.stem+'_'+digest(destination)[:12]+destination.suffix))
        table(out/'tables/support_edges.csv',graph['edges']);table(out/'tables/support_roots.csv',graph['roots']);table(out/'tables/support_component_removal.csv',best['removal']['rows']);table(out/'tables/rejected_support_routes.csv',graph['rejected_routes'])
        final_tips=[t for t in graph['tips'] if t['root_connected']];selected_rows=[]
        preserved={r['feature_id'] for r in (best['slice'].get('features') or {}).get('rows',[]) if r.get('preserved')}
        for c in contacts:
            matched=[t for t in final_tips if c['contact_id'] in t['contact_ids']]
            selected_rows.append({**c,'tip_diameter':best['tip_width_mm'],'generated_tip_ids':[t['tip_id'] for t in matched],
                'actual_tip_xyz':[t['xyz'] for t in matched],'measured_gap_mm':[t['measured_top_gap_mm'] for t in matched],
                'geometry_generated':bool(matched),'slicer_preserved':bool(matched) and all(t['tip_id'] in preserved for t in matched)})
        table(out/'tables/selected_support_tips.csv',selected_rows);table(out/'tables/actual_support_tips.csv',final_tips)
        missing=set((best.get('actual_coverage') or best['planned_coverage'])['uncovered_point_ids']);interval_rows=[]
        for interval in det['critical_intervals']:
            n=len(interval['point_ids']);failed=len(set(interval['point_ids'])&missing)
            interval_rows.append({**interval,'coverage_percent':100*(n-failed)/n,'uncovered_sample_count':failed,'coverage_status':'COVERED' if not failed else 'UNCOVERED_CRITICAL'})
        table(out/'tables/critical_support_intervals.csv',interval_rows)
        comparison=[]
        for name in ['BAMBU_NATIVE_BASELINE','ORCA_NATIVE_BASELINE','GUIDED_ORCA']:
            r=baseline.get(name,{})
            evidence=Path(r.get('folder',''))/'actual_contact_coverage.csv';recomputed=None
            if evidence.is_file():
                ids={x['point_id'] for x in csv.DictReader(evidence.open()) if x['actual_support_underneath']=='True'}
                recomputed=coverage(det['raw'],[p for p in det['raw'] if p['point_id'] in ids],cfg)
            comparison.append(dict(method=name,role='BENCHMARK_ONLY',tip_count=r.get('tip_count'),tip_count_definition='estimated arc contact intervals',contact_area_mm2=r.get('contact_area_mm2'),coverage_percent=(recomputed or {}).get('percent'),trapped_support_count=r.get('trapped_support_count'),removal_definition='legacy sampled probe, not whole-component removal',support_volume_mm3=r.get('estimated_support_volume_mm3'),source_geometry='same frozen mold; coverage reevaluated against custom fused requirements'))
        comparison.append(dict(method='CUSTOM_LITERATURE_AWARE_SUPPORT',role='PRODUCTION_CANDIDATE',tip_count=best['slice'].get('actual_preserved_tip_count',best['planned_tip_count']),tip_count_definition='explicit geometry tips confirmed in toolpaths when sliced',contact_area_mm2=best['slice'].get('contact_area_mm2'),coverage_percent=(best.get('actual_coverage') or best['planned_coverage'])['percent'],trapped_support_count=best['removal']['trapped_support_count'],support_volume_mm3=best['support_volume_mm3']))
        table(out/'tables/baseline_comparison.csv',comparison);save(out/'custom_search.json',candidates)
        status='READY_FOR_HUMAN_REVIEW_AND_COUPON_TEST' if feasible else 'NEEDS_ADJUSTMENT'
        summary=dict(status=status,architecture='CUSTOM_LITERATURE_AWARE_SUPPORT',scale_factor=1.,production_selected=best['candidate_id'] if feasible else None,diagnostic_candidate=best['candidate_id'],
            source_protection=verify(before),input_hashes=before,profile=profile,config=cfg,detection=dict(critical_intervals=len(det['critical_intervals']),longest_span_mm=det['longest_dangerous_span_mm'],freefloating_raw_points=det['freefloating_raw_count'],raw_candidates=len(det['raw']),selected_contacts=len(contacts)),
            selected=best,comparison=comparison,benchmarks_role='BENCHMARK_ONLY',runtime_seconds=time.monotonic()-began,detector_runtime=json.loads((out/'custom_detection_cache.json').read_text()),
            coupon=str(out/'coupon'),physical_print_validated=False,printer_upload=False,printer_start=False)
        summary['core_geometry_search_seconds_including_geometry_and_removal_qc']=summary['detector_runtime']['seconds']+sum(r['graph_runtime_seconds'] for r in candidates)
        summary['critical_coverage_by_domain']={}
        for label,is_port in [('anatomical_vascular',False),('artificial_ports',True)]:
            points=[p for p in det['raw'] if p['risk_level']=='CRITICAL' and p['is_port']==is_port];uncovered=sum(p['point_id'] in missing for p in points)
            summary['critical_coverage_by_domain'][label]=dict(required_samples=len(points),uncovered_samples=uncovered,coverage_percent=100*(len(points)-uncovered)/len(points) if points else 100.)
        save(out/'vascular_support_summary.json',summary)
        from .support_custom_review import generate
        generate(out,mesh,box,det,contacts,graph,best,summary)
        for name in ['custom_pipeline_failure.json','custom_pipeline_failure.png','tables/custom_pipeline_failure.csv']:
            old=out/name
            if old.exists():
                recovered=out/'history/recovered_execution_errors';recovered.mkdir(parents=True,exist_ok=True)
                old.rename(recovered/(old.stem+'_'+digest(old)[:12]+old.suffix))
        from .support_artifacts import write_manifests
        write_manifests(out,summary)
        print('CUSTOM RESULT',status,best['candidate_id'],best['hard_gate_failures'],flush=True);return summary
    except Exception as exc:
        save(out/'custom_pipeline_failure.json',dict(error=str(exc),traceback=traceback.format_exc()));save(out/'vascular_support_summary.json',dict(status='NEEDS_ADJUSTMENT',architecture='CUSTOM_LITERATURE_AWARE_SUPPORT',error=str(exc)))
        table(out/'tables/custom_pipeline_failure.csv',[dict(error=str(exc))]);(out/'vascular_support_report.md').write_text('# 显式支撑生成未完成\n\n'+str(exc)+'\n')
        from .support_aware_review import scene
        scene([(mesh,'#399d63',1),(box,'#99aabb',.12)],out/'custom_pipeline_failure.png','显式支撑未通过：保留输入与失败证据')
        raise
    finally:save(out/'source_protection.json',verify(before))
