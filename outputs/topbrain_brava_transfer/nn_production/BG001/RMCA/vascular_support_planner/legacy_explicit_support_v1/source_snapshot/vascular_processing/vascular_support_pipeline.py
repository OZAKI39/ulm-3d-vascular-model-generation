"""Frozen-model support planning, bounded external search and reproducible review."""
from __future__ import annotations
import hashlib,json,shutil,time,traceback
from pathlib import Path
import numpy as np
import trimesh
import yaml
from .support_inputs import ROOT,discover,save,table,verify,digest,resolve_profile,transform,swc_branches,port_branches
from . import support_detection as detector
from . import support_contact_planner as contact
from . import support_slicer_adapter as slicer
from . import support_qc as qc

def cache_key(inputs,profile,cfg):
    content=dict(hashes=inputs['hashes'],orientation=inputs['print_transform'],layer_height=profile['layer_height'],config=cfg)
    return hashlib.sha256(json.dumps(content,sort_keys=True).encode()).hexdigest()

def sliced(engine,profiles,source,expected,out,timeout):
    key=hashlib.sha256((digest(source)+json.dumps(profiles,sort_keys=True)+str(engine.get('binary_sha256'))).encode()).hexdigest()
    result_path=out/'slice_result.json'
    if result_path.exists() and (out/'cache_key.txt').exists() and (out/'cache_key.txt').read_text()==key:
        r=json.loads(result_path.read_text())
    else:
        r=slicer.slice_local(engine,profiles,source,out,timeout);(out/'cache_key.txt').write_text(key)
    if r.get('slicer_success'):
        geometry=slicer.production_geometry(r['archive'],expected);r['production_geometry']=geometry
        if not geometry['passed']:r['slicer_success']=False;r['status']='SLICED_GEOMETRY_CHANGED'
        save(result_path,r)
    return r

def parameter_grid(cfg):
    # 3 tip diameters x 4 coupled stable-trunk options = 12, not a Cartesian
    # explosion. Each parameter remains explicit in candidate records.
    combinations=[]
    for tip in cfg['tip_diameters_mm']:
        for density,angle,diameter in [(20,15,2.),(30,25,2.5),(40,25,3.),(30,15,3.)]:
            combinations.append(dict(support_type='tree(manual)',tree_support_tip_diameter=str(tip),tree_support_top_rate=f'{density}%',
                tree_support_branch_diameter_organic=str(diameter),tree_support_angle_slow=str(angle)))
    return combinations[:cfg['maximum_slicer_candidates']]

def research(out,contacts,mesh,box,inner,outer,profile,cfg):
    from .support_tree_planner import build,illegal_collision
    from .support_scaffold import build as bridge_build,self_printability
    from .support_accessibility import component_access
    support,graph=build(contacts,mesh,box,inner,profile,cfg);support.export(out/'custom_clever_support.stl');save(out/'support_graph.json',graph)
    audit=illegal_collision(support,mesh,box);audit['access']=component_access(support.split(),mesh,box,outer);audit['self_printability']=self_printability(support,profile);audit['failed_routes']=graph['failed']
    save(out/'custom_support_qc.json',audit)
    rows=[];folder=out/'research';folder.mkdir(exist_ok=True)
    for d in cfg['pillar_diameters_mm']:
        for span in cfg['bridge_sensitivity_mm']:
            model,record=bridge_build(contacts,mesh,box,inner,profile,d,span);model.export(folder/f'bridge_pillar_d{d:g}_span{span:g}.stl')
            record['self_printability']=self_printability(model,profile,span);record['material_mm3']=float(model.volume);rows.append(record)
    save(out/'bridge_sensitivity.json',rows);table(out/'tables/bridge_sensitivity.csv',[dict(pillar_diameter_mm=r['diameter_mm'],bridge_limit_mm=r['max_bridge_span_mm'],bridge_count=len(r['bridges']),pillar_collision_count=sum(not p['passed'] for p in r['pillars']),unsupported_layer_regions=r['self_printability']['unsupported_layer_region_count'],material_mm3=r['material_mm3']) for r in rows])
    return dict(status='RESEARCH_CANDIDATE',slicer_success=False,source_unchanged=True,coverage=None,tip_count=None,planned_tip_count=len(graph['tips']),contact_area_mm2=None,
        planned_nominal_tip_disk_area_mm2=sum(np.pi*c['tip_diameter']**2/4 for c in contacts),
        support_near_bifurcation_count=sum(c['keepout_penalty'] for c in contacts),trapped_support_count=audit['access']['trapped_support_count'],top_accessible_fraction=audit['access']['top_accessible_fraction'],
        estimated_support_mass_g=float(support.volume)*.00104,estimated_print_time_seconds=None,
        hard_gate_failures=['CUSTOM_RESEARCH_NOT_PRODUCTION','SELF_SUPPORT_NOT_CONFIRMED']+(['CUSTOM_TIP_COLLISION'] if graph['failed'] else []))

def run(source=None,output=None,stage='all',config=None,mode='custom'):
    if mode=='benchmark':
        if stage not in ['all','plan','research','slice','report']:raise ValueError('INVALID_BENCHMARK_STAGE')
        return run_benchmark(source,output,stage,config)
    from .support_custom_pipeline import run as run_custom
    return run_custom(source,output,stage,config)

def run_benchmark(source=None,output=None,stage='all',config=None):
    began=time.monotonic();cfg=yaml.safe_load(Path(config or ROOT/'config/vascular_support_planner.yaml').read_text());inputs=discover(source=source)
    out=Path(output).resolve() if output else inputs['root'].parent/'vascular_support_planner/benchmark/recomputed';out.mkdir(parents=True,exist_ok=True)
    snapshot=out/'protected_before.json'
    before=json.loads(snapshot.read_text()) if snapshot.exists() else inputs['hashes'].copy();verify(before)
    before.update(inputs['hashes']);save(snapshot,before);save(out/'source_discovery.json',inputs);save(out/'resolved_config.json',cfg)
    profile,discovery=resolve_profile(inputs,out/'profiles');save(out/'profile_runtime.json',profile)
    mesh=trimesh.load_mesh(inputs['core']);mesh.apply_transform(inputs['print_transform']);box=trimesh.load_mesh(inputs['box']);box.apply_transform(inputs['print_transform'])
    inner=transform(inputs['inner_bounds'],inputs['print_transform']);outer=transform(inputs['outer_bounds'],inputs['print_transform'])
    try:
        if stage in ['all','plan']:
            core_start=time.monotonic();casting=np.array(json.loads(inputs['casting_transform'].read_text())['transform_4x4'])
            branches,topology=swc_branches(inputs['smoothed_swc'],np.array(inputs['print_transform'])@casting,inputs['branch_mapping'],cfg['sample_step_mm']);branches+=port_branches(inputs,cfg['sample_step_mm'])
            det,surface,layer,floating=detector.detect(mesh,branches,topology,profile,cfg,out,box)
            contacts,cov=contact.prune(det['raw'],cfg);contact.write_contacts(out,contacts,cov,cfg);save(out/'branches.json',branches);save(out/'topology.json',topology)
            enforcer,blocker=contact.masks(mesh,contacts);enforcer.export(out/'support_enforcer.stl');blocker.export(out/'support_blocker.stl')
            save(out/'support_requirement_field.json',dict(input_cache_key=cache_key(inputs,profile,cfg),triangle_field='support_requirement_surface.npz',candidate_field='tables/all_candidate_support_points.csv',detectors=['SURFACE_OVERHANG','LAYER_UNSUPPORTED','FREEFLOATING','SWC_LONG_SPAN','LOCAL_MINIMUM']))
            sensitivity=[]
            for span in cfg['span_sensitivity_mm']:
                altered={**cfg,'critical_max_span_mm':span};raw=[]
                for point in det['raw']:
                    p=point.copy();p['risk_level']='CRITICAL' if p['unsupported_span']>span or 'LOCAL_MINIMUM' in p['detectors_triggered'] else 'ATTENTION';raw.append(p)
                cc,cv=contact.prune(raw,altered);sensitivity.append(dict(maximum_span_mm=span,critical_region_count=sum(r['length_mm']>span for r in det['regions']),contact_count=len(cc),coverage_percent=cv['percent']))
            table(out/'tables/span_sensitivity.csv',sensitivity)
            save(out/'planning_runtime.json',dict(seconds=time.monotonic()-core_start,target_seconds=300,cache_key=cache_key(inputs,profile,cfg)))
            print('PLAN',det['critical_region_count'],'critical regions,',len(contacts),'contacts',flush=True)
        else:
            prior=json.loads((out/'support_requirement_field.json').read_text())
            if prior['input_cache_key']!=cache_key(inputs,profile,cfg):
                raise ValueError('PLANNING_CACHE_MISMATCH_RUN_STAGE_PLAN')
            det=json.loads((out/'detection.json').read_text());contacts=json.loads((out/'contact_plan.json').read_text())['contacts'];topology=json.loads((out/'topology.json').read_text())
        if stage in ['all','research']:
            try:r=research(out,contacts,mesh,box,inner,outer,profile,cfg)
            except Exception as exc:
                r=dict(status='CUSTOM_RESEARCH_FAILED',slicer_success=False,source_unchanged=True,contact_area_mm2=None,tip_count=None,coverage=None,
                    hard_gate_failures=['CUSTOM_RESEARCH_FAILED'],error=str(exc))
                save(out/'custom_failure.json',dict(error=str(exc),traceback=traceback.format_exc()));table(out/'tables/custom_failure.csv',[dict(error=str(exc))])
                from .support_aware_review import scene
                scene([(mesh,'#399d63',1),(box,'#99aabb',.15)],out/'custom_failure.png','研究支撑失败，主生产引导路线继续')
            else:
                if (out/'custom_failure.json').exists():
                    archived=out/'history/recovered_custom_failure';archived.mkdir(parents=True,exist_ok=True)
                    for name in ['custom_failure.json','custom_failure.png']:
                        if (out/name).exists():shutil.move(out/name,archived/name)
            save(out/'custom_comparison.json',r)
        records={};engines={}
        if stage in ['all','slice']:
            raw=det['raw']
            for name in ['Bambu','Orca']:
                engines[name]=slicer.discover_engine(name,discovery,out)
                folder=out/'slicer'/(name.lower()+'_native_current');r=sliced(engines[name],slicer.engine_profiles(discovery,name),inputs['oriented_mold'],inputs['oriented_mold'],folder,cfg['slice_timeout_seconds'])
                q=qc.audit_slice(r,mesh,box,outer,raw,[],profile,cfg,folder);q['folder']=str(folder);records[name.upper()+'_NATIVE_BASELINE']=q
                print(name,r['status'],'coverage',(q.get('coverage') or {}).get('percent'),flush=True)
            grid=parameter_grid(cfg);search=[];model=trimesh.load_mesh(inputs['oriented_mold'])
            for i,params in enumerate(grid):
                folder=out/'slicer'/f'guided_{i:02}';folder.mkdir(exist_ok=True,parents=True);tip=float(params['tree_support_tip_diameter']);cc=[{**c,'tip_diameter':tip} for c in contacts]
                qc.nozzle_gate(tip,profile)
                enforcer,blocker=contact.masks(mesh,cc,tip);profiles=slicer.engine_profiles(discovery,'Orca',params);project=folder/'guided_project.3mf'
                slicer.write_project(project,model,enforcer,blocker,profiles);r=sliced(engines['Orca'],profiles,project,inputs['oriented_mold'],folder,cfg['slice_timeout_seconds'])
                q=qc.audit_slice(r,mesh,box,outer,raw,cc,profile,cfg,folder);q.update(folder=str(folder),parameters=params,index=i,project=str(project));search.append(q)
                save(out/'orca_search.json',search);print('GUIDED',i,r['status'],(q.get('coverage') or {}).get('percent'),q.get('contact_area_mm2'),q.get('hard_gate_failures'),flush=True)
            feasible=[q for q in search if not q.get('hard_gate_failures',['FAILED'])]
            ranked=sorted(feasible or search,key=lambda q:(len(q.get('hard_gate_failures',['FAILED'])),q.get('contact_area_mm2') or 1e9,q.get('support_near_bifurcation_count',999),q.get('tip_count') or 999,-(q.get('top_accessible_fraction') or 0),q.get('estimated_support_mass_g') or 1e9,q.get('estimated_print_time_seconds') or 1e9,q['index']))
            best=ranked[0];records['GUIDED_ORCA']=best
            shutil.copy2(best['project'],out/'BG001_RMCA_60x50_guided_support.3mf')
            cc=[{**c,'tip_diameter':float(best['parameters']['tree_support_tip_diameter'])} for c in contacts];enforcer,blocker=contact.masks(mesh,cc,cc[0]['tip_diameter']);enforcer.export(out/'support_enforcer.stl');blocker.export(out/'support_blocker.stl')
            if (Path(best['folder'])/'sliced.3mf').exists():shutil.copy2(Path(best['folder'])/'sliced.3mf',out/'BG001_RMCA_60x50_guided_support_sliced.3mf')
            save(out/'baseline_metrics.json',records)
        elif (out/'baseline_metrics.json').exists():
            records=json.loads((out/'baseline_metrics.json').read_text())
            engines={name:json.loads((out/(name.lower()+'_engine.json')).read_text()) for name in ['Bambu','Orca'] if (out/(name.lower()+'_engine.json')).exists()}
        if stage=='plan' or stage=='research':
            save(out/'source_protection.json',verify(before));return dict(output=str(out),stage=stage)
        if (out/'custom_comparison.json').exists():records['CUSTOM_CLEVER_TREE']=json.loads((out/'custom_comparison.json').read_text())
        comparisons,selected=qc.compare(records);table(out/'tables/baseline_comparison.csv',comparisons)
        status='READY_FOR_HUMAN_REVIEW' if selected and selected=='GUIDED_ORCA' else 'NEEDS_ADJUSTMENT'
        old_summary=json.loads((out/'vascular_support_summary.json').read_text()) if (out/'vascular_support_summary.json').exists() else {}
        summary=dict(status=status,input_hashes=before,source_protection=verify(before),open_source_sources=json.loads((ROOT/'third_party/support_sources.json').read_text()),profile=profile,
            detection=dict(critical_regions=det['critical_region_count'],local_minimum_anchor_points=sum('LOCAL_MINIMUM' in p['detectors_triggered'] for p in det['raw']),longest_span_mm=det['longest_dangerous_span_mm'],raw_points=len(det['raw']),selected_points=len(contacts)),
            coverage=records.get('GUIDED_ORCA',{}).get('coverage'),contact=records.get('GUIDED_ORCA',{}),baseline_comparison=comparisons,selected_candidate=selected,
            orca=engines.get('Orca'),bambu=engines.get('Bambu'),total_seconds=time.monotonic()-began,printer_upload=False,printer_start=False,physical_print_validated=False)
        g=records.get('GUIDED_ORCA',{})
        summary['support']=dict(topology='EXTERNAL_ORCA_ORGANIC_TREE',total_length_mm=g.get('support_path_length_mm'),material_estimate_g=g.get('estimated_support_mass_g'),
            trapped_components=g.get('trapped_support_count'),top_accessible_fraction=g.get('top_accessible_fraction'),chosen_parameters=g.get('parameters'),
            artifact_directory=g.get('folder'),representation='Actual extruded toolpaths; no native separable support solid exported')
        summary['core_planning_runtime']=json.loads((out/'planning_runtime.json').read_text())
        summary['execution_stage']=stage
        summary['architecture_role']='BENCHMARK_ONLY'
        summary['full_run_seconds']=summary['total_seconds'] if stage=='all' else old_summary.get('full_run_seconds',old_summary.get('total_seconds'))
        base=records.get('BAMBU_NATIVE_BASELINE',{})
        reduction=100*(1-g['contact_area_mm2']/base['contact_area_mm2']) if g.get('contact_area_mm2') is not None and base.get('contact_area_mm2') else None
        summary['target_improvements']=dict(contact_area_reduction_percent=reduction,
            at_least_25_percent_area_reduction=reduction is not None and reduction>=25,
            at_least_40_percent_area_reduction=reduction is not None and reduction>=40,
            fewer_contact_intervals=g.get('tip_count') is not None and base.get('tip_count') is not None and g['tip_count']<base['tip_count'],
            zero_bifurcation_neighborhood_contacts=g.get('support_near_bifurcation_count')==0,
            contact_count_definition='Estimated disconnected arc contact intervals, not exact slicer tree-tip count')
        save(out/'vascular_support_summary.json',summary)
        from .support_review import figures,report
        figures(out,mesh,box,det,contacts,topology,records,selected);report(out,summary)
        if (out/'failure.json').exists():
            archived=out/'history/recovered_failure';archived.mkdir(parents=True,exist_ok=True)
            for name in ['failure.json','failure.png']:
                if (out/name).exists():shutil.move(out/name,archived/name)
            if (out/'tables/failure.csv').exists():shutil.move(out/'tables/failure.csv',archived/'failure.csv')
        save(out/'source_protection.json',verify(before));print('RESULT',status,selected,out,flush=True);return summary
    except Exception as exc:
        save(out/'failure.json',dict(error=str(exc),traceback=traceback.format_exc()))
        save(out/'vascular_support_summary.json',dict(status='NEEDS_ADJUSTMENT',error=str(exc),input_hashes=before,physical_print_validated=False))
        table(out/'tables/failure.csv',[dict(error=str(exc))]);(out/'vascular_support_report.md').write_text('# 支撑规划未通过\n\n'+str(exc)+'\n\n已保留中间数据、日志与失败原因。\n')
        from .support_aware_review import scene
        scene([(mesh,'#399d63',1),(box,'#99aabb',.15)],out/'failure.png','支撑规划未通过：保留真实输入与诊断')
        raise
    finally:
        save(out/'source_protection.json',verify(before))
