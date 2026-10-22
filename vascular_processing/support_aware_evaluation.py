"""Fixed-profile local Bambu comparison and top-only removal diagnostics."""
from __future__ import annotations
import copy
import json
import logging
from pathlib import Path
import re
import shutil
import warnings
import xml.etree.ElementTree as ET
import zipfile
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from .support_aware_roi import ROOT,RMCA,DEFAULT_OUTPUT,verify_snapshot
from .support_aware_port_design import configuration,design,load_companion,load_design,ensure_current_design
from .sacrificial_fixture_review import save_json,polydata
from .sacrificial_fixture import export_checked_stl,write_csv
from .topbrain_qc import sha256
from . import support_removal_qc as access
from . import bambu_casting_mold as bambu
from . import bambu_slice_adapter as adapter
from . import bambu_manufacturing as native_bambu
from . import casting_mold_qc as qc
from . import abs_casting_mold as mold
from . import roi_support_metrics as metrics


def archived_warnings(path):
    """Read warning ELEMENTS as well as attributes; keep noncritical warnings."""
    records=[]
    with zipfile.ZipFile(path) as archive:
        for name in archive.namelist():
            if not name.endswith(('.xml','.config')):continue
            try:root=ET.fromstring(archive.read(name))
            except ET.ParseError:continue
            for e in root.iter():
                if 'warning' in e.tag.lower():records.append(dict(source=name,**e.attrib,text=e.text or ''))
                for key,value in e.attrib.items():
                    if 'warning' in key.lower():records.append(dict(source=name,message=value))
    critical=[r for r in records if re.search(r'unsupported|floating|thin[ _-]?feature',json.dumps(r),re.I)]
    return records,critical


def fixed_discovery(out,cfg):
    folder=out/'bambu/evaluation_profile';folder.mkdir(parents=True,exist_ok=True)
    saved=folder/'frozen_discovery.json'
    if saved.exists():
        discovery=json.loads(saved.read_text())
        hashes=json.loads((folder/'sha256.json').read_text())
        if any(sha256(folder/name)!=h for name,h in hashes.items()):raise ValueError('EVALUATION_PROFILE_CHANGED')
        save_json(out/'resolved_evaluation_config.json',cfg)
        return discovery
    discovery=bambu.discover(cfg,folder)
    if discovery['status']=='BAMBU_PROFILE_CONFIRMED':
        config=dict(slicer=dict(enable_support=True))
        paths=native_bambu.prepare_slicer_profiles(discovery,config,folder)
        save_json(folder/'sha256.json',{Path(p).name:sha256(p) for p in paths.values()})
    else:save_json(folder/'sha256.json',{})
    save_json(saved,discovery);save_json(out/'resolved_evaluation_config.json',cfg)
    return discovery


def support_contact(paths,mesh,partition,transform,cfg):
    """Area-weighted triangle-centroid contact envelope around ACTUAL support.

    A G-code file has no separate support solid. This finite-resolution contact
    estimate must not be reported as an exact physical contact area. The mass
    reduction criterion uses measured positive E extrusion instead.
    """
    support=paths['segments'][paths['kind']==1]
    points=bambu.sample_segments(support,.2)
    ids=partition['face_ids'];labels=np.asarray(partition['labels'])
    central=~np.array([str(s).startswith(('I','O','BOX')) for s in labels])
    if not len(points):return dict(vascular_contact_area_mm2=0.,core_and_port_contact_area_mm2=0.,method='ACTUAL_GCODE_SUPPORT_PROXIMITY_TRIANGLE_AREA_ESTIMATE',sampling_step_mm=.2,contact_envelope_mm=.5)
    t=np.asarray(transform);centers=mesh.triangles_center[ids]@t[:3,:3].T+t[:3,3]
    near=cKDTree(points).query(centers)[0]<=.5
    return dict(vascular_contact_area_mm2=float(mesh.area_faces[ids[near&central]].sum()),
        core_and_port_contact_area_mm2=float(mesh.area_faces[ids[near]].sum()),
        method='ACTUAL_GCODE_SUPPORT_PROXIMITY_TRIANGLE_AREA_ESTIMATE',sampling_step_mm=.2,contact_envelope_mm=.5,
        exact_physical_interface_reconstructed=False,triangle_area_weighted=True)


def profile_hashes(record):
    return {p.name:sha256(p) for p in sorted((Path(record['archive']).parent/'profiles').glob('*.json'))}


def slice_cached(discovery,cfg,row,out,logger,mesh,partition,box,vascular_only=False):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);saved=out/'evaluation_slice.json'
    if saved.exists():return json.loads(saved.read_text())
    if discovery['status']!='BAMBU_PROFILE_CONFIRMED':
        result=dict(status='READY_FOR_HUMAN_REVIEW_WITH_BAMBU_PENDING',results=[],printer_job_sent=False,reason=discovery['status'])
        save_json(saved,result);return result
    result=adapter.slice_top_five(discovery,cfg,[row],out,logger)
    for record in result['results']:
        if not record.get('archive'):continue
        if vascular_only:
            paths=bambu.parse_paths(record['archive'])
            folder=Path(record['archive']).parent
            np.savez_compressed(folder/'actual_extrusion_paths.npz',**{k:paths[k] for k in ('segments','kind','extrusion_filament_mm')})
            record['actual_paths_file']=str(folder/'actual_extrusion_paths.npz')
            record['support_quantity']=bambu.support_quantity(paths,row['transform_4x4'],box,discovery)
            record['feature_preservation']=dict(passed=None,status='VASCULAR_ONLY_DIAGNOSTIC_NO_BOX_FEATURE_TEST')
        else:
            bambu.audit_slice(record,cfg,row,mesh,partition,box)
            with np.load(record['actual_paths_file']) as data:paths={k:data[k] for k in data.files}
        record['contact']=support_contact(paths,mesh,partition,row['transform_4x4'],cfg)
        record['warnings'],record['critical_warnings']=archived_warnings(record['archive'])
        record['critical_warnings']+=bambu.native_warning_lines(record)
        record['evaluation_profile_sha256']=profile_hashes(record)
        q=record['support_quantity'];volume=q['total_support_extrusion_filament_mm']*np.pi*(q['filament_diameter_mm']/2)**2
        q.update(total_support_extrusion_mm3=volume,total_support_mass_g=volume/1000*q['filament_density_g_cm3'])
    on=[r for r in result['results'] if r['support_mode']=='ON' and r['status']=='BAMBU_SLICED']
    result['selected']=on[0] if on else None
    save_json(saved,result);return result


def top_visibility(mesh,face,box,cfg):
    from .support_access_qc import line_visible
    gap=cfg['support_removal_probe']['surface_standoff_mm'];p=mesh.triangles_center[face]+mesh.face_normals[face]*gap
    height=box['inner'][1,2]+gap-p[2]
    if height<=0:return False
    for angle in cfg['support']['cone_angles_deg']:
        for az in ([0] if angle==0 else np.arange(0,360,cfg['support']['cone_azimuth_step_deg'])):
            shift=height*np.tan(np.radians(angle));q=p+np.array([shift*np.cos(np.radians(az)),shift*np.sin(np.radians(az)),height])
            if np.all(q[:2]>=box['inner'][0,:2]) and np.all(q[:2]<=box['inner'][1,:2]) and line_visible(box['mesh'],mesh,p,q):return True
    return False


def trapped(visible,probe):return not visible and not probe


def assess_regions(mesh,box,regions,rotation,threshold,cfg,cache,manager):
    """Legacy area/3mm capsule machinery plus explicit small-angle visibility."""
    result=access.assess_regions(mesh,box,regions,rotation,threshold,cfg,cache,manager)
    lookup={r['region_id']:r for r in regions}
    for row in result['rows']:
        visible=row['vertical_visible_to_top']
        if not visible and not row['probe_accessible_to_top']:
            ids=np.array(lookup[row['region_id']]['face_ids']);ids=ids[mesh.face_normals[ids]@rotation[2]<-np.cos(np.radians(threshold['angle_from_horizontal_deg']))]
            for face in access.representative_faces(mesh,ids,cfg['support']['sample_faces_per_region']):
                cached=cache[int(face)] if int(face) in cache else {}
                if 'small_angle_visible' not in cached:cached['small_angle_visible']=top_visibility(mesh,int(face),box,cfg)
                if cached['small_angle_visible']:visible=True;break
        row['visible_to_top_with_small_angle']=visible
        row['status']='TRAPPED_SUPPORT_RISK' if trapped(visible,row['probe_accessible_to_top']) else 'GEOMETRY_PROXY_SUPPORT_REMOVAL_PASS'
    result['trapped_risk_count']=sum(r['status']=='TRAPPED_SUPPORT_RISK' for r in result['rows'])
    result['top_accessible_count']=sum(r['status']!='TRAPPED_SUPPORT_RISK' for r in result['rows'])
    result['top_accessible_fraction']=result['top_accessible_count']/result['candidate_regions'] if result['candidate_regions'] else 1.
    result['status']='TRAPPED_SUPPORT_RISK' if result['trapped_risk_count'] else 'GEOMETRY_PROXY_SUPPORT_REMOVAL_PASS'
    return result


def orientation_rank(row):
    s=row['support'];total=s['internal_support_risk_area_mm2']
    vascular=sum(r['overhang_area_mm2'] for r in s['rows'] if not r['vascular_branch'].startswith('BOX'))
    return [row['trapped_risk_count']!=0,row['trapped_risk_count'],vascular,total,-s['top_accessible_fraction'],
        row['downward_overhang_area_mm2'],row['z_height_mm'],row['bed_footprint_mm2'],row['candidate_id']]


def route_feedback(geometry,assessment,source_sha):
    """Exclude a trapped artificial port's current wall, never alter anatomy."""
    trapped_rows=[r for r in assessment['rows'] if r['status']=='TRAPPED_SUPPORT_RISK']
    ports={r['port_id']:r for r in geometry['ports']}
    previous=geometry.get('route_review') or {}
    allowed=copy.deepcopy(previous.get('allowed_faces_by_endpoint',{}));changed=False
    for port_id in sorted({r['vascular_branch'] for r in trapped_rows}&set(ports)):
        faces=allowed.get(port_id,['-X','+X','-Y','+Y'])
        remaining=[face for face in faces if face!=ports[port_id]['face']]
        if not remaining:return None
        allowed[port_id]=remaining;changed=True
    if not changed:return None
    return dict(vascular_sha256=source_sha,allowed_faces_by_endpoint=allowed,
        reason='TOP_ONLY_SUPPORT_QC_ROUTE_REFINEMENT',previous_upright_trapped_regions=trapped_rows,
        scope='THIS_VASCULAR_CANDIDATE_ONLY',central_geometry_changed=False,
        absolute_minimum_clearance_mm=1.6)


def refine_port_access(assembly,directory,transform,cfg,discovery,logger):
    """At most two new side-wall searches driven by actual top-access failures."""
    directory=Path(directory)
    for attempt in range(3):
        mesh=assembly['union'];regions=access.make_regions(mesh,assembly['partition'],cfg)
        manager=trimesh.collision.CollisionManager();manager.add_object('mold',mesh)
        check=assess_regions(mesh,assembly['box'],regions,np.eye(3),access.overhang_threshold(discovery,cfg),cfg,{},manager)
        save_json(directory/'manufacturing/upright_port_access_qc.json',dict(attempt=attempt,assessment=check))
        if not check['trapped_risk_count'] or attempt==2:return assembly
        feedback=route_feedback(assembly['record'],check,sha256(directory/'vascular_core.stl'))
        if feedback is None:return assembly
        save_json(directory/'port_route_review.json',feedback)
        logger.info('%s rerouting artificial ports after top-access QC: %s',directory.name,feedback['allowed_faces_by_endpoint'])
        assembly=design(directory,transform,cfg,logger)
    return assembly


def evaluate_orientations(design_data,discovery,cfg,out,logger):
    out=Path(out);out.mkdir(exist_ok=True)
    mesh,box,partition=[design_data[k] for k in ('union','box','partition')]
    threshold=access.overhang_threshold(discovery,cfg);regions=access.make_regions(mesh,partition,cfg)
    save_json(out/'support_regions.json',regions)
    rotations=access.rotations(cfg);manager=trimesh.collision.CollisionManager();manager.add_object('mold',mesh)
    cache={};gravity={};rows=[];build=discovery.get('build_volume_mm',cfg['printer']['provisional_build_volume_mm'])
    stored=out/'orientation_candidates.json'
    if stored.exists():rows=json.loads(stored.read_text())
    else:
        for index,(name,rotation) in enumerate(rotations):
            key=tuple(np.round(rotation[2],8))
            if key not in gravity:gravity[key]=assess_regions(mesh,box,regions,rotation,threshold,cfg,cache,manager)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',UserWarning)
                row=access.orientation_row(index,name,rotation,mesh,box,build,gravity[key],cfg)
            row['lexicographic_rank']=orientation_rank(row);rows.append(row)
            logger.info('Orientation %s: traps=%s internal_area=%.3f stable=%s',index,row['trapped_risk_count'],row['internal_support_risk_area_mm2'],row['stable'])
        save_json(stored,rows)
    placed=sorted([r for r in rows if r['geometry_placement_pass']],key=lambda r:r['lexicographic_rank'])
    if not placed:raise ValueError('NO_BUILD_PLACEMENT_ORIENTATION')
    top=placed[:3];upright=rows[0]
    requests=[upright]+[r for r in top if r['candidate_id']!=0]
    records=[]
    for rank,row in enumerate(requests):
        row['rank']=rank;path=out/f"orientation_{row['candidate_id']:02}.stl"
        transformed=mesh.copy();transformed.apply_transform(row['transform_4x4']);export_checked_stl(transformed,path,cfg);row['stl']=str(path)
        response=slice_cached(discovery,cfg,row,out/f"slice_{row['candidate_id']:02}",logger,mesh,partition,box)
        if response['selected']:records.append(dict(orientation=row,slice=response['selected']))
    def rank_actual(r):
        row,record=r['orientation'],r['slice'];q=record['support_quantity'];s=row['support']
        # Zero trapped risk is hard. Among zero-trap actual results, economical
        # internal support precedes fraction/access aesthetics.
        return [not row['passed'],bool(record['critical_warnings']),not record['feature_preservation']['passed'],
            q['estimated_internal_support_mass_g'],record['contact']['vascular_contact_area_mm2'],
            -s['top_accessible_fraction'],q['total_support_mass_g'],row['z_height_mm'],row['bed_footprint_mm2']]
    selected=min(records,key=rank_actual) if records else dict(orientation=top[0],slice=None)
    result=dict(selected=selected,top3=top,actual_slices=records,orientation_count=len(rows),mandatory_upright_first=True,
        proxy_ranking=['zero trapped','vascular triangle contact','internal support area','top access','total support','height','footprint'],
        actual_zero_trap_ranking=['internal support mass','vascular contact','top access','total mass','height','footprint'])
    save_json(out/'orientation_evaluation.json',result)
    write_csv(out/'orientation_candidates.csv',[{k:v for k,v in r.items() if k not in ('support',)} for r in rows])
    return result


def baseline(out=DEFAULT_OUTPUT):
    out=Path(out).resolve();cfg=configuration();discovery=fixed_discovery(out,cfg)
    logger=make_logger(out);folder=out/'baseline';saved=folder/'legacy_70mm_slice_metrics.json'
    if saved.exists():return json.loads(saved.read_text())
    oldcfg=mold.load_config(ROOT/cfg['baseline_config'])
    original=RMCA/'compact_manufacturing_roi'/oldcfg['output']['directory']
    inputs=mold.load_inputs(oldcfg);box=mold.build_intact_box(inputs,oldcfg)
    summary=json.loads((original/'abs_casting_mold_summary.json').read_text())
    source=Path(summary['final_stl']);mesh=trimesh.load_mesh(source)
    partition=qc.partition_vascular(mesh,inputs,box,cfg)
    row=copy.deepcopy(summary['orientation']['selected']);row['rank']=0
    row['stl']=str(Path(summary['final_oriented_stl']))
    result=slice_cached(discovery,cfg,row,folder/'fresh_current_orientation',logger,mesh,partition,box)
    selected=result['selected'];legacy_dims=metrics.dimensions(inputs['central'].vertices)
    record=dict(name='LEGACY_70MM_BASELINE',source_stl=str(source),source_sha256=sha256(source),oriented_stl=row['stl'],
        oriented_sha256=sha256(row['stl']),real_resliced=True,orientation=row,current_box_thickness_mm=box['wall_thickness_mm'],
        vascular_dimensions=legacy_dims,box_dimensions=(box['outer'][1]-box['outer'][0]).tolist(),selected_slice=selected,
        vascular_surface_support=metrics.surface_support(inputs['central'],inputs['companion']['graph'],inputs['companion']['transform'],np.array(row['transform_4x4'])[:3,:3])[0],
        branch_count=5,outlet_count=3,bifurcation_count=2,
        baseline_support_mass_g=selected['support_quantity']['total_support_mass_g'] if selected else None,
        baseline_internal_support_g=selected['support_quantity']['estimated_internal_support_mass_g'] if selected else None,
        baseline_total_filament_g=selected.get('filament_used_g') if selected else None,
        baseline_print_time_min=selected.get('estimated_print_time_seconds')/60 if selected else None,
        baseline_support_extrusion_mm3=selected['support_quantity']['total_support_extrusion_mm3'] if selected else None,
        baseline_warning_count=len(selected['warnings']) if selected else None,
        baseline_support_contact_area=selected['contact']['vascular_contact_area_mm2'] if selected else None)
    save_json(saved,record);return record


def upright_baseline(out=DEFAULT_OUTPUT):
    """Supplemental same-pose comparison, separating pose gains from ROI gains."""
    from .print_orientation import bed_transform
    out=Path(out).resolve();cfg=configuration();saved=out/'baseline/legacy_70mm_upright_metrics.json'
    if saved.exists():return json.loads(saved.read_text())
    discovery=fixed_discovery(out,cfg);legacy=baseline(out)
    oldcfg=mold.load_config(ROOT/cfg['baseline_config']);inputs=mold.load_inputs(oldcfg);box=mold.build_intact_box(inputs,oldcfg)
    mesh=trimesh.load_mesh(legacy['source_stl']);partition=qc.partition_vascular(mesh,inputs,box,cfg)
    mat,_,_,_=bed_transform(mesh.vertices,np.eye(3),discovery.get('build_volume_mm',[256,256,250]),10.)
    placed=mesh.copy();placed.apply_transform(mat);path=out/'baseline/legacy_upright.stl';export_checked_stl(placed,path,cfg)
    row=dict(candidate_id=0,rank=0,stl=str(path),transform_4x4=mat.tolist(),method='CASTING_UPRIGHT',
        downward_overhang_area_mm2=float(mesh.area_faces[mesh.face_normals[:,2]<-np.cos(np.radians(30))].sum()))
    result=slice_cached(discovery,cfg,row,out/'baseline/fresh_upright',make_logger(out),mesh,partition,box)
    record=dict(source_stl=legacy['source_stl'],source_sha256=legacy['source_sha256'],real_resliced=True,
        orientation=row,selected_slice=result['selected'],role='SUPPLEMENTAL_SAME_UPRIGHT_POSE_COMPARISON')
    save_json(saved,record);return record


def make_logger(out):
    logger=logging.getLogger('support_aware');logger.setLevel(logging.INFO)
    if not logger.handlers:
        handler=logging.FileHandler(Path(out)/'run.log',encoding='utf8');handler.setFormatter(logging.Formatter('%(asctime)s %(message)s'));logger.addHandler(handler)
    return logger


def evaluate(out=DEFAULT_OUTPUT):
    out=Path(out).resolve();cfg=configuration();logger=make_logger(out)
    protected=json.loads((out/'protected_before.json').read_text())
    if not verify_snapshot(protected)['all_unchanged']:raise ValueError('PROTECTED_SOURCE_CHANGED')
    discovery=fixed_discovery(out,cfg);legacy=baseline(out)
    t=np.array(json.loads((out/'roi_casting_transform.json').read_text())['transform_4x4'])
    models=json.loads((out/'model_results.json').read_text())
    usable=[r for r in models if r.get('mesh',{}).get('passed')]
    def rank(r):
        return [not r['pdms']['passed'],not r['radius_qc']['passed'],not r['dimensions']['size_target_met'],
            r['dimensions']['size_violation_mm'],r['support']['support_contact_candidate_area_mm2'],
            r['support']['thin_vessel_support_area_mm2'],r['support']['bifurcation_support_area_mm2'],r['candidate_id']]
    ordered=sorted(usable,key=rank);top=ordered[:3]
    # Preserve the requested three-outlet preference through an actual slicing
    # comparison when its proxy cost is not materially (>25%) larger. This is
    # an explicit diversity constraint after the eight-candidate surface stage.
    candidate_info={c['candidate_id']:c for c in json.loads((out/'candidates.json').read_text())}
    three=next((r for r in ordered if candidate_info[r['candidate_id']]['outlet_count']==3),None)
    if three and all(candidate_info[r['candidate_id']]['outlet_count']==2 for r in top) and three['pdms']['passed'] and three['support']['support_contact_candidate_area_mm2']<=1.25*top[0]['support']['support_contact_candidate_area_mm2']:
        top[-1]=three
    save_json(out/'top3_selection.json',[r['candidate_id'] for r in top]);evaluated=[]
    queue=top+[r for r in ordered if r not in top]
    for r in queue:
        if len([v for v in evaluated if v.get('orientations',{}).get('selected',{}).get('slice')])>=3:break
        directory=out/'candidates'/r['candidate_id'];saved=directory/'evaluation.json'
        ensure_current_design(directory,t,cfg)
        if saved.exists():
            previous=json.loads(saved.read_text())
            if 'orientations' in previous:
                previous['geometry']=load_design(directory,t,cfg)['record']
                previous['failures']=sorted(set(r['failures']+previous['geometry']['failures']))
                save_json(saved,previous);evaluated.append(previous);continue
            audit_path=directory/'manufacturing/routing_audit.json'
            if audit_path.exists():
                audit=json.loads(audit_path.read_text())
                if audit.get('expanded_search_performed') and audit['layout']['status']=='NO_COLLISION_FREE_LAYOUT':
                    previous['routing_exhausted']=True;previous['routing_audit']=str(audit_path)
                    save_json(saved,previous);evaluated.append(previous);continue
            attempts=len(list(directory.glob('evaluation_attempt_*.json')))
            saved.rename(directory/f'evaluation_attempt_{attempts+1:02}.json')
        result=dict(candidate_id=r['candidate_id'],failures=list(r['failures']))
        try:
            # Stage 1 real vascular-only slice, independent of any box or ports.
            companion=load_companion(directory,t,cfg,align=False);central=companion['mesh'];empty_box,_=__import__('vascular_processing.support_aware_port_design',fromlist=['new_box']).new_box(central.bounds,cfg)
            samples,labels,_,_,_=metrics.sample_companion(companion['graph'],t)
            face_labels=labels[cKDTree(samples).query(central.triangles_center)[1]]
            part=dict(face_ids=np.arange(len(central.faces)),labels=face_labels,groups={k:np.flatnonzero(face_labels==k) for k in set(face_labels)})
            from .print_orientation import bed_transform
            mat,ext,lo,hi=bed_transform(central.vertices,np.eye(3),discovery.get('build_volume_mm',[256,256,250]),10)
            placed=central.copy();placed.apply_transform(mat);path=directory/'vascular_only_placed.stl';export_checked_stl(placed,path,cfg)
            row=dict(candidate_id=0,rank=0,stl=str(path),transform_4x4=mat.tolist(),downward_overhang_area_mm2=r['support']['overhang_area_mm2'])
            result['vascular_only_slice']=slice_cached(discovery,cfg,row,directory/'bambu_vascular_only',logger,central,part,empty_box,vascular_only=True)
            assembly=design(directory,t,cfg,logger)
            assembly=refine_port_access(assembly,directory,t,cfg,discovery,logger)
            result['geometry']=assembly['record'];result['failures'].extend(assembly['record']['failures'])
            result['orientations']=evaluate_orientations(assembly,discovery,cfg,directory/'orientations',logger)
        except Exception as exc:
            import traceback
            result['failures'].append(type(exc).__name__+': '+str(exc));(directory/'evaluation_error.log').write_text(traceback.format_exc());logger.exception('Candidate %s failed',r['candidate_id'])
        save_json(saved,result);evaluated.append(result)
        print('EVALUATED',r['candidate_id'],result['failures'],flush=True)
    save_json(out/'evaluations.json',evaluated)
    save_json(out/'protected_after_evaluation.json',verify_snapshot(protected))
