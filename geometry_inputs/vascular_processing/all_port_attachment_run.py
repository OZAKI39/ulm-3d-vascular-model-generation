"""Frozen-input four-port alignment and an independent human-review package.

Invoked only through s1-4; no upstream modelling, wall search or box Boolean.
"""
import copy
import json
import logging
from pathlib import Path
import shutil

import numpy as np
import trimesh
import yaml

from . import sacrificial_fixture as f
from . import surface_continuity_qc as continuity
from .port_attachment_alignment import fixed_target_first_bend, clean_boolean_attachment_slivers
from .sacrificial_fixture_review import save_json

PORTS = ('I1', 'O1', 'O2', 'O3')
WALLS = dict(I1='-Y', O1='+Y', O2='-X', O3='+Y')


def saved_routes(old):
    rows = [r for r in f.read_csv(old/'tables/port_face_candidates.csv') if r.get('selected') == 'True']
    if len(rows) != 4 or {r['endpoint_id'] for r in rows} != set(PORTS):
        raise ValueError('FROZEN_ROUTE_PROVENANCE_UNRESOLVED')
    result = {r['endpoint_id']:r for r in rows}
    if any(result[p]['face'] != WALLS[p] for p in PORTS):
        raise ValueError('FROZEN_WALL_ASSIGNMENT_MISMATCH')
    return result


def rebuild_routes(inputs, cfg, old, out, logger):
    source = inputs['mesh']
    box = f.compute_box(source.bounds, cfg)
    collision = f.VesselCollision(source, inputs['graph'], inputs['branches'], inputs['transform'], inputs['endpoints'], cfg)
    rows = saved_routes(old)
    routes, details = [], {}
    for e in inputs['endpoints']:
        row = rows[e.endpoint_id]
        target = np.array([float(row['target_'+k]) for k in 'xyz'])
        kwargs = {}
        if e.endpoint_id == 'O3':
            bend_audit = dict(method='UNCHANGED_VERIFIED_O3_NATURAL_BEND', global_route_search=False,
                wall_target_changed=False, wall_assignment_changed=False, wall_target_mm=target.tolist())
        else:
            kwargs['local_first_bend'], bend_audit = fixed_target_first_bend(e, row['face'], target, cfg)
        route = f.make_route(e, row['face'], target, box, cfg, row['route_type'], float(row['tangent_scale']),
            float(row['guide_bend_radius_mm']) if row['guide_bend_radius_mm'] else None, **kwargs)
        route.vessel_clearance = collision.distance(route)
        route.required_vessel_clearance = 1.6 if e.endpoint_id == 'O3' else 3.
        path = out/'ports'/(e.endpoint_id+'_aligned.stl')
        port_qc = f.export_checked_stl(route.mesh, path, cfg)
        exported = trimesh.load_mesh(path, process=True)
        in_memory = route.mesh
        route.mesh = exported
        section = f.attachment_section_qc(route, cfg)
        # Union uses the freshly generated high precision mesh; exported STL is
        # checked independently, and never loaded from the baseline directory.
        route.mesh = in_memory
        if not section['passed']:
            raise ValueError('ACTUAL_RING_SWEEP_SECTION_FAILED_' + e.endpoint_id)
        if route.vessel_clearance < route.required_vessel_clearance:
            raise ValueError('PORT_CLEARANCE_FAILED_' + e.endpoint_id)
        details[e.endpoint_id] = dict(cap=e.attachment_qc, local_bend=bend_audit,
            port_mesh_qc=port_qc, port_stl=str(path), port_sha256=f.sha256(path),
            section_matches_actual_ring=section, radius_mm=e.radius,
            clearance_mm=route.vessel_clearance, minimum_clearance_mm=route.required_vessel_clearance,
            minimum_bend_radius_mm=route.minimum_bend_radius, wall_target_mm=route.target.tolist(),
            outside_end_mm=route.end.tolist(), original_selected_route=row,
            wall_edge_clearance_mm=route.edge_clearance, radius_changed=False,
            actual_source_ring_used=True, role='GOOD_REFERENCE' if e.endpoint_id == 'O3' else 'REPAIR_TARGET')
        logger.info('%s cap=%s collar=%s bend=%.9f clearance=%.9f',e.endpoint_id,
            e.attachment_qc['cap_vertex_count'],e.attachment_qc['collar_attempts'],route.minimum_bend_radius,route.vessel_clearance)
        routes.append(route)
    pairwise = []
    for j, a in enumerate(routes):
        for b in routes[j+1:]:
            distance = f.surface_distance(a.mesh,b.mesh)
            if distance < cfg['ports']['route_clearance_mm']:
                raise ValueError('PORT_PAIR_CLEARANCE_FAILED_'+a.endpoint.endpoint_id+'_'+b.endpoint.endpoint_id)
            pairwise.append(dict(first=a.endpoint.endpoint_id,second=b.endpoint.endpoint_id,clearance_mm=distance))
    return routes, details, pairwise


def staged_union(source, routes, cfg):
    core = source.copy()
    stages = []
    for route in routes:
        core, qc = f.union_core(core,[route],cfg)
        qc.update(stage='ADD_'+route.endpoint.endpoint_id, euler_characteristic=int(core.euler_number))
        stages.append(qc)
    return core, stages


def run(config_path, output_root=None):
    cfg=f.load_config(config_path)
    if cfg['attachment_alignment']['endpoint_ids'] != list(PORTS):
        raise ValueError('ALL_FOUR_CAP_ALIGNMENTS_REQUIRED')
    raw_cfg=copy.deepcopy(cfg);raw_cfg['attachment_alignment']['endpoint_ids']=[]
    raw=f.load_inputs(raw_cfg)
    old=raw['base']/cfg['all_port_alignment']['baseline_directory']
    out=Path(output_root).resolve() if output_root else raw['base']/cfg['all_port_alignment']['output_directory']
    if out.exists():
        raise FileExistsError('Refusing to overwrite review package: '+str(out))
    for forbidden in ('print_fixture_design_adaptive','print_fixture_design_aligned','print_fixture_design_o3_smooth'):
        if out.is_relative_to(raw['base']/forbidden):
            raise ValueError('PROTECTED_PREVIOUS_OUTPUT_DIRECTORY')
    out.mkdir(parents=True)
    for directory in ('baseline','ports','core','tables','QC'):
        (out/directory).mkdir()
    logger=logging.getLogger('all_port_attachment');logger.setLevel(logging.INFO)
    handler=logging.FileHandler(out/'run.log',encoding='utf-8');logger.addHandler(handler)
    logger.info('NEEDS_ADJUSTMENT_ALL_PORT_ATTACHMENTS')
    snapshot=f.protected_snapshot(raw,out)
    for p in (f.ROOT/'third_party/vascularmd').rglob('*.py'):
        snapshot[str(p.resolve())]=f.sha256(p)
    for name in ('point_predictions.csv','branch_predictions.csv'):
        for p in raw['base'].parent.rglob(name):snapshot[str(p.resolve())]=f.sha256(p)
    save_json(out/'protected_before.json',snapshot)
    (out/'resolved_config.yaml').write_text(yaml.safe_dump(cfg,allow_unicode=True,sort_keys=False),encoding='utf-8')
    summary=dict(initial_status='NEEDS_ADJUSTMENT_ALL_PORT_ATTACHMENTS',status='NEEDS_ADJUSTMENT',
        failures=[],per_port={},output_directory=str(out),source_stl=str(raw['paths']['source_stl']),
        source_sha256=f.sha256(raw['paths']['source_stl']),box_union_performed=False,bambu_slicing_performed=False,
        upstream_recomputed=False,vmtk_installed_or_upgraded=False,production_vmtk_ramp_used=False,
        radius_changed=False,source_vertices_changed=False,global_route_search=False,
        visual_review=dict(automated_images_generated=False,human_review_pending=True))
    try:
        baseline_path=old/'core/BG001_RMCA_BALANCED_core_with_ports.stl'
        baseline=trimesh.load_mesh(baseline_path,process=True)
        shutil.copy2(baseline_path,out/'baseline'/baseline_path.name)
        summary['baseline']=dict(path=str(baseline_path),sha256=f.sha256(baseline_path),
            copy_sha256=f.sha256(out/'baseline'/baseline_path.name),rebuilt=False)
        inputs=f.load_inputs(cfg)
        np.testing.assert_array_equal(raw['mesh'].vertices,inputs['mesh'].vertices)
        np.testing.assert_array_equal(raw['mesh'].faces,inputs['mesh'].faces)
        before_vertices=inputs['mesh'].vertices.copy();before_faces=inputs['mesh'].faces.copy()
        routes,details,pairs=rebuild_routes(inputs,cfg,old,out,logger)
        summary['per_port']=details;summary['pairwise_clearance']=pairs
        core,stages=staged_union(inputs['mesh'],routes,cfg)
        summary['staged_union']=stages
        save_json(out/'staged_boolean_audit.json',stages)
        f.write_csv(out/'tables/staged_boolean_audit.csv',[
            {k:q[k] for k in ('stage','connected_components','watertight','volume_mm3','faces','euler_characteristic','degenerate_faces','passed')}
            for q in stages])
        core,cleanup=clean_boolean_attachment_slivers(core,inputs['mesh'],inputs['endpoints'],cfg)
        summary['numerical_seam_cleanup']=cleanup
        save_json(out/'numerical_seam_cleanup.json',cleanup)
        final_path=out/'core/BG001_RMCA_BALANCED_core_with_ports_all_aligned.stl'
        summary['final_mesh_qc']=f.export_checked_stl(core,final_path,cfg)
        core=trimesh.load_mesh(final_path,process=True)  # QC the delivered file
        summary['final_stl']=str(final_path);summary['final_sha256']=f.sha256(final_path)
        missing=trimesh.boolean.difference([inputs['mesh'],core],engine='manifold')
        lost=abs(float(missing.volume)) if len(missing.faces) else 0.
        summary['source_volume_loss_mm3']=lost
        if lost>cfg['geometry']['volume_tolerance_mm3']:
            raise ValueError('SOURCE_VOLUME_PRESERVATION_FAILED')
        np.testing.assert_array_equal(before_vertices,inputs['mesh'].vertices)
        np.testing.assert_array_equal(before_faces,inputs['mesh'].faces)
        axes=[];rows=[]
        for e in inputs['endpoints']:
            before=continuity.port_continuity(baseline,e,cfg)
            after=continuity.port_continuity(core,e,cfg)
            detail=details[e.endpoint_id]
            gate=continuity.all_port_gate(e.endpoint_id,before,after,detail['clearance_mm'])
            detail.update(baseline_continuity=before,new_continuity=after,gate=gate,
                status='PASS' if gate['passed'] else 'NEEDS_ADJUSTMENT',
                vmtk_status='PENDING_VISUAL_ESCALATION_CHECK' if gate['requires_visual_escalation_check'] else 'VMTK_RAMP_NOT_NEEDED_AT_THIS_STAGE')
            summary['failures'].extend(gate['failures'])
            cap=detail['cap'];old20,new20=before['features'][1],after['features'][1]
            row=dict(port_id=e.endpoint_id,baseline_feature_20_length_mm=old20['total_length_mm'],new_feature_20_length_mm=new20['total_length_mm'],
                baseline_closed_ring_count=old20['circumferential_ring_count'],new_closed_ring_count=new20['circumferential_ring_count'],
                baseline_sharp_circumference_fraction=old20['sharp_feature_circumference_fraction'],new_sharp_circumference_fraction=new20['sharp_feature_circumference_fraction'],
                baseline_normal_p95_deg=before['normal_jumps']['p95_deg'],new_normal_p95_deg=after['normal_jumps']['p95_deg'],
                baseline_normal_max_deg=before['normal_jumps']['maximum_deg'],new_normal_max_deg=after['normal_jumps']['maximum_deg'],
                max_cross_section_area_jump=after['maximum_area_jump_fraction'],cap_tangent_mismatch_deg=cap['before_axis_error_deg'],
                collar_length_mm=cap['collar_length_mm'],collar_outside_volume_mm3=cap['collar_outside_source_volume_mm3'],
                vessel_clearance_mm=detail['clearance_mm'],status=detail['status'])
            rows.append(row)
            axes.append(dict(port_id=e.endpoint_id,**{'mean_tangent_'+k:v for k,v in zip('xyz',cap['reference_swc_mean_tangent'])},
                **{'cap_normal_'+k:v for k,v in zip('xyz',cap['cap_normal'])},angle_deg=cap['before_axis_error_deg']))
            f.write_csv(out/'tables'/(e.endpoint_id+'_cross_section_profile.csv'),
                [dict(version=version,**r) for version,metrics in [('baseline',before),('aligned',after)] for r in metrics['cross_sections']])
            f.write_csv(out/'tables'/(e.endpoint_id+'_feature_edges.csv'),
                [dict(version=version,**r) for version,metrics in [('baseline',before),('aligned',after)] for r in metrics['features']])
            logger.info('%s continuity %s',e.endpoint_id,json.dumps(row))
        f.write_csv(out/'tables/all_port_continuity_summary.csv',rows)
        f.write_csv(out/'tables/port_axis_alignment.csv',axes)
        from .all_port_attachment_review import render_review
        summary['visual_review'].update(render_review(out,baseline,core,inputs,routes,cfg,summary))
        summary['visual_review']['automated_images_generated']=True
        if any(d['gate']['requires_visual_escalation_check'] for d in details.values()):
            summary['failures'].append('VISUAL_ESCALATION_REVIEW_REQUIRED')
        summary['vmtk_status']='VMTK_RAMP_NOT_NEEDED_AT_THIS_STAGE' if not any(d['gate']['requires_visual_escalation_check'] for d in details.values()) else 'PENDING_VISUAL_ESCALATION_CHECK'
    except Exception as exc:
        logger.exception('All-port review run failed')
        summary['failures'].append(type(exc).__name__+': '+str(exc))
    protection=f.verify_snapshot(snapshot)
    save_json(out/'protected_source_hashes.json',protection)
    summary['protected_file_count']=len(snapshot);summary['protected_sources_unchanged']=protection['all_unchanged']
    if not protection['all_unchanged']:summary['failures'].append('PROTECTED_SOURCE_CHANGED')
    if not summary['failures']:summary['status']='READY_FOR_HUMAN_REVIEW'
    save_json(out/'all_port_attachment_summary.json',summary)
    save_json(out/'validation_results.json',dict(status=summary['status'],failures=summary['failures'],
        protected_sources_unchanged=protection['all_unchanged'],geometry=summary.get('final_mesh_qc'),
        per_port={k:v.get('gate') for k,v in summary['per_port'].items()},
        test_results='test_results.xml',human_review_pending=True))
    from .all_port_attachment_review import write_report
    write_report(out,summary)
    logger.info('Final status %s; no box union or slicing',summary['status'])
    handler.close();logger.removeHandler(handler)
    return summary
