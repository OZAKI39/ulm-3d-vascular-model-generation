#!/usr/bin/env python3
"""Design a reviewable two-body casting fixture from frozen BALANCED artifacts."""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import logging
from pathlib import Path
import subprocess
import sys

import cadquery as cq
import numpy as np
import trimesh
import yaml

from vascular_processing.boundary_review_export import protect_inputs
from vascular_processing import sacrificial_fixture as f
from vascular_processing.sacrificial_fixture_review import (
    assembly_export, render_review, save_json, write_report)


def probe_vmtk(cfg):
    command = [cfg['vmtk']['probe_python'], '-c',
               'import sys; from vmtk import vtkvmtk; print(sys.version); print(vtkvmtk.vtkvmtkPolyDataFlowExtensionsFilter)']
    result = dict(command=command, available=False)
    try:
        run = subprocess.run(command, capture_output=True, text=True, timeout=cfg['vmtk']['timeout_seconds'])
        result.update(available=run.returncode == 0, returncode=run.returncode, stdout=run.stdout.strip(), stderr=run.stderr.strip())
    except (OSError, subprocess.TimeoutExpired) as exc:
        result['error'] = str(exc)
    result['decision'] = ('现有 VMTK 可调用，但已验收 STL 没有开放边界，不能直接作 flow extension；采用 CadQuery 局部延长，避免切开或重建原主体'
                          if result['available'] else '当前 VMTK 无法调用；采用 CadQuery 局部延长，不重建原主体')
    return result


def endpoint_rows(endpoints, routes):
    lookup = {r.endpoint.endpoint_id: r for r in routes}
    rows = []
    for e in endpoints:
        route = lookup.get(e.endpoint_id)
        row = dict(endpoint_id=e.endpoint_id, role=e.role, original_swc_id=e.original_swc_id,
                   fitted_swc_id=e.swc_id, branch_id=e.branch_id,
                   x=e.position[0], y=e.position[1], z=e.position[2], radius_mm=e.radius, diameter_mm=2*e.radius,
                   tangent_x=e.tangent[0], tangent_y=e.tangent[1], tangent_z=e.tangent[2],
                   tangent_window_mm=e.tangent_window_mm, sample_count=e.sample_count,
                   tangent_source='frozen_cap_normal' if e.cap_ring is not None else 'swc_arc_length_regression',
                   attachment_axis_correction_deg=e.attachment_qc.get('before_axis_error_deg', 0.0),
                   original_radius_mm=e.original_radius_mm, compensated_radius_mm=e.compensated_radius_mm,
                   radius_source='frozen_native_fitted_swc_and_surface_cap_manifest',
                   direction_meaning='geometric_outward_not_measured_flow',
                   selected_box_face=route.face if route else '',
                   selected_boundary_kind=f.boundary_kind(route.face) if route else '',
                   required_vessel_clearance_mm=route.required_vessel_clearance if route else '',
                   clearance_exception_used=bool(route and route.required_vessel_clearance < route.preferred_vessel_clearance),
                   wall_target_x=route.target[0] if route else '', wall_target_y=route.target[1] if route else '',
                   wall_target_z=route.target[2] if route else '',
                   routing_status='GEOMETRY_PASS' if route else 'PORT_LAYOUT_NEEDS_MANUAL_REVIEW')
        rows.append(row)
    return rows


def run(config_path, output_override=None):
    cfg = f.load_config(config_path)
    inputs = f.load_inputs(cfg)
    out = (Path(output_override).resolve() if output_override else inputs['base'] / cfg['output']['directory'])
    # Never let a CLI path redirect generated geometry into any source folder.
    if out == inputs['base'] or inputs['base'].is_relative_to(out) or any(p.is_relative_to(out) for p in inputs['paths'].values()):
        raise ValueError('Output must be a distinct child directory, never the input directory')
    out.mkdir(parents=True, exist_ok=True)
    for folder in ('core', 'box', 'assembly', 'tables', 'QC'):
        (out / folder).mkdir(exist_ok=True)
    basename = cfg['output']['basename']
    # Remove only this entry's known products, so reruns cannot leave a stale
    # complete-looking STL/STEP/3MF after an incomplete design.
    for suffix in ('', '_PARTIAL_REVIEW_ONLY'):
        for folder, tail in [('core', f'core_with_ports{suffix}.stl'), ('box', f'casting_box{suffix}.stl'),
                             ('box', f'casting_box{suffix}.step')]:
            (out / folder / f'{basename}_{tail}').unlink(missing_ok=True)
    (out / 'assembly' / f'{basename}_assembly.3mf').unlink(missing_ok=True)
    logger = logging.getLogger('sacrificial_fixture')
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    for handler in (logging.StreamHandler(sys.stdout), logging.FileHandler(out / 'design_run.log', mode='w', encoding='utf-8')):
        handler.setFormatter(logging.Formatter('%(message)s'))
        logger.addHandler(handler)
    (out / 'resolved_design_config.yaml').write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True), encoding='utf-8')
    snapshot = f.protected_snapshot(inputs, out)
    save_json(out / 'protected_source_hashes.json', {'before_sha256': snapshot})
    logger.info('[INPUT]\nFrozen source: %s\nProtected files: %s', inputs['paths']['source_stl'], len(snapshot))
    endpoints = inputs['endpoints']
    source = inputs['mesh']
    box = f.compute_box(source.bounds, cfg)
    logger.info('[ENDPOINTS]\n%s inlet, %s outlets', sum(e.role == 'inlet' for e in endpoints), sum(e.role == 'outlet' for e in endpoints))
    logger.info('[BOX]\nInner mm: %s\nOuter mm: %s', box.inner_size, box.outer_size)
    qc = dict(status='NEEDS_ADJUSTMENT', failures=[], warnings=[
        '穿墙孔留有径向 0.20 mm 装配余量，浇注前需要确认定位和防漏方式。',
        '当前固定打印坐标保持不变；两个部件尚未进行新一轮切片或打印误差验证。',
        '端口方向来自中心线或其已确认端面的几何方向，不据此推断生理血流。'],
        inputs={k: str(v) for k, v in inputs['paths'].items()},
        source_sha256={k: snapshot[str(v.resolve())] for k, v in inputs['paths'].items()},
        box=dict(inner_bounds_mm=box.inner, outer_bounds_mm=box.outer, inner_dimensions_mm=box.inner_size,
                 outer_dimensions_mm=box.outer_size, wall_thickness_mm=cfg['box']['wall_thickness_mm'],
                 bottom_thickness_mm=cfg['box']['bottom_thickness_mm'], open_top_design_expected=True,
                 open_mesh_boundary_expected=False),
        source_core=f.mesh_qc(source, cfg),
        inlet_count=sum(e.role == 'inlet' for e in endpoints), outlet_count=sum(e.role == 'outlet' for e in endpoints),
        runtime={p: importlib.metadata.version(p) for p in ('cadquery', 'cadquery-ocp', 'trimesh', 'manifold3d', 'python-fcl', 'vtk', 'pyvista', 'shapely')},
        implementation_sha256={str(p): f.sha256(p) for p in [Path(__file__).resolve(),
            f.ROOT / 'vascular_processing/sacrificial_fixture.py', f.ROOT / 'vascular_processing/sacrificial_fixture_review.py',
            Path(config_path).resolve(), f.ROOT / 'requirements-fixture.txt']},
        vmtk=probe_vmtk(cfg),
        upstream_recomputed=False, global_rotation_changed=False, vascular_radius_redesigned=False,
        allowed_faces=list(f.enabled_faces(cfg)), top_boundary_interpretation='OPEN_TOP_WITHOUT_LID',
        assembly_clearance_definition='RADIAL', distance_method='FCL triangle-to-triangle; mesh tessellation safety allowance applied',
        curvature_method='OCC curvature evaluated on arc-length samples; not a formal continuous curvature proof')
    qc['attachment_alignment'] = {e.endpoint_id: e.attachment_qc for e in endpoints if e.cap_ring is not None}
    routes, candidates, core, box_mesh = [], [], None, None
    basename = cfg['output']['basename']
    with protect_inputs(snapshot):
        try:
            collision = f.VesselCollision(source, inputs['graph'], inputs['branches'], inputs['transform'], endpoints, cfg)
            qc['attachment_exempt_faces'] = collision.exempt_face_counts
            qc['endpoint_start_clearance'] = collision.start_checks
            option_lists = []
            logger.info('[ROUTING]')
            for endpoint in endpoints:
                options, rows = f.generate_options(endpoint, box, cfg, collision)
                option_lists.append(options)
                candidates.extend(rows)
                logger.info('%s: %s candidate routes; %s retained feasible options', endpoint.endpoint_id, len(rows), len(options))
                f.write_csv(out / 'tables' / 'port_face_candidates.csv', candidates)
            routes, qc['strict_layout'] = f.select_layout(option_lists, cfg)
            qc['layout'] = qc['strict_layout']
            if not routes and cfg.get('native_restoration'):
                logger.info('[ORIGINAL VESSEL RESTORATION AUDIT]')
                qc['native_restoration'] = f.audit_native_restoration(inputs, cfg, collision.start_checks)
                save_json(out / 'native_restoration_audit.json', qc['native_restoration'])
                for row in qc['native_restoration']['endpoints']:
                    if row['restoration_needed']:
                        logger.info('%s / raw node %s: %s; available native arc %.4f mm',
                            row['endpoint_id'], row['original_swc_id'], row['status'], row['available_arc_length_mm'])
                routes, qc['clearance_fallback'], rows = f.try_clearance_fallback(
                    endpoints, option_lists, box, cfg, collision, qc['native_restoration'])
                candidates.extend(rows)
                save_json(out / 'clearance_fallback_audit.json', qc['clearance_fallback'])
                if routes:
                    qc['layout'] = qc['clearance_fallback']['layout']
                    for endpoint_id in qc['clearance_fallback']['affected_endpoints']:
                        chosen = next(r for r in routes if r.endpoint.endpoint_id == endpoint_id)
                        qc['warnings'].append(f'{endpoint_id} 已核实没有原始下游血管可补。仅此端口至原血管的间距下限从 '
                            f"{cfg['ports']['route_clearance_mm']:g} mm 逐级放宽至 {chosen.required_vessel_clearance:g} mm；"
                            f'完整延长实体的实测最小间距为 {chosen.vessel_clearance:.4f} mm，未接受相交。')
                    logger.info('Last-resort clearance threshold: %.3f mm; affected endpoints: %s',
                        qc['clearance_fallback']['selected_clearance_mm'], qc['clearance_fallback']['affected_endpoints'])
            if not routes:
                qc['failures'].append('PORT_LAYOUT_NEEDS_MANUAL_REVIEW: ' + qc['layout']['reason'])
                feasible_lists = [options for options in option_lists if options]
                if feasible_lists:
                    routes, qc['partial_layout'] = f.select_layout(feasible_lists, cfg)
                    qc['partial_layout']['review_only'] = True
            if routes:
                selected = {r.option_id: r for r in routes}
                for row in candidates:
                    route = selected.get(row['option_id'])
                    other_routes = [r for r in routes if r.endpoint.endpoint_id != row['endpoint_id']]
                    same_wall = [r for r in other_routes if r.face == row['face']]
                    target = np.array([row['target_x'], row['target_y'], row['target_z']])
                    row['same_wall_spacing_mm'] = min((float(np.linalg.norm(target - r.target)) for r in same_wall), default=None)
                    if route:
                        row.update(selected=True, status='SELECTED_GEOMETRY_PASS', nearest_other_port_mm=route.nearest_other_port,
                                   same_wall_spacing_mm=route.wall_spacing,
                                   required_vessel_clearance_mm=route.required_vessel_clearance)
                for r in routes:
                    logger.info('%s -> %s: %s; %.2f mm; vessel clearance %.3f mm; port clearance %.3f mm',
                                r.endpoint.endpoint_id, r.face, r.kind, r.length, r.vessel_clearance, r.nearest_other_port)
                    if r.endpoint.cap_ring is not None:
                        check = f.attachment_section_qc(r, cfg)
                        qc['attachment_alignment'][r.endpoint.endpoint_id]['generated_mesh_section'] = check
                        if not check['passed']:
                            raise ValueError('ATTACHMENT_ALIGNMENT_FAILED: ' + r.endpoint.endpoint_id)
                        logger.info('%s attachment: axis error %.8f deg; contour error %.8g mm',
                            r.endpoint.endpoint_id, check['axis_error_deg'], check['contour_hausdorff_mm'])
                profiles = f.clearance_profile(routes, collision, cfg)
                f.write_csv(out / 'tables' / 'route_clearance_profile.csv', profiles)
                qc['clearance_profiles'] = {r.endpoint.endpoint_id: dict(
                    required_vessel_clearance_mm=r.required_vessel_clearance,
                    minimum_complete_mesh_clearance_mm=r.vessel_clearance,
                    minimum_outward_sampled_gap_mm=min(p['cross_section_gap_mm'] for p in profiles if p['endpoint_id'] == r.endpoint.endpoint_id),
                    last_outward_sample_below_preferred_mm=max((p['arc_from_original_endpoint_mm'] for p in profiles
                        if p['endpoint_id'] == r.endpoint.endpoint_id and p['below_preferred_clearance']), default=None),
                    sample_spacing_max_mm=cfg['routing']['curvature_step_mm'],
                    note='Sampled cross-sections are diagnostic; acceptance uses full triangle surface including the attachment overlap') for r in routes}
                logger.info('[BOOLEAN]')
                core, qc['core'] = f.union_core(source, routes, cfg)
                suffix = '' if len(routes) == len(endpoints) else '_PARTIAL_REVIEW_ONLY'
                qc['core']['complete_endpoint_set'] = len(routes) == len(endpoints)
                qc['core']['stl'] = str(out / 'core' / f'{basename}_core_with_ports{suffix}.stl')
                qc['core']['exported_stl_qc'] = f.export_checked_stl(core, qc['core']['stl'], cfg)
                if min(r.edge_clearance for r in routes) - cfg['ports']['min_wall_edge_clearance_mm'] < cfg['ports']['assembly_clearance_mm']:
                    qc['warnings'].append(f"入口/出口中最小盒壁边距为 {min(r.edge_clearance for r in routes):.3f} mm，接近 {cfg['ports']['min_wall_edge_clearance_mm']:g} mm 下限，人工审核时需留意实际打印公差。")
                logger.info('Core union PASS: one connected closed positive-volume body')
            body, holes = f.create_box(box, routes, cfg)
            box_mesh = f.shape_mesh(body, cfg)
            top_crossings = [dict(port_id=r.endpoint.endpoint_id, face=r.face, center_mm=r.target.tolist(),
                outside_stub_mm=cfg['ports']['outside_stub_mm'], passage_kind='OPEN_TOP_NO_HOLE')
                for r in routes if f.boundary_kind(r.face) == 'OPEN_TOP']
            material_ports = [r for r in routes if f.boundary_kind(r.face) != 'OPEN_TOP']
            qc['box'].update(mesh=f.mesh_qc(box_mesh, cfg), holes=holes, port_hole_count=len(holes),
                measured_dimensions=f.box_dimension_checks(body, box, cfg),
                expected_port_hole_count=len(material_ports), hole_count_matches_assigned_material_ports=len(holes) == len(material_ports),
                open_top_crossings=top_crossings, open_top_crossing_count=len(top_crossings),
                required_endpoint_count=len(endpoints), boundary_passage_count=len(holes) + len(top_crossings),
                all_endpoints_have_boundary_passage=len(holes) + len(top_crossings) == len(endpoints),
                all_holes_passed=all(h['passed'] for h in holes), top_cavity_open=not body.isInside(cq.Vector(*box.inner.mean(axis=0)[:2], box.inner[1, 2] - cfg['geometry']['coordinate_tolerance_mm'])))
            suffix = '' if len(routes) == len(endpoints) else '_PARTIAL_REVIEW_ONLY'
            qc['box']['step'] = str(out / 'box' / f'{basename}_casting_box{suffix}.step')
            qc['box']['stl'] = str(out / 'box' / f'{basename}_casting_box{suffix}.stl')
            cq.exporters.export(body, qc['box']['step'])
            qc['box']['exported_stl_qc'] = f.export_checked_stl(box_mesh, qc['box']['stl'], cfg)
            if core is not None:
                overlap = trimesh.boolean.intersection([core, box_mesh], engine='manifold')
                overlap_volume = abs(float(overlap.volume)) if len(overlap.faces) else 0.0
                qc['box']['core_box_overlap_volume_mm3'] = overlap_volume
                if overlap_volume > cfg['geometry']['volume_tolerance_mm3']:
                    qc['failures'].append(f'CORE_BOX_COLLISION: intersection {overlap_volume:.6g} mm³')
            if not qc['box']['mesh']['passed'] or not qc['box']['all_holes_passed'] or not qc['box']['top_cavity_open'] or not qc['box']['measured_dimensions']['passed'] or not qc['box']['hole_count_matches_assigned_material_ports']:
                qc['failures'].append('BOX_QC_FAILED')
            if len(routes) != len(endpoints):
                qc['failures'].append('INCOMPLETE_ENDPOINT_EXTENSION')
            logger.info('Box: %s through holes, %s open-top crossings; wall mesh closed=%s', len(holes), len(top_crossings), qc['box']['mesh']['watertight'])
        except Exception as exc:
            logger.exception('Actual design failure')
            qc['failures'].append(type(exc).__name__ + ': ' + str(exc))
        if not qc['failures']:
            qc['status'] = 'READY_FOR_HUMAN_REVIEW'
        f.write_csv(out / 'tables' / 'endpoint_inventory.csv', endpoint_rows(endpoints, routes))
        if candidates:
            f.write_csv(out / 'tables' / 'port_face_candidates.csv', candidates)
        summary = [dict(port_id=r.endpoint.endpoint_id, role=r.endpoint.role, branch_id=r.endpoint.branch_id, wall=r.face,
                   boundary_kind=f.boundary_kind(r.face), boundary_thickness_mm=f.boundary_thickness(r.face, cfg),
                   route_type=r.kind, route_length_mm=r.length, turn_angle_deg=r.turn, minimum_clearance_mm=r.vessel_clearance,
                   required_vessel_clearance_mm=r.required_vessel_clearance,
                   preferred_clearance_met=r.vessel_clearance >= cfg['ports']['route_clearance_mm'],
                   minimum_bend_radius_mm=r.minimum_bend_radius, wall_edge_clearance_mm=r.edge_clearance,
                   nearest_other_port_mm=r.nearest_other_port, same_wall_center_spacing_mm=r.wall_spacing,
                   port_radius_mm=r.endpoint.radius, outside_stub_mm=cfg['ports']['outside_stub_mm'],
                   qc_status=('GEOMETRY_PASS_AUTHORIZED_CLEARANCE_EXCEPTION' if r.required_vessel_clearance < cfg['ports']['route_clearance_mm'] else 'GEOMETRY_PASS')
                   if len(routes) == len(endpoints) else 'GEOMETRY_PASS_PARTIAL_LAYOUT_ONLY') for r in routes]
        for e in endpoints:
            if e.endpoint_id not in {row['port_id'] for row in summary}:
                row = {key: '' for key in (summary[0] if summary else ['port_id', 'role', 'branch_id', 'wall', 'boundary_kind', 'boundary_thickness_mm', 'route_type',
                    'route_length_mm', 'turn_angle_deg', 'minimum_clearance_mm', 'required_vessel_clearance_mm', 'preferred_clearance_met', 'minimum_bend_radius_mm', 'wall_edge_clearance_mm',
                    'nearest_other_port_mm', 'same_wall_center_spacing_mm', 'port_radius_mm', 'outside_stub_mm', 'qc_status'])}
                row.update(port_id=e.endpoint_id, role=e.role, branch_id=e.branch_id, port_radius_mm=e.radius,
                           qc_status='PORT_LAYOUT_NEEDS_MANUAL_REVIEW',
                           minimum_clearance_mm=qc.get('endpoint_start_clearance', {}).get(e.endpoint_id, {}).get('distance_mm', ''))
                summary.append(row)
        f.write_csv(out / 'tables' / 'port_layout_summary.csv', summary)
        qc['port_layout'] = summary
        qc['assembly'] = assembly_export(out / 'assembly', basename, source, routes, box_mesh, core if not qc['failures'] else None)
        try:
            qc['png_files'] = render_review(out / 'QC', source, endpoints, routes, box, box_mesh, core, cfg, qc['status'], qc.get('endpoint_start_clearance', {}))
        except Exception as exc:
            logger.exception('Review figure failure')
            qc['failures'].append('REVIEW_FIGURES_FAILED: ' + str(exc))
            qc['status'] = 'NEEDS_ADJUSTMENT'
    protection = f.verify_snapshot(snapshot)
    save_json(out / 'protected_source_hashes.json', protection)
    qc['protection'] = {k: v for k, v in protection.items() if 'sha256' not in k}
    if not protection['all_unchanged']:
        qc['failures'].append('PROTECTED_SOURCE_CHANGED')
        qc['status'] = 'NEEDS_ADJUSTMENT'
    qc['output_directory'] = str(out)
    save_json(out / 'geometry_qc.json', qc)
    save_json(out / 'attachment_alignment_qc.json', qc['attachment_alignment'])
    write_report(out, inputs, box, endpoints, routes, cfg, qc)
    if qc['failures']:
        (out / 'failure_summary.md').write_text('# 本轮结构设计需要调整\n\n' + '\n'.join('- ' + s for s in qc['failures']) +
            '\n\n保留原模型与可生成的审核图。候选被拒绝的具体原因见 tables/port_face_candidates.csv。\n', encoding='utf-8')
    elif (out / 'failure_summary.md').exists():
        (out / 'failure_summary.md').unlink()
    logger.info('[QC]\n%s\nProtected source hashes unchanged: %s\n[OUTPUT]\n%s', qc['status'], protection['all_unchanged'], out)
    return qc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=f.ROOT / 'config/sacrificial_box_BG001.yaml')
    parser.add_argument('--output-root', type=Path, help='Optional isolated output folder for verification')
    parser.add_argument('--o3-preflight', action='store_true', help='Audit installed VMTK ramp, VTP bridge and the existing O3 seam; does not repair geometry')
    parser.add_argument('--all-ports-aligned', action='store_true', help='Rebuild four attachments against frozen wall targets and generate before/after review')
    args = parser.parse_args()
    use_aligned = args.all_ports_aligned or (not args.o3_preflight and
        f.load_config(args.config).get('all_port_alignment', {}).get('enabled', False))
    if use_aligned:
        if args.o3_preflight:
            parser.error('--all-ports-aligned and --o3-preflight are mutually exclusive')
        from vascular_processing.all_port_attachment_run import run as all_ports
        qc = all_ports(args.config, args.output_root)
        print(qc['status'])
        return 0 if qc['status'] == 'READY_FOR_HUMAN_REVIEW' else 2
    if args.o3_preflight:
        from vascular_processing.o3_smooth_preflight import run as o3_preflight
        qc = o3_preflight(args.config, args.output_root)
        print(qc['status'] + ': ' + qc.get('failure_code', ''))
        return 2
    qc = run(args.config, args.output_root)
    return 0 if qc['status'] == 'READY_FOR_HUMAN_REVIEW' else 2


if __name__ == '__main__':
    raise SystemExit(main())
