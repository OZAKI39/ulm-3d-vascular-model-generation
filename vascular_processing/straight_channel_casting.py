"""A separate straight-channel manufacturing variant; previous Y files stay intact.

All geometry is in mm. The 1.6 mm solid is a sacrificial lumen core, not a
1.6 mm outer-diameter hollow ABS pipe. The horizontal end is buried inside PDMS
for later manual cutting.
The current bottomless fixture uses continuous walls, external dovetails and
raised collar retainers. The historical floor and side-wall gate are retained
only for their archived configurations.
"""
from __future__ import annotations
import hashlib
import json
import math
import shutil
import zipfile
from pathlib import Path
import xml.etree.ElementTree as ET

import cadquery as cq
import numpy as np
import pyvista as pv
import trimesh
import yaml
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
from . import y_channel_casting as shared

ROOT = shared.ROOT
save, digest, box, union = shared.save, shared.digest, shared.box, shared.union
PART_NAMES = ['Straight_lumen_core', 'Straight_casting_box_base_2mm',
              'Straight_retainer_INLET', 'Straight_retainer_OUTLET1']
PREFIX = 'Straight_D1p6_P1S_04_ABS'


def two_top_ports(cfg):
    return cfg.get('channel_topology') == 'two_vertical_ports'


def part_names(cfg):
    if shared.bottomless(cfg):
        from .bottomless_channel_fixture import part_names as fixture_parts
        return fixture_parts('Straight')
    return PART_NAMES+(['Straight_sidewall_gate'] if cfg.get('chamber_mode', 'side_gate') == 'side_gate' else [])


def load_config(path):
    cfg = yaml.safe_load(Path(path).read_text())
    shared.validate_pdms_width_requirement(cfg)
    cfg['output'] = (ROOT/Path(cfg['output'])).resolve()
    for k, v in cfg.items():
        if k == 'base_mm' and shared.bottomless(cfg) and v == 0: continue
        if k.endswith('_mm') and (not np.isfinite(v) or v <= 0):
            raise ValueError('POSITIVE_DIMENSION_REQUIRED: '+k)
    if cfg['diameter_mm'] != 1.6: raise ValueError('THIS_VARIANT_REQUIRES_1P6_MM')
    if 'minimum_pdms_height_mm' in cfg and cfg['minimum_pdms_height_mm'] < 9.:
        raise ValueError('PDMS_HEIGHT_AT_LEAST_9_MM_REQUIRED')
    if cfg['printer']['machine'] != 'Bambu Lab P1S 0.4 nozzle': raise ValueError('P1S_04_REQUIRED')
    if not shared.bottomless(cfg) and cfg['freeboard_mm'] < 2.2: raise ValueError('RETAINER_FEET_IN_PDMS')
    if shared.bottomless(cfg):
        from .bottomless_channel_fixture import validate_config
        validate_config(cfg)
    mode = cfg.get('chamber_mode', 'side_gate')
    topology = cfg.get('channel_topology', 'three_ports')
    if topology not in ['three_ports', 'two_vertical_ports']:
        raise ValueError('UNKNOWN_STRAIGHT_TOPOLOGY')
    if two_top_ports(cfg):
        if not shared.bottomless(cfg) or mode != 'closed_wall_two_vertical_ports':
            raise ValueError('TWO_TOP_PORTS_REQUIRE_CLOSED_BOTTOMLESS_FRAME')
        if any(k in cfg for k in ['horizontal_extension_after_outlet1_mm', 'horizontal_end_clearance_mm', 'manual_trim_overlap_mm']):
            raise ValueError('TWO_TOP_PORTS_MUST_NOT_DECLARE_HORIZONTAL_OUTLET')
        return cfg
    if mode not in ['side_gate', 'closed_wall_internal_end']:
        raise ValueError('UNKNOWN_CHAMBER_MODE')
    if mode == 'side_gate':
        if cfg['side_gate_width_mm'] <= cfg['diameter_mm']+2*cfg['side_gate_fit_clearance_mm']+2:
            raise ValueError('GATE_SIDE_LANDS_TOO_THIN')
    else:
        for key in ['horizontal_extension_after_outlet1_mm', 'horizontal_end_clearance_mm', 'manual_trim_overlap_mm']:
            if key not in cfg: raise ValueError('MISSING_CLOSED_CHAMBER_DIMENSION: '+key)
        if cfg['horizontal_extension_after_outlet1_mm']-cfg['manual_trim_overlap_mm'] <= cfg['diameter_mm']/2+cfg['tee_fillet_mm']:
            raise ValueError('TRIM_PLANE_TOO_CLOSE_TO_JUNCTION')
    return cfg


def cad_mesh(shape, cfg):
    if not shape.isValid() or len(shape.Solids()) != 1: raise ValueError('INVALID_CAD_SOLID')
    vertices, faces = shape.tessellate(cfg['cad_tolerance_mm'], cfg['cad_angular_tolerance'])
    mesh = trimesh.Trimesh([p.toTuple() for p in vertices], faces, process=True)
    before = len(mesh.faces)
    # OCP emits four duplicate facets at the equal-radius tee fillet seams.
    # Remove duplicates/zero-area facets, without moving any coordinates.
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.update_faces(mesh.unique_faces())
    mesh.remove_unreferenced_vertices()
    mesh.metadata['duplicate_or_degenerate_facets_removed'] = before-len(mesh.faces)
    if not mesh.is_volume or len(mesh.split()) != 1: raise ValueError('INVALID_TESSELLATED_SOLID')
    return mesh


def cylinder_x(radius, x0, x1, z):
    mesh = trimesh.creation.cylinder(radius, x1-x0, sections=128)
    m = np.eye(4)
    m[:3, :3] = Rotation.from_euler('y', 90, degrees=True).as_matrix()
    m[:3, 3] = [(x0+x1)/2, 0, z]
    mesh.apply_transform(m)
    return mesh


def design(cfg):
    if two_top_ports(cfg):
        from .straight_two_top_ports import design as two_port_design
        return two_port_design(cfg)
    r = cfg['diameter_mm']/2
    R = r*cfg['bend_radius_to_lumen_radius']
    z = cfg['base_mm']+cfg['pdms_floor_mm']+r
    start = -cfg['parent_length_mm']
    junction = cfg['downstream_length_mm']
    bend_x = start-cfg['terminal_straight_mm']
    inlet_x = bend_x-R
    inner = np.array([[inlet_x-r-cfg['side_clearance_mm'], -r-cfg['side_clearance_mm']],
                      [junction+r+cfg['side_clearance_mm'], r+cfg['side_clearance_mm']]])
    mode = cfg.get('chamber_mode', 'side_gate')
    if mode == 'closed_wall_internal_end':
        horizontal_end = junction+cfg['horizontal_extension_after_outlet1_mm']
        inner[1, 0] = horizontal_end+cfg['horizontal_end_clearance_mm']
    inner = shared.enforce_pdms_width(inner, cfg)
    outer = inner+np.array([[-cfg['wall_mm']]*2, [cfg['wall_mm']]*2])
    # Height is measured from the INNER floor, excluding the ABS base and freeboard.
    fill = max(z+R+cfg['pdms_above_elbow_mm'],
               cfg['base_mm']+shared.minimum_cad(cfg, 'height'))
    rim = fill+cfg['freeboard_mm']
    collar_z = shared.retainer_arm_bottom(cfg, rim)+cfg['retainer_thickness_mm']
    top = collar_z+cfg['port_above_retainer_mm']
    gate_back = None
    if mode == 'side_gate':
        gate_back = outer[1, 0]+cfg['side_gate_fit_clearance_mm']+cfg['side_gate_flange_mm']
        horizontal_end = gate_back+cfg['horizontal_projection_beyond_gate_mm']
    return dict(radius=r, bend_radius=R, center_z=z, start_x=start, junction_x=junction,
        bend_x=bend_x, inlet_x=inlet_x, inner=inner, outer=outer, fill=fill, rim=rim,
        collar_z=collar_z, top=top, gate_back=gate_back, horizontal_end=horizontal_end, chamber_mode=mode)


def build_core(cfg, d):
    """Exact circular CAD sweep plus a filleted tee; no global mesh smoothing."""
    if two_top_ports(cfg):
        from .straight_two_top_ports import build_core as two_port_core
        return two_port_core(cfg,d)
    r, R, z = d['radius'], d['bend_radius'], d['center_z']
    top = cq.Vector(d['inlet_x'], 0, d['top'])
    knee = cq.Vector(d['inlet_x'], 0, z+R)
    mid = cq.Vector(d['bend_x']-R/math.sqrt(2), 0, z+R-R/math.sqrt(2))
    low = cq.Vector(d['bend_x'], 0, z)
    end = cq.Vector(d['horizontal_end'], 0, z)
    path = cq.Wire.assembleEdges([cq.Edge.makeLine(top, knee),
        cq.Edge.makeThreePointArc(knee, mid, low), cq.Edge.makeLine(low, end)])
    circle = cq.Wire.makeCircle(r, top, cq.Vector(0, 0, 1))
    base = cq.Solid.sweep(circle, [], path)
    riser = cq.Solid.makeCylinder(r, d['top']-z, cq.Vector(d['junction_x'], 0, z), cq.Vector(0, 0, 1))
    channel = base.fuse(riser).clean()
    intersection = [e for e in channel.Edges() if e.geomType() == 'ELLIPSE']
    if len(intersection) != 4: raise ValueError('EXPECTED_FOUR_TEE_INTERSECTION_EDGES')
    channel = channel.fillet(cfg['tee_fillet_mm'], intersection)
    if not channel.isValid(): raise ValueError('TEE_FILLET_FAILED')
    core = channel
    for x in [d['inlet_x'], d['junction_x']]:
        collar = cq.Solid.makeCylinder(r+cfg['collar_radial_width_mm'], cfg['collar_thickness_mm'],
                                      cq.Vector(x, 0, d['collar_z']), cq.Vector(0, 0, 1))
        core = core.fuse(collar).clean()
    return core, channel


def build_mold_and_gate(cfg, d):
    inner, outer, z, rim = d['inner'], d['outer'], d['center_z'], d['rim']
    m = trimesh.boolean.difference([box([np.r_[outer[0], 0], np.r_[outer[1], rim]]),
        box([np.r_[inner[0], cfg['base_mm']], np.r_[inner[1], rim+1]])], engine='manifold')
    if d['chamber_mode'] == 'closed_wall_internal_end':
        return m, None
    half = cfg['side_gate_width_mm']/2
    clearance = cfg['side_gate_fit_clearance_mm']
    bore_radius = d['radius']+cfg['side_gate_bore_radial_clearance_mm']
    bore = cylinder_x(bore_radius, inner[1, 0]-.1, d['gate_back']+.1, z)
    notch = box([[inner[1, 0]-.1, -half, z], [outer[1, 0]+.1, half, rim+1]])
    m = trimesh.boolean.difference([m, union([notch, bore])], engine='manifold')
    tongue = box([[inner[1, 0], -half+clearance, z],
                  [outer[1, 0]+clearance+.1, half-clearance, rim]])
    # Exterior overlapping flange covers the two long slot seams. An open
    # lower notch lets the flange pass over the rod while lowering from +Z.
    flange = box([[outer[1, 0]+clearance, -half-1.2, z-bore_radius-1.0],
                  [d['gate_back'], half+1.2, rim+1.2]])
    downslot = box([[inner[1, 0]-.1, -bore_radius, z-4],
                    [d['gate_back']+.1, bore_radius, z]])
    gate = trimesh.boolean.difference([union([tongue, flange]), union([bore, downslot])], engine='manifold')
    return m, gate


def volume_intersection(a, b):
    inter = trimesh.boolean.intersection([a, b], engine='manifold')
    return abs(float(inter.volume)) if len(inter.faces) and min(inter.extents) > 1e-8 else 0.


def section_diameter(mesh, origin, normal):
    section = mesh.section(plane_origin=origin, plane_normal=normal)
    if section is None: raise ValueError('MISSING_LUMEN_SECTION')
    plane, _ = section.to_2D()
    polygons = plane.polygons_full
    if len(polygons) != 1: raise ValueError('SECTION_NOT_SINGLE_CHANNEL')
    area = float(polygons[0].area)
    return 2*np.sqrt(area/np.pi)


def closed_chamber_checks(cfg, d, core, channel, mold):
    """Audit actual wall solids and burial, beyond merely watertight meshes.

    A tray with a side opening is still a watertight STL solid; that alone
    cannot establish a closed casting chamber. Test each entire wall slab.
    """
    inner, outer = d['inner'], d['outer']
    slabs = {'bottom': box([np.r_[outer[0], 0], np.r_[outer[1], cfg['base_mm']]])}
    for axis, name in [(0, 'x'), (1, 'y')]:
        for side in [0, 1]:
            bounds = np.array([np.r_[outer[0], 0], np.r_[outer[1], d['rim']]])
            if side == 0: bounds[1, axis] = inner[0, axis]
            else: bounds[0, axis] = inner[1, axis]
            slabs[name+('_min' if side == 0 else '_max')] = box(bounds)
    missing = {k: max(0., float(m.volume)-volume_intersection(m, mold)) for k, m in slabs.items()}
    # Slice beyond the riser/tee envelope, isolating the horizontal extension.
    extension = trimesh.intersections.slice_mesh_plane(channel, [1, 0, 0],
        [d['junction_x']+d['radius']+cfg['tee_fillet_mm'], 0, 0], cap=True)
    if not extension.is_volume: raise ValueError('HORIZONTAL_EXTENSION_NOT_CLOSED')
    fill = box([np.r_[inner[0], cfg['base_mm']], np.r_[inner[1], d['fill']]])
    lower = extension.bounds[0]-fill.bounds[0]
    upper = fill.bounds[1]-extension.bounds[1]
    outside = max(0., float(extension.volume)-volume_intersection(extension, fill))
    end_clearance = float(inner[1, 0]-extension.bounds[1, 0])
    # The complete core can be lowered through the open top, without a side gate.
    insertion = []
    for lift in np.linspace(0., d['rim']-core.bounds[0, 2]+1, 25):
        moved = core.copy(); moved.apply_translation([0, 0, lift])
        insertion.append(dict(lift_mm=float(lift), core_mold_intersection_mm3=volume_intersection(moved, mold)))
    cut_x = d['horizontal_end']-cfg['manual_trim_overlap_mm']
    return dict(wall_and_bottom_missing_volume_mm3=missing,
        horizontal_extension_bounds_mm=extension.bounds,
        extension_clearance_to_pdms_bounds_mm=dict(lower_xyz=lower, upper_xyz=upper),
        horizontal_extension_outside_pdms_mm3=outside,
        horizontal_end_to_inner_wall_mm=end_clearance, core_insertion_path=insertion,
        manual_cut_reference=dict(plane_x_mm=cut_x, normal=[1, 0, 0],
            pdms_to_remove_from_end_face_mm=float(inner[1, 0]-cut_x),
            overlap_into_horizontal_channel_mm=cfg['manual_trim_overlap_mm'],
            manual_post_cure_operation=True, cut_not_applied_to_printed_geometry=True))


def geometry(cfg, out):
    if shared.bottomless(cfg):
        from .bottomless_channel_fixture import geometry as fixture_geometry
        return fixture_geometry(cfg, out, 'Straight')
    d = design(cfg)
    core_cad, channel_cad = build_core(cfg, d)
    core, channel = cad_mesh(core_cad, cfg), cad_mesh(channel_cad, cfg)
    mold, gate = build_mold_and_gate(cfg, d)
    inlet = shared.retainer(dict(name='INLET', radius_mm=d['radius'], top_center_mm=[d['inlet_x'], 0, d['top']]),
                            d['inner'], cfg, dict(rim_z_mm=d['rim']))
    outlet = shared.retainer(dict(name='OUTLET_A', radius_mm=d['radius'], top_center_mm=[d['junction_x'], 0, d['top']]),
                             d['inner'], cfg, dict(rim_z_mm=d['rim']))
    parts = dict(zip(PART_NAMES, [core, mold, inlet, outlet]))
    if gate is not None: parts['Straight_sidewall_gate'] = gate
    (out/'assembly').mkdir(parents=True, exist_ok=True)
    for name, mesh in parts.items(): mesh.export(out/'assembly'/(name+'.stl'))
    channel.export(out/'assembly/Straight_channel_without_collars_reference.stl')
    cq.exporters.export(core_cad, str(out/'assembly/Straight_lumen_core.step'))
    points = core.triangles_center
    topfaces = np.max(np.abs(core.triangles[:, :, 2]-d['top']), axis=1) < 1e-5
    rightfaces = np.max(np.abs(core.triangles[:, :, 0]-d['horizontal_end']), axis=1) < 1e-5
    port_rows = []
    for name, center, axis, select in [
        ('INLET', [d['inlet_x'], 0, d['top']], [0, 0, 1], topfaces & (points[:, 0] < 0)),
        ('OUTLET1', [d['junction_x'], 0, d['top']], [0, 0, 1], topfaces & (points[:, 0] > 0)),
        ('OUTLET2', [d['horizontal_end'], 0, d['center_z']], [1, 0, 0], rightfaces)]:
        port_rows.append(dict(name=name, center_mm=center, outward_axis=axis,
            measured_cap_diameter_mm=float(2*np.sqrt(core.area_faces[select].sum()/np.pi)),
            minimum_normal_dot=float((core.face_normals[select]@axis).min()), cap_faces=int(select.sum())))
    sampled = [dict(x_mm=x, equivalent_diameter_mm=section_diameter(channel, [x, 0, d['center_z']], [1, 0, 0]))
               for x in [-13., -7., 0., 6., 9., 15., d['horizontal_end']-.5]]
    from .y_uniform_channel import local_section_diameter
    extra_sections = []
    for angle in np.linspace(.12, np.pi/2-.12, 5):
        extra_sections.append(dict(region='inlet_elbow',
            center_mm=[d['bend_x']-d['bend_radius']*np.sin(angle), 0,
                       d['center_z']+d['bend_radius']*(1-np.cos(angle))],
            normal=[-np.cos(angle), 0, np.sin(angle)]))
    for port, x in [('INLET', d['inlet_x']), ('OUTLET1', d['junction_x'])]:
        for height in [d['center_z']+d['bend_radius']+.5, d['fill']-.5, d['top']-.4]:
            extra_sections.append(dict(region=port+'_riser', center_mm=[x, 0, height], normal=[0, 0, 1]))
    for row in extra_sections:
        row['equivalent_diameter_mm'] = local_section_diameter(channel, row['center_mm'], row['normal'])
    names = list(parts)
    collision = {a+' / '+b: volume_intersection(parts[a], parts[b])
                 for i, a in enumerate(names) for b in names[i+1:]}
    # Prove the side gate can actually be lowered into place without crossing
    # the rod or mold; a final zero-volume intersection alone is insufficient.
    gate_path = []
    if gate is not None:
        for lift in np.linspace(0, d['rim']-d['center_z']+3, 25):
            moved = gate.copy()
            moved.apply_translation([0, 0, lift])
            gate_path.append(dict(lift_mm=float(lift), core_intersection_mm3=volume_intersection(moved, core),
                                  mold_intersection_mm3=volume_intersection(moved, mold)))
    fill = box([np.r_[d['inner'][0], cfg['base_mm']], np.r_[d['inner'][1], d['fill']]])
    net = fill.volume-volume_intersection(core, fill)
    mesh_qc = {k: shared.mesh_info(m) for k, m in parts.items()}
    width_qc = shared.pdms_width_qc(mold, d['inner'], cfg, d['fill'])
    height_qc = shared.pdms_height_qc(mold, d['inner'], cfg, d['fill'])
    gates = dict(one_solid_each=all(v['positive_volume'] and v['watertight'] and v['components'] == 1 for v in mesh_qc.values()),
        pdms_width_minimum=width_qc['passed'],
        diameter_1p6_all_three_caps=all(abs(p['measured_cap_diameter_mm']-cfg['diameter_mm']) < .003 for p in port_rows),
        diameter_1p6_straight_sections=all(abs(p['equivalent_diameter_mm']-cfg['diameter_mm']) < .003 for p in sampled),
        diameter_1p6_elbow_and_risers=all(abs(p['equivalent_diameter_mm']-cfg['diameter_mm']) < .003 for p in extra_sections),
        port_directions=all(p['minimum_normal_dot'] > .99999 for p in port_rows),
        pdms_floor=abs(core.bounds[0, 2]-cfg['base_mm']-cfg['pdms_floor_mm']) < .003,
        pdms_height_minimum=height_qc['passed'],
        no_part_penetration=max(collision.values()) < 1e-4,
        retainers_above_fill=min(inlet.bounds[0, 2], outlet.bounds[0, 2])-d['fill'] > .3)
    chamber_qc = None
    if gate is not None:
        gates['gate_insertion_clear'] = max(max(p['core_intersection_mm3'], p['mold_intersection_mm3']) for p in gate_path) < 1e-4
    else:
        chamber_qc = closed_chamber_checks(cfg, d, core, channel, mold)
        bounds_clearance = chamber_qc['extension_clearance_to_pdms_bounds_mm']
        gates.update(continuous_four_walls_and_bottom=max(chamber_qc['wall_and_bottom_missing_volume_mm3'].values()) < 1e-3,
            horizontal_extension_inside_pdms=chamber_qc['horizontal_extension_outside_pdms_mm3'] < 1e-4 and min(*bounds_clearance['lower_xyz'], *bounds_clearance['upper_xyz']) > 0,
            horizontal_end_clearance=abs(chamber_qc['horizontal_end_to_inner_wall_mm']-cfg['horizontal_end_clearance_mm']) < 1e-5,
            core_insertable_from_top=max(p['core_mold_intersection_mm3'] for p in chamber_qc['core_insertion_path']) < 1e-4)
    summary = dict(status='GEOMETRY_PASSED' if all(gates.values()) else 'GEOMETRY_FAILED', units='mm',
        topology='Straight main run, upward inlet, downstream upward outlet 1, coaxial horizontal outlet 2',
        construction='Exact CAD circle sweep and cylinders; four tee edges filleted locally',
        tessellation_cleanup=dict(core_facets_removed=core.metadata['duplicate_or_degenerate_facets_removed'],
            channel_facets_removed=channel.metadata['duplicate_or_degenerate_facets_removed'], vertex_positions_changed=False),
        baseline_straight_length_mm=cfg['parent_length_mm']+cfg['downstream_length_mm'],
        diameter_mm=cfg['diameter_mm'], pdms_minimum_floor_mm=float(core.bounds[0, 2]-cfg['base_mm']),
        mold_outer_size_mm=mold.extents, cavity_xy_size_mm=d['inner'][1]-d['inner'][0], pdms_width_qc=width_qc, pdms_height_qc=height_qc,
        pdms_fill_height_mm=d['fill']-cfg['base_mm'], pdms_net_fill_ml=net/1000,
        wall_mm=cfg['wall_mm'], base_mm=cfg['base_mm'], design=d, ports=port_rows,
        straight_section_checks=sampled, elbow_riser_section_checks=extra_sections, tee_fillet_mm=cfg['tee_fillet_mm'],
        chamber_mode=d['chamber_mode'], production_part_names=names, closed_chamber_qc=chamber_qc,
        horizontal_projection_beyond_gate_mm=cfg.get('horizontal_projection_beyond_gate_mm') if gate is not None else None,
        meshes=mesh_qc, assembly_intersections_mm3=collision, gate_insertion_path=gate_path,
        side_gate_external_seal_required=gate is not None,
        side_gate_radial_gap_mm=cfg.get('side_gate_bore_radial_clearance_mm') if gate is not None else None,
        physical_fit_validated=False, gates=gates)
    summary['assembly_sha256'] = {p.name: digest(p) for p in (out/'assembly').iterdir() if p.is_file()}
    save(out/'geometry_qc.json', summary)
    save(out/'effective_config.json', cfg)
    poly = pv.wrap(channel)
    face_centers = channel.triangles_center
    labels = np.ones(len(channel.faces), dtype=np.int32)
    ztop = np.max(abs(channel.triangles[:, :, 2]-d['top']), axis=1) < 1e-5
    labels[ztop & (face_centers[:, 0] < 0)] = 2
    labels[ztop & (face_centers[:, 0] > 0)] = 3
    labels[np.max(abs(channel.triangles[:, :, 0]-d['horizontal_end']), axis=1) < 1e-5] = 4
    poly.cell_data['ModelFaceID'] = labels
    poly.field_data['coordinate_scale_to_m'] = [.001]
    poly.field_data['ModelFaceNames'] = ['WALL', 'INLET', 'OUTLET1', 'OUTLET2']
    poly.save(out/'labelled_straight_lumen_mm.vtp')
    if not all(gates.values()): raise ValueError('GEOMETRY_GATES_FAILED: '+str(gates))
    return parts, summary


def render(parts, summary, out):
    if summary.get('fixture_mode') == 'bottomless_external_dovetail':
        from .bottomless_channel_fixture import render as fixture_render
        return fixture_render(parts, summary, out)
    (out/'QC').mkdir(exist_ok=True)
    d = summary['design']
    internal = d.get('chamber_mode') == 'closed_wall_internal_end'
    for name in ['assembly', 'top', 'side', 'core', 'exploded', 'internal_end' if internal else 'horizontal_gate']:
        plot = pv.Plotter(off_screen=True, window_size=(1600, 1100))
        plot.set_background('white')
        for key, part in parts.items():
            if name == 'core' and 'core' not in key: continue
            m = part.copy()
            if name == 'exploded':
                if 'gate' in key or 'retainer' in key: m.apply_translation([0, 0, 10])
                elif 'core' in key: m.apply_translation([0, 0, 4])
            color = '#c45142' if 'core' in key else '#84a9bb' if 'box' in key else '#d4a64c'
            plot.add_mesh(pv.wrap(m), color=color, opacity=.20 if 'box' in key else 1., smooth_shading=True)
        if name not in ['core', 'exploded']:
            inner = np.asarray(d['inner'])
            plane = pv.Plane(center=[*inner.mean(0), d['fill']], direction=(0, 0, 1),
                             i_size=np.ptp(inner[:, 0]), j_size=np.ptp(inner[:, 1]))
            plot.add_mesh(plane, color='#64c6d1', opacity=.12)
        if name in ['side', 'core']:
            plot.add_point_labels(np.array([p['center_mm'] for p in summary['ports']]),
                ['INLET (+Z)', 'OUTLET 1 (+Z)', 'OUTLET 2 (internal)' if internal else 'OUTLET 2 (+X)'], font_size=18, text_color='#273541',
                point_size=0, shape_opacity=0, always_visible=True)
        plot.add_text(f"Straight channel | diameter {summary['diameter_mm']:.2f} mm | P1S / 0.4 mm / ABS\n"
                      f"PDMS height {summary['pdms_fill_height_mm']:.2f} mm | floor {summary['pdms_minimum_floor_mm']:.2f} mm | wall / base {summary['wall_mm']:.2f} mm",
                      font_size=14, color='#263238')
        plot.add_axes()
        plot.show_bounds(xtitle='X (mm)', ytitle='Y (mm)', ztitle='Z (mm)', fmt='%.2f', font_size=12, color='#33434b', location='outer')
        if name == 'top': plot.view_xy()
        elif name == 'side': plot.view_xz()
        else: plot.view_isometric()
        plot.camera.zoom(.85)
        if name in ['horizontal_gate', 'internal_end']:
            focus = [d['outer'][1][0], 0, d['center_z']+2]
            plot.camera_position = ([focus[0]+15, -18, focus[2]+13], focus, [0, 0, 1])
        plot.screenshot(out/'QC'/(name+'.png'))
        plot.close()


def arrange(parts, out):
    if any('bottomless_box' in n for n in parts):
        from .bottomless_channel_fixture import arrange as fixture_arrange
        return fixture_arrange(parts, out)
    meshes, records = [], []
    for (name, part), xy in zip(parts.items(), [(80., 85.), (145., 85.), (75., 125.), (105., 125.), (140., 125.)]):
        transform = np.eye(4)
        if 'retainer' in name: transform[:3, :3] = Rotation.from_euler('x', 180, degrees=True).as_matrix()
        elif 'gate' in name: transform[:3, :3] = Rotation.from_euler('y', 90, degrees=True).as_matrix()
        m = part.copy()
        m.apply_transform(transform)
        delta = np.r_[xy, 0]-np.r_[m.bounds.mean(0)[:2], m.bounds[0, 2]]
        m.apply_translation(delta)
        transform[:3, 3] = delta
        records.append(dict(name=name, support_enabled='core' in name, plate_bounds_mm=m.bounds,
            casting_to_plate_transform=transform, plate_to_casting_transform=np.linalg.inv(transform)))
        meshes.append(m)
    envelopes = [m.bounds[:, :2]+[[-7, -7], [7, 7]] for m in meshes]
    for i, a in enumerate(envelopes):
        if not ((a[0] >= [20, 30]).all() and (a[1] <= 250).all()): raise ValueError('OUTSIDE_P1S_BED')
        for b in envelopes[i+1:]:
            if not ((a[1] < b[0]).any() or (b[1] < a[0]).any()): raise ValueError('PLATE_OVERLAP')
    save(out/'print_layout.json', dict(passed=True, reserved_margin_mm=7, records=records))
    (out/'print_parts').mkdir(exist_ok=True)
    for m, r in zip(meshes, records): m.export(out/'print_parts'/(r['name']+'.stl'))
    return meshes, records


def write_project(path, meshes, records, profiles, cfg):
    shared.write_project(path, meshes, records, profiles, cfg)
    with zipfile.ZipFile(path) as z: data = {n: z.read(n) for n in z.namelist()}
    settings = ET.fromstring(data['Metadata/model_settings.config'])
    for row in settings.findall('.//plate/metadata'):
        if row.get('key') == 'plater_name':
            row.set('value', 'Straight D1.6 - P1S 0.4 ABS - '+('two top ports only' if two_top_ports(cfg) else 'two up, one horizontal'))
    data['Metadata/model_settings.config'] = ET.tostring(settings, encoding='utf-8', xml_declaration=True)
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        for name, value in data.items(): z.writestr(name, value)


def slice_plate(cfg, parts, summary, out):
    if cfg.get('delivery_mode')=='manual_slicing':
        raise ValueError('AUTOMATIC_SLICING_DISABLED')
    from .support_slicer_adapter import slice_local, production_geometry
    from .bambu_casting_mold import parse_paths, sample_segments
    from .all_abs_print_plate import read_objects
    meshes, records = arrange(parts, out)
    actual, engine = shared.get_profiles(cfg, out/'profiles')
    prefix = 'Straight_2Top_D1p6_P1S_04_ABS' if two_top_ports(cfg) else PREFIX
    project = out/(prefix+'_all_parts_project.3mf')
    write_project(project, meshes, records, actual, cfg)
    trimesh.util.concatenate(meshes).export(out/'plate_reference.stl')
    result = slice_local(engine, actual, project, out/'bambu', timeout=900)
    if not result['slicer_success']: raise RuntimeError('NATIVE_SLICE_FAILED: '+str(result))
    archive = Path(result['archive'])
    geometry_check = production_geometry(archive, out/'plate_reference.stl')
    paths = parse_paths(archive)
    support = paths['segments'][paths['kind'] == 1]
    model = paths['segments'][paths['kind'] == 0]
    objects, settings, nplates, crc = read_objects(archive)
    lo, hi = meshes[0].bounds[:, :2]
    mid = model.mean(1)
    core_model = model[((mid[:, :2] >= lo-.5) & (mid[:, :2] <= hi+.5)).all(1)]
    if not len(core_model): raise ValueError('CORE_DROPPED_BY_SLICER')
    first_z = float(core_model[:, :, 2].min())
    raft = first_z-float(actual['process']['initial_layer_print_height'])
    if not 0 <= raft <= 1.: raise ValueError('UNEXPECTED_RAFT_LIFT')
    distances = cKDTree(sample_segments(core_model, .15)).query(meshes[0].triangles_center+[0, 0, raft])[0]
    # Check every port tip and the whole horizontal run independently too.
    samples = []
    d = summary['design']
    if two_top_ports(cfg):
        from .straight_two_top_ports import centerline_samples
        samples.extend(centerline_samples(d))
    else:
        for x in np.linspace(d['start_x'], d['horizontal_end']-.2, 80):
            samples.append([x, 0, d['center_z']])
    for x in [d['inlet_x'], d['junction_x']]:
        for z in np.linspace(d['center_z']+d['bend_radius'], d['top']-.2, 30): samples.append([x, 0, z])
    T = np.array(records[0]['casting_to_plate_transform'])
    q = np.array(samples)@T[:3, :3].T+T[:3, 3]+[0, 0, raft]
    center_distances = cKDTree(sample_segments(core_model, .15)).query(q)[0]
    object_settings = [o['metadata'] for o in objects]
    scale_qc = shared.slicer_scale_qc(cfg, settings)
    gates = dict(valid_archive=crc and nplates == 1,
        no_duplicate_dimensional_compensation=scale_qc['passed'],
        geometry_preserved=geometry_check['passed'],
        P1S_04_ABS=settings.get('printer_model') == 'Bambu Lab P1S' and settings.get('nozzle_diameter') == ['0.4'] and settings.get('filament_type') == ['ABS'],
        all_production_objects=len(objects) == len(parts) and {o['metadata'].get('name') for o in objects} == set(parts),
        native_support_present=len(support) > 0,
        whole_core_present=float(distances.max()) < .5,
        all_channel_paths_present=float(center_distances.max()) < d['radius'],
        no_critical_slicer_warnings=not result['warning_audit']['critical'])
    if not two_top_ports(cfg):
        gates['all_three_paths_present'] = gates['all_channel_paths_present']
    audit = dict(gates=gates, actual_settings=result['actual_project_settings'], geometry_check=geometry_check,
        cad_allowance_slicer_scale_qc=scale_qc,
        object_settings=object_settings, expected_object_count=len(parts), actual_object_count=len(objects), support_extrusion_segments=len(support),
        support_filament_mm=float(paths['extrusion_filament_mm'][paths['kind'] == 1].sum()),
        raft_lift_mm=raft, first_core_layer_top_mm=first_z,
        core_surface_max_distance_to_actual_path_mm=float(distances.max()),
        centerline_max_distance_to_actual_path_mm=float(center_distances.max()),
        physical_support_removal_validated=False)
    if shared.bottomless(cfg):
        from .bottomless_channel_fixture import toolpath_qc
        audit['fixture_toolpath_qc'] = toolpath_qc(objects, model, raft)
        gates['all_fixture_toolpaths_present'] = audit['fixture_toolpath_qc']['passed']
    if cfg.get('magnet_slots'):
        from .channel_magnet_slots import native_toolpath_qc
        audit['magnet_toolpath_qc']=native_toolpath_qc(summary,records,model)
        gates['magnet_pockets_preserved']=audit['magnet_toolpath_qc']['passed']
    save(out/'print_qc.json', audit)
    sp = pv.PolyData(support.reshape(-1, 3))
    sp.lines = np.c_[np.full(len(support), 2), np.arange(len(support)*2).reshape(-1, 2)].ravel()
    sp.save(out/'QC/native_print_support_toolpaths_mm.vtp')
    plot = pv.Plotter(off_screen=True, window_size=(1600, 1100))
    plot.set_background('white')
    for obj in objects:
        m = obj['mesh'].copy()
        if 'core' in obj['metadata'].get('name', ''): m.apply_translation([0, 0, raft])
        plot.add_mesh(pv.wrap(m), color='#c45142' if 'core' in obj['metadata'].get('name', '') else '#91acb6', opacity=.75)
    plot.add_mesh(sp, color='#2f7155', line_width=1)
    plot.add_text('Straight D1.6'+(' | two TOP ports only' if two_top_ports(cfg) else '')+' | native P1S / 0.4 ABS slice\nGreen = removable print support', font_size=14, color='black')
    plot.view_isometric()
    plot.add_axes()
    plot.screenshot(out/'QC/native_support_plate.png')
    plot.close()
    sliced = out/(prefix+'_all_parts_sliced.3mf')
    shutil.copy2(archive, sliced)
    final = dict(status='READY_FOR_FIT_AND_PRINT_TRIAL' if all(gates.values()) else 'NEEDS_ADJUSTMENT',
        geometry_status=summary['status'], slice=result, print_qc=audit,
        sidewall_external_seal_required=summary['side_gate_external_seal_required'], physical_casting_validated=False,
        no_printer_upload_or_start=True, files={str(p.relative_to(out)): digest(p)
            for p in [project, sliced, *sorted((out/'assembly').glob('*.stl'))]})
    save(out/'manufacturing_summary.json', final)
    print('STRAIGHT_RESULT', final['status'], result['estimated_print_time_seconds'], result['filament_used_g'], flush=True)
    if not all(gates.values()): raise ValueError('NATIVE_QC_FAILED: '+str(gates))
    return final


def package_delivery(out):
    final = json.loads((out/'manufacturing_summary.json').read_text())
    if final['status'] != 'READY_FOR_FIT_AND_PRINT_TRIAL': raise ValueError('CANNOT_PACKAGE_FAILED_RUN')
    cfg = json.loads((out/'effective_config.json').read_text())
    files = [*sorted(out.glob('*.3mf')), *sorted((out/'QC').glob('*.png')),
             *[out/'assembly'/(n+'.stl') for n in part_names(cfg)],
             out/'assembly/Straight_lumen_core.step', out/'labelled_straight_lumen_mm.vtp',
             *[out/n for n in ['geometry_qc.json', 'print_qc.json', 'manufacturing_summary.json', 'effective_config.json', 'print_layout.json']],
             *[out/'bambu'/n for n in ['command.json', 'stdout.log', 'stderr.log', 'slice_result.json']]]
    if cfg.get('rhino_export'):
        from .rhino_channel_export import verify_exports
        verify_exports(out)
        files += sorted(out.glob('*.3dm')) + [out/'rhino_qc.json']
    members = {str(p.relative_to(out)): p for p in files}
    if cfg.get('magnet_slots'):
        members.update({str(p.relative_to(out)):p for p in [out/'magnet_slots_qc.json',*sorted((out/'reference').iterdir())]})
        members['磁铁槽使用说明.md']=ROOT/'docs/三版模具_外置磁铁槽使用说明.md'
    members['使用说明.md'] = ROOT/'docs'/('直管1.6mm_仅双竖直端口打印说明.md' if two_top_ports(cfg) else '直管1.6mm_双竖直单水平端口打印说明.md')
    hashes = {k: digest(v) for k, v in members.items()}
    dest = ROOT/'deliverables'/('Straight_'+out.name+'.zip')
    dest.parent.mkdir(exist_ok=True)
    with zipfile.ZipFile(dest, 'w', zipfile.ZIP_DEFLATED) as z:
        for name, path in members.items(): z.write(path, 'Straight_channel/'+name)
        z.writestr('Straight_channel/SHA256.json', json.dumps(hashes, indent=2, ensure_ascii=False))
    with zipfile.ZipFile(dest) as z:
        if z.testzip() is not None: raise ValueError('PACKAGE_CRC_FAILED')
        for k, h in hashes.items():
            if hashlib.sha256(z.read('Straight_channel/'+k)).hexdigest() != h: raise ValueError('PACKAGE_MEMBER_HASH_FAILED')
    save(out/'delivery_receipt.json', dict(zip=dest, zip_sha256=digest(dest), zip_bytes=dest.stat().st_size,
        crc_and_hashes_passed=True, members_sha256=hashes))


def run(config_path, stage='all'):
    if stage not in ['all','geometry','project']:
        raise ValueError('AUTOMATIC_SLICING_DISABLED: use geometry, project or all; slice manually in Bambu Studio')
    cfg = load_config(config_path)
    out = cfg['output']
    out.mkdir(parents=True, exist_ok=True)
    print('STRAIGHT_OUTPUT', out, flush=True)
    if stage in ['all', 'geometry']:
        parts, summary = geometry(cfg, out)
        render(parts, summary, out)
        arrange(parts, out)
        if cfg.get('rhino_export'):
            from .rhino_channel_export import export_models
            export_models(cfg, parts, summary, out)
        print('STRAIGHT_GEOMETRY', summary['status'], summary['mold_outer_size_mm'], flush=True)
    else:
        effective = json.loads((out/'effective_config.json').read_text())
        if any(v != effective.get(k) for k, v in cfg.items() if k != 'output'):
            raise ValueError('CONFIG_CHANGED_REGENERATE_GEOMETRY')
        summary = json.loads((out/'geometry_qc.json').read_text())
        if summary['status'] != 'GEOMETRY_PASSED': raise ValueError('GEOMETRY_PASS_REQUIRED')
        for name, expected in summary['assembly_sha256'].items():
            if digest(out/'assembly'/name) != expected: raise ValueError('ASSEMBLY_CHECKSUM_MISMATCH')
        parts = {n: trimesh.load_mesh(out/'assembly'/(n+'.stl')) for n in part_names(cfg)}
    if stage in ['all', 'project']:
        from .channel_manual_delivery import export_delivery
        return export_delivery(cfg,parts,summary,out,'Straight')
    return summary
