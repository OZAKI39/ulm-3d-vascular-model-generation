"""Diameter-controlled manufacturing variant of the immutable labelled Y input.

Source branch axes, terminal positions and lengths are retained, but its wall
is deliberately reconstructed: the old 2.6/2.0 mm mesh is NOT a 1.6 mm lumen.
Circular CAD sweeps and local junction fillets replace that wall. No global
scaling or smoothing is used. All output coordinates are millimetres.
"""
from __future__ import annotations

import cadquery as cq
import numpy as np
import pyvista as pv
import trimesh
from shapely.geometry import Point
from . import y_channel_casting as shared
from .straight_channel_casting import cad_mesh, volume_intersection


def local_section_diameter(mesh, center, normal):
    """Measure the connected section containing the sampled branch centre.

    The same plane can also cut a distant branch or riser. Selecting by the
    known centre avoids adding unrelated section areas to a lumen diameter.
    """
    section = mesh.section(plane_origin=center, plane_normal=normal)
    if section is None:
        raise ValueError('MISSING_BRANCH_SECTION')
    planar, to_world = section.to_2D()
    point = (np.linalg.inv(to_world) @ np.r_[center, 1])[:2]
    polygons = [p for p in planar.polygons_full if p.covers(Point(point))]
    if len(polygons) != 1:
        raise ValueError('AMBIGUOUS_BRANCH_SECTION')
    return float(2*np.sqrt(polygons[0].area/np.pi))


def build(cfg, planes):
    if cfg['diameter_mm'] != 1.6:
        raise ValueError('THIS_VARIANT_REQUIRES_1P6_MM')
    r = cfg['diameter_mm']/2
    R = r*cfg['bend_radius_to_lumen_radius']
    z = cfg['base_mm']+cfg['pdms_floor_mm']+r
    fill = max(z+R+cfg['pdms_above_elbow_mm'],
               cfg['base_mm']+shared.minimum_cad(cfg, 'height'))
    rim = fill+cfg['freeboard_mm']
    collar_z = shared.retainer_arm_bottom(cfg, rim)+cfg['retainer_thickness_mm']
    top = collar_z+cfg['port_above_retainer_mm']
    junction = np.asarray(planes['BIFURCATION']['origin_m'])*1000+[0, 0, z]
    solids, ports, samples = [], [], []
    v = lambda xyz: cq.Vector(*xyz)
    for name, label in shared.PORTS:
        n = np.asarray(planes[name]['normal'], float)
        n /= np.linalg.norm(n)
        origin = np.asarray(planes[name]['origin_m'])*1000+[0, 0, z]
        if abs(n[2]) > 1e-9 or not np.allclose(np.cross(origin-junction, n), 0, atol=1e-8):
            raise ValueError('SOURCE_MUST_BE_STRAIGHT_PLANAR_Y_BRANCHES')
        start = origin+cfg['terminal_straight_mm']*n
        knee = start+R*n+[0, 0, R]
        midpoint = start+R/np.sqrt(2)*n+[0, 0, R*(1-1/np.sqrt(2))]
        cap = np.r_[knee[:2], top]
        wire = cq.Wire.assembleEdges([cq.Edge.makeLine(v(junction), v(start)),
            cq.Edge.makeThreePointArc(v(start), v(midpoint), v(knee)),
            cq.Edge.makeLine(v(knee), v(cap))])
        solids.append(cq.Solid.sweep(cq.Wire.makeCircle(r, v(junction), v(n)), [], wire))
        length = float(np.linalg.norm(origin-junction))
        # Exclude the non-circular merging junction itself. All branches,
        # bends and risers are checked from their actual tessellated surface.
        for s in np.linspace(2.5, length-.2, 6):
            samples.append(dict(port=name, region='body', center_mm=junction+s*n, normal=n))
        for angle in np.linspace(.12, np.pi/2-.12, 5):
            samples.append(dict(port=name, region='elbow',
                center_mm=start+R*np.sin(angle)*n+[0, 0, R*(1-np.cos(angle))],
                normal=np.cos(angle)*n+[0, 0, np.sin(angle)]))
        for h in [knee[2]+.5, fill-.5, top-.4]:
            samples.append(dict(port=name, region='riser', center_mm=np.r_[knee[:2], h], normal=[0, 0, 1]))
        ports.append(dict(name=name, label_id=label, radius_mm=r, nominal_diameter_mm=2*r,
            original_origin_mm=origin, original_normal=n, source_branch_length_mm=length,
            source_original_diameter_mm=planes[name]['radius_m']*2000,
            bend_radius_mm=R, bend_start=start, top_center_mm=cap, collar_bottom_mm=collar_z,
            vertical_axis=[0, 0, 1], analytical_centerline_tangent_jump_deg=0.))
    channel = solids[0].fuse(*solids[1:]).clean()
    # First the two daughter-junction edges, then the four parent edges.
    # The separate fillets avoid OCP's singular simultaneous three-way fillet.
    edges = [e for e in channel.Edges() if e.geomType() == 'ELLIPSE' and e.Length() > 1.5]
    if len(edges) != 2: raise ValueError('EXPECTED_TWO_DAUGHTER_JUNCTION_EDGES')
    channel = channel.fillet(cfg['junction_fillet_mm'], edges)
    edges = [e for e in channel.Edges() if e.geomType() == 'ELLIPSE']
    if len(edges) != 4: raise ValueError('EXPECTED_FOUR_PARENT_JUNCTION_EDGES')
    channel = channel.fillet(cfg['parent_junction_fillet_mm'], edges)
    core = channel
    for p in ports:
        core = core.fuse(cq.Solid.makeCylinder(r+cfg['collar_radial_width_mm'], cfg['collar_thickness_mm'],
            v([*p['top_center_mm'][:2], collar_z]), v([0, 0, 1]))).clean()
    heights = dict(source_z_translation_mm=z, pdms_floor_z_mm=cfg['base_mm'], fill_z_mm=fill,
                   rim_z_mm=rim, port_top_z_mm=top)
    return core, channel, ports, heights, samples


def geometry(cfg, out):
    source, _, planes = shared.import_source(cfg)
    core_cad, channel_cad, ports, heights, samples = build(cfg, planes)
    core, channel = cad_mesh(core_cad, cfg), cad_mesh(channel_cad, cfg)
    inner = channel.bounds[:, :2]+[[-cfg['side_clearance_mm']]*2, [cfg['side_clearance_mm']]*2]
    inner = shared.enforce_pdms_width(inner, cfg)
    outer = inner+[[-cfg['wall_mm']]*2, [cfg['wall_mm']]*2]
    mold = trimesh.boolean.difference([
        shared.box([np.r_[outer[0], 0], np.r_[outer[1], heights['rim_z_mm']]]),
        shared.box([np.r_[inner[0], cfg['base_mm']], np.r_[inner[1], heights['rim_z_mm']+1]])], engine='manifold')
    parts = {'Y_lumen_core': core, 'Y_casting_box_base_2mm': mold}
    for p in ports:
        parts['Y_retainer_'+p['name']] = shared.retainer(p, inner, cfg, heights)
    assembly = out/'assembly'
    assembly.mkdir(parents=True, exist_ok=True)
    for name, mesh in parts.items(): mesh.export(assembly/(name+'.stl'))
    channel.export(assembly/'Y_channel_without_positioning_collars_reference.stl')
    source.apply_translation([0, 0, heights['source_z_translation_mm']])
    source.export(assembly/'Y_original_body_mm_reference.stl')
    cq.exporters.export(core_cad, str(assembly/'Y_lumen_core.step'))
    for row in samples:
        row['equivalent_diameter_mm'] = local_section_diameter(channel, row['center_mm'], row['normal'])
    cap = np.max(abs(channel.triangles[:, :, 2]-heights['port_top_z_mm']), axis=1) < 1e-5
    labels = np.ones(len(channel.faces), dtype=np.int32)
    for p in ports:
        selected = cap & (np.linalg.norm(channel.triangles_center[:, :2]-p['top_center_mm'][:2], axis=1) < p['radius_mm']+.01)
        p['measured_cap_diameter_mm'] = float(2*np.sqrt(channel.area_faces[selected].sum()/np.pi))
        p['minimum_cap_normal_dot'] = float(channel.face_normals[selected, 2].min())
        labels[selected] = p['label_id']
    poly = pv.wrap(channel)
    poly.cell_data['ModelFaceID'] = labels
    poly.field_data['coordinate_scale_to_m'] = [.001]
    poly.field_data['ModelFaceNames'] = ['WALL', 'INLET', 'OUTLET_A', 'OUTLET_B']
    poly.save(out/'labelled_y_D1p6_lumen_mm.vtp')
    keys = list(parts)
    collision = {a+' / '+b: volume_intersection(parts[a], parts[b])
                 for i, a in enumerate(keys) for b in keys[i+1:]}
    fill = shared.box([np.r_[inner[0], cfg['base_mm']], np.r_[inner[1], heights['fill_z_mm']]])
    mesh_qc = {k: shared.mesh_info(m) for k, m in parts.items()}
    floor = float(core.bounds[0, 2]-cfg['base_mm'])
    retainer_clearance = min(m.bounds[0, 2] for k, m in parts.items() if 'retainer' in k)-heights['fill_z_mm']
    width_qc = shared.pdms_width_qc(mold, inner, cfg, heights['fill_z_mm'])
    height_qc = shared.pdms_height_qc(mold, inner, cfg, heights['fill_z_mm'])
    gates = dict(single_closed_parts=all(m.is_volume and len(m.split()) == 1 for m in parts.values()),
        pdms_width_minimum=width_qc['passed'],
        source_file_unchanged=shared.digest(cfg['source']) == cfg['source_sha256'],
        diameter_1p6_all_sampled_sections=all(abs(s['equivalent_diameter_mm']-1.6) < .003 for s in samples),
        diameter_1p6_all_caps=all(abs(p['measured_cap_diameter_mm']-1.6) < .003 for p in ports),
        requested_floor=abs(floor-cfg['pdms_floor_mm']) < .003,
        pdms_height_minimum=height_qc['passed'],
        three_upward_ports=all(p['minimum_cap_normal_dot'] > .99999 for p in ports),
        no_assembly_penetration=max(collision.values()) < 1e-4, retainers_outside_pdms=retainer_clearance > .3)
    summary = dict(status='GEOMETRY_PASSED' if all(gates.values()) else 'GEOMETRY_FAILED', units='mm',
        geometry_mode='uniform_diameter', source_sha256=shared.digest(cfg['source']), planes_sha256=shared.digest(cfg['planes']),
        source_to_mm_scale=1000., source_wall_preserved=False, global_geometric_scaling=False,
        construction='Circular CAD sweeps on original branch axes; local 0.10/0.05 mm junction fillets; uniform diameter manufacturing reconstruction',
        original_body_translation_mm=[0, 0, heights['source_z_translation_mm']],
        diameter_mm=cfg['diameter_mm'], heights=heights, pdms_minimum_floor_mm=floor,
        pdms_fill_height_mm=heights['fill_z_mm']-cfg['base_mm'],
        pdms_fill_volume_ml=(fill.volume-volume_intersection(core, fill))/1000,
        wall_mm=cfg['wall_mm'], base_mm=cfg['base_mm'], cavity_xy_bounds_mm=inner,
        cavity_xy_size_mm=inner[1]-inner[0], mold_outer_size_mm=mold.extents, pdms_width_qc=width_qc, pdms_height_qc=height_qc,
        ports=ports, section_checks=samples, meshes=mesh_qc,
        assembly_pair_intersection_volume_mm3=collision, retainer_minimum_above_fill_mm=retainer_clearance,
        physical_fit_validated=False, optical_quality_validated=False, gates=gates,
        assembly_sha256={p.name: shared.digest(p) for p in sorted(assembly.iterdir())})
    shared.save(out/'geometry_qc.json', summary)
    shared.save(out/'effective_config.json', cfg)
    if not all(gates.values()): raise ValueError('GEOMETRY_GATES_FAILED: '+str(gates))
    return parts, summary
