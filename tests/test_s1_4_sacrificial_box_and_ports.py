"""Retained geometry tests; the tiny fixture contains no BG001 geometry."""
import copy
import json
from pathlib import Path

import cadquery as cq
import numpy as np
import pytest
import trimesh

from vascular_processing import sacrificial_fixture as f
from vascular_processing.boundary_review_export import protect_inputs


@pytest.fixture
def cfg():
    return f.load_config(f.ROOT / 'config/sacrificial_box_BG001.yaml')


@pytest.fixture
def synthetic(cfg):
    path = f.ROOT / 'tests/data/simple_port_tree.swc'
    graph = f.read_source(path).graph
    manifest = json.loads(path.with_suffix('.json').read_text())
    mapping = [dict(swc_id=n, original_swc_id=n, original_radius_mm=.6, manufacturing_radius_mm=.6) for n in graph]
    edges = [dict(parent_id=a, child_id=b, original_branch=name) for name, ids in manifest['branches'].items() for a, b in zip(ids, ids[1:])]
    endpoints, branches = f.decode_endpoints(graph, mapping, edges, np.eye(4), cfg)
    return graph, endpoints, branches, manifest


def endpoint(position=(0, 0, 10), tangent=(1, 0, 0), name='I1', role='inlet'):
    return f.Endpoint(name, role, 1, 1, 'a', np.array(position, float), .6,
                      np.array(tangent, float), 5., 4, .6, .6)


def roomy_box():
    return f.BoxBounds(np.array([[-30, -30, -10], [30, 30, 40]], float),
                       np.array([[-33, -33, -13], [33, 33, 40]], float))


def straight(cfg, position=(0, 0, 10), name='I1', role='inlet'):
    e = endpoint(position, name=name, role=role)
    return f.make_route(e, '+X', np.array([30, position[1], position[2]]), roomy_box(), cfg, 'ROUTE_A_STRAIGHT')


def test_companion_endpoints_and_manifest(synthetic):
    graph, endpoints, _, m = synthetic
    assert [e.swc_id for e in endpoints if e.role == 'inlet'] == m['inlet_ids']
    assert [e.swc_id for e in endpoints if e.role == 'outlet'] == m['outlet_ids']
    assert [n for n in graph if graph.out_degree(n) > 1] == m['bifurcation_ids']
    assert {e.branch_id for e in endpoints} == {'inlet', 'upper', 'lower'}


def test_tangent_stable_and_outward(synthetic, cfg):
    graph, endpoints, _, _ = synthetic
    inlet = endpoints[0]
    np.testing.assert_allclose(inlet.tangent, [-1, 0, 0], atol=1e-9)
    assert inlet.tangent_window_mm == 5 and inlet.sample_count >= 3
    perturbed = copy.deepcopy(graph)
    perturbed.nodes[2]['coords'][1] += .01
    t, *_ = f.endpoint_tangent(perturbed, 1, cfg)
    assert f.angle(t, inlet.tangent) < .1


def test_short_tangent_window_recorded(synthetic, cfg):
    graph, _, _, _ = synthetic
    small = graph.subgraph([1, 2]).copy()
    _, window, count, _ = f.endpoint_tangent(small, 1, cfg)
    assert window == 3 and count == 2


def test_box_margins(cfg):
    b = f.compute_box(np.array([[0, 0, 0], [10, 20, 30]]), cfg)
    np.testing.assert_allclose(b.inner_size, [34, 44, 52])
    np.testing.assert_allclose(b.outer_size, [40, 50, 55])
    np.testing.assert_allclose(b.inner[0], [-12, -12, -10])


def test_open_top_box_dimensions_and_floor(cfg):
    b = roomy_box()
    shape, holes = f.create_box(b, [], cfg)
    assert holes == []
    assert shape.isInside(cq.Vector(0, 0, -11))
    assert not shape.isInside(cq.Vector(0, 0, 39))
    assert shape.isInside(cq.Vector(32, 0, 39))
    assert f.mesh_qc(f.shape_mesh(shape, cfg), cfg)['passed']
    assert shape.BoundingBox().zlen == pytest.approx(53, abs=1e-6)
    measurements = f.box_dimension_checks(shape, b, cfg)
    assert measurements['passed']
    assert measurements['bottom_thickness_mm'] == pytest.approx(3)


def test_straight_route_tangent_and_outside_stub(cfg):
    r = straight(cfg)
    np.testing.assert_allclose(r.target, [30, 0, 10])
    np.testing.assert_allclose(r.end, [43, 0, 10])
    assert r.length == 43 and r.turn == 0
    assert r.mesh.is_volume


def test_straight_route_rejects_direction_kink(cfg):
    with pytest.raises(ValueError, match='STRAIGHT_MUST_FOLLOW'):
        f.make_route(endpoint(), '+X', np.array([30, 5, 10]), roomy_box(), cfg, 'ROUTE_A_STRAIGHT')


def test_curved_route_spline_sweep_and_bend_radius(cfg):
    r = f.make_route(endpoint(tangent=(0, 1, 0)), '+X', np.array([30, 15, 10]), roomy_box(), cfg, 'ROUTE_B_CURVED', 1.5)
    assert r.minimum_bend_radius >= 6
    assert r.shape.isValid() and r.mesh.is_volume
    assert r.turn == 90
    np.testing.assert_allclose(r.end, [43, 15, 10])


def test_wall_edge_clearance():
    assert f.wall_edge_clearance(np.array([30, 25, 10]), '+X', roomy_box()) == 5
    assert f.wall_edge_clearance(np.array([30, 0, 38]), '+X', roomy_box()) == 2


def test_natural_spline_turn_preserves_curvature(cfg):
    e = endpoint(tangent=(0, 1, 0))
    points = f.natural_bend_points(e, '+X', 8, cfg)
    target = points[-1].copy()
    target[0] = 30
    route = f.make_route(e, '+X', target, roomy_box(), cfg, 'ROUTE_B_CURVED', bend_radius=8)
    assert route.minimum_bend_radius > 7.9
    assert route.minimum_bend_radius < 8.1


def test_spacing_uses_stricter_constraint(cfg):
    a, b = endpoint(), endpoint()
    assert f.required_spacing(a, b, cfg) == 8
    a.radius = b.radius = 5
    assert f.required_spacing(a, b, cfg) == 10.4


def test_surface_collision_rejection(cfg):
    route = straight(cfg)
    obstacle = trimesh.creation.box([2, 2, 2], transform=trimesh.transformations.translation_matrix([15, 0, 10]))
    assert f.surface_distance(obstacle, route.mesh) == 0
    obstacle.apply_translation([0, 5, 0])
    assert f.surface_distance(obstacle, route.mesh) > 3


def test_parent_collar_exemption_does_not_hide_other_branches(cfg, synthetic):
    graph, endpoints, branches, _ = synthetic
    inlet = endpoints[0]
    # Two small surfaces: one at parent attachment, one on a sibling endpoint.
    spheres = [trimesh.creation.icosphere(subdivisions=1, radius=.6).apply_translation(e.position) for e in (inlet, endpoints[1])]
    mesh = trimesh.util.concatenate(spheres)
    checker = f.VesselCollision(mesh, graph, branches, np.eye(4), endpoints, cfg)
    assert 0 < checker.exempt_face_counts['I1'] < len(mesh.faces)


def test_face_assignment_deterministic_and_spacing_rejected(cfg):
    a, b, c = straight(cfg), straight(cfg, (0, 12, 10), 'O1', 'outlet'), straight(cfg, (0, -12, 10), 'O2', 'outlet')
    for index, r in enumerate((a, b, c)):
        r.option_id, r.score = str(index), float(index)
    first, _ = f.select_layout([[a], [b], [c]], cfg)
    second, _ = f.select_layout([[a], [b], [c]], cfg)
    assert [r.option_id for r in first] == [r.option_id for r in second] == ['0', '1', '2']
    b.target = a.target.copy()
    chosen, audit = f.select_layout([[a], [b]], cfg)
    assert chosen == [] and audit['status'] == 'PORT_LAYOUT_NEEDS_MANUAL_REVIEW'


def test_port_union_single_connected_and_source_volume_preserved(cfg, tmp_path):
    route = straight(cfg)
    core = f.shape_mesh(cq.Solid.makeCylinder(.6, 5, pnt=cq.Vector(-5, 0, 10), dir=cq.Vector(1, 0, 0)), cfg)
    source = tmp_path / 'source.stl'
    core.export(source)
    before = f.sha256(source)
    result, qc = f.union_core(core, [route], cfg)
    assert qc['connected_components'] == 1 and qc['watertight'] and qc['source_solid_preserved']
    assert result.volume > core.volume
    assert f.sha256(source) == before


def test_complete_synthetic_tree_with_all_three_ports(cfg, synthetic):
    graph, endpoints, _, manifest = synthetic
    pieces = []
    for ids in manifest['branches'].values():
        p, q = [graph.nodes[n]['coords'][:3] for n in (ids[0], ids[-1])]
        delta = q - p
        pieces.append(cq.Solid.makeCylinder(.6, float(np.linalg.norm(delta)),
            pnt=cq.Vector(*p), dir=cq.Vector(*f.unit(delta))))
    pieces.append(cq.Solid.makeSphere(.6, pnt=cq.Vector(0, 0, 10)))
    tree = f.shape_mesh(pieces[0].fuse(*pieces[1:]), cfg)
    assert f.mesh_qc(tree, cfg)['passed']
    routes = []
    for e in endpoints:
        face = '-X' if e.role == 'inlet' else '+X'
        if e.role == 'inlet':
            target = np.array([-30, 0, 10])
            r = f.make_route(e, face, target, roomy_box(), cfg, 'ROUTE_A_STRAIGHT')
        else:
            target = f.natural_bend_points(e, face, 8, cfg)[-1]
            target[0] = 30
            r = f.make_route(e, face, target, roomy_box(), cfg, 'ROUTE_B_CURVED', bend_radius=8)
        routes.append(r)
    core, qc = f.union_core(tree, routes, cfg)
    assert len(routes) == 3 and qc['passed'] and qc['source_solid_preserved']
    box, holes = f.create_box(roomy_box(), routes, cfg)
    assert len(holes) == 3 and all(h['passed'] for h in holes)
    assert f.surface_distance(core, f.shape_mesh(box, cfg)) > .15


def test_start_cross_section_conflict_rejects_every_wall(cfg):
    class BlockedStart:
        start_checks = {'I1': dict(passed=False, distance_mm=1.8)}

        def distance(self, route):
            raise AssertionError('A known infeasible starting disk must not be bypassed')

    cfg['routing']['target_lateral_offsets_mm'] = [0]
    cfg['routing']['target_vertical_offsets_mm'] = [0]
    options, rows = f.generate_options(endpoint(), roomy_box(), cfg, BlockedStart())
    assert not options and {r['face'] for r in rows} == set(f.FACES)
    assert all(r['status'] == 'ENDPOINT_CLEARANCE_CONFLICT' for r in rows)


def test_assembly_layers_and_two_body_3mf_round_trip(cfg, tmp_path):
    from vascular_processing.sacrificial_fixture_review import assembly_export
    import pyvista as pv
    route = straight(cfg)
    source = f.shape_mesh(cq.Solid.makeCylinder(.6, 5, pnt=cq.Vector(-5, 0, 10), dir=cq.Vector(1, 0, 0)), cfg)
    core, _ = f.union_core(source, [route], cfg)
    box, _ = f.create_box(roomy_box(), [route], cfg)
    result = assembly_export(tmp_path, 'synthetic', source, [route], f.shape_mesh(box, cfg), core)
    blocks = pv.read(result['vtm'])
    assert list(blocks.keys()) == ['vascular_core', 'inlet_ports', 'outlet_ports', 'casting_box']
    assert result['three_mf_status'] == 'TWO_BODIES_ROUND_TRIP_VERIFIED'
    assert result['vascular_layers'] == 'DISJOINT_PATCHES_OF_BOOLEAN_EXTERNAL_SURFACE'
    assert not result['internal_attachment_caps_displayed']
    assert sum(blocks[key].n_cells if blocks[key] is not None else 0
        for key in ['vascular_core', 'inlet_ports', 'outlet_ports']) == len(core.faces)
    # Combining the three colored patches must reconstruct one external shell,
    # without duplicated internal source/extension caps.
    vascular = trimesh.util.concatenate([trimesh.Trimesh(blocks[key].points,
        blocks[key].faces.reshape(-1, 4)[:, 1:], process=False)
        for key in ['vascular_core', 'inlet_ports']])
    vascular.merge_vertices()
    assert vascular.is_watertight and vascular.euler_number == core.euler_number


def test_disconnected_core_union_fails(cfg):
    route = straight(cfg)
    core = trimesh.creation.box([1, 1, 1])
    with pytest.raises(ValueError, match='CORE_PORT_UNION_FAILED'):
        f.union_core(core, [route], cfg)


def test_wall_holes_count_clearance_and_boolean(cfg):
    routes = [straight(cfg), straight(cfg, (0, 12, 10), 'O1', 'outlet')]
    box, checks = f.create_box(roomy_box(), routes, cfg)
    assert len(checks) == 2 and all(c['passed'] for c in checks)
    assert all(c['radius_mm'] == pytest.approx(.8) for c in checks)
    assert not box.isInside(cq.Vector(31, 0, 10))
    assert box.isInside(cq.Vector(31, 2, 10))


def test_cad_dimensions_ignore_cached_mesh_deflection(cfg):
    box = roomy_box()
    box.inner += [128, 128, 0]
    box.outer += [128, 128, 0]
    e = endpoint(position=(128, 128, 10))
    route = f.make_route(e, '+X', np.array([158, 128, 10]), box, cfg, 'ROUTE_A_STRAIGHT')
    shape, _ = f.create_box(box, [route], cfg)
    f.shape_mesh(shape, cfg)  # Populate OCC triangulation before the measurement.
    assert f.box_dimension_checks(shape, box, cfg)['passed']


@pytest.mark.parametrize('relative', [
    'tests/data/simple_port_tree.swc',
    's1-2_swc_roi_generate_human.py',
    's1-3_swc_roi_generate_MeVO.py',
])
def test_protected_source_writes_rejected_and_hash_unchanged(relative):
    path = f.ROOT / relative
    snapshot = {str(path): f.sha256(path)}
    with protect_inputs(snapshot), pytest.raises(PermissionError):
        path.write_text('must not be written')
    assert f.verify_snapshot(snapshot)['all_unchanged']


def test_all_face_candidates_are_recorded(cfg):
    # Keep this synthetic enumeration small but use the real CAD implementation.
    cfg['routing']['target_lateral_offsets_mm'] = [0]
    cfg['routing']['target_vertical_offsets_mm'] = [0]
    cfg['routing']['spline_tangent_scales'] = [1.0]
    _, rows = f.generate_options(endpoint(), roomy_box(), cfg)
    assert {r['face'] for r in rows} == set(f.FACES)
    assert all(isinstance(r['score'], float) and r['status'] for r in rows)


@pytest.mark.parametrize('face', ['-Z', '+Z'])
def test_z_faces_are_prohibited_in_configuration_and_geometry(cfg, tmp_path, face):
    import yaml
    cfg['routing']['allowed_faces'] = [face]
    path = tmp_path / 'config.yaml'
    path.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match='Z ports are prohibited'):
        f.load_config(path)
    with pytest.raises(ValueError):
        f.make_route(endpoint(), face, np.array([0, 0, 40]), roomy_box(), cfg, 'ROUTE_A_STRAIGHT')


def test_original_terminal_cannot_be_extended_with_nearby_sibling(synthetic):
    graph, endpoints, _, _ = synthetic
    e = endpoints[-1]
    result = f.inspect_native_continuation(graph, e, set(graph))
    assert result['status'] == 'ORIGINAL_TERMINAL_NO_CONTINUATION'
    assert result['available_outward_node_ids'] == []
    assert result['available_arc_length_mm'] == 0


def test_original_path_length_follows_source_until_branching(synthetic):
    graph, _, _, _ = synthetic
    e = endpoint(role='outlet')
    # Source ID 1 continues along 2,3,4,5; node 5 has two daughters.
    result = f.inspect_native_continuation(graph, e, {1})
    assert result['status'] == 'ORIGINAL_BIFURCATION_AVAILABLE'
    assert len(result['stopping_outward_neighbors']) == 2
    ids = [1] + result['available_outward_node_ids']
    assert all(graph.has_edge(a, b) for a, b in zip(ids, ids[1:]))
    assert result['available_arc_length_mm'] == pytest.approx(sum(np.linalg.norm(
        graph.nodes[b]['coords'][:3] - graph.nodes[a]['coords'][:3]) for a, b in zip(ids, ids[1:])))


def test_last_resort_never_skips_available_original_vessel(cfg, monkeypatch):
    e = endpoint()
    evidence = dict(hash_matches_frozen_provenance=True, all_retained_directed_edges_verified=True,
        endpoints=[dict(endpoint_id=e.endpoint_id, status='ORIGINAL_CONTINUATION_AVAILABLE')])
    def forbidden(*args, **kwargs):
        raise AssertionError('Must restore available native vessels before artificial fallback')
    monkeypatch.setattr(f, 'generate_options', forbidden)
    routes, audit, rows = f.try_clearance_fallback([e], [[]], roomy_box(), cfg, None, evidence)
    assert routes == [] and not audit['used'] and rows == []
    assert 'must be restored' in audit['reason']


def test_last_resort_requires_verified_provenance(cfg):
    with pytest.raises(ValueError, match='verified original anatomy'):
        f.try_clearance_fallback([endpoint()], [[]], roomy_box(), cfg, None,
            dict(endpoints=[], hash_matches_frozen_provenance=False))


def test_last_resort_uses_highest_feasible_grid_value_and_preserves_other_constraints(cfg):
    cfg['routing']['allowed_faces'] = ['+X']
    cfg['routing']['target_lateral_offsets_mm'] = [0]
    cfg['routing']['target_vertical_offsets_mm'] = [0]
    cfg['routing']['spline_tangent_scales'] = [1.0]
    obstacle = trimesh.creation.box([1, 1, 1], transform=trimesh.transformations.translation_matrix([10, 2.6, 10]))
    class Collision:
        start_checks = {'I1': dict(distance_mm=10), 'O1': dict(distance_mm=10)}
        def distance(self, route):
            return f.surface_distance(obstacle, route.mesh)
    a = endpoint()
    b = endpoint((0, 12, 10), name='O1', role='outlet')
    collision = Collision()
    a_options, _ = f.generate_options(a, roomy_box(), cfg, collision)
    b_options, _ = f.generate_options(b, roomy_box(), cfg, collision)
    assert not a_options and b_options
    evidence = dict(hash_matches_frozen_provenance=True, all_retained_directed_edges_verified=True,
        endpoints=[dict(endpoint_id='I1', status='ORIGINAL_TERMINAL_NO_CONTINUATION')])
    routes, audit, _ = f.try_clearance_fallback([a, b], [a_options, b_options], roomy_box(), cfg, collision, evidence)
    assert len(routes) == 2 and audit['used']
    assert audit['affected_endpoints'] == ['I1']
    assert audit['selected_clearance_mm'] == pytest.approx(1.4)
    assert next(r for r in routes if r.endpoint.endpoint_id == 'O1').required_vessel_clearance == 3
    assert min(r.nearest_other_port for r in routes) >= 3
    assert all(r.edge_clearance >= 8 and r.minimum_bend_radius >= 6 for r in routes)
    assert all(r.face in f.SIDE_FACES for r in routes)


def test_zero_mm_fallback_still_rejects_intersection(cfg):
    cfg['routing']['allowed_faces'] = ['+X']
    cfg['routing']['target_lateral_offsets_mm'] = [0]
    cfg['routing']['target_vertical_offsets_mm'] = [0]
    cfg['routing']['spline_tangent_scales'] = [1.0]
    class Collision:
        start_checks = {'I1': dict(distance_mm=10)}
        def distance(self, route):
            return 0.0
    options, rows = f.generate_options(endpoint(), roomy_box(), cfg, Collision(), clearance_override=0)
    assert not options
    assert any(row['status'] == 'PORT_COLLISION_RISK' for row in rows)


def test_relaxation_does_not_disable_port_pair_clearance(cfg):
    a = straight(cfg)
    # Different walls so center spacing does not mask the crossing failure.
    b = f.make_route(endpoint((15, -10, 10), (0, 1, 0), 'O1', 'outlet'), '+Y',
        np.array([15, 30, 10]), roomy_box(), cfg, 'ROUTE_A_STRAIGHT')
    for index, r in enumerate([a, b]):
        r.required_vessel_clearance = 0
        r.option_id, r.score = str(index), .1
    routes, _ = f.select_layout([[a], [b]], cfg)
    assert routes == []


def known_cap_source():
    return trimesh.creation.cylinder(radius=.6, height=5, sections=24,
        transform=trimesh.transformations.translation_matrix([0, 0, -2.5]))


def align_synthetic_cap(cfg):
    e = endpoint(position=(0, 0, 0), tangent=(np.sin(np.radians(42)), 0, np.cos(np.radians(42))),
        name='O3', role='outlet')
    source = known_cap_source()
    f.align_known_cap(e, source, {'vertices': 24}, cfg)
    return e, source


def test_known_cap_corrects_mean_direction_and_keeps_source_unchanged(cfg):
    source = known_cap_source()
    before = source.vertices.copy()
    e = endpoint(position=(0, 0, 0), tangent=(np.sin(np.radians(42)), 0, np.cos(np.radians(42))), role='outlet')
    check = f.align_known_cap(e, source, {'vertices': 24}, cfg)
    assert check['before_axis_error_deg'] == pytest.approx(42)
    np.testing.assert_allclose(e.tangent, [0, 0, 1], atol=1e-12)
    assert e.radius == .6 and len(e.cap_ring) == 24
    np.testing.assert_array_equal(source.vertices, before)
    assert check['collar_contained'] and check['collar_outside_source_volume_mm3'] <= 1e-6


def test_aligned_mesh_section_matches_original_ring_after_stl_round_trip(cfg, tmp_path):
    e, source = align_synthetic_cap(cfg)
    target = f.natural_bend_points(e, '+X', 8, cfg)[-1]
    target[0] = 30
    r = f.make_route(e, '+X', target, roomy_box(), cfg, 'ROUTE_B_CURVED', bend_radius=8)
    path = tmp_path / 'aligned.stl'
    r.mesh.export(path)
    r.mesh = trimesh.load_mesh(path, process=True)
    check = f.attachment_section_qc(r, cfg)
    assert check['passed'] and check['axis_error_deg'] < 1e-6
    assert check['contour_hausdorff_mm'] < 1e-5 and check['center_offset_mm'] < 1e-5
    core, qc = f.union_core(source, [r], cfg)
    assert qc['passed'] and qc['source_solid_preserved'] and core.euler_number == 2


def test_old_slanted_but_closed_port_fails_alignment_check(cfg):
    e, _ = align_synthetic_cap(cfg)
    wrong = endpoint(position=(0, 0, 0), tangent=tuple(e.attachment_qc['reference_swc_mean_tangent']), role='outlet')
    target = f.natural_bend_points(wrong, '+X', 8, cfg)[-1]
    target[0] = 30
    r = f.make_route(wrong, '+X', target, roomy_box(), cfg, 'ROUTE_B_CURVED', bend_radius=8)
    assert r.mesh.is_volume  # The previous checks alone would accept this.
    r.endpoint = e
    check = f.attachment_section_qc(r, cfg)
    assert not check['passed'] and check['axis_error_deg'] == pytest.approx(42)


@pytest.mark.parametrize('corruption', ['center', 'radius', 'vertex_count'])
def test_cap_correspondence_mismatch_is_not_silently_fitted(cfg, corruption):
    e = endpoint(position=(0, 0, 0), tangent=(0, 0, 1), role='outlet')
    record = {'vertices': 24}
    if corruption == 'center':
        e.position += [0, 0, .1]
    elif corruption == 'radius':
        e.radius = .8
    else:
        record['vertices'] = 32
    with pytest.raises(ValueError, match='KNOWN_CAP_'):
        f.align_known_cap(e, known_cap_source(), record, cfg)


def test_real_bg001_o3_cap_provenance_and_alignment(cfg):
    inputs = f.load_inputs(cfg)
    e = next(e for e in inputs['endpoints'] if e.endpoint_id == 'O3')
    assert e.original_swc_id == 2127 and e.attachment_qc['before_axis_error_deg'] > 40
    assert e.radius == pytest.approx(.6311397091056972)
    assert e.attachment_qc['maximum_stl_plane_projection_mm'] < 1e-5
    assert e.attachment_qc['collar_contained']
    assert all(old.cap_ring is None for old in inputs['endpoints'] if old.endpoint_id != 'O3')
