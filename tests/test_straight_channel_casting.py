import copy
import json
from pathlib import Path
import numpy as np
import pyvista as pv
import pytest
import yaml
from scipy.spatial import cKDTree
from vascular_processing import straight_channel_casting as s


@pytest.fixture(scope='module')
def real(tmp_path_factory):
    # Preserve regression coverage for the archived side-gate configuration.
    cfg = s.load_config(s.ROOT/'config/history/before_closed_chamber_20261012/straight_channel_1p6_p1s_abs.yaml')
    out = tmp_path_factory.mktemp('straight_channel')
    parts, result = s.geometry(cfg, out)
    return cfg, parts, result, out


def test_all_parts_one_closed_solid(real):
    _, parts, _, _ = real
    assert len(parts) == 5
    for m in parts.values():
        assert m.is_volume and m.is_watertight and len(m.split()) == 1
        assert np.all(np.bincount(m.edges_unique_inverse) == 2)


def test_nominal_26mm_run_is_straight_and_1p6(real):
    cfg, parts, result, _ = real
    d = result['design']
    assert d['junction_x']-d['start_x'] == 26
    m = parts['Straight_lumen_core']
    for x in [-13.5, -10, -2, 4, 10, 16, 21]:
        assert s.section_diameter(m, [x, 0, d['center_z']], [1, 0, 0]) == pytest.approx(1.6, abs=.003)


def test_two_top_ports_and_one_horizontal_port(real):
    _, parts, report, _ = real
    m = parts['Straight_lumen_core']
    d = report['design']
    top = np.max(abs(m.triangles[:, :, 2]-d['top']), axis=1) < 1e-5
    assert (m.face_normals[top, 2] > .99999).all()
    assert len(np.unique(np.round(m.triangles_center[top, 0], 0))) >= 2
    end = np.max(abs(m.triangles[:, :, 0]-d['horizontal_end']), axis=1) < 1e-5
    assert end.sum() > 0 and (m.face_normals[end, 0] > .99999).all()
    assert 2*np.sqrt(m.area_faces[end].sum()/np.pi) == pytest.approx(1.6, abs=.003)


def test_floor_and_mold_dimensions(real):
    cfg, parts, result, _ = real
    core, mold = parts['Straight_lumen_core'], parts['Straight_casting_box_base_2mm']
    assert core.bounds[0, 2]-cfg['base_mm'] == pytest.approx(1.5, abs=.003)
    assert np.allclose(mold.extents, [40.6, 11.6, 13.4], atol=2e-5)
    assert result['pdms_fill_height_mm'] == pytest.approx(9.)
    assert mold.bounds[1, 2]-cfg['freeboard_mm']-min(intersection[2] for intersection in
        mold.ray.intersects_location([[0, 0, 1]], [[0, 0, 1]])[0]) == pytest.approx(9., abs=2e-5)
    intersections, _, _ = mold.ray.intersects_location([[0, 0, 1]], [[0, 0, 1]])
    assert min(intersections[:, 2]) == pytest.approx(2, abs=1e-6)
    assert np.allclose(mold.bounds[:, 1], [-5.8, 5.8], atol=1e-6)


def test_horizontal_rod_remains_natural_and_projects_4mm(real):
    _, parts, result, _ = real
    core, gate = parts['Straight_lumen_core'], parts['Straight_sidewall_gate']
    assert core.bounds[1, 0]-gate.bounds[1, 0] == pytest.approx(4, abs=2e-5)
    assert result['side_gate_external_seal_required'] is True
    assert result['side_gate_radial_gap_mm'] == .1


def test_gate_can_be_inserted_from_above(real):
    _, parts, result, _ = real
    for lift in [.3, 2, 6.5]:
        gate = parts['Straight_sidewall_gate'].copy()
        gate.apply_translation([0, 0, lift])
        assert s.volume_intersection(gate, parts['Straight_lumen_core']) < 1e-4
        assert s.volume_intersection(gate, parts['Straight_casting_box_base_2mm']) < 1e-4


def test_no_fixture_inside_pdms_or_part_penetration(real):
    _, parts, result, _ = real
    for key in ['Straight_retainer_INLET', 'Straight_retainer_OUTLET1']:
        assert parts[key].bounds[0, 2]-result['design']['fill'] == pytest.approx(.4, abs=2e-5)
    assert max(result['assembly_intersections_mm3'].values()) < 1e-4


def test_vtp_labels_and_explicit_mm_units(real):
    _, _, _, out = real
    m = pv.read(out/'labelled_straight_lumen_mm.vtp')
    assert set(m.cell_data['ModelFaceID']) == {1, 2, 3, 4}
    assert float(m.field_data['coordinate_scale_to_m'][0]) == .001
    assert m.bounds.x_max > 20


def test_print_transform_round_trip_and_bed_contact(real):
    _, parts, _, out = real
    meshes, records = s.arrange(parts, out)
    for original, printed, rec in zip(parts.values(), meshes, records):
        assert abs(printed.bounds[0, 2]) < 1e-7
        back = printed.copy()
        back.apply_transform(rec['plate_to_casting_transform'])
        assert cKDTree(original.vertices).query(back.vertices)[0].max() < 1e-8
    # The gate's thin X direction is vertical during printing, flange down.
    assert meshes[-1].extents[2] < 4


def test_core_only_support_and_five_objects(real):
    from vascular_processing.all_abs_print_plate import read_objects
    cfg, parts, _, out = real
    meshes, records = s.arrange(parts, out)
    profiles = {k: {'name': k} for k in ['machine', 'process', 'filament']}
    path = out/'test.3mf'
    s.write_project(path, meshes, records, profiles, cfg)
    objs, _, plates, crc = read_objects(path)
    assert crc and plates == 1 and len(objs) == 5
    enabled = [o for o in objs if o['metadata']['enable_support'] == '1']
    assert len(enabled) == 1 and enabled[0]['metadata']['name'] == 'Straight_lumen_core'


def test_modified_dimensions_cannot_slice_old_geometry(real, tmp_path):
    cfg, _, _, out = real
    cfg = copy.deepcopy(cfg)
    cfg['output'] = out
    cfg['pdms_floor_mm'] = .5
    file = tmp_path/'changed.yaml'
    file.write_text(yaml.safe_dump({k: str(v) if isinstance(v, Path) else v for k, v in cfg.items()}))
    with pytest.raises(ValueError, match='CONFIG_CHANGED_REGENERATE_GEOMETRY'): s.run(file, 'slice')


def test_adaptive_layout_keeps_horizontal_extension_outside_box(real):
    cfg, _, _, _ = real
    cfg = copy.deepcopy(cfg)
    cfg['side_clearance_mm'] = 4
    cfg['pdms_floor_mm'] = 2
    d = s.design(cfg)
    assert d['horizontal_end']-d['gate_back'] == pytest.approx(4)
    assert d['center_z']-cfg['diameter_mm']/2-cfg['base_mm'] == pytest.approx(2)
    assert d['inner'][1, 0] == pytest.approx(16.8)
