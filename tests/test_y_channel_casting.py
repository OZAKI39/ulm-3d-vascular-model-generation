"""Real labelled lumen checks, including independently computed dimensions."""
import copy
import json
import zipfile
from pathlib import Path
import numpy as np
import pytest
import trimesh
from scipy.spatial import cKDTree
from vascular_processing import y_channel_casting as y


@pytest.fixture(scope='module')
def real(tmp_path_factory):
    # Keep regression coverage for the original wall-preserving historical mode.
    cfg = y.load_config(y.ROOT/'config/history/before_D1p6_H9_20261012/y_channel_p1s_abs.yaml')
    out = tmp_path_factory.mktemp('y_real')
    parts, report = y.geometry(cfg, out)
    return cfg, parts, report, out


def test_source_units_and_actual_labels(real):
    cfg, _, _, _ = real
    src, labels, planes = y.import_source(cfg)
    assert set(labels) == {1, 2, 3, 4}
    assert src.volume == pytest.approx(148.61907072, abs=1e-4)
    assert src.extents[0] == pytest.approx(24.8923059553, abs=1e-5)


def test_no_original_wall_or_bifurcation_reconstruction(real):
    cfg, parts, report, _ = real
    src, labels, _ = y.import_source(cfg)
    wall = src.vertices[np.unique(src.faces[labels == 1])]
    wall += report['original_body_translation_mm']
    distances = cKDTree(parts['Y_lumen_core'].vertices).query(wall)[0]
    assert distances.max() < 1e-5


def test_actual_joint_weld_and_normal_continuity(real):
    _, _, report, _ = real
    assert len(report['ports']) == 3
    for p in report['ports']:
        assert p['ring_vertex_count'] >= 100
        assert p['actual_seam_max_normal_angle_deg'] < 3
        assert p['actual_seam_median_normal_angle_deg'] < 1


def test_ports_actual_top_planes_are_upward_and_diameter_preserved(real):
    cfg, parts, report, _ = real
    core = parts['Y_lumen_core']
    top = report['heights']['port_top_z_mm']
    top_faces = np.max(np.abs(core.triangles[:, :, 2]-top), axis=1) < 2e-5
    assert top_faces.sum() > 3
    assert np.all(core.face_normals[top_faces, 2] > .9999)
    centers = core.triangles_center[top_faces]
    areas = core.area_faces[top_faces]
    for p in report['ports']:
        inside = np.linalg.norm(centers[:, :2]-np.array(p['top_center_mm'])[:2], axis=1) < p['radius_mm']+.01
        equivalent = 2*np.sqrt(areas[inside].sum()/np.pi)
        assert equivalent == pytest.approx(p['nominal_diameter_mm'], abs=.005)


def test_floor_wall_base_and_clearance_measured_from_mesh(real):
    cfg, parts, report, _ = real
    assert parts['Y_lumen_core'].bounds[0, 2]-cfg['base_mm'] == pytest.approx(1.5, abs=1e-5)
    mold = parts['Y_casting_box_base_2mm']
    inner = np.array(report['cavity_xy_bounds_mm'])
    assert np.allclose(mold.bounds[0, :2], inner[0]-2, atol=2e-5)
    assert np.allclose(mold.bounds[1, :2], inner[1]+2, atol=2e-5)
    hits, _, _ = mold.ray.intersects_location([[0, 0, 1]], [[0, 0, 1]])
    assert hits[:, 2].min() == pytest.approx(2., abs=1e-5)


def test_single_closed_solids_no_support_in_pdms(real):
    _, parts, report, _ = real
    assert len(parts) == 5
    for k, m in parts.items():
        assert m.is_volume and len(m.split()) == 1
        if 'retainer' in k:
            assert m.bounds[0, 2] > report['heights']['fill_z_mm']+.39
    assert max(report['assembly_pair_intersection_volume_mm3'].values()) < 1e-4


def test_box_adapts_to_channel_and_floor(real, tmp_path):
    cfg, _, report, _ = real
    alt = copy.deepcopy(cfg)
    alt['pdms_floor_mm'] = 2.5
    alt['side_clearance_mm'] = 4.0
    _, other = y.geometry(alt, tmp_path)
    assert np.allclose(np.array(other['mold_outer_size_mm'])-report['mold_outer_size_mm'], [2, 2, 1], atol=2e-5)


def test_wrong_source_rejected(real):
    cfg, _, _, _ = real
    cfg = copy.deepcopy(cfg)
    cfg['source_sha256'] = '0'*64
    with pytest.raises(ValueError, match='SOURCE_SHA256_MISMATCH'): y.import_source(cfg)


def test_changed_floor_cannot_slice_stale_geometry(real, tmp_path):
    import yaml
    cfg, _, _, out = real
    altered = copy.deepcopy(cfg)
    altered['output'] = out
    altered['pdms_floor_mm'] = .5
    config = tmp_path/'changed.yaml'
    config.write_text(yaml.safe_dump({k: str(v) if isinstance(v, Path) else v for k, v in altered.items()}))
    with pytest.raises(ValueError, match='CONFIG_CHANGED_REGENERATE_GEOMETRY'):
        y.run(config, 'slice')


def test_multiple_cap_boundaries_rejected():
    with pytest.raises(ValueError, match='MULTIPLE_CAP_BOUNDARIES'):
        y.boundary_cycle(np.array([[0, 1, 2], [3, 4, 5]]))


def test_print_parts_and_casting_transforms_roundtrip(real):
    _, parts, _, out = real
    meshes, records = y.arrange(parts, out)
    for original, m, r in zip(parts.values(), meshes, records):
        m = m.copy()
        m.apply_transform(r['plate_to_casting_transform'])
        assert cKDTree(original.vertices).query(m.vertices)[0].max() < 1e-8
        assert r['plate_bounds_mm'][0, 2] == pytest.approx(0, abs=1e-6)


def test_native_project_has_five_independent_parts_and_core_support_only(real):
    from vascular_processing.all_abs_print_plate import read_objects
    cfg, parts, _, out = real
    meshes, records = y.arrange(parts, out)
    profiles = {k: {'name': k} for k in ['machine', 'process', 'filament']}
    path = out/'test.3mf'
    y.write_project(path, meshes, records, profiles, cfg)
    objs, _, nplates, crc = read_objects(path)
    assert crc and len(objs) == 5 and nplates == 1
    enabled = [o for o in objs if o['metadata']['enable_support'] == '1']
    assert len(enabled) == 1
    assert enabled[0]['metadata']['name'] == 'Y_lumen_core'
    assert enabled[0]['metadata']['support_type'] == 'normal(auto)'
    assert enabled[0]['metadata']['raft_layers'] == '2'
