"""Actual mesh measurements for the requested D1.6 / PDMS >=9 mm revision."""
import copy
from pathlib import Path
import numpy as np
import pytest
import trimesh
import yaml
from vascular_processing import y_channel_casting as y
from vascular_processing import straight_channel_casting as straight
from vascular_processing.y_uniform_channel import local_section_diameter, build


@pytest.fixture(scope='module')
def real(tmp_path_factory):
    cfg = y.load_config(y.ROOT/'config/history/before_bottomless_dovetail_20261012/y_channel_p1s_abs.yaml')
    out = tmp_path_factory.mktemp('uniform_y')
    parts, report = y.geometry(cfg, out)
    return cfg, parts, report, out


def test_raw_original_is_unchanged_and_reconstruction_explicit(real):
    cfg, _, report, _ = real
    assert y.digest(cfg['source']) == 'e14b06f1164be3beab261ec052dbe7edad6c2394d80077787cb318cbece15f4b'
    assert not report['source_wall_preserved']
    assert not report['global_geometric_scaling']
    assert [p['source_branch_length_mm'] for p in report['ports']] == pytest.approx([14, 12, 12])


def test_all_branch_bend_riser_sections_and_caps_are_1p6(real):
    _, parts, report, out = real
    mesh = trimesh.load_mesh(out/'assembly/Y_channel_without_positioning_collars_reference.stl')
    assert len(report['section_checks']) == 42
    for row in report['section_checks']:
        assert local_section_diameter(mesh, row['center_mm'], row['normal']) == pytest.approx(1.6, abs=.003)
    core = parts['Y_lumen_core']
    top = np.max(abs(core.triangles[:, :, 2]-report['heights']['port_top_z_mm']), axis=1) < 1e-5
    assert (core.face_normals[top, 2] > .99999).all()
    assert 2*np.sqrt(core.area_faces[top].sum()/(3*np.pi)) == pytest.approx(1.6, abs=.003)


def test_actual_inner_floor_to_fill_has_abs_allowance_and_bottom_cover_is_1p5(real):
    cfg, parts, report, _ = real
    mold = parts['Y_casting_box_base_2mm']
    hits = mold.ray.intersects_location([[0, 0, 1]], [[0, 0, 1]])[0]
    floor = hits[:, 2].min()
    assert floor == pytest.approx(2., abs=1e-5)
    fill = mold.bounds[1, 2]-cfg['freeboard_mm']
    assert fill-floor == pytest.approx(9.6, abs=1e-5)
    assert report['pdms_height_qc']['allowance']['budgeted_finished_lower_bound_mm'] == pytest.approx(9.108, abs=1e-5)
    assert parts['Y_lumen_core'].bounds[0, 2]-floor == pytest.approx(1.5, abs=.003)
    assert report['heights']['port_top_z_mm']-fill == pytest.approx(7.4, abs=1e-5)


def test_five_manifold_parts_and_fixtures_above_fill(real):
    _, parts, report, _ = real
    assert len(parts) == 5 and all(report['gates'].values())
    for name, mesh in parts.items():
        assert mesh.is_volume and len(mesh.split()) == 1
        assert np.all(np.bincount(mesh.edges_unique_inverse) == 2)
        if 'retainer' in name:
            assert mesh.bounds[0, 2]-report['heights']['fill_z_mm'] == pytest.approx(.4, abs=2e-5)


def test_y_pdms_width_already_above_9_is_not_shrunk(real):
    cfg,parts,report,_=real
    m=parts['Y_casting_box_base_2mm']
    hits=m.ray.intersects_location([[0,0,6]],[[0,1,0]])[0]
    assert min(hits[:,1])*2 == pytest.approx(22.59982096,abs=2e-5)
    assert report['pdms_width_qc']['measured_inner_width_mm'] > 9
    assert report['pdms_width_qc']['passed']
    assert np.allclose(m.extents,[41.59037971,26.59982109,14.0],atol=2e-5)


@pytest.mark.parametrize('module,name', [(y,'y_channel_p1s_abs.yaml'),(straight,'straight_channel_1p6_p1s_abs.yaml')])
@pytest.mark.parametrize('invalid', [8.,float('nan')])
def test_minimum_width_cannot_be_lowered_or_nonfinite(module,name,invalid,tmp_path):
    cfg=yaml.safe_load((y.ROOT/'config/history/before_bottomless_dovetail_20261012'/name).read_text());cfg['minimum_pdms_width_mm']=invalid
    p=tmp_path/'bad_width.yaml';p.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError,match='PDMS_WIDTH_AT_LEAST_9_MM_REQUIRED'):module.load_config(p)


def test_width_qc_rejects_outer_width_mistaken_for_pdms_width():
    # Outer width 11.6 mm exceeds 9, but its actual inner gap is only 7.6 mm.
    cfg={'base_mm':2.,'minimum_pdms_width_mm':9.}
    outer=y.box([[-10,-5.8,0],[10,5.8,13.4]])
    cavity=y.box([[-8,-3.8,2],[8,3.8,14]])
    m=trimesh.boolean.difference([outer,cavity],engine='manifold')
    qc=y.pdms_width_qc(m,[[-8,-3.8],[8,3.8]],cfg,11.)
    assert m.extents[1] > 9
    assert qc['measured_inner_width_mm'] == pytest.approx(7.6,abs=1e-5)
    assert not qc['passed']


def test_minimum_height_does_not_limit_larger_needed_clearance(real):
    cfg, _, _, _ = real
    _, _, planes = y.import_source(cfg)
    cfg = copy.deepcopy(cfg)
    cfg['pdms_floor_mm'] = 6.
    *_, heights, samples = build(cfg, planes)
    assert heights['fill_z_mm']-cfg['base_mm'] > 9
    assert heights['fill_z_mm']-(cfg['base_mm']+6+.8+2.4) == pytest.approx(1.5)


@pytest.mark.parametrize('module,name', [(y,'y_channel_p1s_abs.yaml'),(straight,'straight_channel_1p6_p1s_abs.yaml')])
def test_request_cannot_be_silently_lowered_below_9(module, name, tmp_path):
    cfg = yaml.safe_load((y.ROOT/'config/history/before_bottomless_dovetail_20261012'/name).read_text())
    cfg['minimum_pdms_height_mm'] = 8.
    p = tmp_path/'invalid.yaml'; p.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match='PDMS_HEIGHT_AT_LEAST_9_MM_REQUIRED'): module.load_config(p)


def test_new_y_can_be_arranged_as_five_separate_print_objects(real):
    from vascular_processing.all_abs_print_plate import read_objects
    cfg, parts, _, out = real
    meshes, records = y.arrange(parts, out)
    p = out/'new_y.3mf'
    y.write_project(p, meshes, records, {k: {'name':k} for k in ['machine','process','filament']}, cfg)
    objects, _, plates, crc = read_objects(p)
    assert crc and plates == 1 and len(objects) == 5
    enabled = [o for o in objects if o['metadata']['enable_support'] == '1']
    assert len(enabled) == 1 and enabled[0]['metadata']['name'] == 'Y_lumen_core'
