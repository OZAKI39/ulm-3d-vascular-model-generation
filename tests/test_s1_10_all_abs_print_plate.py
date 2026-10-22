"""Real six-object native Bambu project and layout regression checks."""
import copy
import json
import zipfile
import numpy as np
import pytest
import trimesh
from scipy.spatial import cKDTree
from vascular_processing.all_abs_print_plate import (
    FILENAME, default_output, layout, layout_audit, model_output,
    read_objects, validate_archive,
)
from vascular_processing.support_inputs import digest, verify


@pytest.fixture(scope='module')
def plate():
    meshes, records, hashes = layout()
    out = default_output()
    summary = json.loads((out/'single_plate_summary.json').read_text())
    return out, meshes, records, hashes, summary


def test_native_project_has_six_independent_objects_and_one_plate(plate):
    out, meshes, records, _, _ = plate
    assert validate_archive(out/FILENAME, meshes, records)['passed']
    objects, settings, count, crc = read_objects(out/FILENAME)
    assert count == 1 and crc and len(objects) == 6
    assert settings['printer_model'] == 'Bambu Lab X1 Carbon'
    assert settings['nozzle_diameter'] == ['0.2']
    assert settings['filament_type'] == ['ABS']
    assert settings['print_sequence'] == 'by layer'
    assert [o['metadata']['enable_support'] for o in objects] == ['0','1','0','0','0','0']


def test_native_core_matches_latest_assembly_geometry_after_translation(plate):
    _, meshes, _, _, _ = plate
    core = trimesh.load_mesh(model_output()/'vascular_core/assembly_placed_ABS_core.stl')
    core.apply_translation(meshes[1].bounds.mean(0)-core.bounds.mean(0))
    assert cKDTree(core.vertices).query(meshes[1].vertices)[0].max() < 2e-4
    assert cKDTree(meshes[1].vertices).query(core.vertices)[0].max() < 2e-4
    assert len(core.faces) == len(meshes[1].faces)


def test_reserved_brim_support_envelopes_separate_and_fit_bed(plate):
    _, meshes, records, _, _ = plate
    result = layout_audit(meshes, records)
    assert result['passed'] and result['minimum_reserved_gap_mm'] > 8


def test_overlap_is_rejected(plate):
    _, meshes, records, _, _ = plate
    bad = [m.copy() for m in meshes]
    bad[1].apply_translation(np.r_[bad[0].bounds.mean(0)[:2]-bad[1].bounds.mean(0)[:2],0])
    assert not layout_audit(bad, records)['passed']


def test_bed_exclusion_is_rejected(plate):
    _, meshes, records, _, _ = plate
    bad = [m.copy() for m in meshes]
    bad[2].apply_translation(np.r_[[15,15]-bad[2].bounds.mean(0)[:2],0])
    assert not layout_audit(bad, records)['passed']


def test_wrong_object_support_override_is_rejected(plate):
    out, meshes, records, _, _ = plate
    bad = copy.deepcopy(records)
    bad[1]['support_enabled'] = False
    assert not validate_archive(out/FILENAME, meshes, bad)['passed']


def test_export_is_native_slicer_archive_with_single_unmodified_gcode(plate):
    out, _, _, _, summary = plate
    assert digest(out/FILENAME) == summary['sha256'] == digest(summary['slice']['archive'])
    with zipfile.ZipFile(out/FILENAME) as z:
        files = [n for n in z.namelist() if n.endswith('.gcode')]
        assert files == ['Metadata/plate_1.gcode']
        assert z.read(files[0]) == (out/'bambu_all_abs_plate/plate_1.gcode').read_bytes()
    assert summary['slice']['slicer_success']
    assert summary['slice']['exit_code'] == 0
    assert not summary['slice']['warning_audit']['critical']


def test_real_toolpath_preserves_all_parts_and_core_only_support(plate):
    out, _, records, _, summary = plate
    audit = json.loads((out/'actual_toolpath_audit.json').read_text())
    assert audit['passed'] and audit == summary['toolpaths']
    assert audit['model_and_support_inside_bed'] and audit['no_excluded_zone_extrusion']
    for row, record in zip(audit['objects'], records):
        assert row['name'] == record['name']
        assert row['model_segment_count'] > 0
        assert row['sampled_reference_path_coverage'] >= .99
        assert (row['support_segments_near_object'] > 0) == record['support_enabled']


def test_source_models_remain_unchanged_and_no_printer_delivery(plate):
    out, _, _, hashes, summary = plate
    old = json.loads((out/'source_protection_start.json').read_text())
    assert hashes == old and verify(old)['all_unchanged']
    assert summary['source_unchanged'] and summary['no_printer_delivery']
    assert not summary['slice']['printer_upload'] and not summary['slice']['printer_start']
