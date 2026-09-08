"""Synthetic tests only: passing these tests is not real-vessel acceptance."""

import contextlib
import csv
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import yaml

# Direct file execution puts test_code, rather than the project root, on sys.path.
if __name__ == "__main__" and __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from py_scripts.vessel_geometry.export import load_package, run_import
from py_scripts.vessel_geometry.io import (PROJECT_ROOT, load_config, read_json, read_tagged_vtp,
    require_file, resolve_recorded_path, sha256_file, unit_factor)
from py_scripts.vessel_geometry.model import GeometryError
from py_scripts.vessel_geometry.validation import build_geometry, topology_diagnostics, triangle_geometry, validate_arrays
from test_code.review_vessel_geometry import check_html, write_review_html


def write_vtp(path, points, faces, labels, *, polygon_size=3):
    from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    from vtkmodules.vtkCommonCore import vtkPoints
    from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
    from vtkmodules.vtkIOXML import vtkXMLPolyDataWriter
    data = vtkPolyData()
    vertices = vtkPoints(); vertices.SetData(numpy_to_vtk(points, deep=True)); data.SetPoints(vertices)
    cells = vtkCellArray()
    cells.SetData(numpy_to_vtkIdTypeArray(np.arange(0, len(faces) * polygon_size + 1, polygon_size, dtype=np.int64), deep=True),
                  numpy_to_vtkIdTypeArray(np.asarray(faces, dtype=np.int64).ravel(), deep=True))
    data.SetPolys(cells)
    if labels is not None:
        array = numpy_to_vtk(np.asarray(labels), deep=True); array.SetName("CellEntityIds"); data.GetCellData().AddArray(array)
    writer = vtkXMLPolyDataWriter(); writer.SetFileName(str(path)); writer.SetInputData(data)
    assert writer.Write() == 1


class GeometryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="synthetic_vessel_test_", dir=PROJECT_ROOT / "test_code")
        self.root = Path(self.temp.name)
        self.source = self.root / "source"
        self.run = self.source / "known_run"
        self.run.mkdir(parents=True)
        self.points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float32)
        self.faces = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]], dtype=np.int64)
        self.labels = np.array([90, 7, 90, 31], dtype=np.int64)
        self.rows = [dict(vmtk_cap_entity_id="31", port_id="out_b", role="ASSUMED_OUTLET", boundary_origin="TRUE_TERMINAL", boundary_index="12"),
                     dict(vmtk_cap_entity_id="7", port_id="in_a", role="ASSUMED_INLET", boundary_origin="CUT_PORT", boundary_index="3")]
        self.config = {"schema_version": 1, "source_project": str(self.source), "source_run": "known_run",
            "tagged_surface": "surface.vtp", "boundary_manifest": "ports.csv", "source_length_unit": "um",
            "source_unit_config": "source.yaml", "identity_qc": "identity.json",
            "output_root": str(self.root / "data"), "review_output_root": str(self.root / "review")}
        self.config_path = self.root / "import.yaml"
        self.write_fixture()

    def tearDown(self):
        self.temp.cleanup()

    def write_fixture(self):
        write_vtp(self.run / "surface.vtp", self.points, self.faces, self.labels)
        refs = []
        for row in self.rows:
            ids = np.flatnonzero(self.labels == int(row["vmtk_cap_entity_id"]))
            areas, centers, _ = triangle_geometry(self.points, self.faces[ids])
            row["triangle_count"] = str(len(ids))
            row["area_um2"] = str(areas.sum())
            refs.append({**row, "centroid_um": np.average(centers, axis=0, weights=areas).tolist()})
        self.write_rows()
        (self.run / "identity.json").write_text(json.dumps({"boundary_mapping": {"wall_entity_id": 90, "boundaries": refs}}))
        (self.run / "source.yaml").write_text("geometry:\n  input_unit: um\n")
        self.write_config()

    def write_rows(self):
        with (self.run / "ports.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(self.rows[0])); writer.writeheader(); writer.writerows(self.rows)

    def write_config(self):
        self.config_path.write_text(yaml.safe_dump(self.config))

    def build(self, **overrides):
        kwargs = dict(points=self.points, faces=self.faces, labels=self.labels, rows=self.rows,
                      source_length_unit="um", unit_source={"fixture": "SYNTHETIC"}, wall_entity_ids=[90])
        kwargs.update(overrides)
        return build_geometry(**kwargs)

    def test_synthetic_closed_surface_import_and_exact_roundtrip(self):
        result = run_import(self.config_path, progress=None)
        loaded = load_package(result.package_path)
        self.assertEqual(result.manifest["status"], "PASS")
        np.testing.assert_array_equal(loaded.points_source, self.points)
        np.testing.assert_array_equal(loaded.triangles, self.faces)
        np.testing.assert_array_equal(loaded.entity_ids, self.labels)
        self.assertTrue(loaded.checks["original_topology"]["closed"])
        self.assertTrue(loaded.checks["export_roundtrip"]["exact_array_and_identity_checks"]["patches"])

    def test_noncontiguous_unsorted_entity_identity(self):
        geometry = self.build()
        mapping = {p.entity_id: p.port_id for p in geometry.patches}
        self.assertEqual(mapping, {90: None, 7: "in_a", 31: "out_b"})
        self.assertEqual([p.original_role for p in geometry.patches if p.entity_id == 7], ["ASSUMED_INLET"])

    def test_variable_port_counts(self):
        for count in (1, 2, 3):
            with self.subTest(ports=count):
                labels = np.full(4, 90, dtype=np.int64)
                rows = []
                for i in range(count):
                    labels[i] = 11 + i * 17
                    rows.append(dict(entity_id=str(labels[i]), port_id=f"p{i}", role="ASSUMED_INLET" if i % 2 == 0 else "ASSUMED_OUTLET", boundary_origin="CUT_PORT"))
                geometry = self.build(labels=labels, rows=rows)
                self.assertEqual(sum(p.kind != "wall" for p in geometry.patches), count)

    def test_csv_row_order_independent(self):
        first = run_import(self.config_path, progress=None).geometry
        self.rows.reverse(); self.write_rows()
        second = run_import(self.config_path, progress=None).geometry
        self.assertEqual([(p.port_id, p.entity_id, p.face_ids) for p in first.patches],
                         [(p.port_id, p.entity_id, p.face_ids) for p in second.patches])

    def test_missing_cell_labels(self):
        write_vtp(self.run / "surface.vtp", self.points, self.faces, None)
        with self.assertRaisesRegex(GeometryError, "MISSING_LABELS"):
            read_tagged_vtp(self.run / "surface.vtp")

    def test_wrong_label_length(self):
        with self.assertRaisesRegex(GeometryError, "LABEL_LENGTH"):
            self.build(labels=self.labels[:-1])

    def test_illegal_face_indices(self):
        for index in (-1, 999):
            faces = self.faces.copy(); faces[0, 0] = index
            with self.assertRaisesRegex(GeometryError, "FACE_INDEX_OUT_OF_RANGE"):
                self.build(faces=faces)

    def test_noninteger_labels_not_truncated(self):
        for labels in (np.array([90, 7.2, 90, 31]), self.labels.astype(float)):
            with self.assertRaisesRegex(GeometryError, "NONINTEGER_LABELS"):
                self.build(labels=labels)

    def test_nan_and_inf_coordinates(self):
        for value in (float("nan"), float("inf")):
            points = self.points.copy(); points[0, 0] = value
            with self.assertRaisesRegex(GeometryError, "NONFINITE_COORDINATES"):
                self.build(points=points)

    def test_nontriangle_cell_rejected(self):
        write_vtp(self.run / "quad.vtp", self.points, [[0, 1, 2, 3]], [90], polygon_size=4)
        with self.assertRaisesRegex(GeometryError, "NON_TRIANGLE_CELLS"):
            read_tagged_vtp(self.run / "quad.vtp")

    def test_duplicate_mapping(self):
        with self.assertRaisesRegex(GeometryError, "DUPLICATE_MAPPING"):
            self.build(rows=self.rows + [self.rows[0]])
        for key in ("port_id", "boundary_index"):
            rows = [dict(row) for row in self.rows]; rows[1][key] = rows[0][key]
            with self.assertRaisesRegex(GeometryError, "DUPLICATE_MAPPING"):
                self.build(rows=rows)

    def test_conflicting_wall_and_unknown_labels(self):
        with self.assertRaisesRegex(GeometryError, "CONFLICTING_MAPPING"):
            self.build(wall_entity_ids=[7, 90])
        labels = self.labels.copy(); labels[0] = 999
        with self.assertRaisesRegex(GeometryError, "UNMAPPED_ENTITY"):
            self.build(labels=labels)

    def test_missing_wall_evidence(self):
        (self.run / "identity.json").write_text('{}')
        with self.assertRaisesRegex(GeometryError, "WALL_IDENTITY_MISSING"):
            run_import(self.config_path, progress=None)

    def test_unknown_unit_not_guessed(self):
        with self.assertRaisesRegex(GeometryError, "UNKNOWN_UNIT"):
            self.build(source_length_unit="unknown")
        self.config["source_length_unit"] = "pixels"; self.write_config()
        with self.assertRaisesRegex(GeometryError, "UNKNOWN_UNIT"):
            load_config(self.config_path)

    def test_meter_micrometer_coordinates_and_areas(self):
        rows = [{k: v for k, v in row.items() if k not in {"area_um2"}} for row in self.rows]
        microns = self.build(rows=rows)
        meters = self.build(points=self.points.astype(float) * 1e-6, rows=rows, source_length_unit="m")
        np.testing.assert_array_equal(microns.points_m, meters.points_m)
        np.testing.assert_allclose([p.area_m2 for p in microns.patches], [p.area_m2 for p in meters.patches], rtol=1e-14)
        self.assertAlmostEqual(sum(p.area_m2 for p in microns.patches) / 1e-12, 1.5 + np.sqrt(3) / 2)

    def test_unit_evidence_conflict(self):
        (self.run / "source.yaml").write_text("geometry:\n  input_unit: m\n")
        with self.assertRaisesRegex(GeometryError, "UNIT_EVIDENCE_CONFLICT"):
            run_import(self.config_path, progress=None)

    def test_lfs_pointer_recognized_and_blocked(self):
        path = self.run / "surface.vtp"
        path.write_text("version https://git-lfs.github.com/spec/v1\noid sha256:" + "0" * 64 + "\nsize 123\n")
        with self.assertRaisesRegex(GeometryError, "LFS_POINTER") as caught:
            run_import(self.config_path, progress=None)
        self.assertEqual(caught.exception.status, "BLOCKED")
        directory = Path(caught.exception.details["package_path"])
        self.assertFalse((directory / "source_surface.vtp").exists())
        self.assertFalse((directory / "import_manifest.json").exists())

    def test_source_hashes_unchanged(self):
        before = {str(p): sha256_file(p) for p in self.run.iterdir() if p.is_file()}
        result = run_import(self.config_path, progress=None)
        self.assertEqual(before, {str(p): sha256_file(p) for p in self.run.iterdir() if p.is_file()})
        self.assertEqual(sha256_file(result.package_path / "source_surface.vtp"), before[str(self.run / "surface.vtp")])

    def test_wall_open_edges_preserved(self):
        geometry = self.build()
        self.assertTrue(geometry.checks["original_topology"]["closed"])
        self.assertGreater(geometry.checks["patch_topology"]["90"]["boundary_edge_count"], 0)
        self.assertFalse(geometry.checks["surface_processing"]["holes_filled"])
        np.testing.assert_array_equal(geometry.triangles, self.faces)
        self.assertIsNone(next(p for p in geometry.patches if p.kind == "wall").outward_normal)

    def test_missing_input_no_reconstruction_or_success_model(self):
        (self.run / "surface.vtp").unlink()
        with patch("subprocess.run", side_effect=AssertionError("must not invoke external commands")):
            with self.assertRaisesRegex(GeometryError, "INPUT_MISSING") as caught:
                run_import(self.config_path, progress=None)
        directory = Path(caught.exception.details["package_path"])
        self.assertFalse((directory / "vessel_geometry.npz").exists())
        self.assertEqual(read_json(directory / "import_failure.json")["status"], "BLOCKED")

    def test_import_has_no_processing_or_writing_side_effects(self):
        code = '''
import sys
from unittest.mock import patch
with patch('subprocess.run', side_effect=AssertionError('external command')), patch('pathlib.Path.mkdir', side_effect=AssertionError('directory write')), patch('pathlib.Path.write_text', side_effect=AssertionError('file write')):
    import py_scripts.vessel_geometry.model
    import py_scripts.vessel_geometry.io
    import py_scripts.vessel_geometry.validation
    import py_scripts.vessel_geometry.export
    import py_scripts.import_vessel_geometry
assert 'lammps' not in sys.modules
assert not any(name.startswith('utils.cfd') for name in sys.modules)
print('NO_IMPORT_SIDE_EFFECTS')
'''
        result = subprocess.run([sys.executable, "-B", "-c", code], cwd=PROJECT_ROOT, capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("NO_IMPORT_SIDE_EFFECTS", result.stdout)

    def test_both_help_commands_do_no_import(self):
        for module in ("py_scripts.import_vessel_geometry", "test_code.review_vessel_geometry"):
            code = f"import importlib; from unittest.mock import patch\nwith patch('pathlib.Path.mkdir', side_effect=AssertionError('write')):\n m=importlib.import_module('{module}')\n try: m.main(['--help'])\n except SystemExit as e: assert e.code == 0\n"
            result = subprocess.run([sys.executable, "-B", "-c", code], cwd=PROJECT_ROOT, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("--config", result.stdout)

    def test_test_script_runs_directly_from_another_directory(self):
        # Select one existing test so the subprocess cannot recursively run this test.
        result = subprocess.run(
            [sys.executable, "-B", str(Path(__file__).resolve()),
             "GeometryTests.test_noncontiguous_unsorted_entity_identity"],
            cwd=self.root, capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Ran 1 test", result.stderr)
        self.assertIn("OK", result.stderr)

    def test_review_script_runs_directly_from_parent_directory(self):
        result = subprocess.run(
            [sys.executable, "-B", str(PROJECT_ROOT / "test_code" / "review_vessel_geometry.py"),
             "--config", str(self.config_path.relative_to(PROJECT_ROOT))],
            cwd=PROJECT_ROOT.parent, capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        reports = list((self.root / "review").glob("*/review_report.json"))
        self.assertEqual(len(reports), 1)
        report = read_json(reports[0])
        self.assertEqual(report["states"]["automatic_data_check"], "PASS")
        self.assertEqual(report["source_vs_package"]["status"], "PASS")
        self.assertTrue(Path(report["html_path"]).is_file())

    def test_output_symlink_to_source_blocked(self):
        link = self.root / "looks_safe"; link.symlink_to(self.source, target_is_directory=True)
        self.config["output_root"] = str(link / "outputs"); self.write_config()
        with self.assertRaisesRegex(GeometryError, "UNSAFE_OUTPUT"):
            load_config(self.config_path)
        self.assertFalse((self.source / "outputs").exists())

    def test_output_collision_does_not_overwrite(self):
        result = run_import(self.config_path, progress=None)
        digest = sha256_file(result.package_path / "import_manifest.json")
        with patch("py_scripts.vessel_geometry.export.new_run_id", return_value=result.run_id):
            with self.assertRaises(FileExistsError):
                run_import(self.config_path, progress=None)
        self.assertEqual(digest, sha256_file(result.package_path / "import_manifest.json"))

    def test_historical_path_resolution_is_explicit(self):
        (self.run / "boundaries").mkdir()
        path = self.run / "boundaries" / "port.stl"; path.write_bytes(b"synthetic path resolution fixture")
        recorded = r"E:\old\project\known_run\boundaries\port.stl"
        self.assertEqual(resolve_recorded_path(recorded, self.run), path)
        with self.assertRaisesRegex(GeometryError, "HISTORICAL_PATH_AMBIGUOUS"):
            resolve_recorded_path(r"E:\other_run\port.stl", self.run)

    def test_reference_area_count_and_center_failures(self):
        for key, value, code in (("triangle_count", "20", "COUNT_MISMATCH"), ("area_um2", "100", "AREA_MISMATCH")):
            rows = [dict(row) for row in self.rows]; rows[0][key] = value
            with self.assertRaisesRegex(GeometryError, code):
                self.build(rows=rows)
        reference = [{**row, "centroid_um": [99, 99, 99]} for row in self.rows]
        with self.assertRaisesRegex(GeometryError, "CENTER_MISMATCH"):
            self.build(identity_reference=reference)

    def test_open_surface_has_pending_normals_without_repair(self):
        rows = [dict(self.rows[1])]
        geometry = self.build(faces=self.faces[:3], labels=self.labels[:3], rows=rows)
        self.assertFalse(geometry.checks["original_topology"]["closed"])
        self.assertEqual(geometry.checks["original_topology"]["complete_surface_status"], "WARNING")
        self.assertEqual(geometry.patches[1].normal_status, "PENDING")
        self.assertIsNone(geometry.patches[1].outward_normal)
        self.assertEqual(len(geometry.triangles), 3)

    def test_inward_winding_not_flipped(self):
        faces = self.faces[:, ::-1].copy()
        geometry = self.build(faces=faces)
        self.assertLess(geometry.checks["original_topology"]["signed_volume_m3"], 0)
        np.testing.assert_array_equal(geometry.triangles, faces)
        self.assertFalse(geometry.checks["surface_processing"]["face_winding_changed"])

    def test_unused_points_preserved(self):
        points = np.vstack([self.points, [100, 200, 300]]).astype(np.float32)
        geometry = self.build(points=points)
        self.assertEqual(geometry.summary()["unused_point_count"], 1)
        np.testing.assert_array_equal(geometry.points_source, points)

    def test_vtp_cell_identity_conflict(self):
        values = np.array(["", "wrong", "", "out_b"])
        with self.assertRaisesRegex(GeometryError, "CELL_IDENTITY_CONFLICT"):
            self.build(cell_data={"port_id": values})

    def test_package_tampering_detected(self):
        result = run_import(self.config_path, progress=None)
        with (result.package_path / "boundary_patches.json").open("a") as stream:
            stream.write(" ")
        with self.assertRaisesRegex(GeometryError, "PACKAGE_HASH_MISMATCH"):
            load_package(result.package_path)

    def test_html_uses_complete_reloaded_mesh_offline(self):
        result = run_import(self.config_path, progress=None)
        geometry = load_package(result.package_path)
        html = self.root / "review.html"
        write_review_html(geometry, result.manifest, html,
                          {"package_path": str(result.package_path), "states": {"human_review": "PENDING"}})
        checks = check_html(html, geometry)
        self.assertEqual(checks["status"], "PASS")
        self.assertEqual(checks["rendered_triangle_count"], 4)
        self.assertEqual(checks["browser_execution"], "NOT_TESTED")

    def test_nonfinite_csv_error_report_is_valid_json(self):
        self.rows[0]["area_um2"] = "nan"; self.write_rows()
        with self.assertRaisesRegex(GeometryError, "AREA_MISMATCH") as caught:
            run_import(self.config_path, progress=None)
        failure = read_json(Path(caught.exception.details["package_path"]) / "import_failure.json")
        area_check = failure["details"]["comparison"]["checks"][-1]
        self.assertIsNone(area_check["expected"])
        self.assertTrue(failure["input_hashes_unchanged"])

    def test_missing_role_and_origin_fail(self):
        for key, code in (("role", "ROLE_MISSING"), ("boundary_origin", "PORT_IDENTITY_MISSING")):
            rows = [dict(row) for row in self.rows]; del rows[0][key]
            with self.assertRaisesRegex(GeometryError, code):
                self.build(rows=rows)

    def test_current_config_source_conflict_blocks_real_selection(self):
        selected = self.source / "current.yaml"
        selected.write_text("paths:\n  source_surface_run: a_different_run\n")
        self.config["source_selection_config"] = str(selected); self.write_config()
        with self.assertRaisesRegex(GeometryError, "SOURCE_AMBIGUOUS") as caught:
            run_import(self.config_path, progress=None)
        self.assertEqual(caught.exception.status, "BLOCKED")

    def test_duplicate_and_degenerate_diagnostics_do_not_repair(self):
        faces = np.vstack([self.faces, self.faces[0], [0, 0, 1]])
        before = faces.copy()
        checks = topology_diagnostics(self.points.astype(float), faces, 1e-28)
        self.assertEqual(checks["duplicate_triangle_count"], 1)
        self.assertEqual(checks["degenerate_triangle_count"], 1)
        self.assertGreater(checks["nonmanifold_edge_count"], 0)
        np.testing.assert_array_equal(faces, before)

    def test_legacy_and_area_weighted_centers_are_not_conflated(self):
        labels = np.array([7, 90, 90, 7])
        rows = [dict(entity_id="7", port_id="p", role="ASSUMED_INLET", boundary_origin="CUT_PORT")]
        reference = [{"vmtk_cap_entity_id": 7, "port_id": "p", "role": "ASSUMED_INLET", "boundary_origin": "CUT_PORT",
                      "centroid_um": self.points[self.faces[[0, 3]]].mean(axis=1).mean(axis=0).tolist()}]
        geometry = self.build(labels=labels, rows=rows, identity_reference=reference,
                              reference_center_method="triangle_centroid_mean_source_dtype")
        check = geometry.checks["boundary_comparisons"][1]["checks"][-1]
        self.assertEqual(check["absolute_error_m"], 0)
        self.assertGreater(check["different_definition_distance_m"], 0)

    def test_nonfinite_json_rejected(self):
        path = self.root / "bad.json"; path.write_text('{"x": NaN}')
        with self.assertRaisesRegex(GeometryError, "NONFINITE_JSON"):
            read_json(path)


if __name__ == "__main__":
    unittest.main()
