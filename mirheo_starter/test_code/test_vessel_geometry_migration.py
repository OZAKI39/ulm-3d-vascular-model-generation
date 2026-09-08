"""Synthetic package-migration regressions; these do not accept the real vessel."""

from pathlib import Path
import json
import shutil
import subprocess
import sys
import unittest
from unittest.mock import patch

import numpy as np
import yaml

from test_code import test_vessel_geometry as fixtures
from py_scripts.vessel_geometry import migration
from py_scripts.vessel_geometry.export import load_package, run_import
from py_scripts.vessel_geometry.io import PROJECT_ROOT, load_config, read_json, sha256_file
from py_scripts.vessel_geometry.model import GeometryError
from py_scripts.vessel_geometry.validation import assert_same_geometry
from test_code.review_vessel_geometry import ensure_review


class MigrationTests(unittest.TestCase):
    def setUp(self):
        self.f = fixtures.GeometryTests()
        self.f.setUp()
        self.source = run_import(self.f.config_path, progress=None)
        self.path = self.f.root / "migrate.yaml"
        self.config = {
            "schema_version": 1, "mode": "accepted_package",
            "source_package": str(self.source.package_path), "source_run_id": self.source.run_id,
            "source_manifest_sha256": sha256_file(self.source.package_path / "import_manifest.json"),
            "source_acceptance": {"status": "SYNTHETIC_TEST_ONLY", "evidence": "Known tetrahedron fixture"},
            "output_root": str(self.f.root / "migrated"),
            "review_output_root": str(self.f.root / "migrated_review")}
        self.write_config()

    def tearDown(self):
        self.f.tearDown()

    def write_config(self):
        self.path.write_text(yaml.safe_dump(self.config), encoding="utf-8")

    def import_package(self):
        return run_import(self.path, progress=None)

    def test_accepted_package_is_byte_identical_and_roundtrip_exact(self):
        result = self.import_package()
        for source in self.source.package_path.rglob("*"):
            if source.is_file():
                target = result.package_path / source.relative_to(self.source.package_path)
                self.assertEqual(source.read_bytes(), target.read_bytes())
        self.assertEqual(assert_same_geometry(self.source.geometry, result.geometry)["status"], "PASS")
        self.assertEqual(result.manifest["states"]["source_model_historical_acceptance"], "SYNTHETIC_TEST_ONLY")
        self.assertEqual(result.geometry.checks["surface_processing"]["only_transform"].split(':')[0], "NONE")

    def test_repeat_reuses_package_and_html_without_touching_files(self):
        first, report = ensure_review(self.path, progress=None)
        paths = [p for root in (first.package_path, Path(report["html_path"]).parent) for p in root.rglob("*") if p.is_file()]
        before = {str(p): (sha256_file(p), p.stat().st_mtime_ns) for p in paths}
        second, second_report = ensure_review(self.path, progress=None)
        self.assertEqual(first.package_path, second.package_path)
        self.assertEqual(report, second_report)
        self.assertEqual(before, {str(p): (sha256_file(p), p.stat().st_mtime_ns) for p in paths})

    def test_config_change_creates_new_run_and_preserves_old(self):
        first = self.import_package()
        old = sha256_file(first.package_path / "migration_manifest.json")
        self.config["source_acceptance"]["evidence"] += " revised note"
        self.write_config()
        second = self.import_package()
        self.assertNotEqual(first.run_id, second.run_id)
        self.assertEqual(old, sha256_file(first.package_path / "migration_manifest.json"))

    def test_code_change_creates_new_run(self):
        first = self.import_package()
        actual = migration.code_identity()
        with patch.object(migration, "code_identity", return_value={**actual, "test_code_version": "changed"}):
            second = self.import_package()
        self.assertNotEqual(first.run_id, second.run_id)

    def test_loading_relocated_package_never_reads_source_paths(self):
        first = self.import_package()
        relocated = self.f.root / "standalone"
        shutil.copytree(first.package_path, relocated)
        # Remove only this test's synthetic source data, then deny any external read.
        shutil.rmtree(self.source.package_path)
        shutil.rmtree(self.f.source)
        code = '''
import sys
from pathlib import Path
from py_scripts.vessel_geometry.export import load_package
p=Path(sys.argv[1]).resolve()
allowed=(p,Path(sys.base_prefix)/'lib',Path(sys.prefix)/'lib')
def audit(event,args):
    if event=='open' and isinstance(args[0],str):
        path=Path(args[0]).resolve()
        if not any(path.is_relative_to(root) for root in allowed): raise AssertionError('External read: '+str(path))
sys.addaudithook(audit)
g=load_package(p)
assert len(g.triangles)==4
assert g.checks['migration_status']=='PASS'
assert 'mirheo' not in sys.modules
print('STANDALONE_LOAD_PASS')
'''
        process = subprocess.run([sys.executable, "-B", "-c", code, str(relocated)],
            cwd=PROJECT_ROOT, text=True, capture_output=True, timeout=30)
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertIn("STANDALONE_LOAD_PASS", process.stdout)

    def test_original_input_changed_blocks_even_with_cached_result(self):
        self.import_package()
        (self.f.run / "surface.vtp").write_bytes(b"changed source")
        with self.assertRaises(GeometryError) as caught:
            self.import_package()
        self.assertEqual(caught.exception.status, "BLOCKED")

    def test_missing_package_does_not_call_reconstruction_or_solver(self):
        (self.source.package_path / "vessel_geometry.npz").unlink()
        with patch("subprocess.run", side_effect=AssertionError("external process forbidden")), \
             patch.object(migration, "validate_accepted_package", side_effect=AssertionError("cannot fabricate geometry")):
            with self.assertRaises(GeometryError) as caught:
                self.import_package()
        self.assertEqual(caught.exception.status, "BLOCKED")
        failed = Path(caught.exception.details["package_path"])
        self.assertFalse((failed / "vessel_geometry.npz").exists())
        self.assertEqual(read_json(failed / "migration_failure.json")["status"], "BLOCKED")

    def test_lfs_pointer_package_is_blocked(self):
        (self.source.package_path / "vessel_geometry.npz").write_text("version https://git-lfs.github.com/spec/v1\noid sha256:abcd\nsize 30\n")
        with self.assertRaisesRegex(GeometryError, "LFS_POINTER"):
            self.import_package()

    def test_pin_mismatch_and_run_id_mismatch_are_blocked(self):
        self.config["source_manifest_sha256"] = "0" * 64
        self.write_config()
        with self.assertRaisesRegex(GeometryError, "SOURCE_MISMATCH"):
            self.import_package()
        self.config["source_run_id"] = "other_run"
        self.write_config()
        with self.assertRaisesRegex(GeometryError, "SOURCE_MISMATCH"):
            load_config(self.path)

    def test_destination_symlink_to_source_is_blocked(self):
        link = self.f.root / "unsafe"
        link.symlink_to(self.source.package_path, target_is_directory=True)
        self.config["output_root"] = str(link / "output")
        self.write_config()
        with self.assertRaisesRegex(GeometryError, "UNSAFE_OUTPUT"):
            self.import_package()
        self.assertFalse((self.source.package_path / "output").exists())

    def test_target_tampering_never_silently_reused(self):
        first = self.import_package()
        broken = first.package_path / "vessel_geometry.npz"
        broken.write_bytes(b"changed")
        with self.assertRaisesRegex(GeometryError, "PACKAGE_HASH_MISMATCH"):
            load_package(first.package_path)
        second = self.import_package()
        self.assertNotEqual(first.run_id, second.run_id)
        self.assertEqual(broken.read_bytes(), b"changed")

    def test_failed_or_unfinished_migration_cannot_inherit_old_pass(self):
        target = self.f.root / "partial"
        shutil.copytree(self.source.package_path, target)
        (target / "migration_intent.json").write_text('{}')
        with self.assertRaisesRegex(GeometryError, "PACKAGE_NOT_ACCEPTED"):
            load_package(target)
        (target / "migration_failure.json").write_text('{"status":"FAIL"}')
        with self.assertRaisesRegex(GeometryError, "PACKAGE_NOT_ACCEPTED"):
            load_package(target)

    def test_run_collision_never_adds_failure_to_old_result(self):
        first = self.import_package()
        self.config["source_acceptance"]["evidence"] += " changed"
        self.write_config()
        with patch("py_scripts.vessel_geometry.export.new_run_id", return_value=first.run_id):
            with self.assertRaises(FileExistsError):
                self.import_package()
        self.assertFalse((first.package_path / "migration_failure.json").exists())

    def test_multiple_existing_wall_labels_preserved(self):
        self.f.labels[2] = 91
        self.f.write_fixture()
        q = self.f.run / "identity.json"
        data = read_json(q)
        data["boundary_mapping"].pop("wall_entity_id")
        data["boundary_mapping"]["wall_entity_ids"] = [90, 91]
        q.write_text(json.dumps(data))
        source = run_import(self.f.config_path, progress=None)
        self.config.update(source_package=str(source.package_path), source_run_id=source.run_id,
                           source_manifest_sha256=sha256_file(source.package_path / "import_manifest.json"))
        self.write_config()
        result = self.import_package()
        self.assertEqual({p.entity_id for p in result.geometry.patches if p.kind == "wall"}, {90, 91})

    def test_json_exponent_overflow_rejected(self):
        path = self.f.root / "overflow.json"
        path.write_text('{"x":[1e999]}')
        with self.assertRaisesRegex(GeometryError, "NONFINITE_JSON"):
            read_json(path)

    def test_already_meter_package_not_converted_again(self):
        self.f.config["source_length_unit"] = "m"
        self.f.write_config()
        (self.f.run / "source.yaml").write_text("geometry:\n  input_unit: m\n")
        for row in self.f.rows:
            row["area_um2"] = str(float(row["area_um2"]) * 1e12)
        self.f.write_rows()
        identity_path = self.f.run / "identity.json"
        identity = read_json(identity_path)
        for row in identity["boundary_mapping"]["boundaries"]:
            row["area_um2"] = str(float(row["area_um2"]) * 1e12)
            row["centroid_um"] = (np.asarray(row["centroid_um"]) * 1e6).tolist()
        identity_path.write_text(json.dumps(identity))
        source = run_import(self.f.config_path, progress=None)
        self.config.update(source_package=str(source.package_path), source_run_id=source.run_id,
                           source_manifest_sha256=sha256_file(source.package_path / "import_manifest.json"))
        self.write_config()
        loaded = self.import_package().geometry
        self.assertEqual(loaded.to_meter, 1.0)
        np.testing.assert_array_equal(loaded.points_m, self.f.points.astype(float))
        np.testing.assert_array_equal(loaded.points_m, source.geometry.points_m)

    def test_import_and_help_cannot_import_mirheo_cuda_mpi_or_process_data(self):
        code = '''
import sys,runpy
from unittest.mock import patch
class Block:
    def find_spec(self,fullname,*args):
        if fullname.split('.')[0] in {'mirheo','mpi4py','cupy','torch','lammps'}:
            raise AssertionError('Forbidden import '+fullname)
sys.meta_path.insert(0,Block())
with patch('pathlib.Path.mkdir',side_effect=AssertionError('write')), patch('subprocess.run',side_effect=AssertionError('process')):
    import py_scripts.vessel_geometry.migration
    from py_scripts import import_vessel_geometry
    from test_code import review_vessel_geometry
    for module in (import_vessel_geometry,review_vessel_geometry):
        try: module.main(['--help'])
        except SystemExit as e: assert e.code==0
print('IMPORT_HELP_CPU_ONLY_PASS')
'''
        p = subprocess.run([sys.executable, "-B", "-c", code], cwd=PROJECT_ROOT,
                           text=True, capture_output=True, timeout=30)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("IMPORT_HELP_CPU_ONLY_PASS", p.stdout)


if __name__ == "__main__":
    unittest.main()
