from __future__ import annotations

import os
import tempfile
import unittest
from dataclasses import dataclass, replace
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np

from ulm_microbubble_traj_gen_2D import generate_microbubble_trajectories as cli
from ulm_vascular_model_generator.utils.core.models import Vessel
from ulm_vascular_model_generator.utils.io.vessel_transport_export import (
    write_vessel_transport_npz,
)
from ulm_microbubble_traj_gen_2D.utils.workflows.volumetric_runner import (
    VolumetricPipelineUnavailableError,
    run_volumetric_generation,
)
from ulm_microbubble_traj_gen_2D.utils.io.vascular_io import load_physics_input


class GenerationSkipRenderCliTests(unittest.TestCase):
    @staticmethod
    def _write_vascular_result(
        directory: Path,
        *,
        geometry_mode: str = "planar_2d",
        fixed_y_um: float = 1500.0,
        distal_y_um: float | None = None,
    ) -> tuple[Path, Path]:
        directory.mkdir(parents=True)
        swc_path = directory / "tree.swc"
        swc_path.write_text("# test SWC\n", encoding="utf-8")
        vessel_path = directory / "tree.vessels.npz"
        vessel = Vessel(
            vid=0,
            parent_id=-1,
            x_p=np.asarray([0.0, fixed_y_um, 0.0]),
            x_d=np.asarray(
                [
                    100.0,
                    fixed_y_um if distal_y_um is None else distal_y_um,
                    100.0,
                ]
            ),
            radius=10.0,
            flow_rate=1000.0,
            mean_velocity=1.0,
        )
        write_vessel_transport_npz(
            [vessel],
            vessel_path,
            metadata={
                "geometry_mode": geometry_mode,
                "flow_quantity": (
                    "planar_flux_per_unit_depth"
                    if geometry_mode == "planar_2d"
                    else "volume_flow"
                ),
            },
        )
        return swc_path, vessel_path

    def test_windows_compiler_initialization_imports_vcvars_environment(self) -> None:
        vcvars64 = Path(
            r"D:\Microsoft Visual Studio\2022\BuildTools"
            r"\VC\Auxiliary\Build\vcvars64.bat"
        )
        completed = SimpleNamespace(
            returncode=0,
            stdout=(
                b"INCLUDE=C:\\WindowsSDK\\ucrt\r\n"
                b"LIB=C:\\WindowsSDK\\lib\r\n"
                b"PATH=C:\\VisualCpp\\bin\r\n"
            ),
            stderr=b"",
        )
        with (
            mock.patch.object(cli.os, "name", "nt"),
            mock.patch.object(
                cli,
                "_windows_c_compiler_environment_ready",
                return_value=False,
            ),
            mock.patch.object(
                cli,
                "_missing_windows_c_compiler_components",
                return_value=(),
            ),
            mock.patch.object(cli, "_find_vcvars64", return_value=vcvars64),
            mock.patch.object(
                cli.subprocess, "run", return_value=completed
            ) as run,
            mock.patch.dict(cli.os.environ, {}, clear=False),
        ):
            source = cli.initialize_windows_c_compiler_environment()
            imported_include = cli.os.environ["INCLUDE"]
            imported_lib = cli.os.environ["LIB"]
            distutils_use_sdk = cli.os.environ["DISTUTILS_USE_SDK"]
            ms_sdk = cli.os.environ["MSSdk"]

        self.assertEqual(source, vcvars64)
        self.assertEqual(imported_include, r"C:\WindowsSDK\ucrt")
        self.assertEqual(imported_lib, r"C:\WindowsSDK\lib")
        self.assertEqual(distutils_use_sdk, "1")
        self.assertEqual(ms_sdk, "1")
        self.assertIn("vcvars64.bat", run.call_args.args[0])

    def test_ready_windows_environment_is_reused_by_distutils(self) -> None:
        with (
            mock.patch.object(cli.os, "name", "nt"),
            mock.patch.object(
                cli,
                "_windows_c_compiler_environment_ready",
                return_value=True,
            ),
            mock.patch.object(
                cli,
                "_missing_windows_c_compiler_components",
                return_value=(),
            ),
            mock.patch.object(cli.subprocess, "run") as run,
            mock.patch.dict(cli.os.environ, {}, clear=False),
        ):
            source = cli.initialize_windows_c_compiler_environment()
            distutils_use_sdk = cli.os.environ["DISTUTILS_USE_SDK"]
            ms_sdk = cli.os.environ["MSSdk"]

        self.assertIsNone(source)
        self.assertEqual(distutils_use_sdk, "1")
        self.assertEqual(ms_sdk, "1")
        run.assert_not_called()

    def test_windows_environment_reports_all_missing_sdk_tools(self) -> None:
        available = {
            "cl.exe": r"C:\MSVC\cl.exe",
            "link.exe": r"C:\MSVC\link.exe",
            "rc.exe": None,
            "mt.exe": None,
        }
        with (
            mock.patch.object(
                cli.shutil,
                "which",
                side_effect=lambda executable: available[executable],
            ),
            mock.patch.dict(
                cli.os.environ,
                {
                    "INCLUDE": r"C:\WindowsSDK\ucrt",
                    "LIB": r"C:\WindowsSDK\ucrt-lib;C:\WindowsSDK\um-lib",
                },
                clear=False,
            ),
            mock.patch.object(cli.Path, "is_file", return_value=True),
        ):
            missing = cli._missing_windows_c_compiler_components()

        self.assertEqual(missing, ("rc.exe", "mt.exe"))

    def test_windows_compiler_initialization_is_noop_off_windows(self) -> None:
        with (
            mock.patch.object(cli.os, "name", "posix"),
            mock.patch.object(cli.subprocess, "run") as run,
        ):
            source = cli.initialize_windows_c_compiler_environment()

        self.assertIsNone(source)
        run.assert_not_called()

    def test_skip_render_option_is_opt_in(self) -> None:
        with mock.patch("sys.argv", ["generate_microbubble_trajectories.py"]):
            default_args = cli.parse_args()
        with mock.patch(
            "sys.argv",
            ["generate_microbubble_trajectories.py", "--skip-render"],
        ):
            numerical_args = cli.parse_args()

        self.assertFalse(default_args.skip_render)
        self.assertTrue(numerical_args.skip_render)
        self.assertIsNone(default_args.reuse_field_from)

    def test_reuse_field_option_accepts_result_directory_or_npz(self) -> None:
        result_directory = Path("previous-result")
        field_path = Path("previous-result") / "velocity_and_wall_shear_field.npz"

        with mock.patch(
            "sys.argv",
            [
                "generate_microbubble_trajectories.py",
                "--reuse-field-from",
                str(result_directory),
            ],
        ):
            directory_args = cli.parse_args()
        with mock.patch(
            "sys.argv",
            [
                "generate_microbubble_trajectories.py",
                "--reuse-field-from",
                str(field_path),
            ],
        ):
            file_args = cli.parse_args()

        self.assertEqual(directory_args.reuse_field_from, result_directory)
        self.assertEqual(file_args.reuse_field_from, field_path)

    def test_vascular_model_option_accepts_model_dir_alias(self) -> None:
        model_path = Path("generated-planar-model")
        with mock.patch(
            "sys.argv",
            [
                "generate_microbubble_trajectories.py",
                "--model-dir",
                str(model_path),
            ],
        ):
            args = cli.parse_args()

        self.assertEqual(args.vascular_model, model_path)

    def test_generator_output_root_selects_newest_planar_result(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            older = root / "20260729_120000_000000_planar_2d_seed_105"
            newest = root / "20260729_130000_000000_planar_2d_seed_105"
            volumetric = root / "20260729_140000_000000_volumetric_3d_seed_105"
            self._write_vascular_result(older)
            self._write_vascular_result(newest)
            self._write_vascular_result(
                volumetric,
                geometry_mode="volumetric_3d",
                distal_y_um=1600.0,
            )
            os.utime(older, ns=(1_000_000_000, 1_000_000_000))
            os.utime(newest, ns=(2_000_000_000, 2_000_000_000))
            os.utime(volumetric, ns=(3_000_000_000, 3_000_000_000))

            selected = cli.resolve_planar_vascular_model(root)
            newest_any_dimension = cli.resolve_vascular_model(root)

        self.assertEqual(selected, newest.resolve())
        self.assertEqual(newest_any_dimension, volumetric.resolve())

    def test_stable_swc_resolves_matching_timestamped_run(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stable_swc = root / "tree.swc"
            stable_swc.write_text("# stable generator copy\n", encoding="utf-8")
            run_dir = root / "20260729_130000_000000_planar_2d_seed_105"
            self._write_vascular_result(run_dir)

            selected = cli.resolve_planar_vascular_model(stable_swc)

        self.assertEqual(selected, run_dir.resolve())

    def test_explicit_volumetric_result_is_dispatched_by_directory_name(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            model_dir = (
                Path(directory)
                / "20260729_105518_131762_volumetric_3d_seed_105"
            )
            self._write_vascular_result(
                model_dir,
                geometry_mode="volumetric_3d",
                distal_y_um=1600.0,
            )

            selected = cli.resolve_vascular_model(model_dir)
            with self.assertRaisesRegex(ValueError, "No compatible planar_2d"):
                cli.resolve_planar_vascular_model(model_dir)

        self.assertEqual(selected, model_dir.resolve())

    def test_volumetric_placeholder_stops_before_planar_pipeline(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            model_dir = (
                Path(directory)
                / "20260729_105518_131762_volumetric_3d_seed_105"
            )
            self._write_vascular_result(
                model_dir,
                geometry_mode="volumetric_3d",
                distal_y_um=1600.0,
            )

            with self.assertRaisesRegex(
                VolumetricPipelineUnavailableError,
                "not projected into the planar solver",
            ):
                run_volumetric_generation(
                    SimpleNamespace(model_dir=model_dir),
                    render_artifacts=False,
                )

    def test_directory_name_and_transport_metadata_must_agree(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            model_dir = (
                Path(directory)
                / "20260730_005205_527835_planar_2d_seed_105"
            )
            self._write_vascular_result(
                model_dir,
                geometry_mode="volumetric_3d",
                distal_y_um=1600.0,
            )

            with self.assertRaisesRegex(ValueError, "geometry conflict"):
                cli.resolve_vascular_model(model_dir)

    def test_true_planar_flux_uses_configured_extrusion_depth(self) -> None:
        @dataclass(frozen=True)
        class FieldSettings:
            effective_thickness_um: float

        @dataclass(frozen=True)
        class RuntimeSettings:
            raw: dict
            model_dir: Path
            field: FieldSettings

        with tempfile.TemporaryDirectory() as directory:
            model_dir = (
                Path(directory)
                / "20260730_005205_527835_planar_2d_seed_105"
            )
            self._write_vascular_result(model_dir)
            cfg = RuntimeSettings(
                raw={
                    "input": {"model_dir": str(model_dir)},
                    "field": {"effective_thickness_um": 1.0},
                },
                model_dir=model_dir,
                field=FieldSettings(effective_thickness_um=1.0),
            )

            resolved = cli._replace_config_vascular_model(
                cfg,
                model_dir.resolve(),
                "planar_2d",
            )

        self.assertAlmostEqual(
            resolved.field.effective_thickness_um,
            1.0,
        )
        self.assertEqual(
            resolved.raw["_resolved_vascular_geometry_mode"],
            "planar_2d",
        )
        self.assertEqual(
            resolved.raw["_resolved_planar_effective_thickness"]["mode"],
            "configured_extrusion_depth_for_true_planar_flux",
        )

    def test_true_planar_flux_adapter_preserves_solver_flux(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            model_dir = (
                Path(directory)
                / "20260730_005205_527835_planar_2d_seed_105"
            )
            self._write_vascular_result(model_dir)
            physics_input = load_physics_input(
                model_dir,
                planar_extrusion_depth_um=7.5,
            )

        self.assertAlmostEqual(physics_input.vessels[0].flow_rate, 7_500.0)
        self.assertAlmostEqual(
            physics_input.vessels[0].flow_rate / 7.5,
            1_000.0,
        )
        self.assertEqual(
            physics_input.vessel_metadata["source_flow_quantity"],
            "planar_flux_per_unit_depth",
        )

    def test_main_omits_render_outputs_and_viewer_commands(self) -> None:
        output_dir = Path("numerical-result")
        cfg = SimpleNamespace(
            source_path=Path("config.yaml"),
            quick_test=False,
            output_dir=output_dir,
        )
        result = {
            "output_dir": output_dir,
            "field_npz_path": output_dir / "velocity_and_wall_shear_field.npz",
            "trajectories_npz_path": output_dir / "microbubble_field_trajectories.npz",
            "metadata_path": output_dir / "domain_metadata.yaml",
            "molecular_target_npz_path": None,
            "molecular_contact_pilot_paths": (),
            "initial_flow_html_path": None,
            "final_flow_html_path": None,
            "final_wall_shear_html_path": None,
            "red_blood_cell_transport_npz_path": None,
        }
        args = SimpleNamespace(
            config=Path("config.yaml"),
            quick_test=False,
            prepare_target_candidates=False,
            skip_render=True,
            reuse_field_from=None,
        )

        with (
            mock.patch.object(
                cli, "initialize_windows_c_compiler_environment", return_value=None
            ),
            mock.patch.object(cli, "parse_args", return_value=args),
            mock.patch.object(cli, "load_config", return_value=cfg),
            mock.patch.object(cli, "run_generation", return_value=result) as run,
            mock.patch.object(cli, "print_banner"),
            mock.patch.object(cli, "print_success") as print_success,
            mock.patch.object(cli, "print_section") as print_section,
            mock.patch.object(cli, "print_key_values") as print_key_values,
        ):
            cli.main()

        run.assert_called_once_with(
            cfg, render_artifacts=False, reuse_field_from=None
        )
        self.assertIn("rendering was skipped", print_success.call_args.args[0])
        section_names = [call.args[0] for call in print_section.call_args_list]
        self.assertNotIn("Live numerical viewers", section_names)
        saved_items = print_key_values.call_args_list[-1].args[0]
        saved_labels = [label for label, _ in saved_items]
        self.assertNotIn("Initial CFD scene", saved_labels)
        self.assertNotIn("Converged CFD scene", saved_labels)
        self.assertNotIn("Final wall-shear scene", saved_labels)

    def test_main_forwards_explicit_reuse_field_source(self) -> None:
        output_dir = Path("reused-result")
        reuse_source = Path("validated-result")
        cfg = SimpleNamespace(
            source_path=Path("config.yaml"),
            quick_test=False,
            output_dir=output_dir,
        )
        result = {
            "output_dir": output_dir,
            "field_npz_path": output_dir / "velocity_and_wall_shear_field.npz",
            "trajectories_npz_path": output_dir / "microbubble_field_trajectories.npz",
            "metadata_path": output_dir / "domain_metadata.yaml",
            "molecular_target_npz_path": None,
            "molecular_contact_pilot_paths": (),
            "initial_flow_html_path": None,
            "final_flow_html_path": None,
            "final_wall_shear_html_path": None,
            "red_blood_cell_transport_npz_path": None,
        }
        args = SimpleNamespace(
            config=Path("config.yaml"),
            quick_test=False,
            prepare_target_candidates=False,
            skip_render=True,
            reuse_field_from=reuse_source,
        )

        with (
            mock.patch.object(
                cli, "initialize_windows_c_compiler_environment", return_value=None
            ),
            mock.patch.object(cli, "parse_args", return_value=args),
            mock.patch.object(cli, "load_config", return_value=cfg),
            mock.patch.object(cli, "run_generation", return_value=result) as run,
            mock.patch.object(cli, "print_banner"),
            mock.patch.object(cli, "print_success") as print_success,
            mock.patch.object(cli, "print_section"),
            mock.patch.object(cli, "print_key_values"),
        ):
            cli.main()

        run.assert_called_once_with(
            cfg,
            render_artifacts=False,
            reuse_field_from=reuse_source,
        )
        self.assertIn("validated field reuse", print_success.call_args.args[0])

    def test_candidate_preparation_waits_for_selection_and_continues_generation(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output_dir = Path(directory)
            selected_target_path = (
                output_dir / "selected_molecular_target_mask.npz"
            )
            cfg = replace(cli.load_config(), output_dir=output_dir)
            candidate_result = SimpleNamespace(
                output_dir=output_dir,
                field_npz_path=output_dir
                / "velocity_and_wall_shear_field.npz",
                candidate_npz_path=output_dir
                / "molecular_target_candidates.npz",
                candidate_json_path=output_dir
                / "molecular_target_candidates.json",
                final_flow_vti_path=output_dir / "final_flow_field.vti",
                candidate_count=3,
                automatic_target_mask_path=None,
                automatic_selection_report_path=None,
            )
            result = {
                "output_dir": output_dir,
                "field_npz_path": candidate_result.field_npz_path,
                "trajectories_npz_path": output_dir
                / "microbubble_field_trajectories.npz",
                "metadata_path": output_dir / "domain_metadata.yaml",
                "molecular_target_npz_path": output_dir
                / "molecular_target_field.npz",
                "molecular_contact_pilot_paths": (),
                "initial_flow_html_path": None,
                "final_flow_html_path": None,
                "final_wall_shear_html_path": None,
                "red_blood_cell_transport_npz_path": None,
            }
            args = SimpleNamespace(
                config=Path("config.yaml"),
                quick_test=False,
                prepare_target_candidates=True,
                skip_render=True,
                reuse_field_from=None,
            )

            def save_selection(command, check):
                selected_target_path.touch()

            with (
                mock.patch.object(
                    cli,
                    "initialize_windows_c_compiler_environment",
                    return_value=None,
                ),
                mock.patch.object(cli, "parse_args", return_value=args),
                mock.patch.object(cli, "load_config", return_value=cfg),
                mock.patch.object(
                    cli,
                    "run_target_candidate_preparation",
                    return_value=candidate_result,
                ),
                mock.patch.object(
                    cli.subprocess,
                    "run",
                    side_effect=save_selection,
                ) as selector,
                mock.patch.object(
                    cli, "run_generation", return_value=result
                ) as run,
                mock.patch.object(cli, "print_banner"),
                mock.patch.object(cli, "print_command"),
                mock.patch.object(cli, "print_success"),
                mock.patch.object(cli, "print_section"),
                mock.patch.object(cli, "print_key_values"),
            ):
                cli.main()

        selector_args = selector.call_args.args[0]
        self.assertIn("--exit-after-save", selector_args)
        continued_cfg = run.call_args.args[0]
        self.assertTrue(continued_cfg.molecular_target.enabled)
        self.assertEqual(continued_cfg.molecular_target.region_mode, "mask_npz")
        self.assertEqual(
            continued_cfg.molecular_target.mask_npz_path,
            selected_target_path,
        )
        self.assertEqual(
            run.call_args.kwargs["reuse_field_from"],
            candidate_result.field_npz_path,
        )
        self.assertFalse(run.call_args.kwargs["render_artifacts"])


if __name__ == "__main__":
    unittest.main()
