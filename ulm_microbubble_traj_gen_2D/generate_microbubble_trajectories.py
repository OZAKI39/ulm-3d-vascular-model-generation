"""Prepare blood flow, choose a target area, and generate microbubble paths.

The program can generate the paths directly, or it can first divide the vessel
network into target areas. After an area is chosen, the program continues with
the blood flow that was already calculated.
"""

from __future__ import annotations

import argparse
import locale
import os
import shutil
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

# Make the rest of the project available when this file is started on its own.
if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from ulm_microbubble_traj_gen_2D.utils.core.config import DEFAULT_CONFIG_PATH, load_config
from ulm_microbubble_traj_gen_2D.utils.io.vascular_io import (
    VascularGeometryMode,
    geometry_mode_from_model_directory_name,
    load_physics_input,
    planar_inlet_area_equivalent_thickness_um,
    validate_physics_input_geometry,
)
from ulm_microbubble_traj_gen_2D.utils.runtime.console_output import (
    print_banner,
    print_command,
    print_key_values,
    print_section,
    print_success,
    print_warning,
)
from ulm_microbubble_traj_gen_2D.utils.workflows.runner import run_generation
from ulm_microbubble_traj_gen_2D.utils.workflows.target_candidate_runner import (
    run_target_candidate_preparation,
)
from ulm_microbubble_traj_gen_2D.utils.workflows.volumetric_runner import (
    run_volumetric_generation,
)


def initialize_windows_c_compiler_environment() -> Path | None:
    """Prepare a complete x64 Visual C++ environment for DOLFINx JIT builds.

    FFCx compiles generated C code at runtime.  A normal PowerShell or IDE
    launch can find ``cl.exe`` through setuptools while still missing the
    Windows SDK headers or tools.  Import the environment emitted by
    ``vcvars64.bat`` and make setuptools reuse it instead of constructing a
    second, potentially incomplete, compiler environment.
    """

    if os.name != "nt":
        return None

    vcvars64 = None
    if not _windows_c_compiler_environment_ready():
        vcvars64 = _find_vcvars64()
        if vcvars64 is None:
            raise RuntimeError(
                "DOLFINx requires the x64 Visual C++ build environment on "
                "Windows, but vcvars64.bat could not be found. Install Visual "
                "Studio 2022 Build Tools with 'Desktop development with C++' "
                "and a Windows SDK."
            )

        # Pass one raw command line on Windows. ``subprocess`` would otherwise
        # escape the embedded batch-file quotes when serializing an argument
        # list, causing cmd.exe to treat the quotes as part of the filename.
        command = f'cmd.exe /d /c call "{vcvars64}" >nul && set'
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            timeout=60,
        )
        encoding = locale.getpreferredencoding(False) or "utf-8"
        stdout = completed.stdout.decode(encoding, errors="replace")
        stderr = completed.stderr.decode(encoding, errors="replace")
        if completed.returncode != 0:
            detail = (
                " ".join(stderr.split())
                or f"exit code {completed.returncode}"
            )
            raise RuntimeError(
                f"Failed to initialize the x64 Visual C++ environment from "
                f"{vcvars64}: {detail}"
            )

        imported = {}
        for line in stdout.splitlines():
            if "=" not in line or line.startswith("="):
                continue
            key, value = line.split("=", 1)
            if key:
                imported[key] = value
        os.environ.update(imported)

    missing = _missing_windows_c_compiler_components()
    if missing:
        raise RuntimeError(
            "The Visual C++ environment is incomplete; unavailable components: "
            f"{', '.join(missing)}. Repair the Visual Studio C++ Build Tools "
            "and Windows SDK installation."
        )

    # setuptools/distutils normally runs vcvarsall.bat again and replaces PATH
    # for compiler subprocesses.  Reusing this verified environment keeps the
    # Windows SDK bin directory visible when link.exe launches rc.exe/mt.exe.
    os.environ["DISTUTILS_USE_SDK"] = "1"
    os.environ["MSSdk"] = "1"
    return vcvars64


def _windows_c_compiler_environment_ready() -> bool:
    """Return whether the complete MSVC/Windows SDK toolchain is visible."""

    return not _missing_windows_c_compiler_components()


def _missing_windows_c_compiler_components() -> tuple[str, ...]:
    """List compiler or SDK components unavailable to FFCx child builds."""

    missing = [
        executable
        for executable in ("cl.exe", "link.exe", "rc.exe", "mt.exe")
        if shutil.which(executable) is None
    ]
    if not _environment_path_contains_file("INCLUDE", "io.h"):
        missing.append("Windows SDK UCRT header io.h")
    if not _environment_path_contains_file("LIB", "ucrt.lib"):
        missing.append("Windows SDK library ucrt.lib")
    if not _environment_path_contains_file("LIB", "kernel32.lib"):
        missing.append("Windows SDK library kernel32.lib")
    return tuple(missing)


def _environment_path_contains_file(variable: str, filename: str) -> bool:
    """Return whether a semicolon-separated environment path contains a file."""

    return any(
        (Path(entry.strip().strip('"')) / filename).is_file()
        for entry in os.environ.get(variable, "").split(os.pathsep)
        if entry.strip()
    )


def _find_vcvars64() -> Path | None:
    """Locate the newest installed Visual Studio x64 environment script."""

    candidates = []
    vs_install_dir = os.environ.get("VSINSTALLDIR")
    if vs_install_dir:
        candidates.append(
            Path(vs_install_dir) / "VC" / "Auxiliary" / "Build" / "vcvars64.bat"
        )

    program_files_x86 = os.environ.get("ProgramFiles(x86)")
    if program_files_x86:
        vswhere = (
            Path(program_files_x86)
            / "Microsoft Visual Studio"
            / "Installer"
            / "vswhere.exe"
        )
        if vswhere.is_file():
            located = subprocess.run(
                [
                    str(vswhere),
                    "-latest",
                    "-products",
                    "*",
                    "-requires",
                    "Microsoft.VisualStudio.Component.VC.Tools.x86.x64",
                    "-property",
                    "installationPath",
                ],
                check=False,
                capture_output=True,
                text=True,
                timeout=30,
            )
            if located.returncode == 0 and located.stdout.strip():
                candidates.append(
                    Path(located.stdout.strip())
                    / "VC"
                    / "Auxiliary"
                    / "Build"
                    / "vcvars64.bat"
                )

    for root in (
        Path("C:/Program Files/Microsoft Visual Studio/2022"),
        Path("C:/Program Files (x86)/Microsoft Visual Studio/2022"),
        Path("D:/Microsoft Visual Studio/2022"),
    ):
        for edition in ("BuildTools", "Community", "Professional", "Enterprise"):
            candidates.append(
                root / edition / "VC" / "Auxiliary" / "Build" / "vcvars64.bat"
            )

    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    return None


def parse_args() -> argparse.Namespace:
    """
    Read the options chosen when the program starts.
    """

    parser = argparse.ArgumentParser(description="Generate ULM microbubble trajectories from the accepted vascular flow field.")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG_PATH, 
                        help="Path to the field-based YAML configuration file.")
    parser.add_argument(
        "--vascular-model",
        "--model-dir",
        dest="vascular_model",
        type=Path,
        default=None,
        help=(
            "Override input.model_dir. Accepts a vessel-generation run "
            "directory, the generator output root (the newest recognized "
            "planar_2d or volumetric_3d run is selected), or either member of a matching "
            ".swc/.vessels.npz pair."
        ),
    )
    parser.add_argument("--quick-test", action="store_true",
                        help="Use quick-test overrides from the YAML file.")
    parser.add_argument("--prepare-target-candidates", action="store_true", 
                        help=("Create topology-based vessel-bed candidates, complete molecular-target "
                              "selection, and then continue trajectory generation."))
    parser.add_argument("--skip-render", action="store_true",
                        help=("Skip PyVista/VTK visualization rendering in trajectory-generation "
                              "mode. Numerical fields, trajectories, metadata, and configured "
                              "contact-pilot results are still generated."))
    parser.add_argument("--reuse-field-from", type=Path, default=None,
                        help=("Skip the CFD solve by reusing a strictly validated "
                              "velocity_and_wall_shear_field.npz or its result directory."))
    return parser.parse_args()


def resolve_vascular_model(
    source: str | Path,
    *,
    required_mode: VascularGeometryMode | None = None,
) -> Path:
    """Resolve one generator output to a dimension-validated run directory.

    ``vessel_generation.py`` stores each run in a timestamped child directory,
    while also copying only the SWC file to the stable output root.  This
    resolver accepts both layouts and always returns the child directory that
    contains the matching SWC and ``.vessels.npz`` transport export.  Solver
    dimension comes from the run-directory name and must agree with transport
    metadata and coordinate geometry.
    """

    source_path = Path(source).expanduser().resolve()
    if source_path.is_file():
        if (
            source_path.suffix.lower() != ".swc"
            and not source_path.name.lower().endswith(".vessels.npz")
        ):
            raise ValueError(
                "A vascular model file must end in .swc or .vessels.npz: "
                f"{source_path}"
            )
        candidates = _candidate_directories_for_vascular_file(source_path)
    elif source_path.is_dir():
        candidates = [source_path]
        candidates.extend(
            child
            for child in source_path.iterdir()
            if child.is_dir()
        )
    else:
        raise FileNotFoundError(f"Cannot find the vascular model input: {source_path}")

    accepted: list[tuple[Path, VascularGeometryMode]] = []
    rejected: list[str] = []
    for candidate in candidates:
        try:
            mode = _validate_vascular_model_directory(candidate)
        except (FileNotFoundError, ValueError) as exc:
            rejected.append(f"{candidate}: {exc}")
        else:
            if required_mode is None or mode == required_mode:
                accepted.append((candidate.resolve(), mode))
            else:
                rejected.append(
                    f"{candidate}: geometry mode is {mode}, expected {required_mode}."
                )

    accepted = list(dict.fromkeys(accepted))
    if not accepted:
        detail = "\n  ".join(rejected[:8])
        suffix = f"\nChecked candidates:\n  {detail}" if detail else ""
        mode_description = required_mode or "planar_2d or volumetric_3d"
        raise ValueError(
            f"No compatible {mode_description} vascular result was found for "
            f"{source_path}. A compatible result must contain one matching "
            ".swc/.vessels.npz pair and a dimension marker in its directory name."
            f"{suffix}"
        )

    if source_path.is_file():
        matching = [
            (candidate, mode)
            for candidate, mode in accepted
            if _directory_contains_named_vascular_file(candidate, source_path.name)
        ]
        if matching:
            accepted = matching

    selected, _ = max(
        accepted,
        key=lambda item: (item[0].stat().st_mtime_ns, item[0].name),
    )
    return selected


def resolve_planar_vascular_model(source: str | Path) -> Path:
    """Backward-compatible resolver restricted to planar generator runs."""

    return resolve_vascular_model(source, required_mode="planar_2d")


def _candidate_directories_for_vascular_file(source_path: Path) -> list[Path]:
    """Find run directories that may contain a requested SWC or transport NPZ."""

    name = source_path.name
    candidates = [source_path.parent]
    candidates.extend(
        child
        for child in source_path.parent.iterdir()
        if child.is_dir() and _directory_contains_named_vascular_file(child, name)
    )
    return candidates


def _directory_contains_named_vascular_file(directory: Path, name: str) -> bool:
    """Return whether ``directory`` contains the requested file or its pair."""

    if (directory / name).is_file():
        return True
    if name.endswith(".vessels.npz"):
        return (directory / name.removesuffix(".vessels.npz")).with_suffix(
            ".swc"
        ).is_file()
    if name.endswith(".swc"):
        return (directory / name).with_suffix(".vessels.npz").is_file()
    return False


def _validate_vascular_model_directory(
    model_dir: Path,
) -> VascularGeometryMode:
    """Validate one generator pair and return its strict geometry mode."""

    physics_input = load_physics_input(model_dir)
    return validate_physics_input_geometry(model_dir, physics_input)


def _replace_config_vascular_model(
    cfg,
    model_dir: Path,
    geometry_mode: VascularGeometryMode,
):
    """Bind the selected model and its dimension-specific runtime parameters."""

    resolved_raw = dict(cfg.raw)
    input_raw = dict(resolved_raw.get("input", {}))
    input_raw["model_dir"] = str(model_dir)
    resolved_raw["input"] = input_raw
    resolved_raw["_resolved_vascular_geometry_mode"] = geometry_mode
    resolved_raw["_resolved_flow_solver_paradigm"] = (
        "boundary_fitted_dolfinx_stokes_xz_2d"
        if geometry_mode == "planar_2d"
        else "volumetric_3d_placeholder_no_projection"
    )

    if geometry_mode != "planar_2d":
        return replace(cfg, model_dir=model_dir, raw=resolved_raw)

    physics_input = load_physics_input(model_dir)
    validate_physics_input_geometry(
        model_dir,
        physics_input,
        expected_mode="planar_2d",
    )
    configured_thickness_um = float(cfg.field.effective_thickness_um)
    flow_quantity = str(
        physics_input.vessel_metadata.get("flow_quantity", "")
    ).strip().lower()
    if flow_quantity == "planar_flux_per_unit_depth":
        resolved_thickness_um = configured_thickness_um
        thickness_mode = "configured_extrusion_depth_for_true_planar_flux"
    else:
        # Format-v2 planar exports represented circular 3-D tubes embedded in
        # X-Z. Preserve their historical area-equivalent conversion.
        resolved_thickness_um = planar_inlet_area_equivalent_thickness_um(
            physics_input.vessels
        )
        thickness_mode = "legacy_inlet_circular_area_over_planar_width"
    field_raw = dict(resolved_raw.get("field", {}))
    field_raw["effective_thickness_um"] = resolved_thickness_um
    resolved_raw["field"] = field_raw
    resolved_raw["_configured_effective_thickness_um"] = configured_thickness_um
    resolved_raw["_resolved_planar_effective_thickness"] = {
        "mode": thickness_mode,
        "value_um": resolved_thickness_um,
    }
    return replace(
        cfg,
        model_dir=model_dir,
        field=replace(
            cfg.field,
            effective_thickness_um=resolved_thickness_um,
        ),
        raw=resolved_raw,
    )


def main() -> None:
    """
    Prepare the flow, select a target if needed, and generate bubble paths.
    """

    # Read the saved settings and apply the shorter test settings if requested.
    args    = parse_args()
    cfg     = load_config(args.config, quick_test=args.quick_test)
    geometry_mode: VascularGeometryMode = "planar_2d"
    configured_model = getattr(cfg, "model_dir", None)
    if configured_model is not None:
        explicit_model = getattr(args, "vascular_model", None)
        requested_model = explicit_model if explicit_model is not None else configured_model
        try:
            resolved_model = resolve_vascular_model(requested_model)
        except (FileNotFoundError, ValueError):
            if explicit_model is not None:
                raise
            # Older trajectory configurations named one pre-format-v2 run.
            # When that saved run is no longer readable, discover the newest
            # compatible result of the same declared dimension beside it.
            configured_path = Path(configured_model).expanduser().resolve()
            if configured_path.exists():
                # An existing named result that fails dimension/metadata
                # validation is corrupt or mislabeled. Never hide that conflict
                # by silently selecting a different sibling run.
                raise
            if not configured_path.parent.is_dir():
                raise
            try:
                required_mode = geometry_mode_from_model_directory_name(
                    configured_path
                )
            except ValueError:
                required_mode = None
            resolved_model = resolve_vascular_model(
                configured_path.parent,
                required_mode=required_mode,
            )
            print_warning(
                "The configured vascular result is not compatible with the "
                "current transport format. Using the newest compatible "
                f"{required_mode or 'recognized'} result instead: {resolved_model}"
            )
        geometry_mode = geometry_mode_from_model_directory_name(resolved_model)
        cfg = _replace_config_vascular_model(
            cfg,
            resolved_model,
            geometry_mode,
        )

    # FFCx JIT compilation is part of the planar DOLFINx path.  A selected 3-D
    # placeholder must fail without initializing or implying use of the 2-D solver.
    if geometry_mode == "planar_2d":
        compiler_environment = initialize_windows_c_compiler_environment()
        if compiler_environment is not None:
            print(f"[Windows C compiler] Initialized from {compiler_environment}")

    print_banner(
        "ULM microbubble trajectory generation",
        subtitle=("Mode: target-candidate preparation" if args.prepare_target_candidates 
                  else f"Mode: {'quick test' if cfg.quick_test else 'full simulation'}"))
    print_key_values(
        [("Configuration", cfg.source_path),
         *(
             [
                 ("Vascular geometry mode", geometry_mode),
                 ("Vascular model", cfg.model_dir),
                 *(
                     [
                         (
                            "Planar extrusion depth",
                             f"{float(cfg.field.effective_thickness_um):.9g} um",
                         )
                     ]
                     if geometry_mode == "planar_2d"
                     else []
                 ),
             ]
             if getattr(cfg, "model_dir", None) is not None
             else []
         ),
         ("Result directory", cfg.output_dir),
        ]
    )

    # If requested, calculate the blood flow and divide the vessel network into
    # areas that can be chosen as the molecular target.
    reuse_field_from = args.reuse_field_from
    if geometry_mode == "volumetric_3d":
        # The placeholder validates the 3-D contract and raises before any
        # planar target, CFD, particle, reuse, or rendering operation.
        run_volumetric_generation(
            cfg,
            render_artifacts=not args.skip_render,
            reuse_field_from=reuse_field_from,
        )
        return

    if args.prepare_target_candidates:
        candidate_result = run_target_candidate_preparation(cfg)
        print_banner("Target-candidate preparation complete")
        print_success(
            "The accepted flow and selectable candidate vessel beds are ready. "
            "Particle trajectories will be generated after target selection."
        )
        print_section("Saved artifacts")
        print_key_values(
            [
                ("Output directory", candidate_result.output_dir),
                ("Velocity / wall-shear field", candidate_result.field_npz_path),
                ("Candidate catalog", candidate_result.candidate_npz_path),
                ("Human-readable catalog", candidate_result.candidate_json_path),
                ("Selector CFD field", candidate_result.final_flow_vti_path),
                ("Candidate count", candidate_result.candidate_count),
            ]
        )
        if candidate_result.automatic_target_mask_path is not None:
            print_section("Automatic synthetic molecular target")
            print_key_values(
                [
                    ("Selected target mask", candidate_result.automatic_target_mask_path),
                    (
                        "Selection audit report",
                        candidate_result.automatic_selection_report_path,
                    ),
                ]
            )

        # Use the automatically chosen area when one is available. Otherwise,
        # open the selection page and wait for the user to save an area.
        selected_target_path = candidate_result.automatic_target_mask_path
        if selected_target_path is None:
            selector_args = [
                sys.executable,
                "-m",
                "ulm_microbubble_traj_gen.utils.visualization.apps.trame_molecular_target_selector",
                "--result-dir",
                str(candidate_result.output_dir),
                "--open-browser",
                "--exit-after-save",
            ]
            print_section("Interactive molecular-target selection")
            print_command(subprocess.list2cmdline(selector_args))
            subprocess.run(selector_args, check=True)
            selected_target_path = (
                candidate_result.output_dir / "selected_molecular_target_mask.npz"
            )
        if not selected_target_path.is_file():
            raise RuntimeError(
                "Target selection finished without saving "
                "selected_molecular_target_mask.npz."
            )

        # Record the chosen target in both the current settings and the copy of
        # the settings that will be saved with the results.
        target_raw = dict(cfg.raw.get("molecular_target", {}))
        target_raw.update(
            {
                "enabled": True,
                "region_mode": "mask_npz",
                "mask_npz_path": str(selected_target_path),
            }
        )
        resolved_raw = dict(cfg.raw)
        resolved_raw["molecular_target"] = target_raw
        cfg = replace(
            cfg,
            raw=resolved_raw,
            molecular_target=replace(
                cfg.molecular_target,
                enabled=True,
                region_mode="mask_npz",
                mask_npz_path=selected_target_path,
            ),
        )
        # The blood flow was already calculated while preparing the target
        # areas, so use it again instead of repeating the calculation.
        reuse_field_from = candidate_result.field_npz_path
        print_success(
            "Target selection is complete. Continuing with particle trajectory generation."
        )

    # Move the microbubbles through the vessels and create pictures if requested.
    result = run_generation(
        cfg,
        render_artifacts=not args.skip_render,
        reuse_field_from=reuse_field_from,
    )
    print_banner("Generation complete")

    # Explain whether the blood flow was calculated now or loaded from before.
    if args.skip_render:
        if reuse_field_from is not None:
            print_success(
                "The validated field reuse, particle transport, and numerical "
                "saving completed; visualization rendering was skipped."
            )
        else:
            print_success(
                "The flow solve, particle transport, and numerical saving completed; "
                "visualization rendering was skipped."
            )
    else:
        if reuse_field_from is not None:
            print_success(
                "The validated field reuse, particle transport, saving, and "
                "required rendering all completed."
            )
        else:
            print_success(
                "The flow solve, particle transport, saving, and required "
                "rendering all completed."
            )

    # List only the result files that were actually created.
    output_items = [
        ("Output directory", result["output_dir"]),
        ("Velocity / wall-shear field", result["field_npz_path"]),
        ("Particle trajectories", result["trajectories_npz_path"]),
        ("Domain metadata", result["metadata_path"]),
    ]

    if result["molecular_target_npz_path"] is not None:
        output_items.append(("Molecular target field", result["molecular_target_npz_path"]))

    if result["red_blood_cell_transport_npz_path"] is not None:
        output_items.append(("Red-blood-cell transport", result["red_blood_cell_transport_npz_path"]))

    output_items.extend(("Molecular contact pilot", path) for path in result["molecular_contact_pilot_paths"])
    rendered_outputs = [
        ("Initial CFD scene", result["initial_flow_html_path"]),
        ("Converged CFD scene", result["final_flow_html_path"]),
        ("Final wall-shear scene", result["final_wall_shear_html_path"]),
    ]
    output_items.extend((label, path) for label, path in rendered_outputs if path is not None)
    print_section("Saved artifacts")
    print_key_values(output_items)

    # If flow pictures were created, also show how to explore the final blood
    # flow and wall shear in a separate viewing page.
    if result["final_flow_html_path"] is not None:
        flow_probe_command = (
            f'"{sys.executable}" -m ulm_microbubble_traj_gen.utils.visualization.apps.trame_flow_viewer '
            f'--result-dir "{result["output_dir"]}" --stage final --open-browser'
        )
        wall_shear_probe_command = (
            f'"{sys.executable}" -m ulm_microbubble_traj_gen.utils.visualization.apps.trame_flow_viewer '
            f'--result-dir "{result["output_dir"]}" --stage final --view wall-shear --open-browser'
        )
        print_section("Live numerical viewers")
        print_key_values(
            [("Flow probe", flow_probe_command),
             ("Wall-shear probe", wall_shear_probe_command),
            ]
        )


if __name__ == "__main__":
    main()
