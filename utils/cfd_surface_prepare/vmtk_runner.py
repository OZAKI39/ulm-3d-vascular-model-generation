"""Main-environment orchestration for the official VMTK runner."""

from __future__ import annotations

import json
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .config import VmtkConfig
from .io import SurfacePrepareError


@dataclass(frozen=True, slots=True)
class VmtkExchangePaths:
    open_surface_vtp: Path
    centerlines_vtp: Path
    raw_vtp: Path
    raw_stl: Path
    remeshed_open_vtp: Path
    capped_vtp: Path
    request_json: Path
    promotion_request_json: Path
    result_json: Path
    promotion_result_json: Path
    stdout_log: Path
    stderr_log: Path
    promotion_stdout_log: Path
    promotion_stderr_log: Path


@dataclass(frozen=True, slots=True)
class VmtkInvocationResult:
    command: tuple[str, ...]
    request: dict[str, Any]
    runtime: dict[str, Any]


def _git_value(repository: Path, *arguments: str) -> str | None:
    if not (repository / ".git").exists():
        return None
    completed = subprocess.run(
        ["git", "-C", str(repository), *arguments],
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return completed.stdout.strip() if completed.returncode == 0 else None


def official_source_provenance(config: VmtkConfig) -> dict[str, Any]:
    """Record inspected source and the official commit behind runtime v1.5.0."""

    return {
        "repository": "https://github.com/vmtk/vmtk",
        "inspected_repository_path": str(config.official_repository),
        "inspected_repository_commit": _git_value(
            config.official_repository, "rev-parse", "HEAD"
        ),
        "runtime_release_tag": "v1.5.0",
        "runtime_release_tag_commit": "30d5d7cb8e607d153c208a9d7d39c9feb7985476",
        "runtime_package": "conda-forge vmtk 1.5.0 py310h556579a_14",
        "official_source_files_inspected": [
            "vmtkScripts/vmtkflowextensions.py",
            "vtkVmtk/ComputationalGeometry/vtkvmtkPolyDataFlowExtensionsFilter.cxx",
            "vmtkScripts/vmtksurfaceremeshing.py",
            "vmtkScripts/vmtksurfacecapper.py",
        ],
    }


def parameter_mapping(config: VmtkConfig) -> dict[str, Any]:
    return {
        "project_extension_definition": "5D",
        "verified_vmtk_length_definition": (
            "AdaptiveExtensionLength: extensionLength = meanRadius * ExtensionRatio"
        ),
        "verified_official_source_line": (
            "vtkvmtkPolyDataFlowExtensionsFilter.cxx:437-440"
        ),
        "interpolation_mode": config.interpolation_mode,
        "preserve_cross_section_shape": config.preserve_cross_section_shape,
        "extension_mode": config.extension_mode,
        "sigma": config.sigma,
        "transition_ratio": config.transition_ratio,
        "adaptive_extension_length": config.adaptive_extension_length,
        "extension_ratio": config.extension_ratio,
        "adaptive_extension_radius": config.adaptive_extension_radius,
        "adaptive_boundary_points": config.adaptive_boundary_points,
        "remesh_after_extension": config.remesh_after_extension,
        "automatic_fallback": False,
        "parameter_sweep": False,
        "custom_tps_implementation": False,
    }


def _vmtk_process_environment(config: VmtkConfig) -> dict[str, str]:
    """Expose the pinned VMTK binaries only to the pmp child process.

    VMTK 1.5.0 on Windows is built against VTK 9.2.6 and an HDF5 stack that
    cannot be solved into the existing pmp/FEniCS environment.  A process-local
    runtime overlay lets the pmp interpreter load that verified binary stack
    without replacing pmp's packages or activation settings.
    """

    prefix = config.runtime_prefix.resolve()
    site_packages = prefix / "Lib" / "site-packages"
    native_directories = (
        prefix,
        prefix / "Library" / "mingw-w64" / "bin",
        prefix / "Library" / "usr" / "bin",
        prefix / "Library" / "bin",
        prefix / "Scripts",
        prefix / "bin",
    )
    if not site_packages.is_dir() or not (prefix / "Library" / "bin").is_dir():
        raise SurfacePrepareError("VMTK_ENVIRONMENT_BLOCKED")
    environment = os.environ.copy()
    existing_pythonpath = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = os.pathsep.join(
        [str(site_packages), *([existing_pythonpath] if existing_pythonpath else [])]
    )
    existing_path = environment.get("PATH", "")
    environment["PATH"] = os.pathsep.join(
        [*(str(path) for path in native_directories), existing_path]
    )
    environment["VMTK_RUNTIME_PREFIX"] = str(prefix)
    return environment


def run_official_vmtk(
    *,
    config: VmtkConfig,
    paths: VmtkExchangePaths,
    tool_script: Path,
) -> VmtkInvocationResult:
    """Execute exactly one official TPS flow-extension candidate."""

    if not config.environment_python.is_file():
        raise SurfacePrepareError("VMTK_ENVIRONMENT_BLOCKED")
    if not config.runtime_prefix.is_dir():
        raise SurfacePrepareError("VMTK_ENVIRONMENT_BLOCKED")
    if not tool_script.is_file():
        raise SurfacePrepareError(f"Missing VMTK runner: {tool_script}")
    if config.extension_mode not in {"centerlinedirection", "boundarynormal"}:
        raise SurfacePrepareError("INVALID_VMTK_EXTENSION_MODE")
    parameters = parameter_mapping(config)
    request: dict[str, Any] = {
        "operation": "extension",
        "input_surface_vtp": str(paths.open_surface_vtp.resolve()),
        "raw_vtp": str(paths.raw_vtp.resolve()),
        "raw_stl": str(paths.raw_stl.resolve()),
        "parameters": parameters,
    }
    if config.extension_mode == "centerlinedirection":
        if not paths.centerlines_vtp.is_file():
            raise SurfacePrepareError("VMTK_ENVIRONMENT_BLOCKED:centerlines_missing")
        request["centerlines_vtp"] = str(paths.centerlines_vtp.resolve())
    paths.request_json.write_text(json.dumps(request, indent=2), encoding="utf-8")
    command = (
        str(config.environment_python.resolve()),
        str(tool_script.resolve()),
        "--request",
        str(paths.request_json.resolve()),
        "--result",
        str(paths.result_json.resolve()),
    )
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        env=_vmtk_process_environment(config),
    )
    paths.stdout_log.write_text(completed.stdout, encoding="utf-8")
    paths.stderr_log.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0 or not paths.result_json.is_file():
        raise SurfacePrepareError("VMTK_TPS_EXTENSION_FAILED")
    runtime = json.loads(paths.result_json.read_text(encoding="utf-8"))
    required = (paths.raw_vtp, paths.raw_stl)
    if runtime.get("status") != "PASS" or not all(path.is_file() for path in required):
        raise SurfacePrepareError("VMTK_TPS_EXTENSION_FAILED")
    return VmtkInvocationResult(command, request, runtime)


def promote_official_vmtk(
    *,
    config: VmtkConfig,
    paths: VmtkExchangePaths,
    tool_script: Path,
    target_edge_length_um: float,
) -> VmtkInvocationResult:
    """Remesh and cap a RAW candidate only after project-side hard QC passes."""

    if not config.environment_python.is_file() or not paths.raw_vtp.is_file():
        raise SurfacePrepareError("VMTK_ENVIRONMENT_BLOCKED")
    request: dict[str, Any] = {
        "operation": "remesh_cap",
        "raw_vtp": str(paths.raw_vtp.resolve()),
        "remeshed_open_vtp": str(paths.remeshed_open_vtp.resolve()),
        "capped_vtp": str(paths.capped_vtp.resolve()),
        "target_edge_length_um": float(target_edge_length_um),
    }
    paths.promotion_request_json.write_text(
        json.dumps(request, indent=2), encoding="utf-8"
    )
    command = (
        str(config.environment_python.resolve()),
        str(tool_script.resolve()),
        "--request",
        str(paths.promotion_request_json.resolve()),
        "--result",
        str(paths.promotion_result_json.resolve()),
    )
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        env=_vmtk_process_environment(config),
    )
    paths.promotion_stdout_log.write_text(completed.stdout, encoding="utf-8")
    paths.promotion_stderr_log.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0 or not paths.promotion_result_json.is_file():
        raise SurfacePrepareError("VMTK_TPS_EXTENSION_FAILED")
    runtime = json.loads(paths.promotion_result_json.read_text(encoding="utf-8"))
    required = (paths.remeshed_open_vtp, paths.capped_vtp)
    if runtime.get("status") != "PASS" or not all(path.is_file() for path in required):
        raise SurfacePrepareError("VMTK_TPS_EXTENSION_FAILED")
    return VmtkInvocationResult(command, request, runtime)


def exchange_paths(
    *,
    input_directory: Path,
    vmtk_directory: Path,
    geometry_directory: Path,
    extension_mode: str,
) -> VmtkExchangePaths:
    mode_label = (
        "boundarynormal" if extension_mode == "boundarynormal" else "centerline"
    )
    return VmtkExchangePaths(
        open_surface_vtp=input_directory / "open_surface_um.vtp",
        centerlines_vtp=input_directory / "centerlines_um.vtp",
        raw_vtp=geometry_directory / f"vmtk_{mode_label}_raw_um.vtp",
        raw_stl=geometry_directory / f"vmtk_{mode_label}_raw_um.stl",
        remeshed_open_vtp=vmtk_directory / f"vmtk_{mode_label}_remeshed_open_um.vtp",
        capped_vtp=vmtk_directory / f"vmtk_{mode_label}_remeshed_capped_um.vtp",
        request_json=vmtk_directory / "request.json",
        promotion_request_json=vmtk_directory / "promotion_request.json",
        result_json=vmtk_directory / "extension_environment.json",
        promotion_result_json=vmtk_directory / "promotion_environment.json",
        stdout_log=vmtk_directory / "stdout.log",
        stderr_log=vmtk_directory / "stderr.log",
        promotion_stdout_log=vmtk_directory / "promotion_stdout.log",
        promotion_stderr_log=vmtk_directory / "promotion_stderr.log",
    )
