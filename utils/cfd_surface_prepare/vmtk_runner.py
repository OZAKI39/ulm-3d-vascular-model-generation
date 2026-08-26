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
    remeshed_open_stl: Path
    capped_vtp: Path
    request_json: Path
    promotion_request_json: Path
    result_json: Path
    promotion_result_json: Path
    stdout_log: Path
    stderr_log: Path
    promotion_stdout_log: Path
    promotion_stderr_log: Path
    entity_remesh_request_json: Path
    entity_remesh_result_json: Path
    entity_remesh_stdout_log: Path
    entity_remesh_stderr_log: Path


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
        "postprocess_mode": config.postprocess_mode,
        "remesh_after_extension": config.remesh_after_extension,
        "global_surface_remeshing_performed": False,
        "entity_aware_extension_remeshing_performed": config.remesh_after_extension,
        "entity_remesh": {
            "enabled": config.entity_remesh.enabled,
            "entity_array_name": config.entity_remesh.entity_array_name,
            "far_core_entity_id": config.entity_remesh.far_core_entity_id,
            "active_entity_id": config.entity_remesh.active_entity_id,
            "expected_entity_ids": [
                config.entity_remesh.far_core_entity_id,
                config.entity_remesh.active_entity_id,
            ],
            "exclude_entity_ids": list(config.entity_remesh.exclude_entity_ids),
            "active_entity_ids": [config.entity_remesh.active_entity_id],
            "core_collar": {
                "mode": config.entity_remesh.core_collar.mode,
                "face_layers": config.entity_remesh.core_collar.face_layers,
            },
            "element_size_mode": config.entity_remesh.element_size_mode,
            "target_edge_length_um": config.entity_remesh.target_edge_length_um,
            "preserve_boundary_edges": config.entity_remesh.preserve_boundary_edges,
        },
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
    """LEGACY_GLOBAL_REMESH_REFERENCE_ONLY: remesh and cap a RAW candidate."""

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


def cap_official_vmtk(
    *,
    config: VmtkConfig,
    paths: VmtkExchangePaths,
    tool_script: Path,
) -> VmtkInvocationResult:
    """Directly cap the selected open candidate with the official VMTK capper."""

    if config.postprocess_mode == "cap_only" and not config.remesh_after_extension:
        source_open_vtp = paths.raw_vtp
    elif (
        config.postprocess_mode == "cross_seam_active_collar_remesh_then_cap"
        and config.remesh_after_extension
    ):
        source_open_vtp = paths.remeshed_open_vtp
    else:
        raise SurfacePrepareError("INVALID_VMTK_POSTPROCESS_CONFIGURATION")
    if not config.environment_python.is_file() or not source_open_vtp.is_file():
        raise SurfacePrepareError("VMTK_ENVIRONMENT_BLOCKED")
    request: dict[str, Any] = {
        "operation": "cap_only",
        "raw_vtp": str(source_open_vtp.resolve()),
        "source_open_vtp": str(source_open_vtp.resolve()),
        "capped_vtp": str(paths.capped_vtp.resolve()),
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
        raise SurfacePrepareError("VMTK_RAW_DIRECT_CAP_FAILED")
    runtime = json.loads(paths.promotion_result_json.read_text(encoding="utf-8"))
    if runtime.get("status") != "PASS" or not paths.capped_vtp.is_file():
        raise SurfacePrepareError("VMTK_RAW_DIRECT_CAP_FAILED")
    if runtime.get("surface_remesher_called") is not False:
        raise SurfacePrepareError("VMTK_RAW_DIRECT_CAP_FAILED:remesher_called")
    return VmtkInvocationResult(command, request, runtime)


def entity_remesh_official_vmtk(
    *,
    config: VmtkConfig,
    paths: VmtkExchangePaths,
    tool_script: Path,
) -> VmtkInvocationResult:
    """Remesh the cross-seam active entity with only FAR_CORE excluded."""

    settings = config.entity_remesh
    valid = (
        config.postprocess_mode == "cross_seam_active_collar_remesh_then_cap"
        and config.remesh_after_extension
        and settings.enabled
        and settings.entity_array_name == "RemeshEntityId"
        and settings.far_core_entity_id == 1
        and settings.active_entity_id == 2
        and settings.exclude_entity_ids == (1,)
        and settings.core_collar.mode == "core_face_adjacency_layers"
        and settings.core_collar.face_layers == 2
        and settings.element_size_mode == "edgelength"
        and settings.target_edge_length_um == 0.25913916380971913
        and settings.preserve_boundary_edges
    )
    if not valid:
        raise SurfacePrepareError("INVALID_VMTK_POSTPROCESS_CONFIGURATION")
    if not config.environment_python.is_file() or not paths.raw_vtp.is_file():
        raise SurfacePrepareError("VMTK_ENVIRONMENT_BLOCKED")
    request = build_entity_remesh_request(config=config, paths=paths)
    paths.entity_remesh_request_json.write_text(
        json.dumps(request, indent=2), encoding="utf-8"
    )
    command = (
        str(config.environment_python.resolve()),
        str(tool_script.resolve()),
        "--request",
        str(paths.entity_remesh_request_json.resolve()),
        "--result",
        str(paths.entity_remesh_result_json.resolve()),
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
    paths.entity_remesh_stdout_log.write_text(completed.stdout, encoding="utf-8")
    paths.entity_remesh_stderr_log.write_text(completed.stderr, encoding="utf-8")
    required = (paths.remeshed_open_vtp, paths.remeshed_open_stl)
    if completed.returncode != 0 or not paths.entity_remesh_result_json.is_file():
        raise SurfacePrepareError("VMTK_ENTITY_REMESH_GEOMETRY_FAILED")
    runtime = json.loads(paths.entity_remesh_result_json.read_text(encoding="utf-8"))
    if runtime.get("status") != "PASS" or not all(path.is_file() for path in required):
        raise SurfacePrepareError("VMTK_ENTITY_REMESH_GEOMETRY_FAILED")
    if (
        runtime.get("surface_remesher_called") is not True
        or runtime.get("global_surface_remeshing_performed") is not False
        or runtime.get("excluded_entity_ids") != [1]
        or runtime.get("active_entity_ids") != [2]
        or runtime.get("input_entity_ids") != [1, 2]
        or runtime.get("output_entity_ids") != [1, 2]
    ):
        raise SurfacePrepareError("VMTK_ENTITY_ASSIGNMENT_FAILED")
    return VmtkInvocationResult(command, request, runtime)


def build_entity_remesh_request(
    *, config: VmtkConfig, paths: VmtkExchangePaths
) -> dict[str, Any]:
    """Build the generalized entity-remesh request without executing VMTK."""

    settings = config.entity_remesh
    expected_entity_ids = sorted(
        (settings.far_core_entity_id, settings.active_entity_id)
    )
    excluded_entity_ids = sorted(settings.exclude_entity_ids)
    active_entity_ids = sorted(
        set(expected_entity_ids) - set(excluded_entity_ids)
    )
    return {
        "operation": "entity_remesh",
        "raw_vtp": str(paths.raw_vtp.resolve()),
        "remeshed_open_vtp": str(paths.remeshed_open_vtp.resolve()),
        "remeshed_open_stl": str(paths.remeshed_open_stl.resolve()),
        "entity_array_name": settings.entity_array_name,
        "expected_entity_ids": expected_entity_ids,
        "excluded_entity_ids": excluded_entity_ids,
        "active_entity_ids": active_entity_ids,
        "element_size_mode": settings.element_size_mode,
        "target_edge_length_um": settings.target_edge_length_um,
        "preserve_boundary_edges": settings.preserve_boundary_edges,
    }


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
        remeshed_open_vtp=(
            geometry_directory
            / f"vmtk_{mode_label}_crossseam_remeshed_open_um.vtp"
        ),
        remeshed_open_stl=(
            geometry_directory
            / f"vmtk_{mode_label}_crossseam_remeshed_open_um.stl"
        ),
        capped_vtp=vmtk_directory / f"vmtk_{mode_label}_capped_um.vtp",
        request_json=vmtk_directory / "request.json",
        promotion_request_json=vmtk_directory / "cap_request.json",
        result_json=vmtk_directory / "extension_environment.json",
        promotion_result_json=vmtk_directory / "cap_environment.json",
        stdout_log=vmtk_directory / "stdout.log",
        stderr_log=vmtk_directory / "stderr.log",
        promotion_stdout_log=vmtk_directory / "cap_stdout.log",
        promotion_stderr_log=vmtk_directory / "cap_stderr.log",
        entity_remesh_request_json=vmtk_directory / "entity_remesh_request.json",
        entity_remesh_result_json=vmtk_directory / "entity_remesh_environment.json",
        entity_remesh_stdout_log=vmtk_directory / "entity_remesh_stdout.log",
        entity_remesh_stderr_log=vmtk_directory / "entity_remesh_stderr.log",
    )
