"""Process-isolated orchestration for the official VMTK production filters."""

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
    raw_vtp: Path
    raw_stl: Path
    remeshed_open_vtp: Path
    remeshed_open_stl: Path
    capped_vtp: Path
    extension_request_json: Path
    extension_result_json: Path
    extension_stdout_log: Path
    extension_stderr_log: Path
    remesh_request_json: Path
    remesh_result_json: Path
    remesh_stdout_log: Path
    remesh_stderr_log: Path
    cap_request_json: Path
    cap_result_json: Path
    cap_stdout_log: Path
    cap_stderr_log: Path


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
    settings = config.entity_remesh
    return {
        "project_extension_definition": "5D",
        "verified_vmtk_length_definition": (
            "AdaptiveExtensionLength: extensionLength = meanRadius * ExtensionRatio"
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
        "entity_remesh": {
            "entity_array_name": settings.entity_array_name,
            "far_core_entity_id": settings.far_core_entity_id,
            "active_entity_id": settings.active_entity_id,
            "exclude_entity_ids": list(settings.exclude_entity_ids),
            "core_collar": {
                "mode": settings.core_collar.mode,
                "face_layers": settings.core_collar.face_layers,
            },
            "element_size_mode": settings.element_size_mode,
            "target_edge_length_um": settings.target_edge_length_um,
            "preserve_boundary_edges": settings.preserve_boundary_edges,
        },
    }


def _vmtk_process_environment(config: VmtkConfig) -> dict[str, str]:
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
    environment["PATH"] = os.pathsep.join(
        [*(str(path) for path in native_directories), environment.get("PATH", "")]
    )
    environment["VMTK_RUNTIME_PREFIX"] = str(prefix)
    return environment


def _invoke(
    *,
    config: VmtkConfig,
    tool_script: Path,
    request: dict[str, Any],
    request_path: Path,
    result_path: Path,
    stdout_path: Path,
    stderr_path: Path,
    required_outputs: tuple[Path, ...],
    failure_status: str,
) -> VmtkInvocationResult:
    if not config.environment_python.is_file() or not config.runtime_prefix.is_dir():
        raise SurfacePrepareError("VMTK_ENVIRONMENT_BLOCKED")
    if not tool_script.is_file():
        raise SurfacePrepareError(f"Missing VMTK runner: {tool_script}")
    request_path.write_text(json.dumps(request, indent=2), encoding="utf-8")
    command = (
        str(config.environment_python.resolve()),
        str(tool_script.resolve()),
        "--request",
        str(request_path.resolve()),
        "--result",
        str(result_path.resolve()),
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
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0 or not result_path.is_file():
        raise SurfacePrepareError(failure_status)
    runtime = json.loads(result_path.read_text(encoding="utf-8"))
    if runtime.get("status") != "PASS" or not all(
        path.is_file() for path in required_outputs
    ):
        raise SurfacePrepareError(failure_status)
    return VmtkInvocationResult(command, request, runtime)


def run_official_vmtk(
    *, config: VmtkConfig, paths: VmtkExchangePaths, tool_script: Path
) -> VmtkInvocationResult:
    """Run the sole supported TPS boundary-normal extension operation."""

    if config.extension_mode != "boundarynormal":
        raise SurfacePrepareError("INVALID_VMTK_POSTPROCESS_CONFIGURATION")
    request = {
        "operation": "extension",
        "input_surface_vtp": str(paths.open_surface_vtp.resolve()),
        "raw_vtp": str(paths.raw_vtp.resolve()),
        "raw_stl": str(paths.raw_stl.resolve()),
        "parameters": parameter_mapping(config),
    }
    return _invoke(
        config=config,
        tool_script=tool_script,
        request=request,
        request_path=paths.extension_request_json,
        result_path=paths.extension_result_json,
        stdout_path=paths.extension_stdout_log,
        stderr_path=paths.extension_stderr_log,
        required_outputs=(paths.raw_vtp, paths.raw_stl),
        failure_status="VMTK_TPS_EXTENSION_FAILED",
    )


def build_entity_remesh_request(
    *, config: VmtkConfig, paths: VmtkExchangePaths
) -> dict[str, Any]:
    settings = config.entity_remesh
    expected = sorted((settings.far_core_entity_id, settings.active_entity_id))
    excluded = sorted(settings.exclude_entity_ids)
    return {
        "operation": "entity_remesh",
        "raw_vtp": str(paths.raw_vtp.resolve()),
        "remeshed_open_vtp": str(paths.remeshed_open_vtp.resolve()),
        "remeshed_open_stl": str(paths.remeshed_open_stl.resolve()),
        "entity_array_name": settings.entity_array_name,
        "expected_entity_ids": expected,
        "excluded_entity_ids": excluded,
        "active_entity_ids": sorted(set(expected) - set(excluded)),
        "element_size_mode": settings.element_size_mode,
        "target_edge_length_um": settings.target_edge_length_um,
        "preserve_boundary_edges": settings.preserve_boundary_edges,
    }


def entity_remesh_official_vmtk(
    *, config: VmtkConfig, paths: VmtkExchangePaths, tool_script: Path
) -> VmtkInvocationResult:
    if not paths.raw_vtp.is_file():
        raise SurfacePrepareError("VMTK_ENVIRONMENT_BLOCKED")
    request = build_entity_remesh_request(config=config, paths=paths)
    result = _invoke(
        config=config,
        tool_script=tool_script,
        request=request,
        request_path=paths.remesh_request_json,
        result_path=paths.remesh_result_json,
        stdout_path=paths.remesh_stdout_log,
        stderr_path=paths.remesh_stderr_log,
        required_outputs=(paths.remeshed_open_vtp, paths.remeshed_open_stl),
        failure_status="VMTK_CROSS_SEAM_REMESH_FAILED",
    )
    runtime = result.runtime
    if (
        runtime.get("surface_remesher_called") is not True
        or runtime.get("global_surface_remeshing_performed") is not False
        or runtime.get("excluded_entity_ids") != [1]
        or runtime.get("active_entity_ids") != [2]
        or runtime.get("input_entity_ids") != [1, 2]
        or runtime.get("output_entity_ids") != [1, 2]
    ):
        raise SurfacePrepareError("VMTK_CROSS_SEAM_ENTITY_ASSIGNMENT_FAILED")
    return result


def cap_official_vmtk(
    *, config: VmtkConfig, paths: VmtkExchangePaths, tool_script: Path
) -> VmtkInvocationResult:
    if not paths.remeshed_open_vtp.is_file():
        raise SurfacePrepareError("VMTK_ENVIRONMENT_BLOCKED")
    request = {
        "operation": "cap",
        "open_vtp": str(paths.remeshed_open_vtp.resolve()),
        "capped_vtp": str(paths.capped_vtp.resolve()),
    }
    result = _invoke(
        config=config,
        tool_script=tool_script,
        request=request,
        request_path=paths.cap_request_json,
        result_path=paths.cap_result_json,
        stdout_path=paths.cap_stdout_log,
        stderr_path=paths.cap_stderr_log,
        required_outputs=(paths.capped_vtp,),
        failure_status="VMTK_SURFACE_CAP_FAILED",
    )
    if result.runtime.get("surface_remesher_called") is not False:
        raise SurfacePrepareError("VMTK_SURFACE_CAP_FAILED:remesher_called")
    return result


def exchange_paths(
    *, input_directory: Path, vmtk_directory: Path, geometry_directory: Path
) -> VmtkExchangePaths:
    return VmtkExchangePaths(
        open_surface_vtp=input_directory / "open_surface_um.vtp",
        raw_vtp=geometry_directory / "vmtk_boundarynormal_raw_um.vtp",
        raw_stl=geometry_directory / "vmtk_boundarynormal_raw_um.stl",
        remeshed_open_vtp=geometry_directory / "vmtk_boundarynormal_crossseam_remeshed_open_um.vtp",
        remeshed_open_stl=geometry_directory / "vmtk_boundarynormal_crossseam_remeshed_open_um.stl",
        capped_vtp=vmtk_directory / "vmtk_boundarynormal_capped_um.vtp",
        extension_request_json=vmtk_directory / "extension_request.json",
        extension_result_json=vmtk_directory / "extension_environment.json",
        extension_stdout_log=vmtk_directory / "extension_stdout.log",
        extension_stderr_log=vmtk_directory / "extension_stderr.log",
        remesh_request_json=vmtk_directory / "remesh_request.json",
        remesh_result_json=vmtk_directory / "remesh_environment.json",
        remesh_stdout_log=vmtk_directory / "remesh_stdout.log",
        remesh_stderr_log=vmtk_directory / "remesh_stderr.log",
        cap_request_json=vmtk_directory / "cap_request.json",
        cap_result_json=vmtk_directory / "cap_environment.json",
        cap_stdout_log=vmtk_directory / "cap_stdout.log",
        cap_stderr_log=vmtk_directory / "cap_stderr.log",
    )
