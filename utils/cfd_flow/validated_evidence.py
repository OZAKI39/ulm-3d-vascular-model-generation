"""Compact loader for immutable accepted CFD evidence.

Paths are repository-relative and results remain in committed JSON files; no
accepted numerical output is duplicated as a Python constant here.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .io import sha256_file


EVIDENCE_PATHS = {
    "base": "outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_base_anchor003274_20260901/qc/reference_scaled_base_final.json",
    "coarse": "outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_grid_convergence_anchor003274_20260901/qc/coarse_steady_acceptance.json",
    "resolution": "outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_grid_convergence_anchor003274_20260901/qc/coarse_base_resolution_sensitivity.json",
    "fine_resource_budget": "outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_grid_convergence_anchor003274_20260901/qc/resource_budget_termination.json",
    "geometry": "outputs/cfd_flow/healthy_mouse_capillary_dimensionless_qvalue_base_preflight_anchor003274_20260830/qc/final_base_geometry_validation.json",
    "plane_v3": "outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_grid_convergence_anchor003274_20260901/qc/physical_port_flux_plane_contract_v3.json",
    "adaptive_root_cause": "outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_grid_convergence_anchor003274_20260901/qc/fine_adaptive_target_375_376_forensics.json",
    "adaptive_fix": "outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_grid_convergence_anchor003274_20260901/qc/fine_adaptive_target_fix_validation.json",
    "fine_5000_safety": "outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_grid_convergence_anchor003274_20260901/fine_postfix/qc/fine_5000_safety_gate.json",
    "full_v2_base": "outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_base_anchor003274_20260901/qc/reference_scaled_base_full_referee.json",
    "full_v2_contract": "outputs/cfd_flow/healthy_mouse_capillary_tau1_base_anchor003274_20260830/qc/tau1_referee_v2_final.json",
}

ARTIFACT_PATHS = {
    "base_restart": "outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_base_anchor003274_20260901/accepted_restart/tau1_reference_scaled_base_598755.lsb",
    "coarse_restart": "outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_grid_convergence_anchor003274_20260901/coarse/accepted_restart/tau1_reference_scaled_coarse_354295.lsb",
    "fine_resource_restart": "outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_grid_convergence_anchor003274_20260901/fine_postfix/resource_budget_termination/restart/tau1_reference_scaled_fine_postfix_70313.lsb",
    "base_mesh_elemlist": "outputs/cfd_flow/healthy_mouse_capillary_dimensionless_qvalue_base_preflight_anchor003274_20260830/seeder/mesh/elemlist.lsb",
    "base_mesh_bnd": "outputs/cfd_flow/healthy_mouse_capillary_dimensionless_qvalue_base_preflight_anchor003274_20260830/seeder/mesh/bnd.lsb",
    "base_mesh_qval": "outputs/cfd_flow/healthy_mouse_capillary_dimensionless_qvalue_base_preflight_anchor003274_20260830/seeder/mesh/qval.lsb",
    "adaptive_patch": "patches/musubi/adaptive_flux_pressure.patch",
}


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"evidence must be a JSON object: {path}")
    return value


def load_validated_evidence(project_root: Path) -> dict[str, dict[str, Any]]:
    root = Path(project_root).resolve()
    result: dict[str, dict[str, Any]] = {}
    for topic, relative in EVIDENCE_PATHS.items():
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        result[topic] = _load_json(path)
    return result


def artifact_manifest(project_root: Path) -> dict[str, dict[str, Any]]:
    """Hash important payloads without decoding or printing binary content."""

    root = Path(project_root).resolve()
    result: dict[str, dict[str, Any]] = {}
    for role, relative in ARTIFACT_PATHS.items():
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        result[role] = {
            "path": relative,
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    return result
