from __future__ import annotations

from pathlib import Path

from utils.cfd_flow.validated_contract import (
    FULL_TIMESTEP_REFEREE_REVISION,
    PHYSICAL_FLUX_REVISION,
    PHYSICAL_PLANE_REVISION,
)
from utils.cfd_flow.validated_evidence import artifact_manifest, load_validated_evidence


ROOT = Path(__file__).resolve().parents[1]


def test_reference_loader_preserves_frozen_scientific_statuses() -> None:
    evidence = load_validated_evidence(ROOT)
    assert evidence["base"]["status"] == "CFD_FLOW_REPAIRED_BASE_TAU1_STEADY_PASS"
    assert evidence["coarse"]["status"] == "PASS"
    assert evidence["fine_resource_budget"]["scientific_failure"] is False
    assert evidence["fine_resource_budget"]["fine_contracts"]["fine_steady_state"] == "NOT_COMPLETED"
    assert evidence["plane_v3"]["revision"] == PHYSICAL_PLANE_REVISION
    assert evidence["plane_v3"]["flux_algorithm_revision"] == PHYSICAL_FLUX_REVISION
    assert evidence["full_v2_base"]["status"] == "PASS"
    assert evidence["full_v2_contract"]["referee_revision"] == FULL_TIMESTEP_REFEREE_REVISION


def test_reference_loader_exposes_results_without_python_result_constants() -> None:
    evidence = load_validated_evidence(ROOT)
    assert evidence["base"]["steady_audit"]["status"] == "PASS_NON_REFEREE"
    assert evidence["resolution"]["status"] == "PASS_TWO_GRID_RESOLUTION_SENSITIVITY"
    assert evidence["resolution"]["formal_asymptotic_grid_convergence"] is False


def test_binary_provenance_hashes_match_accepted_evidence() -> None:
    evidence = load_validated_evidence(ROOT)
    artifacts = artifact_manifest(ROOT)
    assert artifacts["base_restart"]["sha256"] == evidence["base"]["accepted_restart"]["sha256"]
    assert artifacts["coarse_restart"]["sha256"] == evidence["coarse"]["accepted_restart"]["sha256"]
    assert artifacts["fine_resource_restart"]["sha256"] == evidence["fine_resource_budget"]["restart_integrity"]["payload_sha256"]
    assert artifacts["adaptive_patch"]["sha256"] == "90efce400bb4b6ad5ad22ddd6518ccbeb8ba8a8eece2f76f2109a26dae92758f"


def test_mesh_and_plane_hashes_match_base_runtime_contract() -> None:
    evidence = load_validated_evidence(ROOT)
    artifacts = artifact_manifest(ROOT)
    hashes = evidence["base"]["runtime_contract"]["mesh_hashes"]
    assert artifacts["base_mesh_elemlist"]["sha256"] == hashes["elemlist.lsb"]
    assert artifacts["base_mesh_bnd"]["sha256"] == hashes["bnd.lsb"]
    assert artifacts["base_mesh_qval"]["sha256"] == hashes["qval.lsb"]
    assert evidence["plane_v3"]["contract_sha256"] == evidence["base"]["runtime_contract"]["physical_plane_contract_sha256"]
