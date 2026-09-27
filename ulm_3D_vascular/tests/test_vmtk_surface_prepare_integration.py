"""Static and synthetic checks for the single formal VMTK surface path."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from utils.cfd_surface_prepare.config import load_surface_prepare_config
from utils.cfd_surface_prepare.guarded_remesh import edges_form_simple_closed_loop
from utils.cfd_surface_prepare.io import BoundaryInput
from utils.cfd_surface_prepare.vmtk_pipeline import boundary_plane_alignment_pass
from utils.cfd_surface_prepare.vmtk_qc import (
    normalize_interface_diagnostic,
    symmetric_mesh_size_mismatch,
)
from utils.cfd_surface_prepare.vmtk_runner import (
    build_entity_remesh_request,
    exchange_paths,
    parameter_mapping,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = PROJECT_ROOT / "configs" / "cfd_surface_prepare.yaml"


@pytest.fixture(scope="module")
def formal_config():
    return load_surface_prepare_config(CONFIG_PATH, project_root=PROJECT_ROOT)


def _boundary(index: int = 0) -> BoundaryInput:
    return BoundaryInput(
        index=index,
        port_id=f"synthetic__{index}",
        boundary_origin="CUT_PORT",
        role="ASSUMED_INLET" if index == 0 else "ASSUMED_OUTLET",
        global_node_id=None,
        global_edge_id=index,
        center_um=np.zeros(3),
        source_radius_um=1.0,
        pressure_original_pa=100.0,
        expected_flow_m3_s=1.0e-15,
        simulation_tangent=np.asarray((0.0, 0.0, 1.0)),
        outward_normal=np.asarray((0.0, 0.0, 1.0)),
        extension_length_um=10.0,
        extension_end_um=np.asarray((0.0, 0.0, 10.0)),
    )


def test_formal_configuration_freezes_verified_scientific_parameters(formal_config):
    assert formal_config.backend.method == "vmtk_flowextensions"
    assert formal_config.vmtk.interpolation_mode == "thinplatespline"
    assert formal_config.vmtk.extension_mode == "boundarynormal"
    assert formal_config.vmtk.sigma == 1.0
    assert formal_config.vmtk.transition_ratio == 0.5
    assert formal_config.vmtk.extension_ratio == 10.0
    assert formal_config.qc.minimum_normal_dot == 0.999


def test_formal_configuration_freezes_cross_seam_contract(formal_config):
    settings = formal_config.vmtk.entity_remesh
    assert settings.entity_array_name == "RemeshEntityId"
    assert settings.far_core_entity_id == 1
    assert settings.active_entity_id == 2
    assert settings.exclude_entity_ids == (1,)
    assert settings.core_collar.face_layers == 2
    assert settings.target_edge_length_um == 0.25913916380971913


def test_removed_configuration_sections_are_rejected(tmp_path: Path):
    text = CONFIG_PATH.read_text(encoding="utf-8") + "\nruntime:\n  verbose: true\n"
    path = tmp_path / "invalid.yaml"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match="Unknown keys"):
        load_surface_prepare_config(path, project_root=PROJECT_ROOT)


def test_centerline_extension_mode_is_rejected(tmp_path: Path):
    text = CONFIG_PATH.read_text(encoding="utf-8").replace(
        "extension_mode: boundarynormal", "extension_mode: centerlinedirection"
    )
    path = tmp_path / "invalid.yaml"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match="extension_mode"):
        load_surface_prepare_config(path, project_root=PROJECT_ROOT)


def test_parameter_mapping_contains_only_formal_entity_contract(formal_config):
    mapping = parameter_mapping(formal_config.vmtk)
    assert mapping["extension_mode"] == "boundarynormal"
    assert mapping["entity_remesh"]["exclude_entity_ids"] == [1]
    assert "automatic_fallback" not in mapping
    assert "parameter_sweep" not in mapping


def test_exchange_paths_have_no_centerline_or_promotion_names(tmp_path: Path):
    paths = exchange_paths(
        input_directory=tmp_path / "input",
        vmtk_directory=tmp_path / "vmtk",
        geometry_directory=tmp_path / "geometry",
    )
    names = set(paths.__dataclass_fields__)
    assert "centerlines_vtp" not in names
    assert not any("promotion" in name for name in names)
    assert paths.raw_vtp.name == "vmtk_boundarynormal_raw_um.vtp"


def test_entity_remesh_request_is_far_core_excluded(formal_config, tmp_path: Path):
    paths = exchange_paths(
        input_directory=tmp_path / "input",
        vmtk_directory=tmp_path / "vmtk",
        geometry_directory=tmp_path / "geometry",
    )
    request = build_entity_remesh_request(config=formal_config.vmtk, paths=paths)
    assert request["operation"] == "entity_remesh"
    assert request["expected_entity_ids"] == [1, 2]
    assert request["excluded_entity_ids"] == [1]
    assert request["active_entity_ids"] == [2]


def test_runner_has_only_formal_operations():
    source = (PROJECT_ROOT / "tools" / "run_vmtk_flowextension.py").read_text(
        encoding="utf-8"
    )
    assert 'operation == "extension"' in source
    assert 'operation == "entity_remesh"' in source
    assert 'operation == "cap"' in source
    assert "centerlinedirection" not in source
    assert "remesh_cap" not in source
    assert "cap_only" not in source


def test_pipeline_has_no_resume_or_historical_run_dependency():
    source = (
        PROJECT_ROOT / "utils" / "cfd_surface_prepare" / "vmtk_pipeline.py"
    ).read_text(encoding="utf-8")
    assert "resume_crossseam" not in source
    assert "PREVIOUS_" not in source
    assert "GUARDED_ENTITY_REMESH_RUN_ID" not in source
    assert "seam_quality_comparison.csv" not in source


def test_all_open_hard_qc_precedes_cap():
    source = (
        PROJECT_ROOT / "utils" / "cfd_surface_prepare" / "vmtk_pipeline.py"
    ).read_text(encoding="utf-8")
    cap = source.index("cap = cap_official_vmtk")
    assert source.index("active_collar_original_side_distance_qc") < cap
    assert source.index("active_collar_cross_section_fidelity_qc") < cap
    assert source.index("cut_seam_topology_qc") < cap
    assert source.index("cross_seam_intersection_qc") < cap
    assert source.index("extension_geometry_qc") < cap


def test_boundary_alignment_does_not_relax_0999_threshold():
    profile = {
        "boundaries": [
            {"boundary_plane_normal_abs_dot_expected_outward": 0.999},
            {"boundary_plane_normal_abs_dot_expected_outward": 1.0},
        ]
    }
    assert boundary_plane_alignment_pass(profile)
    profile["boundaries"][0][
        "boundary_plane_normal_abs_dot_expected_outward"
    ] = 0.998999
    assert not boundary_plane_alignment_pass(profile)


def test_simple_closed_loop_detector_rejects_open_chain():
    assert edges_form_simple_closed_loop(
        np.asarray(((0, 1), (1, 2), (2, 3), (3, 0)), dtype=np.int64)
    )
    assert not edges_form_simple_closed_loop(
        np.asarray(((0, 1), (1, 2), (2, 3)), dtype=np.int64)
    )


def test_symmetric_mesh_size_mismatch_is_bidirectional():
    assert symmetric_mesh_size_mismatch(1.0) == 1.0
    assert symmetric_mesh_size_mismatch(2.0) == symmetric_mesh_size_mismatch(0.5)


def test_interface_diagnostic_normalization_uses_canonical_schema():
    report = normalize_interface_diagnostic(
        {
            "boundaries": [
                {
                    "port_id": _boundary().port_id,
                    "interface_edge_count": 8,
                    "normal_jump_P50_deg": 1.0,
                    "normal_jump_P95_deg": 2.0,
                    "normal_jump_P99_deg": 3.0,
                    "normal_jump_max_deg": 4.0,
                }
            ]
        }
    )
    row = report["boundaries"][0]
    assert row["normal_jump_P95_deg"] == 2.0
    assert row["interface_edge_count"] == 8
