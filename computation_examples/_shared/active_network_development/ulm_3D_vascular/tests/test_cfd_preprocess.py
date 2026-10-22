"""Numerical and mapping tests for global-to-ROI CFD preprocessing."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import yaml
from scipy.integrate import quad

from utils.cfd_preprocess.config import ReadinessConfig, load_cfd_preprocess_config
from utils.cfd_preprocess.io import (
    GeometryReferenceError,
    InputValidationError,
    validate_model_run,
    verify_global_edge_manifest,
)
from utils.cfd_preprocess.one_d_flow import edge_resistance, solve_global_flow
from utils.cfd_preprocess.pipeline import next_stage_for_status
from utils.cfd_preprocess.port_transfer import (
    PortTransferError,
    classify_cut_port,
    transfer_all_boundaries,
    transfer_all_ports,
    transfer_cut_port,
    transfer_true_terminal,
)
from utils.cfd_preprocess.qc import roi_readiness
from utils.sampling.roi_extraction import global_model_from_swc
from utils.sampling.sampling_types import CutPort, ROIRecord


MU = 3.45312e-3


def _model(
    points: list[tuple[float, float, float]],
    radii: list[float],
    parents: list[int],
):
    swc = SimpleNamespace(
        node_ids=np.arange(len(points), dtype=np.int64),
        points_um=np.asarray(points, dtype=float),
        radius_raw_um=np.asarray(radii, dtype=float),
        parent_ids=np.asarray(parents, dtype=np.int64),
    )
    return global_model_from_swc(swc, source_model_id="model", source_mouse_id="mouse")


def _solve(model, flow_rate: float = 2.0e-15, leaf_pressure_pa: float = 0.0):
    return solve_global_flow(
        model,
        mu_pa_s=MU,
        leaf_pressure_pa=leaf_pressure_pa,
        boundary_type="flow_rate",
        mean_velocity_mm_s=None,
        prescribed_flow_m3_s=flow_rate,
        relative_mass_tolerance=1.0e-10,
        relative_node_residual_tolerance=1.0e-10,
        reverse_flow_tolerance_m3_s=1.0e-20,
    )


def _roi(
    positions: list[tuple[float, float, float]],
    radii: list[float],
    edges: list[tuple[int, int]],
    ports: list[CutPort],
    *,
    terminal_local_ids: tuple[int, ...] = (),
    terminal_global_ids: tuple[int, ...] = (),
) -> ROIRecord:
    local_edges = np.asarray(edges, dtype=np.int64)
    return ROIRecord(
        roi_id="roi",
        source_model_id="model",
        source_mouse_id="mouse",
        anchor_id=0,
        anchor_position_um=positions[0],
        bbox_min_um=(0.0, 0.0, 0.0),
        bbox_max_um=(10.0, 10.0, 10.0),
        bbox_center_um=(5.0, 5.0, 5.0),
        bbox_size_um=(10.0, 10.0, 10.0),
        global_node_ids=(),
        global_edge_ids=tuple(port.global_edge_id for port in ports),
        local_node_ids=np.arange(len(positions), dtype=np.int64),
        local_node_global_ids=np.full(len(positions), -1, dtype=np.int64),
        local_node_positions_um=np.asarray(positions, dtype=float),
        local_node_radius_um=np.asarray(radii, dtype=float),
        local_edges=local_edges,
        local_edge_ids=np.arange(len(edges), dtype=np.int64),
        local_edge_global_ids=np.arange(len(edges), dtype=np.int64),
        local_edge_points_um=np.asarray(
            [[positions[a], positions[b]] for a, b in edges], dtype=float
        ),
        local_edge_radius_um=np.asarray([[radii[a], radii[b]] for a, b in edges]),
        true_terminal_local_ids=terminal_local_ids,
        true_terminal_global_ids=terminal_global_ids,
        cut_ports=tuple(ports),
        raw_component_count=1,
        raw_total_vessel_length_um=1.0,
        retained_component_length_um=1.0,
    )


def test_constant_radius_edge_resistance_matches_poiseuille() -> None:
    length = 12.0e-6
    radius = 2.1e-6
    expected = 8.0 * MU * length / (np.pi * radius**4)
    assert edge_resistance(length, radius, radius, MU) == pytest.approx(
        expected, rel=1e-14
    )


def test_linear_radius_resistance_matches_quadrature() -> None:
    length = 17.0e-6
    radius0 = 1.2e-6
    radius1 = 2.4e-6
    integral = quad(
        lambda distance: 1.0 / (radius0 + (radius1 - radius0) * distance / length) ** 4,
        0.0,
        length,
        epsabs=0.0,
        epsrel=1.0e-13,
    )[0]
    expected = 8.0 * MU / np.pi * integral
    assert edge_resistance(length, radius0, radius1, MU) == pytest.approx(
        expected, rel=1e-12
    )


def test_simple_chain_pressure_drop() -> None:
    model = _model([(0, 0, 0), (5, 0, 0), (12, 0, 0)], [2, 2, 2], [-1, 0, 1])
    result = _solve(model)
    expected = result.root_flow_m3_s * np.sum(result.resistances_pa_s_m3)
    assert result.pressure_by_node_id[0] == pytest.approx(expected)
    assert result.pressure_by_node_id[2] == 0.0


def test_symmetric_y_has_equal_flow_split() -> None:
    model = _model([(0, 0, 0), (10, 5, 0), (10, -5, 0)], [2, 1.5, 1.5], [-1, 0, 0])
    result = _solve(model)
    assert result.flows_m3_s[0] == pytest.approx(result.flows_m3_s[1], rel=1e-13)


def test_asymmetric_y_lower_resistance_has_larger_flow() -> None:
    model = _model([(0, 0, 0), (10, 5, 0), (20, -10, 0)], [2, 1.8, 1.0], [-1, 0, 0])
    result = _solve(model)
    lower_resistance = int(np.argmin(result.resistances_pa_s_m3))
    higher_resistance = 1 - lower_resistance
    assert result.flows_m3_s[lower_resistance] > result.flows_m3_s[higher_resistance]


def test_global_mass_conservation() -> None:
    model = _model(
        [(0, 0, 0), (5, 0, 0), (10, 4, 0), (10, -4, 0)],
        [2, 2, 1.3, 1.1],
        [-1, 0, 1, 1],
    )
    result = _solve(model)
    assert result.relative_mass_error <= 1.0e-12
    assert result.maximum_internal_relative_residual <= 1.0e-12


def test_cut_port_alpha_recovery() -> None:
    model = _model([(0, 0, 0), (10, 0, 0)], [2, 1], [-1, 0])
    flow = _solve(model)
    port = CutPort("port", 0, 0, (4.0, 0.0, 0.0), 1.6, "xmin")
    roi = _roi([(4, 0, 0), (8, 0, 0)], [1.6, 1.2], [(0, 1)], [port])
    result = transfer_cut_port(
        roi,
        port,
        model,
        flow,
        mu_pa_s=MU,
        maximum_position_error_um=1.0e-8,
        maximum_radius_relative_error=1.0e-12,
        inlet_length_diameters=5.0,
        outlet_length_diameters=5.0,
    )
    assert result.alpha_on_global_edge == pytest.approx(0.4)


def test_cut_pressure_uses_partial_resistance() -> None:
    model = _model([(0, 0, 0), (10, 0, 0)], [2, 1], [-1, 0])
    flow = _solve(model)
    port = CutPort("port", 0, 0, (4.0, 0.0, 0.0), 1.6, "xmin")
    roi = _roi([(4, 0, 0), (8, 0, 0)], [1.6, 1.2], [(0, 1)], [port])
    result = transfer_cut_port(
        roi,
        port,
        model,
        flow,
        mu_pa_s=MU,
        maximum_position_error_um=1.0e-8,
        maximum_radius_relative_error=1.0e-12,
        inlet_length_diameters=5.0,
        outlet_length_diameters=5.0,
    )
    partial = edge_resistance(4.0e-6, 2.0e-6, 1.6e-6, MU)
    expected = flow.pressure_by_node_id[0] - flow.flows_m3_s[0] * partial
    assert result.pressure_pa == pytest.approx(expected, rel=1e-13)


def test_inlet_classification_uses_local_edge_orientation() -> None:
    port = CutPort("in", 0, 0, (0, 0, 0), 1, "xmin")
    roi = _roi([(0, 0, 0), (1, 0, 0)], [1, 1], [(0, 1)], [port])
    assert classify_cut_port(roi, port) == "ASSUMED_INLET"


def test_outlet_classification_uses_local_edge_orientation() -> None:
    port = CutPort("out", 1, 0, (1, 0, 0), 1, "xmax")
    roi = _roi([(0, 0, 0), (1, 0, 0)], [1, 1], [(0, 1)], [port])
    assert classify_cut_port(roi, port) == "ASSUMED_OUTLET"


def test_port_mass_conservation_for_complete_y_cut() -> None:
    model = _model(
        [(0, 0, 0), (10, 0, 0), (20, 5, 0), (20, -5, 0)],
        [2, 1.8, 1.2, 1.2],
        [-1, 0, 1, 1],
    )
    flow = _solve(model)
    ports = [
        CutPort("in", 0, 0, (5, 0, 0), 1.9, "xmin"),
        CutPort("out1", 2, 1, (15, 2.5, 0), 1.5, "xmax"),
        CutPort("out2", 3, 2, (15, -2.5, 0), 1.5, "xmax"),
    ]
    roi = _roi(
        [(5, 0, 0), (10, 0, 0), (15, 2.5, 0), (15, -2.5, 0)],
        [1.9, 1.8, 1.5, 1.5],
        [(0, 1), (1, 2), (1, 3)],
        ports,
    )
    _, error = transfer_all_ports(
        roi,
        model,
        flow,
        mu_pa_s=MU,
        maximum_position_error_um=1.0e-8,
        maximum_radius_relative_error=1.0e-12,
        inlet_length_diameters=5.0,
        outlet_length_diameters=5.0,
    )
    assert error <= 1.0e-12


def _terminal_transfer_kwargs() -> dict[str, float]:
    return {
        "maximum_position_error_um": 1.0e-8,
        "maximum_radius_relative_error": 1.0e-12,
        "inlet_length_diameters": 5.0,
        "outlet_length_diameters": 5.0,
    }


def _chain_with_terminal():
    model = _model([(0, 0, 0), (10, 0, 0), (20, 0, 0)], [2, 2, 2], [-1, 0, 1])
    flow = _solve(model)
    inlet = CutPort("in", 0, 0, (5, 0, 0), 2, "xmin")
    roi = _roi(
        [(5, 0, 0), (10, 0, 0), (20, 0, 0)],
        [2, 2, 2],
        [(0, 1), (1, 2)],
        [inlet],
        terminal_local_ids=(2,),
        terminal_global_ids=(2,),
    )
    return model, flow, roi


def test_chain_true_terminal_is_assumed_outlet_and_balances_mass() -> None:
    model, flow, roi = _chain_with_terminal()
    boundaries, error = transfer_all_boundaries(
        roi, model, flow, mu_pa_s=MU, **_terminal_transfer_kwargs()
    )
    assert [item.role for item in boundaries] == ["ASSUMED_INLET", "ASSUMED_OUTLET"]
    assert boundaries[1].boundary_origin == "TRUE_TERMINAL"
    assert error <= 1.0e-12


def test_y_cut_outlet_plus_terminal_balances_mass() -> None:
    model = _model(
        [(0, 0, 0), (10, 0, 0), (20, 5, 0), (20, -5, 0)],
        [2, 2, 2, 2],
        [-1, 0, 1, 1],
    )
    flow = _solve(model)
    ports = [
        CutPort("in", 0, 0, (5, 0, 0), 2, "xmin"),
        CutPort("cut_out", 2, 1, (15, 2.5, 0), 2, "xmax"),
    ]
    roi = _roi(
        [(5, 0, 0), (10, 0, 0), (15, 2.5, 0), (20, -5, 0)],
        [2, 2, 2, 2],
        [(0, 1), (1, 2), (1, 3)],
        ports,
        terminal_local_ids=(3,),
        terminal_global_ids=(3,),
    )
    boundaries, error = transfer_all_boundaries(
        roi, model, flow, mu_pa_s=MU, **_terminal_transfer_kwargs()
    )
    inlet_flow = sum(
        item.role_flow_m3_s for item in boundaries if item.role == "ASSUMED_INLET"
    )
    cut_out = next(item for item in boundaries if item.port_id == "cut_out")
    terminal = next(
        item for item in boundaries if item.boundary_origin == "TRUE_TERMINAL"
    )
    assert inlet_flow == pytest.approx(cut_out.role_flow_m3_s + terminal.role_flow_m3_s)
    assert error <= 1.0e-12


def test_true_terminal_pressure_and_flow_come_from_global_solution() -> None:
    model, _, roi = _chain_with_terminal()
    flow = _solve(model, leaf_pressure_pa=17.0)
    terminal = transfer_true_terminal(
        roi,
        model,
        flow,
        terminal_index=0,
        local_node_id=2,
        global_node_id=2,
        **_terminal_transfer_kwargs(),
    )
    assert terminal.pressure_pa == flow.pressure_by_node_id[2] == 17.0
    assert terminal.role_flow_m3_s == flow.flow_by_edge_id[1]
    assert terminal.global_edge_id == 1
    assert terminal.alpha_on_global_edge == 1.0


def test_true_terminal_geometry_points_from_roi_to_terminal() -> None:
    model, flow, roi = _chain_with_terminal()
    terminal = transfer_true_terminal(
        roi,
        model,
        flow,
        terminal_index=0,
        local_node_id=2,
        global_node_id=2,
        **_terminal_transfer_kwargs(),
    )
    np.testing.assert_allclose(terminal.geometry.simulation_tangent, [1, 0, 0])
    np.testing.assert_allclose(terminal.geometry.outward_normal, [1, 0, 0])
    assert np.linalg.norm(terminal.geometry.outward_normal) == pytest.approx(1.0)
    assert terminal.geometry.extension_length_um == 20.0


def test_true_terminal_wrong_global_node_is_rejected() -> None:
    model, flow, roi = _chain_with_terminal()
    with pytest.raises(PortTransferError, match="TRUE_TERMINAL_GLOBAL_MAPPING_INVALID"):
        transfer_true_terminal(
            roi,
            model,
            flow,
            terminal_index=0,
            local_node_id=2,
            global_node_id=99,
            **_terminal_transfer_kwargs(),
        )


def test_true_terminal_nonleaf_global_node_is_rejected() -> None:
    model, flow, roi = _chain_with_terminal()
    with pytest.raises(PortTransferError, match="TRUE_TERMINAL_GLOBAL_MAPPING_INVALID"):
        transfer_true_terminal(
            roi,
            model,
            flow,
            terminal_index=0,
            local_node_id=1,
            global_node_id=1,
            **_terminal_transfer_kwargs(),
        )


def test_true_terminal_position_mismatch_is_rejected() -> None:
    model, flow, roi = _chain_with_terminal()
    roi.local_node_positions_um[2, 0] += 0.1
    with pytest.raises(PortTransferError, match="TRUE_TERMINAL_POSITION_MISMATCH"):
        transfer_true_terminal(
            roi,
            model,
            flow,
            terminal_index=0,
            local_node_id=2,
            global_node_id=2,
            **_terminal_transfer_kwargs(),
        )


def test_true_terminal_radius_mismatch_is_rejected() -> None:
    model, flow, roi = _chain_with_terminal()
    roi.local_node_radius_um[2] += 0.1
    with pytest.raises(PortTransferError, match="TRUE_TERMINAL_RADIUS_MISMATCH"):
        transfer_true_terminal(
            roi,
            model,
            flow,
            terminal_index=0,
            local_node_id=2,
            global_node_id=2,
            **_terminal_transfer_kwargs(),
        )


def test_true_terminal_local_degree_must_be_one() -> None:
    model, flow, _ = _chain_with_terminal()
    roi = _roi(
        [(5, 0, 0), (20, 0, 0), (15, 0, 0)],
        [2, 2, 2],
        [(0, 1), (1, 2)],
        [],
        terminal_local_ids=(1,),
        terminal_global_ids=(2,),
    )
    with pytest.raises(PortTransferError, match="TRUE_TERMINAL_LOCAL_TOPOLOGY_INVALID"):
        transfer_true_terminal(
            roi,
            model,
            flow,
            terminal_index=0,
            local_node_id=1,
            global_node_id=2,
            **_terminal_transfer_kwargs(),
        )


def _readiness_config() -> ReadinessConfig:
    return ReadinessConfig(
        minimum_boundary_count=2,
        true_terminal_policy="assumed_outlet",
        required_assumed_inlet_count=1,
        minimum_assumed_outlet_count=1,
        require_connected_roi=True,
        require_cycle_rank_zero=True,
    )


@pytest.mark.parametrize(
    ("field", "value", "check", "expected"),
    [
        ("position_error_um", 0.5e-5, "all_boundary_positions_match", True),
        ("position_error_um", 2.0e-5, "all_boundary_positions_match", False),
        ("radius_relative_error", 0.5e-6, "all_boundary_radii_match", True),
        ("radius_relative_error", 2.0e-6, "all_boundary_radii_match", False),
    ],
)
def test_readiness_uses_real_mapping_tolerances(
    field: str, value: float, check: str, expected: bool
) -> None:
    model, flow, roi = _chain_with_terminal()
    boundaries, error = transfer_all_boundaries(
        roi, model, flow, mu_pa_s=MU, **_terminal_transfer_kwargs()
    )
    boundaries[0] = replace(boundaries[0], **{field: value})
    result = roi_readiness(
        roi,
        boundaries,
        error,
        _readiness_config(),
        boundary_mass_tolerance=1.0e-8,
        maximum_position_error_um=1.0e-5,
        maximum_radius_relative_error=1.0e-6,
    )
    assert result["checks"][check] is expected


def test_next_stage_depends_on_final_status() -> None:
    assert (
        next_stage_for_status("CFD_PREPROCESS_BASELINE_PASS")
        == "CFD_SURFACE_VOLUME_MESH_PREPARATION"
    )
    assert (
        next_stage_for_status("CFD_ROI_NOT_READY")
        != "CFD_SURFACE_VOLUME_MESH_PREPARATION"
    )
    assert next_stage_for_status("GLOBAL_1D_FLOW_FAILED") == "REVIEW_GLOBAL_1D_FLOW"
    assert next_stage_for_status("GLOBAL_TO_ROI_TRANSFER_FAILED") == (
        "REVIEW_GLOBAL_TO_ROI_BOUNDARY_TRANSFER"
    )


def test_global_edge_manifest_mismatch_is_rejected(tmp_path: Path) -> None:
    model = _model([(0, 0, 0), (1, 0, 0)], [1, 1], [-1, 0])
    manifest = tmp_path / "global_edges.csv"
    manifest.write_text(
        "source_model_id,global_edge_id,upstream_node_id,downstream_node_id\n"
        "model,0,99,1\n",
        encoding="utf-8",
    )
    with pytest.raises(InputValidationError, match="GLOBAL_EDGE_MAPPING_MISMATCH"):
        verify_global_edge_manifest(model, manifest)


def test_model_run_roi_mismatch_is_rejected(tmp_path: Path) -> None:
    qc = tmp_path / "qc"
    qc.mkdir()
    (qc / "run_summary.json").write_text(
        json.dumps({"status": "PASS", "roi_id": "wrong", "radius_scale": 0.91}),
        encoding="utf-8",
    )
    (qc / "surface_qc.json").write_text(
        json.dumps({"status": "PASS"}), encoding="utf-8"
    )
    (qc / "radius_fidelity.json").write_text(
        json.dumps({"p95_absolute_relative_error": 0.04}), encoding="utf-8"
    )
    with pytest.raises(GeometryReferenceError, match="ROI mismatch"):
        validate_model_run(tmp_path, roi_id="expected")


def test_yaml_unknown_key_is_rejected(tmp_path: Path) -> None:
    project = Path(__file__).resolve().parents[1]
    payload = yaml.safe_load((project / "configs" / "cfd_preprocess.yaml").read_text())
    payload["fluid"]["silent_fallback"] = 1
    path = tmp_path / "invalid.yaml"
    path.write_text(yaml.safe_dump(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="Unknown keys in fluid"):
        load_cfd_preprocess_config(path, project_root=project)


def test_yaml_requires_exactly_one_inlet_mode(tmp_path: Path) -> None:
    project = Path(__file__).resolve().parents[1]
    payload = yaml.safe_load((project / "configs" / "cfd_preprocess.yaml").read_text())
    payload["global_1d"]["inlet"]["flow_rate_m3_s"] = 1.0e-15
    path = tmp_path / "invalid_inlet.yaml"
    path.write_text(yaml.safe_dump(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="requires only mean_velocity_mm_s"):
        load_cfd_preprocess_config(path, project_root=project)
