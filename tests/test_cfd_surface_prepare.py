"""Focused tests for local-only CFD boundary surgery and BC correction."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
import trimesh

from utils.cfd_surface_prepare.config import (
    LocalCutConfig,
    SurfaceQCConfig,
    load_surface_prepare_config,
)
from utils.cfd_surface_prepare.export import export_geometry, meter_scale_qc
from utils.cfd_surface_prepare.extension import extrude_and_cap
from utils.cfd_surface_prepare.io import (
    BoundaryInput,
    SurfacePrepareError,
    sha256_file,
)
from utils.cfd_surface_prepare.local_cut import local_plane_cut
from utils.cfd_surface_prepare.pressure_correction import (
    CORRECTION_ROLE,
    build_extended_boundary_conditions,
    calculate_pressure_corrections,
)
from utils.cfd_surface_prepare.qc import (
    boundary_geometry_qc,
    core_surface_preservation_qc,
    extension_collision_qc,
    surface_topology_qc,
)
from utils.cfd_surface_prepare.types import BoundarySurfaceResult, TaggedSurface


@pytest.fixture
def boundary() -> BoundaryInput:
    return BoundaryInput(
        index=0,
        port_id="synthetic__cut_000",
        boundary_origin="CUT_PORT",
        role="ASSUMED_INLET",
        global_node_id=None,
        global_edge_id=10,
        center_um=np.asarray((0.0, 0.0, 1.7)),
        source_radius_um=1.0,
        pressure_original_pa=20.0,
        expected_flow_m3_s=1.0e-15,
        simulation_tangent=np.asarray((0.0, 0.0, 1.0)),
        outward_normal=np.asarray((0.0, 0.0, 1.0)),
        extension_length_um=5.0,
        extension_end_um=np.asarray((0.0, 0.0, 6.7)),
    )


@pytest.fixture
def qc_config() -> SurfaceQCConfig:
    return SurfaceQCConfig(
        require_single_component=True,
        require_watertight=True,
        require_zero_boundary_edges=True,
        require_zero_nonmanifold_edges=True,
        require_zero_self_intersections=True,
        require_zero_degenerate_triangles=True,
        expected_boundary_count=1,
        maximum_cap_planarity_error_um=1.0e-5,
        minimum_normal_dot=0.999,
        maximum_extension_length_error_um=1.0e-4,
        maximum_core_surface_distance_um=1.0e-4,
        maximum_core_surface_p95_distance_um=1.0e-5,
        maximum_equivalent_radius_relative_error=0.08,
    )


@pytest.fixture
def local_config() -> LocalCutConfig:
    return LocalCutConfig(3.0, 2.5, 3.0)


@pytest.fixture
def prepared(boundary: BoundaryInput):
    original = trimesh.creation.capsule(radius=1.0, height=4.0, count=[24, 24])
    cut, loop, cut_report = local_plane_cut(
        TaggedSurface.from_mesh(original),
        boundary,
        radial_factor=3.0,
        axial_back_factor=2.5,
        axial_forward_factor=3.0,
    )
    final, result = extrude_and_cap(cut, loop, boundary)
    final.compact()
    return original, final, loop, cut_report, result


def _result(boundary: BoundaryInput, area_um2: float = np.pi) -> BoundarySurfaceResult:
    return BoundarySurfaceResult(
        boundary_index=boundary.index,
        port_id=boundary.port_id,
        boundary_origin=boundary.boundary_origin,
        role=boundary.role,
        source_radius_um=boundary.source_radius_um,
        extension_length_um=boundary.extension_length_um,
        actual_cap_area_um2=area_um2,
        equivalent_radius_um=float(np.sqrt(area_um2 / np.pi)),
        cap_planarity_error_um=0.0,
        minimum_cap_normal_dot=1.0,
        extension_length_error_um=0.0,
        extension_axis_dot=1.0,
        cut_loop_vertex_count=24,
        cap_face_indices=np.arange(22),
        side_face_indices=np.arange(48),
    )


def test_original_surface_sha_is_unchanged(tmp_path: Path):
    path = tmp_path / "immutable.stl"
    path.write_bytes(b"immutable")
    before = sha256_file(path)
    after = sha256_file(path)
    assert before == after


def test_rounded_cap_can_be_cut_flat(prepared):
    _, _, loop, _, _ = prepared
    assert loop.plane_residual_um <= 1.0e-12
    assert loop.area_um2 > 0


def test_local_cut_gets_exactly_one_loop(prepared):
    _, _, _, report, _ = prepared
    assert report["cut_loop_count"] == 1
    assert report["candidate_surface_component_count"] == 1


def test_actual_loop_extrusion_length_is_exact(prepared):
    _, _, _, _, result = prepared
    assert result.extension_length_error_um <= 1.0e-12


def test_distal_cap_is_planar(prepared):
    _, _, _, _, result = prepared
    assert result.cap_planarity_error_um <= 1.0e-12


def test_distal_cap_normal_is_outward(prepared):
    _, _, _, _, result = prepared
    assert result.minimum_cap_normal_dot >= 0.999999
    assert result.extension_axis_dot >= 0.999999


def test_final_synthetic_surface_is_watertight(prepared, qc_config):
    _, final, _, _, _ = prepared
    report, _ = surface_topology_qc(final, qc_config)
    assert report["watertight"]
    assert report["boundary_edge_count"] == 0


def test_final_synthetic_surface_has_zero_nonmanifold_edges(prepared, qc_config):
    _, final, _, _, _ = prepared
    report, _ = surface_topology_qc(final, qc_config)
    assert report["nonmanifold_edge_count"] == 0


def test_core_outside_local_region_is_unchanged(
    prepared, boundary, local_config, qc_config
):
    original, final, _, _, _ = prepared
    report = core_surface_preservation_qc(
        original, final.mesh(), [boundary], local_config, qc_config
    )
    assert report["status"] == "PASS"
    assert report["max_core_surface_distance_um"] <= 1.0e-12


def test_multiple_nearby_surfaces_are_rejected(boundary):
    first = trimesh.creation.capsule(radius=1.0, height=4.0, count=[16, 16])
    second = first.copy()
    second.apply_translation((2.2, 0.0, 0.0))
    combined = trimesh.util.concatenate((first, second))
    with pytest.raises(SurfacePrepareError, match="LOCAL_PORT_CUT_AMBIGUOUS"):
        local_plane_cut(
            TaggedSurface.from_mesh(combined),
            boundary,
            radial_factor=3.0,
            axial_back_factor=2.5,
            axial_forward_factor=3.0,
        )


def test_extension_collision_is_rejected(qc_config):
    first = trimesh.creation.box(extents=(2.0, 2.0, 2.0))
    second = trimesh.creation.box(extents=(2.0, 2.0, 2.0))
    second.apply_translation((1.0, 0.25, 0.25))
    surface = TaggedSurface.from_mesh(trimesh.util.concatenate((first, second)))
    surface.face_kind[len(first.faces) :] = 1
    _, intersections = surface_topology_qc(surface, qc_config)
    report = extension_collision_qc(surface, intersections)
    assert report["status"] == "FAIL"
    assert report["extension_collision_count"] > 0


def test_inlet_flow_rate_is_not_changed(boundary):
    rows = calculate_pressure_corrections(
        [boundary],
        [_result(boundary)],
        dynamic_viscosity_pa_s=0.00345312,
        allow_negative_gauge_pressure=True,
    )
    assert rows[0]["Q_solver_m3_s"] == boundary.expected_flow_m3_s
    assert rows[0]["pressure_correction_applied"] is False


def test_outlet_solver_pressure_is_lower(boundary):
    outlet = replace(boundary, role="ASSUMED_OUTLET")
    row = calculate_pressure_corrections(
        [outlet],
        [_result(outlet)],
        dynamic_viscosity_pa_s=0.00345312,
        allow_negative_gauge_pressure=True,
    )[0]
    assert row["P_solver_boundary_pa"] < outlet.pressure_original_pa
    assert row["extension_pressure_correction_role"] == CORRECTION_ROLE


def test_true_terminal_negative_gauge_is_allowed(boundary):
    terminal = replace(
        boundary,
        role="ASSUMED_OUTLET",
        boundary_origin="TRUE_TERMINAL",
        pressure_original_pa=0.0,
    )
    row = calculate_pressure_corrections(
        [terminal],
        [_result(terminal)],
        dynamic_viscosity_pa_s=0.00345312,
        allow_negative_gauge_pressure=True,
    )[0]
    assert row["P_solver_boundary_pa"] < 0.0


def test_extended_bc_preserves_source_p_and_q(boundary):
    outlet = replace(boundary, index=1, port_id="out", role="ASSUMED_OUTLET")
    rows = calculate_pressure_corrections(
        [boundary, outlet],
        [_result(boundary), _result(outlet)],
        dynamic_viscosity_pa_s=0.00345312,
        allow_negative_gauge_pressure=True,
    )
    original = {
        "method": "saved",
        "fluid": {},
        "flow_model": {},
        "pressure_reference": "gauge",
        "wall": {"type": "NO_SLIP"},
        "inlet": {
            "port_id": boundary.port_id,
            "flow_rate_m3_s": boundary.expected_flow_m3_s,
            "profile": "PARABOLIC",
        },
        "outlets": [{"port_id": outlet.port_id}],
    }
    extended = build_extended_boundary_conditions(
        original, [boundary, outlet], rows
    )
    assert extended["inlet"]["flow_rate_m3_s"] == boundary.expected_flow_m3_s
    assert extended["outlets"][0]["P_original_1D_pa"] == outlet.pressure_original_pa


def test_boundary_tags_have_expected_count(prepared, boundary, qc_config):
    _, final, _, cut_report, result = prepared
    report = boundary_geometry_qc(
        [boundary], [cut_report], [result], qc_config
    )
    assert report["status"] == "PASS"
    assert set(final.boundary_type) == {0, 1}
    assert set(final.boundary_index[final.boundary_type > 0]) == {0}


def test_meter_stl_scale_is_exact(prepared, boundary, tmp_path: Path):
    _, final, _, _, _ = prepared
    from utils.cfd_surface_prepare.export import create_layout

    layout = create_layout(tmp_path, "synthetic")
    paths = export_geometry(
        final, [boundary], layout, create_meter_copy=True
    )
    report = meter_scale_qc(
        paths["cfd_surface_extended_um_stl"],
        paths["cfd_surface_extended_m_stl"],
    )
    assert report["status"] == "PASS"


def test_strict_yaml_rejects_unknown_key(tmp_path: Path):
    source = (
        Path(__file__).resolve().parents[1]
        / "configs"
        / "cfd_surface_prepare.yaml"
    )
    text = source.read_text(encoding="utf-8") + "\nunknown: true\n"
    path = tmp_path / "bad.yaml"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match="Unknown keys"):
        load_surface_prepare_config(path, project_root=source.parents[1])
