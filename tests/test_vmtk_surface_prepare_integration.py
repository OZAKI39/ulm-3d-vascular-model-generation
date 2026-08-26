"""Project-side integration tests for the official VMTK exchange boundary."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pyvista as pv
import pytest

from utils.cfd_surface_prepare.config import (
    VmtkConfig,
    load_surface_prepare_config,
)
from utils.cfd_surface_prepare.io import BoundaryInput, SurfacePrepareError
from utils.cfd_surface_prepare.io import (
    load_original_surface,
    load_surface_inputs,
    sha256_file,
)
from utils.cfd_surface_prepare.local_cut import local_plane_cut
from utils.cfd_surface_prepare.types import TaggedSurface
from utils.cfd_surface_prepare.vmtk_adapter import build_centerline_adapter
from utils.cfd_surface_prepare.vmtk_qc import (
    extension_vector_measurements,
    open_profile_qc,
    polydata_mesh,
    tag_and_export_final_surface,
    topology_qc,
)
from utils.cfd_surface_prepare.vmtk_pipeline import (
    WARNING_STATUS,
    boundary_plane_alignment_pass,
    final_candidate_status,
    raw_geometry_hard_gate_pass,
    should_generate_manual_review_figures,
    should_promote_raw_candidate,
)
from utils.cfd_surface_prepare.vmtk_runner import (
    exchange_paths,
    official_source_provenance,
    parameter_mapping,
    promote_official_vmtk,
    run_official_vmtk,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = PROJECT_ROOT / "configs" / "cfd_surface_prepare.yaml"
PMP_PYTHON = Path("D:/anaconda3/envs/pmp/python.exe")
VMTK_RUNTIME_PREFIX = Path("D:/anaconda3/envs/vmtk-env")


def _boundary(index: int, z: float, outward: float, role: str) -> BoundaryInput:
    center = np.asarray((0.0, 0.0, z))
    normal = np.asarray((0.0, 0.0, outward))
    return BoundaryInput(
        index=index,
        port_id=f"synthetic__{index}",
        boundary_origin="CUT_PORT",
        role=role,
        global_node_id=None,
        global_edge_id=index,
        center_um=center,
        source_radius_um=1.0,
        pressure_original_pa=100.0,
        expected_flow_m3_s=1.0e-15,
        simulation_tangent=np.asarray((0.0, 0.0, 1.0)),
        outward_normal=normal,
        extension_length_um=10.0,
        extension_end_um=center + 10.0 * normal,
    )


def _vmtk_config(extension_mode: str) -> VmtkConfig:
    return VmtkConfig(
        environment_python=PMP_PYTHON,
        runtime_prefix=VMTK_RUNTIME_PREFIX,
        official_repository=PROJECT_ROOT.parent / "external" / "vmtk",
        interpolation_mode="thinplatespline",
        preserve_cross_section_shape=False,
        extension_mode=extension_mode,
        sigma=1.0,
        transition_ratio=0.5,
        adaptive_extension_length=True,
        extension_ratio=10.0,
        adaptive_extension_radius=True,
        adaptive_boundary_points=True,
        remesh_after_extension=True,
    )


@pytest.fixture(scope="module")
def formal_config():
    return load_surface_prepare_config(CONFIG_PATH, project_root=PROJECT_ROOT)


@pytest.fixture(scope="module")
def synthetic_vmtk(tmp_path_factory):
    if not PMP_PYTHON.is_file() or not VMTK_RUNTIME_PREFIX.is_dir():
        pytest.skip("pmp or the pinned VMTK runtime is unavailable")
    root = tmp_path_factory.mktemp("official_vmtk")
    input_dir = root / "input"
    vmtk_dir = root / "vmtk"
    geometry_dir = root / "geometry"
    for directory in (input_dir, vmtk_dir, geometry_dir):
        directory.mkdir()
    paths = exchange_paths(
        input_directory=input_dir,
        vmtk_directory=vmtk_dir,
        geometry_directory=geometry_dir,
        extension_mode="boundarynormal",
    )
    count = 32
    layers = 7
    angle = np.linspace(0.0, 2.0 * np.pi, count, endpoint=False)
    points = np.vstack(
        [
            np.column_stack(
                (1.15 * np.cos(angle), 0.85 * np.sin(angle), np.full(count, z))
            )
            for z in np.linspace(0.0, 5.0, layers)
        ]
    )
    faces: list[list[int]] = []
    for layer in range(layers - 1):
        first = layer * count
        second = (layer + 1) * count
        for point in range(count):
            following = (point + 1) % count
            faces.extend(
                (
                    [first + point, first + following, second + following],
                    [first + point, second + following, second + point],
                )
            )
    vtk_faces = np.column_stack(
        (np.full(len(faces), 3, dtype=np.int64), np.asarray(faces, dtype=np.int64))
    ).ravel()
    pv.PolyData(points, vtk_faces).save(paths.open_surface_vtp, binary=True)
    config = _vmtk_config("boundarynormal")
    before = hashlib.sha256(paths.open_surface_vtp.read_bytes()).hexdigest()
    invocation = run_official_vmtk(
        config=config,
        paths=paths,
        tool_script=PROJECT_ROOT / "tools" / "run_vmtk_flowextension.py",
    )
    promotion = promote_official_vmtk(
        config=config,
        paths=paths,
        tool_script=PROJECT_ROOT / "tools" / "run_vmtk_flowextension.py",
        target_edge_length_um=0.45,
    )
    after = hashlib.sha256(paths.open_surface_vtp.read_bytes()).hexdigest()
    return root, paths, invocation, promotion, before, after


@pytest.fixture(scope="module")
def synthetic_centerline_vmtk(tmp_path_factory, synthetic_vmtk):
    root = tmp_path_factory.mktemp("official_vmtk_centerline")
    input_dir = root / "input"
    vmtk_dir = root / "vmtk"
    geometry_dir = root / "geometry"
    for directory in (input_dir, vmtk_dir, geometry_dir):
        directory.mkdir()
    paths = exchange_paths(
        input_directory=input_dir,
        vmtk_directory=vmtk_dir,
        geometry_directory=geometry_dir,
        extension_mode="centerlinedirection",
    )
    source_paths = synthetic_vmtk[1]
    pv.read(source_paths.open_surface_vtp).save(paths.open_surface_vtp, binary=True)
    centerline = pv.PolyData()
    centerline.points = np.asarray(((0.0, 0.0, -2.0), (0.0, 0.0, 7.0)))
    centerline.lines = np.asarray((2, 0, 1), dtype=np.int64)
    centerline.save(paths.centerlines_vtp, binary=True)
    config = _vmtk_config("centerlinedirection")
    invocation = run_official_vmtk(
        config=config,
        paths=paths,
        tool_script=PROJECT_ROOT / "tools" / "run_vmtk_flowextension.py",
    )
    return paths, invocation


def test_formal_backend_is_official_vmtk(formal_config):
    assert formal_config.backend.method == "vmtk_flowextensions"


def test_formal_yaml_selects_boundarynormal(formal_config):
    assert formal_config.vmtk.extension_mode == "boundarynormal"


def test_formal_interpolation_is_only_thin_plate_spline(formal_config):
    assert formal_config.vmtk.interpolation_mode == "thinplatespline"
    assert formal_config.vmtk.transition_ratio == 0.5


def test_formal_configuration_has_no_shape_preservation_or_fallback(formal_config):
    mapping = parameter_mapping(formal_config.vmtk)
    assert mapping["preserve_cross_section_shape"] is False
    assert mapping["automatic_fallback"] is False
    assert mapping["custom_tps_implementation"] is False


def test_five_d_maps_to_verified_adaptive_ratio_ten(formal_config):
    mapping = parameter_mapping(formal_config.vmtk)
    assert mapping["project_extension_definition"] == "5D"
    assert mapping["extension_ratio"] == 10.0
    assert "meanRadius" in mapping["verified_vmtk_length_definition"]


def test_single_variable_scientific_parameters_are_unchanged(formal_config):
    mapping = parameter_mapping(formal_config.vmtk)
    assert mapping["interpolation_mode"] == "thinplatespline"
    assert mapping["sigma"] == 1.0
    assert mapping["transition_ratio"] == 0.5
    assert mapping["extension_ratio"] == 10.0
    assert mapping["adaptive_extension_length"] is True
    assert mapping["adaptive_extension_radius"] is True
    assert mapping["adaptive_boundary_points"] is True
    assert mapping["preserve_cross_section_shape"] is False
    assert mapping["parameter_sweep"] is False


def test_invalid_yaml_extension_mode_fails(tmp_path: Path):
    invalid = tmp_path / "invalid.yaml"
    invalid.write_text(
        CONFIG_PATH.read_text(encoding="utf-8").replace(
            "extension_mode: boundarynormal", "extension_mode: diagonal"
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="vmtk.extension_mode"):
        load_surface_prepare_config(invalid, project_root=PROJECT_ROOT)


def test_official_source_provenance_records_release_commit(formal_config):
    provenance = official_source_provenance(formal_config.vmtk)
    assert provenance["runtime_release_tag"] == "v1.5.0"
    assert provenance["runtime_release_tag_commit"] == "30d5d7cb8e607d153c208a9d7d39c9feb7985476"


def test_centerline_adapter_uses_saved_real_graph_nodes(tmp_path: Path):
    run = tmp_path / "run" / "global_1d"
    run.mkdir(parents=True)
    (run / "nodes.csv").write_text(
        "node_id,parent_id,x_um,y_um,z_um,radius_um\n"
        "0,-1,0,0,-2,1\n1,0,0,0,0,1\n2,1,0,0,2,1\n3,2,0,0,4,1\n",
        encoding="utf-8",
    )
    (run / "edges.csv").write_text(
        "global_edge_id,parent_node_id,child_node_id\n0,1,2\n",
        encoding="utf-8",
    )
    output = tmp_path / "centerlines.vtp"
    result = build_centerline_adapter(
        tmp_path / "run", [_boundary(0, 1.0, 1.0, "ASSUMED_INLET")], output
    )
    assert result.records[0]["invented_centerline_points"] == 0
    assert set(result.records[0]["source_node_ids"]).issubset({0, 1, 2, 3})
    assert output.is_file()


def test_actual_roi_preflight_has_four_profiles_and_real_centerlines(
    formal_config, tmp_path: Path
):
    inputs = load_surface_inputs(
        formal_config.paths.cfd_preprocess_run,
        expected_boundary_count=4,
    )
    before = {
        inputs.original_surface_um_stl: sha256_file(inputs.original_surface_um_stl),
        inputs.original_surface_um_vtp: sha256_file(inputs.original_surface_um_vtp),
    }
    original = load_original_surface(inputs.original_surface_um_vtp)
    surface = TaggedSurface.from_mesh(original)
    for boundary in inputs.boundaries:
        surface, _, _ = local_plane_cut(
            surface,
            boundary,
            radial_factor=formal_config.local_cut.local_radial_radius_factor,
            axial_back_factor=formal_config.local_cut.local_axial_back_radius_factor,
            axial_forward_factor=formal_config.local_cut.local_axial_forward_radius_factor,
        )
    surface.compact()
    faces = np.column_stack(
        (np.full(len(surface.faces), 3, dtype=np.int64), surface.faces)
    ).ravel()
    open_path = tmp_path / "open_surface_um.vtp"
    pv.PolyData(surface.vertices, faces).save(open_path, binary=True)
    _, open_mesh = polydata_mesh(open_path)
    report, _ = open_profile_qc(open_mesh, inputs.boundaries)
    adapter = build_centerline_adapter(
        inputs.preprocess_run,
        inputs.boundaries,
        tmp_path / "centerlines_um.vtp",
    )
    assert report["status"] == "PASS"
    assert report["profile_count"] == 4
    assert boundary_plane_alignment_pass(report)
    assert all(
        row["boundary_plane_normal_abs_dot_expected_outward"] >= 0.999
        for row in report["boundaries"]
    )
    assert all(row["invented_centerline_points"] == 0 for row in adapter.records)
    assert all(before[path] == sha256_file(path) for path in before)


def test_runner_uses_pmp_environment(synthetic_vmtk):
    _, _, invocation, _, _, _ = synthetic_vmtk
    assert Path(invocation.command[0]).resolve() == PMP_PYTHON.resolve()
    assert Path(invocation.runtime["python_executable"]).resolve() == PMP_PYTHON.resolve()
    assert invocation.runtime["vtk_version"] == "9.2.6"


def test_runner_calls_official_flow_filter(synthetic_vmtk):
    _, _, invocation, _, _, _ = synthetic_vmtk
    assert invocation.runtime["official_flow_filter"] == "vtkvmtkPolyDataFlowExtensionsFilter"
    assert invocation.runtime["custom_tps_implementation"] is False


def test_boundarynormal_runner_calls_official_normal_api(synthetic_vmtk):
    _, _, invocation, _, _, _ = synthetic_vmtk
    assert (
        invocation.runtime["official_direction_api"]
        == "SetExtensionModeToUseNormalToBoundary"
    )
    assert invocation.runtime["extension_mode_effective"] == "boundarynormal"


def test_boundarynormal_runner_does_not_use_centerlines(synthetic_vmtk):
    _, paths, invocation, _, _, _ = synthetic_vmtk
    assert invocation.runtime["centerlines_used_for_extension_direction"] is False
    assert "centerlines_vtp" not in invocation.request
    assert not paths.centerlines_vtp.exists()
    assert (
        invocation.runtime["official_direction_api"]
        != "SetExtensionModeToUseCenterlineDirection"
    )


def test_centerlinedirection_runner_path_remains_available(
    synthetic_centerline_vmtk,
):
    paths, invocation = synthetic_centerline_vmtk
    assert paths.centerlines_vtp.is_file()
    assert "centerlines_vtp" in invocation.request
    assert invocation.runtime["centerlines_used_for_extension_direction"] is True
    assert (
        invocation.runtime["official_direction_api"]
        == "SetExtensionModeToUseCenterlineDirection"
    )


def test_invalid_runner_extension_mode_fails(synthetic_vmtk):
    _, paths, _, _, _, _ = synthetic_vmtk
    with pytest.raises(SurfacePrepareError, match="INVALID_VMTK_EXTENSION_MODE"):
        run_official_vmtk(
            config=_vmtk_config("invalid"),
            paths=paths,
            tool_script=PROJECT_ROOT / "tools" / "run_vmtk_flowextension.py",
        )


def test_vmtk_raw_output_exists(synthetic_vmtk):
    _, paths, _, _, _, _ = synthetic_vmtk
    assert paths.raw_vtp.is_file() and paths.raw_stl.is_file()


def test_vmtk_raw_has_detectable_extensions(synthetic_vmtk):
    _, paths, invocation, _, _, _ = synthetic_vmtk
    assert invocation.runtime["raw_cells"] > invocation.runtime["input_cells"]


def test_vmtk_raw_has_two_synthetic_open_profiles(synthetic_vmtk):
    _, paths, _, _, _, _ = synthetic_vmtk
    _, mesh = polydata_mesh(paths.raw_vtp)
    report, _ = topology_qc(
        mesh, expected_open_profile_count=2, allow_degenerate=True
    )
    assert report["status"] == "PASS"


def test_vmtk_remesher_uses_one_global_target_size(synthetic_vmtk):
    _, _, _, promotion, _, _ = synthetic_vmtk
    assert promotion.runtime["remesh"]["target_size_count"] == 1
    assert promotion.runtime["remesh"]["parameter_sweep_count"] == 0


def test_vmtk_simple_cap_ids_are_present(synthetic_vmtk):
    _, _, _, promotion, _, _ = synthetic_vmtk
    assert promotion.runtime["cap"]["method"] == "simple"
    assert len(promotion.runtime["cap_entity_ids"]) == 2


def test_tagged_caps_keep_project_boundary_metadata(synthetic_vmtk):
    root, paths, _, _, _, _ = synthetic_vmtk
    outputs, report = tag_and_export_final_surface(
        paths.capped_vtp,
        (_boundary(0, 0.0, -1.0, "ASSUMED_INLET"), _boundary(1, 5.0, 1.0, "ASSUMED_OUTLET")),
        root / "final",
        root / "boundaries",
    )
    data = pv.read(outputs["tagged_vtp"])
    assert set(np.unique(data.cell_data["boundary_type_code"])) == {0, 1, 2}
    assert report["distal_boundary_count"] == 2


def test_capped_synthetic_surface_is_watertight(synthetic_vmtk):
    _, paths, _, _, _, _ = synthetic_vmtk
    _, mesh = polydata_mesh(paths.capped_vtp)
    report, _ = topology_qc(mesh, expected_open_profile_count=0)
    assert report["status"] == "PASS"


def test_runner_does_not_modify_exchange_input(synthetic_vmtk):
    _, _, _, _, before, after = synthetic_vmtk
    assert before == after


def test_extension_and_promotion_are_separate_operations(synthetic_vmtk):
    _, _, invocation, promotion, _, _ = synthetic_vmtk
    assert invocation.request["operation"] == "extension"
    assert "remeshed_open_vtp" not in invocation.request
    assert "capped_vtp" not in invocation.request
    assert promotion.request["operation"] == "remesh_cap"


def test_boundarynormal_runtime_keeps_tps_parameters(synthetic_vmtk):
    _, _, invocation, _, _, _ = synthetic_vmtk
    parameters = invocation.runtime["parameters"]
    assert invocation.runtime["official_interpolation_api"] == (
        "SetInterpolationModeToThinPlateSpline"
    )
    assert parameters["sigma"] == 1.0
    assert parameters["transition_ratio"] == 0.5
    assert parameters["extension_ratio"] == 10.0
    assert parameters["adaptive_extension_radius"] is True
    assert parameters["adaptive_boundary_points"] is True


def test_direction_qc_uses_signed_dot():
    measurements = extension_vector_measurements(
        np.zeros(3), np.asarray((0.0, 0.0, -10.0)), np.asarray((0.0, 0.0, 1.0))
    )
    assert measurements["extension_direction_dot"] == -1.0
    assert measurements["actual_axial_length_um"] == -10.0


def test_preboundary_plane_qc_uses_absolute_dot(synthetic_vmtk):
    _, paths, _, _, _, _ = synthetic_vmtk
    _, mesh = polydata_mesh(paths.open_surface_vtp)
    report, _ = open_profile_qc(
        mesh,
        (
            _boundary(0, 0.0, -1.0, "ASSUMED_INLET"),
            _boundary(1, 5.0, 1.0, "ASSUMED_OUTLET"),
        ),
    )
    assert boundary_plane_alignment_pass(report)
    assert all(
        0.999 <= row["boundary_plane_normal_abs_dot_expected_outward"] <= 1.0
        for row in report["boundaries"]
    )


def test_interface_warning_is_not_a_raw_geometry_hard_stop():
    passed = {"status": "PASS"}
    assert raw_geometry_hard_gate_pass(passed, passed, passed)
    assert should_promote_raw_candidate(True)
    assert final_candidate_status(
        final_hard_qc_pass=True, interface_status="WARNING"
    ) == WARNING_STATUS


def test_raw_geometry_failure_blocks_promotion():
    assert not should_promote_raw_candidate(False)


def test_interface_warning_still_requires_diagnostic_figures():
    assert should_generate_manual_review_figures(
        raw_hard_qc_pass=True, interface_status="WARNING"
    )


def test_boundarynormal_output_names_are_explicit(synthetic_vmtk):
    _, paths, _, _, _, _ = synthetic_vmtk
    assert paths.raw_vtp.name == "vmtk_boundarynormal_raw_um.vtp"
    assert paths.raw_stl.name == "vmtk_boundarynormal_raw_um.stl"
