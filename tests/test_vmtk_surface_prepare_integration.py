"""Project-side integration tests for the official VMTK exchange boundary."""

from __future__ import annotations

import hashlib
import ast
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pyvista as pv
import pytest
import trimesh

import utils.cfd_surface_prepare.vmtk_pipeline as vmtk_pipeline

from utils.cfd_surface_prepare.config import (
    CoreCollarConfig,
    EntityRemeshConfig,
    MeshQualityConfig,
    VmtkConfig,
    load_surface_prepare_config,
)
from utils.cfd_surface_prepare.guarded_remesh import (
    active_collar_original_side_distance_qc,
    assign_cross_seam_active_entities,
    core_face_adjacency_layers,
    cross_seam_entity_preservation_qc,
    cross_seam_intersection_qc,
    cut_seam_topology_qc,
    edges_form_simple_closed_loop,
    locked_entity_preservation_qc,
    triangle_pair_intersection_diagnosis,
)
from utils.cfd_surface_prepare.io import BoundaryInput, SurfacePrepareError
from utils.cfd_surface_prepare.io import (
    load_original_surface,
    load_surface_inputs,
    read_json,
    sha256_file,
)
from utils.cfd_surface_prepare.local_cut import local_plane_cut
from utils.cfd_surface_prepare.types import TaggedSurface
from utils.cfd_surface_prepare.vmtk_adapter import build_centerline_adapter
from utils.cfd_surface_prepare.vmtk_qc import (
    active_collar_cross_section_fidelity_qc,
    BoundaryLoop,
    assign_remesh_entities,
    core_symmetric_distance_qc,
    extension_geometry_qc,
    extension_mesh_quality_from_surface,
    extension_vector_measurements,
    meter_scale_qc,
    normalize_interface_diagnostic,
    normal_consistency_qc,
    open_profile_qc,
    polydata_mesh,
    raw_core_exact_copy_qc,
    symmetric_mesh_size_mismatch,
    tag_and_export_final_surface,
    topology_qc,
)
from utils.cfd_surface_prepare.vmtk_pipeline import (
    METER_FAILURE,
    SUCCESS_STATUS,
    TOPOLOGY_FAILURE,
    boundary_plane_alignment_pass,
    final_candidate_status,
    raw_geometry_hard_gate_pass,
    resolve_entity_failure_status,
    resume_open_qc_allows_cap,
    run_crossseam_open_resume,
    should_cap_crossseam_open,
    should_generate_manual_review_figures,
    should_promote_raw_candidate,
)
from utils.cfd_surface_prepare.vmtk_runner import (
    build_entity_remesh_request,
    cap_official_vmtk,
    entity_remesh_official_vmtk,
    exchange_paths,
    official_source_provenance,
    parameter_mapping,
    run_official_vmtk,
)
from utils.cfd_surface_prepare.config import LocalCutConfig


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = PROJECT_ROOT / "configs" / "cfd_surface_prepare.yaml"
PMP_PYTHON = Path("D:/anaconda3/envs/pmp/python.exe")
VMTK_RUNTIME_PREFIX = Path("D:/anaconda3/envs/vmtk-env")
GUARDED_ENTITY_RUN = (
    PROJECT_ROOT
    / "outputs"
    / "cfd_surface_prepare"
    / "vmtk_tps_boundarynormal_guarded_entityremesh_anchor003274_20260826_185752"
)
PREVIOUS_ENTITY_RUN = (
    PROJECT_ROOT
    / "outputs"
    / "cfd_surface_prepare"
    / "vmtk_tps_boundarynormal_entityremesh_anchor003274_20260826_180737"
)


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


def _entity_remesh_config() -> EntityRemeshConfig:
    return EntityRemeshConfig(
        enabled=True,
        entity_array_name="RemeshEntityId",
        far_core_entity_id=1,
        active_entity_id=2,
        exclude_entity_ids=(1,),
        core_collar=CoreCollarConfig(
            mode="core_face_adjacency_layers", face_layers=2
        ),
        element_size_mode="edgelength",
        target_edge_length_um=0.25913916380971913,
        preserve_boundary_edges=True,
    )


def _vmtk_config(
    extension_mode: str, *, entity_remesh: bool = False
) -> VmtkConfig:
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
        postprocess_mode=(
            "cross_seam_active_collar_remesh_then_cap"
            if entity_remesh
            else "cap_only"
        ),
        remesh_after_extension=entity_remesh,
        entity_remesh=_entity_remesh_config(),
    )


@pytest.fixture(scope="module")
def formal_config():
    return load_surface_prepare_config(CONFIG_PATH, project_root=PROJECT_ROOT)


@pytest.fixture(scope="module")
def actual_crossseam_assignment(formal_config, tmp_path_factory):
    """Classify the saved RAW surface without running synthetic VMTK."""

    root = tmp_path_factory.mktemp("saved_raw_crossseam_assignment")
    raw_path = root / "saved_raw.vtp"
    pv.read(
        GUARDED_ENTITY_RUN / "geometry" / "vmtk_boundarynormal_raw_um.vtp"
    ).save(raw_path, binary=True)
    inputs = load_surface_inputs(
        formal_config.paths.cfd_preprocess_run, expected_boundary_count=4
    )
    settings = formal_config.vmtk.entity_remesh
    report, distances = assign_cross_seam_active_entities(
        raw_path,
        inputs.boundaries,
        face_layers=settings.core_collar.face_layers,
        entity_array_name=settings.entity_array_name,
        far_core_entity_id=settings.far_core_entity_id,
        active_entity_id=settings.active_entity_id,
    )
    return raw_path, report, distances, inputs


@pytest.fixture(scope="module")
def crossseam_modified_candidate(actual_crossseam_assignment, tmp_path_factory):
    """Create a deterministic extension-only change without invoking VMTK."""

    raw_path, _, _, inputs = actual_crossseam_assignment
    candidate = pv.read(raw_path).triangulate()
    faces = np.asarray(candidate.faces).reshape((-1, 4))[:, 1:]
    regions = np.asarray(candidate.cell_data["SurfaceRegionId"])
    core_vertices = np.unique(faces[regions == 0])
    extension_vertices = np.unique(faces[regions == 1])
    extension_only = np.setdiff1d(extension_vertices, core_vertices)
    assert len(extension_only) > 0
    candidate.points[int(extension_only[0]), 0] += 1.0e-3
    candidate_path = (
        tmp_path_factory.mktemp("crossseam_modified_candidate")
        / "candidate.vtp"
    )
    candidate.save(candidate_path, binary=True)
    reports = cross_seam_entity_preservation_qc(
        raw_path,
        candidate_path,
        inputs.boundaries,
        restore_region_arrays=False,
    )
    return raw_path, candidate_path, reports


@pytest.fixture()
def four_port_remeshed_open(tmp_path: Path):
    """Controlled four-port distal-open surface for post-remesh QC tests."""

    count = 32
    angle = np.linspace(0.0, 2.0 * np.pi, count, endpoint=False)
    points: list[list[float]] = []
    faces: list[list[int]] = []
    boundaries: list[BoundaryInput] = []
    proximal: dict[int, BoundaryLoop] = {}
    polygon_area = 0.5 * count * np.sin(2.0 * np.pi / count)
    for index, x_center in enumerate((0.0, 4.0, 8.0, 12.0)):
        base = len(points)
        bottom = np.column_stack(
            (
                x_center + np.cos(angle),
                np.sin(angle),
                np.zeros(count),
            )
        )
        top = bottom.copy()
        top[:, 2] = 3.0
        points.extend(bottom.tolist())
        points.extend(top.tolist())
        bottom_center = len(points)
        points.append([x_center, 0.0, 0.0])
        for point in range(count):
            following = (point + 1) % count
            faces.extend(
                (
                    [base + point, base + following, base + count + following],
                    [base + point, base + count + following, base + count + point],
                    [bottom_center, base + following, base + point],
                )
            )
        center = np.asarray((x_center, 0.0, 0.0))
        role = "ASSUMED_INLET" if index == 0 else "ASSUMED_OUTLET"
        boundary = BoundaryInput(
            index=index,
            port_id=f"four_port__{index}",
            boundary_origin="CUT_PORT",
            role=role,
            global_node_id=None,
            global_edge_id=index,
            center_um=center,
            source_radius_um=1.0,
            pressure_original_pa=100.0,
            expected_flow_m3_s=1.0e-15,
            simulation_tangent=np.asarray((0.0, 0.0, 1.0)),
            outward_normal=np.asarray((0.0, 0.0, 1.0)),
            extension_length_um=3.0,
            extension_end_um=center + np.asarray((0.0, 0.0, 3.0)),
        )
        boundaries.append(boundary)
        proximal[index] = BoundaryLoop(
            point_ids=np.arange(count, dtype=np.int64),
            points=bottom,
            center_um=center,
            area_um2=float(polygon_area),
            equivalent_radius_um=float(np.sqrt(polygon_area / np.pi)),
            normal=np.asarray((0.0, 0.0, 1.0)),
            planarity_error_um=0.0,
        )
    vtk_faces = np.column_stack(
        (np.full(len(faces), 3, dtype=np.int64), np.asarray(faces, dtype=np.int64))
    ).ravel()
    path = tmp_path / "four_port_entity_remeshed_open.vtp"
    pv.PolyData(np.asarray(points), vtk_faces).save(path, binary=True)
    _, mesh = polydata_mesh(path)
    interface = {
        "boundaries": [
            {
                "port_id": boundary.port_id,
                "interface_edge_count": count,
                "normal_jump_P50_deg": 0.0,
                "normal_jump_P95_deg": 0.0,
                "normal_jump_P99_deg": 0.0,
                "normal_jump_max_deg": 0.0,
            }
            for boundary in boundaries
        ]
    }
    return mesh, tuple(boundaries), proximal, interface


@pytest.fixture()
def frozen_crossseam_tube(tmp_path: Path):
    """One unchanged frozen tube spanning original CORE and extension regions."""

    count = 48
    levels = np.asarray((-1.25, 0.0, 2.0), dtype=float)
    angle = np.linspace(0.0, 2.0 * np.pi, count, endpoint=False)
    points = np.vstack(
        [
            np.column_stack(
                (np.cos(angle), np.sin(angle), np.full(count, level))
            )
            for level in levels
        ]
    ).astype(np.float32)
    faces: list[list[int]] = []
    regions: list[int] = []
    for layer in range(len(levels) - 1):
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
            regions.extend((layer, layer))
    packed = np.column_stack(
        (np.full(len(faces), 3, dtype=np.int64), np.asarray(faces))
    ).ravel()
    data = pv.PolyData(points, packed)
    data.cell_data["SurfaceRegionId"] = np.asarray(regions, dtype=np.uint8)
    data.cell_data["SurfaceRegion"] = np.where(
        np.asarray(regions) == 0, "CORE", "EXTENSION"
    )
    data.cell_data["RemeshEntityId"] = np.full(
        len(faces), 2, dtype=np.int32
    )
    raw_path = tmp_path / "frozen_raw.vtp"
    open_path = tmp_path / "frozen_open.vtp"
    data.save(raw_path, binary=True)
    data.save(open_path, binary=True)
    boundary = _boundary(0, 0.0, 1.0, "ASSUMED_OUTLET")
    return raw_path, open_path, (boundary,)


@pytest.fixture(scope="module")
def synthetic_vmtk(tmp_path_factory):
    pytest.skip(
        "External synthetic VMTK execution is disabled for the cross-seam task"
    )
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
    # Match the official adaptive boundary-point count so the synthetic seam
    # does not exercise VMTK's unrelated resampling-degeneracy edge case.
    count = 33
    layers = 7
    angle = np.linspace(0.0, 2.0 * np.pi, count, endpoint=False)
    points = np.vstack(
        [
            np.column_stack(
                (
                    (1.15 + 0.035 * np.sin(3.0 * angle)) * np.cos(angle),
                    (0.85 + 0.025 * np.cos(5.0 * angle)) * np.sin(angle),
                    np.full(count, z),
                )
            )
            for z in np.linspace(0.0, 5.0, layers)
        ]
    ).astype(np.float32)
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
    exact_copy = raw_core_exact_copy_qc(paths.open_surface_vtp, paths.raw_vtp)
    promotion = cap_official_vmtk(
        config=config,
        paths=paths,
        tool_script=PROJECT_ROOT / "tools" / "run_vmtk_flowextension.py",
    )
    after = hashlib.sha256(paths.open_surface_vtp.read_bytes()).hexdigest()
    return root, paths, invocation, promotion, before, after, exact_copy


@pytest.fixture(scope="module")
def synthetic_entity_cap(tmp_path_factory):
    """Build a local capped tube without invoking any external VMTK process."""

    root = tmp_path_factory.mktemp("local_entity_cap")
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
    layers = 5
    angle = np.linspace(0.0, 2.0 * np.pi, count, endpoint=False)
    points = np.vstack(
        [
            np.column_stack(
                (1.1 * np.cos(angle), 0.9 * np.sin(angle), np.full(count, z))
            )
            for z in np.linspace(0.0, 4.0, layers)
        ]
    ).astype(np.float32)
    faces: list[list[int]] = []
    regions: list[int] = []
    for layer in range(layers - 1):
        first = layer * count
        second = (layer + 1) * count
        region = 0 if layer < 2 else 1
        for point in range(count):
            following = (point + 1) % count
            faces.extend(
                (
                    [first + point, first + following, second + following],
                    [first + point, second + following, second + point],
                )
            )
            regions.extend((region, region))
    vtk_faces = np.column_stack(
        (np.full(len(faces), 3, dtype=np.int64), np.asarray(faces))
    ).ravel()
    raw = pv.PolyData(points, vtk_faces)
    raw.cell_data["SurfaceRegionId"] = np.asarray(regions, dtype=np.uint8)
    raw.cell_data["SurfaceRegion"] = np.where(
        np.asarray(regions) == 0, "CORE", "EXTENSION"
    )
    raw.save(paths.raw_vtp, binary=True)
    assignment = assign_remesh_entities(paths.raw_vtp)
    pv.read(paths.raw_vtp).save(paths.remeshed_open_vtp, binary=True)

    capped_points = np.vstack(
        (
            points,
            np.asarray(((0.0, 0.0, 0.0), (0.0, 0.0, 4.0)), dtype=np.float32),
        )
    )
    bottom_center = len(points)
    top_center = bottom_center + 1
    top_start = (layers - 1) * count
    capped_faces = list(faces)
    for point in range(count):
        following = (point + 1) % count
        capped_faces.append([bottom_center, following, point])
        capped_faces.append(
            [top_center, top_start + point, top_start + following]
        )
    packed = np.column_stack(
        (
            np.full(len(capped_faces), 3, dtype=np.int64),
            np.asarray(capped_faces),
        )
    ).ravel()
    capped = pv.PolyData(capped_points, packed)
    capped.cell_data["CellEntityIds"] = np.concatenate(
        (
            np.ones(len(faces), dtype=np.int32),
            np.tile(np.asarray((2, 3), dtype=np.int32), count),
        )
    )
    capped.save(paths.capped_vtp, binary=True)
    promotion = SimpleNamespace(
        request={"source_open_vtp": str(paths.remeshed_open_vtp.resolve())},
        runtime={"surface_remesher_called": False},
    )
    return root, paths, assignment, promotion


@pytest.fixture(scope="module")
def clean_open_surface_cap(synthetic_entity_cap):
    """Reuse the local watertight tube; no external capper is executed."""

    _, paths, _, promotion = synthetic_entity_cap
    return paths, promotion


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
    assert (
        formal_config.vmtk.postprocess_mode
        == "cross_seam_active_collar_remesh_then_cap"
    )
    assert formal_config.vmtk.remesh_after_extension is True


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


def test_formal_entity_array_mapping_is_exact(formal_config):
    settings = formal_config.vmtk.entity_remesh
    assert settings.entity_array_name == "RemeshEntityId"
    assert settings.far_core_entity_id == 1
    assert settings.active_entity_id == 2
    assert settings.exclude_entity_ids == (1,)
    assert settings.core_collar.face_layers == 2


def test_formal_entity_remesher_parameters_are_frozen(formal_config):
    settings = formal_config.vmtk.entity_remesh
    assert settings.element_size_mode == "edgelength"
    assert settings.target_edge_length_um == 0.25913916380971913
    assert settings.preserve_boundary_edges is True


def test_parameter_mapping_prohibits_global_remesh(formal_config):
    mapping = parameter_mapping(formal_config.vmtk)
    assert mapping["global_surface_remeshing_performed"] is False
    assert mapping["entity_aware_extension_remeshing_performed"] is True
    assert mapping["entity_remesh"]["expected_entity_ids"] == [1, 2]
    assert mapping["entity_remesh"]["exclude_entity_ids"] == [1]
    assert mapping["entity_remesh"]["active_entity_ids"] == [2]


def test_core_seam_detection(actual_crossseam_assignment):
    raw_path, report, _, _ = actual_crossseam_assignment
    data = pv.read(raw_path).triangulate()
    faces = np.asarray(data.faces).reshape((-1, 4))[:, 1:]
    distances, seeds, seam_edges = core_face_adjacency_layers(
        faces, np.asarray(data.cell_data["SurfaceRegionId"])
    )
    assert report["core_collar_seed_face_count"] == len(seeds)
    assert report["original_cut_seam_edge_count"] == len(seam_edges)
    assert len(seeds) == np.count_nonzero(distances == 0) > 0


def test_core_collar_two_layer_expansion(actual_crossseam_assignment):
    raw_path, report, distances, _ = actual_crossseam_assignment
    data = pv.read(raw_path)
    regions = np.asarray(data.cell_data["SurfaceRegionId"])
    entities = np.asarray(data.cell_data["RemeshEntityId"])
    collar = (regions == 0) & (entities == 2)
    assert report["core_collar_face_layers"] == 2
    assert np.all(np.isin(distances[collar], (0, 1)))
    assert np.count_nonzero(distances == 1) > 0


def test_all_extension_faces_are_active(actual_crossseam_assignment):
    raw_path, _, _, _ = actual_crossseam_assignment
    data = pv.read(raw_path)
    regions = np.asarray(data.cell_data["SurfaceRegionId"])
    entities = np.asarray(data.cell_data["RemeshEntityId"])
    assert np.all(entities[regions == 1] == 2)


def test_far_core_is_excluded_entity_one(actual_crossseam_assignment):
    raw_path, report, _, _ = actual_crossseam_assignment
    data = pv.read(raw_path)
    regions = np.asarray(data.cell_data["SurfaceRegionId"])
    entities = np.asarray(data.cell_data["RemeshEntityId"])
    assert report["far_core_face_count"] > 0
    assert np.all(regions[entities == 1] == 0)


def test_crossseam_assignment_has_only_two_entities(actual_crossseam_assignment):
    _, report, _, _ = actual_crossseam_assignment
    assert report["status"] == "PASS"
    assert report["entity_ids"] == [1, 2]
    assert all(
        row["active_core_collar_face_count"] > 0
        and row["extension_face_count"] > 0
        for row in report["per_port"]
    )


def test_original_cut_seam_is_not_entity_boundary(actual_crossseam_assignment):
    _, report, _, _ = actual_crossseam_assignment
    assert report["original_cut_seam_edge_count"] > 0
    assert report["original_cut_seam_edges_between_different_remesh_entities"] == 0


def test_far_core_active_boundary_is_inside_original_core(
    actual_crossseam_assignment,
):
    _, report, _, _ = actual_crossseam_assignment
    assert report["far_core_active_boundary_edge_count"] > 0
    assert report["far_core_active_boundary_edges_on_original_cut_seam"] == 0


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


def test_entity_remesh_disabled_fails_configuration(tmp_path: Path):
    invalid = tmp_path / "invalid_postprocess.yaml"
    invalid.write_text(
        CONFIG_PATH.read_text(encoding="utf-8").replace(
            "remesh_after_extension: true", "remesh_after_extension: false"
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="vmtk.remesh_after_extension"):
        load_surface_prepare_config(invalid, project_root=PROJECT_ROOT)


def test_entity_target_edge_change_fails_configuration(tmp_path: Path):
    invalid = tmp_path / "invalid_target.yaml"
    invalid.write_text(
        CONFIG_PATH.read_text(encoding="utf-8").replace(
            "target_edge_length_um: 0.25913916380971913",
            "target_edge_length_um: 0.26",
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="target_edge_length_um"):
        load_surface_prepare_config(invalid, project_root=PROJECT_ROOT)


def test_official_source_provenance_records_release_commit(formal_config):
    provenance = official_source_provenance(formal_config.vmtk)
    assert provenance["runtime_release_tag"] == "v1.5.0"
    assert provenance["runtime_release_tag_commit"] == "30d5d7cb8e607d153c208a9d7d39c9feb7985476"


def test_entity_remesh_request_has_exact_two_entity_contract(tmp_path: Path):
    paths = exchange_paths(
        input_directory=tmp_path / "input",
        vmtk_directory=tmp_path / "vmtk",
        geometry_directory=tmp_path / "geometry",
        extension_mode="boundarynormal",
    )
    request = build_entity_remesh_request(
        config=_vmtk_config("boundarynormal", entity_remesh=True), paths=paths
    )
    assert request["expected_entity_ids"] == [1, 2]
    assert request["excluded_entity_ids"] == [1]
    assert request["active_entity_ids"] == [2]
    assert request["preserve_boundary_edges"] is True
    assert request["target_edge_length_um"] == 0.25913916380971913


def test_entity_runner_rejects_configuration_without_crossseam_mode(
    tmp_path: Path,
):
    paths = exchange_paths(
        input_directory=tmp_path / "input",
        vmtk_directory=tmp_path / "vmtk",
        geometry_directory=tmp_path / "geometry",
        extension_mode="boundarynormal",
    )
    with pytest.raises(
        SurfacePrepareError, match="INVALID_VMTK_POSTPROCESS_CONFIGURATION"
    ):
        entity_remesh_official_vmtk(
            config=_vmtk_config("boundarynormal", entity_remesh=False),
            paths=paths,
            tool_script=PROJECT_ROOT / "tools" / "run_vmtk_flowextension.py",
        )


def test_crossseam_preservation_qc_accepts_extension_only_change(
    crossseam_modified_candidate,
):
    _, _, reports = crossseam_modified_candidate
    far_core, boundary, collar, extension = reports
    assert far_core["status"] == "PASS"
    assert boundary["status"] == "PASS"
    assert collar["status"] == "PASS"
    assert extension["status"] == "PASS"
    assert extension["extension_remesh_effect_detected"] is True


def test_crossseam_qc_rejects_unchanged_extension(actual_crossseam_assignment):
    raw_path, _, _, inputs = actual_crossseam_assignment
    far_core, boundary, collar, extension = cross_seam_entity_preservation_qc(
        raw_path, raw_path, inputs.boundaries, restore_region_arrays=False
    )
    assert far_core["status"] == "PASS"
    assert boundary["status"] == "PASS"
    assert collar["status"] == "PASS"
    assert extension["status"] == "FAIL"
    assert extension["extension_remesh_effect_detected"] is False


def test_far_core_exact_qc_catches_vertex_motion(
    actual_crossseam_assignment, tmp_path: Path
):
    raw_path, _, _, _ = actual_crossseam_assignment
    data = pv.read(raw_path).triangulate()
    faces = np.asarray(data.faces).reshape((-1, 4))[:, 1:]
    entities = np.asarray(data.cell_data["RemeshEntityId"])
    far_vertices = np.unique(faces[entities == 1])
    active_vertices = np.unique(faces[entities == 2])
    far_only = np.setdiff1d(far_vertices, active_vertices)
    data.points[int(far_only[0]), 0] += 1.0e-3
    candidate = tmp_path / "far_core_moved.vtp"
    data.save(candidate, binary=True)
    report = locked_entity_preservation_qc(raw_path, candidate)
    assert report["status"] == "FAIL"
    assert not report["checks"]["far_core_vertices_exact_after_output_dtype_cast"]


def test_far_core_exact_qc_catches_connectivity_change(
    actual_crossseam_assignment, tmp_path: Path
):
    raw_path, _, _, _ = actual_crossseam_assignment
    data = pv.read(raw_path).triangulate()
    packed = np.asarray(data.faces).reshape((-1, 4)).copy()
    faces = packed[:, 1:]
    entities = np.asarray(data.cell_data["RemeshEntityId"])
    face_id = int(np.flatnonzero(entities == 1)[0])
    alternatives = np.setdiff1d(np.unique(faces[entities == 1]), faces[face_id])
    packed[face_id, 2] = int(alternatives[0])
    data.faces = packed.ravel()
    candidate = tmp_path / "far_core_connectivity_changed.vtp"
    data.save(candidate, binary=True)
    report = locked_entity_preservation_qc(raw_path, candidate)
    assert report["status"] == "FAIL"
    assert not report["checks"]["far_core_connectivity_unchanged"]


def test_far_core_active_boundary_qc_catches_motion(
    actual_crossseam_assignment, tmp_path: Path
):
    raw_path, _, _, inputs = actual_crossseam_assignment
    data = pv.read(raw_path).triangulate()
    faces = np.asarray(data.faces).reshape((-1, 4))[:, 1:]
    entities = np.asarray(data.cell_data["RemeshEntityId"])
    shared = np.intersect1d(
        np.unique(faces[entities == 1]), np.unique(faces[entities == 2])
    )
    data.points[int(shared[0]), 1] += 1.0e-3
    candidate = tmp_path / "far_active_boundary_moved.vtp"
    data.save(candidate, binary=True)
    _, boundary, _, _ = cross_seam_entity_preservation_qc(
        raw_path, candidate, inputs.boundaries, restore_region_arrays=False
    )
    assert boundary["status"] == "FAIL"
    assert boundary["far_core_active_max_motion_um"] > boundary[
        "machine_precision_tolerance_um"
    ]


def test_simple_closed_loop_detector_rejects_open_chain():
    assert edges_form_simple_closed_loop(
        np.asarray(((0, 1), (1, 2), (2, 3), (3, 0)))
    )
    assert not edges_form_simple_closed_loop(
        np.asarray(((0, 1), (1, 2), (2, 3)))
    )


def test_cut_seam_qc_rejects_an_unchanged_original_ring(
    actual_crossseam_assignment,
):
    raw_path, _, _, inputs = actual_crossseam_assignment
    report = cut_seam_topology_qc(raw_path, raw_path, inputs.boundaries)
    assert report["status"] == "FAIL"
    assert report["original_cut_seam_closed_loop_survives"] is True


def test_crossseam_intersection_qc_accepts_disjoint_triangles(tmp_path: Path):
    data = pv.PolyData(
        np.asarray(
            (
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                (0.0, 0.0, 2.0),
                (1.0, 0.0, 2.0),
                (0.0, 1.0, 2.0),
            )
        ),
        np.asarray((3, 0, 1, 2, 3, 3, 4, 5)),
    )
    data.cell_data["RemeshEntityId"] = np.asarray((1, 2), dtype=np.int32)
    path = tmp_path / "disjoint_crossseam_entities.vtp"
    data.save(path, binary=True)
    report, records = cross_seam_intersection_qc(path)
    assert report["status"] == "PASS"
    assert report["true_self_intersection_count"] == 0
    assert records == []


def test_original_side_collar_distance_uses_complete_frozen_surface(
    frozen_crossseam_tube,
):
    raw_path, open_path, boundaries = frozen_crossseam_tube
    report, hotspots = active_collar_original_side_distance_qc(
        raw_path,
        open_path,
        boundaries,
        maximum_p95_distance_um=0.05,
        cross_seam_neighborhood_um=0.25913916380971913,
    )
    assert report["status"] == "PASS"
    assert report["hard_gate"] is True
    assert report["P50_um"] <= 1.0e-12
    assert report["P95_um"] <= 1.0e-12
    assert report["P99_um"] <= 1.0e-12
    assert report["max_um"] <= 1.0e-12
    assert len(hotspots) == 10
    assert all(
        "distance_from_original_seam_um" in row
        and "within_cross_seam_triangle_neighborhood" in row
        for row in hotspots
    )


def test_three_original_side_cross_sections_preserve_area_and_radius(
    frozen_crossseam_tube,
):
    raw_path, open_path, boundaries = frozen_crossseam_tube
    report = active_collar_cross_section_fidelity_qc(
        raw_path,
        open_path,
        boundaries,
        maximum_equivalent_radius_relative_error=0.08,
    )
    assert report["status"] == "PASS"
    assert report["station_count_per_port"] == 3
    assert report["actual_total_station_count"] == 3
    assert [row["station_offset_in_source_radius"] for row in report["stations"]] == [
        -0.25,
        -0.5,
        -0.75,
    ]
    assert all(row["area_relative_error"] <= 1.0e-12 for row in report["stations"])
    assert all(
        row["radius_relative_error"] <= 1.0e-12 for row in report["stations"]
    )


def test_resume_pipeline_has_no_geometry_regeneration_calls():
    names = set(vmtk_pipeline._run_crossseam_open_resume_impl.__code__.co_names)
    assert "local_plane_cut" not in names
    assert "run_official_vmtk" not in names
    assert "entity_remesh_official_vmtk" not in names
    assert "assign_cross_seam_active_entities" not in names
    assert "cap_official_vmtk" in names


def test_resume_pipeline_reuses_frozen_vtp_and_keeps_legacy_metric_diagnostic():
    source = (
        PROJECT_ROOT / "utils" / "cfd_surface_prepare" / "vmtk_pipeline.py"
    ).read_text(encoding="utf-8")
    function_source = source[
        source.index("def _run_crossseam_open_resume_impl(") :
        source.index("def _post_qc_internal_error_summary(")
    ]
    assert "existing_open_geometry_reused" in function_source
    assert '"geometry_regenerated": False' in function_source
    assert "legacy_reclassified_collar_bidirectional_distance" in function_source
    assert '"hard_gate": False' in function_source


def test_successful_resume_open_qc_is_the_only_route_to_cap():
    passed = {"status": "PASS"}
    failed = {"status": "FAIL"}
    assert resume_open_qc_allows_cap(passed, passed, passed)
    assert not resume_open_qc_allows_cap(passed, failed, passed)
    assert not resume_open_qc_allows_cap()


def test_resume_pipeline_keeps_final_radius_pressure_and_meter_path():
    names = set(vmtk_pipeline._run_crossseam_open_resume_impl.__code__.co_names)
    assert "_radius_qc" in names
    assert "geometry_pressure_correction" in names
    assert "meter_scale_qc" in names
    assert "tag_and_export_final_surface" in names


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


def test_four_open_profiles_are_preserved_after_remesh_qc(
    four_port_remeshed_open,
):
    mesh, boundaries, _, _ = four_port_remeshed_open
    report, mapping = open_profile_qc(mesh, boundaries, distal=True)
    assert report["status"] == "PASS"
    assert report["profile_count"] == 4
    assert report["profile_location"] == "DISTAL_EXTENSION_END"
    assert len(mapping) == 4


def test_direction_length_and_area_qc_survive_remeshed_open_surface(
    four_port_remeshed_open,
):
    mesh, boundaries, proximal, interface = four_port_remeshed_open
    report, distal = extension_geometry_qc(
        mesh, boundaries, proximal, interface
    )
    assert report["status"] == "PASS"
    assert len(distal) == 4
    assert all(
        row["extension_direction_dot"] >= 0.999
        and row["extension_length_relative_error"] <= 0.02
        and row["distal_area_relative_error"] <= 0.05
        for row in report["boundaries"]
    )


def test_crossseam_diagnostic_adapter_preserves_real_count_and_metrics(
    four_port_remeshed_open,
):
    mesh, boundaries, proximal, _ = four_port_remeshed_open
    seam_report = {
        "method": "adjacent face normal jump near original cut plane",
        "boundaries": [
            {
                "port_id": boundary.port_id,
                "adjacent_edge_count": 17 + boundary.index,
                "normal_jump_P50_deg": 1.25 + boundary.index,
                "normal_jump_P95_deg": 2.5 + boundary.index,
                "normal_jump_P99_deg": 3.75 + boundary.index,
                "normal_jump_max_deg": 5.0 + boundary.index,
            }
            for boundary in boundaries
        ],
    }
    normalized = normalize_interface_diagnostic(seam_report)
    report, _ = extension_geometry_qc(mesh, boundaries, proximal, seam_report)
    assert normalized["diagnostic_source"] == "cross_seam_local_normal_jump"
    assert normalized["comparable_to_raw_core_extension_interface"] is False
    assert [row["interface_edge_count"] for row in normalized["boundaries"]] == [
        17,
        18,
        19,
        20,
    ]
    assert [row["normal_jump_P95_deg"] for row in normalized["boundaries"]] == [
        2.5,
        3.5,
        4.5,
        5.5,
    ]
    assert report["status"] == "PASS"
    assert all(
        row["interface_edge_count"] > 0
        and row["interface_diagnostic_available"] is True
        and row["comparable_to_raw_core_extension_interface"] is False
        for row in report["boundaries"]
    )


def test_resume_internal_error_writes_post_qc_failure_summary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    source_run = tmp_path / "source_anchor003274_accepted"
    frozen_open_vtp = (
        source_run
        / "geometry"
        / "vmtk_boundarynormal_crossseam_remeshed_open_um.vtp"
    )
    frozen_open_vtp.parent.mkdir(parents=True)
    frozen_open_vtp.write_bytes(b"immutable frozen OPEN geometry")
    expected_sha = hashlib.sha256(frozen_open_vtp.read_bytes()).hexdigest()
    output_root = tmp_path / "outputs"
    config = SimpleNamespace(paths=SimpleNamespace(output_root=output_root))

    def fail_with_current_schema_bug(*args, **kwargs):
        raise KeyError("interface_edge_count")

    monkeypatch.setattr(
        vmtk_pipeline,
        "_run_crossseam_open_resume_impl",
        fail_with_current_schema_bug,
    )
    with pytest.raises(
        SurfacePrepareError, match="POST_QC_CONTINUATION_INTERNAL_ERROR"
    ):
        run_crossseam_open_resume(
            config,
            project_root=tmp_path,
            source_run=source_run,
            run_id="recovery_internal_error_test",
        )
    summary = read_json(
        output_root / "recovery_internal_error_test" / "qc" / "run_summary.json"
    )
    assert summary["status"] == "POST_QC_CONTINUATION_INTERNAL_ERROR"
    assert summary["stage"] == "post_qc_continuation"
    assert summary["geometry_regenerated"] is False
    assert summary["frozen_source_sha256"] == expected_sha
    assert summary["exception_type"] == "KeyError"
    assert "interface_edge_count" in summary["exception_message"]


def test_runner_uses_pmp_environment(synthetic_vmtk):
    _, _, invocation, _, _, _, _ = synthetic_vmtk
    assert Path(invocation.command[0]).resolve() == PMP_PYTHON.resolve()
    assert Path(invocation.runtime["python_executable"]).resolve() == PMP_PYTHON.resolve()
    assert invocation.runtime["vtk_version"] == "9.2.6"


def test_runner_calls_official_flow_filter(synthetic_vmtk):
    _, _, invocation, _, _, _, _ = synthetic_vmtk
    assert invocation.runtime["official_flow_filter"] == "vtkvmtkPolyDataFlowExtensionsFilter"
    assert invocation.runtime["custom_tps_implementation"] is False


def test_boundarynormal_runner_calls_official_normal_api(synthetic_vmtk):
    _, _, invocation, _, _, _, _ = synthetic_vmtk
    assert (
        invocation.runtime["official_direction_api"]
        == "SetExtensionModeToUseNormalToBoundary"
    )
    assert invocation.runtime["extension_mode_effective"] == "boundarynormal"


def test_boundarynormal_runner_does_not_use_centerlines(synthetic_vmtk):
    _, paths, invocation, _, _, _, _ = synthetic_vmtk
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
    _, paths, _, _, _, _, _ = synthetic_vmtk
    with pytest.raises(SurfacePrepareError, match="INVALID_VMTK_EXTENSION_MODE"):
        run_official_vmtk(
            config=_vmtk_config("invalid"),
            paths=paths,
            tool_script=PROJECT_ROOT / "tools" / "run_vmtk_flowextension.py",
        )


def test_vmtk_raw_output_exists(synthetic_vmtk):
    _, paths, _, _, _, _, _ = synthetic_vmtk
    assert paths.raw_vtp.is_file() and paths.raw_stl.is_file()


def test_vmtk_raw_has_detectable_extensions(synthetic_vmtk):
    _, paths, invocation, _, _, _, _ = synthetic_vmtk
    assert invocation.runtime["raw_cells"] > invocation.runtime["input_cells"]


def test_vmtk_raw_has_two_synthetic_open_profiles(synthetic_vmtk):
    _, paths, _, _, _, _, _ = synthetic_vmtk
    _, mesh = polydata_mesh(paths.raw_vtp)
    report, _ = topology_qc(
        mesh, expected_open_profile_count=2, allow_degenerate=True
    )
    assert report["status"] == "PASS"


def test_cap_only_runtime_proves_global_remesher_was_not_called(synthetic_vmtk):
    _, _, _, promotion, _, _, _ = synthetic_vmtk
    assert promotion.runtime["operation"] == "cap_only"
    assert promotion.runtime["surface_remesher_called"] is False
    assert promotion.runtime["global_surface_remeshing_performed"] is False


def test_vmtk_simple_cap_ids_are_present(synthetic_vmtk):
    _, _, _, promotion, _, _, _ = synthetic_vmtk
    assert promotion.runtime["cap"]["method"] == "simple"
    assert len(promotion.runtime["cap_entity_ids"]) == 2


def test_tagged_caps_keep_project_boundary_metadata(synthetic_vmtk):
    root, paths, _, _, _, _, _ = synthetic_vmtk
    outputs, report = tag_and_export_final_surface(
        paths.capped_vtp,
        (_boundary(0, 0.0, -1.0, "ASSUMED_INLET"), _boundary(1, 5.0, 1.0, "ASSUMED_OUTLET")),
        root / "final",
        root / "boundaries",
        raw_vtp=paths.raw_vtp,
    )
    data = pv.read(outputs["tagged_vtp"])
    assert set(np.unique(data.cell_data["boundary_type_code"])) == {0, 1, 2}
    assert report["distal_boundary_count"] == 2


def test_capped_synthetic_surface_is_watertight(clean_open_surface_cap):
    paths, _ = clean_open_surface_cap
    _, mesh = polydata_mesh(paths.capped_vtp)
    report, _ = topology_qc(mesh, expected_open_profile_count=0)
    assert report["status"] == "PASS"


def test_runner_does_not_modify_exchange_input(synthetic_vmtk):
    _, _, _, _, before, after, _ = synthetic_vmtk
    assert before == after


def test_extension_and_promotion_are_separate_operations(synthetic_vmtk):
    _, _, invocation, promotion, _, _, _ = synthetic_vmtk
    assert invocation.request["operation"] == "extension"
    assert "remeshed_open_vtp" not in invocation.request
    assert "capped_vtp" not in invocation.request
    assert promotion.request["operation"] == "cap_only"
    assert "remeshed_open_vtp" not in promotion.request


def test_boundarynormal_runtime_keeps_tps_parameters(synthetic_vmtk):
    _, _, invocation, _, _, _, _ = synthetic_vmtk
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


def test_intersection_detector_accepts_legal_shared_edge_contact():
    vertices = np.asarray(
        ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0))
    )
    faces = np.asarray(((0, 1, 2), (1, 0, 3)), dtype=np.int64)
    report = triangle_pair_intersection_diagnosis(vertices, faces, 0, 1)
    assert report["shared_vertex_count"] == 2
    assert report["shared_edge_count"] == 1
    assert report["intersection_type"] == "LEGAL_SHARED_EDGE_CONTACT"
    assert report["true_triangle_triangle_intersection"] is False


def test_intersection_segment_calculation_detects_real_penetration():
    vertices = np.asarray(
        (
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.25, -0.25, -1.0),
            (0.25, 0.75, 1.0),
            (0.25, 0.75, -1.0),
        )
    )
    faces = np.asarray(((0, 1, 2), (3, 4, 5)), dtype=np.int64)
    report = triangle_pair_intersection_diagnosis(vertices, faces, 0, 1)
    assert report["shared_vertex_count"] == 0
    assert report["shared_edge_count"] == 0
    assert report["true_triangle_triangle_intersection"] is True
    assert report["intersection_segment_length_um"] > 0.0


def test_saved_previous_penetration_is_raw_layer_zero():
    diagnosis = read_json(
        PREVIOUS_ENTITY_RUN
        / "qc"
        / "previous_entityremesh_intersection_diagnosis.json"
    )
    layer = read_json(
        PREVIOUS_ENTITY_RUN / "qc" / "previous_collision_layer_diagnosis.json"
    )
    assert diagnosis["REAL_GEOMETRIC_PENETRATION"] == "YES"
    assert diagnosis["shared_edge_count"] == 0
    assert diagnosis["intersection_segment_length_um"] > 0.0
    assert layer["port_is_cut_002"] is True
    assert layer["face_adjacency_layer"] == 0


def test_preboundary_plane_qc_uses_absolute_dot(synthetic_vmtk):
    _, paths, _, _, _, _, _ = synthetic_vmtk
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
    ) == SUCCESS_STATUS


def test_raw_geometry_failure_blocks_promotion():
    assert not should_promote_raw_candidate(False)


def test_crossseam_topology_failure_stops_before_cap():
    passed = {"status": "PASS"}
    assert not should_cap_crossseam_open({"status": "FAIL"}, passed, passed)


def test_surviving_original_seam_ring_stops_before_cap():
    passed = {"status": "PASS"}
    assert not should_cap_crossseam_open(passed, passed, {"status": "FAIL"})


def test_clean_crossseam_open_surface_continues_to_cap():
    passed = {"status": "PASS"}
    assert should_cap_crossseam_open(passed, passed, passed)


def test_interface_warning_still_requires_diagnostic_figures():
    assert should_generate_manual_review_figures(
        raw_hard_qc_pass=True, interface_status="WARNING"
    )


def test_boundarynormal_output_names_are_explicit(tmp_path: Path):
    paths = exchange_paths(
        input_directory=tmp_path / "input",
        vmtk_directory=tmp_path / "vmtk",
        geometry_directory=tmp_path / "geometry",
        extension_mode="boundarynormal",
    )
    assert paths.raw_vtp.name == "vmtk_boundarynormal_raw_um.vtp"
    assert paths.raw_stl.name == "vmtk_boundarynormal_raw_um.stl"
    assert (
        paths.remeshed_open_vtp.name
        == "vmtk_boundarynormal_crossseam_remeshed_open_um.vtp"
    )
    assert "remeshed" not in paths.capped_vtp.name


def test_raw_core_points_and_connectivity_are_exact_after_output_cast(synthetic_vmtk):
    _, _, _, _, _, _, report = synthetic_vmtk
    assert report["status"] == "PASS"
    assert report["retained_points_exact_after_output_dtype_cast"] is True
    assert report["original_input_cell_connectivity_changed_count"] == 0
    assert (
        report["retained_input_point_max_motion_um"]
        <= report["machine_precision_tolerance_um"]
    )


def test_raw_surface_region_tagging_is_core_then_extension(synthetic_vmtk):
    _, paths, invocation, _, _, _, _ = synthetic_vmtk
    data = pv.read(paths.raw_vtp).triangulate()
    regions = np.asarray(data.cell_data["SurfaceRegionId"])
    input_cells = invocation.runtime["input_cells"]
    assert np.all(regions[:input_cells] == 0)
    assert np.all(regions[input_cells:] == 1)
    assert set(np.unique(data.cell_data["SurfaceRegion"])) == {"CORE", "EXTENSION"}


def test_cap_only_source_never_references_a_remesher():
    source = (PROJECT_ROOT / "tools" / "run_vmtk_flowextension.py").read_text(
        encoding="utf-8"
    )
    tree = ast.parse(source)
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "_run_cap_only"
    )
    cap_source = ast.get_source_segment(source, function)
    assert cap_source is not None
    assert "vmtkSurfaceRemeshing" not in cap_source
    assert "vtkvmtkPolyDataSurfaceRemeshing" not in cap_source


def test_cap_only_preserves_every_raw_noncap_triangle(synthetic_vmtk, tmp_path: Path):
    _, paths, _, _, _, _, _ = synthetic_vmtk
    _, report = tag_and_export_final_surface(
        paths.capped_vtp,
        (
            _boundary(0, 0.0, -1.0, "ASSUMED_INLET"),
            _boundary(1, 5.0, 1.0, "ASSUMED_OUTLET"),
        ),
        tmp_path / "geometry",
        tmp_path / "boundaries",
        raw_vtp=paths.raw_vtp,
    )
    raw = pv.read(paths.raw_vtp).triangulate()
    assert report["raw_noncap_triangle_count"] == raw.n_cells
    assert report["raw_noncap_triangle_missing_count"] == 0


def test_caponly_final_has_core_extension_and_cap_regions(synthetic_vmtk, tmp_path: Path):
    _, paths, _, _, _, _, _ = synthetic_vmtk
    outputs, _ = tag_and_export_final_surface(
        paths.capped_vtp,
        (
            _boundary(0, 0.0, -1.0, "ASSUMED_INLET"),
            _boundary(1, 5.0, 1.0, "ASSUMED_OUTLET"),
        ),
        tmp_path / "geometry",
        tmp_path / "boundaries",
        raw_vtp=paths.raw_vtp,
    )
    data = pv.read(outputs["tagged_vtp"])
    assert set(np.unique(data.cell_data["SurfaceRegionId"])) == {0, 1, 2}
    assert set(np.unique(data.cell_data["SurfaceRegion"])) == {
        "CORE",
        "EXTENSION",
        "CAP",
    }


def test_formal_boundary_mapping_requires_one_inlet_three_outlets(formal_config):
    inputs = load_surface_inputs(formal_config.paths.cfd_preprocess_run, expected_boundary_count=4)
    assert len(inputs.boundaries) == 4
    assert sum(boundary.role == "ASSUMED_INLET" for boundary in inputs.boundaries) == 1
    assert sum(boundary.role == "ASSUMED_OUTLET" for boundary in inputs.boundaries) == 3


def test_final_winding_consistency_is_a_hard_gate(clean_open_surface_cap):
    paths, _ = clean_open_surface_cap
    _, mesh = polydata_mesh(paths.capped_vtp)
    passed, _ = topology_qc(
        mesh, expected_open_profile_count=0, require_winding_consistent=True
    )
    assert passed["status"] == "PASS"
    broken_faces = np.asarray(mesh.faces).copy()
    broken_faces[0] = broken_faces[0][::-1]
    broken = trimesh.Trimesh(
        vertices=np.asarray(mesh.vertices).copy(),
        faces=broken_faces,
        process=False,
    )
    failed, _ = topology_qc(
        broken, expected_open_profile_count=0, require_winding_consistent=True
    )
    assert failed["checks"]["winding_consistent"] is False
    assert failed["status"] == "FAIL"


def test_bidirectional_distance_exposes_a_final_only_spike(tmp_path: Path):
    points = np.asarray(
        ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0))
    )
    faces = np.asarray(((0, 1, 2), (0, 2, 3)), dtype=np.int64)
    original = trimesh.Trimesh(vertices=points, faces=faces, process=False)
    final_points = np.vstack((points, np.asarray((0.5, 0.0, 2.0))))
    final_faces = np.vstack((faces, np.asarray((0, 1, 4))))
    vtk_faces = np.column_stack(
        (np.full(len(final_faces), 3, dtype=np.int64), final_faces)
    ).ravel()
    final = pv.PolyData(final_points, vtk_faces)
    final.cell_data["SurfaceRegionId"] = np.zeros(len(final_faces), dtype=np.uint8)
    path = tmp_path / "spike.vtp"
    final.save(path, binary=True)
    report = core_symmetric_distance_qc(
        original,
        path,
        (),
        LocalCutConfig(3.0, 2.5, 3.0),
    )
    assert report["original_core_to_caponly_final_core"]["max_um"] == pytest.approx(0.0)
    assert report["caponly_final_core_to_original_core"]["max_um"] > 1.0
    assert report["status"] == "FAIL"


def test_core_normal_evaluator_accepts_exact_caponly_core(synthetic_vmtk, tmp_path: Path):
    _, paths, _, _, _, _, _ = synthetic_vmtk
    outputs, _ = tag_and_export_final_surface(
        paths.capped_vtp,
        (
            _boundary(0, 0.0, -1.0, "ASSUMED_INLET"),
            _boundary(1, 5.0, 1.0, "ASSUMED_OUTLET"),
        ),
        tmp_path / "geometry",
        tmp_path / "boundaries",
        raw_vtp=paths.raw_vtp,
    )
    _, original = polydata_mesh(paths.open_surface_vtp)
    report = normal_consistency_qc(
        original,
        Path(outputs["tagged_vtp"]),
        (),
        LocalCutConfig(3.0, 2.5, 3.0),
    )
    assert report["winding_consistent"] is True
    assert report["core_face_correspondence_missing_count"] == 0
    assert report["status"] == "PASS"


def test_meter_copy_is_exact_for_caponly_export(synthetic_vmtk, tmp_path: Path):
    _, paths, _, _, _, _, _ = synthetic_vmtk
    outputs, _ = tag_and_export_final_surface(
        paths.capped_vtp,
        (
            _boundary(0, 0.0, -1.0, "ASSUMED_INLET"),
            _boundary(1, 5.0, 1.0, "ASSUMED_OUTLET"),
        ),
        tmp_path / "geometry",
        tmp_path / "boundaries",
        raw_vtp=paths.raw_vtp,
    )
    assert meter_scale_qc(
        Path(outputs["manual_review_stl"]), Path(outputs["meter_stl"])
    )["status"] == "PASS"


def test_symmetric_mesh_size_mismatch_is_bidirectional():
    assert symmetric_mesh_size_mismatch(2.0) == 2.0
    assert symmetric_mesh_size_mismatch(0.5) == 2.0


def test_symmetric_mesh_size_mismatch_identity_and_invalid():
    assert symmetric_mesh_size_mismatch(1.0) == 1.0
    assert np.isinf(symmetric_mesh_size_mismatch(0.0))


def test_tail_triangle_counters_detect_catastrophic_sliver(tmp_path: Path):
    points = np.asarray(
        (
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.001, 0.00001, 0.0),
            (0.0, 1.0, 0.0),
        )
    )
    faces = np.asarray(((0, 1, 2), (0, 2, 3)), dtype=np.int64)
    vtk_faces = np.column_stack(
        (np.full(len(faces), 3, dtype=np.int64), faces)
    ).ravel()
    data = pv.PolyData(points, vtk_faces)
    data.cell_data["SurfaceRegionId"] = np.ones(len(faces), dtype=np.uint8)
    path = tmp_path / "sliver_extension.vtp"
    data.save(path, binary=True)
    report = extension_mesh_quality_from_surface(
        path,
        (_boundary(0, 0.0, 1.0, "ASSUMED_INLET"),),
        local_target_edge_um={0: 1.0},
        quality=MeshQualityConfig(
            minimum_triangle_angle_deg=20.0,
            maximum_aspect_ratio=5.0,
            maximum_edge_length_to_local_target_ratio=2.0,
            maximum_neighbor_area_ratio=1.0e9,
            maximum_interface_edge_length_ratio=2.0,
            maximum_bad_triangle_fraction=1.0,
        ),
        method="TAIL_COUNTER_TEST",
        require_symmetric_size_match=False,
    )
    row = report["boundaries"][0]
    assert row["triangle_count_angle_below_5deg"] > 0
    assert row["triangle_count_aspect_above_20"] > 0
    assert report["tail_status"] == "MESH_TAIL_WARNING"


def test_official_cap_preserves_every_remeshed_noncap_triangle(
    synthetic_entity_cap, tmp_path: Path
):
    _, paths, _, promotion = synthetic_entity_cap
    assert Path(promotion.request["source_open_vtp"]) == paths.remeshed_open_vtp.resolve()
    outputs, report = tag_and_export_final_surface(
        paths.capped_vtp,
        (
            _boundary(0, 0.0, -1.0, "ASSUMED_INLET"),
            _boundary(1, 4.0, 1.0, "ASSUMED_OUTLET"),
        ),
        tmp_path / "geometry",
        tmp_path / "boundaries",
        output_stem="synthetic_entityremesh",
        raw_vtp=paths.remeshed_open_vtp,
    )
    remeshed = pv.read(paths.remeshed_open_vtp).triangulate()
    assert report["source_open_noncap_triangle_count"] == remeshed.n_cells
    assert report["source_open_noncap_triangle_missing_count"] == 0
    assert Path(outputs["tagged_vtp"]).is_file()


def test_entityremesh_final_tags_are_far_core_active_and_cap(
    synthetic_entity_cap, tmp_path: Path
):
    _, paths, _, _ = synthetic_entity_cap
    outputs, report = tag_and_export_final_surface(
        paths.capped_vtp,
        (
            _boundary(0, 0.0, -1.0, "ASSUMED_INLET"),
            _boundary(1, 4.0, 1.0, "ASSUMED_OUTLET"),
        ),
        tmp_path / "geometry",
        tmp_path / "boundaries",
        output_stem="synthetic_entityremesh_tags",
        raw_vtp=paths.remeshed_open_vtp,
    )
    data = pv.read(outputs["tagged_vtp"])
    assert report["status"] == "PASS"
    assert set(np.unique(data.cell_data["SurfaceRegionId"])) == {0, 1, 2}
    assert set(np.unique(data.cell_data["RemeshEntityId"])) == {0, 1, 2}
    assert report["remesh_entity_codes"] == {
        "CAP": 0,
        "FAR_CORE": 1,
        "CROSS_SEAM_ACTIVE": 2,
    }


@pytest.mark.parametrize("entity_id", [1, 2])
def test_cap_preserves_each_crossseam_noncap_entity(
    synthetic_entity_cap, tmp_path: Path, entity_id: int
):
    _, paths, _, _ = synthetic_entity_cap
    outputs, _ = tag_and_export_final_surface(
        paths.capped_vtp,
        (
            _boundary(0, 0.0, -1.0, "ASSUMED_INLET"),
            _boundary(1, 4.0, 1.0, "ASSUMED_OUTLET"),
        ),
        tmp_path / f"geometry_{entity_id}",
        tmp_path / f"boundaries_{entity_id}",
        output_stem=f"crossseam_cap_entity_{entity_id}",
        raw_vtp=paths.remeshed_open_vtp,
    )
    source = pv.read(paths.remeshed_open_vtp).triangulate()
    final = pv.read(outputs["tagged_vtp"]).triangulate()
    assert np.count_nonzero(source.cell_data["RemeshEntityId"] == entity_id) == (
        np.count_nonzero(final.cell_data["RemeshEntityId"] == entity_id)
    )


def test_crossseam_final_boundary_tags_are_present(
    synthetic_entity_cap, tmp_path: Path
):
    _, paths, _, _ = synthetic_entity_cap
    outputs, _ = tag_and_export_final_surface(
        paths.capped_vtp,
        (
            _boundary(0, 0.0, -1.0, "ASSUMED_INLET"),
            _boundary(1, 4.0, 1.0, "ASSUMED_OUTLET"),
        ),
        tmp_path / "geometry",
        tmp_path / "boundaries",
        output_stem="crossseam_boundary_tags",
        raw_vtp=paths.remeshed_open_vtp,
    )
    final = pv.read(outputs["tagged_vtp"])
    for name in (
        "boundary_type_code",
        "boundary_index",
        "boundary_origin",
        "port_id",
    ):
        assert name in final.cell_data


def test_meter_scale_dtype_aware_exact_after_cast_passes(
    synthetic_entity_cap, tmp_path: Path
):
    _, paths, _, _ = synthetic_entity_cap
    outputs, _ = tag_and_export_final_surface(
        paths.capped_vtp,
        (
            _boundary(0, 0.0, -1.0, "ASSUMED_INLET"),
            _boundary(1, 4.0, 1.0, "ASSUMED_OUTLET"),
        ),
        tmp_path / "geometry",
        tmp_path / "boundaries",
        output_stem="synthetic_entityremesh_meter",
        raw_vtp=paths.remeshed_open_vtp,
    )
    report = meter_scale_qc(
        Path(outputs["manual_review_stl"]), Path(outputs["meter_stl"])
    )
    assert report["status"] == "PASS"
    assert report["exact_after_float32_cast"] is True


def test_meter_scale_genuine_scaling_error_fails(
    synthetic_entity_cap, tmp_path: Path
):
    _, paths, _, _ = synthetic_entity_cap
    outputs, _ = tag_and_export_final_surface(
        paths.capped_vtp,
        (
            _boundary(0, 0.0, -1.0, "ASSUMED_INLET"),
            _boundary(1, 4.0, 1.0, "ASSUMED_OUTLET"),
        ),
        tmp_path / "geometry",
        tmp_path / "boundaries",
        output_stem="synthetic_entityremesh_bad_meter",
        raw_vtp=paths.remeshed_open_vtp,
    )
    corrupt_path = tmp_path / "incorrect_meter.stl"
    corrupt = trimesh.load_mesh(outputs["meter_stl"], process=False)
    corrupt.vertices = np.asarray(corrupt.vertices) * 1.01
    corrupt.export(corrupt_path)
    assert meter_scale_qc(
        Path(outputs["manual_review_stl"]), corrupt_path
    )["status"] == "FAIL"


def _failure_report(status: str = "PASS") -> dict[str, str]:
    return {"status": status}


def test_meter_failure_maps_to_dedicated_status():
    status = resolve_entity_failure_status(
        topology=_failure_report(),
        boundary_mapping=_failure_report(),
        core_exact=_failure_report(),
        core_distance=_failure_report(),
        normal=_failure_report(),
        radius=_failure_report(),
        meter=_failure_report("FAIL"),
        pressure=_failure_report(),
        integrity=_failure_report(),
    )
    assert status == METER_FAILURE


def test_topology_pass_is_not_renamed_topology_failure_by_meter_qc():
    status = resolve_entity_failure_status(
        topology=_failure_report(),
        boundary_mapping=_failure_report(),
        core_exact=_failure_report(),
        core_distance=_failure_report(),
        normal=_failure_report(),
        radius=_failure_report(),
        meter=_failure_report("FAIL"),
        pressure=_failure_report(),
        integrity=_failure_report(),
    )
    assert status == METER_FAILURE
    assert status != TOPOLOGY_FAILURE


def test_pipeline_recomputes_radius_and_pressure_from_guarded_final():
    source = (PROJECT_ROOT / "utils" / "cfd_surface_prepare" / "vmtk_pipeline.py").read_text(
        encoding="utf-8"
    )
    assert "final_mesh=final_mesh" in source
    assert "geometry_pressure_correction(\n            final_mesh" in source
    assert (
        "extension_pressure_correction_vmtk_boundarynormal_guarded_entityremesh.csv"
        in source
    )


def test_previous_global_remesh_reference_is_read_only():
    root = (
        PROJECT_ROOT
        / "outputs"
        / "cfd_surface_prepare"
        / "vmtk_tps_boundarynormal_anchor003274_20260826_153709"
    )
    paths = (
        root / "input" / "open_surface_um.vtp",
        root / "geometry" / "vmtk_boundarynormal_raw_um.vtp",
        root / "geometry" / "cfd_surface_vmtk_tps_boundarynormal_um.vtp",
    )
    expected = {
        paths[0]: "c132b4e6be39558c6954a332224f58aaff1cf927d1d2ea3cdde0544312d36991",
        paths[1]: "eab3f8a825e7e1413cba63161902e03333de6b176d8ea71071973a4274529323",
    }
    assert all(path.is_file() for path in paths)
    assert all(sha256_file(path) == digest for path, digest in expected.items())
