"""Wrapper-level tests for the one production Seeder/Musubi flow stage."""

from __future__ import annotations

import math
import subprocess
from datetime import datetime
from pathlib import Path

import numpy as np
import pyvista as pv
import pytest

from utils.cfd_flow.apes import (
    bulk_viscosity_from_kinematic,
    compute_lattice_scaling,
    diffusive_time_step,
    generate_musubi_lua,
    generate_seeder_lua,
    load_boundary_conditions,
    parse_mesh_header,
)
from utils.cfd_flow.config import METHOD, load_cfd_flow_config
from utils.cfd_flow.geometry import (
    BOUNDARY_LABELS,
    cells_across_diameter,
    compute_bounding_cube,
    find_seed_point,
    parabolic_velocity,
    partition_surface,
)
from utils.cfd_flow.io import create_run_layout, load_flow_inputs
from utils.cfd_flow.qc import (
    evaluate_mass_conservation,
    validate_and_convert_flow_vtu,
    write_proteus_metadata,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = PROJECT_ROOT / "configs" / "cfd_flow.yaml"


@pytest.fixture(scope="module")
def config():
    return load_cfd_flow_config(CONFIG_PATH, project_root=PROJECT_ROOT)


@pytest.fixture(scope="module")
def inputs(config):
    return load_flow_inputs(config.paths.source_surface_run)


@pytest.fixture(scope="module")
def partition(inputs, tmp_path_factory):
    return partition_surface(inputs, tmp_path_factory.mktemp("solver_patches"))


@pytest.fixture(scope="module")
def boundary_conditions(inputs):
    return load_boundary_conditions(inputs.boundary_conditions)


@pytest.fixture(scope="module")
def scaling(config, boundary_conditions, partition):
    area = partition.patch("inlet").area_um2 * 1.0e-12
    return compute_lattice_scaling(config, boundary_conditions, area)


def test_01_cfd_flow_config_parse(config):
    assert config.method == METHOD
    assert config.mesh.dx_target_um == 0.20


def test_02_final_surface_path_resolution(inputs):
    assert inputs.tagged_surface_vtp.name.endswith("_um.vtp")
    assert inputs.meter_surface_stl.name.endswith("_m.stl")
    assert inputs.tagged_surface_vtp.is_file()


def test_03_cell_entity_ids_map_to_five_patches(partition):
    assert tuple(item.label for item in partition.patches) == BOUNDARY_LABELS
    assert {item.entity_id for item in partition.patches} == {1, 2, 3, 4, 5}


def test_04_patch_union_exact_reconstruction(partition):
    assert partition.qc["status"] == "PASS"
    assert partition.qc["missing_triangles"] == 0
    assert partition.qc["duplicate_triangles"] == 0
    assert partition.qc["patch_triangle_count_sum"] == len(partition.faces)


def test_05_um_to_m_exact_scaling(partition):
    inlet = partition.patch("inlet")
    loaded = pv.read(inlet.path_m)
    source_points = partition.points_um[np.unique(partition.faces[inlet.face_indices])]
    assert np.allclose(np.sort(loaded.bounds), np.sort(pv.PolyData(source_points * 1.0e-6).bounds), rtol=0.0, atol=2e-12)
    assert partition.qc["translation_applied"] is False


def test_06_seed_point_inside_lumen(partition):
    seed = find_seed_point(partition)
    assert seed.seed_inside_lumen is True
    assert seed.candidate_offset_radius in {0.5, 1.0, 2.0}


def test_07_boundary_labels_are_exact():
    assert BOUNDARY_LABELS == ("wall", "inlet", "outlet_01", "outlet_02", "outlet_03")


def test_08_flow_rate_to_mean_velocity(boundary_conditions, partition, scaling):
    expected = boundary_conditions.inlet_flow_m3_s / (partition.patch("inlet").area_um2 * 1e-12)
    assert scaling.velocity_mean_m_s == pytest.approx(expected)


def test_09_flow_rate_to_mass_flow(boundary_conditions):
    expected = 1056.0 * 7.693508475538942e-16
    assert boundary_conditions.density_kg_m3 * boundary_conditions.inlet_flow_m3_s == pytest.approx(expected)


def test_10_parabolic_profile_calculation():
    velocity = parabolic_velocity(
        np.array([0.5, 0.0, 0.0]),
        np.zeros(3),
        np.array([0.0, 0.0, 1.0]),
        1.0,
        2.0,
    )
    assert velocity == pytest.approx([0.0, 0.0, 1.5])
    outside = parabolic_velocity(np.array([2.0, 0.0, 0.0]), np.zeros(3), np.array([0.0, 0.0, 1.0]), 1.0, 2.0)
    assert outside == pytest.approx([0.0, 0.0, 0.0])


def test_11_arbitrary_3d_inlet_normal():
    normal = np.array([1.0, 2.0, 3.0])
    normal /= np.linalg.norm(normal)
    velocity = parabolic_velocity(np.zeros(3), np.zeros(3), normal, 1.0, 4.0)
    assert velocity == pytest.approx(normal * 4.0)
    assert np.cross(velocity, normal) == pytest.approx(np.zeros(3), abs=1e-12)


def test_12_proteus_diffusive_dt_scaling(config):
    dt = diffusive_time_step(2.0e-7, config.physics.reference_dx_m, config.physics.reference_dt_s)
    assert dt == pytest.approx(2.44140625e-8)


def test_13_lattice_nu_tau_omega(scaling):
    assert scaling.nu_lattice == pytest.approx(1.995849609375)
    assert scaling.tau == pytest.approx(scaling.nu_lattice / (1 / 3) + 0.5)
    assert scaling.omega == pytest.approx(1.0 / scaling.tau)
    assert 0.0 < scaling.omega < 2.0


def test_14_expected_mach_calculation(scaling):
    expected = scaling.velocity_max_expected_m_s * scaling.dt_s / scaling.dx_m / math.sqrt(1 / 3)
    assert scaling.mach_max_expected == pytest.approx(expected)
    assert scaling.mach_max_expected < 0.05


def test_15_gauge_pressure_common_offset(scaling, boundary_conditions):
    offsets = np.asarray(scaling.outlet_absolute_pressures_pa) - np.asarray(boundary_conditions.outlet_gauge_pressures_pa)
    assert offsets == pytest.approx(np.full(3, scaling.pressure_reference_pa))


def test_16_pressure_differences_preserved(scaling, boundary_conditions):
    gauge = np.asarray(boundary_conditions.outlet_gauge_pressures_pa)
    absolute = np.asarray(scaling.outlet_absolute_pressures_pa)
    assert gauge[:, None] - gauge[None, :] == pytest.approx(absolute[:, None] - absolute[None, :], abs=1e-10)


def test_17_negative_gauge_pressure_is_safe(scaling, boundary_conditions):
    assert min(boundary_conditions.outlet_gauge_pressures_pa) < 0.0
    assert min(scaling.outlet_absolute_pressures_pa) > 0.0
    assert min(scaling.outlet_lattice_densities) > 0.0


def test_18_musubi_lua_generation(config, partition, boundary_conditions, scaling):
    text = generate_musubi_lua(config, partition, boundary_conditions, scaling)
    assert "kind = 'fluid'" in text
    assert "layout = 'd3q19'" in text
    assert "relaxation = 'bgk'" in text
    assert "kind = 'wall_libb'" in text
    assert "kind = 'mfr_eq'" in text
    assert "mass_flowrate =" in text
    assert "massflowrate" not in text
    assert text.count("kind = 'pressure_eq'") == 3
    assert "bulk_viscosity_phy = (2.0 / 3.0) * nu_phy" in text
    assert "bulk_viscosity = bulk_viscosity_phy" in text


def test_19_seeder_lua_generation(config, partition):
    seed = find_seed_point(partition)
    points = partition.points_um
    bounds = (points[:, 0].min(), points[:, 0].max(), points[:, 1].min(), points[:, 1].max(), points[:, 2].min(), points[:, 2].max())
    cube = compute_bounding_cube(bounds, config.mesh.dx_target_m, 4)
    text = generate_seeder_lua(partition, seed, cube)
    for label in BOUNDARY_LABELS:
        assert f"label = '{label}'" in text
        assert f"{label}.stl' }}}}" in text
    assert text.count("calc_dist = true") == 1
    assert "kind = 'seed'" in text
    assert text.count("{") == text.count("}")


def test_20_vtu_velocity_phy_field_validation(tmp_path):
    grid = pv.ImageData(dimensions=(3, 3, 3), spacing=(2e-7, 2e-7, 2e-7)).cast_to_unstructured_grid()
    grid.cell_data["velocity_phy"] = np.tile([1e-4, 0.0, 0.0], (grid.n_cells, 1))
    grid.cell_data["pressure_phy"] = np.full(grid.n_cells, 100.0)
    source = tmp_path / "source.vtu"
    output = tmp_path / "flow_field.vtu"
    grid.save(source)
    _, result = validate_and_convert_flow_vtu(source, output, pressure_reference_pa=90.0)
    assert result["velocity_phy_components"] == 3
    assert output.is_file()


def test_21_proteus_metadata(tmp_path):
    flow = tmp_path / "flow.vtu"
    flow.touch()
    metadata = write_proteus_metadata(tmp_path / "proteus.json", inlet_equivalent_diameter_m=3e-6, source_flow_vtu=flow)
    assert metadata["lengthUnit"] == 1.0
    assert metadata["velocityUnit"] == 1.0
    assert metadata["velocityField"] == "velocity_phy"
    assert metadata["inletNormal"] is None


def test_22_mass_conservation_evaluator():
    result = evaluate_mass_conservation(10.0, (2.0, 3.0, 5.0))
    assert result["relative_error"] == 0.0
    assert result["flow_signs_pass"] is True


def test_23_no_automatic_resolution_sweep(config):
    assert config.mesh.automatic_resolution_sweep is False
    assert config.mesh.dx_target_um == 0.20


def test_24_no_geometry_regeneration():
    entry = (PROJECT_ROOT / "cfd_flow.py").read_text(encoding="utf-8")
    pipeline = (PROJECT_ROOT / "utils" / "cfd_flow" / "pipeline.py").read_text(encoding="utf-8")
    forbidden_calls = ("run_cfd_preprocess(", "run_vmtk_surface_prepare(", "ultraVessMorpho2Mesh")
    assert not any(value in entry + pipeline for value in forbidden_calls)
    assert "find_seed_point" not in pipeline


def test_25_boundary_manifest_triangle_counts(partition):
    assert {item.label: item.triangle_count for item in partition.patches} == {
        "wall": 67071,
        "inlet": 56,
        "outlet_01": 44,
        "outlet_02": 42,
        "outlet_03": 49,
    }


def test_26_bounding_cube_is_power_of_two_and_safe(config, partition):
    points = partition.points_um
    bounds = (points[:, 0].min(), points[:, 0].max(), points[:, 1].min(), points[:, 1].max(), points[:, 2].min(), points[:, 2].max())
    cube = compute_bounding_cube(bounds, config.mesh.dx_target_m, 4)
    assert cube.cells_per_axis == 2**cube.level
    assert cube.margin_cells_minimum >= 4.0


def test_27_cells_across_diameter_report(partition):
    result = cells_across_diameter(partition, 0.20)
    assert len(result["per_port"]) == 4
    assert result["minimum"] > 10.0


def test_28_corrected_outlet_pressure_source(boundary_conditions):
    assert boundary_conditions.outlet_gauge_pressures_pa == pytest.approx(
        (14.544978101274268, 132.20454922317552, -13.700626673311461)
    )


def test_29_official_boundary_cell_count_parser(tmp_path):
    (tmp_path / "header.lua").write_text(
        """nElems = 12
minLevel = 9
maxLevel = 9
property = {
  { label = 'has boundaries', bitpos = 0, nElems = 4 }
}
""",
        encoding="utf-8",
    )
    (tmp_path / "bnd.lua").write_text(
        """nSides = 3
nBCtypes = 5
bclabel = { 'wall', 'inlet', 'outlet_01', 'outlet_02', 'outlet_03' }
""",
        encoding="utf-8",
    )
    boundary_ids = np.asarray(
        (
            (1, 2, 0),
            (1, 3, 0),
            (4, 0, 0),
            (5, 0, 0),
        ),
        dtype="<i8",
    )
    (tmp_path / "bnd.lsb").write_bytes(boundary_ids.tobytes())
    result = parse_mesh_header(tmp_path)
    assert result["fluid_element_count"] == 12
    assert result["minimum_level"] == result["maximum_level"] == 9
    assert result["boundary_cell_counts"] == {
        "wall": 2,
        "inlet": 1,
        "outlet_01": 1,
        "outlet_02": 1,
        "outlet_03": 1,
    }


def test_30_recovery_run_directory_is_separate(tmp_path):
    layout = create_run_layout(
        tmp_path,
        timestamp=datetime(2026, 8, 28, 17, 0, 0),
        recovery=True,
    )
    assert layout.root.name == "musubi_recovery_anchor003274_20260828_170000"


def test_31_bulk_viscosity_policy_does_not_change_lattice_scaling(config, scaling):
    before = (scaling.dt_s, scaling.nu_lattice, scaling.tau, scaling.omega)
    bulk = bulk_viscosity_from_kinematic(3.27e-6)
    after = (scaling.dt_s, scaling.nu_lattice, scaling.tau, scaling.omega)
    assert bulk == pytest.approx(2.18e-6, rel=0.0, abs=1.0e-18)
    assert before == after
    assert config.physics.bulk_viscosity_source == "MUSUBI_D3Q19_REQUIRED_EXPLICIT_PARAMETER"
    assert config.physics.bulk_viscosity_policy == "OFFICIAL_MUSUBI_BASELINE_TWO_THIRDS_KINEMATIC_VISCOSITY"


def test_32_generated_musubi_lua_passes_luac(config, partition, boundary_conditions, scaling):
    text = generate_musubi_lua(config, partition, boundary_conditions, scaling)
    result = subprocess.run(
        ["wsl.exe", "-d", config.apes.wsl_distribution, "--", "/home/lzy/.local/bin/luac", "-p", "-"],
        input=text,
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert result.returncode == 0, result.stderr


def test_33_solver_recovery_directory_is_separate(tmp_path):
    layout = create_run_layout(
        tmp_path,
        timestamp=datetime(2026, 8, 28, 18, 0, 0),
        solver_recovery=True,
    )
    assert layout.root.name == "musubi_solver_recovery_anchor003274_20260828_180000"
    assert not layout.seeder_mesh.exists()
