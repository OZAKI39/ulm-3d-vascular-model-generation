"""Production Tau1 contract, replay, provenance, and visual-package tests."""

from __future__ import annotations

import dataclasses
import json
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest
import yaml

from utils.cfd_flow.apes import (
    compute_lattice_scaling,
    diffusive_time_step,
    generate_musubi_lua,
    inspect_apes_environment,
)
from utils.cfd_flow.config import METHOD, SCHEMA_VERSION, load_cfd_flow_config
from utils.cfd_flow.io import sha256_file
from utils.cfd_flow.physical_port_flux import (
    FLUX_ALGORITHM_REVISION,
    FLUX_DEFINITION,
    evaluate_physical_port_fluxes,
    mesh_origin_dx,
)
from utils.cfd_flow.production import (
    ACCEPTED_ITERATION,
    ACCEPTED_RESTART_SHA256,
    parse_controller_records,
    validate_full_v2,
    validate_local_artifacts,
    validated_scaling_record,
)
from utils.cfd_flow.restart_decode import (
    parse_restart_header,
    read_restart_pdf,
    read_treelm_elemlist,
    reconstruct_macroscopic_field,
    tree_ids_to_ijk,
)
from utils.cfd_flow.visualization import _validate_png, _write_html


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "cfd_flow.yaml"
REGRESSION_CONFIG = PROJECT_ROOT / "configs" / "cfd_flow_promotion_regression.yaml"


@pytest.fixture(scope="module")
def config():
    return load_cfd_flow_config(REGRESSION_CONFIG, project_root=PROJECT_ROOT)


@pytest.fixture(scope="module")
def plane_contract(config):
    return json.loads(config.paths.physical_plane_contract.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def scaling(config, plane_contract):
    area = plane_contract["ports"]["inlet"]["planes"]["central"]["aperture_physical_area_m2"]
    return compute_lattice_scaling(config, None, area)


@pytest.fixture(scope="module")
def lua(config, scaling):
    return generate_musubi_lua(
        config, None, None, scaling, mesh_path="/validated/base/mesh/", maximum_iterations=5000
    )


@pytest.fixture(scope="module")
def replay(config, scaling, plane_contract):
    header = parse_restart_header(config.paths.accepted_base_restart_header)
    pdf = read_restart_pdf(
        config.paths.accepted_base_restart_binary,
        n_elems=header.n_elems,
        n_components=header.n_components,
    )
    field = reconstruct_macroscopic_field(
        pdf,
        dx_m=scaling.dx_m,
        dt_s=scaling.dt_s,
        rho0_kg_m3=config.physics.density_kg_m3,
    )
    tree_ids, _, _ = read_treelm_elemlist(
        config.paths.frozen_base_mesh / "elemlist.lsb", n_elems=header.n_elems
    )
    origin, dx = mesh_origin_dx(config.paths.frozen_base_mesh)
    points = origin + (tree_ids_to_ijk(tree_ids) + 0.5) * dx
    flux = evaluate_physical_port_fluxes(
        plane_contract, points, field.velocity_phy, field.density_lattice, dx_m=dx
    )
    return header, pdf, field, flux


def test_01_schema_v2_load_and_modes(config):
    default = load_cfd_flow_config(DEFAULT_CONFIG, project_root=PROJECT_ROOT)
    assert config.schema_version == SCHEMA_VERSION == 2
    assert config.execution.mode == "VALIDATED_BASE_PROMOTION_REPLAY"
    assert default.execution.mode == "FRESH_STEADY"


def test_02_old_schema_cannot_silently_migrate(tmp_path):
    value = yaml.safe_load(REGRESSION_CONFIG.read_text(encoding="utf-8"))
    value["schema_version"] = 1
    value["method"] = "PROTEUS_COMPATIBLE_SEEDER_MUSUBI_STEADY_LBM_BASELINE"
    path = tmp_path / "old.yaml"
    path.write_text(yaml.safe_dump(value), encoding="utf-8")
    with pytest.raises(ValueError):
        load_cfd_flow_config(path, project_root=PROJECT_ROOT)


def test_03_method_and_base_dx(config):
    assert config.method == METHOD
    assert config.mesh.dx_m == 2.0e-7


def test_04_dt_direct_tau1_formula(config, scaling):
    assert diffusive_time_step(config.mesh.dx_m, config.physics.kinematic_viscosity_m2_s) == pytest.approx(
        config.mesh.dx_m**2 / (6 * config.physics.kinematic_viscosity_m2_s), rel=0, abs=0
    )
    assert scaling.dt_s == pytest.approx(2.038735983690112e-9, rel=1e-15)


def test_05_nu_tau_omega_exact(scaling):
    assert scaling.nu_lattice == pytest.approx(1 / 6, rel=1e-15)
    assert scaling.tau == 1.0
    assert scaling.omega == 1.0


def test_06_dynamic_pressure_reference(config, scaling):
    expected = config.physics.density_kg_m3 * config.physics.lattice_cs_squared * (
        config.mesh.dx_m / scaling.dt_s
    ) ** 2
    assert scaling.pressure_reference_pa == pytest.approx(expected, rel=1e-15)
    assert scaling.pressure_reference_pa == pytest.approx(3387510.7199999993)


def test_07_gauge_pressures_preserved_and_absolute_derived(config, scaling):
    absolute = np.asarray(scaling.outlet_absolute_pressures_pa)
    gauge = np.asarray(config.boundary.outlet_gauge_pressures_pa)
    assert absolute - scaling.pressure_reference_pa == pytest.approx(gauge, abs=3e-10)
    assert absolute == pytest.approx(
        (3387525.2649781005, 3387642.9245492225, 3387497.019373326)
    )


def test_08_target_lattice_formula(config, scaling):
    record = validated_scaling_record(config, scaling)
    expected = config.boundary.target_mass_flow_kg_s / 1056 * scaling.dt_s / config.mesh.dx_m**3
    assert record["target_lattice_flux"] == pytest.approx(expected, rel=0, abs=0)
    assert expected == pytest.approx(0.0006974804380964758)


def test_09_lua_uses_validated_boundaries(lua):
    assert "kind = 'wall_libb'" in lua
    assert "kind = 'adaptive_flux_pressure'" in lua
    assert lua.count("kind = 'pressure_eq'") == 3
    assert "mfr_eq" not in lua


def test_10_lua_fresh_initialization(lua, scaling):
    assert f"pressure_reference_phy = {scaling.pressure_reference_pa:.17g}" in lua
    assert "initial_condition = { pressure = pressure_reference_phy" in lua
    assert all(name + " = 0.0" in lua for name in ("velocityX", "velocityY", "velocityZ"))
    assert "read =" not in lua


def test_11_lua_has_exact_physics(lua, scaling):
    assert f"dx = {scaling.dx_m:.17g}" in lua
    assert f"dt = {scaling.dt_s:.17g}" in lua
    assert "rho0_phy = 1056" in lua
    assert "nu_phy = 3.27e-06" in lua
    assert "bulk_viscosity_phy = 2.1799999999999999e-06" in lua


def test_12_lua_passes_luac(config, lua):
    result = subprocess.run(
        ["wsl.exe", "-d", config.apes.wsl_distribution, "--", "/home/lzy/.local/bin/luac", "-p", "-"],
        input=lua, capture_output=True, text=True, check=False, timeout=30,
    )
    assert result.returncode == 0, result.stderr


def test_13_corrected_binary_sha_is_enforced(config):
    environment = inspect_apes_environment(config.apes)
    assert environment.status == "PASS"
    checks = next(
        record["binary_sha256_checks"]
        for record in environment.commands if "binary_sha256_checks" in record
    )
    assert checks["musubi"]["actual_sha256"] == config.apes.musubi_expected_sha256
    bad = dataclasses.replace(config.apes, musubi_expected_sha256="0" * 64)
    assert inspect_apes_environment(bad).status != "PASS"


def test_14_mesh_and_restart_hash_enforcement(config):
    result = validate_local_artifacts(config)
    assert result["status"] == "PASS"
    assert result["accepted_restart"]["sha256"] == ACCEPTED_RESTART_SHA256
    assert all(value["status"] == "PASS" for value in result["mesh_hashes"].values())


def test_15_physical_plane_contract(plane_contract):
    assert plane_contract["contract_sha256"] == "ffaa49bdb6e43fb7208ff29df07a90d4e92ef9bfa4b96ca4f997d4f453a7f005"
    assert plane_contract["revision"] == "STANDARDIZED_INTERIOR_PHYSICAL_PORT_PLANES_V3"


def test_16_accepted_restart_decode_statistics(replay):
    header, pdf, field, _ = replay
    assert header.iteration == ACCEPTED_ITERATION
    assert np.mean(field.density_lattice) == pytest.approx(1.0000492996793002, rel=1e-14)
    assert np.percentile(field.density_lattice, 1) == pytest.approx(0.9999980621262552, rel=1e-14)
    assert np.percentile(field.density_lattice, 99) == pytest.approx(1.0001547035004035, rel=1e-14)
    assert np.min(pdf) == pytest.approx(0.027777432554529028, rel=1e-14)


def test_17_physical_flux_reproduction(replay):
    _, _, _, flux = replay
    assert flux["flux_definition"] == FLUX_DEFINITION
    assert flux["algorithm_revision"] == FLUX_ALGORITHM_REVISION
    assert flux["Qin_m3_s"] == pytest.approx(2.728393297831303e-15, rel=1e-13)
    assert flux["Qout_m3_s"] == pytest.approx(2.724180907077963e-15, rel=1e-13)
    assert flux["closure"] == pytest.approx(0.0015439089213012042, rel=1e-12)


def test_18_outlet_fluxes_and_fractions(replay):
    flux = replay[3]
    expected_q = (1.417617879220853e-16, 1.9936341084531973e-15, 5.887850107026804e-16)
    expected_f = (0.052038316381177195, 0.7318288235826557, 0.21613286003616722)
    for index, (q, fraction) in enumerate(zip(expected_q, expected_f), start=1):
        label = f"outlet_{index:02d}"
        assert flux["ports"][label]["physical_q_m3_s"] == pytest.approx(q, rel=1e-13)
        assert flux["flow_fractions"][label] == pytest.approx(fraction, rel=1e-13)


def test_19_gauge_pressure_conversion_from_accepted_controller(config, scaling):
    history = json.loads(config.paths.accepted_base_checkpoint_history.read_text(encoding="utf-8"))
    record = next(item for item in history["checkpoint_history"] if item["iteration"] == ACCEPTED_ITERATION)
    inlet_gauge = record["controller"]["pressure_pa"] - scaling.pressure_reference_pa
    assert inlet_gauge == pytest.approx(531.431946845226, abs=1e-9)
    drops = [inlet_gauge - value for value in config.boundary.outlet_gauge_pressures_pa]
    assert drops == pytest.approx((516.8869687439517, 399.2273976220505, 545.1325735185375), abs=1e-9)


def test_20_steady_gates_are_accepted_artifact(config):
    value = json.loads(config.paths.accepted_base_qc.read_text(encoding="utf-8"))
    steady = value["steady_audit"]
    assert steady["status"] == "PASS_NON_REFEREE"
    assert not steady["failed_gates"]
    assert all(steady["gates"].values())


def test_21_full_timestep_v2_lineage(config):
    result = validate_full_v2(config)
    assert result["status"] == "PASS"
    assert result["residual"] == pytest.approx(7.913943402747673e-10)
    assert result["new_solver_calls"] == 0


def test_22_controller_log_parser():
    text = "ADAPTIVE_FLUX_PRESSURE iter=5000 target_lattice= 6.9748043809647912E-004 controlled_lattice= 6.9748043810104361E-004 relative_error= 6.5442568663590096E-012 rho_boundary= 1.0000890779025819E+000 pressure_pa= 3.3878124723499105E+006 max_lattice_velocity= 8.2216649264649707E-006 minimum_pdf= 2.7779565620380503E-002 globBC_count=223"
    record = parse_controller_records(text)[0]
    assert record["iteration"] == 5000
    assert record["active_global_boundary_count"] == 223
    assert record["relative_error"] < 1e-8


def test_23_visual_nonblank_hard_gate(tmp_path, config):
    path = tmp_path / "figure.png"
    figure, axis = plt.subplots(
        figsize=(config.visualization.width_px / config.visualization.dpi, config.visualization.height_px / config.visualization.dpi),
        dpi=config.visualization.dpi,
    )
    axis.plot([0, 1], [0, 1])
    axis.set_title("finite nonblank CFD test")
    figure.savefig(path)
    plt.close(figure)
    assert _validate_png(path, config)["status"] == "PASS"


def test_24_offline_html_generation_and_disclaimers(tmp_path, config):
    metrics = {
        "Qin_m3_s": 2.728e-15, "physical_volume_closure": 0.0015,
        "inlet_gauge_pressure_pa": 531.4, "rho_mean": 1.00005,
        "flow_fractions": {"outlet_01": 0.052, "outlet_02": 0.732, "outlet_03": 0.216},
    }
    steady = {"R_mass_short": 0.0015, "R_mass_long": 0.0015, "R_velocity": 1e-10, "R_pressure": 1e-10, "R_inlet": 0.0031}
    full = {"status": "PASS", "residual": 7.9e-10, "gate": 1e-8}
    coarse = {"maximum_absolute_percent_difference": 1.45}
    path = tmp_path / "production_review.html"
    _write_html(path, config, metrics, steady, full, coarse, "PASS")
    text = path.read_text(encoding="utf-8")
    assert "Formal three-grid GCI was not completed because Fine steady computation was terminated under resource-budget constraints." in text
    assert "DEFERRED_TO_POST_GRID_PRODUCTION_VALIDATION" in text
    assert "https://" not in text and "http://" not in text


def test_25_visual_manifest_contract_is_implemented():
    source = (PROJECT_ROOT / "utils/cfd_flow/visualization.py").read_text(encoding="utf-8")
    for key in ("filename", "purpose", "source_data", "units", "raw_min", "raw_max", "display_min", "display_max", "sha256", "width_px", "height_px", "status"):
        assert key in source
    assert "visual_manifest.json" in source


def test_26_fine_transient_cannot_be_steady_source():
    source = (PROJECT_ROOT / "utils/cfd_flow/pipeline.py").read_text(encoding="utf-8")
    assert '"steady_solution_source": "VALIDATED_RESEARCH_BASE_ACCEPTED_RESTART"' in source
    assert '"fresh_full_production_steady_solve": False' in source
    assert '"fine_transient_used": False' in (PROJECT_ROOT / "utils/cfd_flow/visualization.py").read_text(encoding="utf-8")


def test_27_three_grid_cannot_be_false_pass(config):
    evidence = json.loads(config.paths.coarse_base_grid_evidence.read_text(encoding="utf-8"))
    assert evidence["formal_asymptotic_grid_convergence"] is False
    assert evidence["three_grid_metrics"] == "NOT_AVAILABLE_FINE_STEADY_NOT_COMPLETED"
    assert "grid independent proven" in evidence["claims_not_made"]


def test_28_wss_cannot_be_false_validated():
    source = (PROJECT_ROOT / "utils/cfd_flow/pipeline.py").read_text(encoding="utf-8")
    assert '"WSS_status": "DEFERRED_TO_POST_GRID_PRODUCTION_VALIDATION"' in source


def test_29_corrected_source_and_patch_provenance(config):
    evidence = json.loads(
        (PROJECT_ROOT / "outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_grid_convergence_anchor003274_20260901/qc/fine_adaptive_target_fix_validation.json").read_text(encoding="utf-8")
    )
    assert evidence["candidate_binary_sha256"] == config.apes.musubi_expected_sha256
    assert evidence["candidate_patch_sha256"] == config.apes.musubi_patch_sha256
    assert evidence["candidate_source_sha256"] == config.apes.musubi_patched_source_sha256


def test_30_production_csv_and_vtu_contract_are_implemented():
    production = (PROJECT_ROOT / "utils/cfd_flow/production.py").read_text(encoding="utf-8")
    assert "production_primary_metrics.csv" in (PROJECT_ROOT / "utils/cfd_flow/pipeline.py").read_text(encoding="utf-8")
    for field in ("velocity_phy", "velocity_magnitude_m_s", "velocity_magnitude_mm_s", "pressure_gauge_pa", "pressure_absolute_solver_pa", "rho_lattice"):
        assert field in production


def test_31_no_old_production_contract_tokens(lua):
    combined = lua + DEFAULT_CONFIG.read_text(encoding="utf-8") + REGRESSION_CONFIG.read_text(encoding="utf-8")
    assert "MUSUBI_ONLY_RECOVERY" not in combined
    assert "reference_dx_m" not in combined
    assert "reference_dt_s" not in combined
    assert "mfr_eq" not in combined


def test_32_mesh_hash_values_are_exact(config):
    assert sha256_file(config.paths.frozen_base_mesh / "elemlist.lsb") == config.mesh.elemlist_sha256
    assert sha256_file(config.paths.frozen_base_mesh / "bnd.lsb") == config.mesh.bnd_sha256
    assert sha256_file(config.paths.frozen_base_mesh / "qval.lsb") == config.mesh.qval_sha256
