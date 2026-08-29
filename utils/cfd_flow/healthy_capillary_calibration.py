"""Isolated healthy mouse cortical-capillary adaptive-flux calibration.

This module intentionally does not alter the production CFD pipeline.  It
continues the accepted axis-aligned mesh from the frozen iteration-159326
restart, runs the pinned adaptive-flux Musubi binary at one pre-registered
healthy continuum target, and independently audits the final restart PDFs.
"""

from __future__ import annotations

import math
import shutil
from datetime import datetime
from pathlib import Path, PurePosixPath
from typing import Any, Iterable

import numpy as np

from .adaptive_flux_pressure_audit import (
    musubi_pressure_flux_affine_coefficients,
    physical_mass_factor,
    reconstruct_musubi_boundary_state,
    solve_boundary_density,
)
from .adaptive_flux_steady import (
    AXIS_MESH_RUN,
    BULK_NU_M2_S,
    EXPECTED_FLUID_CELLS,
    EXPECTED_INLET_GLOBBC,
    EXPECTED_RESTART_BYTES,
    MPI_RANKS,
    NU_M2_S,
    PRESSURE_REFERENCE_PA,
    RUNTIME_BASE_WSL,
    SHORT_SIMULATION_NAME,
    WRAPPER_MARGIN_S,
    CheckpointTracker,
    _archive_runtime,
    _git,
    _maximum_runtime_output_path_length,
    _one_wsl_health_check,
    _predicted_maximum_runtime_path_length,
    _prepare_short_tracking_root,
    _run_luac,
    _run_monitored_wsl,
    _wsl_path_to_unc,
    summarize_controller_csv,
)
from .adaptive_flux_validation import (
    BINARY_WSL,
    EXPECTED_BINARY_SHA256,
    EXPECTED_DT_S,
    MPIRUN_WSL,
)
from .apes import parse_mesh_header, windows_to_wsl
from .config import load_cfd_flow_config
from .diagnostics import parse_official_steady_termination
from .exact_link_flux import (
    CS2,
    EXPECTED_DX_M,
    _file_manifest,
    _link_rows,
    _write_link_csv,
    build_coordinate_lookup,
    equilibrium_pdf,
    pull_fetch_pdfs,
    reconstruct_boundary,
    signed_mass_balance,
    velocity_from_pdf,
)
from .io import FlowError, sha256_file, write_json
from .mcclure_adaptive_flux_reference import physical_volume_flux_to_lattice
from .port_flux_audit import (
    PORT_LABELS,
    extract_boundary_property_indices,
    parse_bnd_header,
    parse_boundary_property_header,
    read_boundary_ids,
)
from .restart_decode import (
    D3Q19_DIRECTIONS,
    parse_restart_header,
    read_restart_pdf,
    read_treelm_elemlist,
    restart_binary_size_contract,
    tree_ids_to_ijk,
    tree_levels,
)


RUN_PREFIX = "healthy_mouse_capillary_calibration_anchor003274"
SOURCE_RUN = "mcclure_adaptive_flux_resume_0p5_anchor003274_20260829_170926"
RESUME_ITERATION = 159_326
RESUME_BINARY_NAME = "a3274_3.890E-03.lsb"
EXPECTED_RESUME_SHA256 = "c911dd8085fea971590758bb0732cc2245f39f9c3720a283215237de0280d3d0"
EXPECTED_FIRST_CONTROLLER_ITERATION = RESUME_ITERATION + 1

REFERENCE_DENSITY_KG_M3 = 1056.0
SMOOTH_INLET_AREA_M2 = 7.819752111687344e-12
PRIMARY_MEAN_VELOCITY_M_S = 3.5e-4
OLD_TARGET_SCALE_FACTOR = 3.557431889225085
TARGET_BRACKET_M_S = {"LOW": 3.0e-4, "PRIMARY": 3.5e-4, "HIGH": 4.0e-4}
TARGET_Q_M3_S = PRIMARY_MEAN_VELOCITY_M_S * SMOOTH_INLET_AREA_M2
TARGET_MASS_FLOW_KG_S = REFERENCE_DENSITY_KG_M3 * TARGET_Q_M3_S

CONVERGENCE_INTERVAL = 100
CONVERGENCE_NVALS = 100
EFFECTIVE_CONVERGENCE_WINDOW_ITERATIONS = 10_000
EARLIEST_FINAL_STEADY_ITERATION = RESUME_ITERATION + EFFECTIVE_CONVERGENCE_WINDOW_ITERATIONS
VELOCITY_THRESHOLD_M_S = 1.75e-6
PRESSURE_THRESHOLD_PA = 0.08
CHECKPOINT_INTERVAL = 10_000
WALLCLOCK_LIMIT_S = 3_600
MAXIMUM_ITERATIONS = 1_000_000

OUTLET_GAUGE_PRESSURE_PA = {
    "outlet_01": 14.544978101274268,
    "outlet_02": 132.20454922317552,
    "outlet_03": -13.700626673311461,
}
CORRECTED_RUNTIME_NORMALS = {
    "outlet_01": np.asarray((0, -1, -1), dtype=np.int8),
    "outlet_02": np.asarray((1, 0, 1), dtype=np.int8),
    "outlet_03": np.asarray((-1, 0, 0), dtype=np.int8),
}
EXPECTED_PRESSURE_COUNTS = {
    "outlet_01": (178, 178, 0),
    "outlet_02": (179, 161, 18),
    "outlet_03": (193, 193, 0),
}
INLET_PRESSURE_NEIGHBOR_DIRECTION = np.asarray((0, 0, -1), dtype=np.int8)

INLET_TARGET_HARD_GATE = 0.01
INLET_TARGET_PREFERRED_GATE = 0.001
MASS_ACCUMULATION_HARD_GATE = 0.01
MASS_ACCUMULATION_PREFERRED_GATE = 0.001
PRESSURE_CROSSCHECK_GATE = 0.001
BOUNDARY_BALANCE_HARD_GATE = 0.01
BOUNDARY_BALANCE_PREFERRED_GATE = 0.001
OUTLET_02_BACKFLOW_FRACTION_GATE = 0.05
HEALTHY_UMEAN_RANGE_M_S = (3.0e-4, 4.0e-4)
HEALTHY_WSS_RANGE_DYN_CM2 = (17.0, 97.0)

STATUS_PASS = "CFD_FLOW_HEALTHY_MOUSE_CAPILLARY_CALIBRATION_PASS"
STATUS_MASS_INCOMPLETE = "CFD_FLOW_HEALTHY_CALIBRATION_MASS_STEADY_INCOMPLETE"
STATUS_OUTLETS_UNRESOLVED = "CFD_FLOW_HEALTHY_CALIBRATION_OUTLET_PRESSURES_UNRESOLVED"
STATUS_SAFETY_FAILED = "CFD_FLOW_HEALTHY_CALIBRATION_NUMERICAL_SAFETY_FAILED"
STATUS_FAILED = "CFD_FLOW_HEALTHY_MOUSE_CAPILLARY_CALIBRATION_FAILED"
NEXT_GRID = "RUN HEALTHY ADAPTIVE-FLUX GRID CONVERGENCE"
NEXT_RESUME = "RESUME HEALTHY FLOW FROM LATEST CHECKPOINT"
NEXT_OUTLETS = "RECALIBRATE OUTLET PRESSURE SPLIT"
NEXT_SAFETY = "REVIEW HEALTHY TARGET LATTICE SCALING"
NEXT_REVIEW = "REVIEW HEALTHY CALIBRATION AUDIT FAILURE"


class HealthyCheckpointTracker(CheckpointTracker):
    """Checkpoint tracker whose evidence records the healthy 10k interval."""

    def _write(self) -> None:
        write_json(
            self.manifest_path,
            {
                "status": "PASS",
                "checkpoint_interval_iterations": CHECKPOINT_INTERVAL,
                "runtime_root_wsl": self.runtime_root_wsl,
                "records": self.records,
            },
        )


def healthy_literature_reference_contract() -> dict[str, Any]:
    """Return only values verified against the listed primary articles."""

    studies = [
        {
            "title": "Quantitative imaging of red blood cell velocity in vivo using optical coherence Doppler tomography",
            "authors": ["Hugang Ren", "Congwu Du", "Kicheon Park", "Nora D. Volkow", "Yingtian Pan"],
            "year": 2012,
            "doi": "10.1063/1.4726115",
            "url": "https://pubmed.ncbi.nlm.nih.gov/22904572/",
            "exact_value_used": {
                "condition": "healthy mouse cerebral cortex under normocapnia",
                "rbc_velocity_mean_mm_s": 0.72,
                "rbc_velocity_sd_mm_s": 0.15,
            },
            "role": "experimental RBC-velocity reference",
        },
        {
            "title": "More homogeneous capillary flow and oxygenation in deeper cortical layers correlate with increased oxygen extraction",
            "authors": [
                "Baoqiang Li", "Tatiana V. Esipova", "Ikbal Sencan", "Kivılcım Kılıç",
                "Buyin Fu", "Michele Desjardins", "Mohammad Moeini", "Sreekanth Kura",
                "Mohammad A. Yaseen", "Frederic Lesage", "Leif Østergaard", "Anna Devor",
                "David A. Boas", "Sergei A. Vinogradov", "Sava Sakadžić",
            ],
            "year": 2019,
            "doi": "10.7554/eLife.42299",
            "url": "https://elifesciences.org/articles/42299",
            "exact_value_used": {
                "condition": "awake mouse cortex",
                "line_scan_mean_rbc_speed_mm_s": 0.61,
                "point_rbc_passage_mean_speed_mm_s": 0.69,
            },
            "role": "experimental RBC-velocity reference; modalities retained separately",
        },
        {
            "title": "A Pilot Study Investigating Changes in Capillary Hemodynamics and Its Modulation by Exercise in the APP-PS1 Alzheimer Mouse Model",
            "authors": [
                "Xuecong Lu", "Mohammad Moeini", "Baoqiang Li", "Yuankang Lu",
                "Rafat Damseh", "Philippe Pouliot", "Éric Thorin", "Frédéric Lesage",
            ],
            "year": 2019,
            "doi": "10.3389/fnins.2019.01261",
            "url": "https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2019.01261/full",
            "exact_value_used": {
                "condition": "wild-type control",
                "capillary_diameter_mean_um": 4.67,
                "capillary_diameter_sem_um": 0.06,
                "rbc_velocity_mean_mm_s": 0.87,
                "rbc_velocity_sem_mm_s": 0.03,
            },
            "role": "experimental WT diameter and RBC-velocity reference",
        },
        {
            "title": "Cortical microvascular blood flow velocity mapping by combining dynamic light scattering optical coherence tomography and two-photon microscopy",
            "authors": [
                "Qi Pian", "Mohammed Alfadhel", "Jianbo Tang", "Grace V. Lee",
                "Baoqiang Li", "Buyin Fu", "Yagmur Ayata", "Mohammad Abbas Yaseen",
                "David A. Boas", "Timothy W. Secomb", "Sava Sakadžić",
            ],
            "year": 2023,
            "doi": "10.1117/1.JBO.28.7.076003",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10362155/",
            "exact_value_used": {
                "condition": "average capillary center-line velocity across three mice",
                "mean_mm_s": 1.0,
                "between_mouse_sd_mm_s": 0.2,
                "individual_mouse_means_mm_s": [1.1, 1.0, 0.8],
            },
            "role": "experimental center-line blood-velocity reference",
        },
        {
            "title": "Simulation of oxygen transport and estimation of tissue perfusion in extensive microvascular networks: Application to cerebral cortex",
            "authors": [
                "Jose T. Celaya-Alcala", "Grace V. Lee", "Amy F. Smith", "Bohan Li",
                "Sava Sakadžić", "David A. Boas", "Timothy W. Secomb",
            ],
            "year": 2021,
            "doi": "10.1177/0271678X20927100",
            "url": "https://pubmed.ncbi.nlm.nih.gov/32501155/",
            "exact_value_used": {
                "mouse_cortical_network_segment_flow_nl_min": {"minimum": "<0.1", "maximum": ">100"},
                "nodal_pressure_distribution_peak_mmHg": 31.0,
                "pressure_dependent_target_wss_range_dyn_cm2": [17.0, 97.0],
            },
            "role": "experimentally constrained computational network plausibility reference",
        },
    ]
    required = ("title", "authors", "year", "doi", "url", "exact_value_used")
    checks = {
        "five_primary_articles": len(studies) == 5,
        "all_required_metadata_present": all(all(item.get(key) for key in required) for item in studies),
        "rbc_and_continuum_not_equated": True,
        "network_values_not_misrepresented_as_direct_project_measurements": True,
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "studies": studies,
        "healthy_rbc_reference_range_mm_s": [0.6, 1.0],
        "continuum_mapping": {
            "status": "MODEL_CALIBRATION_ASSUMPTION",
            "assumption": "Newtonian Poiseuille proxy U_center approximately equals 2*U_mean",
            "selected_mean_fluid_velocity_mm_s": 0.35,
            "corresponding_centerline_proxy_mm_s": 0.70,
            "not_present_in_solver": [
                "discrete RBCs", "hematocrit transport", "Fahraeus effect", "Fahraeus-Lindqvist model"
            ],
        },
    }


def calculate_healthy_targets() -> dict[str, Any]:
    diameter = math.sqrt(4.0 * SMOOTH_INLET_AREA_M2 / math.pi)
    dynamic_viscosity = REFERENCE_DENSITY_KG_M3 * NU_M2_S
    candidates: dict[str, Any] = {}
    for label, velocity in TARGET_BRACKET_M_S.items():
        q = velocity * SMOOTH_INLET_AREA_M2
        mass = REFERENCE_DENSITY_KG_M3 * q
        tau_pa = 8.0 * dynamic_viscosity * velocity / diameter
        candidates[label] = {
            "mean_velocity_m_s": velocity,
            "mean_velocity_mm_s": velocity * 1.0e3,
            "volume_flow_m3_s": q,
            "mass_flow_kg_s": mass,
            "volume_flow_nl_min": q * 1.0e12 * 60.0,
            "centerline_proxy_mm_s": 2.0 * velocity * 1.0e3,
            "wss_proxy_pa": tau_pa,
            "wss_proxy_dyn_cm2": tau_pa * 10.0,
            "reynolds_number": velocity * diameter / NU_M2_S,
        }
    return {
        "status": "PASS",
        "density_kg_m3": REFERENCE_DENSITY_KG_M3,
        "kinematic_viscosity_m2_s": NU_M2_S,
        "dynamic_viscosity_pa_s": dynamic_viscosity,
        "smooth_inlet_area_m2": SMOOTH_INLET_AREA_M2,
        "equivalent_inlet_diameter_m": diameter,
        "discrete_inlet_cell_count": EXPECTED_INLET_GLOBBC,
        "discrete_cell_area_proxy_m2": EXPECTED_INLET_GLOBBC * EXPECTED_DX_M**2,
        "discrete_cell_area_proxy_is_physical_area": False,
        "candidates": candidates,
        "selected": "PRIMARY",
        "selected_target": candidates["PRIMARY"],
        "selected_to_old_target_ratio": OLD_TARGET_SCALE_FACTOR,
        "solver_sweep_count": 1,
    }


def generate_healthy_calibration_lua(
    *, mesh_wsl: str, resume_header_wsl: str, maximum_iterations: int = MAXIMUM_ITERATIONS
) -> str:
    checkpoint_start = RESUME_ITERATION + CHECKPOINT_INTERVAL
    return f"""-- Isolated healthy mouse cortical-capillary adaptive-flux calibration.
simulation_name = '{SHORT_SIMULATION_NAME}'
printRuntimeInfo = true
timing_file = 'timing/timing.res'
mesh = '{mesh_wsl}/'
scaling = 'diffusive'
logging = {{ level = 5 }}

dx = {EXPECTED_DX_M:.17g}
dt = {EXPECTED_DT_S:.17g}
rho0_phy = {REFERENCE_DENSITY_KG_M3:.17g}
nu_phy = {NU_M2_S:.17g}
bulk_viscosity_phy = {BULK_NU_M2_S:.17g}
pressure_reference_phy = {PRESSURE_REFERENCE_PA:.17g}
maximum_iterations = {int(maximum_iterations)}

function outlet_01_pressure(x,y,z,t) return {PRESSURE_REFERENCE_PA + OUTLET_GAUGE_PRESSURE_PA['outlet_01']:.17g} end
function outlet_02_pressure(x,y,z,t) return {PRESSURE_REFERENCE_PA + OUTLET_GAUGE_PRESSURE_PA['outlet_02']:.17g} end
function outlet_03_pressure(x,y,z,t) return {PRESSURE_REFERENCE_PA + OUTLET_GAUGE_PRESSURE_PA['outlet_03']:.17g} end

sim_control = {{
  time_control = {{
    max = {{ iter = maximum_iterations, clock = {WALLCLOCK_LIMIT_S} }},
    interval = {{ iter = {CONVERGENCE_INTERVAL} }}
  }},
  abort_criteria = {{
    stop_file = 'stop',
    steady_state = true,
    convergence = {{
      variable = {{ 'pressure_phy', 'vel_mag_phy' }},
      shape = {{ kind = 'all' }},
      reduction = {{ 'average', 'average' }},
      time_control = {{ min = {{ iter = {RESUME_ITERATION} }}, max = {{ iter = maximum_iterations }}, interval = {{ iter = {CONVERGENCE_INTERVAL} }} }},
      norm = 'average',
      nvals = {CONVERGENCE_NVALS},
      absolute = true,
      condition = {{
        {{ threshold = {PRESSURE_THRESHOLD_PA:.17g}, operator = '<=' }},
        {{ threshold = {VELOCITY_THRESHOLD_M_S:.17g}, operator = '<=' }}
      }}
    }}
  }}
}}

physics = {{ dt = dt, rho0 = rho0_phy }}
identify = {{ label = 'ROI003274', kind = 'fluid', layout = 'd3q19', relaxation = 'bgk' }}
fluid = {{ kinematic_viscosity = nu_phy, bulk_viscosity = bulk_viscosity_phy }}
initial_condition = {{ pressure = pressure_reference_phy, velocityX = 0.0, velocityY = 0.0, velocityZ = 0.0 }}

boundary_condition = {{
  {{ label = 'wall', kind = 'wall_libb' }},
  {{ label = 'inlet', kind = 'adaptive_flux_pressure', mass_flowrate = {TARGET_MASS_FLOW_KG_S:.17g} }},
  {{ label = 'outlet_01', kind = 'pressure_eq', pressure = outlet_01_pressure }},
  {{ label = 'outlet_02', kind = 'pressure_eq', pressure = outlet_02_pressure }},
  {{ label = 'outlet_03', kind = 'pressure_eq', pressure = outlet_03_pressure }}
}}

tracking = {{
  {{
    label = 'p', folder = 'tracking/p/',
    variable = {{ 'pressure_phy' }}, shape = {{ kind = 'all' }}, reduction = {{ 'average' }},
    time_control = {{ min = {{ iter = {RESUME_ITERATION} }}, max = {{ iter = maximum_iterations }}, interval = {{ iter = {CONVERGENCE_INTERVAL} }} }},
    output = {{ format = 'ascii' }}
  }},
  {{
    label = 'u', folder = 'tracking/u/',
    variable = {{ 'vel_mag_phy' }}, shape = {{ kind = 'all' }}, reduction = {{ 'average' }},
    time_control = {{ min = {{ iter = {RESUME_ITERATION} }}, max = {{ iter = maximum_iterations }}, interval = {{ iter = {CONVERGENCE_INTERVAL} }} }},
    output = {{ format = 'ascii' }}
  }}
}}

restart = {{
  read = '{resume_header_wsl}',
  write = 'restart/',
  time_control = {{
    min = {{ iter = {checkpoint_start} }},
    max = {{ iter = maximum_iterations }},
    interval = {{ iter = {CHECKPOINT_INTERVAL} }}
  }}
}}
"""


def healthy_lua_contract(text: str, *, resume_header_wsl: str) -> dict[str, Any]:
    checks = {
        "restart_read": f"read = '{resume_header_wsl}'" in text,
        "not_fresh": "read =" in text,
        "healthy_mass_target": f"mass_flowrate = {TARGET_MASS_FLOW_KG_S:.17g}" in text,
        "fixed_pressure_reference": f"pressure_reference_phy = {PRESSURE_REFERENCE_PA:.17g}" in text,
        "fixed_outlet_pressures": all(
            f"return {PRESSURE_REFERENCE_PA + gauge:.17g}" in text
            for gauge in OUTLET_GAUGE_PRESSURE_PA.values()
        ),
        "fixed_physics": all(
            token in text
            for token in (
                f"dx = {EXPECTED_DX_M:.17g}", f"dt = {EXPECTED_DT_S:.17g}",
                f"rho0_phy = {REFERENCE_DENSITY_KG_M3:.17g}", f"nu_phy = {NU_M2_S:.17g}",
                f"bulk_viscosity_phy = {BULK_NU_M2_S:.17g}", "layout = 'd3q19'",
                "relaxation = 'bgk'", "kind = 'wall_libb'",
            )
        ),
        "adaptive_flux": "kind = 'adaptive_flux_pressure'" in text,
        "three_pressure_eq_outlets": text.count("kind = 'pressure_eq'") == 3,
        "criterion_interval": f"interval = {{ iter = {CONVERGENCE_INTERVAL} }}" in text,
        "criterion_nvals": f"nvals = {CONVERGENCE_NVALS}" in text,
        "criterion_min_resume": f"min = {{ iter = {RESUME_ITERATION} }}" in text,
        "pressure_threshold": f"threshold = {PRESSURE_THRESHOLD_PA:.17g}" in text,
        "velocity_threshold": f"threshold = {VELOCITY_THRESHOLD_M_S:.17g}" in text,
        "checkpoint_interval": f"interval = {{ iter = {CHECKPOINT_INTERVAL} }}" in text,
        "wallclock": f"clock = {WALLCLOCK_LIMIT_S}" in text,
        "no_harvester": "harvest" not in text.lower(),
        # The accepted frozen mesh path legitimately contains a ``seeder``
        # directory name.  Guard executable/config invocation, not provenance.
        "no_seeder_call": "seeder.lua" not in text.lower() and "seed.lua" not in text.lower(),
    }
    return {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}


def _pressure_neighbor_indices(
    cell_indices: np.ndarray,
    cell_ijk: np.ndarray,
    lookup: dict[tuple[int, int, int], int],
    direction: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    normal = np.asarray(direction, dtype=np.int64).reshape(3)
    neighbors1 = np.full(len(cell_indices), -1, dtype=np.int64)
    neighbors2 = np.full(len(cell_indices), -1, dtype=np.int64)
    for row, cell_index in enumerate(np.asarray(cell_indices, dtype=np.int64)):
        coordinate = cell_ijk[cell_index]
        neighbors1[row] = lookup.get(tuple(int(v) for v in coordinate + normal), -1)
        neighbors2[row] = lookup.get(tuple(int(v) for v in coordinate + 2 * normal), -1)
    valid = (neighbors1 >= 0) & (neighbors2 >= 0)
    return valid, neighbors1, neighbors2


def _mesh_boundary_context(mesh: Path) -> dict[str, Any]:
    property_header = parse_boundary_property_header((mesh / "header.lua").read_text(encoding="utf-8"))
    boundary_header = parse_bnd_header((mesh / "bnd.lua").read_text(encoding="utf-8"))
    tree_ids, property_bits, elemlist_contract = read_treelm_elemlist(
        mesh / "elemlist.lsb", n_elems=EXPECTED_FLUID_CELLS
    )
    levels = tree_levels(tree_ids)
    if not np.all(levels == 9):
        raise FlowError(STATUS_FAILED, "Frozen adaptive mesh is not uniform root level 9")
    cell_ijk = tree_ids_to_ijk(tree_ids, levels)
    property_indices = extract_boundary_property_indices(property_bits, property_header.bit_position)
    boundary_ids = read_boundary_ids(
        mesh / "bnd.lsb", element_count=property_header.element_count,
        side_count=boundary_header.side_count,
    )
    label_to_id = {label: index for index, label in enumerate(boundary_header.labels, start=1)}
    boundaries = {
        label: reconstruct_boundary(
            boundary_ids, property_indices, label=label, boundary_id=label_to_id[label]
        )
        for label in PORT_LABELS
    }
    return {
        "tree_ids": tree_ids,
        "cell_ijk": cell_ijk,
        "lookup": build_coordinate_lookup(cell_ijk),
        "boundaries": boundaries,
        "elemlist_contract": elemlist_contract,
    }


def source_only_pressure_outlet_contract(mesh: Path) -> dict[str, Any]:
    context = _mesh_boundary_context(mesh)
    records: dict[str, Any] = {}
    checks: list[bool] = []
    for label, direction in CORRECTED_RUNTIME_NORMALS.items():
        boundary = context["boundaries"][label]
        valid, _, _ = _pressure_neighbor_indices(
            boundary.cell_indices, context["cell_ijk"], context["lookup"], direction
        )
        observed = (len(boundary.cell_indices), int(np.count_nonzero(valid)), int(np.count_nonzero(~valid)))
        expected = EXPECTED_PRESSURE_COUNTS[label]
        records[label] = {
            "globbc_total": observed[0], "valid_pressure_cells": observed[1],
            "removed_cells": observed[2], "runtime_normal": direction.astype(int).tolist(),
            "curved": False, "neighbors": ["cell+n", "cell+2n"],
            "expected": list(expected), "matches_expected": observed == expected,
        }
        checks.append(observed == expected)
    return {"status": "PASS" if all(checks) else "FAIL", "outlets": records}


def _parse_pressure_tracking(paths: Iterable[Path]) -> dict[int, float]:
    records: dict[int, float] = {}
    for path in paths:
        for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) < 2:
                continue
            time_s, pressure_pa = float(fields[0]), float(fields[1])
            iteration = int(round(time_s / EXPECTED_DT_S))
            if iteration in records and not math.isclose(records[iteration], pressure_pa, rel_tol=0.0, abs_tol=1e-8):
                raise FlowError(STATUS_FAILED, f"Conflicting pressure tracking at iteration {iteration}")
            records[iteration] = pressure_pa
    if not records:
        raise FlowError(STATUS_FAILED, "No global-average pressure tracking records found")
    return records


def global_mass_accumulation_audit(
    *, start_binary: Path, end_binary: Path, start_iteration: int,
    end_iteration: int, pressure_records: dict[int, float],
) -> dict[str, Any]:
    difference = end_iteration - start_iteration
    if difference < EFFECTIVE_CONVERGENCE_WINDOW_ITERATIONS:
        raise FlowError(STATUS_FAILED, f"Mass window is only {difference} iterations")
    start_pdf = read_restart_pdf(start_binary, n_elems=EXPECTED_FLUID_CELLS, n_components=19)
    end_pdf = read_restart_pdf(end_binary, n_elems=EXPECTED_FLUID_CELLS, n_components=19)
    start_mass = float(np.sum(start_pdf, dtype=np.float64))
    end_mass = float(np.sum(end_pdf, dtype=np.float64))
    delta = end_mass - start_mass
    factor = REFERENCE_DENSITY_KG_M3 * EXPECTED_DX_M**3 / EXPECTED_DT_S
    accumulation = (delta / difference) * factor
    if start_iteration not in pressure_records or end_iteration not in pressure_records:
        missing = [i for i in (start_iteration, end_iteration) if i not in pressure_records]
        raise FlowError(STATUS_FAILED, f"Pressure tracking missing mass-window iterations: {missing}")
    p_start = pressure_records[start_iteration]
    p_end = pressure_records[end_iteration]
    lattice_pressure_scale = REFERENCE_DENSITY_KG_M3 * EXPECTED_DX_M**2 / EXPECTED_DT_S**2 * CS2
    pressure_delta_sum = EXPECTED_FLUID_CELLS * (p_end - p_start) / lattice_pressure_scale
    pressure_accumulation = (pressure_delta_sum / difference) * factor
    direction_same = (
        accumulation == 0.0 or pressure_accumulation == 0.0
        or math.copysign(1.0, accumulation) == math.copysign(1.0, pressure_accumulation)
    )
    discrepancy = abs(accumulation - pressure_accumulation) / max(
        abs(accumulation), abs(pressure_accumulation), 1.0e-30
    )
    ratio = abs(accumulation) / TARGET_MASS_FLOW_KG_S
    return {
        "status": "PASS" if ratio <= MASS_ACCUMULATION_HARD_GATE and direction_same and discrepancy <= PRESSURE_CROSSCHECK_GATE else "FAIL",
        "start_iteration": start_iteration, "end_iteration": end_iteration,
        "iteration_difference": difference, "start_binary": str(start_binary),
        "end_binary": str(end_binary), "start_total_pdf_mass": start_mass,
        "end_total_pdf_mass": end_mass, "total_pdf_mass_change": delta,
        "mass_accumulation_kg_s": accumulation,
        "accumulation_to_healthy_target_ratio": ratio,
        "hard_gate": MASS_ACCUMULATION_HARD_GATE,
        "preferred_gate": MASS_ACCUMULATION_PREFERRED_GATE,
        "pressure_start_pa": p_start, "pressure_end_pa": p_end,
        "pressure_derived_total_pdf_mass_change": pressure_delta_sum,
        "pressure_derived_accumulation_kg_s": pressure_accumulation,
        "pdf_pressure_direction_same": direction_same,
        "pdf_pressure_relative_discrepancy": discrepancy,
        "pdf_pressure_crosscheck_gate": PRESSURE_CROSSCHECK_GATE,
        "pdf_pressure_crosscheck_status": "PASS" if direction_same and discrepancy <= PRESSURE_CROSSCHECK_GATE else "FAIL",
    }


def exact_boundary_flux_audit(
    *, mesh: Path, restart_binary: Path, link_csv: Path
) -> dict[str, Any]:
    context = _mesh_boundary_context(mesh)
    boundaries = context["boundaries"]
    tree_ids = context["tree_ids"]
    cell_ijk = context["cell_ijk"]
    lookup = context["lookup"]
    if len(boundaries["inlet"].cell_indices) != EXPECTED_INLET_GLOBBC:
        raise FlowError(STATUS_FAILED, "Final inlet globBC count is not 287")
    pdf = read_restart_pdf(restart_binary, n_elems=EXPECTED_FLUID_CELLS, n_components=19)
    values = np.asarray(pdf)
    all_finite = bool(np.all(np.isfinite(values)))
    minimum_pdf = float(np.min(values))
    density = np.sum(values, axis=1, dtype=np.float64)
    if not all_finite or np.any(density <= 0.0):
        raise FlowError(STATUS_SAFETY_FAILED, "Final restart PDFs are non-finite or have non-positive density")
    velocity = (values @ D3Q19_DIRECTIONS.astype(np.float64)) / density[:, None]
    maximum_lattice_speed = float(np.max(np.linalg.norm(velocity, axis=1)))

    inlet = boundaries["inlet"]
    valid, neighbor1, neighbor2 = _pressure_neighbor_indices(
        inlet.cell_indices, cell_ijk, lookup, INLET_PRESSURE_NEIGHBOR_DIRECTION
    )
    if not np.all(valid):
        raise FlowError(STATUS_FAILED, f"Adaptive inlet lacks two neighbors for {np.count_nonzero(~valid)} cells")
    fetched1 = pull_fetch_pdfs(pdf, cell_ijk, neighbor1, coordinate_lookup=lookup)
    fetched2 = pull_fetch_pdfs(pdf, cell_ijk, neighbor2, coordinate_lookup=lookup)
    extrapolated_velocity = 1.5 * velocity_from_pdf(fetched1) - 0.5 * velocity_from_pdf(fetched2)
    stored = np.asarray(pdf[inlet.cell_indices], dtype=np.float64)
    alpha, beta = musubi_pressure_flux_affine_coefficients(
        stored_boundary_pdfs=stored, incoming_masks=inlet.incoming_masks,
        extrapolated_velocity=extrapolated_velocity,
    )
    target_flux_lattice = physical_volume_flux_to_lattice(
        TARGET_Q_M3_S, dx_m=EXPECTED_DX_M, dt_s=EXPECTED_DT_S
    )
    rho_boundary = solve_boundary_density(target_flux_lattice, alpha, beta)
    fetched_boundary = pull_fetch_pdfs(pdf, cell_ijk, inlet.cell_indices, coordinate_lookup=lookup)
    reconstructed_inlet = reconstruct_musubi_boundary_state(
        rho_boundary, fetched_boundary_pdfs=fetched_boundary,
        incoming_masks=inlet.incoming_masks, extrapolated_velocity=extrapolated_velocity,
    )
    mass_factor = physical_mass_factor(
        density_kg_m3=REFERENCE_DENSITY_KG_M3, dx_m=EXPECTED_DX_M, dt_s=EXPECTED_DT_S
    )
    link_records, inlet_mass = _link_rows(
        boundary=inlet, selected_rows=np.arange(len(inlet.cell_indices), dtype=np.int64),
        tree_ids=tree_ids, old_pdf=pdf, new_pdf=reconstructed_inlet,
        mass_factor=mass_factor, outward_sign=False,
    )
    inlet_error = abs(inlet_mass - TARGET_MASS_FLOW_KG_S) / TARGET_MASS_FLOW_KG_S
    pressure_scale = REFERENCE_DENSITY_KG_M3 * EXPECTED_DX_M**2 / EXPECTED_DT_S**2 * CS2
    inlet_absolute_pressure = rho_boundary * pressure_scale
    inlet_gauge_pressure = inlet_absolute_pressure - pressure_scale

    outlets: dict[str, Any] = {}
    outlet_mass: dict[str, float] = {}
    for label, direction in CORRECTED_RUNTIME_NORMALS.items():
        boundary = boundaries[label]
        valid, neighbor1, neighbor2 = _pressure_neighbor_indices(
            boundary.cell_indices, cell_ijk, lookup, direction
        )
        selected = np.flatnonzero(valid).astype(np.int64)
        observed = (len(boundary.cell_indices), len(selected), int(np.count_nonzero(~valid)))
        if observed != EXPECTED_PRESSURE_COUNTS[label]:
            raise FlowError(STATUS_FAILED, f"{label} pressure cell contract changed: {observed}")
        fetched1 = pull_fetch_pdfs(pdf, cell_ijk, neighbor1[selected], coordinate_lookup=lookup)
        fetched2 = pull_fetch_pdfs(pdf, cell_ijk, neighbor2[selected], coordinate_lookup=lookup)
        outlet_velocity = 1.5 * velocity_from_pdf(fetched1) - 0.5 * velocity_from_pdf(fetched2)
        absolute_pressure = PRESSURE_REFERENCE_PA + OUTLET_GAUGE_PRESSURE_PA[label]
        replacement = equilibrium_pdf(absolute_pressure / pressure_scale, outlet_velocity)
        records, mass = _link_rows(
            boundary=boundary, selected_rows=selected, tree_ids=tree_ids,
            old_pdf=pdf, new_pdf=replacement, mass_factor=mass_factor, outward_sign=True,
        )
        link_records.extend(records)
        outlet_mass[label] = mass
        outlets[label] = {
            "globbc_total": observed[0], "valid_pressure_cells": observed[1],
            "removed_cells": observed[2], "runtime_normal": direction.astype(int).tolist(),
            "curved": False, "neighbor_semantics": ["cell+n", "cell+2n"],
            "velocity_semantics": "1.5*u1 - 0.5*u2", "gauge_pressure_pa": OUTLET_GAUGE_PRESSURE_PA[label],
            "absolute_pressure_pa": absolute_pressure, "exact_signed_mass_flow_kg_s": mass,
        }
    _write_link_csv(link_csv, link_records)
    balance = signed_mass_balance(inlet_mass, outlet_mass.values())
    output2 = outlet_mass["outlet_02"]
    significant_backflow = output2 < 0.0 and abs(output2) > OUTLET_02_BACKFLOW_FRACTION_GATE * TARGET_MASS_FLOW_KG_S
    exact_q = inlet_mass / REFERENCE_DENSITY_KG_M3
    umean = exact_q / SMOOTH_INLET_AREA_M2
    diameter = math.sqrt(4.0 * SMOOTH_INLET_AREA_M2 / math.pi)
    mu = REFERENCE_DENSITY_KG_M3 * NU_M2_S
    wss_pa = 8.0 * mu * umean / diameter
    return {
        "status": "PASS" if (
            inlet_error <= INLET_TARGET_HARD_GATE
            and balance["relative_error"] <= BOUNDARY_BALANCE_HARD_GATE
            and not significant_backflow and all_finite
            and minimum_pdf > 0.0 and maximum_lattice_speed < 0.05
        ) else "FAIL",
        "method": "final-restart PULL/FETCH incoming-link replacement audit",
        "controller_used_as_flux_source": False,
        "inlet": {
            "globbc_count": len(inlet.cell_indices),
            "pressure_neighbor_direction": INLET_PRESSURE_NEIGHBOR_DIRECTION.astype(int).tolist(),
            "target_mass_flow_kg_s": TARGET_MASS_FLOW_KG_S,
            "exact_mass_flow_kg_s": inlet_mass, "target_relative_error": inlet_error,
            "preferred_gate": INLET_TARGET_PREFERRED_GATE, "hard_gate": INLET_TARGET_HARD_GATE,
            "alpha": alpha, "beta": beta, "rho_boundary": rho_boundary,
            "absolute_pressure_pa": inlet_absolute_pressure, "gauge_pressure_pa": inlet_gauge_pressure,
        },
        "outlets": outlets, "outlet_signed_sum_kg_s": balance["outlet_signed_sum_kg_s"],
        "instantaneous_mass_balance_relative_error": balance["relative_error"],
        "mass_balance_preferred_gate": BOUNDARY_BALANCE_PREFERRED_GATE,
        "mass_balance_hard_gate": BOUNDARY_BALANCE_HARD_GATE,
        "significant_outlet_02_backflow": significant_backflow,
        "outlet_02_to_healthy_target_ratio": output2 / TARGET_MASS_FLOW_KG_S,
        "all_finite": all_finite, "minimum_pdf": minimum_pdf,
        "maximum_lattice_speed": maximum_lattice_speed,
        "exact_volume_flow_m3_s": exact_q, "smooth_area_mean_velocity_m_s": umean,
        "centerline_proxy_mm_s": 2.0 * umean * 1.0e3,
        "wss_proxy_pa": wss_pa, "wss_proxy_dyn_cm2": wss_pa * 10.0,
        "reynolds_number": umean * diameter / NU_M2_S,
        "healthy_velocity_compatibility": "PASS" if HEALTHY_UMEAN_RANGE_M_S[0] <= umean <= HEALTHY_UMEAN_RANGE_M_S[1] else "FAIL",
        "healthy_wss_proxy_compatibility": "PASS" if HEALTHY_WSS_RANGE_DYN_CM2[0] <= wss_pa * 10.0 <= HEALTHY_WSS_RANGE_DYN_CM2[1] else "FAIL",
        "per_link_csv": str(link_csv), "elemlist_contract": context["elemlist_contract"],
    }


def classify_healthy_calibration(
    *, official_steady: bool, additional_iterations: int,
    exact: dict[str, Any], mass: dict[str, Any], source_unchanged: bool,
) -> tuple[str, str, str | None]:
    safety_failures = []
    if not exact["all_finite"]:
        safety_failures.append("non-finite final PDFs")
    if exact["minimum_pdf"] <= 0.0:
        safety_failures.append("minimum PDF <= 0")
    if exact["maximum_lattice_speed"] >= 0.05:
        safety_failures.append("maximum lattice speed >= 0.05")
    if safety_failures:
        return STATUS_SAFETY_FAILED, NEXT_SAFETY, "; ".join(safety_failures)
    if not source_unchanged:
        return STATUS_FAILED, NEXT_REVIEW, "source frozen files changed"
    if exact["inlet"]["target_relative_error"] > INLET_TARGET_HARD_GATE:
        return STATUS_FAILED, NEXT_REVIEW, "healthy adaptive inlet target error > 1%"
    if exact["healthy_velocity_compatibility"] != "PASS":
        return STATUS_FAILED, NEXT_REVIEW, "smooth-area mean velocity outside 0.30-0.40 mm/s"
    mass_ok = (
        official_steady and additional_iterations >= EFFECTIVE_CONVERGENCE_WINDOW_ITERATIONS
        and mass["accumulation_to_healthy_target_ratio"] <= MASS_ACCUMULATION_HARD_GATE
        and mass["pdf_pressure_crosscheck_status"] == "PASS"
    )
    if not mass_ok:
        reason = (
            "p/u official steady not reached" if not official_steady
            else "global mass accumulation/crosscheck gate failed"
        )
        return STATUS_MASS_INCOMPLETE, NEXT_RESUME, reason
    boundary_ok = exact["instantaneous_mass_balance_relative_error"] <= BOUNDARY_BALANCE_HARD_GATE
    if not boundary_ok:
        return STATUS_OUTLETS_UNRESOLVED, NEXT_OUTLETS, "instantaneous exact boundary balance > 1%"
    if exact["significant_outlet_02_backflow"]:
        return STATUS_OUTLETS_UNRESOLVED, NEXT_OUTLETS, "significant outlet_02 backflow"
    if exact["healthy_wss_proxy_compatibility"] != "PASS":
        return STATUS_FAILED, NEXT_REVIEW, "Poiseuille WSS proxy outside 17-97 dyn/cm2"
    return STATUS_PASS, NEXT_GRID, None


def _archive_checkpoint_manifest(
    tracker: HealthyCheckpointTracker, runtime_root_wsl: str, run_root: Path
) -> None:
    for checkpoint in tracker.records:
        hrel = PurePosixPath(checkpoint["header_path_wsl"]).relative_to(PurePosixPath(runtime_root_wsl))
        brel = PurePosixPath(checkpoint["binary_path_wsl"]).relative_to(PurePosixPath(runtime_root_wsl))
        ah = run_root.joinpath(*hrel.parts)
        ab = run_root.joinpath(*brel.parts)
        archive_sha = sha256_file(ab)
        checkpoint.update({
            "archived_header_path": str(ah), "archived_binary_path": str(ab),
            "archived_binary_sha256": archive_sha,
            "archive_sha256_match": archive_sha == checkpoint["sha256"],
        })
        if not checkpoint["archive_sha256_match"]:
            raise FlowError(STATUS_FAILED, f"Archived checkpoint SHA mismatch at {checkpoint['iteration']}")
    tracker._write()


def _select_mass_window(
    *, final_iteration: int, final_binary: Path, source_binary: Path,
    checkpoints: list[dict[str, Any]],
) -> tuple[int, Path, int, Path]:
    candidates = [(RESUME_ITERATION, source_binary)]
    candidates.extend(
        (int(item["iteration"]), Path(item["archived_binary_path"]))
        for item in checkpoints
        if int(item["iteration"]) <= final_iteration - EFFECTIVE_CONVERGENCE_WINDOW_ITERATIONS
    )
    usable = [item for item in candidates if item[0] <= final_iteration - EFFECTIVE_CONVERGENCE_WINDOW_ITERATIONS]
    if not usable:
        raise FlowError(STATUS_FAILED, "No complete restart at least 10000 iterations before final")
    start_iteration, start_binary = max(usable, key=lambda item: item[0])
    return start_iteration, start_binary, final_iteration, final_binary


def run_healthy_capillary_calibration(project_root: Path) -> dict[str, Any]:
    """Run exactly one isolated PRIMARY healthy-target Musubi continuation."""

    root = Path(project_root).resolve()
    head = _git(root, "rev-parse", "HEAD")
    branch = _git(root, "branch", "--show-current")
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    output_root = config.paths.output_root
    mesh = output_root / AXIS_MESH_RUN / "seeder" / "mesh"
    inlet_rim = output_root / AXIS_MESH_RUN / "qc" / "inlet_rim_audit.json"
    source_root = output_root / SOURCE_RUN
    source_header = source_root / "restart" / f"{SHORT_SIMULATION_NAME}_lastHeader.lua"
    source_binary = source_root / "restart" / RESUME_BINARY_NAME
    binary_windows = Path(r"\\wsl.localhost\Ubuntu\home\lzy\apes-worktrees\musubi_mcclure_adaptive_flux_20260829_1300\build\musubi_adaptive_flux")

    production_paths = (
        root / "cfd_flow.py", root / "configs" / "cfd_flow.yaml",
        root / "utils" / "cfd_flow" / "pipeline.py",
    )
    production_diff = _git(root, "diff", "--", *(str(p.relative_to(root)) for p in production_paths))
    if production_diff:
        raise FlowError(STATUS_FAILED, "Production pipeline has local modifications")
    if config.mesh.dx_target_m != EXPECTED_DX_M or config.solver.wallclock_limit_s != WALLCLOCK_LIMIT_S:
        raise FlowError(STATUS_FAILED, "Frozen dx or wallclock setting changed")
    if parse_mesh_header(mesh)["fluid_element_count"] != EXPECTED_FLUID_CELLS:
        raise FlowError(STATUS_FAILED, "Frozen mesh fluid cell count changed")
    source = parse_restart_header(source_header)
    source_contract = restart_binary_size_contract(
        source_binary, n_elems=source.n_elems, n_components=source.n_components, n_dofs=source.n_dofs
    )
    source_sha = sha256_file(source_binary)
    binary_sha = sha256_file(binary_windows)
    if source.iteration != RESUME_ITERATION or source.n_elems != EXPECTED_FLUID_CELLS or source.n_components != 19 or source.n_dofs != 1 or source_contract["status"] != "PASS":
        raise FlowError(STATUS_FAILED, f"Healthy source restart contract failed: {source}")
    if source_sha != EXPECTED_RESUME_SHA256:
        raise FlowError(STATUS_FAILED, f"Healthy source restart SHA changed: {source_sha}")
    if binary_sha != EXPECTED_BINARY_SHA256:
        raise FlowError(STATUS_FAILED, f"Pinned adaptive binary SHA changed: {binary_sha}")
    outlet_preflight = source_only_pressure_outlet_contract(mesh)
    if outlet_preflight["status"] != "PASS":
        raise FlowError(STATUS_FAILED, f"Source-only pressure outlet contract failed: {outlet_preflight}")

    frozen_paths = (
        *production_paths, *sorted(path for path in mesh.iterdir() if path.is_file()),
        inlet_rim, source_header, source_binary, binary_windows,
    )
    frozen_before = _file_manifest(frozen_paths)
    local_frozen_before = _file_manifest(frozen_paths[:-1])
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_root = output_root / f"{RUN_PREFIX}_{stamp}"
    runtime_root_wsl = f"{RUNTIME_BASE_WSL}/h_{stamp[-6:]}"
    runtime_root = _wsl_path_to_unc(config.apes.wsl_distribution, runtime_root_wsl)
    qc = run_root / "qc"
    tracking = run_root / "tracking"
    restart = run_root / "restart"
    qc.mkdir(parents=True, exist_ok=False)
    tracking.mkdir(parents=True, exist_ok=False)
    manifest_path = qc / "healthy_calibration_manifest.json"
    literature = healthy_literature_reference_contract()
    targets = calculate_healthy_targets()
    write_json(qc / "healthy_mouse_capillary_reference.json", literature)
    write_json(qc / "healthy_flow_target_calculation.json", targets)
    write_json(qc / "healthy_pressure_outlet_source_preflight.json", outlet_preflight)

    summary: dict[str, Any] = {
        "status": STATUS_FAILED, "next": NEXT_REVIEW, "run_root": str(run_root),
        "actual_head": head, "branch": branch, "production_pipeline_modified": False,
        "resume_iteration": RESUME_ITERATION, "resume_header": str(source_header),
        "resume_binary": str(source_binary), "resume_binary_sha256": source_sha,
        "adaptive_binary_sha256": binary_sha, "adaptive_binary_rebuilt": False,
        "mesh": str(mesh), "fluid_cells": EXPECTED_FLUID_CELLS,
        "seeder_calls": 0, "musubi_calls": 0, "harvester_calls": 0,
        "mpi_ranks": MPI_RANKS, "grid_convergence": "NOT_RUN",
        "dx_m": EXPECTED_DX_M, "dt_s": EXPECTED_DT_S,
        "target_q_m3_s": TARGET_Q_M3_S, "target_mass_flow_kg_s": TARGET_MASS_FLOW_KG_S,
        "outlet_gauge_pressures_pa": OUTLET_GAUGE_PRESSURE_PA,
        "pressure_reference_pa": PRESSURE_REFERENCE_PA,
        "convergence": {"interval": CONVERGENCE_INTERVAL, "nvals": CONVERGENCE_NVALS,
                        "effective_window": EFFECTIVE_CONVERGENCE_WINDOW_ITERATIONS,
                        "velocity_threshold_m_s": VELOCITY_THRESHOLD_M_S,
                        "pressure_threshold_pa": PRESSURE_THRESHOLD_PA,
                        "earliest_final_steady_iteration": EARLIEST_FINAL_STEADY_ITERATION},
        "checkpoint_interval": CHECKPOINT_INTERVAL, "wallclock_limit_s": WALLCLOCK_LIMIT_S,
        "runtime_root_wsl": runtime_root_wsl, "frozen_files_before": frozen_before,
        "literature_reference_contract": literature["status"], "started_at": datetime.now().isoformat(),
    }
    write_json(manifest_path, summary)
    tracker: HealthyCheckpointTracker | None = None
    try:
        health = _one_wsl_health_check(config.apes.wsl_distribution)
        write_json(qc / "wsl_health_preflight.json", health)
        if health["status"] != "PASS":
            raise FlowError(STATUS_FAILED, "WSL health preflight failed")
        if runtime_root.exists():
            raise FlowError(STATUS_FAILED, f"Short WSL runtime root already exists: {runtime_root_wsl}")
        _prepare_short_tracking_root(distribution=config.apes.wsl_distribution, root_wsl=runtime_root_wsl)
        staged_header = runtime_root / "restart" / source_header.name
        staged_binary = runtime_root / "restart" / source_binary.name
        shutil.copy2(source_header, staged_header)
        shutil.copy2(source_binary, staged_binary)
        if sha256_file(staged_binary) != source_sha or sha256_file(staged_header) != sha256_file(source_header):
            raise FlowError(STATUS_FAILED, "Staged source checkpoint SHA mismatch")
        resume_header_wsl = f"restart/{source_header.name}"
        lua = generate_healthy_calibration_lua(
            mesh_wsl=windows_to_wsl(mesh, config.apes.wsl_distribution),
            resume_header_wsl=resume_header_wsl,
        )
        lua_path = runtime_root / "musubi.lua"
        lua_path.write_text(lua, encoding="utf-8")
        lua_contract = healthy_lua_contract(lua, resume_header_wsl=resume_header_wsl)
        write_json(qc / "healthy_lua_contract.json", lua_contract)
        if lua_contract["status"] != "PASS":
            raise FlowError(STATUS_FAILED, f"Healthy Lua contract failed: {lua_contract}")
        luac_rc = _run_luac(
            distribution=config.apes.wsl_distribution, workdir_wsl=runtime_root_wsl,
            stdout_path=tracking / "luac_stdout.log", stderr_path=tracking / "luac_stderr.log",
        )
        preflight_checks = {
            "literature_contract": literature["status"] == "PASS",
            "target_calculation": targets["status"] == "PASS",
            "source_restart": source_contract["status"] == "PASS",
            "source_restart_sha": source_sha == EXPECTED_RESUME_SHA256,
            "adaptive_binary_sha": binary_sha == EXPECTED_BINARY_SHA256,
            "outlet_source_contract": outlet_preflight["status"] == "PASS",
            "lua_contract": lua_contract["status"] == "PASS", "lua_syntax": luac_rc == 0,
            "persistent_wsl_ext4": runtime_root_wsl.startswith(f"{RUNTIME_BASE_WSL}/"),
            "maximum_path_length_at_most_80": _predicted_maximum_runtime_path_length(runtime_root_wsl) <= 80,
            "production_pipeline_unchanged": not production_diff,
        }
        preflight = {"status": "PASS" if all(preflight_checks.values()) else "FAIL", "checks": preflight_checks}
        write_json(qc / "healthy_static_preflight.json", preflight)
        if preflight["status"] != "PASS":
            raise FlowError(STATUS_FAILED, f"Healthy static preflight failed: {preflight}")

        tracker = HealthyCheckpointTracker(
            runtime_root=runtime_root, runtime_root_wsl=runtime_root_wsl,
            manifest_path=run_root / "checkpoint_manifest.json",
        )
        summary.update({"status": "CFD_FLOW_HEALTHY_CALIBRATION_RUNNING", "musubi_calls": 1})
        write_json(manifest_path, summary)
        run = _run_monitored_wsl(
            distribution=config.apes.wsl_distribution, workdir_wsl=runtime_root_wsl,
            command=["env", "OMPI_ALLOW_RUN_AS_ROOT=1", "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1",
                     MPIRUN_WSL, "-np", str(MPI_RANKS), BINARY_WSL, "musubi.lua"],
            stdout_path=tracking / "musubi_stdout.log", stderr_path=tracking / "musubi_stderr.log",
            controller_csv_path=tracking / "controller_records.csv",
            timeout_s=WALLCLOCK_LIMIT_S + WRAPPER_MARGIN_S, checkpoint_tracker=tracker,
            minimum_controller_iteration_exclusive=RESUME_ITERATION,
        )
        summary.update({
            "musubi_returncode": run.returncode, "musubi_wall_time_s": run.wall_time_s,
            "live_safety_failure": run.safety_failure, "continuity_failure": run.continuity_failure,
            "wsl_infrastructure_failure": run.infrastructure_failure, "wrapper_timeout": run.wrapper_timeout,
            "controller_record_count": run.controller_record_count,
            "last_controller_iteration": run.latest_controller_iteration,
            "checkpoint_iterations": [int(item["iteration"]) for item in tracker.records],
        })
        write_json(manifest_path, summary)
        if run.safety_failure:
            raise FlowError(STATUS_SAFETY_FAILED, run.safety_failure)
        if run.continuity_failure or run.infrastructure_failure or run.wrapper_timeout or run.returncode != 0:
            raise FlowError(STATUS_FAILED, run.continuity_failure or run.infrastructure_failure or ("wrapper timeout" if run.wrapper_timeout else f"Musubi returned {run.returncode}"))
        maximum_path_length = _maximum_runtime_output_path_length(runtime_root_wsl, runtime_root)
        if maximum_path_length > 80:
            raise FlowError(STATUS_FAILED, f"Musubi-visible path length {maximum_path_length} > 80")
        _archive_runtime(runtime_root, run_root)
        _archive_checkpoint_manifest(tracker, runtime_root_wsl, run_root)

        controller = summarize_controller_csv(tracking / "controller_records.csv")
        if controller["first_iteration"] != EXPECTED_FIRST_CONTROLLER_ITERATION:
            raise FlowError(STATUS_FAILED, f"First controller iteration {controller['first_iteration']} != {EXPECTED_FIRST_CONTROLLER_ITERATION}")
        write_json(qc / "healthy_runtime_controller.json", controller)
        combined_log = (tracking / "musubi_stdout.log").read_text(encoding="utf-8", errors="replace") + "\n" + (tracking / "musubi_stderr.log").read_text(encoding="utf-8", errors="replace")
        steady = parse_official_steady_termination(combined_log)
        official = bool(steady["official_steady_termination"])
        final_iteration = int(controller["final_iteration"])
        additional = final_iteration - RESUME_ITERATION
        if official and int(steady["confirmation_iteration"]) != final_iteration:
            raise FlowError(STATUS_FAILED, "Official p/u confirmation does not match final controller iteration")
        if official and final_iteration < EARLIEST_FINAL_STEADY_ITERATION:
            raise FlowError(STATUS_FAILED, "Official p/u steady occurred before the mandatory 10000-iteration window")
        write_json(qc / "healthy_pu_steady_termination.json", steady)

        final_header = restart / f"{SHORT_SIMULATION_NAME}_lastHeader.lua"
        final = parse_restart_header(final_header)
        final_binary = restart / final.binary_path.name
        final_contract = restart_binary_size_contract(
            final_binary, n_elems=final.n_elems, n_components=final.n_components, n_dofs=final.n_dofs
        )
        if final.iteration != final_iteration or final.n_elems != EXPECTED_FLUID_CELLS or final.n_components != 19 or final.n_dofs != 1 or final_contract["status"] != "PASS":
            raise FlowError(STATUS_FAILED, f"Final restart contract failed: {final}")
        pressure_records = _parse_pressure_tracking((tracking / "p").glob("*.res"))
        start_iter, start_bin, end_iter, end_bin = _select_mass_window(
            final_iteration=final_iteration, final_binary=final_binary,
            source_binary=source_binary, checkpoints=tracker.records,
        )
        mass = global_mass_accumulation_audit(
            start_binary=start_bin, end_binary=end_bin, start_iteration=start_iter,
            end_iteration=end_iter, pressure_records=pressure_records,
        )
        write_json(qc / "healthy_global_mass_accumulation.json", mass)
        exact = exact_boundary_flux_audit(
            mesh=mesh, restart_binary=final_binary,
            link_csv=qc / "healthy_exact_boundary_links.csv",
        )
        write_json(qc / "healthy_exact_boundary_flux.json", exact)
        physiology = {
            "status": "PASS" if exact["healthy_velocity_compatibility"] == "PASS" and exact["healthy_wss_proxy_compatibility"] == "PASS" else "FAIL",
            "rbc_reference_range_mm_s": literature["healthy_rbc_reference_range_mm_s"],
            "continuum_mapping": literature["continuum_mapping"],
            "exact_smooth_area_mean_velocity_mm_s": exact["smooth_area_mean_velocity_m_s"] * 1.0e3,
            "poiseuille_centerline_proxy_mm_s": exact["centerline_proxy_mm_s"],
            "poiseuille_wss_proxy_pa": exact["wss_proxy_pa"],
            "poiseuille_wss_proxy_dyn_cm2": exact["wss_proxy_dyn_cm2"],
            "healthy_wss_reference_dyn_cm2": list(HEALTHY_WSS_RANGE_DYN_CM2),
            "reynolds_number": exact["reynolds_number"],
            "actual_cfd_wall_shear": "NOT_COMPUTED_NO_SOURCE_PROVEN_RESTART_METHOD",
        }
        write_json(qc / "healthy_physiology_comparison.json", physiology)
        frozen_after = _file_manifest(frozen_paths)
        source_unchanged = frozen_before == frozen_after
        write_json(qc / "source_frozen_files_unchanged_qc.json", {
            "status": "PASS" if source_unchanged else "FAIL",
            "source_frozen_files_unchanged": source_unchanged,
            "before": frozen_before, "after": frozen_after,
        })
        status, next_step, first_failure = classify_healthy_calibration(
            official_steady=official, additional_iterations=additional,
            exact=exact, mass=mass, source_unchanged=source_unchanged,
        )
        pressure_drops = {
            label: {
                "pa": exact["inlet"]["gauge_pressure_pa"] - gauge,
                "mmHg": (exact["inlet"]["gauge_pressure_pa"] - gauge) / 133.322387415,
            }
            for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items()
        }
        summary.update({
            "status": status, "next": next_step, "first_failure": first_failure,
            "additional_iterations": additional, "final_absolute_iteration": final_iteration,
            "official_pu_steady": official, "steady_confirmation_iteration": steady.get("confirmation_iteration"),
            "runtime_controller": controller, "maximum_runtime_output_path_length": maximum_path_length,
            "final_restart_header": str(final_header), "final_restart_binary": str(final_binary),
            "final_restart_binary_sha256": sha256_file(final_binary),
            "mass_accumulation": mass, "exact_boundary_flux": exact,
            "pressure_drops": pressure_drops, "physiology": physiology,
            "source_frozen_files_unchanged": source_unchanged,
            "completed_at": datetime.now().isoformat(),
        })
        write_json(manifest_path, summary)
        return summary
    except Exception as error:
        status = error.status if isinstance(error, FlowError) else STATUS_FAILED
        next_step = NEXT_SAFETY if status == STATUS_SAFETY_FAILED else NEXT_REVIEW
        local_frozen_after = _file_manifest(frozen_paths[:-1])
        full_unchanged = False
        try:
            full_unchanged = frozen_before == _file_manifest(frozen_paths)
        except (OSError, ValueError):
            pass
        summary.update({
            "status": status, "next": next_step, "first_failure": str(error),
            "source_frozen_files_unchanged": full_unchanged or local_frozen_before == local_frozen_after,
            "completed_at": datetime.now().isoformat(),
        })
        if tracker is not None:
            summary["checkpoint_iterations"] = [int(item["iteration"]) for item in tracker.records]
            summary["latest_checkpoint_path"] = tracker.records[-1]["binary_path_wsl"] if tracker.records else None
        write_json(manifest_path, summary)
        return summary
