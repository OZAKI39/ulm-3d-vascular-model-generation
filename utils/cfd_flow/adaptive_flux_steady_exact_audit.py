"""Independent exact PDF-link audit of an adaptive-flux steady restart."""

from __future__ import annotations

import math
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .adaptive_flux_pressure_audit import (
    physical_mass_factor,
    musubi_pressure_flux_affine_coefficients,
    reconstruct_musubi_boundary_state,
    solve_boundary_density,
)
from .adaptive_flux_steady import (
    AXIS_MESH_RUN,
    BULK_NU_M2_S,
    EXPECTED_FLUID_CELLS,
    EXPECTED_INLET_GLOBBC,
    NU_M2_S,
    PRESSURE_REFERENCE_PA,
    STEADY_PENDING_AUDIT,
    STEADY_PREFIX,
)
from .adaptive_flux_validation import EXPECTED_DT_S
from .exact_link_flux import (
    CS2,
    EXPECTED_DX_M,
    GLOBAL_NORMALS,
    REFERENCE_DENSITY_KG_M3,
    TARGET_MASS_FLOW_KG_S,
    TARGET_Q_M3_S,
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
from .io import FlowError, read_json, write_json
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
    tree_ids_to_ijk,
    tree_levels,
)


STEADY_BASELINE_PASS = "CFD_FLOW_MCCLURE_ADAPTIVE_FLUX_STEADY_BASELINE_PASS"
RUNTIME_TARGET_MISMATCH = "CFD_FLOW_ADAPTIVE_FLUX_RUNTIME_TARGET_MISMATCH"
MASS_BALANCE_FAILED = "CFD_FLOW_ADAPTIVE_FLUX_MASS_BALANCE_FAILED"
EXACT_AUDIT_FAILED = "CFD_FLOW_MCCLURE_ADAPTIVE_FLUX_STEADY_EXACT_AUDIT_FAILED"
NEXT_GRID = "RUN ADAPTIVE-FLUX GRID CONVERGENCE"
NEXT_TARGET_REVIEW = "REVIEW DIFFERENCE BETWEEN CONTROLLER RESIDUAL AND EXACT PDF-LINK FLUX"
NEXT_OUTLET_REVIEW = "REVIEW PRESSURE OUTLET EXACT FLUX AT STEADY STATE"
NEXT_AUDIT_REVIEW = "REVIEW INDEPENDENT FINAL-RESTART PDF-LINK AUDIT"

MAXIMUM_RELATIVE_ERROR = 0.01
PREFERRED_RELATIVE_ERROR = 0.001
OLD_MFR_EQ_AREA_RATIO = 1.4680772275174907
OLD_MFR_EQ_VELOCITY_RATIO = 0.6811630759309535
INLET_PRESSURE_NEIGHBOR_DIRECTION = np.asarray((0, 0, -1), dtype=np.int8)


def _git(root: Path, *arguments: str) -> str:
    return subprocess.run(
        ["git", "-C", str(root), *arguments],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    ).stdout.strip()


def _pressure_neighbor_indices(
    cell_indices: np.ndarray,
    cell_ijk: np.ndarray,
    lookup: dict[tuple[int, int, int], int],
    direction: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply the existing pressure_eq two-neighbor rule for one known direction."""

    normal = np.asarray(direction, dtype=np.int64).reshape(3)
    neighbors1 = np.full(len(cell_indices), -1, dtype=np.int64)
    neighbors2 = np.full(len(cell_indices), -1, dtype=np.int64)
    for row, cell_index in enumerate(np.asarray(cell_indices, dtype=np.int64)):
        coordinate = cell_ijk[cell_index]
        neighbors1[row] = lookup.get(
            tuple(int(value) for value in coordinate + normal), -1
        )
        neighbors2[row] = lookup.get(
            tuple(int(value) for value in coordinate + 2 * normal), -1
        )
    valid = (neighbors1 >= 0) & (neighbors2 >= 0)
    return valid, neighbors1, neighbors2


def classify_steady_exact_audit(
    *,
    inlet_relative_error: float,
    mass_balance_relative_error: float,
    all_pdfs_finite: bool,
    maximum_lattice_speed: float,
    minimum_pdf: float,
) -> tuple[str, str]:
    safety_pass = (
        all_pdfs_finite
        and math.isfinite(maximum_lattice_speed)
        and math.isfinite(minimum_pdf)
        and maximum_lattice_speed < 0.05
        and minimum_pdf > 0.0
    )
    if not safety_pass:
        return EXACT_AUDIT_FAILED, NEXT_AUDIT_REVIEW
    if inlet_relative_error > MAXIMUM_RELATIVE_ERROR:
        return RUNTIME_TARGET_MISMATCH, NEXT_TARGET_REVIEW
    if mass_balance_relative_error > MAXIMUM_RELATIVE_ERROR:
        return MASS_BALANCE_FAILED, NEXT_OUTLET_REVIEW
    return STEADY_BASELINE_PASS, NEXT_GRID


def _latest_pending_steady(output_root: Path) -> Path:
    candidates = sorted(
        path for path in output_root.glob(f"{STEADY_PREFIX}_*") if path.is_dir()
    )
    for run_root in reversed(candidates):
        manifest = run_root / "qc" / "adaptive_flux_steady_manifest.json"
        if manifest.is_file() and read_json(manifest).get("status") == STEADY_PENDING_AUDIT:
            return run_root
    raise FlowError(EXACT_AUDIT_FAILED, "No steady run is pending exact audit")


def _outlet_pressures() -> dict[str, float]:
    return {
        "outlet_01": PRESSURE_REFERENCE_PA + 14.544978101274268,
        "outlet_02": PRESSURE_REFERENCE_PA + 132.20454922317552,
        "outlet_03": PRESSURE_REFERENCE_PA - 13.700626673311461,
    }


def run_adaptive_flux_steady_exact_audit(
    project_root: Path, steady_run_root: Path | None = None
) -> dict[str, Any]:
    """Read final PDFs directly; controller logs are never a flux data source."""

    root = Path(project_root).resolve()
    head = _git(root, "rev-parse", "HEAD")
    output_root = root / "outputs" / "cfd_flow"
    steady_root = (
        Path(steady_run_root).resolve()
        if steady_run_root is not None
        else _latest_pending_steady(output_root)
    )
    steady_manifest_path = steady_root / "qc" / "adaptive_flux_steady_manifest.json"
    steady = read_json(steady_manifest_path)
    if steady.get("status") != STEADY_PENDING_AUDIT:
        raise FlowError(EXACT_AUDIT_FAILED, "Steady criterion did not pass before exact audit")
    if steady.get("actual_head") != head:
        raise FlowError(EXACT_AUDIT_FAILED, "Project HEAD changed after the steady run")
    mesh = output_root / AXIS_MESH_RUN / "seeder" / "mesh"
    restart_header_path = Path(steady["final_restart_header"])
    restart_header = parse_restart_header(restart_header_path)
    restart_binary = Path(steady["final_restart_binary"])

    frozen_paths = tuple(Path(path) for path in steady["frozen_files_before"])
    critical_paths = (
        *frozen_paths,
        steady_root / "musubi.lua",
        steady_manifest_path,
        restart_header_path,
        restart_binary,
    )
    before = _file_manifest(critical_paths)
    original_frozen_unchanged = steady["frozen_files_before"] == _file_manifest(
        frozen_paths
    )

    audit_root = steady_root / "exact_audit"
    audit_root.mkdir(parents=True, exist_ok=False)
    manifest_path = audit_root / "adaptive_flux_steady_exact_audit_manifest.json"
    summary: dict[str, Any] = {
        "status": EXACT_AUDIT_FAILED,
        "next": NEXT_AUDIT_REVIEW,
        "run_root": str(steady_root),
        "audit_root": str(audit_root),
        "actual_head": head,
        "production_pipeline_modified": False,
        "mesh_path": str(mesh),
        "restart_header": str(restart_header_path),
        "restart_binary": str(restart_binary),
        "seeder_calls": 0,
        "musubi_calls": 0,
        "harvester_calls": 0,
        "grid_convergence": "NOT_RUN",
        "controller_output_used_as_flux_source": False,
        "third_treelm_parser_written": False,
        "started_at": datetime.now().isoformat(),
    }
    write_json(manifest_path, summary)

    try:
        if not original_frozen_unchanged:
            raise FlowError(EXACT_AUDIT_FAILED, "Frozen source changed before exact audit")
        if (
            restart_header.n_elems != EXPECTED_FLUID_CELLS
            or restart_header.n_components != 19
            or restart_header.n_dofs != 1
        ):
            raise FlowError(EXACT_AUDIT_FAILED, "Final restart dimensions changed")

        property_header = parse_boundary_property_header(
            (mesh / "header.lua").read_text(encoding="utf-8")
        )
        boundary_header = parse_bnd_header(
            (mesh / "bnd.lua").read_text(encoding="utf-8")
        )
        tree_ids, property_bits, elemlist_contract = read_treelm_elemlist(
            mesh / "elemlist.lsb", n_elems=EXPECTED_FLUID_CELLS
        )
        levels = tree_levels(tree_ids)
        if not np.all(levels == 9):
            raise FlowError(EXACT_AUDIT_FAILED, "Adaptive mesh is not uniform level 9")
        cell_ijk = tree_ids_to_ijk(tree_ids, levels)
        property_indices = extract_boundary_property_indices(
            property_bits, property_header.bit_position
        )
        boundary_ids = read_boundary_ids(
            mesh / "bnd.lsb",
            element_count=property_header.element_count,
            side_count=boundary_header.side_count,
        )
        label_to_id = {
            label: index
            for index, label in enumerate(boundary_header.labels, start=1)
        }
        boundaries = {
            label: reconstruct_boundary(
                boundary_ids,
                property_indices,
                label=label,
                boundary_id=label_to_id[label],
            )
            for label in PORT_LABELS
        }
        boundary_counts = {
            label: len(boundary.cell_indices)
            for label, boundary in boundaries.items()
        }
        if boundary_counts["inlet"] != EXPECTED_INLET_GLOBBC:
            raise FlowError(
                EXACT_AUDIT_FAILED,
                f"Adaptive inlet globBC changed: {boundary_counts['inlet']}",
            )

        pdf = read_restart_pdf(
            restart_binary,
            n_elems=EXPECTED_FLUID_CELLS,
            n_components=19,
        )
        values = np.asarray(pdf)
        all_finite = bool(np.all(np.isfinite(values)))
        minimum_pdf = float(np.min(values))
        density = np.sum(values, axis=1, dtype=np.float64)
        if not all_finite or np.any(density <= 0.0):
            raise FlowError(EXACT_AUDIT_FAILED, "Final restart PDFs are non-finite or non-positive density")
        velocity = (values @ D3Q19_DIRECTIONS.astype(np.float64)) / density[:, None]
        maximum_lattice_speed = float(np.max(np.linalg.norm(velocity, axis=1)))
        lookup = build_coordinate_lookup(cell_ijk)

        inlet = boundaries["inlet"]
        inlet_valid, inlet_neighbor1, inlet_neighbor2 = _pressure_neighbor_indices(
            inlet.cell_indices,
            cell_ijk,
            lookup,
            INLET_PRESSURE_NEIGHBOR_DIRECTION,
        )
        if not np.all(inlet_valid):
            raise FlowError(
                EXACT_AUDIT_FAILED,
                f"Adaptive inlet lacks two pressure neighbors for {np.count_nonzero(~inlet_valid)} cells",
            )
        fetched1 = pull_fetch_pdfs(
            pdf, cell_ijk, inlet_neighbor1, coordinate_lookup=lookup
        )
        fetched2 = pull_fetch_pdfs(
            pdf, cell_ijk, inlet_neighbor2, coordinate_lookup=lookup
        )
        extrapolated_velocity = 1.5 * velocity_from_pdf(fetched1) - 0.5 * velocity_from_pdf(
            fetched2
        )
        stored = np.asarray(pdf[inlet.cell_indices], dtype=np.float64)
        alpha, beta = musubi_pressure_flux_affine_coefficients(
            stored_boundary_pdfs=stored,
            incoming_masks=inlet.incoming_masks,
            extrapolated_velocity=extrapolated_velocity,
        )
        target_flux_lattice = physical_volume_flux_to_lattice(
            TARGET_Q_M3_S, dx_m=EXPECTED_DX_M, dt_s=EXPECTED_DT_S
        )
        rho_boundary = solve_boundary_density(target_flux_lattice, alpha, beta)
        fetched_boundary = pull_fetch_pdfs(
            pdf, cell_ijk, inlet.cell_indices, coordinate_lookup=lookup
        )
        reconstructed_inlet = reconstruct_musubi_boundary_state(
            rho_boundary,
            fetched_boundary_pdfs=fetched_boundary,
            incoming_masks=inlet.incoming_masks,
            extrapolated_velocity=extrapolated_velocity,
        )
        mass_factor = physical_mass_factor(
            density_kg_m3=REFERENCE_DENSITY_KG_M3,
            dx_m=EXPECTED_DX_M,
            dt_s=EXPECTED_DT_S,
        )
        inlet_rows = np.arange(len(inlet.cell_indices), dtype=np.int64)
        link_records, exact_inlet_mass = _link_rows(
            boundary=inlet,
            selected_rows=inlet_rows,
            tree_ids=tree_ids,
            old_pdf=pdf,
            new_pdf=reconstructed_inlet,
            mass_factor=mass_factor,
            outward_sign=False,
        )
        inlet_relative_error = (
            abs(exact_inlet_mass - TARGET_MASS_FLOW_KG_S)
            / TARGET_MASS_FLOW_KG_S
        )

        pressure_factor = (
            REFERENCE_DENSITY_KG_M3 * EXPECTED_DX_M**2 / EXPECTED_DT_S**2
        )
        pressure_reference = pressure_factor * CS2
        final_absolute_pressure = rho_boundary * pressure_reference
        final_gauge_pressure = final_absolute_pressure - pressure_reference
        outlet_mass: dict[str, float] = {}
        outlet_qc: dict[str, Any] = {}
        for label, absolute_pressure in _outlet_pressures().items():
            boundary = boundaries[label]
            valid, neighbor1, neighbor2 = _pressure_neighbor_indices(
                boundary.cell_indices,
                cell_ijk,
                lookup,
                GLOBAL_NORMALS[label],
            )
            selected = np.flatnonzero(valid).astype(np.int64)
            if len(selected) == 0:
                raise FlowError(EXACT_AUDIT_FAILED, f"No valid pressure cells for {label}")
            fetched1 = pull_fetch_pdfs(
                pdf, cell_ijk, neighbor1[selected], coordinate_lookup=lookup
            )
            fetched2 = pull_fetch_pdfs(
                pdf, cell_ijk, neighbor2[selected], coordinate_lookup=lookup
            )
            outlet_velocity = 1.5 * velocity_from_pdf(fetched1) - 0.5 * velocity_from_pdf(
                fetched2
            )
            outlet_density = absolute_pressure / pressure_reference
            replacement = equilibrium_pdf(outlet_density, outlet_velocity)
            records, mass = _link_rows(
                boundary=boundary,
                selected_rows=selected,
                tree_ids=tree_ids,
                old_pdf=pdf,
                new_pdf=replacement,
                mass_factor=mass_factor,
                outward_sign=True,
            )
            link_records.extend(records)
            outlet_mass[label] = mass
            outlet_qc[label] = {
                "globbc_count": len(boundary.cell_indices),
                "valid_pressure_count": len(selected),
                "removed_solid_neighbor_count": int(np.count_nonzero(~valid)),
                "pressure_neighbor_direction": GLOBAL_NORMALS[label].astype(int).tolist(),
                "absolute_pressure_pa": absolute_pressure,
                "exact_signed_mass_flow_kg_s": mass,
            }

        link_csv = audit_root / "exact_boundary_link_replacements.csv"
        _write_link_csv(link_csv, link_records)
        balance = signed_mass_balance(exact_inlet_mass, outlet_mass.values())
        final_status, next_step = classify_steady_exact_audit(
            inlet_relative_error=inlet_relative_error,
            mass_balance_relative_error=balance["relative_error"],
            all_pdfs_finite=all_finite,
            maximum_lattice_speed=maximum_lattice_speed,
            minimum_pdf=minimum_pdf,
        )

        after = _file_manifest(critical_paths)
        unchanged = before == after and original_frozen_unchanged
        if not unchanged:
            raise FlowError(EXACT_AUDIT_FAILED, "Source or final restart changed during exact audit")
        exact_qc = {
            "status": "PASS" if final_status == STEADY_BASELINE_PASS else "FAIL",
            "method": "independent final-restart PULL/FETCH PDF replacement link deltas",
            "runtime_controller_output_used": False,
            "inlet": {
                "pressure_neighbor_direction": INLET_PRESSURE_NEIGHBOR_DIRECTION.astype(int).tolist(),
                "globbc_count": len(inlet.cell_indices),
                "alpha_recomputed_from_final_pdfs": alpha,
                "beta_recomputed_from_final_pdfs": beta,
                "rho_boundary_recomputed_from_final_pdfs": rho_boundary,
                "replacement_pdfs_constructed": True,
                "flux_sum_expression": "sum(new incoming PDF - old outgoing PDF) over boundary links",
                "target_mass_flow_kg_s": TARGET_MASS_FLOW_KG_S,
                "exact_mass_flow_kg_s": exact_inlet_mass,
                "target_relative_error": inlet_relative_error,
                "preferred_gate": PREFERRED_RELATIVE_ERROR,
                "hard_gate": MAXIMUM_RELATIVE_ERROR,
            },
            "outlets": outlet_qc,
            "outlet_signed_sum_kg_s": balance["outlet_signed_sum_kg_s"],
            "mass_balance_relative_error": balance["relative_error"],
            "mass_balance_preferred_gate": PREFERRED_RELATIVE_ERROR,
            "mass_balance_hard_gate": MAXIMUM_RELATIVE_ERROR,
            "all_pdfs_finite": all_finite,
            "maximum_lattice_speed": maximum_lattice_speed,
            "minimum_pdf": minimum_pdf,
            "final_rho_boundary": rho_boundary,
            "final_inlet_absolute_pressure_pa": final_absolute_pressure,
            "final_inlet_gauge_pressure_pa": final_gauge_pressure,
            "per_link_csv": str(link_csv),
        }
        write_json(audit_root / "adaptive_flux_steady_exact_flux_qc.json", exact_qc)
        write_json(
            audit_root / "source_frozen_files_unchanged_qc.json",
            {
                "status": "PASS",
                "source_frozen_files_unchanged": True,
                "before": before,
                "after": after,
            },
        )
        summary.update(
            {
                "status": final_status,
                "next": next_step,
                "steady_status": steady["status"],
                "steady_iterations": steady["total_steady_iterations"],
                "physics": {
                    "dx_m": EXPECTED_DX_M,
                    "dt_s": EXPECTED_DT_S,
                    "density_kg_m3": REFERENCE_DENSITY_KG_M3,
                    "kinematic_viscosity_m2_s": NU_M2_S,
                    "bulk_viscosity_m2_s": BULK_NU_M2_S,
                },
                "boundary_counts": boundary_counts,
                "target_q_m3_s": TARGET_Q_M3_S,
                "target_mass_flow_kg_s": TARGET_MASS_FLOW_KG_S,
                "final_rho_boundary": rho_boundary,
                "final_inlet_absolute_pressure_pa": final_absolute_pressure,
                "final_inlet_gauge_pressure_pa": final_gauge_pressure,
                "independent_exact_inlet_mass_flow_kg_s": exact_inlet_mass,
                "independent_inlet_target_relative_error": inlet_relative_error,
                "exact_outlet_mass_flow_kg_s": outlet_mass,
                "outlet_signed_sum_kg_s": balance["outlet_signed_sum_kg_s"],
                "exact_mass_balance_relative_error": balance["relative_error"],
                "all_pdfs_finite": all_finite,
                "maximum_lattice_speed": maximum_lattice_speed,
                "minimum_pdf": minimum_pdf,
                "old_mfr_eq_area_ratio": OLD_MFR_EQ_AREA_RATIO,
                "old_mfr_eq_velocity_ratio": OLD_MFR_EQ_VELOCITY_RATIO,
                "source_frozen_files_unchanged": True,
                "elemlist_contract": elemlist_contract,
                "link_csv": str(link_csv),
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(manifest_path, summary)
        return summary
    except Exception as error:
        after = _file_manifest(critical_paths)
        summary.update(
            {
                "status": EXACT_AUDIT_FAILED,
                "next": NEXT_AUDIT_REVIEW,
                "first_failure": str(error),
                "source_frozen_files_unchanged": before == after
                and original_frozen_unchanged,
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(manifest_path, summary)
        return summary
