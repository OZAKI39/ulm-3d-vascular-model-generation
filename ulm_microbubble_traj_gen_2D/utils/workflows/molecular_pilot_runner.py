"""Connect saved no-bond trajectories to the dimensionless contact-pilot study."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..core.config import BindingScenarioSweepConfig, MolecularBindingConfig
from ..molecular.molecular_binding import reaction_disk_radius_um, surface_slip_velocity_um_s
from ..molecular.molecular_binding_scenarios import build_da_on_scenarios
from ..molecular.molecular_contact_pilot import (
    _analyze_prepared_molecular_contact_pilot,
    _prepare_molecular_contact_records,
    save_contact_pilot_yaml,
)
from ..molecular.molecular_target_field import MolecularTargetField
from ..particles.particle_trajectory_kinematics import realized_center_velocities_um_s
from ..core.types import ParticleTrajectories


_NUMERICAL_WALL_LOCK_SPEED_ATOL_UM_S = 1.0e-12


def _internal_substep_exposure_summary(
    trajectories: ParticleTrajectories,
    *,
    requested_capture_distance_um: float,
    transported_capture_distance_um: float,
) -> dict[str, object] | None:
    """Summarize authoritative internal-step exposure at its native distance.

    A configured scenario sweep may ask about other hypothetical capture
    distances.  Those continue to use the saved-frame pilot, but the ratio
    matching the transported molecular model must expose the internal-step
    registry result so the YAML cannot silently present the coarser estimator
    as authoritative.
    """

    requested = float(requested_capture_distance_um)
    transported = float(transported_capture_distance_um)
    scale = max(abs(requested), abs(transported), 1.0)
    if (
        transported <= 0.0
        or abs(requested - transported)
        > 256.0 * np.finfo(np.float64).eps * scale
    ):
        return None
    required = (
        trajectories.registry_target_exposure_time_s,
        trajectories.registry_target_exposure_event_count,
        trajectories.registry_target_reaction_area_time_um2_s,
        trajectories.registry_target_exposure_right_censored,
    )
    if any(value is None for value in required):
        return None
    count = int(trajectories.registry_bubble_id.size)
    exposure = np.asarray(required[0], dtype=np.float64).reshape(-1)
    events = np.asarray(required[1], dtype=np.int64).reshape(-1)
    area_time = np.asarray(required[2], dtype=np.float64).reshape(-1)
    censored = np.asarray(required[3], dtype=bool).reshape(-1)
    if not all(value.size == count for value in (exposure, events, area_time, censored)):
        raise ValueError(
            "Internal target-exposure registry arrays must match registry_bubble_id."
        )
    applicability = (
        None
        if trajectories.registry_target_exposure_quantitative_applicability_fraction
        is None
        else np.asarray(
            trajectories.registry_target_exposure_quantitative_applicability_fraction,
            dtype=np.float64,
        ).reshape(-1)
    )
    if applicability is not None and applicability.size != count:
        raise ValueError(
            "Internal target-exposure applicability must match registry_bubble_id."
        )
    total_exposure = float(np.sum(exposure, dtype=np.float64))
    weighted_applicability = (
        None
        if applicability is None or total_exposure <= 0.0
        else float(np.dot(exposure, applicability) / total_exposure)
    )
    return {
        "capture_distance_um": transported,
        "encountered_bubble_count_N_enc": int(np.count_nonzero(events > 0)),
        "independent_exposure_event_count": int(np.sum(events, dtype=np.int64)),
        "total_exposure_time_s": total_exposure,
        "total_reaction_area_time_E_T_um2_s": float(
            np.sum(area_time, dtype=np.float64)
        ),
        "right_censored_bubble_count": int(np.count_nonzero(censored)),
        "exposure_weighted_rbc_quantitative_applicability_fraction": (
            weighted_applicability
        ),
        "sampling": "accepted_internal_physical_substeps_and_zero_time_random_sweep",
    }


def _numerical_wall_lock_mask(
    trajectories: ParticleTrajectories,
    realized_velocities_um_s: np.ndarray,
) -> tuple[np.ndarray, str]:
    """Identify non-physical wall locks without depending on the retired slide solver.

    In a current trajectory, ``velocities_um_s`` stores the last accepted
    internal position-level translational velocity for each output record.  A
    contact-active record is therefore suspicious when this velocity is clearly
    non-zero but the same permanent bubble ID makes essentially no saved-frame
    displacement.  This saved-frame mask is only a secondary reporting aid; the
    authoritative v15 consistency check integrates position and velocity over
    accepted internal contact intervals.  The deliberately tiny absolute
    tolerance only treats numerical zero as zero; it is not a physical adhesion
    or low-flow threshold.

    Historical file compatibility belongs to the raw result-loading layer.  A
    current runtime trajectory intentionally contains no state from the retired
    wall-sliding algorithm.
    """

    count = int(trajectories.bubble_id.size)
    constraint_active = trajectories.contact_constraint_active
    if constraint_active is not None:
        active = np.asarray(constraint_active, dtype=bool)
        if active.shape != (count,):
            raise ValueError(
                "Contact-constraint activity must match the trajectory record count."
            )
        constrained_rhs_speed = np.linalg.norm(
            np.asarray(trajectories.velocities_um_s, dtype=np.float64)[:, (0, 2)],
            axis=1,
        )
        realized_speed = np.linalg.norm(
            np.asarray(realized_velocities_um_s, dtype=np.float64)[:, (0, 2)],
            axis=1,
        )
        finite = np.isfinite(constrained_rhs_speed) & np.isfinite(realized_speed)
        lock = (
            active
            & finite
            & (constrained_rhs_speed > _NUMERICAL_WALL_LOCK_SPEED_ATOL_UM_S)
            & (realized_speed <= _NUMERICAL_WALL_LOCK_SPEED_ATOL_UM_S)
        )
        return (
            lock,
            "v16_accepted_internal_velocity_vs_same_id_realized_velocity",
        )

    return np.zeros(count, dtype=bool), "unavailable"


def run_configured_contact_pilot(
    output_dir: Path,
    trajectories: ParticleTrajectories,
    target_field: MolecularTargetField,
    binding_cfg: MolecularBindingConfig,
    sweep_cfg: BindingScenarioSweepConfig,
) -> tuple[Path, ...]:
    """Write one no-bond contact report for each configured capture-distance ratio.

    Contact opportunity changes with capture distance, so every ratio receives
    its own exposure report.  The scenario table is instead built from the
    reference time fixed in configuration before transport; observed exposure
    never changes ``k_on``.  The routine never modifies trajectories or targets.
    """

    count = int(trajectories.bubble_id.size)
    if trajectories.positions_um.shape != (count, 3):
        raise ValueError("Trajectory positions must have shape (record_count, 3).")
    if trajectories.velocities_um_s.shape != (count, 3):
        raise ValueError("Trajectory velocities must have shape (record_count, 3).")
    if trajectories.wall_normal_xz.shape != (count, 2):
        raise ValueError("Trajectory wall normals must have shape (record_count, 2).")

    radii_um = 0.5 * np.asarray(trajectories.diameter_um, dtype=np.float64)
    gaps_um = np.asarray(trajectories.wall_gap_um, dtype=np.float64)
    normals = np.asarray(trajectories.wall_normal_xz, dtype=np.float64)
    normal_norm = np.linalg.norm(normals, axis=1)
    valid_normals = np.isfinite(normal_norm) & (normal_norm > 1.0e-12)
    safe_normals = np.zeros_like(normals)
    safe_normals[valid_normals] = normals[valid_normals] / normal_norm[valid_normals, None]
    # A centreline/medial-axis record may have no unique nearest-wall normal.
    # It is valid trajectory data but cannot define molecular wall contact, so
    # retain it in population statistics while assigning zero reaction area.
    safe_normals[~valid_normals, 1] = 1.0
    centres_xz_um = np.asarray(trajectories.positions_um[:, (0, 2)], dtype=np.float64)
    output_dt_s = float(trajectories.metadata.get("output_dt_s", 0.0))
    prepared_contact_records = _prepare_molecular_contact_records(
        trajectories.frame_offsets,
        trajectories.bubble_id,
        output_dt_s,
    )
    if trajectories.realized_velocities_um_s is None:
        realized_velocities = realized_center_velocities_um_s(
            trajectories.frame_offsets,
            trajectories.bubble_id,
            trajectories.positions_um,
            output_dt_s,
        )
    else:
        realized_velocities = np.asarray(
            trajectories.realized_velocities_um_s,
            dtype=np.float64,
        )
    if realized_velocities.shape != (count, 3):
        raise ValueError(
            "Realized trajectory velocities must have shape (record_count, 3)."
        )
    particle_velocity_xz_um_s = realized_velocities[:, (0, 2)]
    tangents = np.column_stack((-safe_normals[:, 1], safe_normals[:, 0]))
    tangential_velocity_um_s = np.sum(particle_velocity_xz_um_s * tangents, axis=1)
    slip_um_s = surface_slip_velocity_um_s(
        tangential_velocity_um_s,
        radii_um,
        np.asarray(trajectories.angular_velocity_rad_s, dtype=np.float64),
    )
    numerical_wall_lock, numerical_wall_lock_source = _numerical_wall_lock_mask(
        trajectories,
        realized_velocities,
    )

    output_paths: list[Path] = []
    for ratio_index, ratio in enumerate(sweep_cfg.capture_distance_to_rest_length_ratios):
        capture_distance_um = float(ratio) * float(binding_cfg.rest_length_um)
        internal_exposure = _internal_substep_exposure_summary(
            trajectories,
            requested_capture_distance_um=capture_distance_um,
            transported_capture_distance_um=float(binding_cfg.capture_distance_um),
        )
        reaction_radius_um = reaction_disk_radius_um(
            radii_um,
            gaps_um,
            capture_distance_um,
        )
        reaction_radius_um[~valid_normals] = 0.0
        # Molecular geometry is needed only inside the capture layer.  In a
        # long no-bond trajectory this is normally a tiny fraction of all saved
        # records; sending the dense population to the exact nearest-face query
        # used to allocate one tied-face object per record and dominate the
        # complete pilot runtime.
        reaction_area_um2 = np.zeros(count, dtype=np.float64)
        reaction_lanes = np.flatnonzero(reaction_radius_um > 0.0)
        if reaction_lanes.size:
            active_area = np.asarray(
                target_field.reaction_area_um2(
                    centres_xz_um[reaction_lanes],
                    tangents[reaction_lanes],
                    reaction_radius_um[reaction_lanes],
                ),
                dtype=np.float64,
            )
            if active_area.shape != (reaction_lanes.size,):
                raise ValueError(
                    "Molecular target reaction areas must match the capture-eligible "
                    "trajectory record count."
                )
            reaction_area_um2[reaction_lanes] = active_area
        pilot = _analyze_prepared_molecular_contact_pilot(
            prepared_contact_records,
            reaction_area_um2,
            slip_um_s,
            numerical_wall_lock=numerical_wall_lock,
        )
        scenarios = build_da_on_scenarios(
            da_on_reference_time_s=sweep_cfg.da_on_reference_time_s,
            da_on_levels=sweep_cfg.da_on_levels,
            target_density_molecules_per_um2_levels=(
                sweep_cfg.target_density_molecules_per_um2_levels
            ),
            ligand_density_molecules_per_um2_levels=(
                sweep_cfg.ligand_density_molecules_per_um2_levels
            ),
            capture_distance_to_rest_length_ratios=(float(ratio),),
            rest_length_um=float(binding_cfg.rest_length_um),
        )
        ratio_label = f"{float(ratio):.6g}".replace("-", "m").replace(".", "p")
        path = output_dir / (
            f"molecular_contact_pilot_{ratio_index:02d}_capture_ratio_{ratio_label}.yaml"
        )
        output_paths.append(
            save_contact_pilot_yaml(
                path,
                pilot,
                scenarios,
                study_context={
                    "capture_distance_to_rest_length_ratio": float(ratio),
                    "rest_length_um": float(binding_cfg.rest_length_um),
                    "capture_distance_um": capture_distance_um,
                    "exposure_authority": (
                        "internal_physical_substep_registry"
                        if internal_exposure is not None
                        else "saved_frame_endpoint_pilot_for_hypothetical_capture_distance"
                    ),
                    "internal_physical_substep_exposure": internal_exposure,
                    "saved_frame_pilot_scope": (
                        "secondary_compatibility_estimator; it may miss contacts "
                        "between saved frame endpoints"
                    ),
                    "reaction_area_integration": (
                        "analytic_target_roi_intersected_with_eligible_closed_wall_intervals"
                    ),
                    "target_region_mode": str(
                        getattr(target_field, "region_mode", "test_or_external_field")
                    ),
                    "target_density_molecules_per_m2_in_pilot": float(
                        getattr(target_field, "target_density_molecules_per_m2", 0.0)
                    ),
                    "eligible_wall_site_count": int(
                        np.count_nonzero(
                            getattr(target_field, "solid_wall_mask", np.empty(0, dtype=bool))
                        )
                    ),
                    "target_wall_site_count": int(
                        np.count_nonzero(
                            getattr(target_field, "target_wall_mask", np.empty(0, dtype=bool))
                        )
                    ),
                    "source_mask_npz_path": (
                        str(getattr(target_field, "source_mask_npz_path"))
                        if getattr(target_field, "source_mask_npz_path", None) is not None
                        else None
                    ),
                    "source_trajectory_binding_enabled": False,
                    "da_on_reference_time_s": float(
                        sweep_cfg.da_on_reference_time_s
                    ),
                    "association_rate_parameter_source": (
                        "predeclared_dimensionless_scenario_independent_of_observed_exposure"
                    ),
                    "source_trajectory_schema": str(
                        trajectories.metadata.get("trajectory_schema", "unknown")
                    ),
                    "source_trajectory_record_count": int(trajectories.bubble_id.size),
                    "centre_velocity_source": (
                        "same_id_geometry_accepted_saved_frame_displacement"
                    ),
                    "numerical_wall_lock_record_count": int(
                        np.count_nonzero(numerical_wall_lock)
                    ),
                    "numerical_wall_lock_detection": numerical_wall_lock_source,
                    "numerical_wall_lock_speed_atol_um_s": float(
                        _NUMERICAL_WALL_LOCK_SPEED_ATOL_UM_S
                    ),
                    "source_trajectory_registry_count": int(
                        trajectories.registry_bubble_id.size
                    ),
                    "records_without_unique_wall_normal": int(
                        np.count_nonzero(~valid_normals)
                    ),
                    "scenario_execution": (
                        "definitions_only; derived cases are not automatically rerun"
                    ),
                },
            )
        )
    return tuple(output_paths)
