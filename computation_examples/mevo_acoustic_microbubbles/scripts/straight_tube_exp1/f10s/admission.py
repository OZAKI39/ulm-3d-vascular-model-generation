"""F10P gates with explicit requested-window coverage; thresholds unchanged."""
import numpy as np

def layered(record, dt, cache, topologies, policy, duration_s):
    complete = record["status"] == "PASS" and abs(record["completed_s"] - duration_s) < 1e-13
    available = dt is not None and cache is not None
    comparisons = [dt, cache] if available else []
    full = available and all(c["complete_window"] for c in comparisons)
    position = complete and full and all(c["COM_m"] <= policy["group_COM_tolerance_m"] for c in comparisons)
    signal = float(np.linalg.norm(record["displacement_m"]))
    noise = max([c["COM_m"] for c in comparisons], default=float("inf"))
    signal_verified = signal > noise
    # OFF is a numerically admitted null reference, not demonstrated actuation.
    position = position and (signal_verified or record["source"] == "OFF")
    shape = (
        complete
        and full
        and all(
            c["Rg_m"] <= policy["group_Rg_tolerance_m"]
            and c["directional_xy"] <= policy["group_xy_tolerance"]
            and c["Q_relative"] <= policy["group_Q_relative"]
            and c["aspect_relative"] <= policy["group_aspect_relative"]
            and c["matched_RMS_over_R"] <= policy["distribution_RMS_over_R"]
            and (c["axis_angle_deg"] is None or c["axis_angle_deg"] <= policy["principal_axis_angle_deg"])
            for c in comparisons
        )
    )
    particle = (
        complete and full and all(c["fixed_id_m"] <= policy["fixed_id_tolerance_m"] for c in comparisons)
    )
    # Conservative gate: this initial implementation admits survival of the initial
    # single network only. Split/merge cases retain ledgers but require explicit review.
    topo = (
        complete
        and full
        and len(topologies) == 3
        and all(t["initial_network_persists_entire_window"] for t in topologies)
    )
    return dict(
        BASIC_NUMERICS_PASS="PASS" if complete else "FAIL",
        GROUP_POSITION_NUMERICALLY_ADMITTED="PASS" if position else ("FAIL" if available else "NOT_RUN"),
        GROUP_SHAPE_NUMERICALLY_ADMITTED="PASS" if shape else ("FAIL" if available else "NOT_RUN"),
        GROUP_TOPOLOGY_NUMERICALLY_ADMITTED="PASS" if topo else ("FAIL" if available else "NOT_RUN"),
        PARTICLE_ID_NUMERICALLY_ADMITTED="PASS" if particle else ("FAIL" if available else "NOT_RUN"),
        particle_status="PARTICLE_ID_NUMERICALLY_ADMITTED" if particle else "PARTICLE_ID_NOT_ADMITTED",
        simulation_group_position_eligible=bool(position),
        simulation_group_shape_eligible=bool(shape),
        simulation_topology_eligible=bool(topo),
        simulation_particle_id_eligible=bool(particle),
        control_signal_m=signal,
        COM_noise_m=noise if np.isfinite(noise) else None,
        control_signal_exceeds_numerical_difference=bool(signal_verified),
        topology_scope="Initial 32-member network survival for complete requested window; later events may be right-censored. No physical fusion.",
        EXPERIMENTALLY_CALIBRATED="NOT_VERIFIED",
        TRAINING_READY=False,
        training_ready=False,
    )

