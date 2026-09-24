"""Permanent Particle-1 development cases. Every dt here is VALIDATION_ONLY."""
from dataclasses import dataclass
from itertools import combinations
import numpy as np
from .field import FlowSample
from .microbubble import MicrobubbleState, equilibrium_state
from .integrator import advance_single_microbubble, euler_position
from .validation_boundary import ValidationBoundaryClassifier

EPS = np.finfo(np.float64).eps
SYNTHETIC_ATOL = 256 * EPS  # fixed O(1) angular rate float64 budget
SYNTHETIC_RTOL = 256 * EPS
VALIDATION_ROLE = "VALIDATION_ONLY"


@dataclass(frozen=True)
class AffineValidationField:
    offset_m_s: np.ndarray
    gradient_s_inv: np.ndarray

    def sample(self, position_m):
        g = self.gradient_s_inv
        curl = np.array([g[2, 1] - g[1, 2], g[0, 2] - g[2, 0], g[1, 0] - g[0, 1]])
        return FlowSample(self.offset_m_s + g @ np.asarray(position_m), 10., g.copy(), curl,
                          .5 * (g + g.T), True, 0)


def uniform_case():
    velocity = np.array([2e-4, -1e-4, .5e-4])
    initial = np.array([1e-5, 2e-5, -1e-5])
    field = AffineValidationField(velocity, np.zeros((3, 3)))
    duration = .05  # VALIDATION_ONLY known straight-line case
    dts = np.array([.0025, .00125, .000625])
    rows = []
    for dt in dts:
        state = MicrobubbleState(0, initial, 1e-6, np.zeros(3), np.zeros(3))
        state = equilibrium_state(state, field.sample(state.position_m))
        count = int(round(duration / dt))
        # Summation roundoff grows with step count and coordinate magnitude.
        bound = 16 * (count + 1) * EPS * np.max(np.abs(initial) + np.abs(velocity) * duration)
        for step in range(count + 1):
            t = step * dt
            exact = initial + velocity * t
            row = dict(dt_s=float(dt), timestep_role=VALIDATION_ROLE, step=step, time_s=t,
                       position_error_m=float(np.max(np.abs(state.position_m - exact))), position_error_bound_m=bound,
                       velocity_error_m_s=float(np.max(np.abs(state.velocity_m_s - velocity))),
                       angular_velocity_norm_s_inv=float(np.linalg.norm(state.angular_velocity_s_inv)))
            row.update({f"{axis}_m": state.position_m[i] for i, axis in enumerate("xyz")})
            row.update({f"exact_{axis}_m": exact[i] for i, axis in enumerate("xyz")})
            row.update({f"velocity_{axis}_m_s": state.velocity_m_s[i] for i, axis in enumerate("xyz")})
            rows.append(row)
            if step < count: state = advance_single_microbubble(state, field, dt)
    return rows


def rotation_case():
    omega = np.array([1.25, -2., 3.])
    x, y, z = omega
    g = np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])
    field = AffineValidationField(np.zeros(3), g)
    positions = np.array([[a, b, 0.] for a in np.linspace(-1e-5, 1e-5, 7) for b in np.linspace(-1e-5, 1e-5, 7)])
    rows = []
    for idx, position in enumerate(positions):
        sample = field.sample(position)
        state = equilibrium_state(MicrobubbleState(0, position, 1e-6, np.zeros(3), np.zeros(3)), sample)
        row = dict(query_id=idx, x_m=position[0], y_m=position[1], z_m=position[2],
                   angular_velocity_error_s_inv=float(np.max(np.abs(state.angular_velocity_s_inv - omega))),
                   strain_norm_s_inv=float(np.linalg.norm(sample.strain_rate_s_inv)))
        for i, axis in enumerate("xyz"):
            row.update({f"analytic_omega_{axis}_s_inv": omega[i], f"sampled_omega_{axis}_s_inv": state.angular_velocity_s_inv[i],
                        f"vorticity_{axis}_s_inv": sample.vorticity_s_inv[i], f"velocity_{axis}_m_s": sample.velocity_m_s[i]})
        rows.append(row)
    return rows


def single_step_case():
    g = np.array([[2., -3., .5], [4., .25, -2.], [-1.5, 3., 1.]])
    field = AffineValidationField(np.array([2e-4, -1e-4, 1.5e-4]), g)
    initial = MicrobubbleState(0, [1e-5, 2e-5, 3e-5], 1e-6, np.zeros(3), np.zeros(3))
    dt = .001  # VALIDATION_ONLY, not a default argument
    sample = field.sample(initial.position_m)
    state = advance_single_microbubble(initial, field, dt)
    return dict(timestep_role=VALIDATION_ROLE, dt_s=dt, old_position_m=initial.position_m.tolist(),
                old_sampled_velocity_m_s=sample.velocity_m_s.tolist(), old_vorticity_s_inv=sample.vorticity_s_inv.tolist(),
                old_omega_s_inv=(.5 * sample.vorticity_s_inv).tolist(), new_position_m=state.position_m.tolist(),
                exact_euler_position_m=(initial.position_m + dt * sample.velocity_m_s).tolist(),
                endpoint_velocity_m_s=state.velocity_m_s.tolist(), endpoint_omega_s_inv=state.angular_velocity_s_inv.tolist(),
                input_position_after_call_m=initial.position_m.tolist(), radius_m=float(state.radius_m),
                production_particle_timestep_frozen=False)


def select_inlet_centroid(field, inlet):
    """Max inward centroid velocity among ALL inlet-adjacent canonical tetra;
    exact score ties resolve to the smallest canonical ID. No trajectory trial.
    """
    faces = inlet.faces.reshape(-1, 4)[:, 1:]
    xyz = np.asarray(inlet.points)[faces]
    area_vectors = .5 * np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
    inward = -area_vectors.sum(axis=0); inward /= np.linalg.norm(inward)
    ids = np.unique(np.asarray(inlet.cell_data["GlobalElementID"], dtype=np.int64) - 1)
    centroids = field.points_m[field.tetra[ids]].mean(axis=1)
    # A tetra centroid has four equal weights; evaluate candidate scores directly
    # from the already-validated P0 nodal representation, then confirm selection.
    velocities = field.velocity_nodes_m_s[field.tetra[ids]].mean(axis=1)
    scores = velocities @ inward
    index = np.lexsort((ids, -scores))[0]
    chosen = int(ids[index]); position = centroids[index]
    sample = field.sample(position)
    if not sample.inside_lumen or sample.tetra_id != chosen or scores[index] <= 0:
        raise ValueError("selected inlet-adjacent centroid is not an inward-flowing interior point")
    vertices = field.points_m[field.tetra[chosen]]
    shortest_edge = min(np.linalg.norm(vertices[a] - vertices[b]) for a, b in combinations(range(4), 2))
    speed = float(np.linalg.norm(sample.velocity_m_s))
    dt0 = float(.25 * shortest_edge / speed)
    horizon = float(4 * np.linalg.norm(np.ptp(field.points_m, axis=0)) / speed)
    record = dict(selection_rule="maximum dot(centroid nodal-mean velocity, area-weighted INLET inward normal); exact tie -> minimum canonical tetra_id",
                  candidate_count=len(ids), initial_tetra_id=chosen, initial_position_m=position.tolist(),
                  initial_velocity_m_s=sample.velocity_m_s.tolist(), inward_normal=inward.tolist(), inward_speed_m_s=float(scores[index]),
                  local_shortest_edge_m=float(shortest_edge), initial_speed_m_s=speed,
                  validation_dt_rule="dt0=0.25*initial_tetra_shortest_edge/initial_speed; dt0/2; dt0/4; fixed before trajectories",
                  validation_timesteps_s=[dt0, dt0 / 2, dt0 / 4], timestep_role=VALIDATION_ROLE,
                  horizon_s=horizon, horizon_rule="4*domain_bounding_box_diagonal/initial_speed; validation stop budget, not production",
                  production_particle_timestep_frozen=False)
    candidates = [dict(tetra_id=int(cell), inward_speed_m_s=float(score), **{f"{a}_m": float(p[i]) for i, a in enumerate("xyz")})
                  for cell, score, p in zip(ids, scores, centroids)]
    return record, candidates


def trajectory_row(state, sample, time_s, dt_s, index, event="ACTIVE"):
    row = dict(step=index, time_s=float(time_s), validation_dt_s=float(dt_s), timestep_role=VALIDATION_ROLE,
               particle_id=state.particle_id, radius_m=float(state.radius_m), tetra_id=int(sample.tetra_id),
               speed_m_s=float(np.linalg.norm(state.velocity_m_s)), pressure_pa=float(sample.pressure_pa),
               vorticity_norm_s_inv=float(np.linalg.norm(sample.vorticity_s_inv)), omega_norm_s_inv=float(np.linalg.norm(state.angular_velocity_s_inv)),
               boundary_event=event)
    for i, axis in enumerate("xyz"):
        row.update({f"{axis}_m": float(state.position_m[i]), f"u_{axis}_m_s": float(sample.velocity_m_s[i]),
                    f"V_{axis}_m_s": float(state.velocity_m_s[i]), f"vorticity_{axis}_s_inv": float(sample.vorticity_s_inv[i]),
                    f"Omega_{axis}_s_inv": float(state.angular_velocity_s_inv[i])})
    return row


def real_trajectory(field, classifier, initialization, radius_m, dt_s, progress=None):
    """Single trajectory, fixed dt; terminate at first boundary contact.

    Event position/time are the straight Euler segment/surface intersection,
    not a corrected step that restarts. All original endpoints are retained in
    the event record. WALL/INLET -> FAIL, no reflection or retry with a new dt.
    """
    state = MicrobubbleState(0, initialization["initial_position_m"], radius_m, np.zeros(3), np.zeros(3))
    sample = field.sample(state.position_m); state = equilibrium_state(state, sample)
    rows = [trajectory_row(state, sample, 0., dt_s, 0)]
    length, max_step, event_record = 0., 0., None
    max_steps = int(np.ceil(initialization["horizon_s"] / dt_s))
    segment_count = 0
    for step in range(1, max_steps + 1):
        predicted = euler_position(state.position_m, state.velocity_m_s, dt_s)
        event = classifier.first_event(state.position_m, predicted)
        segment_count += 1
        step_length = float(np.linalg.norm(predicted - state.position_m))
        max_step = max(max_step, step_length)
        if event is not None:
            time_s = ((step - 1) + event.segment_fraction) * dt_s
            event_record = dict(role=event.role, exit_time_s=float(time_s), position_m=event.position_m.tolist(),
                                segment_index=step - 1, segment_fraction=event.segment_fraction,
                                segment_start_m=state.position_m.tolist(), unmodified_trial_endpoint_m=predicted.tolist(),
                                triangle_id=event.triangle_id, role_triangle_id=event.role_triangle_id,
                                simultaneous_roles=list(event.simultaneous_roles))
            sample = field.sample(event.position_m)
            terminal = MicrobubbleState(0, event.position_m, radius_m, state.velocity_m_s, state.angular_velocity_s_inv)
            terminal = equilibrium_state(terminal, sample)
            rows.append(trajectory_row(terminal, sample, time_s, dt_s, step, event.role))
            length += event.segment_fraction * step_length
            break
        new_state = advance_single_microbubble(state, field, dt_s)
        # Separately sampled output supports row-wise V=u and Omega=curl/2 checks.
        sample = field.sample(new_state.position_m)
        rows.append(trajectory_row(new_state, sample, step * dt_s, dt_s, step))
        length += step_length
        state = new_state
        if progress and step % 500 == 0: progress(step, step * dt_s)
    numeric = [[v for v in row.values() if isinstance(v, (int, float))] for row in rows]
    finite = all(np.isfinite(values).all() for values in numeric)
    velocity_error = max(abs(r[f"V_{a}_m_s"] - r[f"u_{a}_m_s"]) for r in rows for a in "xyz")
    omega_error = max(abs(r[f"Omega_{a}_s_inv"] - .5 * r[f"vorticity_{a}_s_inv"]) for r in rows for a in "xyz")
    role = event_record["role"] if event_record else "HORIZON_REACHED"
    passed = finite and velocity_error == 0 and omega_error == 0 and role in ["OUTLET_01", "OUTLET_02", "OUTLET_03"]
    summary = dict(passed=passed, timestep_role=VALIDATION_ROLE, validation_dt_s=float(dt_s), active_and_terminal_finite=finite,
                   wall_crossing=role == "WALL", inlet_crossing=role == "INLET", exit_boundary=role,
                   exit_time_s=None if event_record is None else event_record["exit_time_s"],
                   trajectory_length_m=length, max_trial_step_length_m=max_step, row_count=len(rows),
                   every_segment_checked=segment_count == len(rows) - 1, checked_segment_count=segment_count,
                   max_velocity_relation_error_m_s=velocity_error, max_angular_relation_error_s_inv=omega_error,
                   event=event_record, production_particle_timestep_frozen=False,
                   visual_step_jump_review="PENDING_USER_REVIEW")
    return rows, summary
