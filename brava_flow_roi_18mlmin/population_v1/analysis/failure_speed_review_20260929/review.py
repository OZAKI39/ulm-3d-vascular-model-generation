"""Read saved BraVa data; write derived diagnostics here; never integrate a particle."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
POP = OUT.parent.parent
BASE = POP.parent
TRACKS = POP / 'analysis/realization_0000_contactfix/tracks'
RUN = POP / 'runs/realization_0000_contactfix'
INPUTS = {}


def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def register(path):
    INPUTS[str(path)] = digest(path)
    return path


def read_json(path):
    return json.loads(register(path).read_text())


def save_json(name, data):
    (OUT / name).write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def track(pid):
    with np.load(register(TRACKS / f'mb_{pid:06d}.npz')) as z:
        return {k: z[k] for k in z.files}


flow_path = register(BASE / 'cases/balanced_pressure_final/frozen_flow/steady_flow.vtu')
identity = read_json(RUN / 'IDENTITY.json')
assert INPUTS[str(flow_path)] == identity['flow_sha256']
with np.load(register(BASE / 'cases/balanced_pressure_final/frozen_flow/flow_arrays_si.npz')) as z:
    points, tetra, velocity = z['points_m'], z['tetra'], z['velocity_m_s']
active = read_json(BASE / 'reports/ACTIVE_FLOW.json')
assert INPUTS[str(BASE / 'cases/balanced_pressure_final/frozen_flow/flow_arrays_si.npz')] == active['postprocess_data_sha256']['frozen_flow/flow_arrays_si.npz']


def geometry(cid):
    q = points[tetra[cid]]
    a = (q[1:] - q[0]).T
    inv = np.linalg.inv(a)
    uv = velocity[tetra[cid]]
    grad = (uv[1:] - uv[0]).T @ inv
    scale = np.linalg.norm(a, np.inf)
    tau = 64 * np.finfo(float).eps * scale * np.linalg.norm(inv, np.inf) * (1 + np.abs(q).max() / scale)
    return q, uv, inv, grad, tau


def sample(cid, x):
    q, uv, inv, grad, tau = geometry(cid)
    tail = inv @ (x - q[0])
    weights = np.r_[1 - tail.sum(), tail]
    return weights, weights @ uv, grad, tau


diagnosed = read_json(POP / 'analysis/realization_0000_contactfix/failure_diagnostics.json')
failures = []
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
fig, axes = plt.subplots(2, 2, figsize=(12, 7.4), constrained_layout=True)
for col, pid in enumerate([440, 1505]):
    z = track(pid)
    t, x, v, cell = z['physical_time_s'], z['x_m'], z['velocity_m_s'], z['tetra_id']
    first = int(np.flatnonzero(cell == 138584)[0])
    prior = int(cell[first - 1])
    weights, terminal_u, grad, tau = sample(138584, x[-1])
    prior_w, prior_u, prior_g, prior_tau = sample(prior, x[-1])
    q, uv, inv, _, _ = geometry(prior)
    nonzero = np.flatnonzero(np.linalg.norm(uv, axis=1) > 0)
    assert len(nonzero) == 1 and prior == 578387
    k = int(nonzero[0])
    shape_grad = np.vstack([-inv.sum(axis=0), inv])[k]
    interface_distance = abs(prior_w[k]) / np.linalg.norm(shape_grad)
    diag = next(r for r in diagnosed if r['particle_id'] == pid)
    row = dict(particle_id=pid, terminal_time_s=float(t[-1]),
               first_zero_cell_assignment_s=float(t[first]), prior_cell=prior,
               zero_cell=138584, zero_cell_sample_count=int(np.count_nonzero(cell == 138584)),
               zero_velocity_sample_count=int(np.count_nonzero(np.linalg.norm(v, axis=1) == 0)),
               position_change_after_zero_cell_assignment_m=float(np.linalg.norm(x[-1] - x[first])),
               terminal_gap_um=float(z['wall_gap_m'][-1] * 1e6),
               contact_rank=diag['contact_ranks'], trigger=diag['trigger'],
               zero_cell_terminal_weights=weights.tolist(), zero_cell_weight_tolerance=float(tau),
               within_zero_cell_tolerance=bool(weights.min() >= -tau and weights.max() <= 1 + tau),
               prior_cell_terminal_weights=prior_w.tolist(),
               distance_to_shared_zero_velocity_face_m=float(interface_distance),
               terminal_speed_using_prior_cell_m_s=float(np.linalg.norm(prior_u)),
               prior_cell_divergence_s_inv=float(np.trace(prior_g)),
               prior_cell_linear_decay_timescale_s=float(-1 / np.trace(prior_g)),
               interpretation='Approach to an internal zero-velocity face; final cell assignment is within roundoff tolerance, not macroscopic penetration')
    failures.append(row)
    rel = (t[:-1] - t[first]) * 1000
    speed = np.linalg.norm(v[1:], axis=1) * 1000  # held velocity evaluated at preceding position
    keep = (rel >= -65) & (rel <= 18)
    ax = axes[0, col]
    ax.plot(rel[keep], speed[keep], color='#246ca6', lw=1.4)
    ax.set_yscale('symlog', linthresh=1e-10)
    ax.axvline(0, color='#ba4a38', ls='--', label='First assignment to cell 138584')
    ax.set(title=f'Failure ID {pid}: recorded speed', ylabel='Speed (mm/s)', xlabel='Time from zero-cell assignment (ms)')
    ax.grid(alpha=.2)
    ax.legend(fontsize=8)
    ax = axes[1, col]
    ax.plot(rel[keep], z['wall_gap_m'][:-1][keep] * 1e6, color='#246ca6')
    ax.axvline(0, color='#ba4a38', ls='--')
    ax.axhline(.002, color='#777777', ls=':', label='2 nm handoff scale')
    ax.set(title='Finite distance from the actual wall', ylabel='Sphere-surface wall gap (micrometers)', xlabel='Time from zero-cell assignment (ms)')
    ax.grid(alpha=.2)
    ax.legend(fontsize=8)
fig.suptitle('Two numerical failures: stagnation near an internal zero-velocity face\nSaved historical trajectories; no new integration', fontsize=13)
fig.savefig(OUT / 'failures.png', dpi=300)
fig.savefig(OUT / 'failures.pdf')
plt.close(fig)
save_json('failure_details.json', failures)

current_differences = []
for name, expected in identity['code_files_sha256'].items():
    path = BASE / 'code/particle_core' / name
    if digest(path) != expected:
        current_differences.append(name)

scan = read_json(POP / 'analysis/acceleration_diagnosis/SCAN_SUMMARY.json')
diagnosis = read_json(POP / 'analysis/acceleration_diagnosis/DIAGNOSIS_SUMMARY.json')
candidate_path = register(POP / 'analysis/acceleration_diagnosis/speedup_diagnosis_all_candidates.csv')
with candidate_path.open() as f:
    candidates = list(csv.DictReader(f))
overlap_path = register(POP / 'cfd_wall_audit/acceleration_overlap.csv')
with overlap_path.open() as f:
    overlaps = list(csv.DictReader(f))
animation_path = register(POP / 'candidate_0/views_0to1p2s_large_axes/az240_el25p0/data/animation_samples.npz')
with np.load(animation_path) as z:
    ft, ids, pos, means, displayed, prior = (z[k] for k in ['frame_time_s', 'particle_id', 'position_m', 'mean_speed_m_s', 'displayed', 'prior_exclusion'])
alive = np.isfinite(pos).all(axis=2)
assert np.array_equal(displayed, alive & ~prior[None, :] & (means >= .002))
on = ~displayed[:-1] & displayed[1:] & alive[:-1] & alive[1:]
off = displayed[:-1] & ~displayed[1:] & alive[:-1] & alive[1:]
early = []
for r in candidates:
    if float(r['physical_time_s']) > 1.2:
        continue
    pid = int(r['particle_id'])
    j = int(np.flatnonzero(ids == pid)[0])
    early.append(dict(particle_id=pid, physical_time_s=float(r['physical_time_s']),
                      before_mm_s=float(r['before_m_s']) * 1000,
                      after_mm_s=float(r['after_m_s']) * 1000,
                      before_gap_um=float(r['before_gap_median_m']) * 1e6,
                      after_gap_um=float(r['after_gap_median_m']) * 1e6,
                      pair_rows_in_window=int(r['pair_rows_in_window']),
                      permanently_excluded=bool(prior[j]),
                      threshold_reappearance_times_s=ft[1:][on[:, j]].tolist()))
visibility = dict(frame_count=len(ft), physical_frame_interval_s=float(np.diff(ft).mean()),
                  threshold_reappearance_transitions_while_alive=int(on.sum()),
                  threshold_disappearance_transitions_while_alive=int(off.sum()),
                  affected_particle_count=int(np.count_nonzero(on.any(axis=0) | off.any(axis=0))),
                  exact_mask_reconstruction_pass=True,
                  interpretation='Transitions exclude birth/death between frames and count visibility changes, not physical jumps',
                  early_existing_candidates=early)
save_json('display_review.json', visibility)

fig, axes = plt.subplots(2, 2, figsize=(12, 7.4), constrained_layout=True)
for col, pid in enumerate([1189, 642]):
    z = track(pid)
    event = next(r for r in candidates if int(r['particle_id']) == pid)
    center = float(event['physical_time_s'])
    ts = z['physical_time_s'][:-1]
    keep = (ts >= center - .12) & (ts <= center + .12)
    indices = np.flatnonzero(keep)
    field_speeds = np.array([np.linalg.norm(sample(int(z['tetra_id'][i]), z['x_m'][i])[1]) for i in indices]) * 1000
    ax = axes[0, col]
    ax.plot(ts[keep], np.linalg.norm(z['velocity_m_s'][1:][keep], axis=1) * 1000, label='Recorded bubble speed', color='#246ca6')
    ax.plot(ts[keep], field_speeds, label='P1 fluid speed at same start position', color='#e59831', ls='--')
    ax.set_yscale('log')
    ax.set(title=f'ID {pid}: near-wall motion and release', xlabel='Physical time (s)', ylabel='Speed (mm/s)')
    ax.grid(alpha=.2)
    ax.legend(fontsize=8)
    ax = axes[1, col]
    if pid == 1189:
        ax.plot(ts[keep], z['wall_gap_m'][:-1][keep] * 1e6, color='#246ca6')
        ax.set_yscale('log')
        ax.axhline(float(event['radius_m']) * 1e6, color='#777777', ls=':', label='One particle radius')
        ax.set(title='Wall gap grows during release', xlabel='Physical time (s)', ylabel='Wall gap (micrometers)')
    else:
        j = int(np.flatnonzero(ids == pid)[0])
        fk = (ft >= center - .15) & (ft <= min(1.2, center + .12))
        ax.plot(ft[fk], means[fk, j] * 1000, color='#246ca6', label='Rolling 0.1 s mean speed')
        ax.axhline(2, color='#ba4a38', ls='--', label='Display threshold: 2 mm/s')
        ax.scatter(ft[fk & displayed[:, j]], means[fk & displayed[:, j], j] * 1000, s=13, color='#246ca6', label='Shown frames', zorder=3)
        ax.scatter(ft[fk & ~displayed[:, j]], means[fk & ~displayed[:, j], j] * 1000, s=16, marker='x', color='#ba4a38', label='Hidden frames', zorder=3)
        ax.set_yscale('log')
        ax.set(title='Latest animation: threshold controls visibility', xlabel='Physical time (s)', ylabel='Rolling mean speed (mm/s)')
    ax.grid(alpha=.2)
    ax.legend(fontsize=8)
fig.suptitle('Speed changes exist in saved trajectories; display filtering adds appearance changes\nID 1189 is excluded from the latest display; ID 642 is retained', fontsize=12)
fig.savefig(OUT / 'speed_and_visibility.png', dpi=300)
fig.savefig(OUT / 'speed_and_visibility.pdf')
plt.close(fig)

# Locate the strongest increments in the existing accepted-state data; no state updates.
extremes = []
double_jump_tracks = 0
for n, path in enumerate(sorted(TRACKS.glob('mb_*.npz'))):
    with np.load(register(path)) as z:
        s = np.linalg.norm(z['velocity_m_s'], axis=1)
        delta = np.diff(s)
        if len(delta) < 2:
            continue
        qualifying = (delta > .005) & (s[1:] > 2 * np.maximum(s[:-1], 1e-14))
        qualifying[0] = False  # birth row has different velocity semantics
        double_jump_tracks += int(qualifying.any())
        i = int(np.argmax(delta[1:])) + 2
        before_cell, after_cell = int(z['tetra_id'][i - 2]), int(z['tetra_id'][i - 1])
        extremes.append(dict(particle_id=int(path.stem[3:]),
            evaluation_time_s=float(z['physical_time_s'][i - 1]),
            evaluation_dt_s=float(z['physical_time_s'][i - 1] - z['physical_time_s'][i - 2]),
            before_mm_s=float(s[i - 1] * 1000), after_mm_s=float(s[i] * 1000),
            increase_mm_s=float(delta[i - 1] * 1000),
            before_tetra=before_cell, after_tetra=after_cell,
            before_gap_um=float(z['wall_gap_m'][i - 2] * 1e6),
            after_gap_um=float(z['wall_gap_m'][i - 1] * 1e6)))
    if n % 400 == 0:
        print(f'Read {n} existing tracks', flush=True)
extremes.sort(key=lambda r: -r['increase_mm_s'])
strongest = extremes[0]
z = track(strongest['particle_id'])
i = int(np.flatnonzero(z['physical_time_s'] == strongest['evaluation_time_s'])[0])
context = []
for j in [i - 1, i]:
    cid = int(z['tetra_id'][j])
    _, u, grad, _ = sample(cid, z['x_m'][j])
    context.append(dict(time_s=float(z['physical_time_s'][j]), cell=cid,
        bubble_speed_mm_s=float(np.linalg.norm(z['velocity_m_s'][j + 1]) * 1000),
        FEM_speed_mm_s=float(np.linalg.norm(u) * 1000),
        velocity_error_m_s=float(np.max(abs(z['velocity_m_s'][j + 1] - u))),
        divergence_s_inv=float(np.trace(grad)), gradient_norm_s_inv=float(np.linalg.norm(grad))))
save_json('largest_increment_context.json', dict(particle_id=strongest['particle_id'], states=context,
    displacement_um=float(np.linalg.norm(z['x_m'][i] - z['x_m'][i - 1]) * 1e6),
    fractional_speed_increase=context[1]['bubble_speed_mm_s'] / context[0]['bubble_speed_mm_s'] - 1,
    interpretation='Fast bulk-flow interval; held bubble velocities match local P1 samples. This does not certify time accuracy.'))
save_json('largest_step_increments.json', dict(track_count=len(extremes),
    tracks_with_increment_gt5_mm_s_and_speed_more_than_doubled=double_jump_tracks,
    top_per_track=extremes[:20],
    limitation='Threshold-specific screen, not a proof of temporal accuracy; held velocities aligned to their evaluation positions'))

save_json('SUMMARY.json', dict(
    run='realization_0000_contactfix', protocol='historical_outcome_dependent_stop',
    flow_sha256=identity['flow_sha256'], current_vs_run_source_differences=current_differences,
    failures=failures, visibility=visibility,
    previous_acceleration_screen_count=len(candidates),
    previous_screen_no_pair_in_window=int(diagnosis['no_pair_during_speedup']),
    previous_screen_high_divergence_overlap=sum(r['window_any_eta_gt_0p5'] == 'True' for r in overlaps),
    temporal_screen_tracks=len(extremes), temporal_screen_double_jump_tracks=double_jump_tracks,
    strongest_step_increment=extremes[0],
    cfd_solved=False, trajectories_integrated=False, scientific_inputs_modified=False,
    scope='Postprocessing stored data and inspecting source; no convergence study or new physical validation'))
assert all(digest(Path(path)) == expected for path, expected in INPUTS.items())
save_json('INPUT_HASHES.json', INPUTS)
print(json.dumps(dict(output=str(OUT), inputs_checked=len(INPUTS), largest_increment=extremes[0], display_transitions=visibility['affected_particle_count']), ensure_ascii=False), flush=True)
