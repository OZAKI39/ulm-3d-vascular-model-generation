"""Passive float64 GPU ledger of the unchanged F8 M0 integration path.

Production kernels are called unchanged. A read-only repeated wall query records
metadata, validated against the actual before/after projection vector. Buffers
are bounded, sparse wall events are compacted on GPU, and transfer occurs only
at accepted save intervals. No step-level CPU download or trajectory correction.
"""

from time import perf_counter
import numpy as np
import warp as wp
from .engine import FastSimulator
from .f8_gate import GateSimulator, GateFeedbackSimulator
from . import kernels as k
from . import vessel_kernels as vk

TERMS = [
    "flow",
    "primary",
    "secondary",
    "mobility",
    "pair_contact",
    "wall",
    "boundary",
    "accepted",
    "residual",
]
EVENT_COLUMNS = [
    "step",
    "start_s",
    "dt_s",
    "particle_id",
    "after_pair_round",
    "before_x_m",
    "before_y_m",
    "before_z_m",
    "after_x_m",
    "after_y_m",
    "after_z_m",
    "dx_m",
    "dy_m",
    "dz_m",
    "length_m",
    "gap_before_m",
    "normal_x",
    "normal_y",
    "normal_z",
    "wall_face_index",
    "preceding_pair_dx_m",
    "preceding_pair_dy_m",
    "preceding_pair_dz_m",
    "projection_reconstruction_error_m",
]


@wp.kernel
def begin_step(
    x: wp.array(dtype=wp.vec3d),
    old: wp.array(dtype=wp.vec3d),
    work: wp.array2d(dtype=wp.vec3d),
    norms: wp.array2d(dtype=wp.float64),
    pair_block: wp.array(dtype=wp.vec3d),
):
    i = wp.tid()
    old[i] = x[i]
    pair_block[i] = wp.vec3d(wp.float64(0.0))
    for j in range(9):
        work[i, j] = wp.vec3d(wp.float64(0.0))
    for j in range(3):
        norms[i, j] = wp.float64(0.0)


@wp.kernel
def free_terms(
    primary: wp.array(dtype=wp.vec3d),
    total: wp.array(dtype=wp.vec3d),
    flow: wp.array(dtype=wp.vec3d),
    free: wp.array(dtype=wp.vec3d),
    velocity: wp.array(dtype=wp.vec3d),
    parent: wp.array(dtype=int),
    delta: wp.array(dtype=wp.vec3d),
    count: wp.array(dtype=int),
    active: wp.array(dtype=int),
    drag: wp.float64,
    dt: wp.float64,
    flow_enabled: int,
    work: wp.array2d(dtype=wp.vec3d),
):
    i = wp.tid()
    if active[i] == 0:
        return
    if flow_enabled == 1:
        work[i, 0] = dt * flow[i]
    work[i, 1] = dt * (primary[i] / drag)
    work[i, 2] = dt * ((total[i] - primary[i]) / drag)
    root = k.root_of(parent, i)
    balanced = velocity[i] + delta[root] / wp.float64(count[root])
    work[i, 3] = dt * (balanced - free[i])


@wp.kernel
def pair_delta(
    before: wp.array(dtype=wp.vec3d),
    after: wp.array(dtype=wp.vec3d),
    work: wp.array2d(dtype=wp.vec3d),
    norms: wp.array2d(dtype=wp.float64),
    pair_block: wp.array(dtype=wp.vec3d),
):
    i = wp.tid()
    d = after[i] - before[i]
    work[i, 4] = work[i, 4] + d
    norms[i, 0] = norms[i, 0] + wp.length(d)
    pair_block[i] = pair_block[i] + d


@wp.kernel
def wall_metadata(
    closed: wp.uint64,
    wall: wp.uint64,
    vertices: wp.array(dtype=wp.vec3d),
    faces: wp.array2d(dtype=int),
    tags: wp.array(dtype=int),
    wall_faces: wp.array2d(dtype=int),
    origin: wp.vec3d,
    padding: wp.float64,
    x: wp.array(dtype=wp.vec3d),
    radius: wp.float64,
    active: wp.array(dtype=int),
    cache_x: wp.array(dtype=wp.vec3d),
    cache_d: wp.array(dtype=wp.float64),
    error: wp.array(dtype=int),
    gap: wp.array(dtype=wp.float64),
    normal: wp.array(dtype=wp.vec3d),
    face: wp.array(dtype=int),
):
    i = wp.tid()
    face[i] = -1
    gap[i] = wp.float64(0.0)
    normal[i] = wp.vec3d(wp.float64(0.0))
    if active[i] == 0 or error[0] != 0:
        return
    lower = cache_d[i] - wp.length(x[i] - cache_x[i])
    if lower >= radius + padding:
        return
    h = vk.wall_query(closed, wall, vertices, faces, tags, wall_faces, origin, padding, x[i])
    if h.valid == 1:
        face[i] = h.face
        gap[i] = h.distance - radius
        normal[i] = h.normal


@wp.kernel
def wall_delta(
    before: wp.array(dtype=wp.vec3d),
    after: wp.array(dtype=wp.vec3d),
    work: wp.array2d(dtype=wp.vec3d),
    norms: wp.array2d(dtype=wp.float64),
    pair_block: wp.array(dtype=wp.vec3d),
    gap: wp.array(dtype=wp.float64),
    normal: wp.array(dtype=wp.vec3d),
    face: wp.array(dtype=int),
    index: wp.array(dtype=int),
    offset: wp.array(dtype=int),
    clock: wp.array(dtype=int),
    origin: wp.array(dtype=wp.float64),
    dt: wp.float64,
    round_: int,
    event_count: wp.array(dtype=int),
    events: wp.array2d(dtype=wp.float64),
    audit_error: wp.array(dtype=int),
):
    i = wp.tid()
    d = after[i] - before[i]
    length = wp.length(d)
    work[i, 5] = work[i, 5] + d
    norms[i, 1] = norms[i, 1] + length
    reconstruction = wp.length(d - wp.min(gap[i], wp.float64(0.0)) * normal[i])
    if length > wp.float64(0.0):
        norms[i, 2] = wp.max(norms[i, 2], reconstruction)
        row = wp.atomic_add(event_count, 0, 1)
        if row >= events.shape[0]:
            wp.atomic_or(audit_error, 0, 1)
        else:
            events[row, 0] = wp.float64(offset[0] + index[0])
            events[row, 1] = origin[0] + wp.float64(clock[0]) * dt
            events[row, 2] = dt
            events[row, 3] = wp.float64(i)
            events[row, 4] = wp.float64(round_)
            for axis in range(3):
                events[row, 5 + axis] = before[i][axis]
                events[row, 8 + axis] = after[i][axis]
                events[row, 11 + axis] = d[axis]
                events[row, 16 + axis] = normal[i][axis]
                events[row, 20 + axis] = pair_block[i][axis]
            events[row, 14] = length
            events[row, 15] = gap[i]
            events[row, 19] = wp.float64(face[i])
            events[row, 23] = reconstruction
    pair_block[i] = wp.vec3d(wp.float64(0.0))


@wp.kernel
def finish_step(
    old: wp.array(dtype=wp.vec3d),
    proposal: wp.array(dtype=wp.vec3d),
    accepted: wp.array(dtype=wp.vec3d),
    work: wp.array2d(dtype=wp.vec3d),
    norms: wp.array2d(dtype=wp.float64),
    n: int,
    clock: wp.array(dtype=int),
    origin: wp.array(dtype=wp.float64),
    dt: wp.float64,
    error: wp.array(dtype=int),
    diag: wp.array2d(dtype=wp.float64),
    radius: wp.float64,
    index: wp.array(dtype=int),
    rows: wp.array2d(dtype=wp.vec3d),
    meta: wp.array2d(dtype=wp.float64),
    cumulative: wp.array2d(dtype=wp.vec3d),
    audit_error: wp.array(dtype=int),
):
    row = index[0]
    if row >= rows.shape[0]:
        wp.atomic_or(audit_error, 0, 2)
        return
    valid = int(error[0] == 0)
    for i in range(n):
        if diag[i, 3] > wp.float64(0.01) * radius or diag[i, 6] > wp.float64(0.1) * radius:
            valid = 0
    meta[row, 0] = origin[0] + wp.float64(clock[0]) * dt
    meta[row, 1] = dt
    meta[row, 2] = wp.float64(valid)
    for j in range(3, 7):
        meta[row, j] = wp.float64(0.0)
    for j in range(9):
        rows[row, j] = wp.vec3d(wp.float64(0.0))
    if valid == 1:
        for i in range(n):
            work[i, 6] = accepted[i] - proposal[i]
            work[i, 7] = accepted[i] - old[i]
            s = wp.vec3d(wp.float64(0.0))
            for j in range(7):
                s = s + work[i, j]
            work[i, 8] = work[i, 7] - s
            for j in range(9):
                rows[row, j] = rows[row, j] + work[i, j] / wp.float64(n)
                cumulative[i, j] = cumulative[i, j] + work[i, j]
            meta[row, 3] = wp.max(meta[row, 3], wp.length(work[i, 8]))
            meta[row, 4] = meta[row, 4] + norms[i, 0] / wp.float64(n)
            meta[row, 5] = meta[row, 5] + norms[i, 1] / wp.float64(n)
            meta[row, 6] = wp.max(meta[row, 6], norms[i, 2])
    index[0] = row + 1


class AuditMixin:
    def __init__(self, *args, audit_enabled=False, audit_capacity=256, **kwargs):
        self.audit_ready = False
        self.audit_enabled = bool(audit_enabled)
        super().__init__(*args, **kwargs)
        if not self.audit_enabled:
            return
        if self.batch != 1 or self.mobility_mode != "f5_local" or self.unit_velocity is not None:
            raise ValueError("F9 audit is admitted only for original M0 single-group force/flow dynamics")
        self.audit_capacity = int(audit_capacity)
        self.audit_old = wp.empty_like(self.x)
        self.audit_primary = wp.empty_like(self.force)
        self.audit_pair_before = wp.empty_like(self.x)
        self.audit_wall_before = wp.empty_like(self.x)
        self.audit_proposal = wp.empty_like(self.x)
        self.audit_pair_block = wp.zeros_like(self.x)
        self.audit_work = wp.zeros((self.n, 9), dtype=wp.vec3d, device=self.device)
        self.audit_norms = wp.zeros((self.n, 3), dtype=wp.float64, device=self.device)
        self.audit_gap = wp.zeros(self.n, dtype=wp.float64, device=self.device)
        self.audit_normal = wp.zeros_like(self.x)
        self.audit_face = wp.zeros(self.n, dtype=int, device=self.device)
        self.audit_index = wp.zeros(1, dtype=int, device=self.device)
        self.audit_offset = wp.zeros(1, dtype=int, device=self.device)
        self.audit_error = wp.zeros(1, dtype=int, device=self.device)
        self.audit_event_count = wp.zeros(1, dtype=int, device=self.device)
        self.audit_rows = wp.zeros((self.audit_capacity, 9), dtype=wp.vec3d, device=self.device)
        self.audit_meta = wp.zeros((self.audit_capacity, 7), dtype=wp.float64, device=self.device)
        self.audit_events = wp.zeros(
            (self.audit_capacity * self.n * 5, len(EVENT_COLUMNS)), dtype=wp.float64, device=self.device
        )
        self.audit_cumulative = wp.zeros((self.n, 9), dtype=wp.vec3d, device=self.device)
        self.audit_chunks, self.audit_event_chunks, self.audit_rejected = [], [], []
        self.audit_transfer_s = 0.0
        self.audit_steps = 0
        self.audit_ready = True
        wp.load_module(module=begin_step.module, device=self.device)

    def reset(self):
        super().reset()
        if self.audit_ready:
            for a in [
                self.audit_index,
                self.audit_offset,
                self.audit_error,
                self.audit_event_count,
                self.audit_cumulative,
            ]:
                a.zero_()
            self.audit_chunks, self.audit_event_chunks, self.audit_rejected = [], [], []
            self.audit_steps = 0
            self.audit_transfer_s = 0.0

    def _audit_launch(self, kernel, dim, args):
        FastSimulator.launch(self, kernel, dim, args)

    def fields(self):
        super().fields()
        if self.audit_ready:
            wp.copy(self.audit_primary, self.force)

    def step(self, dt):
        if self.audit_ready:
            self._audit_launch(
                begin_step,
                self.n,
                [self.x, self.audit_old, self.audit_work, self.audit_norms, self.audit_pair_block],
            )
        return super().step(dt)

    def _record_pair(self, out):
        self._audit_launch(
            pair_delta,
            self.n,
            [self.audit_pair_before, out, self.audit_work, self.audit_norms, self.audit_pair_block],
        )
        self._pair_recorded = True

    def wall_pass(self, x):
        if not self.audit_ready:
            return super().wall_pass(x)
        if self._audit_stage == "pair":
            self._record_pair(x)
        wp.copy(self.audit_wall_before, x)
        self._audit_launch(
            wall_metadata,
            self.n,
            [
                *self.geometry.args(),
                x,
                self.p.radius_m,
                self.active,
                self.cache_x,
                self.cache_d,
                self.error,
                self.audit_gap,
                self.audit_normal,
                self.audit_face,
            ],
        )
        result = super().wall_pass(x)
        self._audit_launch(
            wall_delta,
            self.n,
            [
                self.audit_wall_before,
                x,
                self.audit_work,
                self.audit_norms,
                self.audit_pair_block,
                self.audit_gap,
                self.audit_normal,
                self.audit_face,
                self.audit_index,
                self.audit_offset,
                self.step_count,
                self.time_origin,
                self._current_dt,
                self._pair_contact_pass,
                self.audit_event_count,
                self.audit_events,
                self.audit_error,
            ],
        )
        return result

    def launch(self, kernel, dim, args):
        if not self.audit_ready:
            return super().launch(kernel, dim, args)
        if kernel is k.predict:
            self._audit_stage = "predict"
            self._audit_launch(
                free_terms,
                self.n,
                [
                    self.audit_primary,
                    self.force,
                    self.flow_query.velocity,
                    self.free,
                    args[1],
                    self.component_parent,
                    self.component_delta,
                    self.component_count,
                    self.active,
                    6 * np.pi * self.p.viscosity_Pa_s * self.p.radius_m,
                    args[3],
                    int(self.flow_enabled),
                    self.audit_work,
                ],
            )
        elif kernel is k.contact:
            self._audit_stage = "pair"
            self._pair_recorded = False
            wp.copy(self.audit_pair_before, args[0])
        elif kernel is k.accept:
            self._audit_stage = "accept"
            wp.copy(self.audit_proposal, args[0])
        result = super().launch(kernel, dim, args)
        if kernel is k.contact and not self._pair_recorded:
            self._record_pair(args[1])
        elif kernel is k.accept:
            self._audit_launch(
                finish_step,
                1,
                [
                    self.audit_old,
                    self.audit_proposal,
                    self.x,
                    self.audit_work,
                    self.audit_norms,
                    self.n,
                    self.step_count,
                    self.time_origin,
                    self._current_dt,
                    self.error,
                    self.diag,
                    self.p.radius_m,
                    self.audit_index,
                    self.audit_rows,
                    self.audit_meta,
                    self.audit_cumulative,
                    self.audit_error,
                ],
            )
        return result

    def _flush_audit(self, rejected=False):
        begin = perf_counter()
        count = int(self.audit_index.numpy()[0])
        event_count = int(self.audit_event_count.numpy()[0])
        code = int(self.audit_error.numpy()[0])
        if code or count > self.audit_capacity or event_count > self.audit_events.shape[0]:
            raise RuntimeError("F9_AUDIT_BUFFER_OVERFLOW: no silent truncation")
        rows = self.audit_rows.numpy()[:count].copy()
        meta = self.audit_meta.numpy()[:count].copy()
        host = wp.empty((event_count, len(EVENT_COLUMNS)), dtype=wp.float64, device="cpu", pinned=True)
        if event_count:
            wp.copy(host, self.audit_events, count=event_count * len(EVENT_COLUMNS))
            wp.synchronize_device(self.device)
        events = host.numpy().copy()
        if len(events):
            events = events[np.lexsort((events[:, 3], events[:, 4], events[:, 0]))]
        record = (rows, meta)
        if rejected or not np.all(meta[:, 2] == 1):
            self.audit_rejected.append(dict(rows=rows, meta=meta, events=events))
        else:
            self.audit_chunks.append(record)
            self.audit_event_chunks.append(events)
        self.audit_steps += count
        self.audit_offset.assign(np.array([self.audit_steps], np.int32))
        self.audit_index.zero_()
        self.audit_event_count.zero_()
        self.audit_transfer_s += perf_counter() - begin

    def advance(self, dt, count, graphs=True):
        if not self.audit_ready:
            return super().advance(dt, count, graphs)
        if count > self.audit_capacity:
            raise ValueError("Audit save interval exceeds registered bounded capacity")
        try:
            result = super().advance(dt, count, graphs)
        except Exception:
            self._flush_audit(rejected=True)
            raise
        self._flush_audit()
        return result

    def run(self, *args, **kwargs):
        if self.audit_ready and kwargs.get("adaptive", False):
            raise ValueError("Adaptive rollback not admitted for this fixed-step audit")
        result = super().run(*args, **kwargs)
        if self.audit_ready:
            result["motion_ledger"] = dict(
                displacement_m=np.concatenate([a[0] for a in self.audit_chunks])
                if self.audit_chunks
                else np.empty((0, 9, 3)),
                step_meta=np.concatenate([a[1] for a in self.audit_chunks])
                if self.audit_chunks
                else np.empty((0, 7)),
                wall_contact_events=np.concatenate(self.audit_event_chunks)
                if self.audit_event_chunks
                else np.empty((0, len(EVENT_COLUMNS))),
                per_ID_cumulative_m=self.audit_cumulative.numpy(),
                rejected_intervals=self.audit_rejected,
            )
            result["summary"]["audit"] = dict(
                enabled=True,
                transfer_and_reduction_readback_s=self.audit_transfer_s,
                accepted_steps=len(result["motion_ledger"]["step_meta"]),
                rejected_intervals=len(self.audit_rejected),
                step_buffer_capacity=self.audit_capacity,
                event_buffer_capacity=self.audit_events.shape[0],
                scope="Exact signed before/after production updates; a same-trajectory numerical identity, not a no-wall counterfactual. GPU reductions and sparse events; batch transfer at save intervals only.",
            )
        else:
            result["summary"]["audit"] = dict(enabled=False)
        return result


class AuditSimulator(AuditMixin, GateSimulator):
    pass


class AuditFeedbackSimulator(AuditMixin, GateFeedbackSimulator):
    pass
