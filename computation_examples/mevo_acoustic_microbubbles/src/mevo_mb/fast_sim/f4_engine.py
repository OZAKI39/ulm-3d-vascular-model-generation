"""Opt-in F4 vessel constraints; disabled vessel delegates the archived F3 path."""

from time import perf_counter
import hashlib
from pathlib import Path
import numpy as np
import warp as wp
from .f3_engine import F3Simulator
from . import kernels as k
from . import f2_kernels as f2
from . import f3_kernels as f3
from . import vessel_kernels as vk
from . import vessel_active_kernels as masked
from .vessel_geometry import ORIGINAL_SHA256


class GPUVessel:
    def __init__(self, prepared, original, device="cuda:0"):
        begin = perf_counter()
        if hashlib.sha256(Path(original).read_bytes()).hexdigest() != ORIGINAL_SHA256:
            raise ValueError("Original labelled Y mesh identity mismatch")
        a = np.load(prepared)
        assert str(a["original_sha256"]) == ORIGINAL_SHA256
        self.device = device
        self.origin = wp.vec3d(*a["local_origin_m"])
        self.padding = float(a["broadphase_padding_m"])
        self.vertices = wp.array(a["vertices_m"], dtype=wp.vec3d, device=device)
        self.faces = wp.array(a["faces"], dtype=wp.int32, device=device)
        self.tags = wp.array(a["labels"], dtype=wp.int32, device=device)
        wf = a["faces"][a["labels"] == 1]
        self.wall_faces = wp.array(wf, dtype=wp.int32, device=device)
        self.points32 = wp.array(
            (a["vertices_m"] - a["local_origin_m"]).astype(np.float32), dtype=wp.vec3, device=device
        )
        self.indices = wp.array(a["faces"].ravel(), dtype=wp.int32, device=device)
        self.wall_indices = wp.array(wf.ravel(), dtype=wp.int32, device=device)
        self.closed = wp.Mesh(points=self.points32, indices=self.indices)
        self.wall = wp.Mesh(points=self.points32, indices=self.wall_indices)
        self.cap_origins = wp.array(a["cap_origins_m"], dtype=wp.vec3d, device=device)
        self.cap_normals = wp.array(a["cap_normals"], dtype=wp.vec3d, device=device)
        self.build_s = perf_counter() - begin
        self.metadata = dict(
            original_sha256=ORIGINAL_SHA256,
            prepared_sha256=hashlib.sha256(Path(prepared).read_bytes()).hexdigest(),
            mesh_build_s=self.build_s,
            static_BVH=True,
            wall_triangle_count=len(wf),
            units="m",
            broadphase_padding_m=self.padding,
        )

    def args(self):
        return [
            self.closed.id,
            self.wall.id,
            self.vertices,
            self.faces,
            self.tags,
            self.wall_faces,
            self.origin,
            self.padding,
        ]

    def query(self, x):
        # Read-only initialization/offline adapter. Physical kernels use persistent buffers.
        gpu = wp.array(np.asarray(x).reshape(-1, 3), dtype=wp.vec3d, device=self.device)
        d = wp.empty(len(x), dtype=wp.float64, device=self.device)
        n = wp.empty(len(x), dtype=wp.vec3d, device=self.device)
        f = wp.empty(len(x), dtype=wp.int32, device=self.device)
        inside = wp.empty(len(x), dtype=wp.int32, device=self.device)
        wp.launch(
            vk.query_points, dim=len(x), inputs=[*self.args(), gpu, d, n, f, inside], device=self.device
        )
        result = dict(distance_m=d.numpy(), normal=n.numpy(), face=f.numpy(), inside=inside.numpy())
        if (result["face"] < 0).any():
            raise RuntimeError("GPU wall query failed")
        return result


class F4VesselSimulator(F3Simulator):
    def __init__(
        self,
        *args,
        vessel=None,
        vessel_enabled=True,
        wall_contact_enabled=True,
        geometry_unit_test_velocity=None,
        **kwargs,
    ):
        self.vessel_enabled = bool(vessel_enabled)
        self.wall_contact_enabled = bool(wall_contact_enabled)
        self.geometry = vessel
        self._vessel_ready = False
        self._wall_frames = None
        self._active_frames = None
        self._wall_stats_frames = None
        self.unit_velocity = None
        super().__init__(*args, **kwargs)
        if not self.vessel_enabled:
            return
        if vessel is None:
            raise ValueError("Vessel enabled requires original GPU geometry")
        self.active = wp.ones(self.total, dtype=wp.int32, device=self.device)
        self.pending = wp.zeros(self.total, dtype=wp.int32, device=self.device)
        self.fractions = wp.zeros(self.total, dtype=wp.float64, device=self.device)
        self.exit_label = wp.zeros(self.total, dtype=wp.int32, device=self.device)
        self.exit_time = wp.zeros(self.total, dtype=wp.float64, device=self.device)
        self.wall_stats = wp.zeros((self.total, 12), dtype=wp.float64, device=self.device)
        self.cache_x = wp.empty_like(self.x)
        self.cache_d = wp.empty(self.total, dtype=wp.float64, device=self.device)
        self.wall_gap = wp.empty(self.total, dtype=wp.float64, device=self.device)
        self.query_d = wp.empty_like(self.wall_gap)
        self.query_n = wp.empty_like(self.x)
        self.query_face = wp.empty(self.total, dtype=wp.int32, device=self.device)
        self.query_inside = wp.empty_like(self.query_face)
        self.field_enabled = int(self.basis is not None)
        if self.field_enabled:
            margin = self.p.radius_m + 2e-8
            self.field_lower = wp.vec3d(*(self.basis.origin + margin))
            self.field_upper = wp.vec3d(
                *(self.basis.origin + self.basis.spacing * (np.array(self.basis.shape) - 1) - margin)
            )
            x = self.initial.reshape(-1, 3)
            if (x < np.array(self.field_lower)).any() or (x > np.array(self.field_upper)).any():
                raise ValueError("Initial particle outside radius-inset stored field")
        else:
            self.field_lower = wp.vec3d(-1.0)
            self.field_upper = wp.vec3d(1.0)
        if geometry_unit_test_velocity is not None:
            self.unit_velocity = wp.array(
                np.broadcast_to(geometry_unit_test_velocity, (self.total, 3)).copy(),
                dtype=wp.vec3d,
                device=self.device,
            )
        self._vessel_ready = True
        self.reset()
        begin = perf_counter()
        wp.load_module(module=vk, device=self.device)
        wp.load_module(module=masked, device=self.device)
        wp.synchronize_device(self.device)
        self.compile_s += perf_counter() - begin

    def reset(self):
        super().reset()
        if not self._vessel_ready:
            return
        self.active.fill_(1)
        for a in [self.pending, self.exit_label, self.exit_time, self.wall_stats]:
            a.zero_()
        q = self.geometry.query(self.initial.reshape(-1, 3))
        if not q["inside"].all() or (q["distance_m"] < self.p.radius_m * (1 - 0.01)).any():
            raise ValueError("Initial state violates original lumen/wall radius budget")
        self.cache_x.assign(self.initial.reshape(-1, 3))
        self.cache_d.assign(q["distance_m"])
        self.wall_gap.assign(q["distance_m"] - self.p.radius_m)
        d = np.zeros((self.total, 12))
        d[:, 0] = q["distance_m"] - self.p.radius_m
        self.wall_stats.assign(d)

    def wall_pass(self, x):
        super().launch(
            vk.wall_project,
            self.total,
            [
                *self.geometry.args(),
                x,
                self.p.radius_m,
                self.active,
                self.cache_x,
                self.cache_d,
                self.wall_stats,
                self.error,
            ],
        )

    def step(self, dt):
        self._pair_contact_pass = 0
        return super().step(dt)

    def launch(self, kernel, dim, args):
        if not self.vessel_enabled or not self._vessel_ready:
            return super().launch(kernel, dim, args)
        if kernel is k.accept:
            g = self.geometry
            x = args[0]
            super().launch(
                vk.clip_boundaries,
                self.total,
                [
                    g.closed.id,
                    g.tags,
                    g.origin,
                    g.padding,
                    self.x,
                    x,
                    self.active,
                    self.error,
                    g.cap_origins,
                    g.cap_normals,
                    self.field_lower,
                    self.field_upper,
                    self.field_enabled,
                    self.pending,
                    self.fractions,
                ],
            )
            super().launch(
                vk.final_geometry_check,
                self.total,
                [
                    *g.args(),
                    self.x,
                    x,
                    self.n,
                    self.p.radius_m,
                    self.active,
                    self.cache_x,
                    self.cache_d,
                    self.wall_stats,
                    self.error,
                    g.cap_origins,
                    g.cap_normals,
                    self.field_lower,
                    self.field_upper,
                    self.field_enabled,
                    int(self.wall_contact_enabled),
                    self.pending,
                    self.fractions,
                    self.wall_gap,
                ],
            )
            super().launch(
                vk.commit_events,
                self.total,
                [
                    self.active,
                    self.pending,
                    self.fractions,
                    self.exit_label,
                    self.exit_time,
                    self.wall_stats,
                    self.error,
                    self.step_count,
                    self.time_origin,
                    self._current_dt,
                ],
            )
        mapping = {
            k.field_response: masked.field_response,
            k.local_mobility: masked.local_mobility,
            k.component_links: masked.component_links,
            k.contact: masked.contact,
            k.accept: masked.accept,
            f3.pair_buffer: masked.pair_buffer,
            f2.pair_force: masked.pair_force,
        }
        if kernel is k.predict and self.unit_velocity is not None:
            super().launch(
                vk.prescribed_prediction,
                self.total,
                [self.x, self.unit_velocity, args[2], self._current_dt, self.active],
            )
        elif kernel in mapping:
            super().launch(mapping[kernel], dim, [*args, self.active])
        else:
            super().launch(kernel, dim, args)
        if self.wall_contact_enabled:
            if kernel is k.predict:
                self.wall_pass(args[2])
            if kernel is k.contact:
                self._pair_contact_pass += 1
                if self._pair_contact_pass % 3 == 0:
                    self.wall_pass(args[1])

    def state(self):
        x = super().state()
        if self._vessel_ready and self._wall_frames is not None:
            super().launch(
                vk.query_points,
                self.total,
                [
                    *self.geometry.args(),
                    self.x,
                    self.query_d,
                    self.query_n,
                    self.query_face,
                    self.query_inside,
                ],
            )
            self._wall_frames.append((self.query_d.numpy() - self.p.radius_m).reshape(self.batch, self.n))
            self._active_frames.append(self.active.numpy().reshape(self.batch, self.n))
            self._wall_stats_frames.append(self.wall_stats.numpy().reshape(self.batch, self.n, 12))
        return x

    def quality(self, sample=None):
        q = super().quality(sample)
        if self._vessel_ready:
            d = self.wall_stats.numpy()
            if sample is not None:
                d = d[sample * self.n : (sample + 1) * self.n]
            q["wall"] = dict(
                minimum_queried_gap_m=float(d[:, 0].min()),
                max_predicted_penetration_m=float(d[:, 1].max()),
                max_pass_correction_m=float(d[:, 2].max()),
                cumulative_correction_m=float(d[:, 3].sum()),
                refined_queries=int(d[:, 4].sum()),
                point_triangle_tests=int(d[:, 5].sum()),
                swept_triangle_tests=int(d[:, 6].sum()),
                correction_passes=int(d[:, 7].sum()),
                particle_contact_steps=int(d[:, 8].sum()),
                lipschitz_fast_paths=int(d[:, 9].sum()),
                terminal_exits=int(d[:, 10].sum()),
                field_domain_exits=int(d[:, 11].sum()),
            )
        return q

    def run(self, *args, **kwargs):
        if not self.vessel_enabled:
            return super().run(*args, **kwargs)
        self._wall_frames = []
        self._active_frames = []
        self._wall_stats_frames = []
        result = super().run(*args, **kwargs)
        assert len(result["time_s"]) == len(self._wall_frames)
        result["wall_gap_m"] = np.asarray(self._wall_frames)
        result["active"] = np.asarray(self._active_frames)
        result["wall_stats"] = np.asarray(self._wall_stats_frames)
        labels = self.exit_label.numpy()
        times = self.exit_time.numpy()
        result["boundary_events"] = [
            dict(
                particle_id=int(i % self.n),
                batch=int(i // self.n),
                label=int(labels[i]),
                event="FIELD_DOMAIN_EXIT" if labels[i] == 100 else "TERMINAL_EXIT",
                time_s=float(times[i]),
            )
            for i in np.flatnonzero(labels)
        ]
        s = result["summary"]
        s.update(
            vessel_enabled=True,
            wall_contact_enabled=self.wall_contact_enabled,
            wall_physics="GEOMETRIC_WALL_ONLY",
            geometry=self.geometry.metadata,
            actual_additional_substeps=0,
            wall_passes_per_step=5 if self.wall_contact_enabled else 0,
            scenario_scope="GEOMETRY_UNIT_TEST" if self.unit_velocity is not None else "ACOUSTIC_SIMULATION",
            terminal_exit_count=int(np.sum((labels > 0) & (labels < 100))),
            field_exit_count=int(np.sum(labels == 100)),
        )
        s["model"]["wall"] = (
            "geometric nonpenetration only; open terminals"
            if self.wall_contact_enabled
            else "disabled matched vessel comparator"
        )
        s["model"]["flow"] = "zero"
        self._wall_frames = None
        self._active_frames = None
        self._wall_stats_frames = None
        return result
