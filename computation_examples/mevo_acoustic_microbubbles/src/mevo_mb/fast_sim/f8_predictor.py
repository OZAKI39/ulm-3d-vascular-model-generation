"""GPU force/RPY candidate evaluator, separate from all historical integrators."""

import numpy as np
import warp as wp
from .consistent_field import query_numpy
from .engine import complex_vectors
from . import f7_kernels as fk, f3_kernels as f3, rpy
from .f5_engine import flow_free_velocity


class MatchedPredictor:
    def __init__(self, basis, response, physics, n=32, device="cuda:0", mode="rpy_free_space"):
        if mode not in ("rpy_free_space", "diagonal_stokes"):
            raise ValueError("Explicit candidate mobility required")
        self.mode = mode
        self.basis, self.response, self.p, self.n = basis, response, physics, n
        self.device = wp.get_device(device)

        def arr(a, dtype):
            return wp.array(np.ascontiguousarray(a), dtype=dtype, device=self.device)

        def zeros(shape, dtype):
            return wp.zeros(shape, dtype=dtype, device=self.device)

        self.x = zeros(n, wp.vec3d)
        self.uf = zeros(n, wp.vec3d)
        self.samples = zeros((n, 4), wp.vec2d)
        self.gradients = zeros((n, 3, 4), wp.vec2d)
        self.controls = zeros(4, wp.vec2d)
        self.force = zeros(n, wp.vec3d)
        self.primary = zeros(n, wp.vec3d)
        self.velocity = zeros(n, wp.vec3d)
        self.pressure = zeros(n, wp.vec2d)
        self.harmonics = zeros((n, 8), wp.vec2d)
        self.pairs = zeros((n, n), wp.vec3d)
        self.diag = zeros((n, 11), wp.float64)
        self.near = zeros((n, 3), wp.float64)
        self.active = arr(np.ones(n, np.int32), wp.int32)
        self.path = zeros(n, wp.float64)
        self.error = zeros(1, wp.int32)
        self.knots = arr(response.knots, wp.float64)
        self.values = arr(complex_vectors(response.values), wp.vec2d)
        self.valid = arr(response.valid, wp.int32)
        self.kv = arr(complex_vectors(response.knot_values), wp.vec2d)
        self.kval = arr(response.knot_valid, wp.int32)
        self.chi = response.query(np.array([response.knots[0] + 0j]))[0, 0] / response.knots[0]
        for mod in (fk, f3, rpy):
            wp.load_module(module=mod, device=self.device)
        with wp.ScopedCapture(device=self.device) as cap:
            self.launch()
        self.graph = cap.graph

    def launch(self):
        p = self.p
        n = self.n

        def run(kernel, dim, args):
            wp.launch(kernel, dim=dim, inputs=args, device=self.device)

        self.error.zero_()
        run(
            fk.cached_response,
            n,
            [
                self.samples,
                self.gradients,
                self.controls,
                self.knots,
                self.values,
                self.valid,
                self.kv,
                self.kval,
                self.harmonics,
                self.force,
                self.pressure,
                self.diag,
                self.error,
            ],
        )
        wp.copy(self.primary, self.force)
        run(
            f3.pair_buffer,
            (1, n, n),
            [
                self.x,
                self.harmonics,
                self.pressure,
                self.pairs,
                n,
                p.radius_m,
                2 * np.pi * p.frequency_Hz,
                p.density_kg_m3,
                p.transition_over_sum,
                p.soft_core_over_sum,
                p.feedback_end_over_sum,
                wp.vec2d(self.chi.real, self.chi.imag),
                float(self.response.knots[0]),
                float(self.response.knots[-1]),
                self.near,
                self.error,
            ],
        )
        run(f3.pair_sum, n, [self.pairs, self.force, n])
        if self.mode == "diagonal_stokes":
            run(
                flow_free_velocity,
                n,
                [
                    self.force,
                    self.uf,
                    self.velocity,
                    self.active,
                    6 * np.pi * p.viscosity_Pa_s * p.radius_m,
                    1,
                    0.0,
                    self.path,
                ],
            )
        else:
            run(
                rpy.velocity,
                n,
                [
                    self.x,
                    self.force,
                    self.uf,
                    self.active,
                    self.velocity,
                    n,
                    p.radius_m,
                    p.viscosity_Pa_s,
                    1,
                    0.0,
                    self.path,
                    self.error,
                ],
            )

    def prepare(self, positions, flow):
        x = np.asarray(positions, float)
        uf = np.asarray(flow, float)
        if x.shape != (self.n, 3) or uf.shape != x.shape or not np.isfinite(uf).all():
            raise ValueError("Finite actual per-bubble positions and flow required")
        samples = [
            query_numpy(v, self.basis.origin, self.basis.spacing, x) for v in self.basis.db.pressure_basis
        ]
        self.samples.assign(complex_vectors(np.stack([s[0] for s in samples], axis=-1)))
        self.gradients.assign(complex_vectors(np.stack([s[1] for s in samples], axis=-1)))
        self.x.assign(x)
        self.uf.assign(uf)
        self.positions = x.copy()

    def evaluate(self, amplitudes, phases, enabled=None, details=False):
        a = np.asarray(amplitudes)
        ph = np.asarray(phases)
        if a.shape != (4,) or ph.shape != (4,) or np.any(a < 0) or not np.isfinite(a + ph).all():
            raise ValueError("Four finite nonnegative controls required")
        u = a * np.exp(1j * ph)
        if enabled is not None:
            u = u * np.asarray(enabled)
        self.controls.assign(complex_vectors(u))
        wp.capture_launch(self.graph)
        velocity = self.velocity.numpy()
        error = int(self.error.numpy()[0])
        if error:
            raise ValueError("Candidate outside pressure/closure domain: " + str(error))
        if not np.isfinite(velocity).all():
            raise ValueError("Nonfinite candidate velocity")
        if details:
            return dict(
                velocity_m_s=velocity,
                primary_N=self.primary.numpy(),
                secondary_N=self.force.numpy() - self.primary.numpy(),
                pressure_complex_Pa=self.pressure.numpy(),
            )
        return velocity
