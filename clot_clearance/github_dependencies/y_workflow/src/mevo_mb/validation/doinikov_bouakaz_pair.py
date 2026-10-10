"""N=2 Doinikov–Bouakaz 2015 potential-flow reference (PRE92,043001).

Finite Legendre projection of Eqs8–13, kinetic energy Eq41, Lagrange Eq39–43,
gas pressure Eq67 and drag Eq68. This evaluates the paper's boundary system
directly, rather than copying its d^-9 symbolic series. Multipole-order
convergence is mandatory. No encapsulation; no N-body superposition.
Author Mathematica supplement exists but its download was inaccessible.
"""

from dataclasses import dataclass
from math import comb

import numpy as np

from .zhang2011 import solve_pair


@dataclass(frozen=True)
class DoinikovBouakazParameters:
    radius1_m: float = 1.5e-6
    radius2_m: float = 2e-6
    distance_m: float = 14e-6
    frequency_Hz: float = 2.95e6
    pressure_Pa: float = 1000.0
    density_kg_m3: float = 998.0
    viscosity_Pa_s: float = 0.001
    surface_tension_N_m: float = 0.072
    ambient_pressure_Pa: float = 101300.0
    sound_speed_m_s: float = 1482.0
    gas_exponent: float = 1.4
    pulse_cycles: int = 60
    # Paper states Gaussian/cycles, but no explicit width convention is given.
    # ASSUMPTION: full duration = 6 sigma. No tuning to reproduce a curve.
    gaussian_duration_sigmas: float = 6.0

    def __post_init__(self):
        if self.distance_m <= self.radius1_m + self.radius2_m:
            raise ValueError("Initially separated spheres required")
        if min(self.radius1_m, self.radius2_m, self.frequency_Hz, self.gaussian_duration_sigmas) <= 0:
            raise ValueError("Positive parameters required")


class DoinikovBouakazPairReferenceModel:
    """Experimental reference, not validated arbitrary-near-contact physics."""

    def __init__(self, parameters=DoinikovBouakazParameters(), multipole_order=6):
        if not 1 <= multipole_order <= 24:
            raise ValueError("Use multipole order 1..24 and check convergence")
        self.p, self.order = parameters, multipole_order
        self.n = np.arange(multipole_order + 1)
        self.binomial = np.array([[comb(int(n + m), int(n)) for m in self.n] for n in self.n])

    def mass_matrix(self, q):
        """Exact symmetric kinetic-energy form for the finite projection, Eq41."""
        r, s, x, z = q
        d = z - x
        n, count = self.n, len(self.n)
        a, b = r / d, s / d
        c12 = self.binomial * (-1.0) ** n[None, :] * b ** (n[None, :] + 1)
        c21 = self.binomial * (-1.0) ** n[:, None] * a ** (n[None, :] + 1)
        k12 = n[:, None] / (n[:, None] + 1) * a ** n[:, None] * c12
        k21 = n[:, None] / (n[:, None] + 1) * b ** n[:, None] * c21
        mat = np.eye(2 * count, dtype=np.result_type(q))
        mat[:count, count:] = -k12
        mat[count:, :count] = -k21
        rhs = np.zeros((2 * count, 4), dtype=mat.dtype)
        rhs[0, 0], rhs[count, 1] = -r, -s
        rhs[1, 2], rhs[count + 1, 3] = -r / 2, -s / 2
        coefficients = np.linalg.solve(mat, rhs)
        c1, c2 = coefficients[:count], coefficients[count:]
        b2, b1 = c12 @ c2, c21 @ c1
        B = np.array(
            [
                r * r * (c1[0] + b2[0]),
                s * s * (c2[0] + b1[0]),
                r * r / 3 * (c1[1] + a * b2[1]),
                s * s / 3 * (c2[1] + b * b1[1]),
            ]
        )
        return -2 * np.pi * self.p.density_kg_m3 * (B + B.T)

    def rhs(self, t, y):
        q, vel = y[:4], y[4:]
        r, s = q[:2]
        p = self.p
        if min(r, s) <= 0 or q[3] <= q[2]:
            raise RuntimeError("Reference left spherical positive-radius domain")
        mass = self.mass_matrix(q)
        dm = np.empty((4, 4, 4))  # dm[k,i,j] = dM_ij/dq_k
        for k in range(4):
            qc = q.astype(complex)
            qc[k] += 1e-25j
            dm[k] = self.mass_matrix(qc).imag / 1e-25
        coriolis = np.einsum("kij,j,k->i", dm, vel, vel) - 0.5 * np.einsum("ijk,j,k->i", dm, vel, vel)
        duration = p.pulse_cycles / p.frequency_Hz
        envelope = np.exp(-0.5 * ((t - duration / 2) / (duration / p.gaussian_duration_sigmas)) ** 2)
        drive = -p.pressure_Pa * envelope * np.sin(2 * np.pi * p.frequency_Hz * t)
        radii = q[:2]
        r0 = np.array([p.radius1_m, p.radius2_m])
        pressure = (
            (p.ambient_pressure_Pa + 2 * p.surface_tension_N_m / r0)
            * (r0 / radii) ** (3 * p.gas_exponent)
            * (1 - 3 * p.gas_exponent * vel[:2] / p.sound_speed_m_s)
            - 4 * p.viscosity_Pa_s * vel[:2] / radii
            - 2 * p.surface_tension_N_m / radii
            - p.ambient_pressure_Pa
            - drive
        )
        force = np.r_[4 * np.pi * radii**2 * pressure, -12 * np.pi * p.viscosity_Pa_s * radii * vel[2:]]
        return np.r_[vel, np.linalg.solve(mass, force - coriolis)]

    def solve(self, cycles=None, samples_per_cycle=48, rtol=2e-8, max_steps_per_cycle=32):
        return solve_pair(self, cycles or self.p.pulse_cycles, samples_per_cycle, rtol, max_steps_per_cycle)
