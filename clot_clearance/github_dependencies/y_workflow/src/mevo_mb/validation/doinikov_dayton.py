"""N=1 reference only: Doinikov & Dayton 2006, DOI 10.1121/1.2215228.

Equations 6,31,32,33,36,40; thick viscoelastic solid shell. This is NOT the
SonoVue backend. No public author implementation was found (search evidence
under results/stage2/translation). SciPy integrates these four particular ODEs.
No history force, shell buckling, dissolution, streaming or finite-Re drag.
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class DoinikovDaytonParameters:
    outer_radius_m: float = 2e-6
    shell_thickness_m: float = 15e-9
    frequency_Hz: float = 2.5e6
    pressure_Pa: float = 150e3
    density_liquid_kg_m3: float = 1000.0
    density_shell_kg_m3: float = 1100.0
    viscosity_liquid_Pa_s: float = 0.001
    viscosity_shell_Pa_s: float = 1.77
    shear_modulus_Pa: float = 88.8e6
    sigma_inner_N_m: float = 0.04
    sigma_outer_N_m: float = 0.005
    ambient_pressure_Pa: float = 101300.0
    sound_speed_m_s: float = 1500.0
    gas_exponent: float = 1.4
    gas_density_kg_m3: float = 1.2

    def __post_init__(self):
        if not 0 <= self.shell_thickness_m < self.outer_radius_m:
            raise ValueError("Shell thickness must be smaller than outer radius")
        if self.shear_modulus_Pa <= 0 or self.frequency_Hz <= 0:
            raise ValueError("Positive modulus and frequency required")


class DoinikovDaytonSingleBubbleReference:
    def __init__(self, parameters: DoinikovDaytonParameters):
        self.p = p = parameters
        self.r0 = p.outer_radius_m - p.shell_thickness_m
        self.delta = p.outer_radius_m**3 - self.r0**3
        self.q = (p.density_liquid_kg_m3 - p.density_shell_kg_m3) / p.density_shell_kg_m3
        self.omega = 2 * np.pi * p.frequency_Hz
        self.k = self.omega / p.sound_speed_m_s
        self.mb = 4 * np.pi / 3 * (p.density_shell_kg_m3 * self.delta + p.gas_density_kg_m3 * self.r0**3)
        self.re = self.r0 * (
            1
            - (p.ambient_pressure_Pa + 2 * p.sigma_outer_N_m / p.outer_radius_m)
            * p.outer_radius_m**3
            / (4 * p.shear_modulus_Pa * self.r0**3)
        )  # Eq40
        d = 1 - self.re / self.r0
        self.pg0 = (
            p.ambient_pressure_Pa
            + 2 * p.sigma_inner_N_m / self.r0
            + 2 * p.sigma_outer_N_m / p.outer_radius_m
            + 4
            * p.shear_modulus_Pa
            * self.delta
            / p.outer_radius_m**3
            * d
            * (1 + 0.5 * d * (1 - 3 * self.r0**3 / p.outer_radius_m**3))
        )  # Eq33

    def outer(self, r):
        return (r**3 + self.delta) ** (1 / 3)

    def pressure(self, x, t):
        return self.p.pressure_Pa * np.sin(self.omega * t - self.k * x)

    def gradient(self, x, t):
        return -self.k * self.p.pressure_Pa * np.cos(self.omega * t - self.k * x)

    def G(self, r, u, x, v, t):
        """Right hand side of published Eq31, in m²/s²."""
        p = self.p
        s = self.outer(r)
        strain = 1 - self.re / r
        return (
            p.density_liquid_kg_m3 / p.density_shell_kg_m3 * v * v / 4
            + self.pg0 / p.density_shell_kg_m3 * (self.r0 / r) ** (3 * p.gas_exponent)
            - 2 * p.sigma_inner_N_m / (p.density_shell_kg_m3 * r)
            - 2 * p.sigma_outer_N_m / (p.density_shell_kg_m3 * s)
            - (p.ambient_pressure_Pa + self.pressure(x, t)) / p.density_shell_kg_m3
            - 4
            * u
            / (p.density_shell_kg_m3 * r * s**3)
            * (p.viscosity_liquid_Pa_s * r**3 + p.viscosity_shell_Pa_s * self.delta)
            - 4
            * p.shear_modulus_Pa
            * self.delta
            / (p.density_shell_kg_m3 * s**3)
            * strain
            * (1 + 0.5 * strain * (1 - 3 * r**3 / s**3))
        )

    def derivatives_G(self, r, u, x, v, t):
        # Complex-step differentiation of a smooth algebraic expression avoids
        # subtractive cancellation. No finite-difference ODE or hidden fitting.
        args = np.array([r, u, x, v, t], dtype=complex)
        partials = []
        for i in range(5):
            perturbed = args.copy()
            perturbed[i] += 1e-25j
            partials.append(float(np.imag(self.G(*perturbed)) / 1e-25))
        return partials

    def rhs(self, t, y):
        r, u, x, v = y
        p = self.p
        if r <= 0:
            raise RuntimeError("Nonpositive inner radius in reference solver")
        s = self.outer(r)
        alpha = 1 + self.q * r / s
        beta = 1.5 + self.q * r / s * (4 * s**3 - r**3) / (2 * s**3)
        A = 1 + self.q * r**4 / s**4
        C = 1 + self.q * r**4 * (2 * s**3 - r**3) / s**7
        mass = self.mb + 2 * np.pi / 3 * p.density_liquid_kg_m3 * s**3
        vdot = (
            -4 * np.pi / 3 * s**3 * self.gradient(x, t)
            - 6 * np.pi * p.viscosity_liquid_Pa_s * s * v
            - 2 * np.pi * p.density_liquid_kg_m3 * r * r * u * v
        ) / mass  # Eq32, Eq6
        gr, gu, gx, gv, gt = self.derivatives_G(r, u, x, v, t)
        compress = p.density_liquid_kg_m3 / (p.sound_speed_m_s * p.density_shell_kg_m3 * alpha)
        # Substitute Eq36 into Eq31, collect all Rddot terms. Gdot contains
        # gu*Rddot; gv*vdot is known from Eq32. No third derivative is integrated.
        numerator = (
            self.G(r, u, x, v, t)
            - beta * u * u
            + compress * (r * (gr * u + gx * v + gv * vdot + gt) + 2 * u**3 * C)
        )
        denominator = alpha * r - compress * (r * gu + 2 * r * u * A)
        if denominator <= 0:
            raise RuntimeError("Compressible reference denominator outside validity")
        return np.array([u, numerator / denominator, v, vdot])

    def solve(self, cycles=20, samples_per_cycle=256, rtol=1e-9, max_steps_per_cycle=128):
        scale = np.array([self.r0, self.r0 * self.omega, self.r0, self.r0 * self.omega])

        def fun(tau, y):
            return self.rhs(tau / self.omega, y * scale) / (scale * self.omega)

        tau = np.linspace(0, 2 * np.pi * cycles, cycles * samples_per_cycle + 1)
        result = solve_ivp(
            fun,
            (0.0, tau[-1]),
            [1.0, 0.0, 0.0, 0.0],
            t_eval=tau,
            method="DOP853",
            rtol=rtol,
            atol=rtol * 0.01,
            max_step=2 * np.pi / max_steps_per_cycle,
        )
        if not result.success or len(result.t) != len(tau):
            raise RuntimeError(result.message)
        y = result.y.T * scale
        return dict(
            time_s=result.t / self.omega,
            inner_radius_m=y[:, 0],
            outer_radius_m=self.outer(y[:, 0]),
            inner_speed_m_s=y[:, 1],
            displacement_m=y[:, 2],
            velocity_m_s=y[:, 3],
            nfev=result.nfev,
        )

    def linear_force(self):
        """Published Eq63 with Eqs44–49; independent analytic small-signal check."""
        p, r, s = self.p, self.r0, self.p.outer_radius_m
        alpha = 1 + self.q * r / s
        delta_eta = (
            4
            * (p.viscosity_liquid_Pa_s * r**3 + p.viscosity_shell_Pa_s * self.delta)
            / (alpha * p.density_shell_kg_m3 * r * r * s**3)
        )
        delta_rad = (
            p.density_liquid_kg_m3 * r * self.omega**2 / (alpha * p.sound_speed_m_s * p.density_shell_kg_m3)
        )
        D = 1 - self.re / r
        omega0sq = (
            3 * p.gas_exponent * self.pg0
            - 2 * p.sigma_inner_N_m / r
            - 2 * p.sigma_outer_N_m * r**3 / s**4
            + 4 * p.shear_modulus_Pa * self.delta / s**3 * (1 - D * (1 + 3 * r**3 / s**3))
        ) / (alpha * p.density_shell_kg_m3 * r * r)
        damping = delta_eta + delta_rad
        return (
            2
            * np.pi
            * p.pressure_Pa**2
            * r
            * self.omega
            * damping
            * self.k
            / (alpha * p.density_shell_kg_m3 * ((self.omega**2 - omega0sq) ** 2 + self.omega**2 * damping**2))
        )
