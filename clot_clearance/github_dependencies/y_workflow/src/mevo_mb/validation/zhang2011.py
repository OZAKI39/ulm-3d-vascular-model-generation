"""N=2 reference: Zhang et al. CPB 20,114302 (2011), Eqs 3–11.

DOI 10.1088/1674-1056/20/11/114302. Elastic-shell modified Herring model,
NOT full buckling/rupture Marmottant and NOT a general radial backend.
SI units. Uniform acoustic excitation; no primary pressure gradient.
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Zhang2011Parameters:
    radius1_m: float = 3e-6
    radius2_m: float = 6e-6
    distance_m: float = 50e-6
    frequency_Hz: float = 2.25e6
    pressure_Pa: float = 1e5
    density_kg_m3: float = 998.0
    sound_speed_m_s: float = 1500.0
    ambient_pressure_Pa: float = 1.01e5
    surface_tension_N_m: float = 0.051
    gas_exponent: float = 1.07
    viscosity_Pa_s: float = 0.001
    elasticity_N_m: float = 1.0
    shell_viscosity_Pa_s: float = 1.27
    shell_thickness_m: float = 2e-9

    def __post_init__(self):
        if min(self.radius1_m, self.radius2_m) <= self.shell_thickness_m:
            raise ValueError("Radius must exceed shell thickness")
        if self.distance_m <= self.radius1_m + self.radius2_m or self.frequency_Hz <= 0:
            raise ValueError("Separated bubbles and positive frequency required")


class Zhang2011EncapsulatedPairReference:
    def __init__(self, parameters=Zhang2011Parameters()):
        self.p = parameters

    def wall_pressure(self, r, u, r0):
        """Eq9 and analytic partials dP/dR, dP/dU for Eqs7,8."""
        p = self.p
        gas = (p.ambient_pressure_Pa + 2 * p.surface_tension_N_m / r0) * (r0 / r) ** (3 * p.gas_exponent)
        e, m = p.shell_thickness_m, 12 * p.shell_viscosity_Pa_s * p.shell_thickness_m
        pressure = (
            gas
            - 4 * p.viscosity_Pa_s * u / r
            - 2 * p.surface_tension_N_m / r
            - 4 * p.elasticity_N_m * (1 / r0 - 1 / r)
            - m * u / (r * (r - e))
        )
        pr = (
            -3 * p.gas_exponent * gas / r
            + 4 * p.viscosity_Pa_s * u / r**2
            + 2 * p.surface_tension_N_m / r**2
            - 4 * p.elasticity_N_m / r**2
            + m * u * (2 * r - e) / (r * r * (r - e) ** 2)
        )
        pu = -4 * p.viscosity_Pa_s / r - m / (r * (r - e))
        return pressure, pr, pu

    def rhs(self, t, y):
        # q=[R1,R2,x1,x2], followed by their four velocities.
        r, s, x, z, u, w, v, h = y
        p, d = self.p, z - x
        if min(r, s) <= p.shell_thickness_m or d <= 0:
            raise RuntimeError("Reference left spherical positive-radius domain")
        rho = p.density_kg_m3
        p1, pr1, pu1 = self.wall_pressure(r, u, p.radius1_m)
        p2, pr2, pu2 = self.wall_pressure(s, w, p.radius2_m)
        drive = -p.pressure_Pa * np.sin(2 * np.pi * p.frequency_Hz * t)
        drag1 = -12 * np.pi * p.viscosity_Pa_s * r * (v - (-s * s * w / d**2 + s**3 * h / d**3))
        drag2 = -12 * np.pi * p.viscosity_Pa_s * s * (h - (r * r * u / d**2 + r**3 * v / d**3))
        # Collect all four accelerations on the LHS; Pdot includes P_U*Rddot.
        matrix = np.array(
            [
                [r - r * pu1 / (rho * p.sound_speed_m_s), s * s / d, 0.0, -(s**3) / (2 * d * d)],
                [r * r / d, s - s * pu2 / (rho * p.sound_speed_m_s), r**3 / (2 * d * d), 0.0],
                [0.0, r * s * s / d**2, r / 3, -r * s**3 / d**3],
                [-s * r * r / d**2, 0.0, -s * r**3 / d**3, s / 3],
            ]
        )
        right = np.array(
            [
                (p1 + r / p.sound_speed_m_s * pr1 * u - p.ambient_pressure_Pa - drive) / rho
                - 1.5 * u * u
                + v * v / 4
                - 2 * s * w * w / d
                + s * s * (v * w + 5 * w * h) / (2 * d * d)
                - s**3 * h * (v + 2 * h) / (2 * d**3),
                (p2 + s / p.sound_speed_m_s * pr2 * w - p.ambient_pressure_Pa - drive) / rho
                - 1.5 * w * w
                + h * h / 4
                - 2 * r * u * u / d
                - r * r * (h * u + 5 * u * v) / (2 * d * d)
                - r**3 * v * (h + 2 * v) / (2 * d**3),
                drag1 / (2 * np.pi * rho * r * r)
                - u * v
                - (u * s * s * w + 2 * r * s * w * w) / d**2
                + s * s * (s * u * h + 5 * r * w * h) / d**3,
                drag2 / (2 * np.pi * rho * s * s)
                - w * h
                + (w * r * r * u + 2 * s * r * u * u) / d**2
                + r * r * (r * w * v + 5 * s * u * v) / d**3,
            ]
        )
        return np.r_[y[4:], np.linalg.solve(matrix, right)]

    def solve(self, cycles=64, samples_per_cycle=64, rtol=2e-8, max_steps_per_cycle=48):
        return solve_pair(self, cycles, samples_per_cycle, rtol, max_steps_per_cycle)


def solve_pair(model, cycles, samples_per_cycle, rtol, max_steps_per_cycle):
    """Dimensionless integration of these specific four coupled second-order ODEs."""
    p = model.p
    radius = min(p.radius1_m, p.radius2_m)
    omega = 2 * np.pi * p.frequency_Hz
    scale = np.r_[np.full(4, radius), np.full(4, radius * omega)]
    init = np.array([p.radius1_m, p.radius2_m, 0.0, p.distance_m, 0.0, 0.0, 0.0, 0.0]) / scale

    def fun(tau, y):
        return model.rhs(tau / omega, y * scale) / (scale * omega)

    def contact(tau, y):
        return y[3] - y[2] - y[0] - y[1]

    contact.terminal = True
    contact.direction = -1
    tau = np.linspace(0, 2 * np.pi * cycles, int(round(cycles * samples_per_cycle)) + 1)
    z = solve_ivp(
        fun,
        (0, tau[-1]),
        init,
        method="DOP853",
        t_eval=tau,
        rtol=rtol,
        atol=rtol * 0.01,
        max_step=2 * np.pi / max_steps_per_cycle,
        events=contact,
    )
    if not z.success:
        raise RuntimeError(z.message)
    y = z.y.T * scale
    return dict(
        time_s=z.t / omega,
        radius_m=y[:, :2],
        position_m=y[:, 2:4],
        radial_speed_m_s=y[:, 4:6],
        velocity_m_s=y[:, 6:8],
        distance_m=y[:, 3] - y[:, 2],
        contact_reached=bool(z.t_events[0].size),
        contact_time_s=float(z.t_events[0][0] / omega) if z.t_events[0].size else None,
        nfev=z.nfev,
    )


def secondary_coefficient(z, density, final_cycles=8, frequency=2.25e6):
    """Eq11, cycle average in N*m². Positive denotes attraction."""
    t = z["time_s"]
    use = t >= t[-1] - final_cycles / frequency - 1e-16
    dotv = 4 * np.pi * z["radius_m"] ** 2 * z["radial_speed_m_s"]
    return float(
        density / (4 * np.pi) * np.trapezoid(np.prod(dotv[use], axis=1), t[use]) / (t[use][-1] - t[use][0])
    )
