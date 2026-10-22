"""Uncalibrated lossless uniform baffled-piston control basis in water.

FNM single integral, Kelly & McGough 2006 Eq. (9), PMC2117898.
H = p/(rho*c*v_n), Re[P exp(+i*omega*t)]; source_Pa means rho*c*v_n,
not a measured face pressure or voltage. Each piston radiates independently;
no PDMS, mutual scattering, elastic waves or transducer transfer function.
"""

import json
import numpy as np
import h5py
from numpy.polynomial.legendre import leggauss
from .basis import ComplexAcousticField


class WaterDiscs:
    def __init__(self, elements, frequency_Hz=2.5e6, sound_speed_m_s=1500.0):
        self.elements = elements
        self.frequency = float(frequency_Hz)
        self.k = 2 * np.pi * self.frequency / sound_speed_m_s

    def evaluate(self, points, order=48):
        """Pressure and analytic derivative of the SAME converged integral."""
        x = np.atleast_2d(points)
        nodes, weights = leggauss(order)
        cs = np.cos((nodes + 1) * np.pi / 2)
        weights = weights / 2
        pressures, gradients = [], []
        for e in self.elements:
            n = np.asarray(e["normal_local"], float)
            delta = x - e["center_local_m"]
            z = delta @ n
            lateral = delta - z[:, None] * n
            r = np.linalg.norm(lateral, axis=1)
            a = e["active_diameter_m"] / 2
            if np.any(z <= 0) or not np.isfinite(x).all():
                raise ValueError("Only finite points in front of each piston are admitted")
            rr, zz = r[:, None], z[:, None]
            D = rr**2 + a * a - 2 * a * rr * cs
            N = a * (a - rr * cs)
            W = N / D
            d = np.sqrt(zz**2 + D)
            E = np.exp(-1j * self.k * d)
            E0 = np.exp(-1j * self.k * zz)
            Wr = (-a * cs * D - N * (2 * rr - 2 * a * cs)) / D**2
            p = ((W * (E0 - E)) * weights).sum(1)
            pr = ((Wr * (E0 - E) + W * 1j * self.k * E * (rr - a * cs) / d) * weights).sum(1)
            pz = ((W * (-1j * self.k * E0 + 1j * self.k * E * zz / d)) * weights).sum(1)
            direction = np.divide(lateral, r[:, None], out=np.zeros_like(lateral), where=r[:, None] > 0)
            pressures.append(p)
            gradients.append(pz[:, None] * n + pr[:, None] * direction)
        return np.asarray(pressures), np.asarray(gradients)

    def rayleigh_reference(self, points, radial_order=64, angular_order=192):
        """Independent two-dimensional area quadrature including its derivative."""
        x = np.atleast_2d(points)
        t, w = leggauss(radial_order)
        angle = np.arange(angular_order) * 2 * np.pi / angular_order
        pp, gg = [], []
        for e in self.elements:
            n = np.asarray(e["normal_local"], float)
            u = np.cross(n, [0, 0, 1])
            u /= np.linalg.norm(u)
            v = np.cross(n, u)
            a = e["active_diameter_m"] / 2
            r = (t + 1) * a / 2
            src = np.asarray(e["center_local_m"]) + (
                r[:, None, None] * (np.cos(angle)[None, :, None] * u + np.sin(angle)[None, :, None] * v)
            ).reshape(-1, 3)
            weights = np.repeat(w * a / 2 * r * 2 * np.pi / angular_order, angular_order)
            dvec = x[:, None, :] - src
            d = np.linalg.norm(dvec, axis=-1)
            wave = np.exp(-1j * self.k * d)
            pp.append(1j * self.k / (2 * np.pi) * np.sum(wave / d * weights, axis=1))
            gg.append(
                1j
                * self.k
                / (2 * np.pi)
                * np.sum((wave * (-1j * self.k * d - 1) / d**3 * weights)[..., None] * dvec, axis=1)
            )
        return np.asarray(pp), np.asarray(gg)


def cubic_weights(t):
    w = np.stack(
        [
            -0.5 * t + t * t - 0.5 * t**3,
            1 - 2.5 * t * t + 1.5 * t**3,
            0.5 * t + 2 * t * t - 1.5 * t**3,
            -0.5 * t * t + 0.5 * t**3,
        ],
        axis=-1,
    )
    d = np.stack(
        [-0.5 + 2 * t - 1.5 * t * t, -5 * t + 4.5 * t * t, 0.5 + 4 * t - 4.5 * t * t, -t + 1.5 * t * t],
        axis=-1,
    )
    return w, d


def cubic_query(pressure, origin, spacing, points):
    """C1 Catmull-Rom interpolant and exact polynomial gradient, no extrapolation."""
    x = np.atleast_2d(points)
    q = (x - origin) / spacing
    shape = np.array(pressure.shape)
    if not np.isfinite(q).all() or np.any(q < 1) or np.any(q >= shape - 2):
        raise ValueError("Outside cubic cache interior; no extrapolation")
    cell = np.floor(q).astype(int)
    w, d = cubic_weights(q - cell)
    p = np.zeros(len(x), complex)
    g = np.zeros((len(x), 3), complex)
    for a in range(4):
        for b in range(4):
            for c in range(4):
                v = pressure[tuple((cell + np.array([a, b, c]) - 1).T)]
                p += w[:, 0, a] * w[:, 1, b] * w[:, 2, c] * v
                g[:, 0] += d[:, 0, a] * w[:, 1, b] * w[:, 2, c] * v / spacing[0]
                g[:, 1] += w[:, 0, a] * d[:, 1, b] * w[:, 2, c] * v / spacing[1]
                g[:, 2] += w[:, 0, a] * w[:, 1, b] * d[:, 2, c] * v / spacing[2]
    return p, g


class CachedWaterBasis:
    def __init__(self, path):
        with h5py.File(path) as h:
            self.axes = tuple(h[f"axes/{s}_m"][:] for s in "xyz")
            self.pressures = h["pressure_basis"][:]
            self.metadata = json.loads(h.attrs["metadata_json"])
        self.origin = np.array([a[0] for a in self.axes])
        self.spacing = np.array([a[1] - a[0] for a in self.axes])
        self.shape = tuple(len(a) for a in self.axes)
        self.combine_count = 0

    def vector(self, amplitude_Pa, phase_rad, enabled=None):
        a, p = np.asarray(amplitude_Pa, float), np.asarray(phase_rad, float)
        enabled = np.ones(4, bool) if enabled is None else np.asarray(enabled, bool)
        if (
            a.shape != (4,)
            or p.shape != (4,)
            or enabled.shape != (4,)
            or not np.isfinite(a).all()
            or not np.isfinite(p).all()
            or np.any(a < 0)
        ):
            raise ValueError("Four finite independent amplitudes/phases/enables required")
        return a * enabled * np.exp(1j * p)

    def combine(self, amplitude_Pa, phase_rad, enabled=None):
        u = self.vector(amplitude_Pa, phase_rad, enabled)
        p = np.einsum("i,ixyz->xyz", u, self.pressures)
        # Nodal gradients are informational only. GPU forces use the derivative
        # of the SAME cubic pressure, never interpolation of these gradients.
        g = np.stack(np.gradient(p, *self.spacing, edge_order=2), axis=-1)
        field = ComplexAcousticField(self.axes, p, g, self.metadata["frequency_Hz"])
        packed = np.concatenate([p[..., None], g], axis=-1)
        self.combine_count += 1
        return field, np.ascontiguousarray(np.stack([packed.real, packed.imag], axis=-1)), u

    def query_basis(self, x):
        rows = [cubic_query(p, self.origin, self.spacing, x) for p in self.pressures]
        return np.array([r[0] for r in rows]), np.array([r[1] for r in rows])
