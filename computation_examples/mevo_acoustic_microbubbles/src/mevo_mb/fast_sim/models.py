"""Bounded reuse of the actual acoustic and APECSS databases; NumPy references."""

from dataclasses import dataclass, asdict
from pathlib import Path
import hashlib
import numpy as np

from mevo_mb.acoustics.basis import AcousticBasisDatabase
from mevo_mb.bubbles.adaptive_table import AdaptiveResponseTableBubbleModel, CORNERS


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@dataclass(frozen=True)
class Physics:
    radius_m: float = 1.78e-6
    frequency_Hz: float = 2.5e6
    density_kg_m3: float = 1000.0
    viscosity_Pa_s: float = 0.001
    radial_effective_viscosity_Pa_s: float = 0.002
    secondary: bool = True
    near_field: str = "smooth_near_field_exploratory"
    transition_over_sum: float = 2.0
    soft_core_over_sum: float = 1.0
    local_lubrication: bool = True
    lubrication_activation: float = 0.1
    lubrication_floor: float = 0.001
    mobility_iterations: int = 2
    contact_iterations: int = 4
    contact_relaxation: float = 0.5

    def __post_init__(self):
        if self.near_field not in ("classical_far_field", "smooth_near_field_exploratory"):
            raise ValueError("Unknown near-field model")
        if not (self.radius_m > 0 and self.frequency_Hz > 0 and self.viscosity_Pa_s > 0):
            raise ValueError("Positive physical parameters required")
        if not 0 < self.lubrication_floor < self.lubrication_activation:
            raise ValueError("Invalid lubrication range")
        if not 1 <= self.mobility_iterations <= 4 or not 1 <= self.contact_iterations <= 12:
            raise ValueError("Bounded local iteration counts required")
        if not 0 < self.contact_relaxation <= 1 or self.transition_over_sum <= 1:
            raise ValueError("Invalid contact relaxation / near-field transition")

    def metadata(self):
        return dict(
            asdict(self),
            physics_scope="EXPLORATORY",
            harmonic_convention="Re[P exp(+i omega t)]",
            force_direction="e_ij=(x_j-x_i)/d; positive correlation attracts",
            near_field_formula="F=C*dx/[d^2+(s*Rsum*(1-q)^2)^2]^(3/2), q=clip(d/(transition*Rsum),0,1); classical above transition",
            near_field_calibration="none; signed smooth exploratory regularization, not a Doinikov solution",
            mobility="Two relaxed block-Jacobi local normal-resistance passes, initialized at free Stokes velocity; no global solve",
            contact="Jacobi PBD in diagonal Stokes mobility metric, weight mu_i/(mu_i+mu_j); monodisperse F1 weights 1/2",
            contact_vs_hydrodynamics="PBD uses baseline mobility, not the full coupled hydrodynamic metric",
            wall="none",
            flow="zero",
            cluster_locking=False,
            adhesion=False,
            coalescence=False,
        )


class ResponseSection:
    """Exact restriction of valid adaptive trilinear cells to one fixed R and f.

    Intervals preserve the original max-leaf-index tie rule. Values within an
    interval are linear in pressure; invalid gaps stay invalid. Zero drive is
    the exact quiescent solution, not extrapolation below the table minimum.
    """

    def __init__(self, path, physics):
        self.model = m = AdaptiveResponseTableBubbleModel(path)
        m.validate_medium(asdict(physics))
        r, f = physics.radius_m, physics.frequency_Hz
        contains = (m.lower[:, 0] <= r) & (r <= m.upper[:, 0]) & (m.lower[:, 1] <= f) & (f <= m.upper[:, 1])
        boundaries = np.unique(np.r_[m.lower[contains, 2], m.upper[contains, 2]])
        if len(boundaries) < 2:
            raise ValueError("Radius/frequency outside actual response table")
        lower, upper, values, valid, leaf = [], [], [], [], []
        for lo, hi in zip(boundaries[:-1], boundaries[1:]):
            p = (lo + hi) / 2
            chosen = np.flatnonzero(contains & m.valid & (m.lower[:, 2] <= p) & (p <= m.upper[:, 2]))
            cell = int(chosen[-1]) if len(chosen) else -1
            vals = np.zeros((2, 8), complex)
            if cell >= 0:
                for end, pressure in enumerate([lo, hi]):
                    w = (np.array([r, f, pressure]) - m.lower[cell]) / (m.upper[cell] - m.lower[cell])
                    coeff = np.prod(np.where(CORNERS, w, 1 - w), axis=1)
                    vals[end] = coeff @ m.harmonics[cell]
            lower.append(lo)
            upper.append(hi)
            values.append(vals)
            valid.append(cell >= 0)
            leaf.append(cell)
        self.lower, self.upper = np.array(lower), np.array(upper)
        self.values, self.valid = np.array(values), np.array(valid, np.int32)
        self.knots = boundaries
        self.knot_values = np.zeros((len(boundaries), 8), complex)
        self.knot_valid = np.zeros(len(boundaries), np.int32)
        for k, pressure in enumerate(boundaries):
            try:
                self.knot_values[k] = m.query(r, f, pressure)["harmonics_m3"]
                self.knot_valid[k] = 1
            except ValueError:
                pass
        self.metadata = dict(
            path=str(Path(path).resolve()),
            sha256=sha(path),
            original_metadata=m.metadata,
            fixed_radius_m=r,
            fixed_frequency_Hz=f,
            intervals_Pa=np.c_[lower, upper].tolist(),
            valid=self.valid.tolist(),
            original_leaf_indices=leaf,
            zero_pressure="exact zero acoustic response",
            extrapolation=False,
        )

    def query(self, p):
        p = np.asarray(p)
        out = np.zeros(p.shape + (8,), complex)
        nonzero = np.abs(p) > 0
        if np.any(nonzero):
            out[nonzero] = self.model.harmonic_response(
                self.metadata["fixed_radius_m"], self.metadata["fixed_frequency_Hz"], p[nonzero]
            )["harmonics_m3"]
        return out


class Basis:
    def __init__(self, path, frequency_Hz):
        self.db = AcousticBasisDatabase.load_basis(path)
        if self.db.metadata["frequency_Hz"] != frequency_Hz:
            raise ValueError("Acoustic basis is frequency-specific")
        self.axes = self.db.axes_m
        self.origin = np.array([a[0] for a in self.axes])
        self.spacing = np.array([a[1] - a[0] for a in self.axes])
        if any(
            not np.allclose(np.diff(a), d, rtol=1e-10, atol=1e-15) for a, d in zip(self.axes, self.spacing)
        ):
            raise ValueError("Warp F1 requires a uniform stored grid")
        self.shape = tuple(map(len, self.axes))
        self.metadata = dict(
            path=str(Path(path).resolve()),
            sha256=sha(path),
            basis_metadata=self.db.metadata,
            grid_shape=self.shape,
            axes_bounds_m=[[float(a[0]), float(a[-1])] for a in self.axes],
        )

    def combine(self, amplitude_Pa, phase_rad, enabled=None):
        a, phase = np.asarray(amplitude_Pa), np.asarray(phase_rad)
        if a.shape != (len(self.db.pressure_basis),) or phase.shape != a.shape or np.any(a < 0):
            raise ValueError("One nonnegative amplitude and phase per actual PZT")
        enabled = np.ones(a.shape, bool) if enabled is None else np.asarray(enabled, bool)
        if enabled.shape != a.shape or not np.isfinite(a).all() or not np.isfinite(phase).all():
            raise ValueError("Invalid source controls")
        u = a * enabled * np.exp(1j * phase)
        field = self.db.combine(u)
        packed = np.concatenate([field.pressure_Pa[..., None], field.gradient_Pa_m], axis=-1)
        return field, np.ascontiguousarray(np.stack([packed.real, packed.imag], axis=-1)), u


def pair_forces_numpy(x, harmonics, p):
    delta = x[None, :, :] - x[:, None, :]
    d2 = np.sum(delta**2, axis=-1)
    if np.any((d2 == 0) & ~np.eye(len(x), dtype=bool)):
        raise ValueError("Coincident centers")
    d = np.sqrt(d2)
    if p.near_field == "smooth_near_field_exploratory":
        q = np.clip(d / (p.transition_over_sum * 2 * p.radius_m), 0, 1)
        core = p.soft_core_over_sum * 2 * p.radius_m * (1 - q) ** 2
        d2 = d2 + core**2
    np.fill_diagonal(d2, np.inf)
    dotv = 2j * np.pi * p.frequency_Hz * harmonics * np.arange(1, 9)
    corr = 0.5 * np.real(dotv @ dotv.conj().T)
    return p.density_kg_m3 / (4 * np.pi) * np.sum(delta * (corr / d2**1.5)[..., None], axis=1)


def random_cloud(n, radius, cloud_radius, seed, center=(0, 0, 0)):
    rng = np.random.default_rng(seed)
    points = []
    for _ in range(100000):
        x = rng.uniform(-cloud_radius, cloud_radius, 3)
        if np.linalg.norm(x) > cloud_radius:
            continue
        if points and np.min(np.linalg.norm(np.asarray(points) - x, axis=1)) < 2.2 * radius:
            continue
        points.append(x)
        if len(points) == n:
            return np.asarray(points) + center
    raise ValueError("Requested concentration cannot be initialized without overlap")
