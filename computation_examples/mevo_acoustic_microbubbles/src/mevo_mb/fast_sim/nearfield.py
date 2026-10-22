"""F2 pairwise monopole-feedback hypothesis and its independent NumPy reference."""

from dataclasses import dataclass, asdict
import numpy as np
from .models import Physics


@dataclass(frozen=True)
class F2Physics(Physics):
    contact_iterations: int = 12
    force_model: str = "f2_near_field_candidate"
    contact_model: str = "f1_pbd_contact"
    feedback_end_over_sum: float = 4.0

    def __post_init__(self):
        super().__post_init__()
        if self.force_model not in ("f1_smooth_near_field", "f2_near_field_candidate"):
            raise ValueError("Unknown force model")
        if self.contact_model not in ("f1_pbd_contact", "f2_mobility_aware_contact"):
            raise ValueError("Unknown contact model")
        if self.near_field != "smooth_near_field_exploratory" or self.feedback_end_over_sum <= 1:
            raise ValueError("F2 comparison requires the frozen F1 geometric regularizer")

    def metadata(self):
        m = super().metadata()
        m.update(asdict(self))
        m.update(
            model_version="fast_pzt_f2_pair_feedback_tensor_contact_v1",
            near_field_candidate="Pairwise linear fundamental monopole feedback using local table secant compliance; cubic signed correlation blend; harmonics 2-8 unchanged",
            near_field_limits="Not a collective radial response, multipole theory or calibrated near-field solution; Primary remains the uncoupled table Primary",
            feedback_window="Cubic smoothstep: 1 at d<=Ri+Rj; 0 at d>=4(Ri+Rj)",
            zero_drive_compliance="Lowest valid table-pressure secant; exact zero total drive still has zero response",
            contact_candidate="Frozen local normal-mobility tensors; symmetric pair H=(Mi+Mj)/2; equal/opposite correction overlap*H*n/(2*n.H*n), with fixed Jacobi relaxation",
            contact_vs_hydrodynamics="Local symmetric relative-mobility metric, not a full coupled Stokes contact projection",
        )
        return m


def feedback_weight(d, radius, end=4.0):
    s = np.clip((np.asarray(d) / (2 * radius) - 1) / (end - 1), 0, 1)
    return 1 - s * s * (3 - 2 * s)


def coupled_fundamental(v_i, v_j, p_i, p_j, d, physics, chi_zero):
    chi_i = v_i / p_i if abs(p_i) > 0 else chi_zero
    chi_j = v_j / p_j if abs(p_j) > 0 else chi_zero
    k = physics.density_kg_m3 * (2 * np.pi * physics.frequency_Hz) ** 2 / (4 * np.pi * d)
    a, b = k * chi_i, k * chi_j
    determinant = 1 - a * b
    if abs(determinant) < 1e-6:
        raise ValueError("Singular pair closure; no force clipping")
    u_i, u_j = (v_i - a * v_j) / determinant, (v_j - b * v_i) / determinant
    return u_i, u_j, p_i - k * u_j, p_j - k * u_i, determinant


def pair_reference(x, harmonics, pressure, physics, chi_zero):
    f = np.zeros_like(x, dtype=float)
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            delta = x[j] - x[i]
            d = float(np.linalg.norm(delta))
            if d <= 0:
                raise ValueError("Coincident bubbles")
            corr = np.sum(np.arange(1, 9) ** 2 * np.real(harmonics[i] * np.conj(harmonics[j])))
            if physics.force_model == "f2_near_field_candidate":
                w = feedback_weight(d, physics.radius_m, physics.feedback_end_over_sum)
                if w > 0:
                    vi, vj, _, _, _ = coupled_fundamental(
                        harmonics[i, 0], harmonics[j, 0], pressure[i], pressure[j], d, physics, chi_zero
                    )
                    corr += w * (
                        np.real(vi * np.conj(vj)) - np.real(harmonics[i, 0] * np.conj(harmonics[j, 0]))
                    )
            q = min(d / (physics.transition_over_sum * 2 * physics.radius_m), 1)
            core = physics.soft_core_over_sum * 2 * physics.radius_m * (1 - q) ** 2
            pair = (
                physics.density_kg_m3
                * (2 * np.pi * physics.frequency_Hz) ** 2
                / (8 * np.pi)
                * (corr * delta)
                / (d * d + core * core) ** 1.5
            )
            f[i] += pair
            f[j] -= pair
    return f
