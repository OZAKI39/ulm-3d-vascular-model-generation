"""SI instantaneous force/velocity diagnostics, independent of trajectory shape inference."""

import numpy as np


def rms_vectors(x):
    return float(np.sqrt(np.mean(np.sum(np.asarray(x) ** 2, axis=-1))))


def affine_velocity(x, u, max_condition=1e8):
    relative = x - x.mean(0)
    velocity = u - u.mean(0)
    condition = float(np.linalg.cond(relative))
    rank = int(np.linalg.matrix_rank(relative))
    if rank < 3 or not np.isfinite(condition) or condition > max_condition:
        return dict(
            status="DEGENERATE",
            condition_number=condition if np.isfinite(condition) else None,
            rank=rank,
            mean_velocity_m_s=u.mean(0),
            relative_velocity_RMS_m_s=rms_vectors(velocity),
        )
    coef, _, _, _ = np.linalg.lstsq(relative, velocity, rcond=None)
    L = coef.T
    E = (L + L.T) / 2
    W = (L - L.T) / 2
    dev = E - np.eye(3) * np.trace(E) / 3
    residual = velocity - relative @ coef
    return dict(
        status="IDENTIFIABLE",
        condition_number=condition,
        rank=rank,
        L_per_s=L,
        E_per_s=E,
        W_per_s=W,
        deviatoric_E_per_s=dev,
        deviatoric_strain_norm_per_s=float(np.linalg.norm(dev)),
        rotation_norm_per_s=float(np.linalg.norm(W)),
        dilation_per_s=float(np.trace(E)),
        mean_velocity_m_s=u.mean(0),
        relative_velocity_RMS_m_s=rms_vectors(velocity),
        fit_residual_RMS_m_s=rms_vectors(residual),
        fit_relative_residual=rms_vectors(residual) / max(rms_vectors(velocity), 1e-300),
    )


def force_summary(primary, secondary, force_scale):
    delta = primary - primary.mean(0)
    dp = rms_vectors(delta)
    ds = rms_vectors(secondary)
    return dict(
        mean_primary_N=primary.mean(0),
        D_P_N=dp,
        D_S_N=ds,
        force_scale_N=force_scale,
        D=dp / max(ds, force_scale),
        primary_RMS_N=rms_vectors(primary),
        primary_nonuniformity=dp / max(np.linalg.norm(primary.mean(0)), force_scale),
        secondary_net_N=secondary.sum(0),
        total_RMS_N=rms_vectors(primary + secondary),
    )
