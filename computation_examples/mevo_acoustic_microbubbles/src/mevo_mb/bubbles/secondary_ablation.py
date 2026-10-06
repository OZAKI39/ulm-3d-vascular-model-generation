"""Stage 4M sensitivity policies; these are NOT corrected near-field physics."""

import numpy as np

BANNER = "Sensitivity-only smooth removal of unresolved classical secondary interaction inside the lubrication zone."
POLICIES = ("A0", "A1", "A2")
FLOOR = 0.001
ACTIVATION = 0.1


def near_contact_secondary_ablation(gap_m, radius_m, policy="A0", xp=np):
    """One symmetric gap weight per unordered pair, with fixed Stage 4L scales."""
    if policy not in POLICIES:
        raise ValueError(f"Unknown secondary ablation policy: {policy}")
    h = xp.asarray(gap_m)
    if policy == "A0":
        return xp.ones_like(h, dtype=float)
    if policy == "A1":
        return (h > FLOOR * radius_m).astype(float)
    s = xp.clip((h - FLOOR * radius_m) / ((ACTIVATION - FLOOR) * radius_m), 0, 1)
    # Factorization avoids powers; explicit endpoint branches give exact 0/1.
    w = s * s * s * (10 + s * (-15 + 6 * s))
    return xp.where(h <= FLOOR * radius_m, 0.0, xp.where(h >= ACTIVATION * radius_m, 1.0, w))


def pair_norm_diagnostics(d2, correlation, weights, gaps, radius_m, density, xp=np):
    """Sum magnitudes once per unordered pair, never interpret as net force."""
    upper = xp.triu(xp.ones(d2.shape, dtype=bool), k=1)
    classical = density / (4 * np.pi) * xp.abs(correlation[upper]) / d2[upper]
    w = weights[upper]
    total = float(xp.sum(classical))
    effective = float(xp.sum(classical * w))
    count = int(w.size)
    full, zero = int(xp.count_nonzero(w == 1)), int(xp.count_nonzero(w == 0))
    taper = count - full - zero
    ratio = effective / total if total else None
    return dict(
        S_classical_N=total, S_effective_N=effective,
        effective_classical_ratio=ratio,
        removed_classical_force_norm_fraction=None if ratio is None else 1 - ratio,
        unordered_pair_count=count, full_weight_pair_count=full,
        taper_zone_pair_count=taper, zeroed_pair_count=zero,
        full_weight_pair_fraction=full / count if count else None,
        taper_zone_pair_fraction=taper / count if count else None,
        zeroed_pair_fraction=zero / count if count else None,
        subfloor_pair_count=int(xp.count_nonzero(gaps[upper] <= FLOOR * radius_m)),
        semantics="Unordered pair-force magnitude sums, not net physical force",
    )
