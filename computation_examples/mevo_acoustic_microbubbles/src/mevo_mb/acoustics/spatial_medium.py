"""Strict spatial-medium input bridge for audited k-wave-python 0.6.2 Python backend."""

import hashlib
from pathlib import Path
import numpy as np
import h5py


def medium_arrays(m, shape):
    if "maps_h5" in m:
        path = Path(m["maps_h5"])
        if hashlib.sha256(path.read_bytes()).hexdigest() != m["maps_sha256"]:
            raise ValueError("Material-map checksum mismatch")
        with h5py.File(path) as h:
            c = h["sound_speed_map"][:]
            rho = h["density_map"][:]
            alpha = h["alpha_coeff_map"][:]
    else:
        c = np.asarray(m.get("sound_speed_map", m.get("sound_speed_m_s")), float)
        rho = np.asarray(m.get("density_map", m.get("density_kg_m3")), float)
        alpha = np.asarray(m.get("alpha_coeff_map", m.get("alpha_coeff", 0.0)), float)
    for name, v, positive in [
        ("sound_speed", c, True),
        ("density", rho, True),
        ("alpha_coeff", alpha, False),
    ]:
        if v.ndim and v.shape != tuple(shape):
            raise ValueError(name + " map shape must exactly equal grid shape")
        if not np.isfinite(v).all() or np.any(v <= 0 if positive else v < 0):
            raise ValueError("Invalid " + name)
    power = np.asarray(m.get("alpha_power", 1.5))
    if power.ndim or not np.isfinite(power) or not 0 < float(power) < 3:
        raise ValueError("Installed backend requires ONE global scalar alpha_power")
    mode = m.get("alpha_mode")
    if mode not in (None, "no_dispersion", "no_absorption"):
        raise ValueError("Unverified absorption mode")
    if np.any(alpha) and abs(float(power) - 1) < 1e-6 and mode != "no_dispersion":
        raise ValueError("Dispersive alpha_power=1 is singular")
    return c, rho, alpha, float(power), mode


def make_medium(m, shape):
    from kwave.kmedium import kWaveMedium

    c, rho, alpha, power, mode = medium_arrays(m, shape)
    # Preserve historical lossless scalar pathway exactly.
    kwargs = dict(sound_speed=c, density=rho)
    if np.any(alpha):
        kwargs.update(alpha_coeff=alpha, alpha_power=power, alpha_mode=mode)
    return kWaveMedium(**kwargs)


def speed_bounds(m, shape):
    c, _, _, _, _ = medium_arrays(m, shape)
    return float(np.min(c)), float(np.max(c))
