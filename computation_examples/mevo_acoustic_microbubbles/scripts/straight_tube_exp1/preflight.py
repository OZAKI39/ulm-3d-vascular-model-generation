"""Deterministic apparatus, source-grid, response identity, and budget admission."""

import shutil
import numpy as np
from scipy.fft import next_fast_len
from scipy.spatial.distance import pdist
from mevo_mb.coupling.straight_tube import apparatus, geometry_audit, disc_points, ZeroFlow
from mevo_mb.coupling.geometry import grid_axes
from mevo_mb.fast_sim.models import random_cloud, ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics
from .common import BASE, ROOT, configs, write, sha, identity


def build():
    c = configs()
    g, a, b = c["demo"]["geometry"], c["acoustic"], c["bubbles"]
    if c["demo"]["flow"] != {"mode": "zero"}:
        raise ValueError("Explicit zero flow only")
    tube, elements = apparatus(g, a["acoustics"]["frequency_Hz"])
    audit = geometry_audit(tube, elements, g["gap_m"])
    if not all(e["gap_matches_registered"] for e in audit["pzt"]):
        raise ValueError("Source geometry disagrees with registered gaps")
    cloud = random_cloud(b["count"], b["radius_m"], b["cloud_radius_m"], b["seed"], b["cloud_center_local_m"])
    x = tube.world(cloud)
    if b["radius_distribution"] != "monodisperse":
        raise ValueError("M0 currently supports equal radii only; no silent distribution replacement")
    assert tube.contains(x, b["radius_m"]).all()
    assert pdist(x).min() >= 2.2 * b["radius_m"] * (1 - 1e-12)
    surface = tube.surface().reshape(-1, 3)
    allpoints = np.vstack([surface, *[disc_points(e) for e in elements]])
    h = a["grid_design"]["spacing_m"]
    pad = a["acoustics"]["pml_size"]
    bound = np.max(abs(allpoints), axis=0)
    required = 2 * (np.ceil(bound / h).astype(int) + a["grid_design"]["bli_margin_cells"] + 1)
    # Choose even FFT-friendly *expanded* dimensions, then remove the outside PML.
    full = []
    for n in required + 2 * pad:
        m = next_fast_len(int(n))
        while m % 2:
            m = next_fast_len(m + 1)
        full.append(m)
    shape = np.asarray(full) - 2 * pad
    roi = (np.ceil(np.max(abs(surface), axis=0) / h).astype(int) + 2) * h
    cfg = dict(
        name=c["demo"]["name"],
        grid=dict(shape=shape.tolist(), spacing_m=[h] * 3),
        medium=a["medium"],
        acoustics=dict(a["acoustics"], roi_half_width_m=roi.tolist()),
        pzt=elements,
    )
    axes = grid_axes(shape.tolist(), [h] * 3)
    samples_per_period = int(
        np.ceil(a["medium"]["sound_speed_m_s"] / (a["acoustics"]["CFL"] * h * a["acoustics"]["frequency_Hz"]))
    )
    transit = int(
        np.ceil(
            a["acoustics"]["transit_multiples"]
            * np.linalg.norm(shape * h)
            / 1500
            * cfg["acoustics"]["frequency_Hz"]
        )
    )
    nt = (a["acoustics"]["ramp_periods"] + transit + a["acoustics"]["recording_periods"]) * samples_per_period
    voxels = int(np.prod(full))
    # Historical 4090 measurement: 5.25 GiB, 25.91 s at 288x160x160, Nt=2712.
    ratio = voxels / (288 * 160 * 160)
    estimate_gib = 5.25 * ratio * 1.3
    estimate_s = 3 * 25.913232 * ratio * (nt / 2712) * 1.4 + 60
    if estimate_gib > a["grid_design"]["maximum_estimated_vram_GiB"] or estimate_s > 600:
        raise ValueError("Registered acoustic budget exceeded; do not silently coarsen grid")
    checks = dict(
        geometry=audit,
        background_velocity_exact_zero=bool(np.all(ZeroFlow().query(x)[0] == 0)),
        finite_radius_initial_inside=bool(tube.contains(x, b["radius_m"]).all()),
        minimum_initial_pair_distance_m=float(pdist(x).min()),
        tube_bounds_m=[[float(surface[:, i].min()), float(surface[:, i].max())] for i in range(3)],
        source_envelope_bounds_m=[[-float(v), float(v)] for v in bound],
        inner_diameter_m=2 * tube.radius_m,
        outer_diameter_m=2 * tube.outer_radius_m,
        length_m=tube.length_m,
        wall_thickness_m=tube.thickness_m,
        experimental_wall_thickness_m=None,
        ends=list(tube.ends),
        grid=cfg["grid"],
        pml_outside=True,
        expanded_shape=full,
        wavelength_m=1500 / cfg["acoustics"]["frequency_Hz"],
        points_per_wavelength=1500 / cfg["acoustics"]["frequency_Hz"] / h,
        grid_convergence="NOT_ESTABLISHED",
        samples_per_period=samples_per_period,
        Nt=nt,
        recording_periods=4,
        estimated_peak_vram_GiB=estimate_gib,
        estimated_total_acoustic_s=estimate_s,
        estimate_basis="Historical Stage 3 RTX 4090 runtime and memory scaled by voxel count and Nt; includes 40% timing margin, not a guarantee",
        grid_bounds_m=[[float(axis[0]), float(axis[-1])] for axis in axes],
        injection_interference="UNKNOWN",
        identity=identity(),
    )
    BASE.mkdir(parents=True, exist_ok=True)
    write(BASE / "CONFIG_SNAPSHOT.json", c)
    write(BASE / "RESOLVED_ACOUSTIC_CONFIG.json", cfg)
    write(BASE / "analysis/PREFLIGHT.json", checks)
    np.savez_compressed(
        BASE / "initial_state.npz", positions_m=x, particle_id=np.arange(b["count"]), radius_m=b["radius_m"]
    )
    dest = ROOT / b["response_path"]
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / "results/fast_pzt/f5/inputs/f4/response.h5", dest)
    assert sha(dest) == b["response_sha256"]
    response = ResponseSection(
        dest, F2Physics(radius_m=b["radius_m"], frequency_Hz=cfg["acoustics"]["frequency_Hz"])
    )
    write(BASE / "inputs/RESPONSE_IDENTITY.json", response.metadata)
    write(
        BASE / "INPUT_IDENTITY.json",
        dict(
            identity(), response_sha256=sha(dest), historical_basis_not_reused=True, old_CFD_not_required=True
        ),
    )
    print(
        {
            k: checks[k]
            for k in ["grid", "expanded_shape", "estimated_peak_vram_GiB", "estimated_total_acoustic_s"]
        }
    )
    return cfg


if __name__ == "__main__":
    build()
