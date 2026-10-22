"""Geometry, documented literature conversion and resource admission, no GPU solve."""

import itertools
import shutil
import numpy as np
import h5py
from scipy.fft import next_fast_len
from scipy.spatial.distance import pdist
from mevo_mb.coupling.pdms_block import apparatus, audit, material_labels, top_port_audit
from mevo_mb.coupling.straight_tube import disc_points, ZeroFlow
from mevo_mb.coupling.geometry import grid_axes
from mevo_mb.fast_sim.models import random_cloud, ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics
from .common import ROOT, BASE, configs, read, write, sha, identity


def reference_material(m, f):
    lo, hi = m["measured_frequency_range_Hz"]
    if not lo <= f <= hi:
        raise ValueError("PDMS literature frequency outside measured 3–7 MHz; no silent extrapolation")
    p = m["pdms"]
    ref = m["reference_frequency_Hz"]
    y = p["alpha_power"]
    alpha = p["alpha_reference_Np_m"]
    term = alpha * np.tan(np.pi * y / 2) / (2 * np.pi * ref)
    c0 = 1 / (1 / p["phase_speed_at_reference_m_s"] - term)
    cf = 1 / (1 / p["phase_speed_at_reference_m_s"] + term * ((f / ref) ** (y - 1) - 1))
    coeff = alpha * 20 / np.log(10) / 100 / (ref / 1e6) ** y
    if m["shared_alpha_power"] != y:
        raise ValueError("Incompatible regional power laws")
    for name in ["channel_liquid", "exterior"]:
        if m[name]["alpha_coeff"] != 0:
            raise ValueError("Nonzero liquid loss needs a compatible registered global power law")
    return dict(
        frequency_Hz=f,
        phase_speed_m_s=cf,
        zero_frequency_sound_speed_m_s=c0,
        density_kg_m3=p["density_kg_m3"],
        alpha_Np_m=alpha * (f / ref) ** y,
        alpha_coeff_dB_MHz_power_cm=coeff,
        alpha_power=y,
        alpha_mode=m["alpha_mode"],
        conversion="20/log(10), divide by 100 cm/m and by f_ref_MHz**y; c0 from Eq18 to avoid double dispersion",
        measured_frequency_range_Hz=[lo, hi],
        experimental_calibration=False,
    )


def build():
    c = configs()
    g, m, p = c["geometry"], c["materials"], c["protocol"]
    b = p["bubbles"]
    BASE.mkdir(parents=True, exist_ok=True)
    if (BASE / "CONFIG_SNAPSHOT.json").exists():
        assert read(BASE / "CONFIG_SNAPSHOT.json") == c, "Changed config requires separate results"
    lumen, elements = apparatus(g)
    geo = audit(g)
    geo["top_port"] = top_port_audit(g)
    ref = reference_material(m, p["full_field_candidate_frequency_Hz"])
    x = random_cloud(b["count"], b["radius_m"], b["cloud_radius_m"], b["seed"], b["cloud_center_m"])
    assert lumen.contains(x, b["radius_m"]).all() and pdist(x).min() > 2 * b["radius_m"]
    assert p["flow"] == {"mode": "zero"} and np.count_nonzero(ZeroFlow().query(x)[0]) == 0
    corners = np.array(list(itertools.product(*[[-a / 2, a / 2] for a in g["block_dimensions_m"]])))
    bounds = np.max(abs(np.vstack([corners, *[disc_points(e) for e in elements]])), axis=0)
    assert p["BLI_margin_cells"] >= int(np.ceil(1 / (np.pi * p["BLI_tolerance"]))) + 2
    h = p["grid_spacing_m"]
    pad = p["pml_size"]
    required = 2 * (np.ceil(bounds / h).astype(int) + p["BLI_margin_cells"] + 1)
    full = []
    for n in required + 2 * pad:
        v = next_fast_len(int(n))
        while v % 2:
            v = next_fast_len(v + 1)
        full.append(v)
    shape = np.array(full) - 2 * pad
    axes = grid_axes(shape.tolist(), [h] * 3)
    ppw = ref["zero_frequency_sound_speed_m_s"] / p["full_field_candidate_frequency_Hz"] / h
    assert ppw >= p["minimum_points_per_wavelength"]
    ntperiod = int(np.ceil(1500 / (p["CFL"] * h * p["full_field_candidate_frequency_Hz"])))
    nt = int(
        (
            12
            + np.ceil(
                2
                * np.linalg.norm(shape * h)
                / ref["zero_frequency_sound_speed_m_s"]
                * p["full_field_candidate_frequency_Hz"]
            )
            + 4
        )
        * ntperiod
    )
    voxels = int(np.prod(full))
    old = 378 * 192 * 192
    memory = 15100 / 1024 * voxels / old * 1.2  # actual F10L peak plus heterogeneous margin
    seconds = 373.70579146797536 * voxels / old * (nt / 3504) * 1.5
    reasons = []
    if memory > p["memory_budget_GiB"]:
        reasons.append("ESTIMATED_VRAM_EXCEEDS_20_GiB_BUDGET_AND_24_GiB_DEVICE")
    if seconds > p["full_field_budget_s"]:
        reasons.append("ESTIMATED_THREE_BASIS_TIME_EXCEEDS_600_S")
    reasons.extend(
        ["SOLID_END_PZT_ACOUSTICS_UNCALIBRATED", "TOP_PORT_AXIAL_REFERENCE_AND_BLIND_DEPTH_PROVISIONAL"]
    )
    if not geo["end_mounting_confirmed"]:
        reasons.append("END_PZT_MOUNTING_UNCONFIRMED")
    if not 2.499e6 <= p["full_field_candidate_frequency_Hz"] <= 2.501e6:
        reasons.append("NO_APECSS_RESPONSE_AT_CANDIDATE_FREQUENCY_FOR_ON")
    grid = dict(
        shape=shape.tolist(),
        spacing_m=[h] * 3,
        expanded_shape=full,
        axes_bounds_m=[[v[0], v[-1]] for v in axes],
        pml_outside=True,
        pml_cells=pad,
        bli_margin_cells=p["BLI_margin_cells"],
        bli_tolerance=p["BLI_tolerance"],
        source_BLI_support="BOUNDED_ANALYTIC_PREFLIGHT; full kWaveArray weights NOT_RUN after memory rejection",
        source_centers_and_rims_inside=True,
        minimum_ppw=ppw,
        dt_s=1 / (ntperiod * p["full_field_candidate_frequency_Hz"]),
        samples_per_period=ntperiod,
        Nt=nt,
        CFL_requested=p["CFL"],
        estimated_peak_GiB=memory,
        estimated_three_basis_seconds=seconds,
        estimate_method="Measured F10L15100MiB scaled by expanded voxels, 20% heterogeneity allowance; time scales by Nt with 50% allowance. Estimate not measured F10M apparatus usage.",
        frequency_Hz=p["full_field_candidate_frequency_Hz"],
    )
    out = BASE / "materials"
    out.mkdir(exist_ok=True)
    path = out / "material_maps.h5"
    if not path.exists():
        with h5py.File(path, "w") as h5:
            h5.attrs["material_ids"] = "0 exterior-water DEMO; 1 PDMS reference; 2 channel-liquid-water DEMO"
            for name, v in zip("xyz", axes):
                h5.create_dataset(name + "_m", data=v)
            datasets = {
                name: h5.create_dataset(
                    name,
                    shape=tuple(shape),
                    dtype=dtype,
                    compression="gzip",
                    shuffle=True,
                    chunks=(1, int(shape[1]), int(shape[2])),
                )
                for name, dtype in [
                    ("material_id", "u1"),
                    ("sound_speed_map", "f4"),
                    ("density_map", "f4"),
                    ("alpha_coeff_map", "f4"),
                ]
            }
            yy, zz = np.meshgrid(axes[1], axes[2], indexing="ij")
            for i, xx in enumerate(axes[0]):
                labels = material_labels(np.stack([np.full_like(yy, xx), yy, zz], axis=-1), g)
                datasets["material_id"][i] = labels
                datasets["sound_speed_map"][i] = np.where(
                    labels == 1, ref["zero_frequency_sound_speed_m_s"], 1500
                )
                datasets["density_map"][i] = np.where(labels == 1, ref["density_kg_m3"], 1000)
                datasets["alpha_coeff_map"][i] = np.where(labels == 1, ref["alpha_coeff_dB_MHz_power_cm"], 0)
            h5.attrs["alpha_power"] = ref["alpha_power"]
            h5.attrs["frequency_Hz"] = grid["frequency_Hz"]
            h5.attrs["source"] = "Xu2020 Tables II/III, 10:1 85C; liquids lossless DEMO"
    inputs = BASE / "inputs"
    inputs.mkdir(exist_ok=True)
    response = inputs / "response.h5"
    if not response.exists():
        shutil.copy2(ROOT / "results/straight_tube_exp1/f10l/inputs/response.h5", response)
    assert sha(response) == b["response_sha256"]
    ResponseSection(response, F2Physics(radius_m=b["radius_m"], frequency_Hz=b["frequency_Hz"]))
    np.savez_compressed(
        BASE / "initial_state.npz", positions_m=x, particle_id=np.arange(len(x)), radius_m=b["radius_m"]
    )
    write(BASE / "CONFIG_SNAPSHOT.json", c)
    write(BASE / "analysis/GEOMETRY.json", geo)
    write(BASE / "analysis/RESOURCE_PREFLIGHT.json", grid)
    write(BASE / "materials/REFERENCE.json", ref)
    write(
        BASE / "materials/IDENTITY.json",
        dict(path=str(path.relative_to(ROOT)), sha256=sha(path), grid_shape=shape),
    )
    write(
        BASE / "analysis/ACOUSTIC_ADMISSION.json",
        dict(
            status="NOT_RUN" if reasons else "ADMITTED",
            reasons=reasons,
            actual_apparatus_basis_count=0,
            old_F10L_fields_reused=False,
            training_ready=False,
            identity=identity(),
        ),
    )
    write(BASE / "INPUT_IDENTITY.json", identity())
    print(
        dict(shape=shape, expanded=full, ppw=ppw, estimated_GiB=memory, estimated_s=seconds, reasons=reasons),
        flush=True,
    )


if __name__ == "__main__":
    build()
