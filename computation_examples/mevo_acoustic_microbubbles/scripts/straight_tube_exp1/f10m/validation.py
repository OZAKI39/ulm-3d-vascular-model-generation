"""Bounded real 3-D GPU propagation tests, distinct from full apparatus fields."""

import os

os.environ["KWAVE_F10L_PYTHON_ONLY"] = "1"
import time
import importlib.metadata
import numpy as np
from .common import BASE, configs, write, identity
from .preflight import reference_material
from mevo_mb.acoustics.spatial_medium import make_medium, medium_arrays
from mevo_mb.acoustics.kwave_backend import GPUTelemetry, extract_complex


def solve(folder, label, med):
    from kwave.kgrid import kWaveGrid
    from kwave.ksource import kSource
    from kwave.ksensor import kSensor
    from kwave.kspaceFirstOrder import kspaceFirstOrder

    shape = (192, 16, 16)
    h = 25e-6
    f = 5e6
    spp = 40
    dt = 1 / (spp * f)
    nt = 2400
    kg = kWaveGrid(shape, [h] * 3)
    kg.setTime(nt, dt)
    source = kSource()
    source.p_mask = np.zeros(shape, bool)
    source.p_mask[24, :, :] = True
    times = np.arange(nt) * dt
    ramp = 0.5 * (1 - np.cos(np.pi * np.minimum(times * f / 8, 1)))
    source.p = (100 * np.cos(2 * np.pi * f * times) * ramp)[None, :].astype("f4")
    source.p_mode = "additive"
    sensor = kSensor()
    sensor.mask = np.zeros(shape, bool)
    for x in [64, 96, 144]:
        sensor.mask[x, 8, 8] = True
    sensor.record = ["p"]
    sensor.record_start_index = nt - 8 * spp + 1
    start = time.perf_counter()
    result = kspaceFirstOrder(
        kg,
        make_medium(med, shape),
        source,
        sensor,
        backend="python",
        device="gpu",
        dtype="float32",
        pml_size=[20, 0, 0],
        pml_inside=False,
        pml_alpha=2,
        quiet=True,
        smooth_p0=False,
    )
    raw = result["p"]
    raw = raw.get() if hasattr(raw, "get") else np.asarray(raw)
    t = times[-8 * spp :]
    P = extract_complex(raw, t, f)
    delta = float(
        np.linalg.norm(
            extract_complex(raw[:, : 4 * spp], t[: 4 * spp], f)
            - extract_complex(raw[:, 4 * spp :], t[4 * spp :], f)
        )
        / np.linalg.norm(P)
    )
    np.savez_compressed(
        folder / (label + ".npz"),
        time_s=t,
        pressure_Pa=raw,
        P_complex=P,
        sensor_x_m=np.array([64, 96, 144]) * h,
    )
    return dict(phasor=P, elapsed_s=time.perf_counter() - start, half_window_relative=delta)


def main():
    folder = BASE / "adapter_validation"
    folder.mkdir(exist_ok=False)
    c = configs()
    ref = reference_material(c["materials"], 5e6)
    shape = (192, 16, 16)
    records = {}
    checks = {}
    start = time.perf_counter()
    with GPUTelemetry() as tel:
        records["scalar_water"] = solve(
            folder, "scalar_water", dict(sound_speed_m_s=1500.0, density_kg_m3=1000.0)
        )
        records["map_water"] = solve(
            folder,
            "map_water",
            dict(sound_speed_map=np.full(shape, 1500.0), density_map=np.full(shape, 1000.0)),
        )
        cm = np.full(shape, 1500.0)
        rho = np.full(shape, 1000.0)
        cm[120:] = ref["zero_frequency_sound_speed_m_s"]
        rho[120:] = ref["density_kg_m3"]
        records["lossless_interface"] = solve(
            folder, "lossless_interface", dict(sound_speed_map=cm, density_map=rho)
        )
        records["pdms_lossless"] = solve(
            folder,
            "pdms_lossless",
            dict(
                sound_speed_map=np.full(shape, ref["zero_frequency_sound_speed_m_s"]),
                density_map=np.full(shape, ref["density_kg_m3"]),
            ),
        )
        records["pdms_absorbing"] = solve(
            folder,
            "pdms_absorbing",
            dict(
                sound_speed_map=np.full(shape, ref["zero_frequency_sound_speed_m_s"]),
                density_map=np.full(shape, ref["density_kg_m3"]),
                alpha_coeff_map=np.full(shape, ref["alpha_coeff_dB_MHz_power_cm"]),
                alpha_power=ref["alpha_power"],
            ),
        )
        # True spatial attenuation map plus density / speed interface, not only constant maps.
        alpha = np.zeros(shape)
        alpha[120:] = ref["alpha_coeff_dB_MHz_power_cm"]
        records["absorbing_interface"] = solve(
            folder,
            "absorbing_interface",
            dict(sound_speed_map=cm, density_map=rho, alpha_coeff_map=alpha, alpha_power=ref["alpha_power"]),
        )
    scalar = records["scalar_water"]["phasor"]
    mapped = records["map_water"]["phasor"]
    scalar_error = float(np.linalg.norm(scalar - mapped) / np.linalg.norm(scalar))
    expected_R = (ref["zero_frequency_sound_speed_m_s"] * ref["density_kg_m3"] - 1.5e6) / (
        ref["zero_frequency_sound_speed_m_s"] * ref["density_kg_m3"] + 1.5e6
    )
    measured_R = (records["lossless_interface"]["phasor"][0] - scalar[0]) / scalar[0]
    dx = (96 - 64) * 25e-6
    lossless = records["pdms_lossless"]["phasor"]
    absorbing = records["pdms_absorbing"]["phasor"]
    loss_ratio = float(abs((absorbing[1] / absorbing[0]) / (lossless[1] / lossless[0])))
    expected_loss = float(np.exp(-ref["alpha_Np_m"] * dx))
    phase_error = float(
        np.angle(absorbing[1] / absorbing[0] * np.exp(2j * np.pi * 5e6 * dx / ref["phase_speed_m_s"]))
    )
    water_phase_error = float(np.angle(scalar[1] / scalar[0] * np.exp(2j * np.pi * 5e6 * dx / 1500)))
    checks.update(
        scalar_constant_map_relative_error_lt_1e_6=scalar_error < 1e-6,
        reflection_magnitude_absolute_error_lt_0p02=abs(abs(measured_R) - abs(expected_R)) < 0.02,
        power_law_amplitude_ratio_absolute_error_lt_0p025=abs(loss_ratio - expected_loss) < 0.025,
        positive_time_phasor_water_phase_error_lt_0p05_rad=abs(water_phase_error) < 0.05,
        pdms_reference_phase_error_lt_0p12_rad=abs(phase_error) < 0.12,
        all_outputs_finite=all(np.isfinite(r["phasor"]).all() for r in records.values()),
        steady_half_window_relative_lt_0p01=max(r["half_window_relative"] for r in records.values()) < 0.01,
    )
    # Errors must be rejected before the backend silently reads only the first exponent.
    for name, med in [
        (
            "spatial_exponent",
            dict(sound_speed_m_s=1500, density_kg_m3=1000, alpha_power=np.full(shape, 1.494)),
        ),
        ("wrong_map_shape", dict(sound_speed_map=np.ones((4, 4, 4)) * 1500, density_kg_m3=1000)),
        ("negative_density", dict(sound_speed_m_s=1500, density_kg_m3=-1)),
        (
            "unit_exponent_dispersion",
            dict(sound_speed_m_s=1500, density_kg_m3=1000, alpha_coeff=0.1, alpha_power=1),
        ),
    ]:
        try:
            medium_arrays(med, shape)
        except ValueError:
            checks[name + "_rejected"] = True
        else:
            checks[name + "_rejected"] = False
    summary = dict(
        status="PASS" if all(checks.values()) else "FAIL",
        scope="3-D GPU planar verification, transverse periodic with zero transverse PML; NOT the PDMS apparatus field",
        frequency_Hz=5e6,
        grid_shape=shape,
        spacing_m=25e-6,
        PML=[20, 0, 0],
        samples_per_period=40,
        Nt=2400,
        checks=checks,
        scalar_map_relative_error=scalar_error,
        reflection_expected_pressure_coefficient=expected_R,
        reflection_measured_magnitude=abs(measured_R),
        attenuation_ratio_expected=expected_loss,
        attenuation_ratio_measured=loss_ratio,
        pdms_phase_error_rad=phase_error,
        water_phase_error_rad=water_phase_error,
        actual_wall_s=time.perf_counter() - start,
        sampled_device_peak_MiB=max(s["used_MiB"] for s in tel.samples),
        runs={
            k: dict(
                elapsed_s=v["elapsed_s"],
                half_window_relative=v["half_window_relative"],
                P_real=v["phasor"].real,
                P_imag=v["phasor"].imag,
            )
            for k, v in records.items()
        },
        environment={p: importlib.metadata.version(p) for p in ["k-wave-python", "cupy-cuda12x", "numpy"]},
        identity=identity(),
    )
    write(folder / "VALIDATION.json", summary)
    write(folder / "GPU_TELEMETRY.json", dict(samples=tel.samples, errors=tel.errors))
    print(summary, flush=True)
    if not all(checks.values()):
        raise RuntimeError("Bounded medium validation failed; preserve outputs")


if __name__ == "__main__":
    main()
