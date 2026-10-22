"""One bounded miniature three-disc run tests CW map integration, not apparatus physics."""

import os

os.environ["KWAVE_F10L_PYTHON_ONLY"] = "1"
import time
import numpy as np
import h5py
from mevo_mb.acoustics.kwave_backend import KWaveCW
from mevo_mb.fast_sim.consistent_field import query_numpy
from .common import BASE, configs, write, sha, identity
from .preflight import reference_material


def main():
    folder = BASE / "source_bridge"
    folder.mkdir(exist_ok=False)
    ref = reference_material(configs()["materials"], 5e6)
    shape = [112] * 3
    h = 25e-6
    cm = np.full(shape, 1500, dtype="f4")
    rho = np.full(shape, 1000, dtype="f4")
    alpha = np.zeros(shape, dtype="f4")
    cm[56:] = ref["zero_frequency_sound_speed_m_s"]
    rho[56:] = ref["density_kg_m3"]
    alpha[56:] = ref["alpha_coeff_dB_MHz_power_cm"]
    path = folder / "validation_materials.h5"
    with h5py.File(path, "w") as f:
        for k, v in [("sound_speed_map", cm), ("density_map", rho), ("alpha_coeff_map", alpha)]:
            f.create_dataset(k, data=v, compression="gzip")
    elements = []
    for name, pos, n in [
        ("test_y_minus", [0, -0.00025, 0], [0, 1, 0]),
        ("test_y_plus", [0, 0.00025, 0], [0, -1, 0]),
        ("test_x_plus", [0.00025, 0, 0], [-1, 0, 0]),
    ]:
        elements.append(
            dict(
                id=name,
                position_m=pos,
                normal=n,
                geometry="disc",
                size_m=[0.0002],
                amplitude_Pa=0.0,
                phase_rad=0.0,
                enabled=True,
            )
        )
    cfg = dict(
        grid=dict(shape=shape, spacing_m=[h] * 3),
        medium=dict(maps_h5=str(path), maps_sha256=sha(path), alpha_power=ref["alpha_power"]),
        pzt=elements,
        acoustics=dict(
            frequency_Hz=5e6,
            CFL=0.2,
            transit_multiples=2.0,
            ramp_periods=8,
            recording_periods=4,
            roi_half_width_m=[0.0003] * 3,
            bli_tolerance=0.01,
            upsampling_rate=4,
            pml_size=16,
            pml_inside=False,
            pml_alpha=2,
            backend="python",
            device="gpu",
            BonA=None,
            allow_python_gpu_fallback=False,
        ),
    )
    write(folder / "CONFIG.json", cfg)
    t = time.perf_counter()
    adapter = KWaveCW(cfg, folder)
    occupied = np.argwhere(adapter.source_mask)
    clear = np.minimum(occupied.min(0), np.array(shape) - 1 - occupied.max(0))
    assert min(clear) >= 1, "Truncated BLI support"
    u = np.array([100, 50j, -30], complex)
    p, _, meta = adapter.run(u, "mixed_control")
    points = np.array([[1e-5, 2e-5, 3e-5], [-1.1e-5, 2.2e-5, 3.3e-5]])
    origin = [v[0] for v in adapter.roi_axes]
    spacing = [h] * 3
    pp, gg = query_numpy(p, origin, spacing, points)
    eps = 1e-9
    fd = np.stack(
        [
            (
                query_numpy(p, origin, spacing, points + np.eye(3)[j] * eps)[0]
                - query_numpy(p, origin, spacing, points - np.eye(3)[j] * eps)[0]
            )
            / (2 * eps)
            for j in range(3)
        ],
        axis=-1,
    )
    err = float(np.linalg.norm(fd - gg) / np.linalg.norm(gg))
    gate = dict(
        status="PASS",
        scope="MINIATURE_NUMERICAL_ADAPTER_TEST_ONLY: 0.2 mm discs, not 7 mm apparatus",
        source_count=3,
        actual_source_control_real_Pa=u.real,
        actual_source_control_imag_Pa=u.imag,
        edge_clearance_cells=clear,
        source_weights_sha256=adapter.weights_sha256,
        medium_sha256=sha(path),
        F6_consistent_gradient_relative_error=err,
        wall_s=time.perf_counter() - t,
        metadata=meta,
        identity=identity(),
    )
    if err >= 1e-6 or meta["recording_half_window_relative_L2"] >= 0.01:
        gate["status"] = "FAIL"
    write(folder / "GATE.json", gate)
    assert gate["status"] == "PASS"
    print(dict(status=gate["status"], wall_s=gate["wall_s"], edge=clear, gradient_error=err), flush=True)


if __name__ == "__main__":
    main()
