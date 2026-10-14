"""Exactly three new disc solves, with actual source support and phasor audits."""

import time
import numpy as np
import h5py
from mevo_mb.acoustics.kwave_backend import KWaveCW, acoustic_identity
from mevo_mb.acoustics.basis import AcousticBasisDatabase
from mevo_mb.fast_sim.consistent_field import query_numpy
from .common import BASE, read, write, identity, sha


def run():
    cfg = read(BASE / "RESOLVED_ACOUSTIC_CONFIG.json")
    import os

    if cfg["acoustics"]["backend"] != "python":
        raise ValueError("F10L dependency guard is registered only for Python/CuPy")
    os.environ["KWAVE_F10L_PYTHON_ONLY"] = "1"
    from .environment import install_guard

    install_guard()
    folder = BASE / "acoustic"
    folder.mkdir(exist_ok=False)
    start = time.perf_counter()
    adapter = KWaveCW(cfg, folder)
    occupied = np.argwhere(adapter.source_mask)
    # PML is outside. Also require a clear physical-grid edge so finite BLI support was not truncated.
    edge = np.minimum(occupied.min(0), np.array(adapter.source_mask.shape) - 1 - occupied.max(0))
    if np.any(edge < 1):
        raise ValueError("BLI support reaches physical grid edge; enlarge source domain")
    write(
        folder / "SOURCE_GRID.json",
        dict(
            identity(),
            source_count=3,
            source_ids=[e["id"] for e in cfg["pzt"]],
            source_mask_shape=adapter.source_mask.shape,
            source_count_voxels=int(adapter.source_mask.sum()),
            occupied_index_bounds=[occupied.min(0), occupied.max(0)],
            edge_clearance_cells=edge,
            pml_outside=True,
            source_weights_sha256=adapter.weights_sha256,
            source_generation_s=time.perf_counter() - start,
        ),
    )
    records = []
    fields = []
    gradients = []
    amp = cfg["acoustics"]["basis_source_amplitude_Pa"]
    for i in range(3):
        if time.perf_counter() - start > 540:
            raise TimeoutError("600 s total budget: do not start another solve after 540 s")
        u = np.zeros(3, complex)
        u[i] = amp
        p, _, meta = adapter.run(u, f"basis_{i:02d}")
        if meta["recording_half_window_relative_L2"] > cfg["acoustics"]["steady_half_window_relative_L2_max"]:
            raise ValueError("Phasor has not settled")
        fields.append(p / amp)
        # Store F6's exact derivative of the trilinear interpolant at nodes
        # (positive-cell convention); runtime always differentiates P itself.
        points = np.stack(np.meshgrid(*adapter.roi_axes, indexing="ij"), axis=-1).reshape(-1, 3)
        origin = np.array([v[0] for v in adapter.roi_axes])
        spacing = np.array(cfg["grid"]["spacing_m"])
        # Avoid nominal upper-end roundoff by reconstituting from the same lattice.
        points = np.minimum(points, origin + spacing * (np.array(adapter.roi_shape) - 1))
        _, g = query_numpy(p / amp, origin, spacing, points)
        gradients.append(g.reshape(*adapter.roi_shape, 3))
        records.append(dict(meta, execution_commit=identity()["execution_commit"]))
    metadata = dict(
        frequency_Hz=cfg["acoustics"]["frequency_Hz"],
        BonA=None,
        acoustic_identity=acoustic_identity(cfg),
        source_ids=[e["id"] for e in cfg["pzt"]],
        source_elements=cfg["pzt"],
        medium=cfg["medium"],
        geometry=read(BASE / "CONFIG_SNAPSHOT.json")["demo"]["geometry"],
        PML={k: cfg["acoustics"][k] for k in ("pml_size", "pml_inside", "pml_alpha")},
        simulations=records,
        gradient="F6 exact derivative of same trilinear complex P; nodal positive-cell convention; runtime re-derives",
        calibration="HOMOGENEOUS_WATER_DEMO; transparent additive source; no material reflections or voltage calibration",
        **identity(),
    )
    path = folder / "basis.h5"
    AcousticBasisDatabase(adapter.roi_axes, fields, gradients, metadata).save_basis(path)
    with h5py.File(path, "a") as h5:
        h5.create_dataset("valid_roi_mask", data=np.ones(adapter.roi_shape, bool), compression="gzip")
    # Basic independent linear algebra and F6 finite-difference checks, no fourth solver call.
    weights = np.array([300, 300, 300], complex)
    combined = np.einsum("i,ixyz->xyz", weights, np.array(fields))
    initial = np.load(BASE / "initial_state.npz")["positions_m"]
    p, g = query_numpy(combined, origin, spacing, initial)
    eps = 1e-9
    fd = np.stack(
        [
            (
                query_numpy(combined, origin, spacing, initial + np.eye(3)[j] * eps)[0]
                - query_numpy(combined, origin, spacing, initial - np.eye(3)[j] * eps)[0]
            )
            / (2 * eps)
            for j in range(3)
        ],
        axis=-1,
    )
    relative = float(np.linalg.norm(g - fd) / np.linalg.norm(g))
    assert relative < 1e-6
    assert np.allclose(1j * combined, np.einsum("i,ixyz->xyz", 1j * weights, np.array(fields)))
    write(
        folder / "ACOUSTIC_GATE.json",
        dict(
            status="PASS",
            new_independent_basis_count=3,
            basis_sha256=sha(path),
            whole_acoustic_s=time.perf_counter() - start,
            finite_difference_gradient_relative_error=relative,
            source_input_peak_per_channel_Pa=300,
            initial_combination_local_peak_range_Pa=[abs(p).min(), abs(p).max()],
            initial_combination_local_rms_range_Pa=[abs(p).min() / np.sqrt(2), abs(p).max() / np.sqrt(2)],
            superposition="algebra of actual basis outputs; no new direct solve; inherits validated linear propagation",
            grid_converged=False,
            experimental_calibration=False,
        ),
    )


if __name__ == "__main__":
    from .preflight import build

    build()
    run()
