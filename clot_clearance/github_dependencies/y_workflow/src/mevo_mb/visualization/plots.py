"""Static scientific figures with PNG/PDF artifacts; original values, no smoothing."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def save(fig, path: Path, dpi: int = 300):
    import json
    from mevo_mb.validation.gates import validation_watermark

    run = path.parent.parent
    metadata = {}
    for candidate in (run / "metadata.json", run / "particles/SUMMARY.json"):
        if candidate.exists():
            metadata = json.loads(candidate.read_text())
            break
    status = dict(
        metadata.get("physics_validity", {}), **{k: v for k, v in metadata.items() if k != "physics_validity"}
    )
    engine = fig.get_layout_engine()
    if engine is not None and type(engine).__name__ == "ConstrainedLayoutEngine":
        engine.set(rect=(0, 0.07, 1, 0.93))
    else:
        fig.subplots_adjust(bottom=0.12)
    fig.text(
        0.5, 0.012, validation_watermark(status), ha="center", color="darkred", fontsize=10, weight="bold"
    )
    fig.savefig(path.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_superposition(run, axes, basis, direct, metrics, config):
    mid = np.argmin(np.abs(axes[0]))
    extent = [axes[2][0] * 1e3, axes[2][-1] * 1e3, axes[1][0] * 1e3, axes[1][-1] * 1e3]
    fig, axs = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    images = [np.abs(basis[mid]) / 1e3, np.abs(direct[mid]) / 1e3, np.abs(basis[mid] - direct[mid])]
    for ax, values, title in zip(
        axs, images, ["Basis superposition (kPa)", "Direct k-Wave (kPa)", "Absolute error (Pa)"], strict=True
    ):
        im = ax.imshow(values, origin="lower", extent=extent, aspect="equal")
        ax.set(xlabel="z (mm)", ylabel="y (mm)", title=title)
        fig.colorbar(im, ax=ax, shrink=0.75)
    fig.suptitle(f"Linear superposition: relative L2 = {metrics['relative_L2']:.3g}")
    save(fig, run / "figures/superposition", config["output"]["figure_dpi"])


def plot_experiment(run, field, config):
    import h5py
    import pandas as pd
    from mevo_mb.bubbles.forces import primary_bjerknes
    from mevo_mb.bubbles.response import LinearResonantBubbleModel

    axes = field.axes_m
    ix = np.argmin(np.abs(axes[0]))
    extent = [axes[2][0] * 1e3, axes[2][-1] * 1e3, axes[1][0] * 1e3, axes[1][-1] * 1e3]
    response = config["bubble_response"]
    if not config["particles"]["primary_enabled"]:
        volume = np.zeros(field.pressure_Pa[ix].shape, dtype=complex)
    elif response["model"] == "linear_surrogate":
        model = LinearResonantBubbleModel(response["resonance_frequency_Hz"], response["damping_ratio"])
        volume = model.response(
            config["particles"]["size_distribution"]["radius_m"],
            field.frequency_Hz,
            field.pressure_Pa[ix],
            config["medium"],
        )
    else:
        from mevo_mb.bubbles.response_table import ResponseTableBubbleModel, ResponseTableBoundsError
        from mevo_mb.bubbles.adaptive_table import AdaptiveResponseTableBubbleModel

        cls = (
            AdaptiveResponseTableBubbleModel
            if response["model"] == "adaptive_response_table"
            else ResponseTableBubbleModel
        )
        model = cls(response["path"])
        model.validate_medium(config["medium"])
        volume = np.full(field.pressure_Pa[ix].shape, np.nan + 1j * np.nan)
        for index in np.ndindex(volume.shape):
            try:
                volume[index] = model.response(
                    config["particles"]["size_distribution"]["radius_m"],
                    field.frequency_Hz,
                    field.pressure_Pa[ix][index],
                    config["medium"],
                )
            except ResponseTableBoundsError:
                pass  # Display a masked pixel; never extrapolate a force for plotting.
    force = np.linalg.norm(primary_bjerknes(volume, field.gradient_Pa_m[ix]), axis=-1)
    fig, axs = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    for ax, values, title, cmap in zip(
        axs,
        [np.abs(field.pressure_Pa[ix]) / 1e3, np.angle(field.pressure_Pa[ix]), force * 1e12],
        ["Local pressure amplitude (kPa)", "Complex pressure phase (rad)", "Primary force magnitude (pN)"],
        ["viridis", "twilight", "magma"],
        strict=True,
    ):
        im = ax.imshow(values, origin="lower", extent=extent, cmap=cmap, interpolation="nearest")
        ax.set(xlabel="z (mm)", ylabel="y (mm)", title=title)
        ax.add_patch(
            plt.Circle((0, 0), config["vessel"]["radius_m"] * 1e3, fill=False, color="white", lw=0.8)
        )
        fig.colorbar(im, ax=ax, shrink=0.8)
    save(fig, run / "figures/acoustic_and_primary_force", config["output"]["figure_dpi"])
    with h5py.File(run / "particles/trajectories.h5") as f:
        x = f["positions_m"][:]
        labels = f["cluster_labels"][-1]
        times = f["time_s"][:]
    fig = plt.figure(figsize=(13, 4), layout="constrained")
    ax = fig.add_subplot(131, projection="3d")
    take = np.unique(np.linspace(0, len(times) - 1, min(200, len(times))).astype(int))
    for i in range(x.shape[1]):
        ax.plot(*(x[take, i] * 1e3).T, lw=0.5, alpha=0.65)
    ax.set(xlabel="x (mm)", ylabel="y (mm)", zlabel="z (mm)", title="Accepted trajectories")
    for panel, data, title in [(132, x[0], "Initial distribution"), (133, x[-1], "Final spatial clusters")]:
        ax = fig.add_subplot(panel)
        ax.scatter(
            data[:, 1] * 1e6,
            data[:, 2] * 1e6,
            c=labels if panel == 133 else "tab:blue",
            s=9,
            cmap="tab20" if panel == 133 else None,
        )
        ax.set(xlabel="y (µm)", ylabel="z (µm)", title=title, aspect="equal")
    save(fig, run / "figures/trajectories_and_distribution", config["output"]["figure_dpi"])
    metrics = pd.read_csv(run / "metrics.csv")
    fig, axs = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    t_ms = metrics.time_s * 1e3
    axs[0, 0].plot(t_ms, metrics.mean_nearest_neighbor_m * 1e6, label="Mean")
    axs[0, 0].plot(t_ms, metrics.median_nearest_neighbor_m * 1e6, label="Median")
    axs[0, 0].set(ylabel="Nearest neighbor (µm)")
    axs[0, 0].legend()
    axs[0, 1].plot(t_ms, metrics.cluster_count)
    axs[0, 1].set(ylabel="Cluster count")
    axs[1, 0].plot(t_ms, metrics.largest_cluster_fraction, label="Largest cluster")
    axs[1, 0].plot(t_ms, metrics.near_contact_fraction, label="Near contact")
    axs[1, 0].set(ylabel="Population fraction", ylim=(-0.02, 1.02))
    axs[1, 0].legend()
    axs[1, 1].plot(t_ms, metrics.mean_lateral_displacement_m * 1e6)
    axs[1, 1].set(ylabel="Mean lateral displacement (µm)")
    for ax in axs.flat:
        ax.set(xlabel="Exposure time (ms)")
        ax.grid(alpha=0.25)
    fig.suptitle("Spatial clustering metrics; no fusion model")
    save(fig, run / "figures/metrics", config["output"]["figure_dpi"])
