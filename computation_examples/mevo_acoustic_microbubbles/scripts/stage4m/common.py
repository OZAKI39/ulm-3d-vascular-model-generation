"""Stage 4M immutable references and read-only analysis helpers."""

from pathlib import Path
import json
import h5py
import numpy as np
from mevo_mb.bubbles.metrics import aggregation_metrics
from scripts.stage4l.dense_probe import sha, write

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/stage4/secondary_ablation"
DENSE_A0 = dict(
    B2="20261005T214817062274Z_B2_lubrication_coupled_5ms",
    B3="20261005T214817162616Z_B3_lubrication_coupled_5ms",
)
PARTIAL_A0 = dict(B2="20261005T223817_397200Z_B2_dtfactor1", B3="20261005T223818_038339Z_B3_dtfactor1")


def refs(legacy, stage4l, case):
    return dict(
        L0=legacy / "results/stage4/performance_refactor/gold_baseline/cases" / case,
        A0=stage4l / "results/stage4/lubrication_upgrade/secondary_short" / PARTIAL_A0[case],
        dense_A0=stage4l / "results/stage4/lubrication_upgrade/dense_probe" / DENSE_A0[case],
    )


def formation_metrics(x, radii, ids, config):
    if len(x) < 2:
        raise ValueError("Formation analysis requires at least two active particles")
    row, labels, _, _ = aggregation_metrics(x, radii, config)
    sizes = np.bincount(labels)
    largest = labels == np.argmax(sizes)
    centroid = x[largest].mean(axis=0)
    row.update(
        active_particles=len(x),
        cluster_size_histogram={str(k): int(np.sum(sizes == k)) for k in np.unique(sizes)},
        largest_cluster_particle_ids=ids[largest].tolist(),
        largest_cluster_centroid_m=centroid.tolist(),
        largest_cluster_radius_of_gyration_m=float(
            np.sqrt(np.mean(np.sum((x[largest] - centroid) ** 2, axis=1)))
        ),
    )
    for k in (3, 5, 10):
        n = int(np.sum(sizes[sizes >= k]))
        row[f"particles_in_clusters_ge_{k}"] = n
        row[f"fraction_in_clusters_ge_{k}"] = n / len(x)
    return row


def exact_frame(path, time_s):
    with h5py.File(path, "r") as h:
        times = h["time_s"][:]
        indices = np.flatnonzero(abs(times - time_s) < 1e-10)
        if len(indices) != 1:
            raise ValueError(f"No unique saved frame at {time_s}: {path}")
        k = indices[0]
        active = h["active"][k]
        return (
            h["positions_m"][k][active],
            h["radii_m"][:][active],
            h["particle_id"][:][active],
            float(times[k]),
        )


def recorded_pair_diagnostics(tracker, n):
    if n == 0:
        from mevo_mb.bubbles.secondary_ablation import pair_norm_diagnostics

        empty = np.empty((0, 0))
        return pair_norm_diagnostics(
            empty, empty, empty, empty, 1.78e-6, tracker.config["medium"]["density_kg_m3"]
        )
    pair = tracker.pair
    if hasattr(pair, "cpu"):
        pair = pair.cpu if n < pair.crossover or pair.gpu_disabled else pair.gpu
    return pair.pair_diagnostics()


def reference_analysis(legacy, stage4l):
    out = BASE / "l0_vs_a0"
    out.mkdir(parents=True, exist_ok=True)
    rows = {}
    identities = {}
    for case in ("B2", "B3"):
        paths = refs(legacy, stage4l, case)
        config = json.loads((paths["L0"] / "RUN_IDENTITY.json").read_text())["config"]["metrics"]
        aconfig = json.loads((paths["A0"] / "RUN_IDENTITY.json").read_text())["config"]["metrics"]
        assert config == aconfig
        rows[case] = {}
        initial = []
        for policy in ("L0", "A0"):
            trajectory = paths[policy] / "particles/trajectories.h5"
            x, r, ids, saved_t = exact_frame(trajectory, 0.08)
            rows[case][policy] = dict(time_s=saved_t, **formation_metrics(x, r, ids, config))
            initial.append(exact_frame(trajectory, 0)[:3])
            identities[f"{case}_{policy}"] = dict(
                trajectory=str(trajectory), sha256=sha(trajectory), metrics_config=config
            )
        for a, b in zip(*initial):
            np.testing.assert_array_equal(a, b)
    note = "This comparison isolates the effect of the Stage 4L translation change before any secondary ablation."
    result = dict(time_s=0.08, statement=note, rows=rows, inputs=identities, matching_initial_population=True)
    write(out / "LUBRICATION_ONLY_FORMATION_EFFECT.json", result)
    md = [
        "# Lubrication-only formation effect",
        "",
        note,
        "",
        "Exact shared HDF5 saved frame: 0.08 s. No simulation was run. Initial positions, IDs and radii match bitwise.",
        "",
        "| Case | Model | Largest cluster | Cluster count | Near-contact fraction | Median NN (µm) |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for case, policies in rows.items():
        for policy, row in policies.items():
            md.append(
                f"| {case} | {policy} | {row['largest_cluster_size']} | {row['cluster_count']} | {row['near_contact_fraction']:.4f} | {row['median_nearest_neighbor_m'] * 1e6:.3f} |"
            )
    md += [
        "",
        "Cluster threshold: 10 µm center distance; near-contact particle fraction: minimum surface gap ≤0.5 µm, unchanged from the stored Stage 4 configuration. Dense same-ID diagnostics separately use the inherited 2.1R graph definition. Full distributions, centroids, gyration radii and membership counts are in the JSON.",
        "",
        "L0→A0 is the translation effect. A0→A2 will be the incremental ablation effect. L0→A2 combines both changes.",
    ]
    (out / "LUBRICATION_ONLY_FORMATION_EFFECT.md").write_text("\n".join(md) + "\n")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(9, 3.1), constrained_layout=True)
    for ax, key, label in zip(
        axes,
        ["largest_cluster_size", "cluster_count", "near_contact_fraction"],
        ["Largest cluster (bubbles)", "Number of clusters", "Near-contact fraction"],
    ):
        for j, (policy, color) in enumerate([("L0", "#777777"), ("A0", "#6b4c9a")]):
            ax.bar(
                np.arange(2) + (j - 0.5) * 0.34,
                [rows[c][policy][key] for c in ("B2", "B3")],
                0.32,
                label=policy,
                color=color,
            )
        ax.set_xticks([0, 1], ["B2", "B3"])
        ax.set_ylabel(label)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylim(0, 26)
    axes[1].set_ylim(0, 95)
    axes[2].set_ylim(0, 1.05)
    axes[0].legend(frameon=False, ncol=2, loc="upper right")
    fig.suptitle("Lubrication-only formation effect at 0.08 s")
    for suffix in ("png", "pdf"):
        fig.savefig(out / f"Figure_M1_lubrication_only.{suffix}", dpi=300)
    plt.close(fig)
    print(json.dumps(rows), flush=True)
    return result


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--legacy-root", type=Path, required=True)
    p.add_argument("--stage4l-root", type=Path, required=True)
    a = p.parse_args()
    reference_analysis(a.legacy_root, a.stage4l_root)
