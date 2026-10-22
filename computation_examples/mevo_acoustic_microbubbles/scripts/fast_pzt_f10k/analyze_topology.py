"""Reuse registered F10H morphology gates and independently rebuild every edge."""

from time import perf_counter
import h5py
import numpy as np
from scipy import linalg
from mevo_mb.fast_sim.f10g_contact import resistance, normals
from scripts.fast_pzt_f10h.analyze import compare_paths, geometric_series, gz
from .common import ROOT, BASE, identity, js, write, sha


def main():
    started = perf_counter()
    paths = [BASE / "topology_dt10_v1/trial", BASE / "topology_dt5_v1/trial"]
    comparison = compare_paths(*paths)
    audits = []
    saved = {}
    for tag, path in zip(["DT10", "DT5"], paths):
        data = geometric_series(path)
        rows = gz(path / "SOFT_STEPS.json.gz")
        i, j = np.triu_indices(32, 1)
        d = np.linalg.norm(data["x"][:, j] - data["x"][:, i], axis=-1)
        linked = d <= 2.1 * 1.78e-6
        expected = sorted(
            (round(float(data["times"][f + 1]), 12), int(i[p]), int(j[p]), int(linked[f + 1, p]))
            for f, p in np.argwhere(linked[1:] != linked[:-1])
        )
        actual = sorted(
            (round(float(t), 12), int(i), int(j), int(on))
            for t, (_, i, j, on) in zip(data["edge_times"], data["edge_ids"])
        )
        assert expected == actual
        residuals = []
        mismatch = []
        overlaps = []
        corrections = []
        for row in rows:
            assert row["accepted"] and not row["error_bits"]
            x = np.array(row["initial_positions_m"])
            r = row["result"]
            y = np.array(r["positions_m"])
            np.testing.assert_array_equal(y, row["actual_accepted_positions_m"])
            a, b, dist, n = normals(y)
            ov = np.maximum(3.56e-6 - dist, 0)
            force = np.zeros_like(y)
            pair = 0.04641319518378658 * ov[:, None] * n
            np.add.at(force, a, -pair)
            np.add.at(force, b, pair)
            A = resistance(x, 1.78e-6)["A"]
            v = linalg.solve(
                A,
                np.array(r["free_m_s"]).ravel() + force.ravel() / (6 * np.pi * 0.001 * 1.78e-6),
                assume_a="pos",
            )
            residuals.append(float(np.linalg.norm(y - x - r["dt_s"] * v.reshape(-1, 3), axis=1).max()))
            mismatch.append(float(np.max(abs(force - np.array(r["contact_force_N"])))))
            overlaps.append(float(ov.max()))
            corrections.append(r["max_single_iteration_m"])
        with h5py.File(path / "trajectory.h5") as h:
            wall = float(h["trajectory/wall_stats"][:, :, :, 3].max())
        audit = dict(
            tag=tag,
            path=str(path.relative_to(ROOT)),
            trajectory_sha256=sha(path / "trajectory.h5"),
            accepted_steps=len(rows),
            reconstructed_edges=len(expected),
            edge_ledger_exact=True,
            max_equation_residual_m=max(residuals),
            max_force_geometry_difference_N=max(mismatch),
            max_overlap_m=max(overlaps),
            max_single_correction_m=max(corrections),
            contact_quality_pass=bool(
                max(residuals) <= 1.01e-13
                and max(overlaps) <= 17.8e-9
                and max(corrections) <= 178e-9
                and max(mismatch) < 1e-18
            ),
            wall_stat_column3_max=wall,
            topology=data["topology"],
            timing=js(path / "TIMING.json"),
            source_identity=js(path.parent / "IDENTITY.json"),
        )
        assert audit["contact_quality_pass"]
        audits.append(audit)
        for key in ["times", "x", "labels", "edge_times", "edge_ids"]:
            saved[f"{tag}_{key}"] = data[key]
        saved[f"{tag}_equation_residual_m"] = residuals
    np.savez_compressed(BASE / "analysis/TOPOLOGY_SERIES.npz", **saved)
    write(
        BASE / "analysis/B2_TOPOLOGY_COMPARISON.json",
        dict(
            comparison=comparison,
            audits=audits,
            identity=identity(),
            analysis_wall_s=perf_counter() - started,
            statement="One cross-model initialized interior scenario. Geometric graph merging is not physical bubble fusion. OldF10H/I strict12/29FAIL remains unchanged.",
            training_ready=False,
        ),
    )
    print(comparison)


if __name__ == "__main__":
    main()
