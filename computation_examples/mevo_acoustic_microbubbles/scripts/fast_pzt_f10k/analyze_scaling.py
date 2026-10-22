"""Aggregate only actual completed same-quality measurements, retain startup failures."""

import numpy as np
from .common import BASE, ROOT, js, sha, write, identity
from .scaling import state


def main():
    result = js(BASE / "scaling_v3/RESULTS.json")
    for row in result["rows"]:
        n = row["N"]
        folder = BASE / f"scaling_v3/N{n}"
        if not folder.exists():
            assert row["status"] == "NOT_RUN"
            continue
        meta = js(folder / "INPUT.json")
        saved = np.load(folder / "INPUT.npz")
        x, free, reconstructed = state(n)
        np.testing.assert_array_equal(saved["positions_m"], x)
        np.testing.assert_array_equal(saved["free_velocity_m_s"], free)
        assert meta["contact_pairs"] == reconstructed["contact_pairs"]
        assert meta["npz_sha256"] == sha(folder / "INPUT.npz")
        row.update(
            input_contact_pairs=meta["contact_pairs"],
            input_label=meta["label"],
            construction=js(folder / "CONSTRUCTION.json"),
            input_sha256=sha(folder / "INPUT.npz"),
        )
        assert row["status"] == "PASS" and row["comparison"]["pass"]
        for m in row["measurements"]:
            assert m["quality_pass"] and m["converged"]
            assert m["metrics"]["equation_residual_m"] <= 1e-14
            assert m["metrics"]["max_overlap_m"] <= 17.8e-9
            assert m["metrics"]["single_correction_m"] <= 178e-9
            if m["method"] == "gpu":
                a = np.array(m["components"]["linear_audit"])
                assert np.all(a[:, :2] <= a[:, 2, None])
            else:
                c = m["components"]
                assert c["base_true_residual"] <= c["base_target"]
                assert all(v["pass_"] for v in c["direction_linear_audit"])
        for m in ["cpu", "gpu"]:
            trials = [r for r in row["measurements"] if r["method"] == m]
            assert sum(t["kind"] == "warmup" for t in trials) == 1
            assert sum(t["kind"] == "formal" for t in trials) == 2
            assert all(
                b["num_threads"] == 1 for b in row["construction"]["actual_BLAS"] if b["user_api"] == "blas"
            )
        c = np.load(folder / "cpu_latest.npz")
        g = np.load(folder / "gpu_latest.npz")
        assert np.linalg.norm(c["positions_m"] - g["positions_m"], axis=1).max() <= 1e-14
        assert np.max(abs(c["contact_force_N"] - g["contact_force_N"])) <= 1e-15
        np.testing.assert_array_equal(c["active_pairs"], g["active_pairs"])
    proc = js(BASE / "scaling_v3_PROCESS.json")
    result["full_supervised_process_s"] = proc["full_process_wall_s"]
    samples = [
        int(line.split(",")[-1].strip().split()[0])
        for r in proc["resource_samples"]
        for line in r["processes"].splitlines()
        if line and line.split(",")[-1].strip().split()[0].isdigit()
    ]
    result["sampled_GPU_process_peak_MiB"] = max(samples, default=0)
    result["failed_startup_process_s"] = {
        p.stem: js(p)["full_process_wall_s"] for p in BASE.glob("scaling_v[12]_PROCESS.json")
    }
    result["source"] = "results/fast_pzt/f10k/scaling_v3/RESULTS.json"
    result["source_sha256"] = sha(ROOT / result["source"])
    result["execution_commit"] = js(BASE / "scaling_v3/IDENTITY.json")["execution_commit"]
    result["analysis_identity"] = identity()
    result["complexity_limit"] = (
        "Unchanged CPU all-pair contact B has O(N^3) storage; regular symmetric512 has only4 PCG iterations, unlike real32 needing132. No universal crossover inferred."
    )
    result["initialization_regression"] = (
        "32 opt-in GPU matches archived F10J bitwise in all3 calls; original OptimizedB2(128) still rejects."
    )
    result["largest_measured_N"] = max(r["N"] for r in result["rows"] if r["status"] == "PASS")
    assert result["solve_and_warmup_s"] < 120 and result["stage_wall_s"] < 300
    write(BASE / "analysis/SCALING.json", result)
    print("Three sizes,18 solves: quality and budget audit PASS")


if __name__ == "__main__":
    main()
