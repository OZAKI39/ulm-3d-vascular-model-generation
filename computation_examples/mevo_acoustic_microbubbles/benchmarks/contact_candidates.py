"""Candidate-only benchmark; production contact projection is unchanged."""

import itertools
import json
from pathlib import Path
import sys
import time
import tracemalloc

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/stage2_5"))
from common import resumable_run, finish, atomic_json, watermark


def brute(x, cutoff):
    pairs = []
    for i in range(0, len(x), 128):
        d = x[i : i + 128, None, :] - x[None, :, :]
        a, b = np.nonzero(np.sum(d * d, axis=-1) <= cutoff**2)
        a += i
        good = a < b
        pairs.extend(zip(a[good].tolist(), b[good].tolist()))
    return np.asarray(pairs, dtype=int).reshape(-1, 2)


def bins(x, cutoff):
    cells = {}
    for i, loc in enumerate(np.floor(x / cutoff).astype(int)):
        cells.setdefault(tuple(loc), []).append(i)
    return cells


def cell_query(x, cutoff, cells):
    out = []
    offsets = list(itertools.product([-1, 0, 1], repeat=3))
    for key, indices in cells.items():
        for d in offsets:
            other = tuple(k + v for k, v in zip(key, d))
            if other < key or other not in cells:
                continue
            a, b = np.meshgrid(indices, cells[other], indexing="ij")
            a = a.ravel()
            b = b.ravel()
            good = (a < b) if other == key else np.ones(len(a), bool)
            good &= np.sum((x[a] - x[b]) ** 2, axis=1) <= cutoff**2
            out.extend(zip(np.minimum(a[good], b[good]).tolist(), np.maximum(a[good], b[good]).tolist()))
    return np.asarray(out, int).reshape(-1, 2)


def run_case(n, density, seed):
    rng = np.random.default_rng(seed)
    r = rng.uniform(0.8e-6, 1.2e-6, n)
    x = rng.uniform(0, density * 1e-6 * n ** (1 / 3), size=(n, 3))
    cutoff = 2 * r.max()
    rows = []
    reference = None
    for method in ["brute", "ckdtree", "uniform_cell"]:
        tracemalloc.start()
        start = time.perf_counter()
        structure = (
            cKDTree(x) if method == "ckdtree" else (bins(x, cutoff) if method == "uniform_cell" else None)
        )
        built = time.perf_counter()
        candidates = (
            brute(x, cutoff)
            if method == "brute"
            else (
                structure.query_pairs(cutoff, output_type="ndarray")
                if method == "ckdtree"
                else cell_query(x, cutoff, structure)
            )
        )
        end = time.perf_counter()
        peak = tracemalloc.get_traced_memory()[1]
        tracemalloc.stop()
        a, b = candidates.T
        exact = candidates[np.linalg.norm(x[a] - x[b], axis=1) <= r[a] + r[b]]
        actual = set(map(tuple, exact.tolist()))
        if reference is None:
            reference = actual
        missed = len(reference - actual)
        extra = len(actual - reference)
        row = dict(
            N=n,
            spacing_scale_um=density,
            seed=seed,
            method=method,
            build_s=built - start,
            query_s=end - built,
            candidate_pairs=len(candidates),
            exact_contacts=len(actual),
            missed_contacts=missed,
            extra_contacts=extra,
            traced_peak_bytes=peak,
            memory_note="tracemalloc Python/NumPy tracked allocations, not full process RSS",
            exact=missed == 0 and extra == 0,
        )
        rows.append(row)
        print(n, density, method, row["query_s"], missed, flush=True)
    return rows


def main():
    start = time.perf_counter()
    out = resumable_run(
        "contact_benchmark",
        dict(N=[128, 256, 512, 1000, 2000, 5000, 10000], spacing_scale_um=[2.1, 4.0, 8.0], seed=20261002),
        ["Geometric contact candidate exactness; no physical adhesion model"],
    )
    rows = json.loads((out / "CHECKPOINT.json").read_text()) if (out / "CHECKPOINT.json").exists() else []
    for n in [128, 256, 512, 1000, 2000, 5000, 10000]:
        for density in [2.1, 4.0, 8.0]:
            if any(r["N"] == n and r["spacing_scale_um"] == density for r in rows):
                continue
            rows.extend(run_case(n, density, 20261002 + n))
            atomic_json(out / "CHECKPOINT.json", rows)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    for kind in ["brute", "ckdtree", "uniform_cell"]:
        rr = [r for r in rows if r["method"] == kind and r["spacing_scale_um"] == 2.1]
        ax.loglog([r["N"] for r in rr], [r["build_s"] + r["query_s"] for r in rr], "o-", label=kind)
    ax.set(xlabel="Particles", ylabel="Build + query (s)", title="Contact candidates, dense random fixture")
    ax.legend()
    watermark(fig, "SOFTWARE BENCHMARK — ZERO MISSED CONTACTS REQUIRED")
    fig.savefig(out / "CONTACT_SEARCH.png", dpi=180)
    finish(
        out,
        dict(
            status="PASS_EXACTNESS" if all(r["exact"] for r in rows) else "FAILED",
            rows=rows,
            runtime_s=time.perf_counter() - start,
            production_solver_changed=False,
        ),
    )


if __name__ == "__main__":
    main()
