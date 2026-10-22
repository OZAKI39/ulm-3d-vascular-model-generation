"""Build/resume the focused adaptive database; retain independent holdout errors."""

import argparse
from pathlib import Path
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mevo_mb.bubbles.adaptive_table import (
    AdaptiveResponseTableBuilder,
    AdaptiveResponseTableBubbleModel,
    ResponseNodeCache,
    response_errors,
)
from common import new_run, finish, atomic_json, watermark


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--resume")
    p.add_argument("--profile", default="production")
    p.add_argument("--workers", type=int, default=6)
    a = p.parse_args()
    out = (
        Path(a.resume)
        if a.resume
        else new_run(
            "response_table",
            vars(a),
            ["APECSS v1.8; Tu2009 Fig4 mapping; Marmottant2005 reversible specialization"],
        )
    )
    b = AdaptiveResponseTableBuilder(
        out,
        "third_party/build/apecss/ultrasound_apecss",
        "configs/stage2/sonovue_parameters.json",
        profile=a.profile,
        workers=a.workers,
        cache="results/stage2_5/response_table/node_cache",
    )
    result = b.build()
    model = AdaptiveResponseTableBubbleModel(out / "adaptive_response.h5")
    rng = np.random.default_rng(20261003)
    cache = ResponseNodeCache(b.cache, b.exe)
    rows = []
    start = time.perf_counter()
    for point in rng.uniform(b.bounds[:, 0], b.bounds[:, 1], size=(128, 3)):
        actual = cache.get(b.case(point))
        if actual["status"] != "COMPLETE":
            rows.append(dict(point=point.tolist(), status="FAILED_NODE"))
            continue
        try:
            pred = model.query(*point)
            errors = response_errors(
                pred["harmonics_m3"],
                actual,
                point[1],
                float(pred["max_R_over_R0"]),
                float(pred["min_R_over_R0"]),
            )
            rows.append(
                dict(
                    point=point.tolist(),
                    status="PASS" if max(errors.values()) <= b.config["tolerance"] else "FAIL",
                    errors={k: float(v) for k, v in errors.items()},
                )
            )
        except ValueError as error:
            rows.append(dict(point=point.tolist(), status="REJECTED", error=str(error)))
    atomic_json(
        out / "HOLDOUTS.json",
        dict(seed=20261003, independent_uniform_queries=rows, runtime_s=time.perf_counter() - start),
    )
    failed = sum(r["status"] in ("FAIL", "FAILED_NODE") for r in rows)
    rejected = sum(r["status"] == "REJECTED" for r in rows)
    result.update(
        independent_holdout_failures=failed,
        independent_holdout_rejections=rejected,
        holdout_count=len(rows),
        physics_status="partially validated",
        table=str(out / "adaptive_response.h5"),
    )
    if failed:
        result["status"] = "FAILED_INDEPENDENT_HOLDOUT"
        for row in rows:
            if row["status"] not in ("FAIL", "FAILED_NODE"):
                continue
            point = np.asarray(row["point"])
            for cell in b.leaves:
                if np.all(point >= cell["lo"]) and np.all(point <= cell["hi"]):
                    cell.update(status="INVALID_CELL", reason="independent_holdout_failure")
        b.checkpoint()
        b.export()
        result["valid_leaves"] = sum(c["status"] == "VALID_CELL" for c in b.leaves)
        result["coverage_fraction"] = float(
            sum(np.prod(np.asarray(c["hi"]) - c["lo"]) for c in b.leaves if c["status"] == "VALID_CELL")
            / np.prod(b.bounds[:, 1] - b.bounds[:, 0])
        )
    before = time.perf_counter()
    again = AdaptiveResponseTableBuilder(
        out,
        b.exe,
        "configs/stage2/sonovue_parameters.json",
        profile=a.profile,
        workers=a.workers,
        cache=b.cache,
    ).build()
    result["resume_check"] = dict(
        new_apecss_calls=again["new_apecss_calls"], runtime_s=time.perf_counter() - before
    )
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
    leaves = b.leaves
    ax[0].hist([c["depth"] for c in leaves], bins=np.arange(14) - 0.5)
    ax[0].set(xlabel="Binary refinement depth", ylabel="Leaf cells")
    good = [r for r in rows if "errors" in r]
    for k in good[0]["errors"] if good else []:
        ax[1].plot([r["errors"][k] for r in good], ".", label=k)
    ax[1].axhline(b.config["tolerance"], color="k", ls="--")
    ax[1].set(xlabel="Independent holdout query", ylabel="Error (phase in rad)")
    ax[1].legend(fontsize=6)
    watermark(fig, "PARTIALLY VALIDATED — INTERPOLATION CHECK, NOT SONOVUE CALIBRATION")
    fig.savefig(out / "ADAPTIVE_COVERAGE.png", dpi=180)
    plt.close(fig)
    finish(out, result)


if __name__ == "__main__":
    main()
