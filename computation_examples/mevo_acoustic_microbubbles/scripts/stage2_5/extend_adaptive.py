"""Extend an existing adaptive tree without repeating completed APECSS nodes.

Old run is immutable. New budgets and old holdouts used for refinement are recorded;
a fresh random set is held out after the extension.
"""

import argparse
import json
from pathlib import Path
import time

import numpy as np
from mevo_mb.bubbles.adaptive_table import (
    AdaptiveResponseTableBuilder,
    AdaptiveResponseTableBubbleModel,
    ResponseNodeCache,
    CORNERS,
    response_errors,
)
from common import new_run, finish, atomic_json


def main():
    p = argparse.ArgumentParser()
    p.add_argument("parent")
    p.add_argument("--resume", type=Path)
    p.add_argument("--max-nodes", type=int, default=24000)
    p.add_argument("--max-depth", type=int, default=20)
    p.add_argument("--seed", type=int, default=20261004)
    a = p.parse_args()
    parent = Path(a.parent)
    tree = json.loads((parent / "CHECKPOINT.json").read_text())
    old_nodes = len(tree["nodes"])
    if a.resume:
        out = a.resume
        config = json.loads((out / "CONFIG.json").read_text())["config"]
        expected = dict(parent=str(parent), max_nodes=a.max_nodes, max_depth=a.max_depth, holdout_seed=a.seed)
        if config != expected:
            raise ValueError("Extension resume configuration mismatch")
        tree = json.loads((out / "CHECKPOINT.json").read_text())
        extension = out / "EXTENSION.json"
        reopened = json.loads((extension if extension.exists() else out / "SUMMARY.json").read_text())[
            "reopened_leaves"
        ]
    else:
        out = new_run(
            "response_table",
            dict(parent=str(parent), max_nodes=a.max_nodes, max_depth=a.max_depth, holdout_seed=a.seed),
            ["APECSS v1.8; unchanged Tu2009 parameter manifold"],
        )
        tree["config"].update(max_nodes=a.max_nodes, max_depth=a.max_depth)
        hold = json.loads((parent / "HOLDOUTS.json").read_text())
        failures = [r["point"] for r in hold["independent_uniform_queries"] if r["status"] == "FAIL"]
        leaves = []
        reopened = 0
        for cell in tree["leaves"]:
            lo, hi = np.asarray(cell["lo"]), np.asarray(cell["hi"])
            failed = any(np.all(np.array(q) >= lo) & np.all(np.array(q) <= hi) for q in failures)
            if failed or cell.get("reason") == "independent_holdout_failure":
                for corner in CORNERS:
                    aa = lo + corner * (hi - lo) / 2
                    tree["queue"].append(
                        dict(lo=aa.tolist(), hi=(aa + (hi - lo) / 2).tolist(), depth=cell["depth"] + 3)
                    )
                reopened += 1
            elif cell.get("reason") in ["node_budget", "max_depth"]:
                tree["queue"].append({k: cell[k] for k in ["lo", "hi", "depth"]})
                reopened += 1
            else:
                leaves.append(cell)
        tree["leaves"] = leaves
        atomic_json(out / "CHECKPOINT.json", tree)
        atomic_json(out / "EXTENSION.json", dict(reopened_leaves=reopened))
    b = AdaptiveResponseTableBuilder(
        out,
        "third_party/build/apecss/ultrasound_apecss",
        "configs/stage2/sonovue_parameters.json",
        profile=tree["config"]["profile"],
        max_nodes=a.max_nodes,
        max_depth=a.max_depth,
        workers=7,
        cache="results/stage2_5/response_table/node_cache",
    )
    result = b.build()
    model = AdaptiveResponseTableBubbleModel(out / "adaptive_response.h5")
    cache = ResponseNodeCache(b.cache, b.exe)
    rng = np.random.default_rng(a.seed)
    rows = []
    start = time.perf_counter()
    for point in rng.uniform(b.bounds[:, 0], b.bounds[:, 1], (256, 3)):
        actual = cache.get(b.case(point))
        if actual["status"] != "COMPLETE":
            rows.append(dict(point=point.tolist(), status="FAILED_NODE"))
            continue
        try:
            q = model.query(*point)
            e = response_errors(
                q["harmonics_m3"], actual, point[1], float(q["max_R_over_R0"]), float(q["min_R_over_R0"])
            )
            rows.append(
                dict(
                    point=point.tolist(),
                    status="PASS" if max(e.values()) <= b.config["tolerance"] else "FAIL",
                    errors={k: float(v) for k, v in e.items()},
                )
            )
        except ValueError as error:
            rows.append(dict(point=point.tolist(), status="REJECTED", error=str(error)))
    atomic_json(
        out / "HOLDOUTS.json",
        dict(seed=a.seed, independent_uniform_queries=rows, runtime_s=time.perf_counter() - start),
    )
    failed = [r for r in rows if r["status"] in ["FAIL", "FAILED_NODE"]]
    if failed:
        for row in failed:
            point = np.array(row["point"])
            for cell in b.leaves:
                if np.all(point >= cell["lo"]) and np.all(point <= cell["hi"]):
                    cell.update(status="INVALID_CELL", reason="independent_holdout_failure")
        b.checkpoint()
        b.export()
        result["status"] = "FAILED_INDEPENDENT_HOLDOUT"
    result.update(
        parent=str(parent),
        reopened_leaves=reopened,
        reused_existing_tree_nodes=old_nodes,
        independent_holdout_failures=len(failed),
        independent_holdout_rejections=sum(r["status"] == "REJECTED" for r in rows),
        holdout_count=len(rows),
        table=str(out / "adaptive_response.h5"),
        physics_status="partially validated",
    )
    result["coverage_fraction"] = float(
        sum(np.prod(np.array(c["hi"]) - c["lo"]) for c in b.leaves if c["status"] == "VALID_CELL")
        / np.prod(b.bounds[:, 1] - b.bounds[:, 0])
    )
    result["valid_leaves"] = sum(c["status"] == "VALID_CELL" for c in b.leaves)
    resumed = AdaptiveResponseTableBuilder(
        out,
        b.exe,
        "configs/stage2/sonovue_parameters.json",
        profile=tree["config"]["profile"],
        max_nodes=a.max_nodes,
        max_depth=a.max_depth,
        workers=7,
        cache=b.cache,
    ).build()
    result["resume_check"] = {
        "new_apecss_calls": resumed["new_apecss_calls"],
        "runtime_s": resumed["runtime_s"],
    }
    finish(out, result)


if __name__ == "__main__":
    main()
