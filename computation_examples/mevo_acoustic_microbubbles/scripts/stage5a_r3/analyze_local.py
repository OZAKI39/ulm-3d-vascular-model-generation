"""Same-ID local refinement and exclusive dynamic-profile summaries."""

from pathlib import Path
import pstats
import h5py
import numpy as np
from scripts.stage5a_r3.common import BASE, R2, read, write, sha
from scripts.stage5a_r3.lineage import sample_all, metrics, summarize


def analyze():
    refinements = []
    for path in sorted((BASE / "diagnosis").glob("LOCAL_REFINEMENT_B*.json")):
        rows = read(path)
        if len(rows) < 3:
            continue
        trajectories = [Path(r["run"]) / "trajectories.h5" for r in rows]
        limits = []
        for p in trajectories:
            with h5py.File(p) as h:
                limits.append((float(h["time_s"][0]), float(h["time_s"][-1])))
        start = limits[0][0]
        end = min(b for a, b in limits)
        if "initial" in path.name:
            end = min(end, start + 1e-5)
        a = [sample_all(p, np.array([start, end])) for p in trajectories]
        assert all(
            np.array_equal(x["ids"], a[0]["ids"]) and np.array_equal(x["positions"][0], a[0]["positions"][0])
            for x in a
        )
        case = rows[0]["summary"]["case"]
        group = np.load(R2 / "comparisons" / f"{case}_0.001_v4a_legacy_reference.npz")["particle_ids"]
        errors = []
        for left, right in [(0, 1), (1, 2), (0, 2)]:
            groups = {}
            for name, mask in [
                ("all_active", np.ones(len(a[0]["ids"]), bool)),
                ("original_group", np.isin(a[0]["ids"], group)),
                ("complement", ~np.isin(a[0]["ids"], group)),
            ]:
                norms = metrics(a[left]["positions"][:, mask], a[right]["positions"][:, mask])
                groups[name] = {mode: summarize(v[-1]) for mode, v in norms.items()}
                groups[name]["worst_particle_id"] = int(a[0]["ids"][mask][np.argmax(norms["raw"][-1])])
            errors.append(
                dict(
                    left=left,
                    right=right,
                    groups=groups,
                    sampling_bound_m=float(a[left]["uncertainty_m"][-1] + a[right]["uncertainty_m"][-1]),
                )
            )
        old, new = [v["groups"]["all_active"]["raw"]["rms"] for v in errors[:2]]
        order = float(np.log2(old / new)) if new > 0 and old > 0 else None
        report = dict(
            name=path.stem,
            case=case,
            relative_duration_s=end - start,
            actual_end_times_s=[float(x["time"][-1]) for x in a],
            comparisons=errors,
            observed_rms_refinement_order=order,
            accuracy_status="NOT_EVALUABLE",
            reason="Observed three-level local refinement is not an independent global error bound",
            runs=rows,
        )
        refinements.append(report)
        save = {"particle_id": a[0]["ids"], "clock_s": np.array([start, end])}
        save.update({f"level_{i}_positions_m": x["positions"] for i, x in enumerate(a)})
        np.savez_compressed(BASE / "comparisons" / (path.stem + ".npz"), **save)
    profile = []
    for path in sorted((BASE / "runs").glob("*/DYNAMIC_PROFILE.pstats")):
        data = pstats.Stats(str(path))
        functions = []
        for (file, line, name), (primitive, calls, self_s, total_s, callers) in data.stats.items():
            functions.append(
                dict(
                    file=file,
                    line=line,
                    function=name,
                    primitive_calls=primitive,
                    calls=calls,
                    exclusive_self_s=self_s,
                    inclusive_s=total_s,
                )
            )
        profile.append(
            dict(
                run=str(path.parent),
                source_sha256=sha(path),
                total_function_calls=data.total_calls,
                total_profile_s=data.total_tt,
                functions=sorted(functions, key=lambda x: x["exclusive_self_s"], reverse=True),
                timing_scope="cProfile changes overhead; diagnostic only. Exclusive self times partition profiler elapsed time, inclusive times must not be added",
            )
        )
    methods = []
    for path in sorted((BASE / "runs").glob("*/SUMMARY.json")):
        summary = read(path)
        accepted = read(path.parent / "ACCEPTED.jsonl")
        for method in ["S0", "S1"]:
            rows = [r for r in accepted if r["substep_method"] == method]
            methods.append(
                dict(
                    run=str(path.parent),
                    case=summary["case"],
                    method=method,
                    accepted_substeps=len(rows),
                    physical_duration_s=sum(r["dt_s"] for r in rows),
                    dt_quantiles_s=np.quantile([r["dt_s"] for r in rows], [0, 0.5, 1]).tolist()
                    if rows
                    else None,
                    first_time_s=rows[0]["time_s"] if rows else None,
                    last_time_s=rows[-1]["time_s"] if rows else None,
                    error_max=max((r["error_ratio"] for r in rows), default=None),
                    nonlinear_kkt_status="NOT_APPLICABLE" if method == "S0" else "EVALUATED",
                    nonlinear_iterations=sum(r["nonlinear_iterations"] for r in rows),
                    qp_iterations=sum(r["qp_iterations"] for r in rows),
                    qp_cost_s={
                        key: sum(r["qp_timings"][key] for r in rows)
                        for key in ["setup_s", "update_s", "solve_s", "polish_s"]
                    },
                )
            )
    write(
        BASE / "final/S0_S1_USAGE_AND_REFINEMENT.json",
        dict(
            local_refinement=refinements,
            usage=methods,
            policy="S0 caps independently halved 50/25/12.5 ns; S0 and S1 remain different methods; no extra recovery bursts",
        ),
    )
    write(
        BASE / "final/DYNAMIC_COST_PROFILE.json",
        dict(
            profiles=profile,
            accepted_method_costs=methods,
            other_costs="Each SUMMARY.json includes total wall time, force/assembly/QP/certificate, I/O/checkpoint and monitoring inclusive fields; do not sum overlapping fields",
        ),
    )
    print(
        [
            (
                r["name"],
                r["observed_rms_refinement_order"],
                r["comparisons"][1]["groups"]["all_active"]["raw"],
            )
            for r in refinements
        ],
        flush=True,
    )


if __name__ == "__main__":
    analyze()
