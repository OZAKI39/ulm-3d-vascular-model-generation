import argparse
from time import perf_counter
from mevo_mb.fast_sim.f10e_isolation import isolate
from mevo_mb.fast_sim.f10d_metrics import process
from mevo_mb.fast_sim.f10c_families import control_families
from scripts.fast_pzt_f10c.runner import run
from .common import BASE, CAMPAIGN, identity, original_records, objects, write, js
from .prefix import compare


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    args = ap.parse_args()
    out = BASE / args.label
    out.mkdir(parents=True, exist_ok=False)
    source = identity()
    write(out / "IDENTITY.json", source)
    started = perf_counter()
    obj = objects()
    static = perf_counter() - started
    trees = []
    for control in ["S_1000", "O_1000"]:
        records = original_records(control)

        def execute(group, parent):
            family = control_families(group, max_batch=64)[0]
            name = control + "_" + "_".join(r["split_group"].split("_")[-1] for r in group)
            path = out / name
            begin = perf_counter()
            timing = run(
                path,
                family["positions_m"],
                family["events"],
                obj,
                sample_ids=family["sample_ids"],
                postprocess=process,
                identity_provider=lambda: source,
                dataset_metadata=dict(
                    initial_state_family=[r["initial_state_family"] for r in group],
                    split_groups=[r["split_group"] for r in group],
                    split_group=group[0]["split_group"],
                    labels=["EXPLORATORY", "PHYSICS_UNCALIBRATED"],
                    independent_scientific_initial_states=len(group),
                    initial_state_sources=[r["initial_state_source"] for r in group],
                    original_physics_mode="M0_DT10_PBD12",
                    postprocessing_mode="F10D_optimized",
                ),
            )
            timing["complete_run_wall_s"] = perf_counter() - begin
            write(path / "PIPELINE_TIMING.json", timing)
            summary = js(path / "SUMMARY.json")
            assert len(group) == 1  # Actual F10D parents each have two members.
            regression = compare(
                CAMPAIGN / control / "trajectory.h5",
                path / "trajectory.h5",
                group[0]["initial_state_source"]["parent_member"],
            )
            write(path / "PREFIX_REGRESSION.json", regression)
            result = dict(
                attempt_id=name,
                status=summary["status"],
                reason=summary["reason"],
                timing=timing,
                per_sample_quality=summary["per_sample_quality"],
                accepted_end_s=regression["common_accepted_end_s"],
                prefix_regression=regression["status"],
            )
            if regression["status"] != "PASS":
                result["status"] = "UNRESOLVED"
            write(
                path / "RETRY_PROVENANCE.json",
                dict(result, parent=parent, initial_state_sources=[r["initial_state_source"] for r in group]),
            )
            print(result, flush=True)
            return result

        trees.append(isolate(records, "F10D/campaign_v2_diversity/" + control, execute))
        write(
            out / "RESULTS.json",
            dict(
                trees=trees,
                static_prepare_s=static,
                process_body_s=perf_counter() - started,
                training_ready=False,
            ),
        )


if __name__ == "__main__":
    main()
