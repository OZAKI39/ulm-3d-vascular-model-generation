"""Read-only scientific summary/catalog, including every rejected-interval record."""

from collections import Counter
import json
import h5py
import numpy as np
from .common import BASE, ROOT, js, sha, write, identity


def main():
    out = BASE / "analysis"
    out.mkdir(exist_ok=True)
    catalog = []
    failures = []
    for folder in ["campaign_v2_benchmark", "campaign_v2_families", "campaign_v2_diversity"]:
        for path in sorted((BASE / folder).glob("*/trajectory.h5")):
            with h5py.File(path) as h:
                meta = json.loads(h.attrs["metadata_json"])
                assert not h.attrs["training_ready"]
                for group in [
                    "source",
                    "medium",
                    "flow",
                    "bubbles",
                    "control",
                    "trajectory",
                    "morphology",
                    "quality",
                ]:
                    assert group in h, (path, group)
                np.testing.assert_array_equal(h["particle_id"][:], np.arange(32))
                assert h["positions_m"].dtype == np.dtype("float64")
                ids = list(h["source/sample_id"].asstr()[:])
                assert ids == meta["sample_id"]
                config = js(path.parent / "CONFIG.json")
                times = h["time_s"][:]
                timing = js(path.parent / "PIPELINE_TIMING.json")
                for i, sid in enumerate(ids):
                    row = dict(
                        sample_id=sid,
                        batch_id=meta["batch_id"],
                        member=i,
                        file=str(path.relative_to(ROOT)),
                        sha256=sha(path),
                        source_commit=meta["source"]["execution_commit"],
                        control_family_id=meta["control_family_id"],
                        initial_state_family=meta["initial_state_family"][i],
                        split_group=meta["split_groups"][i],
                        initial_state_source=meta["initial_state_sources"][i],
                        original_physics_mode=meta["original_physics_mode"],
                        sample_class="EXPLORATORY_NUMERICAL_SAMPLE"
                        if folder.endswith("diversity")
                        else "BENCHMARK_ONLY",
                        execution_status=meta["status"],
                        numerical_quality="NUMERICAL_PASS"
                        if meta["status"] == "COMPLETE"
                        else "ABORTED_BATCH",
                        physical_duration_s=float(times[-1]),
                        requested_duration_s=config["duration_s"],
                        stored_frames=len(times),
                        particle_count=32,
                        preparation_time_s=timing["model_prepare_s"],
                        integration_time_s=timing["integration_quality_wall_s"],
                        morphology_time_s=timing["morphology_quality_output_prepare_s"],
                        output_time_s=timing["hdf5_output_s"] + timing["serialization_and_tail_s"],
                        total_elapsed_time_s=timing["complete_run_wall_s"],
                        valid_samples_per_gpu_hour=timing["valid_samples_per_gpu_hour"],
                        timing_scope="Per batch, not per member; static/process overhead recorded in parent RESULTS",
                        physics_scope="PHYSICS_UNCALIBRATED",
                        training_ready=False,
                    )
                    assert (
                        row["split_group"] == row["initial_state_family"]
                        or row["sample_class"] == "BENCHMARK_ONLY"
                    )
                    catalog.append(row)
                if meta["status"] == "FAILED":
                    assert "rejected_interval" in h and meta["valid_sample_count"] == 0
                    group = h["rejected_interval"]
                    proposal = group["proposal_m"][:]
                    evidence = []
                    for i, sid in enumerate(ids):
                        d = np.linalg.norm(proposal[i, :, None, :] - proposal[i, None, :, :], axis=-1)
                        pairs = np.argwhere(np.triu(d < 1.99 * h["radii_m"][0], 1))
                        evidence.append(
                            dict(
                                sample_id=sid,
                                violating_pairs=pairs.tolist(),
                                gap_m=[float(d[a, b] - 2 * h["radii_m"][0]) for a, b in pairs],
                                attribution="OBSERVED_REJECTED_PROPOSAL_CONTACT"
                                if len(pairs)
                                else "COLLATERAL_ABORT_NO_LOCAL_VIOLATION_OBSERVED",
                            )
                        )
                    failures.append(
                        dict(
                            path=str(path.relative_to(ROOT)),
                            reason=meta["reason"],
                            original_failures=meta["failures"],
                            evidence=evidence,
                            first_error_time_bracket_s=[
                                float(group.attrs["start_s"]),
                                float(group.attrs["attempted_end_s"]),
                            ],
                            exact_first_step_time_known=False,
                            scope="Offline rejected-proposal evidence, not an accepted state or proof of future validity; GPU global bits do not record the first emitting lane.",
                        )
                    )
    write(
        out / "SAMPLE_CATALOG.json",
        dict(
            records=catalog,
            timing_scope="All timing fields describe the containing batch",
            training_ready=False,
        ),
    )
    write(
        out / "FAILURE_AUDIT.json",
        dict(
            status="RECORDED_NOT_FIXED",
            error_2048_definition="vessel_kernels.final_geometry_check: distance < 1.99 * radius",
            source_sha256=sha(ROOT / "src/mevo_mb/fast_sim/vessel_kernels.py"),
            rows=failures,
            affected_trajectories=sum(len(f["evidence"]) for f in failures),
            collateral_trajectories_without_observed_local_violation=sum(
                e["attribution"].startswith("COLLATERAL") for f in failures for e in f["evidence"]
            ),
            limitation="50 us host-check intervals; no fabricated substep time; original quality/status is never overwritten; no retry or discarded failure.",
            training_ready=False,
        ),
    )

    val = js(BASE / "validate_v1/RESULTS.json")
    bench = js(BASE / "campaign_v2_benchmark/RESULTS.json")
    family = js(BASE / "campaign_v2_families/RESULTS.json")
    diversity = js(BASE / "campaign_v2_diversity/RESULTS.json")
    checks = [c for r in val["rows"] for c in r["checks"]]
    checks += [c for r in diversity["comparisons"] for c in r["checks"]]
    checks += [c for r in bench["rows"] for c in r.get("comparison", {}).get("checks", [])]
    checks += [c for r in family["comparisons"] for c in r.get("checks", [])]
    maxima = {}
    for c in checks:
        assert c["status"] == "PASS", c
        if not c.get("bitwise", True):
            maxima[c["key"]] = max(maxima.get(c["key"], 0), c.get("max_absolute_error", 0))
    timing_rows = []
    for row in bench["rows"]:
        key = "full_process_wall_s" if row["startup"] == "cold" else "complete_run_wall_s"
        old, new = [row["modes"][m][key] for m in ["legacy", "optimized"]]
        timing_rows.append(
            dict(
                B=row["B"],
                startup=row["startup"],
                legacy_s=old,
                optimized_s=new,
                speedup=old / new,
                legacy_samples_per_hour=3600 * row["B"] / old,
                optimized_samples_per_hour=3600 * row["B"] / new,
            )
        )
    top = js(BASE / "validate_v1/F5_TOPOLOGY.json")[0]
    f5 = np.load(BASE / "validate_v1/F5_RELEASE_METRICS.npz")
    write(
        out / "SUMMARY.json",
        dict(
            engineering_status="PASS",
            diverse_full_control_coverage="INCOMPLETE_NUMERICAL_FAILURES",
            identity=identity(),
            timings=timing_rows,
            postprocess_speedup_B64_warm=val["rows"][1]["speedup"],
            comparison_count=len(checks),
            non_bitwise_maxima=maxima,
            families=dict(
                grouping_s=family["grouping_s"],
                static_s=family["static_prepare_s"],
                grouped_s=family["grouped_execution_s"],
                grouped_with_static_s=family["grouped_total_with_static_s"],
                serial_s=family["serial_total_s"],
                resident_speedup=family["actual_grouping_speedup_vs_serial"],
                static_inclusive_speedup=(family["serial_total_s"] + family["static_prepare_s"])
                / family["grouped_total_with_static_s"],
                valid=family["valid_sample_count"],
                samples_per_hour=family["total_valid_samples_per_gpu_hour"],
                independent_initial_families=0,
            ),
            diversity={
                k: diversity[k]
                for k in [
                    "requested_seeds",
                    "prepared_valid",
                    "valid_control_trajectories",
                    "valid_initial_families",
                    "total_generation_s",
                    "unique_initial_families_per_gpu_hour",
                    "scope",
                ]
            },
            complete_three_control_initial_families=0,
            diversity_cold_workflow_s=js(BASE / "campaign_v2/COMPLETE.json")["jobs"][-1][
                "full_process_wall_s"
            ],
            f5=dict(
                frames=len(f5["time_s"]),
                raw_edges=top["raw_edge_changes"],
                topology_counts=top["counts"],
                initial_largest=int(f5["largest_size"][0, 0]),
                at_release_largest=int(f5["largest_size"][600, 0]),
                final_largest=int(f5["largest_size"][-1, 0]),
            ),
            catalog_count=len(catalog),
            catalog_status_counts=dict(Counter(r["numerical_quality"] for r in catalog)),
            physics_scope="EXPLORATORY / PHYSICS_UNCALIBRATED",
            training_ready=False,
        ),
    )
    print("Catalog", len(catalog), "failure batches", len(failures), "non-bitwise maxima", maxima)


if __name__ == "__main__":
    main()
