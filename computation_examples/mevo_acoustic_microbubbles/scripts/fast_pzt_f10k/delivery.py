"""Consolidate measured evidence without running or modifying scientific trajectories."""

import h5py
import numpy as np
from .common import ROOT, BASE, BASELINE, identity, js, sha, write


def main():
    analysis = BASE / "analysis"
    topology = js(analysis / "B2_TOPOLOGY_COMPARISON.json")
    comparison = topology["comparison"]
    scaling = js(analysis / "SCALING.json")
    navigation = js(analysis / "NAVIGATION.json")
    selected = js(analysis / "ALTERNATIVE_SCREEN.json")["selected"]
    boundary = []
    for path in [
        ROOT / selected["source"],
        *[ROOT / a["path"] / "trajectory.h5" for a in topology["audits"]],
    ]:
        with h5py.File(path) as h:
            labels = h["boundary/label"][:]
            times = h["boundary/time_s"][:]
            active = h["trajectory/active"][:]
            row = dict(
                source=str(path.relative_to(ROOT)),
                sha256=sha(path),
                boundary_labels=labels.tolist(),
                boundary_times_s=times.tolist(),
                nonzero_boundary_labels=int(np.count_nonzero(labels)),
                all_saved_IDs_active=bool(active.all()),
                minimum_saved_wall_gap_m=float(h["trajectory/wall_gap_lower_bound_m"][:].min()),
                wall_stat_max_by_column=h["trajectory/wall_stats"][:].max(axis=(0, 1, 2)).tolist(),
            )
            assert row["nonzero_boundary_labels"] == 0 and row["all_saved_IDs_active"]
            boundary.append(row)
    write(
        analysis / "BOUNDARY_AUDIT.json",
        dict(
            status="PASS",
            rows=boundary,
            note="Zero boundary labels mean no stored exit events; time=0 is inactive sentinel, not an exit at zero. Includes entire selected historical file and both new windows.",
        ),
    )

    baseline_scope = js(analysis / "SCIENTIFIC_SCOPE.json")
    new_rows = [
        dict(
            model="B2",
            task="One interior persistent graph merge",
            status="PERSISTENT_MERGE_CHECKED_ONE_CASE",
            evidence="Same IDs and 0.290 ms event at dt10/5 us, observed for2.710 ms; original2 ms persistence criterion.",
            limitation="Cross-model initialization, one 31+1 case; not physical bubble fusion or arbitrary topology validation.",
            source="B2_TOPOLOGY_COMPARISON.json",
        ),
        dict(
            model="B2",
            task="Persistent split",
            status="PERSISTENT_TOPOLOGY_UNTESTED",
            evidence="No split in either new3 ms window; near-wall F5 candidate not admitted.",
            limitation="Missing positive persistent splitting and independent physical calibration.",
            source="F5_ADMISSION.json",
        ),
        dict(
            model="B2",
            task="Internal morphology and instantaneous edges",
            status="MORPHOLOGY_LIMITED",
            evidence="Positions, Rg, aspect, tensor and connected members pass existing gates; raw46/47 edges and up to60 us fail.",
            limitation="Old12/29 strict timing FAIL remains; no threshold relaxation.",
            source="B2_TOPOLOGY_COMPARISON.json",
        ),
        dict(
            model="M0",
            task="Actual interior downstream-plane control",
            status="GROUP_MOTION_CHECKED",
            evidence="Actual O trajectories cross20 um at2.02433 ms; O-OFF20.32/20.53 nm, below original50 nm useful-effect gate.",
            limitation="Navigation to prospective+1 um Z target is not established; 10 ms field support is right-censored.",
            source="NAVIGATION.json",
        ),
        dict(
            model="B2",
            task="B=1 fixed-state size scaling",
            status="NUMERICAL_CHECKED",
            evidence="32 real and128/512 synthetic pass both solvers; GPU advantage only measured at128/512;1024 NOT_RUN.",
            limitation="Regular synthetic networks, no long trajectories, no B64 production or universal crossover.",
            source="SCALING.json",
        ),
    ]
    for row in new_rows:
        row.update(
            source_sha256=sha(analysis / row["source"]),
            training_ready=False,
            physical_calibration="PHYSICS_UNCALIBRATED",
        )
    write(
        analysis / "SCIENTIFIC_SCOPE_FINAL.json",
        dict(
            baseline_through="F10J",
            baseline_rows=baseline_scope["rows"],
            baseline_scope_sha256=sha(analysis / "SCIENTIFIC_SCOPE.json"),
            F10K_additions=new_rows,
            labels=["EXPLORATORY", "PHYSICS_UNCALIBRATED"],
            training_ready=False,
            interpretation="New rows add limited evidence; historical model results remain descriptive, not cross-model truth.",
        ),
    )

    jobs = {
        p.stem.removesuffix("_PROCESS"): dict(
            full_process_wall_s=js(p)["full_process_wall_s"], exit_code=js(p)["exit_code"], sha256=sha(p)
        )
        for p in sorted(BASE.glob("*_PROCESS.json"))
    }
    failures = dict(
        topology_dt10_v1="300 scientific steps and durable files completed, then metadata KeyError protocol_sha256. Original exit1 and files retained; recovery writes sidecars only, no science rerun.",
        scaling_v1="Missing threadpoolctl at import, no numerical solve.",
        scaling_v2="Missing threadpoolctl at import, no numerical solve.",
        package_install="Environment has no pip; retained failed pip log, installed only threadpoolctl==3.6.0 with uv --no-deps and recorded official wheel source.",
        navigation_v1="Offline analysis float equality selection failed. Integer-derived dt matching fixed; original trajectories unchanged; subsequent analysis passed.",
    )
    write(
        analysis / "COST_AND_FAILURES.json",
        dict(
            jobs=jobs,
            full_supervised_process_wall_sum_s=sum(p["full_process_wall_s"] for p in jobs.values()),
            topology_integration_quality_sum_s=sum(
                a["timing"]["integration_quality_wall_s"] for a in topology["audits"]
            ),
            scaling_solve_warmup_s=scaling["solve_and_warmup_s"],
            scaling_stage_s=scaling["stage_wall_s"],
            startup_attempts_s=sum(scaling["failed_startup_process_s"].values()),
            sampled_scaling_GPU_process_peak_MiB=scaling["sampled_GPU_process_peak_MiB"],
            scaling_peak_RSS_bytes=max(
                max(r.get("peak_RSS_bytes", 0), r.get("controller_sampled_peak_RSS_bytes", 0))
                for r in scaling["rows"]
            ),
            failures=failures,
            scope="Process sum is not total project/human elapsed time; local read-only analysis, plotting, transfer and dependency installation are separate.",
        ),
    )

    write(
        analysis / "FINAL_GATE.json",
        dict(
            delivery_status="COMPLETED_EXPLORATORY",
            baseline_commit=BASELINE,
            engineering="PASS_WITH_RETAINED_FAILURES_AND_VERIFIED_RECOVERY",
            F5_admission="B2_TOPOLOGY_CASE_NOT_ADMITTED",
            alternative_initialization="ADMITTED_INTERIOR_CROSS_MODEL_COUNTERFACTUAL",
            numerical_contact_quality="PASS_BOTH_SHORT_WINDOWS",
            group_motion="CHECKED_LIMITED_SCENARIOS",
            individual_shape_gates={k: v for k, v in comparison["gates"].items() if k != "edge_events"},
            morphology=comparison["status"],
            instantaneous_edge_gate="FAIL",
            historical_strict_edge_failure="UNCHANGED_12_OF_29_FAIL_MAX_130_US",
            persistent_merge="ONE_POSITIVE_SAME_MODEL_DT_PAIR",
            persistent_split="PERSISTENT_TOPOLOGY_UNTESTED",
            navigation="NOT_ESTABLISHED",
            original50nm_effect_gate_at_plane=navigation["original50nm_effect_gate_at_plane"],
            fixed_state_scaling={str(r["N"]): r["status"] for r in scaling["rows"]},
            largest_measured_N=scaling["largest_measured_N"],
            general_GPU_crossover="NOT_ESTABLISHED_REGULAR_SYNTHETIC_CASES_ONLY",
            near_wall="NEAR_WALL_UNCALIBRATED",
            physics="PHYSICS_UNCALIBRATED",
            labels=["EXPLORATORY", "PHYSICS_UNCALIBRATED"],
            training_ready=False,
            scope="Graph merging is not bubble fusion. CPU/GPU correctness does not calibrate contact/near-field/wall physics. No new production admission.",
            analysis_sources={
                str(p.relative_to(ROOT)): sha(p)
                for p in sorted(analysis.glob("*.json"))
                if p.name != "FINAL_GATE.json"
            },
            identity=identity(),
        ),
    )
    print("Layered final gate written; morphology/navigation/AI limitations retained.")


if __name__ == "__main__":
    main()
