"""Assemble capability-specific gates and a hash index after all runs/rendering close."""

import hashlib
import json
import subprocess
from .common import BASE, ROOT, write


def get(name):
    return json.loads((BASE / name).read_text())


def main():
    identity = get("balanced_morphodynamics_v1/controlled/IDENTITY.json")
    for key in ["response", "basis"]:
        assert hashlib.sha256((BASE / "inputs" / (key + ".h5")).read_bytes()).hexdigest() == identity[key]["sha256"]
    checks = get("balanced_core_v1/GPU_CHECKS.json")
    analysis = get("ANALYSIS_CHECKS.json")
    source = get("FINAL_SOURCE_CHECK.json")
    performance = get("balanced_performance_v1/RESULT.json")
    throughput = get("balanced_throughput_v1/RESULT.json")
    morph = get("balanced_morphodynamics_v1/RESULT.json")
    video = get("visualization/VIDEO_VALIDATION.json")
    step = get("balanced_steps_v1/RESULT.json")
    selected = next(r for r in step["rows"] if r["name"] == step["selected"])
    final_runs = {}
    for path in sorted(BASE.glob("balanced*/**/SUMMARY.json")):
        summary = json.loads(path.read_text())
        final_runs[str(path.parent.relative_to(BASE))] = dict(
            status=summary["status"],
            numerical_quality=summary["numerical_quality"],
            sample_quality=summary["sample_numerical_quality"],
            max_overlap_m=summary["quality"]["max_accepted_penetration_m"],
            physical_duration_s=summary["completed_physical_s"],
            macro_step_count=summary["macro_step_count"],
            actual_substep_count=summary["actual_substep_count"],
        )
    core_pass = checks["status"] == "PASS" and source["status"] == "PASS"
    numerical_pass = all(r["numerical_quality"] == "PASS" for r in final_runs.values())
    aggregate_pass = (
        analysis["matched_OFF_max_displacement_m"] == 0
        and selected["initial_largest_size"] == [1]
        and selected["final_largest_size"] == [32]
    )
    short_pass = performance["short_target_10s"] == "PASS"
    long_pass = (
        throughput["long_run"]["status"] == "COMPLETE"
        and throughput["long_run"]["numerical_quality"] == "PASS"
        and throughput["long_run"]["process_body_wall_s"] < 60
    )
    batch_pass = (
        throughput["batch_certified_valid_samples"] == 8
        and throughput["serial_valid_samples"] == 8
        and throughput["raw_sample_throughput_speedup"] > 1
    )
    assert video["status"] == "PASS" and analysis["status"] == "PASS"
    for stem in [
        "Fig1_aggregation_formation",
        "Fig2_Secondary_ON_OFF",
        "Fig3_PZT_morphodynamics",
        "Fig4_step_performance",
    ]:
        for extension in ["png", "pdf"]:
            assert (BASE / f"figures/{stem}.{extension}").is_file()
    gate = dict(
        schema="fast_pzt_f1_capabilities_v1",
        model="fast_pzt_f1_v2_balanced_mobility",
        branch="dev/fast-pzt-aggregation-f1-20261007",
        reviewed_base="5c0a0f10b19b376d9321510149cb6db50bcf4331",
        scientific_execution_commit="9829156893339efff713483f3125848b84adbc98",
        matched_controls_execution_commit="3721bebd597164262691af766f723d245e043116",
        FAST_GPU_CORE="PASS" if core_pass else "FAIL",
        BASIC_AGGREGATION="PASS" if aggregate_pass else "FAIL",
        PZT_MORPHODYNAMICS=morph["pzt_morphodynamics"],
        NUMERICAL_QUALITY="PASS" if numerical_pass else "FAIL",
        THROUGHPUT_TARGET="PASS" if short_pass and long_pass and batch_pass else "FAIL",
        PHYSICS_SCOPE="EXPLORATORY",
        quantitative_physics_calibration=False,
        selected_macro_step_s=selected["nominal_dt_s"],
        physical_substeps_per_macro=1,
        local_contact_iterations=12,
        mobility_iterations=2,
        effective_steps_per_physical_second=50000,
        small_GPU_checks=len(checks["checks"]),
        quality_protocol="configs/fast_pzt/f1_protocol.json (unchanged)",
        selected_step_statistics=selected["comparison"],
        selected_runs=final_runs,
        pzt_same_state_comparison=morph["comparison"],
        throughput=dict(
            N200_5ms_status="PASS" if short_pass else "FAIL",
            N200_5ms_simulation_wall_s=performance["short_run"]["wall_s"],
            N200_5ms_process_body_wall_s=performance["short_run"]["process_body_wall_s"],
            N32_1s_status="PASS" if long_pass else "FAIL",
            N32_1s_simulation_wall_s=throughput["long_run"]["wall_s"],
            N32_1s_process_body_wall_s=throughput["long_run"]["process_body_wall_s"],
            N200_1s_status="NOT_RUN",
            batch_status="PASS" if batch_pass else "FAIL",
            valid_batch_samples=throughput["batch_certified_valid_samples"],
            valid_serial_samples=throughput["serial_valid_samples"],
            valid_samples_per_GPU_hour=throughput["batch_valid_samples_per_GPU_hour"],
            throughput_speedup=throughput["raw_sample_throughput_speedup"],
            throughput_scope="8 independent 32-bubble / 50-ms samples; includes sampled metrics/HDF5, excludes one-time startup/rendering/archive transfer",
        ),
        figures="figures/",
        movie=video,
        superseded_evidence="PROTOTYPE_RETIREMENT.json",
        limitations=[
            "Uncalibrated smooth near-field force and local hydrodynamics",
            "Independent radial table; no collective acoustic/radial scattering",
            "Packed contact correction is substantial and can influence morphology",
            "Graph splits/merges are threshold-defined topology, not bubble breakup/fusion",
            "Atomic accumulation is not bitwise deterministic; no individual-trajectory accuracy claim",
            "No walls, flow or polydisperse batch implemented",
        ],
    )
    write(BASE / "FINAL_GATE.json", gate)
    # This index excludes itself to avoid circular hashes. Both roots contain the same artifacts after sync.
    rows = []
    for path in sorted(BASE.rglob("*")):
        if not path.is_file() or path.name == "ARTIFACT_MANIFEST.json":
            continue
        rel = path.relative_to(BASE)
        size = path.stat().st_size
        keep = size <= 4 * 1024**2
        if path.suffix == ".h5":
            keep = str(rel) == "inputs/response.h5" or (
                str(rel).startswith(
                    ("balanced_steps_v1/", "balanced_controls_v1/", "balanced_morphodynamics_v1/")
                )
                and size <= 4 * 1024**2
            )
        rows.append(
            dict(
                relative_path=str(rel),
                bytes=size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                local_path=str(path),
                server_path="/workspace/mevo_fast_pzt_f1/results/fast_pzt/f1/" + str(rel),
                git_archive=keep,
                role="retired_diagnostic"
                if str(rel)
                .split("/")[0]
                .startswith(
                    (
                        "smoke",
                        "basic_v1",
                        "time_step_v1",
                        "contact8_v1",
                        "contact12_v1",
                        "short_performance",
                        "morphodynamics_v1",
                        "core_layout_v2",
                    )
                )
                else "final_evidence_or_input",
            )
        )
    write(
        BASE / "ARTIFACT_MANIFEST.json",
        dict(
            schema="fast_pzt_f1_artifacts_v1",
            artifacts=rows,
            archive_policy="Small metadata/logs/figures and necessary scene HDF5 <=4 MiB; other raw HDF5 and files >4 MiB remain in both indexed roots",
            reproduction="scripts/fast_pzt/README.md",
            table_basis_identity="balanced_morphodynamics_v1/controlled/IDENTITY.json",
            source_tree_commit_at_assembly=subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
            ).strip(),
            physics_scope="EXPLORATORY",
        ),
    )
    print(
        json.dumps(
            {
                k: gate[k]
                for k in [
                    "FAST_GPU_CORE",
                    "BASIC_AGGREGATION",
                    "PZT_MORPHODYNAMICS",
                    "NUMERICAL_QUALITY",
                    "THROUGHPUT_TARGET",
                    "PHYSICS_SCOPE",
                ]
            }
        )
    )


if __name__ == "__main__":
    main()
