"""Build separate F3 gates and content-addressed delivery indices from closed evidence."""

import hashlib
import subprocess
from .common import ROOT, BASE, F2, PROTOCOL, write, js

SERVER = "/workspace/mevo_fast_pzt_f3/results/fast_pzt/f3"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    force = js(BASE / "CONTROL_FORCE_DIAGNOSTICS.json")
    repeat = js(BASE / "NUMERICAL_REPRODUCIBILITY.json")
    shape = js(BASE / "MORPHODYNAMICS_COMPARISON.json")
    perf = js(BASE / "PERFORMANCE.json")
    checks = js(BASE / "ANALYSIS_CHECKS.json")
    preservation = js(BASE / "SERVER_PRESERVATION.json")
    assert checks["status"] == "PASS"
    assert all(not v["missing"] and not v["mismatched"] for v in preservation.values())
    assert js(BASE / "visualization/VIDEO_VALIDATION.json")["status"] == "PASS"
    statuses = [v["status"] for v in shape["effects"].values()]
    status = (
        "ROBUST_SHAPE_RESPONSE"
        if "ROBUST_SHAPE_RESPONSE" in statuses
        else (
            "TRANSIENT_SHAPE_RESPONSE"
            if "TRANSIENT_SHAPE_RESPONSE" in statuses
            else "NO_ROBUST_SHAPE_RESPONSE"
        )
    )
    gate = dict(
        schema="fast_pzt_f3_separate_gates_v1",
        PZT_DIFFERENTIAL_DRIVING=dict(
            status=force["differential_driving_status"],
            measurable_nonzero_differential_primary=True,
            control_screening=force["screening_status"],
            meaning="Direct spatial Primary drive is small relative to Secondary; amplitude control changes Total by 4.30% of baseline DS, predominantly by modulation of the existing Secondary pattern. This diagnostic is not a validated shape threshold.",
            force_scale_N=force["force_scale_N"],
            evidence="CONTROL_FORCE_DIAGNOSTICS.json",
        ),
        NUMERICAL_REPRODUCIBILITY=dict(
            status=repeat["status"],
            repeats_per_mode=3,
            stable_positions_bitwise_equal=repeat["modes"]["stable_reduction"]["position_bitwise_equal"],
            stable_sorted_events_identical=repeat["modes"]["stable_reduction"][
                "raw_sorted_edge_events_identical"
            ],
            scope=repeat["repeat_scope"],
            evidence="NUMERICAL_REPRODUCIBILITY.json",
        ),
        PZT_TRANSLATION_CONTROL=dict(
            status="PASS"
            if all(v["translation_resolved"] for v in shape["effects"].values())
            else "NOT_RESOLVED",
            endpoint_distance_m={n: v["translation_endpoint_m"] for n, v in shape["effects"].items()},
            meaning="Control-minus-unchanged centroid contrast; not net displacement from the mature state",
        ),
        PZT_SHAPE_CONTROL=dict(
            status=status,
            endpoint_status={n: v["endpoint_status"] for n, v in shape["effects"].items()},
            longest_eigenvalue_window_s={
                n: v["metrics"]["eigenvalues_relative"]["above_repeat_model_initial_windows"][
                    "longest_sampled_span_s"
                ]
                for n, v in shape["effects"].items()
            },
            sustained_requirement_s=PROTOCOL["observation"]["sustained_window_s"],
            topology=shape["topology_status"],
            controlled_rigid_rotation_established=False,
            sustained_elongation_established=False,
            uncertainty_scope=shape["envelope_scope"],
            evidence="MORPHODYNAMICS_COMPARISON.json",
        ),
        NUMERICAL_QUALITY=dict(
            status="PASS" if all(v["quality"] == "PASS" for v in checks["trajectories"].values()) else "FAIL",
            closed_trajectories=len(checks["trajectories"]),
            max_accepted_overlap_nm=max(
                v["max_accepted_overlap_nm"] for v in checks["trajectories"].values()
            ),
            max_contact_pass_correction_nm=max(
                v["max_contact_pass_correction_nm"] for v in checks["trajectories"].values()
            ),
            limits=PROTOCOL["quality"],
            dt_s=PROTOCOL["dt_s"],
            actual_substeps_per_macro=1,
            meaning="Inherited finite/domain/clock/overlap/per-pass criteria only; does not validate local physical accuracy or dense-cluster morphology",
            dense_contact_correction_over_mobility_path={
                n: v["correction_over_mobility_path"] for n, v in shape["quality"].items()
            },
        ),
        PERFORMANCE=dict(
            status="MEASURED_NO_MATERIAL_REGRESSION",
            N=200,
            physical_duration_s=0.005,
            simulation_wall_s={n: v["wall_s"] for n, v in perf["rows"].items()},
            dataset_wall_s={n: v["dataset_wall_s"] for n, v in perf["rows"].items()},
            stable_to_atomic_wall_ratio=perf["stable_to_atomic_wall_ratio"],
            scope=perf["scope"],
            evidence="PERFORMANCE.json",
        ),
        PHYSICS_SCOPE="EXPLORATORY",
        limitation_flags=[
            "FIELD_CONTROL_LIMITED",
            "FIELD_DIFFERENTIAL_DRIVE_LIMITED",
            "MORPHOLOGY_NUMERICALLY_SENSITIVE",
        ],
        unresolved=[
            "24.2–38.0% stored-gradient versus interpolated-pressure derivative mismatch across eligible controls",
            "Truncated local mobility does not enforce angular balance; instantaneous residual recorded separately",
            "Dense contact corrections comparable to summed mobility-predicted path",
            "Actual 1 pm/1 nm initial perturbations amplify despite same-input determinism",
            "Exploratory pair feedback and near field; no new physical or experimental calibration",
        ],
        source_archive_base=PROTOCOL["f2_archive"],
        delivery_source_commit=commit,
        no_universal_pass=True,
    )
    write(BASE / "FINAL_GATE.json", gate)

    inputs = {}
    old = js(F2 / "INPUT_MANIFEST.json")["inputs"]
    for name in ["basis", "response"]:
        path = BASE / "inputs" / (name + ".h5")
        assert sha(path) == old[name]["sha256"]
        inputs[name] = dict(
            sha256=sha(path),
            bytes=path.stat().st_size,
            local_path=str(path.resolve()),
            server_source=old[name]["server_source"],
            reused=True,
            git_duplicate=False,
        )
    for name, path in {
        "mature_state": F2 / "mature_v1/MATURE_STATE.npz",
        "historical_F1_near_field_contrast": F2 / "mature_v1/F1_near_field/CONTRAST.npz",
        "historical_mobility4_contrast": F2 / "mature_v1/mobility4/CONTRAST.npz",
    }.items():
        inputs[name] = dict(
            sha256=sha(path),
            bytes=path.stat().st_size,
            relative_path=str(path.relative_to(ROOT)),
            local_path=str(path.resolve()),
            server_source="/workspace/mevo_fast_pzt_f2/" + str(path.relative_to(ROOT)),
            reused=True,
            git_duplicate=False,
        )
    write(BASE / "INPUT_MANIFEST.json", dict(inputs=inputs, f2_archive=PROTOCOL["f2_archive"]))

    jobs = {}
    for name in ["checks_v2", "driving_v3", "repeats_v1", "controls_v1", "performance_v1"]:
        row = js(BASE / name / "RESULT.json")
        jobs[name] = row["identity"]
    # Executed core files must still have exactly the bytes recorded by the closed jobs.
    for job in jobs.values():
        for name, digest in job["source_sha256"].items():
            if name.startswith("src/"):
                assert sha(ROOT / name) == digest, name
    write(
        BASE / "EXECUTION_INDEX.json",
        dict(
            f2_archive=PROTOCOL["f2_archive"],
            f2_scientific_executions=dict(
                core="31e66ce73c26e421101c335009f20a9462179a40",
                sensitivity="cd7d6f5ca4c8a3298b5e7aedf2f64642f9c85316",
                mature="bfd228ab69fe7fe7a3ba29c43b06cf53d0f4f3c2",
                performance="8feabe08c664a190cdea719dbd43aee96e3a4856",
                phase="03655ec93f2cb518a428d9db496a5312258d1575",
            ),
            f3_scientific_jobs=jobs,
            delivery_source_commit=commit,
            offline_render_identity="OFFLINE_RENDER_IDENTITY.json",
            archive_identity="The final branch commit containing this index; see git log. Archive SHA is not substituted for scientific execution SHAs.",
            failed_attempts=[
                dict(
                    label="driving_v1",
                    execution_commit="facc9ca8d816628b312177a74d10f4aa334316a2",
                    reason="Missing write re-export; no simulation",
                ),
                dict(
                    label="checks_v1",
                    execution_commit="2626672f787e312c2d6ee99de2bdf3276091475b",
                    reason="Missing write re-export; no simulation",
                ),
                dict(
                    label="driving_v2",
                    execution_commit="8cf42b5817f9f545cc6016b2546a59e8f122e90f",
                    reason="Prescribed Z-pair ROI-RMS equation infeasible; no trajectory; slot subsequently marked ineligible",
                ),
            ],
        ),
    )
    write(
        BASE / "DELIVERY_SOURCE.json",
        dict(
            commit=commit,
            sha256={
                str(p.relative_to(ROOT)): sha(p)
                for folder in [ROOT / "scripts/fast_pzt_f3", ROOT / "configs/fast_pzt/f3"]
                for p in sorted(folder.rglob("*"))
                if p.is_file() and "__pycache__" not in p.parts
            },
        ),
    )
    artifacts = []
    for path in sorted(BASE.rglob("*")):
        rel = path.relative_to(BASE)
        if (
            not path.is_file()
            or path.is_symlink()
            or rel.parts[0] == "inputs"
            or path.name == "ARTIFACT_MANIFEST.json"
        ):
            continue
        artifacts.append(
            dict(
                relative_path=str(rel),
                bytes=path.stat().st_size,
                sha256=sha(path),
                git_archive=path.stat().st_size <= 4 * 1024**2,
                local_path=str(path),
                server_path=SERVER + "/" + str(rel),
                role="F3_evidence",
            )
        )
    write(
        BASE / "ARTIFACT_MANIFEST.json",
        dict(
            schema="fast_pzt_f3_artifacts_v1",
            artifacts=artifacts,
            note="Closed F3 evidence only; immutable inherited inputs are indexed separately; manifest excludes itself",
        ),
    )
    print(
        {
            k: gate[k]["status"]
            for k in [
                "PZT_DIFFERENTIAL_DRIVING",
                "NUMERICAL_REPRODUCIBILITY",
                "PZT_TRANSLATION_CONTROL",
                "PZT_SHAPE_CONTROL",
                "NUMERICAL_QUALITY",
                "PERFORMANCE",
            ]
        }
    )
    print("F3 indexed files:", len(artifacts), "bytes:", sum(a["bytes"] for a in artifacts))


if __name__ == "__main__":
    main()
