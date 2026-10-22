"""Capability-specific gates and small-file archive index; run after artifacts close."""

import hashlib
import subprocess
from .common import BASE, ROOT, write
from .analyze import js, brief, BRANCHES


def main():
    identity = js(BRANCHES["phase"] / "IDENTITY.json")
    inputs = {}
    for name in ["response", "basis"]:
        path = BASE / "inputs" / (name + ".h5")
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        assert sha == identity[name]["sha256"]
        inputs[name] = dict(
            sha256=sha,
            bytes=path.stat().st_size,
            local_path=str(path.resolve()),
            server_source=identity[name]["path"],
            metadata=identity[name],
        )
    write(
        BASE / "INPUT_MANIFEST.json",
        dict(
            inputs=inputs,
            policy="Reuse F1 response and basis unchanged; response already archived in F1; large basis retained externally by hash",
        ),
    )
    unchanged = {}
    for name in ["nearfield.py", "contact_models.py", "f2_engine.py", "f2_kernels.py", "morphodynamics.py"]:
        rel = "src/mevo_mb/fast_sim/" + name
        sha = hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()
        assert sha == identity["source_sha256"][rel]
        unchanged[rel] = sha
    write(
        BASE / "FINAL_SOURCE_CHECK.json",
        dict(
            status="PASS",
            compared_execution=identity["execution_commit"],
            unchanged_scientific_modules=unchanged,
            F1_preservation="F1_PRESERVATION.json",
            note="Offline analysis/report scripts are newer; executed source hashes and commit identities remain immutable",
        ),
    )
    checks = js(BASE / "core_v1/GPU_CHECKS.json")
    assert checks["status"] == "PASS" and js(BASE / "ANALYSIS_CHECKS.json")["status"] == "PASS"
    pzt = js(BASE / "PZT_MORPHODYNAMICS.json")
    perf = js(BASE / "PERFORMANCE.json")
    video = js(BASE / "visualization/VIDEO_VALIDATION.json")
    assert video["status"] == "PASS"
    paths = list(BRANCHES.values()) + [
        BASE / "mature_v1/formation",
        BASE / "sensitivity_v1/near_10us",
        BASE / "performance_v1/N200_5ms",
    ]
    selected = {str(p.relative_to(BASE)): brief(p) for p in paths}
    assert all(r["quality"] == "PASS" for r in selected.values())
    failed = {
        str(p.parent.relative_to(BASE)): brief(p.parent)
        for p in BASE.glob("sensitivity_v1/*/SUMMARY.json")
        if js(p)["numerical_quality"] != "PASS"
    }
    for name in [
        "Fig_F2_1_near_field",
        "Fig_F2_2_contact_sensitivity",
        "Fig_F2_3_PZT_morphodynamics",
        "Fig_F2_4_3D_cluster_shapes",
    ]:
        for ext in ["png", "pdf"]:
            assert (BASE / "figures" / f"{name}.{ext}").is_file()
    executions = {
        str(p.parent.relative_to(BASE)): js(p)["execution_commit"]
        for p in sorted(BASE.rglob("IDENTITY.json"))
    }
    gate = dict(
        schema="fast_pzt_f2_capabilities_v1",
        branch="dev/fast-pzt-aggregation-f2-20261008",
        NEAR_FIELD_IMPLEMENTATION="PASS",
        CONTACT_SENSITIVITY="CONTACT_MODEL_SENSITIVE",
        contact_candidate_disposition="REJECTED_AS_DEFAULT_ANGULAR_BIAS",
        PZT_MORPHODYNAMICS=pzt["status"],
        pzt_status_scope="120 ms endpoint; transient eigenvalue responses are separately observed, not hidden by this endpoint label",
        pzt_experiment_execution="PASS_SAME_MATURE_STATE",
        transient_shape_response="DESCRIPTIVE_WINDOWS_OBSERVED_CONTACT_ROBUSTNESS_NOT_ESTABLISHED",
        NUMERICAL_QUALITY="PASS_SELECTED_BACKEND_ONLY",
        PERFORMANCE="SAME_ORDER_2_060X_F1_WALL_COST",
        PHYSICS_SCOPE="EXPLORATORY",
        calibrated=False,
        selected_backend=js(BASE / "adoption_v1/SELECTED.json"),
        dt_s=1e-5,
        physical_substeps_per_macro=1,
        mobility_iterations=2,
        contact_iterations=12,
        quality_limits=dict(max_accepted_overlap_nm=17.8, max_contact_pass_correction_nm=178),
        selected_run_quality=selected,
        failed_candidate_runs_preserved=failed,
        N200_5ms_simulation_wall_s=perf["run"]["wall_s"],
        N200_5ms_F1_ratio=perf["wall_ratio_to_F1"],
        GPU_batching="PRESERVED_TINY_INDEPENDENCE_CHECK_PASSED_NO_NEW_BATCH8_BENCHMARK",
        execution_commits=executions,
        F1_preservation="F1_PRESERVATION.json",
        source_check="FINAL_SOURCE_CHECK.json",
        video=video,
        limitations=[
            "Pairwise fundamental closure is not a collective radial/multipole theory; Primary remains uncoupled",
            "Contact candidate has spurious angular corrections; F1 central PBD remains selected",
            "Packing and ID-level rearrangements depend on force/contact/mobility and atomic order",
            "Single seed and one repeat provide descriptive sensitivity, not confidence intervals or calibrated prediction",
            "No lasting controlled elongation, isolated rigid rotation, member exchange or persistent component split/merge established",
            "No flow, walls, polydispersity, arbitrary frequency/size or physical bubble fusion/breakup",
        ],
    )
    write(BASE / "FINAL_GATE.json", gate)
    artifacts = []
    for p in sorted(BASE.rglob("*")):
        if not p.is_file() or p.name == "ARTIFACT_MANIFEST.json":
            continue
        rel = str(p.relative_to(BASE))
        size = p.stat().st_size
        keep = size <= 4 * 1024**2 and not rel.startswith("inputs/")
        artifacts.append(
            dict(
                relative_path=rel,
                bytes=size,
                sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                git_archive=keep,
                local_path=str(p),
                server_path="/workspace/mevo_fast_pzt_f2/results/fast_pzt/f2/" + rel,
                role="external_unchanged_F1_input"
                if rel.startswith("inputs/")
                else "diagnostic_or_failed_candidate"
                if rel.startswith(("sensitivity_v1/contact_", "sensitivity_v1/near_20us"))
                else "F2_evidence",
            )
        )
    write(
        BASE / "ARTIFACT_MANIFEST.json",
        dict(
            schema="fast_pzt_f2_artifacts_v1",
            artifacts=artifacts,
            archive_policy="Necessary raw coordinates, metrics, configs, logs, figures and movie <=4 MiB; larger artifacts and unchanged input files indexed in both roots. Manifest excludes itself.",
            source_tree_commit_at_assembly=subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
            ).strip(),
            execution_identity="Per-run IDENTITY.json; final archive commit is distinct and available from Git",
            reproduction="scripts/fast_pzt_f2/README.md",
        ),
    )
    print(
        {
            k: gate[k]
            for k in [
                "NEAR_FIELD_IMPLEMENTATION",
                "CONTACT_SENSITIVITY",
                "PZT_MORPHODYNAMICS",
                "NUMERICAL_QUALITY",
                "PERFORMANCE",
                "PHYSICS_SCOPE",
            ]
        }
    )
    print(
        "Archive files",
        sum(r["git_archive"] for r in artifacts),
        "bytes",
        sum(r["bytes"] for r in artifacts if r["git_archive"]),
    )


if __name__ == "__main__":
    main()
