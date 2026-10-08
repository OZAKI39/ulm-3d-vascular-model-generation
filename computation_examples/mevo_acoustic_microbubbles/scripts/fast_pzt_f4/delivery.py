"""Separate scientific gates and content-addressed F4 archive indices."""

import hashlib
import subprocess
from .common import ROOT, BASE, js, write

SERVER = "/workspace/mevo_fast_pzt_f4/results/fast_pzt/f4"
BRANCH = "dev/fast-pzt-f4-y-geometry-20261008"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    for name in ["ANALYSIS_CHECKS", "ARRAY_VALIDATION"]:
        assert js(BASE / (name + ".json"))["status"] == "PASS"
    assert js(BASE / "visualization/VIDEO_VALIDATION.json")["status"] == "PASS"
    preservation = js(BASE / "PRESERVATION.json")
    assert all(
        not r["missing"] and not r["mismatched"]
        for side in ["local", "server"]
        for r in preservation[side].values()
    )
    quality = js(BASE / "ANALYSIS_CHECKS.json")
    control = js(BASE / "PZT_WALL_CONTROL.json")
    perf = js(BASE / "PERFORMANCE.json")
    gate = dict(
        VESSEL_GEOMETRY=dict(status="PASS_ORIGINAL_MESH", evidence="GEOMETRY_IDENTITY.json"),
        FIELD_GEOMETRY_COMPATIBILITY=dict(
            status="PASS_PARTIAL_OVERLAP_ONLY",
            whole_Y_supported=False,
            evidence="GEOMETRY_FIELD_OVERLAP.json",
        ),
        WALL_CONTACT=dict(
            status="PASS_WITHIN_FIXED_BUDGET",
            evidence="WALL_CONTACT_VALIDATION.json",
            pair_and_wall_budget_m=17.8e-9,
            per_pass_limit_m=178e-9,
        ),
        NEAR_WALL_AGGREGATION=dict(
            status="PASS_BOTH_REGIONS",
            count=32,
            duration_s=0.6,
            wall_contact_observed=True,
            evidence="NEAR_WALL_AGGREGATION.json",
        ),
        PZT_WITH_VESSEL=dict(
            status="WEAK_COMPUTED_RESPONSE",
            robust_shape_control=False,
            additional_wall_contact=False,
            centroid_contrast_m=control["endpoint_contrasts"]["wall_on_switched_vs_wall_on_baseline"][
                "centroid_distance_m"
            ],
            evidence="PZT_WALL_CONTROL.json",
        ),
        PERFORMANCE=dict(
            status="PASS_BOUNDED_MEASURED_OVERHEAD",
            N200_5ms_wall_on_s=perf["N200"]["rows"]["wall_on"]["wall_s"],
            N200_wall_ratio=perf["N200"]["N200_wall_ratio"],
            evidence="PERFORMANCE.json",
        ),
        PHYSICS_SCOPE="EXPLORATORY",
        wall_physics="GEOMETRIC_WALL_ONLY",
        flow_enabled=False,
        numerical_quality=quality,
        no_universal_physical_validation=True,
        limitations=[
            "Stored field covers only a subset of the original Y vessel",
            "Original finite mesh facets and 100 um acoustic grid are preserved",
            "F3 angular-balance, stored-gradient and dense rearrangement sensitivity remain",
            "No wall hydrodynamics, scattering, friction or adhesion",
            "No new F4 repeat/initial-perturbation uncertainty envelope; small PZT contrast is not robust shape control",
        ],
        delivery_source_commit=commit,
        branch=BRANCH,
    )
    write(BASE / "FINAL_GATE.json", gate)
    inputs = {}
    sources = {
        "basis.h5": "/workspace/mevo_acoustic_microbubbles/results/stage3/acoustic/20261002T114813_234502Z_water_2p5MHz/acoustic/basis.h5",
        "response.h5": "/workspace/mevo_acoustic_microbubbles/results/stage4/acoustic/response_extension/20261002T144654140341Z/adaptive_response.h5",
        "labelled_lumen_m.vtp": "/workspace/mevo_acoustic_microbubbles/results/stage4/geometry/20261002T135639758501Z/labelled_lumen_m.vtp",
        "PLANES.json": "/workspace/mevo_acoustic_microbubbles/results/stage4/geometry/20261002T135639758501Z/PLANES.json",
    }
    for p in sorted((BASE / "inputs").iterdir()):
        if p.is_file():
            inputs[p.name] = dict(
                sha256=sha(p),
                bytes=p.stat().st_size,
                local_path=str(p.resolve()),
                server_source=sources.get(p.name, SERVER + "/inputs/" + p.name),
                git_archive=not p.is_symlink(),
                role="Original reused input"
                if p.name in sources
                else "Exact geometry preprocessing / initial state / field-only control preflight",
            )
    write(BASE / "INPUT_MANIFEST.json", dict(inputs=inputs))
    jobs = {}
    for name in ["checks_v2", "edge_checks_v1", "aggregation_v1", "controls_v1", "performance_v1"]:
        jobs[name] = js(BASE / name / "RESULT.json")["identity"]
        for rel, digest in jobs[name]["source_sha256"].items():
            if rel.startswith("src/"):
                assert sha(ROOT / rel) == digest, (name, rel)
    write(
        BASE / "EXECUTION_INDEX.json",
        dict(
            base_f3="3eeb94eb885f5df4f0e214ba9b7904dbdecec744",
            scientific_and_unit_test_jobs=jobs,
            offline_render_identity=js(BASE / "OFFLINE_RENDER_IDENTITY.json"),
            failed_attempts=[
                dict(
                    label="checks_v1",
                    execution_commit="15e8ce4de38d9eeb3432ba8cca4f5ea428adbf08",
                    reason="Pair geometry unit test initialized too close to facet edge; no acoustic trajectory",
                    log_recovery="LOG_RECOVERY.json",
                )
            ],
            source_commit=commit,
            archive_identity="Final branch commit containing this index; archive SHA is distinct from execution SHA",
        ),
    )
    sourcefiles = []
    for dirname in ["src/mevo_mb/fast_sim", "scripts/fast_pzt_f4", "configs/fast_pzt/f4"]:
        sourcefiles += [
            p for p in (ROOT / dirname).rglob("*") if p.is_file() and "__pycache__" not in p.parts
        ]
    sourcefiles += [ROOT / "docs/FAST_PZT_F4_Y_GEOMETRY_DEVELOPMENT.md", ROOT / "docs/FAST_PZT_F4_PLAN.md"]
    write(
        BASE / "DELIVERY_SOURCE.json",
        dict(commit=commit, source_sha256={str(p.relative_to(ROOT)): sha(p) for p in sorted(sourcefiles)}),
    )
    artifacts = []
    for p in sorted(BASE.rglob("*")):
        rel = p.relative_to(BASE)
        if not p.is_file() or p.is_symlink() or p.name == "ARTIFACT_MANIFEST.json":
            continue
        archive = p.stat().st_size <= 4 * 1024**2 or (p.suffix == ".mp4" and p.stat().st_size < 50 * 1024**2)
        # Intermediate preview renders are available locally but not part of final Git evidence.
        if p.name.startswith("preview_"):
            archive = False
        artifacts.append(
            dict(
                relative_path=str(rel),
                sha256=sha(p),
                bytes=p.stat().st_size,
                git_archive=archive,
                local_path=str(p),
                server_path=SERVER + "/" + str(rel),
                role="F4 evidence",
            )
        )
    write(
        BASE / "ARTIFACT_MANIFEST.json",
        dict(
            artifacts=artifacts,
            archive_policy="Files <=4 MiB plus explicitly requested final MP4 up to 50 MiB; reused large inputs indexed only. Full trajectories remain on server/local disk; compact NPZ retains every saved numeric array and event losslessly. Manifest excludes itself.",
            source_commit=commit,
            branch=BRANCH,
        ),
    )
    print("F4 gates and archive index written:", len(artifacts), "files")


if __name__ == "__main__":
    main()
