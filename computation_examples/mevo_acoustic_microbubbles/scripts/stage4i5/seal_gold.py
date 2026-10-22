"""Seal the two naturally completed Stage 4I cases without rerunning them."""
import argparse
import hashlib
import importlib.metadata
import json
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BATCH = "results/stage4/primary_secondary/20261003T132343245938Z"
DEST = "results/stage4/performance_refactor/gold_baseline"


def read(path):
    return json.loads(Path(path).read_text())


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def validate_case(row):
    folder = ROOT / row["run"]
    summary = read(folder / "particles/SUMMARY.json")
    identity = read(folder / "RUN_IDENTITY.json")
    sources = read(folder / "SOURCE_IDENTITY.json")
    branch = read(folder / "BRANCH_METRICS.json")
    assert summary["status"] == "COMPLETE" and summary["duration_s"] == 6.0
    assert not (folder / "EXECUTION_FAILURE.json").exists()
    assert not list(folder.rglob("*.tmp"))
    with zipfile.ZipFile(folder / "logs/source_snapshot.zip") as archive:
        for name, expected in sources["source_sha256"].items():
            assert hashlib.sha256(archive.read(name)).hexdigest() == expected, name
    with np.load(folder / "INITIAL_POPULATION.npz") as initial:
        initial_hash = hashlib.sha256(initial["positions_m"].tobytes()
            + initial["radii_m"].tobytes() + initial["particle_id"].tobytes()).hexdigest()
        assert initial_hash == identity["initial_population_sha256"]
        initial_positions = initial["positions_m"].copy()
        ids = initial["particle_id"].copy()
    assert np.array_equal(ids, np.arange(256))
    assert sum(branch[k] for k in ["N_A", "N_B", "N_unassigned", "N_wall_loss"]) == len(ids)
    crossings = branch["crossings"]
    assert len({e["particle_id"] for e in crossings}) == len(crossings)
    with h5py.File(folder / "particles/trajectories.h5", "r") as h:
        assert bool(h.attrs["complete"])
        times, active = h["time_s"][:], h["active"][:]
        positions = h["positions_m"][:]
        assert times[0] == 0.0 and times[-1] == 6.0 and np.all(np.diff(times) > 0)
        assert np.array_equal(h["particle_id"][:], ids)
        assert np.array_equal(positions[0], initial_positions)
        assert np.isfinite(positions).all()
        assert not np.any(active[1:] & ~active[:-1])
        assert int(active[-1].sum()) == branch["N_unassigned"]
        finite_fields = []
        for name in h:
            ds = h[name]
            if ds.shape[:2] == active.shape and ds.dtype.kind in "fc":
                assert np.isfinite(ds[:][active]).all(), name
                finite_fields.append(name)
        for event in crossings:
            pid = event["particle_id"]
            assert branch["assignment"][str(pid)] == event["branch"]
            assert not active[-1, pid]
            assert np.array_equal(positions[-1, pid], event["position_m"])
            assert 0 <= event["crossing_time_s"] <= 6
        saved_dt = h["accepted_dt_s"][:]
        slip = h["total_velocity_m_s"][:] - h["flow_velocity_m_s"][:]
        medium = identity["config"]["medium"]
        reynolds = (2 * medium["density_kg_m3"] * h["radii_m"][:][None, :]
                    * np.linalg.norm(slip, axis=-1) / medium["viscosity_Pa_s"])
        frame, particle = np.unravel_index(np.where(active, reynolds, -np.inf).argmax(), active.shape)
        trajectory = dict(frames=len(times),initial_population_sha256=initial_hash,
            finite_active_fields=finite_fields,all_positions_finite=True,
            complete=True,time_start_s=float(times[0]),time_end_s=float(times[-1]),
            saved_positive_dt_min_s=float(saved_dt[saved_dt > 0].min()),
            saved_dt_max_s=float(saved_dt.max()),
            saved_constrained_slip_Re=dict(maximum=float(reynolds[frame,particle]),
                time_s=float(times[frame]),particle_id=int(ids[particle]),
                active_particle_time_samples=int(active.sum()),
                samples_Re_at_least_1=int(((reynolds >= 1) & active).sum()),
                scope="Saved active states, not all accepted steps; inactive NaNs are intentional."))
    return summary, identity, sources, trajectory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--boundary-evidence", type=Path, required=True)
    args = parser.parse_args()
    batch = read(ROOT / BATCH / "SUMMARY.json")
    state = read(ROOT / "results/stage4/primary_secondary/STATE.json")
    assert batch == state and state["status"] == "PASS_WITH_SCOPE"
    assert state["run"] == BATCH and {r["case"] for r in state["rows"]} == {"B2", "B3"}
    evidence = read(args.boundary_evidence)
    assert evidence["boundary_finalized"]["orchestrator_exists"] is False
    assert evidence["boundary_finalized"]["convergence_batches"] == []
    assert not list((ROOT / "results/stage4/convergence").glob("*/IDENTITY.json"))
    verified = {r["case"]:validate_case(r) for r in state["rows"]}
    out = ROOT / DEST
    out.mkdir(parents=True, exist_ok=False)
    write(out / "BOUNDARY_OPERATIONS.json", evidence)
    records = {}

    def preserve(source, relative):
        source = Path(source)
        target = out / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        before = source.stat()
        expected = digest(source)
        shutil.copy2(source, target)
        assert digest(target) == expected
        after = source.stat()
        assert (before.st_size,before.st_mtime_ns,before.st_ctime_ns) == (after.st_size,after.st_mtime_ns,after.st_ctime_ns)
        records[str(relative)] = dict(source=str(source.relative_to(ROOT)),bytes=before.st_size,sha256=expected)
        print(json.dumps(dict(preserved=str(relative),bytes=before.st_size,sha256=expected)),flush=True)

    cases = {}
    inputs = {}
    for row in state["rows"]:
        label = row["case"]
        summary, identity, sources, trajectory = verified[label]
        folder = ROOT / row["run"]
        for source in sorted(folder.rglob("*")):
            if source.is_file():
                preserve(source, Path("cases") / label / source.relative_to(folder))
        config = identity["config"]
        for name, path, expected in [
            ("wall_surface",config["geometry"]["mesh_surface"],None),
            ("flow",config["flow"]["snapshot"],identity["flow_sha256"]),
            ("acoustic_field",identity["acoustic_field_path"],identity["acoustic_field_sha256"]),
            ("response_table",identity["response_table_path"],identity["response_table_sha256"])]:
            p = Path(path) if Path(path).is_absolute() else ROOT / path
            inputs[str(p.relative_to(ROOT))] = dict(kind=name,expected_sha256=expected)
        process = evidence["boundary_hold"]["worker_final_samples"][label]
        cases[label] = dict(result_id=row["run"],git_commit_at_run=sources["git_commit"],
            source_sha256=sources["source_sha256"],config_sha256=digest(folder/"config.yaml"),
            identity=identity,branch_counts={k:row[k] for k in ["N_A","N_B","N_unassigned","N_wall_loss","f_A","f_B"]},
            initial_particle_file=f"cases/{label}/INITIAL_POPULATION.npz",
            branch_assignments_and_crossings=f"cases/{label}/BRANCH_METRICS.json",
            trajectory_file=f"cases/{label}/particles/trajectories.h5",trajectory_checks=trajectory,
            accepted_steps=summary["accepted_steps"],nominal_macro_steps=12000,
            additional_accepted_substeps=summary["accepted_steps"]-12000,
            rejected_steps=dict(value=None,status="NOT_RECORDED",scope="No error-estimator rejection/retry loop in this implementation; failures would abort, not silently retry."),
            effective_dt=dict(minimum_s=summary["minimum_dt_s"],mean_s=6/summary["accepted_steps"],
                maximum_s=None,maximum_status="NOT_RECORDED_FOR_ALL_STEPS",upper_bound_s=.0005,
                maximum_saved_s=trajectory["saved_dt_max_s"],scope="Only 301 visualization frames stored. Small macro-boundary remainder steps can determine the recorded minimum."),
            contact_pair_events=summary["contact_pair_events"],wall_hit_events=summary["wall_hit_events"],
            wall_hit_particle_count=summary["wall_hit_particle_count"],
            wall_projection_call_count=dict(value=None,status="NOT_RECORDED"),
            geometry_query_counters=summary["geometry_cache"],
            pair_force_calls=dict(value=summary["physics_validity"]["checks"],status="DERIVED_FROM_CODE",scope="Each successful evaluate call records one validity check and one secondary-force call; not a separately instrumented counter."),
            validity=summary["physics_validity"],runtime_s=summary["elapsed_s"],timings=summary["timings"],
            worker_process_cpu=dict(user_s=process["user_cpu_s"],system_s=process["system_cpu_s"],sampled_at=process["time"],scope="Process-lifetime /proc sample at completed-case boundary; includes initialization and negligible pool cleanup."),
            peak_RAM=dict(VmHWM_kB=process["peak_rss_kB"],sampled_at=process["time"],scope="Linux process high-water mark observed after scientific completion, before natural worker exit."),
            gpu_time=dict(value=None,status="NOT_INSTRUMENTED"),
            sampled_VRAM=dict(evidence="BOUNDARY_OPERATIONS.json",scope="Device-wide nvidia-smi samples; not a per-case lifetime maximum or a GPU-time measurement."))
    for kind, key, identity_key in [("acoustic","basis","acoustic_basis_sha256"),("simvascular/mesh","mesh","mesh_sha256")]:
        saved = read(ROOT / "results/stage4" / kind / "STATE.json")
        inputs[saved[key]] = dict(kind=kind,expected_sha256=state["identity"][identity_key])
    geometry = read(ROOT / "results/stage4/geometry/STATE.json")
    for name in ["PLANES.json","IDENTITY.json","SUMMARY.json"]:
        inputs[str(Path(geometry["run"])/name)] = dict(kind="geometry_metadata",expected_sha256=None)
    for name, spec in inputs.items():
        preserve(ROOT / name, Path("inputs") / name)
        spec["sha256"] = records[str(Path("inputs") / name)]["sha256"]
        assert spec["expected_sha256"] is None or spec["sha256"] == spec["expected_sha256"], name
    for name in ["IDENTITY.json","source_snapshot.zip","ROWS.json","SUMMARY.json"]:
        preserve(ROOT / BATCH / name, Path("batch") / name)
    versions = {}
    for package in ["numpy","scipy","h5py","vtk","pyvista","cupy-cuda12x","cupy-cuda13x"]:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            pass
    manifest = dict(schema="stage4i5_gold_v1",status="SEALED",created_UTC=datetime.now(timezone.utc).isoformat(),
        original_batch=BATCH,cases=cases,input_snapshots=inputs,files=records,
        boundary_operations_sha256=digest(out/"BOUNDARY_OPERATIONS.json"),
        repository_commit_at_sealing=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        sealer_sha256=digest(__file__),environment_observed_at_sealing=dict(python=sys.version,packages=versions),
        preservation="Independent byte copies, verified against unchanged originals; all sealed files and directories made read-only. Original results are not rerun or rewritten.",
        scientific_scope="GOLD is the numerical comparison reference, not physical truth. Existing near-contact, wall and finite-Re limitations remain. Missing lifetime instrumentation is explicit, not reconstructed as measured data.")
    lines = ["# Stage 4I gold baseline", "", "Both original cases completed naturally before the orchestration-only boundary stop. J was never started.", "",
        "| Case | A | B | Observation cutoff | Accepted steps | Runtime (h) | Mean dt (us) | Saved max constrained Re |",
        "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for label, c in cases.items():
        b = c["branch_counts"]
        lines.append(f"| {label} | {b['N_A']} | {b['N_B']} | {b['N_unassigned']} | {c['accepted_steps']} | {c['runtime_s']/3600:.6f} | {c['effective_dt']['mean_s']*1e6:.6f} | {c['trajectory_checks']['saved_constrained_slip_Re']['maximum']:.6f} |")
    lines += ["", "All 256 IDs are accounted for per case. Both trajectories span 0–6 s and contain 301 saved states. Active fields are finite; inactive force/velocity NaNs are intentional.", "",
        "The run source commit is `bb6261902dbdc4a4c3a6121144e1cad52181b3f5`; the repository commit at sealing is recorded separately. Each archived source file was checked against its original source manifest.", "",
        "Full-step maximum dt, rejected-step counter, wall-projection call counter and GPU elapsed time were not recorded. The manifest distinguishes missing measurements, code-derived counters, saved-frame diagnostics and process samples. The minimum dt can be a macro-step remainder and is not by itself a physical stiffness measure.", "",
        "Contact-event counts include repeated geometric corrections across steps, not distinct colliding pairs or fusion events. Geometry cache misses count point queries, not projection calls. Timing components overlap and must not be summed as independent costs.", "",
        "The complete source results, including raw event logs, are preserved as independent byte copies. Input snapshots and SHA256 checks prevent accidental reuse of changed geometry, fields or response tables. No 6 s baseline rerun is permitted.", "",
        "Gold status does not remove model warnings. Small-Re Stokes translation and classical far-field secondary forces are not validated in the observed finite-Re/near-contact region.", ""]
    (out / "STAGE4I_BASELINE_SUMMARY.md").write_text("\n".join(lines))
    manifest["baseline_summary_sha256"] = digest(out/"STAGE4I_BASELINE_SUMMARY.md")
    write(out / "GOLD_BASELINE_MANIFEST.json", manifest)
    (out / "MANIFEST_SHA256.txt").write_text(digest(out/"GOLD_BASELINE_MANIFEST.json") + "  GOLD_BASELINE_MANIFEST.json\n")
    for path in out.rglob("*"):
        if path.is_file():
            path.chmod(0o444)
    for path in sorted((p for p in out.rglob("*") if p.is_dir()),key=lambda p:len(p.parts),reverse=True):
        path.chmod(0o555)
    out.chmod(0o555)
    print(json.dumps(dict(status="SEALED",path=DEST,manifest_sha256=digest(out/"GOLD_BASELINE_MANIFEST.json"))),flush=True)


if __name__ == "__main__":
    main()
