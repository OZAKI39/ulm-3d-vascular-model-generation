"""Immutable paired Y-control runs; production forces and FlowField API are reused."""
import argparse
import copy
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context
import numpy as np
import yaml
from mevo_mb.acoustics.basis import ComplexAcousticField
from mevo_mb.bubbles.inlet_sampling import initialize_inlet
from mevo_mb.coupling.bifurcation_experiment import BifurcationExperimentTracker
from mevo_mb.coupling.geometry import geometry_from_config
from mevo_mb.flow.fields import VTKFlowField
from mevo_mb.io.artifacts import new_run as particle_run
from mevo_mb.utils.config import load_config
from scripts.stage4.common import ROOT, new_run, require, sha, write, finish
from scripts.stage4.synthetic import AcousticOff


def setup():
    geometry = yaml.safe_load((ROOT / "configs/stage4/bifurcation.yaml").read_text())
    protocol = yaml.safe_load((ROOT / "configs/stage4/experiment.yaml").read_text())
    source = require("simvascular/steady")
    imported = require("simvascular/import")
    mesh = require("simvascular/mesh")
    acoustic = require("acoustic")
    controls = require("phase_screening")
    config = load_config(ROOT / "configs/stage4/acoustic.yaml")
    config.pop("vessel",None)
    config["geometry"] = dict(geometry["geometry"],kind="bifurcation",
        mesh_surface=str((ROOT / mesh["mesh"]).parent / "mesh-complete.exterior.vtp"))
    config["medium"] = geometry["medium"]
    config["particles"].update(protocol["particles"])
    config["output"].update(protocol["output"])
    config["classification_tolerance_m"] = geometry["geometry"]["classification_tolerance_m"]
    config["bubble_response"]["path"] = controls["response_table"]
    config["acoustics"]["geometry_hash"] = source["geometry_hash"]
    config["flow"] = dict(backend="vtk",snapshot=source["snapshot"],velocity_array=source["velocity_array"],
                         coordinate_scale_to_m=source["coordinate_scale_to_m"],
                         velocity_scale_to_m_s=source["velocity_scale_to_m_s"],spatial_method=imported["spatial_method"],time_unit="s")
    config["validation_gates"] = dict(stage="Stage4",baseline_tag="stage3-minimal-validated",
        geometry="results/stage4/geometry/STATE.json",flow_import="results/stage4/simvascular/import/STATE.json",
        acoustic="results/stage4/acoustic/STATE.json",controls="results/stage4/phase_screening/STATE.json")
    config["experiment_scope"] = "Stage4 finite initial ensemble, homogeneous water, geometric walls; no clinical concentration or patient calibration"
    flow = VTKFlowField(ROOT / source["snapshot"],velocity_array=source["velocity_array"],
                        coordinate_scale_to_m=source["coordinate_scale_to_m"],
                        velocity_scale_to_m_s=source["velocity_scale_to_m_s"],spatial_method=imported["spatial_method"])
    identity = dict(geometry_hash=source["geometry_hash"],mesh_sha256=source["mesh_sha256"],
                    flow_sha256=source["snapshot_sha256"],acoustic_basis_sha256=acoustic["basis_sha256"],
                    response_table_sha256=sha(ROOT / controls["response_table"]))
    return config,geometry,protocol,flow,controls,identity


def run_case(parent,label,seed,*,secondary=False,dt_factor=1.,policy=None,count=None,
             field_override=None,tolerance=None,duration=None,provided=None,wall_certificate=False,
             performance_backend=None,displacement_factor=1.,termination_objective=None):
    config,geometry,protocol,flow,controls,identity = provided if provided is not None else setup()
    config = copy.deepcopy(config)
    if wall_certificate:
        config["performance_execution"] = dict(wall_projection="per_particle_interior_certificate_v1",
            distance_roundoff_margin_m=1e-10, changed_force_or_timestep=False,
            helper_sha256=sha(ROOT/"scripts/stage4/geometry_certificate.py"),
            equivalence_evidence=require("performance/later_particle_geometry_certificate")["run"])
    if performance_backend is not None:
        decision=require("performance_refactor/backend_decision")
        if performance_backend!=decision["recommended_backend_for_JK"]:
            raise ValueError("Formal J/K backend must match the reviewed I.5 recommendation")
        config.setdefault("performance_execution",{}).update(
            backend=performance_backend,backend_decision=decision["run"],
            backend_decision_sha256=sha(ROOT/decision["run"]/"SUMMARY.json"))
        config["particles"]["pair_benchmark_path"]="results/stage4/performance_refactor/protocol/PAIR_BACKEND_STATE.json"
    if not np.isfinite(displacement_factor) or displacement_factor<=0:
        raise ValueError("Positive finite displacement factor required")
    config["particles"]["max_displacement_per_radius"]*=displacement_factor
    if termination_objective is not None:
        from dataclasses import asdict
        from mevo_mb.coupling.termination import SimulationTerminationPolicy
        termination=SimulationTerminationPolicy(objective=termination_objective)
        config["termination_policy"]=asdict(termination)
    sampling = dict(geometry["particles"],seed=int(seed))
    if policy is not None:
        sampling["sampling_policy"] = policy
    if count is not None:
        sampling["count"] = int(count)
    config["particles"].update(primary_enabled=label != "B0",secondary_enabled=secondary,
                              dt_s=config["particles"]["dt_s"]*dt_factor,count=sampling["count"])
    if duration is not None:
        config["particles"]["duration_s"] = duration
    if tolerance is not None:
        config["classification_tolerance_m"] = tolerance
    config["seed"],config["sampling"] = int(seed),sampling
    control = controls["cases"][label]
    config["pzt_control"] = control
    for element,amplitude,phase,enabled in zip(config["pzt"],control["amplitudes_Pa"],control["phases_rad"],control["enabled"],strict=True):
        element.update(amplitude_Pa=amplitude,phase_rad=phase,enabled=enabled)
    load_start = time.perf_counter()
    field = AcousticOff() if label == "B0" else ComplexAcousticField.load(ROOT / (field_override or control["field"]))
    load_seconds = time.perf_counter()-load_start
    config["effective_acoustic_field"] = dict(path=(field_override or control["field"]) if label!="B0" else None,
        basis_override=field_override is not None,
        roi_grid_shape=[len(axis) for axis in field.axes_m] if label!="B0" else None,
        roi_grid_spacing_m=[float(axis[1]-axis[0]) for axis in field.axes_m] if label!="B0" else None,
        scope="Loaded field is authoritative for interpolation; base grid config identifies the original basis. Overrides retain their own field grid and direct-solver metadata.")
    g = geometry_from_config(config)
    population = initialize_inlet(g,flow,sampling)
    out = particle_run(parent,label+"_seed"+str(seed),config)
    initial = out / "INITIAL_POPULATION.npz"
    np.savez_compressed(initial,positions_m=population.positions_m,radii_m=population.radii_m,particle_id=population.ids)
    import hashlib
    initial_hash = hashlib.sha256(population.positions_m.tobytes()+population.radii_m.tobytes()+population.ids.tobytes()).hexdigest()
    write(out / "RUN_IDENTITY.json",dict(identity,particle_seed=int(seed),initial_population_sha256=initial_hash,
           acoustic_field_sha256=sha(ROOT / (field_override or control["field"])) if label != "B0" else None,
           acoustic_field_path=config["effective_acoustic_field"]["path"],response_table_path=controls["response_table"],
           acoustic_OFF=label == "B0",control_case=label,secondary_enabled=secondary,config=config))
    if performance_backend is None:
        tracker = BifurcationExperimentTracker(config,field,flow)
    else:
        from mevo_mb.coupling.integrators import GlobalSynchronousIntegrator
        tracker=GlobalSynchronousIntegrator(config,field,flow)
        if performance_backend!="P0":
            from scripts.stage4i5.variants import install_variant,validate_runtime_identity
            backend_runtime=install_variant(tracker,performance_backend)
            try:
                validate_runtime_identity(backend_runtime,decision)
            except Exception:
                if hasattr(tracker,"_i5_cleanup"):
                    tracker._i5_cleanup()
                raise
            write(out/"PERFORMANCE_BACKEND_RUNTIME.json",backend_runtime)
    certificate = None
    if wall_certificate:
        from scripts.stage4.geometry_certificate import install
        certificate = install(tracker.cylinder,per_particle=True)
    try:
        summary = tracker.run(population,out)
    finally:
        if hasattr(tracker,"_i5_cleanup"):
            tracker._i5_cleanup()
    if termination_objective is not None:
        reason=termination.reason(time_s=summary["duration_s"],observation_limit_s=config["particles"]["duration_s"],
            total_particles=len(population.ids),terminally_assigned=summary["branch_metrics"]["N_assigned"])
        if reason is None:
            write(out/"TERMINATION_POLICY_VIOLATION.json",dict(status="FAIL",
                reason="Legacy absorbing population became empty before the objective observation limit; physical terminality not proven",
                policy=config["termination_policy"]))
            raise RuntimeError("Run does not satisfy its objective-driven termination policy")
        summary["termination_reason"]=reason
        summary["termination_policy"]=config["termination_policy"]
    summary["timings"]["acoustic_loading_s"] = load_seconds
    if certificate is not None:
        summary["wall_projection_certificate"] = certificate
    write(out / "particles/SUMMARY.json",summary)
    row = dict(case=label,seed=int(seed),secondary=secondary,run=str(out.relative_to(ROOT)),
                N_A=summary["branch_metrics"]["N_A"],N_B=summary["branch_metrics"]["N_B"],
                N_unassigned=summary["branch_metrics"]["N_unassigned"],N_wall_loss=summary["branch_metrics"]["N_wall_loss"],
                f_A=summary["branch_metrics"]["f_A"],f_B=summary["branch_metrics"]["f_B"],
                initial_population_sha256=initial_hash,elapsed_s=summary["elapsed_s"],timings=summary["timings"])
    print(json.dumps(row),flush=True)
    return row


def _run_worker(arguments):
    return run_case(*arguments)


def _run_secondary_worker(arguments):
    return run_case(*arguments[:3],secondary=True,wall_certificate=arguments[3] if len(arguments)>3 else False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stage",choices=["primary","seeds","secondary"])
    parser.add_argument("--reuse-completed",help="Existing ROWS.json containing only fully completed, identity-matched cases")
    parser.add_argument("--certified-wall",action="store_true",help="Use the verified exact interior-wall certificate for secondary cases")
    args = parser.parse_args()
    if args.certified_wall and args.stage != "secondary":
        parser.error("--certified-wall is limited to the secondary comparison")
    provided = setup()
    config,geometry,protocol,flow,controls,identity = provided
    kind = {"primary":"primary_only","seeds":"statistics/paired_runs","secondary":"primary_secondary"}[args.stage]
    prereq = {"primary":["phase_screening"],"seeds":["primary_only"],"secondary":["statistics"]}[args.stage]
    out = new_run(kind,dict(config=config,protocol=protocol,identity=identity,
                           wall_projection_certificate=args.certified_wall),prereq)
    seeds = geometry["particles"]["seeds"] if args.stage == "seeds" else [geometry["particles"]["seed"]]
    labels = ["B2","B3"] if args.stage == "secondary" else ["B0","B1","B2","B3"]
    rows = []
    # Seed 0 primary trials are reused by identity; no duplicated scientific run.
    existing = require("primary_only")["rows"] if args.stage == "seeds" else []
    if args.reuse_completed:
        assert args.stage == "primary"
        existing = json.loads((ROOT / args.reuse_completed).read_text())
        for row in existing:
            folder = ROOT / row["run"]
            recorded = json.loads((folder / "RUN_IDENTITY.json").read_text())
            assert all(recorded[key] == value for key,value in identity.items())
            assert recorded["config"]["pzt_control"] == controls["cases"][row["case"]]
            assert not recorded["secondary_enabled"]
            for key in ["dt_s","duration_s","max_displacement_per_radius","integrator"]:
                assert recorded["config"]["particles"][key] == config["particles"][key]
            import h5py
            with h5py.File(folder / "particles/trajectories.h5") as h:
                assert h.attrs["complete"]
        write(out / "REUSED_COMPLETED_CASES.json",dict(source=args.reuse_completed,rows=existing,
              justification="Scientific identities/protocol checked; exact geometry-cache equivalence verified separately"))
    if args.stage == "seeds":
        rows.extend(existing)
        tasks = [(out,label,seed) for seed in seeds for label in labels
                 if not any(r["seed"] == seed and r["case"] == label for r in existing)]
        failures = []
        with ProcessPoolExecutor(max_workers=protocol["statistics"]["worker_processes"],mp_context=get_context("spawn")) as pool:
            futures = {pool.submit(_run_worker,task):task for task in tasks}
            for future in as_completed(futures):
                try:
                    rows.append(future.result())
                except Exception as error:
                    _,label,seed = futures[future]
                    failures.append(dict(case=label,seed=seed,error=repr(error)))
                    write(out / "FAILED_CASES.json",failures)
                    print(json.dumps(dict(incomplete_case=label,seed=seed,error=repr(error))),flush=True)
                write(out / "ROWS.json",rows)
        if failures:
            finish(out,dict(status="FAIL",rows=rows,identity=identity,failures=failures,
                            scope="Incomplete cases excluded from completed list; no statistical acceptance of a partial seed batch"))
            raise RuntimeError("Paired batch has incomplete cases; preserved completed runs and failure identities")
        for seed in seeds:
            assert len({r["initial_population_sha256"] for r in rows if r["seed"] == seed}) == 1
        finish(out,dict(status="PASS_WITH_SCOPE",rows=rows,identity=identity,
                        scope="Paired-seed executions complete; effect assessed separately; independent cases use finite worker pool"))
        return
    if args.stage == "secondary":
        with ProcessPoolExecutor(max_workers=min(2,protocol["statistics"]["worker_processes"]),mp_context=get_context("spawn")) as pool:
            futures = [pool.submit(_run_secondary_worker,(out,label,seeds[0],args.certified_wall)) for label in labels]
            for future in as_completed(futures):
                rows.append(future.result())
                write(out / "ROWS.json",rows)
        rows.sort(key=lambda row:row["case"])
        assert len({r["initial_population_sha256"] for r in rows})==1
        finish(out,dict(status="PASS_WITH_SCOPE",rows=rows,identity=identity,
                        scope="Secondary executions complete on one shared seed; effect and validity quantified separately"))
        return
    for seed in seeds:
        for label in labels:
            previous = next((r for r in existing if r["seed"] == seed and r["case"] == label),None)
            row = previous if previous else run_case(out,label,seed,secondary=args.stage == "secondary",provided=provided)
            rows.append(row)
            write(out / "ROWS.json",rows)
        hashes = {r["initial_population_sha256"] for r in rows if r["seed"] == seed}
        if len(hashes) != 1:
            raise RuntimeError("Paired controls do not share exact initial IDs/positions")
    finish(out,dict(status="PASS_WITH_SCOPE",rows=rows,identity=identity,
                    scope="Execution and paired identity passed; effect and robustness are assessed separately, unassigned outcomes retained"))


if __name__ == "__main__":
    main()
