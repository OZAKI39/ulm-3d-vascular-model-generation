"""Paired timestep, acoustic-grid, classification and inlet-sampling diagnostics."""
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context
import h5py
import numpy as np
from mevo_mb.acoustics.basis import ComplexAcousticField
from mevo_mb.acoustics.convergence import run_refined, load_refined
from mevo_mb.bubbles.adaptive_table import AdaptiveResponseTableBubbleModel
from mevo_mb.bubbles.forces import primary_bjerknes
from scripts.stage4.acoustic import configuration
from scripts.stage4.common import ROOT, new_run, require, sha, write, finish
from scripts.stage4.experiments import setup, run_case


def particle_worker(task):
    out,label,seed,options = task
    return run_case(out,label,seed,**options)


def compare(reference,candidate):
    a = json.loads((ROOT / reference["run"] / "BRANCH_METRICS.json").read_text())
    b = json.loads((ROOT / candidate["run"] / "BRANCH_METRICS.json").read_text())
    assert reference["initial_population_sha256"] == candidate["initial_population_sha256"]
    changed = [int(pid) for pid in a["assignment"] if a["assignment"][pid] != b["assignment"][pid]]
    ea = {e["particle_id"]:e for e in a["crossings"]}
    eb = {e["particle_id"]:e for e in b["crossings"]}
    shared = sorted(set(ea)&set(eb))
    times = np.array([abs(ea[i]["crossing_time_s"]-eb[i]["crossing_time_s"]) for i in shared])
    with h5py.File(ROOT / reference["run"] / "particles/trajectories.h5") as ha, h5py.File(ROOT / candidate["run"] / "particles/trajectories.h5") as hb:
        distance = np.linalg.norm(ha["positions_m"][-1]-hb["positions_m"][-1],axis=1)
    outcomes=["BRANCH_A","BRANCH_B","UNASSIGNED","WALL_LOSS"]
    matrix={origin:{destination:sum(a["assignment"][pid]==origin and b["assignment"][pid]==destination
        for pid in a["assignment"]) for destination in outcomes} for origin in outcomes}
    early={}
    metrics_a=json.loads((ROOT/reference["run"]/"metrics.json").read_text())["rows"]
    metrics_b=json.loads((ROOT/candidate["run"]/"metrics.json").read_text())["rows"]
    grid=np.array([r["time_s"] for r in metrics_a if r["time_s"]<=.1+1e-12])
    for key in ["mean_nearest_neighbor_m","near_contact_fraction","largest_cluster_fraction","cluster_count"]:
        ta,va=zip(*[(r["time_s"],r[key]) for r in metrics_a if r.get(key) is not None])
        tb,vb=zip(*[(r["time_s"],r[key]) for r in metrics_b if r.get(key) is not None])
        xa,xb=np.interp(grid,ta,va),np.interp(grid,tb,vb)
        early[key]=dict(maximum_absolute_difference=float(np.max(np.abs(xa-xb))),
            relative_L2_difference=float(np.linalg.norm(xa-xb)/max(np.linalg.norm(xa),1e-30)),
            reference=xa.tolist(),candidate=xb.tolist())
    na=len(a["assignment"])
    timestep={}
    for label,source in [("reference",reference),("candidate",candidate)]:
        summary=json.loads((ROOT/source["run"]/"particles/SUMMARY.json").read_text())
        timestep[label]=dict(accepted_steps=summary["accepted_steps"],
            nominal_dt_s=summary["nominal_dt_s"],minimum_effective_dt_s=summary["minimum_dt_s"],
            mean_effective_dt_s=summary["duration_s"]/summary["accepted_steps"],
            maximum_effective_dt_s=None,maximum_dt_scope="Not instrumented in retained production loop",
            observed_duration_s=summary["duration_s"])
    return dict(reference=reference["run"],candidate=candidate["run"],changed_particle_ids=changed,
                effective_timestep=timestep,
                assignment_agreement=1-len(changed)/len(a["assignment"]),
                f_A_difference=b["f_A"]-a["f_A"] if a["f_A"] is not None and b["f_A"] is not None else None,
                N_unassigned_difference=b["N_unassigned"]-a["N_unassigned"],
                A_over_initial_difference=(b["N_A"]-a["N_A"])/na,
                B_over_initial_difference=(b["N_B"]-a["N_B"])/na,
                branch_switch_matrix=matrix,early_cluster_common_times_s=grid.tolist(),early_cluster_metrics=early,
                shared_outlet_crossings=len(shared),mean_absolute_crossing_time_difference_s=float(times.mean()) if len(times) else None,
                maximum_absolute_crossing_time_difference_s=float(times.max()) if len(times) else None,
                final_or_crossing_position_RMS_difference_m=float(np.sqrt(np.mean(distance**2))),
                final_or_crossing_position_maximum_difference_m=float(distance.max()))


def main():
    decision=require("performance_refactor/backend_decision")
    backend=decision["recommended_backend_for_JK"]
    provided = setup()
    config,geometry,protocol,flow,controls,identity = provided
    out = new_run("convergence",dict(dt_factor=.5,displacement_factor=.5,grid_spacing_factor=.8,grid_count_factor=1.25,
          performance_backend=backend,backend_decision=decision["run"],
          termination_objective="branch_crossing_and_early_cluster_sensitivity_at_fixed_observation",duration_s=6.,
          diagnostic_gates=dict(minimum_per_particle_assignment_agreement=.95,
              maximum_initial_population_arrival_fraction_difference=.05,
              maximum_early_cluster_relative_L2=.1),
          classification_tolerance_m=1e-7,alternate_sampling="uniform_area",optional_N512_runtime_limit_s=120.),["primary_secondary"])
    primary = require("primary_only")["rows"]
    seed = geometry["particles"]["seed"]
    reference = {r["case"]:r for r in primary if r["seed"]==seed}
    secondary_reference={r["case"]:r for r in require("primary_secondary")["rows"]}
    report = dict(status="RUNNING",timestep={},sampling=[],population_size=[],classification={},acoustic_grid={})
    # All tasks belong to Stage 4J; independent cases use a finite spawn pool.
    # Generate the one fine field first, then share only immutable paths/identities.
    c = configuration()
    for element,amplitude,phase in zip(c["pzt"],controls["cases"]["B2"]["amplitudes_Pa"],controls["cases"]["B2"]["phases_rad"],strict=True):
        element.update(amplitude_Pa=amplitude,phase_rad=phase,enabled=amplitude>0)
    (out / "acoustic").mkdir()
    refined = run_refined(c,out)
    fine = load_refined(refined["refined_field"],c["acoustics"]["frequency_Hz"])
    field_path = out / "B2_refined_field.h5"
    fine.save(field_path)
    coarse = ComplexAcousticField.load(ROOT / controls["cases"]["B2"]["field"])
    with h5py.File(ROOT / reference["B2"]["run"] / "particles/trajectories.h5") as h:
        all_positions = h["positions_m"][:][h["active"][:]]
    selection = np.linspace(0,len(all_positions)-1,min(10000,len(all_positions)),dtype=int)
    x = all_positions[selection]
    model = AdaptiveResponseTableBubbleModel(ROOT / controls["response_table"])
    p0,p1 = coarse.interpolate_complex_pressure(x),fine.interpolate_complex_pressure(x)
    g0,g1 = coarse.interpolate_complex_gradient(x),fine.interpolate_complex_gradient(x)
    f0 = primary_bjerknes(model.response(1.78e-6,2.5e6,p0,c["medium"]),g0)
    f1 = primary_bjerknes(model.response(1.78e-6,2.5e6,p1,c["medium"]),g1)
    def rel(a,b):
        return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-30))
    report["acoustic_grid"] = dict(common_sample_count=len(x),pressure_relative_L2=rel(p0,p1),
        gradient_relative_L2=rel(g0,g1),primary_force_relative_L2=rel(f0,f1),
        coarse_pressure_range_Pa=[float(abs(p0).min()),float(abs(p0).max())],
        fine_pressure_range_Pa=[float(abs(p1).min()),float(abs(p1).max())],
        refined_field=str(field_path.relative_to(ROOT)),refined_field_sha256=sha(field_path),
        solver_metadata=refined["solver_metadata"],scope="Identical coordinates and physical domain, one direct fine control; two-grid diagnostic, not proof of grid independence")
    np.savez_compressed(out / "GRID_COMMON_POINTS.npz",positions_m=x,coarse_pressure_Pa=p0,fine_pressure_Pa=p1,
                        coarse_gradient_Pa_m=g0,fine_gradient_Pa_m=g1,coarse_force_N=f0,fine_force_N=f1)
    write(out / "PROGRESS.json",report)
    tasks = [("timestep",label,dict(dt_factor=.5,displacement_factor=.5)) for label in ["B0","B2","B3"]]
    tasks += [("timestep_secondary","B2",dict(dt_factor=.5,displacement_factor=.5,secondary=True))]
    tasks += [("acoustic_grid","B2",dict(field_override=str(field_path.relative_to(ROOT)))),
              ("classification","B2",dict(tolerance=1e-7))]
    tasks += [("sampling",label,dict(policy="uniform_area")) for label in ["B0","B2","B3"]]
    if np.median([r["elapsed_s"] for r in primary]) < 120.:
        tasks += [("population_size",label,dict(count=512)) for label in ["B0","B2","B3"]]
    else:
        report["population_size_scope"] = "N512 NOT_TESTED: measured median primary case runtime exceeded predeclared optional-run cost limit"
    common_options=dict(performance_backend=backend,wall_certificate=True,
        termination_objective="branch_crossing_and_early_cluster_sensitivity_at_fixed_observation")
    tasks=[(kind,label,dict(common_options,**options)) for kind,label,options in tasks]
    write(out/"PLANNED_TASKS.json",[dict(kind=kind,case=label,options=options) for kind,label,options in tasks])
    with ProcessPoolExecutor(max_workers=protocol["statistics"]["worker_processes"],mp_context=get_context("spawn")) as pool:
        futures = {pool.submit(particle_worker,(out,label,seed,options)):(kind,label) for kind,label,options in tasks}
        for future in as_completed(futures):
            kind,label = futures[future]
            row = future.result()
            if kind=="timestep":
                report[kind][label] = compare(reference[label],row)
            elif kind=="timestep_secondary":
                report["timestep"][label+"_primary_secondary"] = compare(secondary_reference[label],row)
            elif kind=="acoustic_grid":
                report[kind]["branch_comparison"] = compare(reference[label],row)
            elif kind=="classification":
                report[kind] = compare(reference[label],row)
            else:
                report[kind].append(row)
            write(out / "PROGRESS.json",report)
    for key in ["sampling","population_size"]:
        report[key].sort(key=lambda row:row["case"])
    stats = require("statistics")
    intervals = [stats["paired_effects"][key]["bootstrap_95CI"] for key in ["B2_minus_B0_f_A","B3_minus_B0_f_B","B2_minus_B3_f_A"]]
    lower_effect = min(ci[0] for ci in intervals) if all(ci is not None for ci in intervals) else None
    def maximum_defined(values):
        return max(abs(v) for v in values) if all(v is not None for v in values) else None
    errors = dict(timestep=maximum_defined([v["f_A_difference"] for v in report["timestep"].values()]),
                  grid=maximum_defined([report["acoustic_grid"]["branch_comparison"]["f_A_difference"]]),
                  classification=maximum_defined([report["classification"]["f_A_difference"]]))
    detailed_gates={key:dict(assignment=value["assignment_agreement"]>=.95,
        arrival_fractions=max(abs(value["A_over_initial_difference"]),abs(value["B_over_initial_difference"]))<=.05,
        early_clusters=all(v["relative_L2_difference"]<=.1 for v in value["early_cluster_metrics"].values()))
        for key,value in report["timestep"].items()}
    report.update(status="PASS_WITH_SCOPE",observed_fraction_uncertainties=errors,
                  backend=backend,backend_decision=decision["run"],
                  displacement_limit_reduced=True,secondary_timestep_case="B2",
                  branch_assignment_sensitivity={key:value["assignment_agreement"] for key,value in report["timestep"].items()},
                  detailed_timestep_gates=detailed_gates,
                  timestep_outcomes_robust=all(all(gates.values()) for gates in detailed_gates.values()),
                  paired_effect_lower_95CI=lower_effect,
                  effect_exceeds_tested_numerical_uncertainties=lower_effect is not None and all(v is not None for v in errors.values()) and lower_effect>max(errors.values()),
                  scope="Execution complete; directional/quantitative robustness is determined from numerical differences and paired intervals, not a prescribed enrichment percentage")
    finish(out,report)


if __name__ == "__main__":
    main()
