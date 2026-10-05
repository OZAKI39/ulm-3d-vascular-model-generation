"""Stage 4K: reproducibly selected critical particles, reference stays outside production."""
import ast
import copy
import json
import numpy as np
from mevo_mb.acoustics.basis import ComplexAcousticField
from mevo_mb.bubbles.population import BubblePopulation
from mevo_mb.coupling.integrators import GlobalSynchronousIntegrator
from mevo_mb.coupling.geometry import geometry_from_config
from mevo_mb.io.artifacts import new_run as particle_run
from validation.history_branch_robustness.closure import build_closure
from validation.history_branch_robustness.moving import run_reference
from scripts.stage4.common import ROOT,new_run,require,write,finish
from scripts.stage4.experiments import setup


def main():
    require("convergence")
    decision=require("performance_refactor/backend_decision")
    backend=decision["recommended_backend_for_JK"]
    config,geometry,protocol,flow,controls,identity = setup()
    primary = {row["case"]:row for row in require("primary_only")["rows"]}
    a = json.loads((ROOT / primary["B0"]["run"] / "BRANCH_METRICS.json").read_text())
    b = json.loads((ROOT / primary["B2"]["run"] / "BRANCH_METRICS.json").read_text())
    junction = {row["particle_id"]:row for row in b["bifurcation_crossings"]}
    candidates = []
    for pid,value in b["assignment"].items():
        if value not in ["BRANCH_A","BRANCH_B"]:
            continue
        switched = a["assignment"][pid] in ["BRANCH_A","BRANCH_B"] and value != a["assignment"][pid]
        y = abs(junction[int(pid)]["position_m"][1]) if int(pid) in junction else float("inf")
        candidates.append(dict(particle_id=int(pid),switch_B0_B2=switched,bifurcation_abs_y_m=y))
    candidates.sort(key=lambda row:(not row["switch_B0_B2"],row["bifurcation_abs_y_m"],row["particle_id"]))
    selected = candidates[:protocol["history_reference"]["critical_particle_count"]]
    if len(selected)<10:
        raise RuntimeError("Fewer than 10 assigned critical production particles; cannot manufacture a history agreement sample")
    out = new_run("history_branch",dict(selection_rule=protocol["history_reference"]["selection"],identity=identity,
                  selected=selected,performance_backend=backend,backend_decision=decision["run"],
                  reference="Adiabatic periodic variable-radius plus slow equilibrium-radius history",dt_s=[.001,.0005]),["convergence"])
    write(out / "SELECTION.json",dict(ordered_candidates=candidates,selected=selected,production_B2=primary["B2"]["run"]))
    forbidden = []
    for path in (ROOT / "src").rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text())):
            names = [node.module or ""] if isinstance(node,ast.ImportFrom) else [x.name for x in node.names] if isinstance(node,ast.Import) else []
            if any("history_branch_robustness" in name or "history_force_sensitivity" in name for name in names):
                forbidden.append(str(path.relative_to(ROOT)))
    assert not forbidden
    ids = np.array([row["particle_id"] for row in selected])
    with np.load(ROOT / primary["B2"]["run"] / "INITIAL_POPULATION.npz") as initial:
        positions,radii = initial["positions_m"][ids],initial["radii_m"][ids]
    field = ComplexAcousticField.load(ROOT / controls["cases"]["B2"]["field"])
    c = copy.deepcopy(config)
    c["particles"].update(primary_enabled=True,secondary_enabled=False,count=len(ids))
    c["selected_global_ids"] = ids.tolist()
    c["pzt_control"] = controls["cases"]["B2"]
    for element,amp,phase,on in zip(c["pzt"],c["pzt_control"]["amplitudes_Pa"],c["pzt_control"]["phases_rad"],c["pzt_control"]["enabled"],strict=True):
        element.update(amplitude_Pa=amp,phase_rad=phase,enabled=on)
    c["performance_execution"]=dict(backend=backend,backend_decision=decision["run"])
    c["particles"]["pair_benchmark_path"]="results/stage4/performance_refactor/protocol/PAIR_BACKEND_STATE.json"
    reduced = particle_run(out,"reduced_population_production_control",c)
    tracker=GlobalSynchronousIntegrator(c,field,flow)
    if backend!="P0":
        from scripts.stage4i5.variants import install_variant,validate_runtime_identity
        backend_runtime=install_variant(tracker,backend)
        try:
            validate_runtime_identity(backend_runtime,decision)
        except Exception:
            if hasattr(tracker,"_i5_cleanup"):
                tracker._i5_cleanup()
            raise
        write(reduced/"PERFORMANCE_BACKEND_RUNTIME.json",backend_runtime)
    from scripts.stage4.geometry_certificate import install
    install(tracker.cylinder,per_particle=True)
    try:
        production = tracker.run(BubblePopulation(positions.copy(),radii,np.arange(len(ids))),reduced)
    finally:
        if hasattr(tracker,"_i5_cleanup"):
            tracker._i5_cleanup()
    control_agreement = [production["branch_metrics"]["assignment"][str(i)] == b["assignment"][str(pid)] for i,pid in enumerate(ids)]
    closure = build_closure(out,ROOT)
    results = []
    for label,enabled,dt in [("B",False,.0005),("C",True,.001),("C",True,.0005)]:
        results.append(run_reference(out,geometry_from_config(config),flow,field,closure,positions,ids,c,label=label,dt=dt,history=enabled))
        write(out / "REFERENCES.json",results)
    reference = results[-1]
    previous = results[-2]
    rows = []
    prod_events = {row["particle_id"]:row for row in b["crossings"]}
    ref_events = {row["particle_id"]:row for row in reference["branch"]["crossings"]}
    for index,pid in enumerate(ids):
        event = ref_events.get(int(pid))
        source = prod_events[int(pid)]
        rows.append(dict(particle_id=int(pid),production_branch=b["assignment"][str(pid)],
                    history_reference_branch=reference["branch"]["assignment"][str(pid)],
                    reduced_population_control_agrees=bool(control_agreement[index]),
                    production_crossing=source,history_reference_crossing=event,
                    displacement_difference_m=float(np.linalg.norm(np.array(reference["final_positions_m"])[index]-source["position_m"])),
                    crossing_time_difference_s=event["crossing_time_s"]-source["crossing_time_s"] if event else None))
    agreement = float(np.mean([r["production_branch"]==r["history_reference_branch"] for r in rows]))
    dt_agreement = float(np.mean([reference["branch"]["assignment"][str(pid)]==previous["branch"]["assignment"][str(pid)] for pid in ids]))
    inertia = results[0]
    history_only_agreement = float(np.mean([reference["branch"]["assignment"][str(pid)]==inertia["branch"]["assignment"][str(pid)] for pid in ids]))
    reference_comparisons = {}
    for name,other in [("C_coarse_vs_C_fine",previous),("B_inertia_vs_C_history",inertia)]:
        events = {e["particle_id"]:e for e in other["branch"]["crossings"]}
        common = sorted(set(events)&set(ref_events))
        distances = np.linalg.norm(np.array(reference["final_positions_m"])-other["final_positions_m"],axis=1)
        reference_comparisons[name] = dict(
            position_RMS_difference_m=float(np.sqrt(np.mean(distances**2))),
            position_maximum_difference_m=float(distances.max()),shared_crossings=len(common),
            mean_absolute_crossing_time_difference_s=float(np.mean([abs(events[i]["crossing_time_s"]-ref_events[i]["crossing_time_s"]) for i in common])) if common else None,
            changed_global_ids=[int(pid) for pid in ids if reference["branch"]["assignment"][str(pid)]!=other["branch"]["assignment"][str(pid)]])
    classification = "robust in tested critical trajectories" if agreement>=.95 else "moderately sensitive" if agreement>=.80 else "branch assignment sensitive to translation model"
    # Any selected-sample nonconvergence or omitted-contact effect remains a warning.
    status = "PASS_WITH_SCOPE" if agreement>=.95 and dt_agreement==1. and all(control_agreement) else "WARNING"
    finish(out,dict(status=status,rows=rows,branch_agreement_rate=agreement,classification=classification,
                   branch_assignment_agreement_rate=agreement,performance_backend=backend,backend_decision=decision["run"],
                   reference_dt_assignment_agreement=dt_agreement,reduced_population_control_agreement=float(np.mean(control_agreement)),
                   B_vs_C_history_only_branch_agreement=history_only_agreement,reference_comparisons=reference_comparisons,
                   selected_count=len(ids),production_history_imports=forbidden,reference_results=results,
                   scope="Only selected B2 critical particles; adiabatic two-scale reference, not a fully microcycle-resolved moving history solution. Production-vs-C includes averaging, inertia and different radial/closure interpolation approximations; B-vs-C tests history within the same reference construction. Wall hydrodynamics remain unvalidated. No quantitative clinical fraction or experimental validation."))


if __name__ == "__main__":
    main()
