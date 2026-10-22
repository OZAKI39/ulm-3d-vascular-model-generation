"""Correctness-first I.5 decision from completed candidate evidence, never preset speed."""
import json

from scripts.stage4i5.common import BASE,GOLD,ROOT,new_run,finish,write,sha


def prerequisite_gates():
    """A short trajectory pass cannot bypass independent geometry/contact tests."""
    import vtk

    wall={}
    neighbor=None
    ordered=None
    for kind in ["p1_vtk_contact","p2_warp"]:
        for path in sorted((BASE/kind).glob("*/SUMMARY.json")):
            result=json.loads(path.read_text())
            if result.get("test")=="gold_projection_shadow":
                ordered=result
            if result.get("point_count",0)>=10000 and result.get("VTK_version")==vtk.vtkVersion.GetVTKVersion():
                wall[result["backend"]]=result
            if result.get("rows") and {r.get("kind") for r in result["rows"]}=={
                    "sparse","dense","wall_near","high_velocity","variable_radius","gold_B2","gold_B3"}:
                neighbor=result
    gates={}
    for backend,result in wall.items():
        projection=result["projection"]
        geometric=(all(r["max_distance_error_m"]<=1e-12
                       and r["max_closest_point_error_m"]<=1e-12
                       and r["max_normal_angle_deg"]<=1e-5
                       and not any(r[key] for key in ["inside_outside_mismatches",
                           "projection_decision_mismatches","velocity_constraint_decision_mismatches"])
                       for r in result["groups"])
                   and projection.get("projection_failures",1)==0
                   and projection.get("max_position_error_m",float("inf"))<=1e-12
                   and projection.get("hit_mismatches",1)==0
                   and projection.get("maximum_penetration_m",float("inf"))<=1e-12)
        gates[backend]=dict(pass_geometry=geometric,run=result["run"],
            mesh_assumptions_pass=result["mesh_audit"]["closed_consistent_outward"],
            query_speedup=result["speedup"],
            max_distance_error_m=max(r["max_distance_error_m"] for r in result["groups"]),
            max_closest_point_error_m=max(r["max_closest_point_error_m"] for r in result["groups"]),
            max_normal_angle_deg=max(r["max_normal_angle_deg"] for r in result["groups"]),
            inside_outside_mismatches=sum(r["inside_outside_mismatches"] for r in result["groups"]),
            projection_decision_mismatches=sum(r["projection_decision_mismatches"] for r in result["groups"]),
            velocity_constraint_decision_mismatches=sum(r["velocity_constraint_decision_mismatches"] for r in result["groups"]),
            full_projection=result["projection"],backend_details=result.get("backend_details",{}))
    gates["neighbor"]=dict(pass_zero_miss=bool(neighbor and all(
        r["missed_true_contacts"]==0 and r["extra_narrow_phase_contacts"]==0
        and r["sorted_pairs_identical"] for r in neighbor["rows"])),run=neighbor["run"] if neighbor else None)
    gates["ordered_contact"]=dict(pass_bitwise_shadow=bool(ordered and ordered["status"]=="PASS"
        and len(ordered["rows"])==20 and all(r["bitwise_and_events_identical"] for r in ordered["rows"])),
        run=ordered["run"] if ordered else None)
    return gates


def main():
    protocols=json.loads((BASE/"protocol/PROTOCOL.json").read_text())
    segment_ids={s["id"] for s in protocols["segments"]}
    comparisons={}
    for path in sorted((BASE/"full_regression").glob("*/SUMMARY.json")):
        row=json.loads(path.read_text())
        comparisons[row["backend"]]=row
    for name in ["P1_exact_cache","P1_neighbor","P1_combined","P1_static","P2_warp","P1_ordered"]:
        if name not in comparisons:
            raise RuntimeError("Missing candidate gate report: "+name)
    limiter=json.loads((BASE/"limiter_analysis/STATE.json").read_text())
    if limiter["status"]!="COMPLETE":
        raise RuntimeError("Limiter analysis incomplete")
    eligible_local=[r for r in limiter["segments"] if r["localized_few_particle_limiter"]]
    if eligible_local:
        raise RuntimeError("Localized limiters detected: explicitly review buffered active-set/P3 feasibility before deciding")
    prerequisites=prerequisite_gates()
    outcomes=[]
    for backend,result in comparisons.items():
        rows=result["rows"]
        complete={r["segment"] for r in rows}==segment_ids
        independent_gate=True
        if backend in {"P1_neighbor","P1_combined"}:
            independent_gate=prerequisites["neighbor"]["pass_zero_miss"]
        if backend=="P1_static":
            independent_gate=prerequisites.get("vtk",{}).get("pass_geometry",False)
        if backend=="P2_warp":
            gate=prerequisites.get("warp",{})
            independent_gate=gate.get("pass_geometry",False) and gate.get("mesh_assumptions_pass",False)
        if backend=="P1_ordered":
            independent_gate=prerequisites["ordered_contact"]["pass_bitwise_shadow"]
        correct=independent_gate and complete and all(r["strict_P1_gate"]=="PASS" for r in rows)
        speed=sum(r["baseline_wall_s"] for r in rows)/sum(r["candidate_wall_s"] for r in rows) if complete else None
        outcomes.append(dict(backend=backend,complete_segments=len(rows),strict_correctness=correct,
            aggregate_speedup=speed,gate_report=result["run"],independent_fixture_gate=independent_gate,
            status="ELIGIBLE" if correct else "REJECTED_OR_INCOMPLETE",regression=rows))
    valid=[r for r in outcomes if r["strict_correctness"] and r["aggregate_speedup"]>1.05]
    timing=None
    if valid:
        from scripts.stage4i5.timing_confirmation import confirm
        timing=confirm(valid)
        evidence={row["backend"]:row for row in timing["candidates"]}
        for row in outcomes:
            row["minimal_timer_confirmation"]=evidence.get(row["backend"])
        valid=[row for row in valid if evidence[row["backend"]]["strict_correctness"]
               and evidence[row["backend"]]["aggregate_speedup"]>1.05]
    selected=max(valid,key=lambda r:r["minimal_timer_confirmation"]["aggregate_speedup"])["backend"] if valid else "P0"
    # Float32 Warp is never adopted through the strict-P1 shortcut; its own
    # geometrical and trajectory gate must be reviewed explicitly if it ever passes.
    if selected=="P2_warp":
        raise RuntimeError("Warp passed short comparison: explicit mesh/normal/penetration and transfer review required")
    out=new_run("backend_decision",dict(correctness_first=True,noise_margin_for_speed_selection=.05))
    source=json.loads((GOLD/"cases/B2/SOURCE_IDENTITY.json").read_text())["source_sha256"]
    physics={k:v for k,v in source.items() if k.startswith("src/")}
    changed=[k for k,v in physics.items() if sha(ROOT/k)!=v]
    if changed:
        raise RuntimeError("Original production physics source changed: "+str(changed))
    p3=dict(status="NOT_IMPLEMENTED_DATA_DO_NOT_SUPPORT_LOCAL_SUBCYCLING",
        segment_evidence=[dict(id=r["segment"]["id"],ids_covering_90_percent=r["particles_covering_90_percent"],
            median_stiff=r["median_particles_requiring_substeps"],active=r["median_active_particles"],
            median_contact_component=r["median_limiter_contact_component"]) for r in limiter["segments"]],
        scope="No experimental local solver or full optimized P3 runs launched. Global all-pair forces and widespread stiffness make a small active subset unsupported by these samples.")
    write(out/"P3_FEASIBILITY.json",p3)
    p3dir=BASE/"p3_local_subcycling"
    p3dir.mkdir(exist_ok=True)
    write(p3dir/"STATE.json",dict(p3,evidence=limiter["run"],decision=str(out.relative_to(ROOT))))
    near=BASE/"near_contact_sensitivity"
    near.mkdir(exist_ok=True)
    write(near/"STATE.json",dict(status="NOT_RUN_OPTIONAL",reason="Near-contact policy unchanged throughout main performance regression and J/K; no force cutoff substituted"))
    lines=["# Stage 4I.5 backend decision","",f"recommended_backend_for_JK: **{selected}**","",
        "Correctness takes priority. Candidates need all four completed S1/S2 regressions and the strict trajectory/contact gate. Aggregate speedup uses summed matched solver times. A component-timed gain greater than 5% is a screening criterion; adoption additionally requires the same four windows with per-component wrappers removed, strict timer-removal/state regression and a confirmed gain greater than 5%. This controls the API-dependent cost of millions of per-point VTK timer calls. The 2× engineering target is not a scientific acceptance threshold.","",
        "| Backend | Completed segments | Strict gate | Aggregate speedup |","|---|---:|---|---:|"]
    for r in outcomes:
        lines.append(f"| {r['backend']} | {r['complete_segments']} | {r['strict_correctness']} | {r['aggregate_speedup'] if r['aggregate_speedup'] is not None else 'Not comparable'} |")
    lines += ["", "## Minimal-timer confirmation", "",
        (f"Evidence: {timing['run']}/SUMMARY.json. Common accepted-step diagnostic recording, dt and classifier timers remain; per-component and per-point wall wrappers are disabled. These short replay timings are not full-production runtime measurements." if timing
         else "Not run: no strict candidate exceeded the predeclared 5% screening margin."),"",
        "| Backend | Strict confirmation | Aggregate speedup with component wrappers removed |",
        "|---|---|---:|"]
    for r in (timing["candidates"] if timing else []):
        lines.append(f"| {r['backend']} | {r['strict_correctness']} | {r['aggregate_speedup']} |")
    lines += ["", "## Independent wall fixture", "",
        "| Backend | Geometry gate | Query speedup | Max distance error [m] | Max normal angle [deg] | Inside mismatches | Projection decisions | Velocity decisions |",
        "|---|---|---:|---:|---:|---:|---:|---:|"]
    wall_notes=[]
    for name in ["vtk","warp"]:
        if name in prerequisites:
            r=prerequisites[name]
            lines.append(f"| {name} | {r['pass_geometry']} | {r['query_speedup']:.5g} | {r['max_distance_error_m']:.5g} | {r['max_normal_angle_deg']:.5g} | {r['inside_outside_mismatches']} | {r['projection_decision_mismatches']} | {r['velocity_constraint_decision_mismatches']} |")
            wall_notes.extend(["",f"{name} full projection: `{json.dumps(r['full_projection'])}`. Evidence: {r['run']}/SUMMARY.json.",""])
    lines.extend(wall_notes)
    lines += ["", "## Original physics and terminal boundaries", "",
        f"All {len(physics)} original production source files match sealed gold byte-for-byte. New optional modules do not change force laws, constraints or physical tolerances.",
        "Existing assigned-particle removal is an absorbing model boundary, not proven terminal-outlet physics. No additional deletion or endpoint shortcut is enabled. J retains the 6 s observation with objective-policy validation; a plateau with unassigned particles is insufficient to stop.","",
        "## Local subcycling", "", p3["scope"], json.dumps(p3["segment_evidence"],indent=2),"",
        "## Scope of gold regression", "",
        "The original B2/B3 are the unique full gold runs and were never rerun. P0 checkpoints are independently reconstructed and compared with gold endpoints. P1 tests compare every accepted state, dt, contact hash, wall flags and assignments with matched P0 short replays. No measured full-6-s optimized branch agreement is claimed when no such run was performed. Most dense windows contain no new outlet crossing; their strongest tests are state/contact equivalence.","",
        "## Warp and SDF", "",
        "Warp remains optional/benchmark-only unless its explicit mesh, projection, trajectory and runtime gates all pass. Float32 geometry is never made acceptable by loosening the original tolerances. SDF is not substituted: the new profile attributes major dense-case cost to contact projection/NNLS, so another approximate wall representation does not address that dominant cost."]
    (out/"STAGE4I5_BACKEND_DECISION.md").write_text("\n".join(lines)+"\n")
    finish(out,dict(status="PASS_WITH_SCOPE",recommended_backend_for_JK=selected,candidates=outcomes,
        original_physics_files_unchanged=physics,gold_manifest_sha256=sha(GOLD/"GOLD_BASELINE_MANIFEST.json"),
        p3=p3,limiter_analysis=limiter["run"],independent_prerequisite_gates=prerequisites,
        minimal_timer_confirmation=timing,new_terminal_deactivation_enabled=False,
        early_termination_enabled=False,full_optimized_P3_runs="NOT_RUN_NOT_ELIGIBLE",
        numerical_gate_scope="Four strict saved-state replay regressions, not new physical validation or a full optimized 6 s rerun"))


if __name__=="__main__":
    main()
