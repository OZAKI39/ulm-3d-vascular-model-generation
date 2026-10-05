"""Strict short-replay gates; incomplete/failed candidates never become production."""
import argparse
import json

import h5py
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree

from mevo_mb.coupling.bifurcation import BifurcationGeometry

from scripts.stage4i5.common import BASE,GOLD,ROOT,new_run,finish


def comparison(reference,candidate):
    a=json.loads((reference/"SUMMARY.json").read_text())
    b=json.loads((candidate/"SUMMARY.json").read_text())
    if a["segment"]!=b["segment"] or a["status"]!="COMPLETE" or b["status"]!="COMPLETE":
        raise ValueError("Only completed matching predeclared segments may pass")
    with h5py.File(reference/"replay.h5") as h:
        x=h["positions_m"][:]
        times=h["rows"][:,0]
        rows=h["rows"][:]
        final=h["final_positions_m"][:]
        ids=h["particle_id"][:]
        initial_active=h["active"][0]
        mask=h["active"][:]
        contacts=h["contact_sha256"][:]
        wall=h["wall_hits"][:]
    with h5py.File(candidate/"replay.h5") as h:
        y=h["positions_m"][:]
        time_b=h["rows"][:,0]
        rows_b=h["rows"][:]
        mask_b=h["active"][:]
        final_b=h["final_positions_m"][:]
        same_ids=bool(np.array_equal(ids,h["particle_id"][:]))
        same_steps=bool(np.array_equal(times,time_b))
        same_dt=bool(same_steps and np.array_equal(rows[:,2],h["rows"][:,2]))
        same_contacts=bool(same_steps and np.array_equal(contacts,h["contact_sha256"][:]))
        same_wall=bool(same_steps and np.array_equal(wall,h["wall_hits"][:]))
        same_masks=bool(same_steps and np.array_equal(mask,h["active"][:]))
        same_flags=bool(same_steps and np.array_equal(rows[:,12:17],h["rows"][:,12:17]))
        radii=h["radii_m"][:]
        final_ids=h["final_active_ids"][:]
    config=json.loads((GOLD/"cases"/b["segment"]["case"]/"RUN_IDENTITY.json").read_text())["config"]
    geometry=BifurcationGeometry(config["geometry"])
    relative=config["geometry"]["mesh_surface"].split("/workspace/mevo_acoustic_microbubbles/")[-1]
    geometry.bind_mesh_surface(pv.read(GOLD/"inputs"/relative))
    maximum_penetration=maximum_overlap=0.
    frames=[(frame[active_frame],radii[active_frame]) for frame,active_frame in zip(y,mask_b)]
    frames.append((final_b[final_ids],radii[final_ids]))
    for xx,rr in frames:
        if not len(rr):
            continue
        maximum_penetration=max(maximum_penetration,float(np.maximum(rr-geometry.wall_distance(xx),0).max()))
        pairs=cKDTree(xx).query_pairs(np.nextafter(2*rr.max(),np.inf),output_type="ndarray")
        if len(pairs):
            i,j=pairs.T
            maximum_overlap=max(maximum_overlap,float(np.max(rr[i]+rr[j]-np.linalg.norm(xx[j]-xx[i],axis=1))))
    wall_safe=maximum_penetration<=1e-12
    contact_safe=maximum_overlap<=config["particles"]["contact_tolerance_m"]
    if same_steps:
        delta=np.linalg.norm(x[:,initial_active]-y[:,initial_active],axis=-1)
        position_method="Identical accepted times, initially active IDs only"
    else:
        extended_times=np.r_[time_b,b["end_time_s"]]
        extended_x=np.concatenate([y,final_b[None]],axis=0)
        interpolated=np.empty_like(x[:,initial_active])
        for j,pid in enumerate(np.flatnonzero(initial_active)):
            for axis in range(3):
                interpolated[:,j,axis]=np.interp(times,extended_times,extended_x[:,pid,axis])
        delta=np.linalg.norm(x[:,initial_active]-interpolated,axis=-1)
        position_method="Linear time alignment of candidate accepted states; initially active IDs only, not bitwise"
    max_position=max(float(delta.max()),float(np.linalg.norm(final[initial_active]-final_b[initial_active],axis=1).max()))
    rms=float(np.sqrt(np.mean(delta**2)))
    branch_a=json.loads((reference/"branch.json").read_text())["assignment"]
    branch_b=json.loads((candidate/"branch.json").read_text())["assignment"]
    changed=[int(pid) for pid in branch_a if branch_a[pid]!=branch_b[pid]]
    crossings_a=json.loads((reference/"crossings.json").read_text())
    crossings_b=json.loads((candidate/"crossings.json").read_text())
    events_a={event["particle_id"]:event for event in crossings_a}
    events_b={event["particle_id"]:event for event in crossings_b}
    common_crossings=sorted(set(events_a)&set(events_b))
    crossing_errors=[abs(events_a[i]["crossing_time_s"]-events_b[i]["crossing_time_s"]) for i in common_crossings]
    validity_keys=["validity_status","physics_status","invalid_reason","within_response_table",
        "NEAR_CONTACT_PHYSICS_APPROXIMATED","invalid_accepted_duration_s"]
    validity_flags=all(a["physics_validity"].get(key)==b["physics_validity"].get(key) for key in validity_keys)
    strict=(wall_safe and contact_safe and same_ids and same_dt and same_contacts and same_wall and same_masks and same_flags and validity_flags
            and not changed and max_position<=1e-12 and crossings_a==crossings_b)
    return dict(segment=a["segment"]["id"],backend=b["backend"],reference=str(reference.relative_to(ROOT)),
        candidate=str(candidate.relative_to(ROOT)),strict_P1_gate="PASS" if strict else "FAIL",
        particle_ids_identical=same_ids,accepted_times_identical=same_steps,dt_identical=same_dt,
        contact_sequence_identical=same_contacts,wall_hits_identical=same_wall,active_masks_identical=same_masks,
        accepted_regime_and_flag_rows_identical=same_flags,new_crossing_events_identical=crossings_a==crossings_b,
        branch_agreement_all_ids=1-len(changed)/len(ids),changed_particle_ids=changed,
        initially_active_particles=int(initial_active.sum()),new_crossings_reference=len(crossings_a),new_crossings_candidate=len(crossings_b),
        common_new_crossings=len(common_crossings),
        maximum_new_crossing_time_error_s=max(crossing_errors) if crossing_errors else None,
        max_position_error_m=max_position,rms_position_error_m=rms,position_comparison=position_method,
        bitwise_positions=bool(same_steps and np.array_equal(x,y) and np.array_equal(final,final_b)),
        contact_pair_event_difference=b["contact_pair_events"]-a["contact_pair_events"],
        maximum_original_mesh_penetration_m=maximum_penetration,wall_tolerance_m=1e-12,
        maximum_remaining_pair_overlap_m=maximum_overlap,contact_tolerance_m=config["particles"]["contact_tolerance_m"],
        wall_feasibility_pass=wall_safe,contact_feasibility_pass=contact_safe,
        baseline_wall_s=a["solver_wall_s"],candidate_wall_s=b["solver_wall_s"],
        speedup=a["solver_wall_s"]/b["solver_wall_s"],accepted_step_difference=len(rows_b)-len(rows),
        speed_measurement=dict(cprofile=False,
            reference_component_timers=a.get("component_timers",True),
            candidate_component_timers=b.get("component_timers",True),
            scope="Runtime ratios across different instrumentation are diagnostic only"),
        validity_summary_exact=a["physics_validity"]==b["physics_validity"],
        validity_flag_summary_equal=validity_flags,
        limitation="No new branch crossing in a segment makes that segment's branch agreement a weak check; trajectory/contact gates still apply.")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--backend",required=True)
    args=parser.parse_args()
    runs={}
    for path in sorted((BASE/"p1_vtk_contact").glob("*/SUMMARY.json")):
        s=json.loads(path.read_text())
        if s.get("status")=="COMPLETE" and not s.get("cprofile"):
            runs[(s["backend"],s["segment"]["id"])]=path.parent
    out=new_run("full_regression",dict(scope="Short-segment regression, not full 6 s optimized results",backend=args.backend))
    results=[]
    for (backend,segment),path in runs.items():
        if backend==args.backend and ("P0",segment) in runs:
            results.append(comparison(runs[("P0",segment)],path))
    status="PASS" if len(results)==4 and all(x["strict_P1_gate"]=="PASS" for x in results) else "FAIL_OR_INCOMPLETE"
    finish(out,dict(status=status,backend=args.backend,rows=results))


if __name__=="__main__":
    main()
