"""Postprocess limiter/validity evidence without perturbing the scientific replay."""
from collections import Counter,defaultdict
import json

import h5py
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree

from mevo_mb.coupling.bifurcation import BifurcationGeometry
from mevo_mb.flow.fields import VTKFlowField
from scripts.stage4i5.common import BASE,GOLD,new_run,finish,write


def component_size(pairs,particle,n):
    parent=np.arange(n)
    def root(i):
        while parent[i]!=i:
            i=parent[i]
        return i
    for a,b in pairs:
        a,b=root(a),root(b)
        parent[b]=a
    target=root(particle)
    return sum(root(i)==target for i in range(n))


def analyze(path):
    summary=json.loads((path/"SUMMARY.json").read_text())
    case=summary["segment"]["case"]
    config=json.loads((GOLD/"cases"/case/"RUN_IDENTITY.json").read_text())["config"]
    branch=json.loads((GOLD/"cases"/case/"BRANCH_METRICS.json").read_text())
    assigned_at={event["particle_id"]:event["crossing_time_s"] for event in branch["crossings"] if "crossing_time_s" in event}
    geom=config["geometry"]
    geometry=BifurcationGeometry(geom)
    relative=geom["mesh_surface"].split("/workspace/mevo_acoustic_microbubbles/")[-1]
    geometry.bind_mesh_surface(pv.read(GOLD/"inputs"/relative))
    with h5py.File(path/"replay.h5") as h:
        rows=h["rows"][:]
        pos=h["positions_m"][:]
        active=h["active"][:]
        velocity=h["total_velocity_m_s"][:]
        forces=h["forces_N"][:]
        radii=h["radii_m"][:]
    limiter=rows[:,5].astype(int)
    points=pos[np.arange(len(pos)),limiter]
    wall=geometry.wall_distance(points)
    flow_config=config["flow"]
    flow=VTKFlowField(GOLD/"inputs"/flow_config["snapshot"],velocity_array=flow_config["velocity_array"],
        coordinate_scale_to_m=flow_config["coordinate_scale_to_m"],
        velocity_scale_to_m_s=flow_config["velocity_scale_to_m_s"],spatial_method=flow_config["spatial_method"])
    # Stage4 flow is steady. Re uses slip relative to that exact imported flow,
    # not lab-frame transport speed (the displacement limiter uses total speed).
    fluid_velocity=flow.velocity(points,0.)
    counts=Counter()
    durations=defaultdict(float)
    categories=Counter()
    details=[]
    stiff_counts=[]
    for k,row in enumerate(rows):
        pid=int(row[5])
        ids=np.flatnonzero(active[k])
        local=int(np.flatnonzero(ids==pid)[0])
        x,r=pos[k,ids],radii[ids]
        pairs=cKDTree(x).query_pairs(np.nextafter(2*r.max()+config["particles"]["contact_tolerance_m"],np.inf),output_type="ndarray")
        if len(pairs):
            i,j=pairs.T
            pairs=pairs[np.linalg.norm(x[j]-x[i],axis=1)<=r[i]+r[j]+config["particles"]["contact_tolerance_m"]]
        size=component_size(pairs,local,len(x))
        ratios=np.linalg.norm(x-x[local],axis=1)/(r+r[local])
        ratios[local]=np.inf
        near=float(ratios.min())<2
        near_wall=wall[k]/radii[pid]<5
        labels=[]
        if size>1:
            labels.append("contact")
        if size>2:
            labels.append("cluster")
        if near_wall:
            labels.append("near_wall")
        if near:
            labels.append("near_contact_validity_region")
        if not labels:
            labels.append("free")
        particle_limited=row[4]==1
        assigned=assigned_at.get(pid,np.inf)<=row[0]
        if assigned:
            raise RuntimeError("A previously assigned gold particle is still limiting the replay")
        if particle_limited:
            counts[pid]+=1
            durations[pid]+=row[2]
            categories.update(labels)
        speed=np.linalg.norm(velocity[k,ids],axis=1)
        dt_i=config["particles"]["max_displacement_per_radius"]*r/np.maximum(speed,1e-30)
        stiff_counts.append(int((dt_i<config["particles"]["dt_s"]).sum()))
        Re=config["medium"]["density_kg_m3"]*2*radii[pid]*np.linalg.norm(velocity[k,pid]-fluid_velocity[k])/config["medium"]["viscosity_Pa_s"]
        details.append(dict(step=k,start_s=float(row[0]),dt_s=float(row[2]),
            reason="particle_displacement" if particle_limited else "macro_boundary_remainder",
            fastest_particle_id=pid,limiter_particle_id=pid if particle_limited else None,
            categories=labels,branch_assigned=bool(assigned),contact_component_size=int(size),
            near_contact_distance_ratio=float(ratios.min()),wall_distance_over_radius=float(wall[k]/radii[pid]),
            position_m=points[k].tolist(),Re=float(Re),
            primary_force_N=float(np.linalg.norm(forces[k,pid,0])),secondary_force_N=float(np.linalg.norm(forces[k,pid,1])),
            particles_requiring_substeps=stiff_counts[-1]))
    total=sum(counts.values())
    top=counts.most_common()
    prefix=np.cumsum([n for _,n in top])
    top90=int(np.searchsorted(prefix,.9*total)+1) if total else 0
    paircounts=rows[:,6]*(rows[:,6]-1)/2
    weighting={}
    for name,weight in [("accepted_physical_time",rows[:,2]),("measured_step_wall_time",rows[:,11]),
                        ("accepted_state_pair_samples",paircounts)]:
        weighting[name]={key:float(np.average(rows[:,12+i],weights=weight))
                        for i,key in enumerate(["FAR_FIELD","TRANSITION","NEAR_CONTACT","CONTACT"])}
        weighting[name]["any_NEAR_CONTACT_or_CONTACT"]=float(np.average(rows[:,14]+rows[:,15]>0,weights=weight))
    return dict(segment=summary["segment"],source=str(path),accepted_steps=len(rows),
        particle_limited_steps=total,macro_boundary_limited_steps=int((rows[:,4]==0).sum()),
        distinct_limiter_particles=len(counts),particles_covering_90_percent=top90,
        top_limiters=[dict(id=i,steps=n,step_fraction=n/max(total,1),duration_s=durations[i]) for i,n in top],
        category_counts=dict(categories),median_particles_requiring_substeps=float(np.median(stiff_counts)),
        median_active_particles=float(np.median(rows[:,6])),
        median_limiter_contact_component=float(np.median([d["contact_component_size"] for d in details])),
        maximum_limiter_contact_component=max(d["contact_component_size"] for d in details),
        localized_few_particle_limiter=bool(0<top90<=5),validity_weighting=weighting,
        original_evaluation_pair_weighted=summary["physics_validity"]["pair_regime"],details=details)


def main():
    out=new_run("limiter_analysis",dict(method="Offline analysis of every accepted-step start; no diagnostic wall queries injected into solver"))
    selected={}
    for path in sorted((BASE/"profiles").glob("*/SUMMARY.json")):
        summary=json.loads(path.read_text())
        if summary.get("status")=="COMPLETE" and summary.get("backend")=="P0":
            selected[summary["segment"]["id"]]=path.parent
    results=[]
    for name,path in selected.items():
        row=analyze(path)
        details=row.pop("details")
        write(out/(name+"_EVERY_ACCEPTED_STEP.json"),details)
        results.append(row)
    lines=["# Timestep limiter analysis", "", "Global particle IDs are retained. Macro-boundary remainders are counted separately from displacement limits.",
        "Categories overlap; a contact/cluster limiter can also be near-wall and outside the far-field approximation.","",
        "| Segment | Steps | Particle-limited | IDs covering 90% | Median stiff/active | Median contact component |",
        "|---|---:|---:|---:|---:|---:|"]
    for r in results:
        lines.append(f"| {r['segment']['id']} | {r['accepted_steps']} | {r['particle_limited_steps']} | {r['particles_covering_90_percent']} | {r['median_particles_requiring_substeps']}/{r['median_active_particles']} | {r['median_limiter_contact_component']} |")
    lines += ["", "Local-subcycling eligibility is evidence-based: few limiting IDs alone are insufficient if the contact/near-contact buffer covers much of the population. Long-range active–inactive forces must still be evaluated; freezing the remaining bubbles is not allowed."]
    (out/"TIMESTEP_LIMITER_ANALYSIS.md").write_text("\n".join(lines)+"\n")
    valid=["# Performance versus physics validity", "", "These are four short checkpoint replays, not an instrumented wall-time census of the full gold runs.",
           "The original evaluation-pair fractions include RK stages. Accepted-time and measured-step wall-time fractions below use the **start** of each accepted interval.",
           "FAR/TRANSITION/NEAR/CONTACT values are pair fractions averaged with the named weights, not a partition of CPU time between force kernels. The any_NEAR_CONTACT_or_CONTACT value measures time/step coverage by at least one such pair.",
           "Pair fractions and the fraction of steps containing any near/contact pair answer different questions; neither establishes local physical accuracy.", ""]
    for r in results:
        valid.append("## "+r['segment']['id'])
        for name,value in r["validity_weighting"].items():
            valid.append("- "+name+": "+json.dumps(value))
        valid.append("")
    valid.append("Large cost in NEAR_CONTACT/CONTACT reflects the numerical stiffness of advancing an approximate far-field model in close configurations. Smaller timesteps do not by themselves make that force model physically valid. Finite-Re and near-wall limitations remain.")
    (out/"PERFORMANCE_VS_PHYSICS_VALIDITY.md").write_text("\n".join(valid)+"\n")
    finish(out,dict(status="COMPLETE" if len(results)==4 else "PARTIAL",segments=results))


if __name__=="__main__":
    main()
