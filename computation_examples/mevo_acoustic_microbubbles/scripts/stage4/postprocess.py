"""Read-only numerical summaries and ParaView exports of completed Stage 4 runs."""
import csv
import json
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET
import h5py
import numpy as np
import pyvista as pv
from scripts.stage4.common import ROOT, require, write, sha


def read(path):
    return json.loads((ROOT / path).read_text())


def stage(kind):
    return read(Path("results/stage4") / kind / "STATE.json")


def saved_slip_reynolds(row):
    """Measured saved-state diagnostics; no finite-Re correction is applied."""
    folder = ROOT/row["run"]
    medium = read(folder/"RUN_IDENTITY.json")["config"]["medium"]
    path = folder/"particles/trajectories.h5"
    with h5py.File(path) as h:
        if not h.attrs["complete"]:
            raise ValueError("Reynolds summary requires a closed, complete diagnostic or formal run")
        active = h["active"][:]
        flow = h["flow_velocity_m_s"][:]
        radii = h["radii_m"][:]
        times = h["time_s"][:]
        ids = h["particle_id"][:]
        result = {}
        for label,key in [("constrained","total_velocity_m_s"),("unconstrained","free_total_velocity_m_s")]:
            slip = h[key][:]-flow
            re = 2*medium["density_kg_m3"]*radii[None,:]*np.linalg.norm(slip,axis=-1)/medium["viscosity_Pa_s"]
            good = active & np.isfinite(re)
            if good.any():
                masked = np.where(good,re,-np.inf)
                frame,particle = np.unravel_index(masked.argmax(),masked.shape)
                result[label] = dict(maximum=float(re[frame,particle]),time_s=float(times[frame]),
                    particle_id=int(ids[particle]),saved_active_sample_count=int(good.sum()),
                    sample_count_Re_at_least_1=int(((re>=1)&good).sum()),
                    sample_fraction_Re_at_least_1=float(((re>=1)&good).sum()/good.sum()))
            else:
                result[label] = dict(maximum=None,saved_active_sample_count=0)
    return dict(run=row["run"],trajectory_sha256=sha(path),formula="Re = 2*rho*R0*|v-u_flow|/mu",
                measurements=result,
                scope="Saved active states only; peaks between saved frames can be larger. Re >= 1 violates the small-Re assumption behind Stokes mobility; values below 1 alone do not establish asymptotic validity. Unconstrained values omit geometric reactions and are reported separately.")


def scalar_summary(rows,key):
    good = [(r["time_s"],r[key]) for r in rows if r.get(key) is not None]
    if not good:
        return dict(first_defined=None,last_defined=None,minimum=None,maximum=None,
                    time_average_over_defined_intervals=None,covered_duration_s=0.,defined_samples=0,total_samples=len(rows))
    t,v = np.array(good).T
    integral,covered = 0.,0.
    for before,after in zip(rows[:-1],rows[1:],strict=True):
        if before.get(key) is not None and after.get(key) is not None:
            dt = after["time_s"]-before["time_s"]
            integral += .5*dt*(before[key]+after[key])
            covered += dt
    return dict(first_defined=float(v[0]),last_defined=float(v[-1]),minimum=float(v.min()),maximum=float(v.max()),
                time_average_over_defined_intervals=float(integral/covered) if covered>0 else None,
                covered_duration_s=float(covered),observation_span_s=float(rows[-1]["time_s"]-rows[0]["time_s"]),
                defined_samples=len(v),total_samples=len(rows))


def aggregation(row):
    folder = ROOT / row["run"]
    metrics = read(folder / "metrics.json")
    rows = metrics["rows"]
    summary = read(folder / "particles/SUMMARY.json")
    names = ["mean_nearest_neighbor_m","cluster_count","largest_cluster_fraction","near_contact_fraction",
             "BRANCH_A_roi_count","BRANCH_B_roi_count","BRANCH_A_normalized_occupancy",
             "BRANCH_B_normalized_occupancy","A_to_B_concentration_ratio","B_to_A_concentration_ratio"]
    result = {key:scalar_summary(rows,key) for key in names}
    early = [r for r in rows if r["time_s"]<=.1+1e-9 and r.get("mean_nearest_neighbor_m") is not None]
    result["early_NN_slope_m_s"] = float(np.polyfit([r["time_s"] for r in early],
                               [r["mean_nearest_neighbor_m"] for r in early],1)[0]) if len(early)>=3 else None
    result["early_NN_slope_window_s"] = [0.,.1]
    result["pair_regime"] = summary["physics_validity"]["pair_regime"]
    result["physics_validity_status"] = summary["physics_validity"]["validity_status"]
    result["wall_hit_particle_count"] = summary["wall_hit_particle_count"]
    result["contact_pair_events"] = summary["contact_pair_events"]
    result["saved_slip_reynolds"] = saved_slip_reynolds(row)
    result["branch_spatial_clusters"] = {}
    for branch in ["BRANCH_A","BRANCH_B"]:
        selected = [dict(r[branch+"_spatial_aggregation"],time_s=r["time_s"])
                    for r in rows if branch+"_spatial_aggregation" in r]
        result["branch_spatial_clusters"][branch] = {key:scalar_summary(selected,key)
            for key in ["cluster_count","largest_cluster_fraction","mean_nearest_neighbor_m","near_contact_fraction"]}
    # Direct simultaneous snapshots identify contacts as particles, separate from pair fractions.
    particle_contact_fractions = []
    with h5py.File(folder / "particles/trajectories.h5") as h:
        radii = h["radii_m"][:]
        from scipy.spatial import cKDTree
        tolerance = read(folder / "RUN_IDENTITY.json")["config"]["particles"]["contact_tolerance_m"]
        for x,active in zip(h["positions_m"],h["active"],strict=True):
            q,r = x[active],radii[active]
            pairs = cKDTree(q).query_pairs(2*r.max()+tolerance,output_type="ndarray") if len(q)>1 else np.empty((0,2),int)
            if len(pairs):
                pairs = pairs[np.linalg.norm(q[pairs[:,0]]-q[pairs[:,1]],axis=1)<=r[pairs].sum(axis=1)+tolerance]
            particle_contact_fractions.append(len(np.unique(pairs))/len(q) if len(q) else 0.)
    result["contacting_particle_fraction"] = dict(maximum=float(max(particle_contact_fractions)),
        final=float(particle_contact_fractions[-1]),snapshot_mean=float(np.mean(particle_contact_fractions)),
        definition="Active particles in at least one geometric contact / active particles, at saved states; zero if none active")
    result["scope"] = "Spatial clusters by distance threshold, never fusion. Survivor population changes over time. Conditional concentration ratios undefined at zero denominator; defined-only temporal summary must not hide empty ROIs. Branch clusters are spatial groups, not assigned-outcome clusters."
    return result


def secondary_comparison(out):
    primary = {r["case"]:r for r in require("primary_only")["rows"]}
    secondary = {r["case"]:r for r in require("primary_secondary")["rows"]}
    comparisons = {}
    for label in ["B2","B3"]:
        a,b = primary[label],secondary[label]
        assert a["initial_population_sha256"] == b["initial_population_sha256"]
        pa,pb = aggregation(a),aggregation(b)
        outcome_a,outcome_b = read(Path(a["run"]) / "BRANCH_METRICS.json"),read(Path(b["run"]) / "BRANCH_METRICS.json")
        assert outcome_a["assignment"].keys() == outcome_b["assignment"].keys()
        initial_count = len(outcome_a["assignment"])
        assert initial_count > 0
        transition = {}
        for pid,value in outcome_a["assignment"].items():
            key = value+"->"+outcome_b["assignment"][pid]
            transition[key] = transition.get(key,0)+1
        target = "f_A" if label=="B2" else "f_B"
        target_count = "N_A" if label=="B2" else "N_B"
        counts = []
        for row in [a,b]:
            counts.append(dict(N_initial=initial_count,
                **{k:row[k] for k in ["N_A","N_B","N_unassigned","N_wall_loss",target]},
                target_fraction_of_initial_population=row[target_count]/initial_count))
        comparisons[label] = dict(primary_run=a["run"],secondary_run=b["run"],
            primary=pa,secondary=pb,primary_counts=counts[0],secondary_counts=counts[1],
            target_fraction_difference=b[target]-a[target] if a[target] is not None and b[target] is not None else None,
            target_arrival_fraction_difference=(b[target_count]-a[target_count])/initial_count,
            changed_assignments=sum(v for k,v in transition.items() if k.split("->")[0]!=k.split("->")[1]),
            transitions=transition)
    result = dict(status="PASS_WITH_SCOPE",cases=comparisons,seed=primary["B2"]["seed"],
                  scope="Secondary comparison on one paired seed per steering direction; descriptive, no multi-seed secondary significance claim. Near-contact force remains approximate.")
    write(out / "SECONDARY_COMPARISON.json",result)
    return result


def performance(out):
    rows = require("statistics/paired_runs")["rows"]+require("primary_secondary")["rows"]
    data = [dict(case=r["case"],seed=r["seed"],secondary=r["secondary"],run=r["run"],
                 elapsed_s=r["elapsed_s"],**r["timings"]) for r in rows]
    with (out / "PERFORMANCE.csv").open("w",newline="") as stream:
        writer = csv.DictWriter(stream,fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)
    summary = {label:dict(median_elapsed_s=float(np.median([r["elapsed_s"] for r in rows if r["case"]==label and not r["secondary"]])),
                         maximum_elapsed_s=float(max(r["elapsed_s"] for r in rows if r["case"]==label and not r["secondary"])))
               for label in ["B0","B1","B2","B3"]}
    diagnostics = {}
    for name in ["profile_B2_secondary", "certified_wall_trial", "component_contact_shadow",
                 "secondary_contact_snapshot", "later_geometry_certificate", "later_particle_geometry_certificate",
                 "certificate_prefix_check", "secondary_numerical_diagnostic"]:
        path = ROOT/"results/stage4/performance"/name/"STATE.json"
        if path.exists():
            saved = read(path)
            evidence = str(Path(saved["run"])/"SUMMARY.json")
            diagnostics[name] = dict(status=saved["status"],summary=evidence,sha256=sha(ROOT/evidence))
    write(out / "PERFORMANCE.json",dict(cases=summary,rows=data,diagnostics=diagnostics,
          timing_scope="Subtimers overlap: integration contains force/flow/contact evaluations. Do not sum all fields as disjoint time components. Concurrent runs include CPU contention."))
    return summary


def paraview(out):
    destination = out / "paraview"
    destination.mkdir()
    geometry,flow = require("geometry"),require("simvascular/steady")
    for source,name in [(geometry["surface"],"vessel_branch_labels.vtp"),
                        (flow["snapshot"],"steady_flow.vtu"),
                        (str(Path(geometry["run"])/"centerlines.vtp"),"centerlines.vtp")]:
        shutil.copy2(ROOT / source,destination / name)
    cfg = read(Path(require("primary_only")["rows"][0]["run"]) / "RUN_IDENTITY.json")["config"]
    pzts = []
    for i,pzt in enumerate(cfg["pzt"]):
        disc = pv.Disc(center=pzt["position_m"],normal=pzt["normal"],inner=0,outer=pzt["size_m"][0]/2,c_res=96)
        disc.cell_data["PZT_ID"] = np.full(disc.n_cells,i)
        pzts.append(disc)
    pv.merge(pzts).save(destination / "pzt_array.vtp")
    planes = read(Path(geometry["run"]) / "PLANES.json")
    for name,plane in planes.items():
        if "radius_m" in plane:
            pv.Disc(center=plane["origin_m"],normal=plane["normal"],inner=0,
                    outer=plane["radius_m"],c_res=96).save(destination / (name.lower()+"_plane.vtp"))
    mapping = {"PARENT":0,"BRANCH_A":1,"BRANCH_B":2,"UNASSIGNED":-1,"WALL_LOSS":-2}
    collections = []
    for row in require("primary_only")["rows"]+require("primary_secondary")["rows"]:
        case = row["case"]+("_secondary" if row["secondary"] else "_primary")
        folder = destination / case
        folder.mkdir()
        branch = read(Path(row["run"]) / "BRANCH_METRICS.json")
        root = ET.Element("VTKFile",type="Collection",version="0.1",byte_order="LittleEndian")
        collection = ET.SubElement(root,"Collection")
        with h5py.File(ROOT / row["run"] / "particles/trajectories.h5") as h:
            assert h.attrs["complete"]
            ids,radii = h["particle_id"][:],h["radii_m"][:]
            assignment = np.array([mapping[branch["assignment"][str(pid)]] for pid in ids],np.int32)
            for k,t in enumerate(h["time_s"]):
                active = h["active"][k]
                points = pv.PolyData(h["positions_m"][k][active])
                points["particle_id"] = ids[active]
                points["radius_m"] = radii[active]
                points["outcome_by_end_of_observation"] = assignment[active]
                points["current_region"] = np.where(points.points[:,0]<0,0,np.where(points.points[:,1]>=0,1,2)).astype(np.int32)
                for name in h:
                    if name in ["particle_id","radii_m","positions_m","active","time_s","accepted_dt_s"] or h[name].shape[:2]!=(len(h["time_s"]),len(ids)):
                        continue
                    values = h[name][k][active]
                    if np.iscomplexobj(values):
                        points[name+"_real"] = values.real
                        points[name+"_imag"] = values.imag
                    else:
                        points[name] = values
                name = f"particles_{k:04d}.vtp"
                points.save(folder / name)
                ET.SubElement(collection,"DataSet",timestep=str(float(t)),group="",part="0",file=case+"/"+name)
        path = destination / (case+".pvd")
        ET.ElementTree(root).write(path,encoding="utf-8",xml_declaration=True)
        collections.append(str(path.relative_to(ROOT)))
    write(destination / "EXPORT_MANIFEST.json",dict(units="metres, seconds, Pa, newtons",branch_ids=mapping,
        collections=collections,source_runs=require("primary_only")["rows"]+require("primary_secondary")["rows"],
        semantics="Only active particles exported per frame; first-exit outcomes stored separately from instantaneous geometric region. Points/glyphs may be enlarged for visibility; radius_m is physical. No interpolation or smoothing of scientific positions.",
        files={str(p.relative_to(destination)):sha(p) for p in destination.rglob("*") if p.is_file()}))
    return collections
