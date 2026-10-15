"""Actual fixed-ID metrics, inherited layered gates, and control/uncertainty separation."""
import time
import numpy as np
import h5py
from scripts.straight_tube_exp1.f10o.analyze import load, metrics, topology_record
from scripts.straight_tube_exp1.f10p.analyze import comparison, layered
from .common import BASE, Q, CONFIG, policy, read, write, identity


def historical_force_summary():
    data=read(BASE/'diagnostics/HISTORICAL_Y_FROZEN_AUDIT.json')
    def norm(x):return np.linalg.norm(np.asarray(x),axis=-1)
    out={}
    for key,row in data['runs'].items():
        f=row['last_force_sample'];frozen=data['frozen'][key]
        out[key]=dict(last_saved_force_time_s=f['time_s'],last_step_contact_pairs=row['last_step_contact_pair_count'],
                      pressure_end_min_max_Pa=row['step_pressure_end_Pa'],
                      last_saved_force_pN={k:dict(mean_norm=float(norm(f[k]).mean()*1e12),max_norm=float(norm(f[k]).max()*1e12))
                                          for k in ['primary_force_N','secondary_force_N','contact_force_N']},
                      endpoint_F2_increment_pN=dict(mean=float(norm(frozen['coarse']['F2_increment_N']).mean()*1e12),
                                                   max=float(norm(frozen['coarse']['F2_increment_N']).max()*1e12)),
                      same_geometry_cache_difference_pN={k:float(norm(np.array(frozen['coarse'][k])-frozen['fine'][k]).max()*1e12)
                                                         for k in ['primary_N','secondary_N','F2_increment_N']})
    base=data['frozen']['DT1']['coarse']
    geometry={key:{k:float(norm(np.array(base[k])-data['frozen'][key]['coarse'][k]).max()*1e12)
                   for k in ['primary_N','secondary_N','F2_increment_N']} for key in ['DT05','FINE']}
    write(BASE/'analysis/HISTORICAL_Y_FORCE_SUMMARY.json',dict(runs=out,different_final_geometry_same_cache_pN=geometry,
        conclusion='Small direct cache differences at fixed states can seed different contact rearrangements; geometry-dependent F2/secondary/contact responses and dt sensitivity coexist. This frozen-state audit does not identify one exclusive causal mechanism.',
        Y_scope='HOMOGENEOUS_WATER_PROXY_SPECIFIC'))


def main():
    start=time.perf_counter();out=BASE/"analysis";out.mkdir(exist_ok=True)
    cfg=read(CONFIG);p=policy();runs={};paths={};costs=[]
    for cp in sorted(BASE.glob("*/CAMPAIGN.json")):
        costs.append(read(cp))
        for rp in sorted((cp.parent/"runs").glob("*/RESULT.json")):
            tag=rp.parent.name;r=read(rp);paths[tag]=rp.parent
            r['original_source_metadata']=r['source']
            r['source']=read(rp.parent/'INPUT.json')['source']
            t,x=load(rp.parent);m=metrics(t,x);top,labs,ets,eids=topology_record(t,x)
            np.savez_compressed(out/(tag+".npz"),**m,labels=labs,edge_times_s=ets,edge_ids=eids)
            with h5py.File(rp.parent/"trajectory.h5") as h:
                rows=h["step_metrics"][:];extra=h["extra_metrics"][:]
                actual_steps=h["step_time_duration_segment"][:]
            assert np.max(np.linalg.norm(x[:,:,1:],axis=-1))+p["radius_m"]<.0008
            assert np.max(abs(x[:,:,0]))+p["radius_m"]<.01625
            switches=read(rp.parent/"SWITCHES.json")
            assert all(abs(e["ideal_command_time_s"]-e["effective_time_s"])<1e-15 for e in switches)
            assert np.all(abs(actual_steps[1:,0]-np.sum(actual_steps[:-1,:2],axis=1))<1e-13)
            r.update(topology=top,rg_initial_final_m=m["rg_m"][[0,-1]],xy_initial_final=m["directional_xy"][[0,-1]],
                     Q_initial_final_m2=m["tensor_m2"][[0,-1]],aspect_initial_final=m["aspect"][[0,-1]],
                     CPU_solver_s=float(rows[:,7].sum()),CPU_resistance_s=float(rows[:,8].sum()),
                     actual_pressure_min_max_Pa=[float(extra[:,3].min()),float(extra[:,4].max())],
                     accepted_max_overlap_m=float(rows[rows[:,0].astype(bool),2].max(initial=0)),
                     accepted_max_true_complementarity_m=float(rows[rows[:,0].astype(bool),17].max(initial=0)),
                     target_error_m={a:float(np.linalg.norm(m["centroid_m"][-1]-target)) for a,target in cfg["targets_m"].items()},
                     Rg_relative_change=float(m["rg_m"][-1]/m["rg_m"][0]-1),
                     switch_count=len(switches),switch_times_s=[e["effective_time_s"] for e in switches])
            runs[tag]=r
    comps={};gates={}
    for name,r in runs.items():
        if r["requested_s"]<.002-1e-13:
            gates[name]=dict(BASIC_NUMERICS_PASS="PASS" if r["status"]=="PASS" else "FAIL",
                            GROUP_POSITION_NUMERICALLY_ADMITTED="NOT_RUN",GROUP_SHAPE_NUMERICALLY_ADMITTED="NOT_RUN",
                            GROUP_TOPOLOGY_NUMERICALLY_ADMITTED="NOT_RUN",PARTICLE_ID_NUMERICALLY_ADMITTED="NOT_RUN",
                            TARGET_REACHED={a:False for a in ["X","Y"]},scope="Short-window screening only, not 2ms admission")
            continue
        if name.endswith(("_DT05","_FINE")):continue
        tops=[r["topology"]];available={}
        for suffix in ["DT05","FINE"]:
            stem=name[:-4] if name.endswith("_DT1") else name
            other=stem+"_"+suffix
            path=paths.get(other)
            if name in ["X_COMPAT_DT1","OFF_DT1"]:
                path=Q/"b_navigation_v1/runs"/(("NAV_X_" if name=="X_COMPAT_DT1" else "OFF_")+suffix)
            if path is not None and path.exists():
                s,a=comparison(paths[name],path,"Same frozen F10Q SQP and full-precision control; "+suffix)
                comps[name+"__"+suffix]=s;available[suffix]=s
                np.savez_compressed(out/("compare_"+name+"__"+suffix+".npz"),**a)
                tt,xx=load(path);tops.append(topology_record(tt,xx)[0])
        g=layered(r,available.get("DT05"),available.get("FINE"),tops,p)
        g["target_error_m"]=r["target_error_m"]
        g["endpoint_within_tolerance_single_resolution"]={a:bool(r["status"]=="PASS" and err<=cfg["target_tolerance_m"]) for a,err in r["target_error_m"].items()}
        g["TARGET_REACHED"]={a:bool(g["simulation_group_position_eligible"] and err<=cfg["target_tolerance_m"]) for a,err in r["target_error_m"].items()}
        g["navigation_shape_cost_acceptable"]=bool(abs(r["Rg_relative_change"])<=.05 and r["topology"]["max_components"]==1)
        g["Y_scope"]="HOMOGENEOUS_WATER_PROXY_SPECIFIC"
        g["scope"]="Full numerical admission only when dt/cache repeats exist. Unrepeated robust-state/phase samples are single-resolution observations."
        if name in ["X_COMPAT_DT1","OFF_DT1"]:
            g["comparison_scope"]="F10Q retained half-step/fine-cache originals, valid because new full-precision interface replay was independently identical; no duplicated scientific cost"
        gates[name]=g
    phase=read(BASE/"diagnostics/X_PHASE_SENSITIVITY.json")
    measured=[]
    for shift in cfg["phase_full_trajectory_deg"]:
        name=f"X_PHASE_{shift:+d}_DT1"
        if name not in runs:continue
        predictor=next(v for v in phase["rows"] if v["delta_Xplus_deg"]==shift)
        actual=np.asarray(runs[name]["endpoint_m"])
        measured.append(dict(shift_deg=shift,trajectory=name,endpoint_m=actual,target_error_m=runs[name]["target_error_m"]["X"],
                             force_only_predicted_endpoint_m=predictor["predicted_F10Q_scaled_endpoint_um"]*1e-6,
                             observed_within_tolerance=gates[name]["endpoint_within_tolerance_single_resolution"]["X"],
                             validated_target=gates[name]["TARGET_REACHED"]["X"]))
    robustness={k:dict(initial=r["initial"],endpoint_m=r["endpoint_m"],error_m=r["target_error_m"]["X"],
                       relative_Rg_change=r["Rg_relative_change"],basic=r["status"],
                       numerical_repeats="NOT_RUN" if "rotate" in k else "F10Q_REFERENCE_REPEATS")
                for k,r in runs.items() if k=="X_COMPAT_DT1" or "rotate" in k}
    historical=read(Q/"analysis/FINAL_GATE.json")["groups"]["NAV_Y"]
    gate=dict(identity=identity(),ENGINEERING_PASS="PENDING_FINAL_CHECKS",labels=cfg["labels"],training_ready=False,
              Y_scope="HOMOGENEOUS_WATER_PROXY_SPECIFIC",groups=gates,historical_Y=historical,
              new_Y_selected=read(BASE/"FULL_SELECTION.json").get("Y_selected") if (BASE/"FULL_SELECTION.json").exists() else None,
              TARGET_DEFINITIONS=cfg["targets_m"],target_tolerance_m=cfg["target_tolerance_m"],
              HARDWARE_EXECUTED=False,ELECTRICAL_TO_ACOUSTIC_CALIBRATED=False,EXPERIMENTALLY_CALIBRATED="NOT_VERIFIED",
              CONTACT_MODEL_DEPENDENCE="NOT_VERIFIED_FOR_NEW_CONTROLS; prior historical dependence remains",
              long_trap_navigation="NOT_RUN",arbitrary_Z_navigation="NOT_DEMONSTRATED",
              longest_completed_physical_s=max((r["completed_s"] for r in runs.values()),default=0),
              actual_failures={k:r["reasons"] for k,r in runs.items() if r["status"]!="PASS"})
    local=sum(read(BASE/n)["wall_s"] for n in ["diagnostics/EXECUTION.json","diagnostics/PEAK_BRANCH_REVIEW.json","CONTROL_REGISTRATION.json"])
    costs=dict(remote_campaign_wall_s=sum(c["wall_s"] for c in costs),local_field_diagnostics_s=local,
               CPU_solver_s=sum(r["CPU_solver_s"] for r in runs.values()),CPU_resistance_s=sum(r["CPU_resistance_s"] for r in runs.values()),
               integration_s=sum(r["integration_s"] for r in runs.values()),
               peak_RSS_MiB=max((c.get("peak_rss_MiB",0) for c in costs),default=0),GPU_kernel_only_s=None,
               timing_scope="Actual hybrid GPU force/field/walls plus CPU Sequential QP; upload/render/human development excluded from scientific budget",
               trajectory_attempts=len(runs),failed_attempts=sum(r["status"]!="PASS" for r in runs.values()))
    write(out/"RUN_ANALYSIS.json",runs);write(out/"COMPARISONS.json",comps);write(out/"FINAL_GATE.json",gate)
    historical_force_summary()
    write(out/"X_PHASE_RESULTS.json",dict(field=phase,actual_2ms=measured,unrun_phase_shifts=sorted(set(cfg["phase_perturbation_deg"])-set(cfg["phase_full_trajectory_deg"])-{0}),scope="No trajectory inferred for field-only phase samples"))
    write(out/"X_INITIAL_ROBUSTNESS.json",robustness)
    costs["analysis_s"]=time.perf_counter()-start
    costs["total_scientific_s"]=costs["remote_campaign_wall_s"]+local+costs["analysis_s"]
    write(out/"COSTS.json",costs)
    print("Scientific seconds",costs["total_scientific_s"],"trajectories",len(runs))
    print({k:{n:g[n] for n in ["BASIC_NUMERICS_PASS","GROUP_POSITION_NUMERICALLY_ADMITTED","GROUP_SHAPE_NUMERICALLY_ADMITTED","PARTICLE_ID_NUMERICALLY_ADMITTED","TARGET_REACHED"]} for k,g in gates.items() if not k.endswith("SHORT")})


if __name__=="__main__":main()
