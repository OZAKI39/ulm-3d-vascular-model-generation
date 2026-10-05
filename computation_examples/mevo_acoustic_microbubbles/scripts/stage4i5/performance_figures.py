"""Six evidence-only performance figures; data values are never edited for display."""
import json
from pathlib import Path
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scripts.stage4.figures import style,save
from scripts.stage4i5.common import ROOT,BASE,GOLD,write


def completed_replays():
    found={}
    for file in sorted((BASE/"p1_vtk_contact").glob("*/SUMMARY.json")):
        row=json.loads(file.read_text())
        if row.get("status")=="COMPLETE" and not row.get("cprofile"):
            found[(row["backend"],row["segment"]["id"])]=row
    return found


def cumulative_figure(out):
    fig,axes=plt.subplots(1,2,figsize=(15,6),layout="constrained")
    for ax,label in zip(axes,["B2","B3"]):
        branch=json.loads((GOLD/"cases"/label/"BRANCH_METRICS.json").read_text())
        events=sorted(branch["crossings"],key=lambda e:e["crossing_time_s"])
        t=np.r_[0.,[e["crossing_time_s"] for e in events],6.]
        a=np.r_[0,np.cumsum([e["branch"]=="BRANCH_A" for e in events]),branch["N_A"]]
        b=np.r_[0,np.cumsum([e["branch"]=="BRANCH_B" for e in events]),branch["N_B"]]
        for values,name,color in [(a,"Branch A","#59dbca"),(b,"Branch B","#54a8ff"),(256-a-b,"Unassigned","#eeb770")]:
            ax.step(t,values,where="post",label=name,color=color,lw=2.5)
        ax.set(xlim=(0,6),ylim=(0,265),xlabel="Simulated time [s]",ylabel="Particles",title=f"{label}: A {branch['N_A']} / B {branch['N_B']} / U {branch['N_unassigned']}")
        ax.grid(alpha=.3)
        ax.legend()
    fig.suptitle("P6  |  Cumulative first branch assignments — immutable gold")
    return save(fig,out,"figure_P6_cumulative_assignments","Exact saved first-crossing event times; 256 initial particles per case. A plateau while particles remain unassigned does not authorize early termination or establish their eventual branch.")


def make_performance_figures(out,decision):
    style()
    out=Path(out)/"performance_figures"
    out.mkdir()
    figures=[]
    runs=completed_replays()
    segments=["B2_S1","B2_S2","B3_S1","B3_S2"]
    chosen=decision["recommended_backend_for_JK"]
    components={"Wall queries/batch":["wall_closest_distance","wall_normal_including_repeated_closest_point","wall_cache_and_batch","wall_combined_closest_and_normal"],
        "Contact position":["contact_position_projection"],"NNLS / velocity":["NNLS","contact_velocity_assembly"],
        "Candidate search":["contact_candidate_search"],"Wall constraint":["wall_projection","wall_velocity_constraint"],
        "Response / pair / flow":["response_lookup","pair_force","flow_interpolation"],
        "Validity / other":["physics_validity","acoustic_interpolation","branch_classification","timestep_computation"]}
    fig,axes=plt.subplots(2,2,figsize=(16,11),layout="constrained")
    colors=plt.get_cmap("tab20")(np.linspace(0,1,len(components)))
    for ax,segment in zip(axes.flat,segments):
        data=[runs[(backend,segment)] for backend in ("P0",chosen)]
        bottom=np.zeros(2)
        for color,(label,keys) in zip(colors,components.items()):
            values=np.array([sum(r["timers_exclusive_s"].get(k,0.) for k in keys) for r in data])
            ax.bar([0,1],values,bottom=bottom,label=label,color=color)
            bottom+=values
        ax.set(xticks=[0,1],xticklabels=["P0",chosen],ylabel="Exclusive timed seconds",title=segment)
    handles,labels=axes.flat[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="outside lower center",ncol=4)
    fig.suptitle("P1  |  Matched short-replay runtime components")
    figures.append(save(fig,out,"figure_P1_runtime_breakdown","Matched explicit timers, without cProfile. Components exclude initialization, HDF output and diagnostic recording. The recorded residual Python/timer overhead is reported in JSON, so bars are not claimed to exhaust total runtime."))

    fig,axes=plt.subplots(1,2,figsize=(15,6),layout="constrained")
    for ax,label in zip(axes,["B2","B3"]):
        for segment,color in [("S1","#59dbca"),("S2","#54a8ff")]:
            with h5py.File(ROOT/runs[("P0",label+"_"+segment)]["run"] / "replay.h5") as h:
                dt=h["rows"][:,2]*1e6
                reason=h["rows"][:,4]
            bins=np.geomspace(max(dt.min(),1e-7),max(dt.max(),1e-6),45)
            ax.hist(dt,bins=bins,histtype="step",lw=2,label=f"{segment}; macro remainders {int((reason==0).sum())}",color=color)
        ax.set(xscale="log",xlabel="Accepted dt [µs]",ylabel="Step count",title=label)
        ax.legend()
        ax.grid(alpha=.3)
    fig.suptitle("P2  |  Effective timestep distribution in complete short replays")
    figures.append(save(fig,out,"figure_P2_effective_dt","All accepted steps of the selected segments. Macro-boundary remainders are included and explicitly counted; the smallest dt is not automatically a stiffness measurement. Not a full gold timestep histogram."))

    limiter=json.loads((BASE/"limiter_analysis/STATE.json").read_text())
    fig,axes=plt.subplots(4,2,figsize=(17,18),layout="constrained")
    for (ax,category_ax),row in zip(axes,limiter["segments"]):
        top=row["top_limiters"][:12]
        ax.bar(np.arange(len(top)),[100*r["step_fraction"] for r in top],color="#54a8ff")
        ax.set(xticks=np.arange(len(top)),xticklabels=[str(r["id"]) for r in top],xlabel="Global particle ID (top 12)",ylabel="Particle-limited steps [%]",
            title=f"{row['segment']['id']}  |  90% needs {row['particles_covering_90_percent']} IDs")
        ax.set_ylim(0,max(100*r["step_fraction"] for r in top)*1.35)
        ax.text(.98,.95,f"Median stiff: {row['median_particles_requiring_substeps']:.0f}/{row['median_active_particles']:.0f}\nMedian contact component: {row['median_limiter_contact_component']:.0f}",ha="right",va="top",transform=ax.transAxes,fontsize=12)
        categories=row["category_counts"]
        names=["free","near_wall","contact","cluster","near_contact_validity_region","branch_assigned"]
        values=[100*categories.get(name,0)/row["particle_limited_steps"] for name in names]
        labels=[name.replace("_"," ").replace("validity region","region") for name in names]
        category_ax.barh(labels,values,color="#59dbca")
        category_ax.set(xlabel="Particle-limited steps [%]; categories overlap",xlim=(0,105),title=row["segment"]["id"]+" categories")
    fig.suptitle("P3  |  Sources of the global timestep limit")
    figures.append(save(fig,out,"figure_P3_limiter_particles","Particle-limit events only; macro-boundary remainders excluded from the denominator. Category percentages overlap and must not be added. Global IDs, not compact active-array indices; every-step details remain in TIMESTEP_LIMITER_ANALYSIS."))

    wall={}
    gold_vtk=json.loads((GOLD/"GOLD_BASELINE_MANIFEST.json").read_text())["environment_observed_at_sealing"]["packages"]["vtk"]
    for kind in ["p1_vtk_contact","p2_warp"]:
        for path in sorted((BASE/kind).glob("*/SUMMARY.json")):
            row=json.loads(path.read_text())
            if row.get("point_count")==10000 and row.get("VTK_version")==gold_vtk:
                wall[row["backend"]]=row
    fig,ax=plt.subplots(figsize=(11,6),layout="constrained")
    names,seconds=[],[]
    for key,row in wall.items():
        if not names:
            names.append("Original VTK")
            seconds.append(row["reference_query_s"])
        gate=decision.get("independent_prerequisite_gates",{}).get(key,{})
        name="Static VTK" if key=="vtk" else "Warp + host transfers"
        if gate.get("pass_geometry") is False:
            name+="\nGeometry gate: FAIL"
        names.append(name)
        seconds.append(row["candidate_query_s"])
    bars=ax.bar(names,np.array(seconds)*100,color=["#8ba6b8","#59dbca","#eeb770"][:len(names)])
    ax.bar_label(bars,fmt="%.3g",padding=5)
    ax.set_ylim(0,max(seconds)*110)
    ax.set(ylabel="Query time [µs / point]",title="P4  |  Same 10,000-point wall fixture on the server")
    figures.append(save(fig,out,"figure_P4_wall_backend","Build/compile excluded after warmup; Warp host adapter includes actual transfers. Query speed alone does not imply projection or trajectory equivalence. Geometry decision errors and adoption gates are reported separately."))

    comparison_reports={}
    for path in sorted((BASE/"full_regression").glob("*/SUMMARY.json")):
        summary=json.loads(path.read_text())
        comparison_reports[summary["backend"]]=summary
    comparisons=[row for summary in comparison_reports.values() for row in summary.get("rows",[])]
    if not comparisons:
        raise RuntimeError("No completed candidate comparisons: refusing an empty speedup/error figure")
    fig,ax=plt.subplots(figsize=(12,7),layout="constrained")
    palette=dict(zip(sorted({r["backend"] for r in comparisons}),plt.get_cmap("tab10").colors))
    seen=set()
    for row in comparisons:
        label=row["backend"]
        ax.scatter(row["speedup"],row["rms_position_error_m"],color=palette[label],s=65,clip_on=False,
            marker="o" if row["strict_P1_gate"]=="PASS" else "x",label=label if label not in seen else None)
        seen.add(label)
    ax.axvline(1.,color="#8ba6b8",ls="--")
    ax.set(yscale="symlog",xlabel="Component-timed short-replay speedup vs P0",ylabel="Position RMS difference [m]",title="P5  |  Speedup and trajectory error")
    ax.set_yscale("symlog",linthresh=1e-12)
    ax.set_ylim(0,max(1e-12,1.5*max(r["rms_position_error_m"] for r in comparisons)))
    if all(r["rms_position_error_m"]==0 for r in comparisons):
        ax.text(.5,.5,"All displayed RMS differences are exactly zero",ha="center",transform=ax.transAxes)
    warp=decision.get("independent_prerequisite_gates",{}).get("warp",{})
    if warp.get("pass_geometry") is False:
        ax.text(.5,.43,"Warp: geometry gate failed; failed dense replays are excluded",ha="center",
            transform=ax.transAxes,fontsize=11,color="#eeb770")
    ax.grid(alpha=.3)
    ax.legend()
    figures.append(save(fig,out,"figure_P5_speedup_error","Component-timed runs only; the separate minimal-timer adoption confirmation is reported in the backend decision and report. Zero errors stay exactly zero on a symmetric-log axis. Circles pass the segment trajectory/contact gate, crosses fail it; a segment pass does not imply backend adoption. Each point is one segment; failed/incomplete trajectories are excluded from speedup claims and retained in failure records."))
    figures.append(cumulative_figure(out))
    write(out/"FIGURE_MANIFEST.json",dict(figures=figures,backend_decision=decision["run"],
        comparison_reports={name:row["run"] for name,row in comparison_reports.items()}))
    return figures
