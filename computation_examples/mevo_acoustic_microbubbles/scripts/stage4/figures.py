"""Eight evidence figures from completed data; no changes to simulation fields."""
from pathlib import Path
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import pyvista as pv
from scripts.stage4.common import ROOT, require, write
from scripts.stage4.postprocess import read, stage
from mevo_mb.acoustics.basis import ComplexAcousticField

COLORS = dict(B0="#8ba6b8",B1="#d2c2ee",B2="#59dbca",B3="#54a8ff",UNASSIGNED="#eeb770")


def style():
    plt.rcParams.update({"figure.facecolor":"#05090f","axes.facecolor":"#05090f",
        "savefig.facecolor":"#05090f","text.color":"#ecf3ff","axes.labelcolor":"#c5d7e8",
        "axes.edgecolor":"#4b6478","xtick.color":"#9aafc2","ytick.color":"#9aafc2",
        "grid.color":"#344654","font.size":13,"axes.titlesize":17,"axes.labelsize":14,
        "legend.frameon":False,"pdf.fonttype":42,"axes.spines.top":False,"axes.spines.right":False})


def save(fig,out,name,caption):
    fig.savefig(out / (name+".png"),dpi=180,bbox_inches="tight")
    fig.savefig(out / (name+".pdf"),bbox_inches="tight")
    plt.close(fig)
    return dict(name=name,png=str((out/(name+".png")).relative_to(ROOT)),caption=caption)


def window(title):
    p = pv.Plotter(off_screen=True,window_size=(1920,1080))
    p.set_background("#05090f")
    p.add_text(title,position="upper_left",font_size=22,color="#ecf3ff")
    return p


def finish3d(p,path):
    p.camera_position=[(24,-32,24),(-1,0,0),(0,0,1)]
    p.enable_parallel_projection()
    p.show_bounds(xtitle="x [mm]",ytitle="y [mm]",ztitle="z [mm]",grid="back",location="outer",
                  color="#9aafc2",font_size=17,all_edges=False,ticks="outside",n_xlabels=5,n_ylabels=5,n_zlabels=3)
    p.camera.zoom(1.12)
    p.screenshot(path)
    p.close()


def make_figures(out,secondary):
    style()
    out = Path(out) / "figures"
    out.mkdir()
    g,cfd = require("geometry"),require("simvascular/steady")
    vessel = pv.read(ROOT / g["surface"])
    vessel_mm = vessel.scale(1000,inplace=False)
    primary = {r["case"]:r for r in require("primary_only")["rows"]}
    config = read(Path(primary["B0"]["run"]) / "RUN_IDENTITY.json")["config"]
    manifest = []
    p = window("01  |  Y-bifurcation and four ideal acoustic elements")
    p.add_mesh(vessel_mm,color="#a9cddd",opacity=.4,smooth_shading=True)
    for i,element in enumerate(config["pzt"]):
        pos = np.array(element["position_m"])*1000
        disc = pv.Disc(center=pos,normal=element["normal"],inner=0,outer=element["size_m"][0]*500,c_res=96)
        p.add_mesh(disc,color="#59dbca",opacity=.7)
        p.add_point_labels(pos[None], [element["id"]],font_size=22,text_color="white",shape=None,
                           show_points=False,always_visible=True)
    p.add_point_labels(np.array([[10.2,6,0],[10.2,-6,0],[-14,0,0]]),["Branch A","Branch B","Parent inlet"],
                       font_size=24,text_color="white",shape=None,show_points=False,always_visible=True)
    name = "figure_01_geometry_pzt"
    finish3d(p,out/(name+".png"))
    manifest.append(dict(name=name,png=str((out/(name+".png")).relative_to(ROOT)),
        caption="Actual labelled lumen; four ideal 6 mm acoustic discs. Display coordinates in mm, source data SI. No electromechanical PZT model."))

    mesh = pv.read(ROOT / cfd["snapshot"])
    mesh["speed_mm_s"] = np.linalg.norm(mesh["Velocity"],axis=1)*1000
    indices = np.arange(160)
    theta = indices*np.pi*(3-np.sqrt(5))
    r = 1.15e-3*np.sqrt((indices+.5)/len(indices))
    seeds = pv.PolyData(np.column_stack([np.full(len(indices),-.0135),r*np.cos(theta),r*np.sin(theta)]))
    lines = mesh.streamlines_from_source(seeds,vectors="Velocity",integration_direction="forward",
              initial_step_length=.1,min_step_length=.01,max_step_length=.25,max_steps=10000,max_length=.06,
              terminal_speed=1e-8,compute_vorticity=False)
    if lines.n_points < 100:
        raise RuntimeError("Too few valid CFD streamlines for the requested figure")
    lines["speed_mm_s"] = np.linalg.norm(lines["Velocity"],axis=1)*1000
    p = window("02  |  SimVascular steady velocity and streamlines")
    p.add_mesh(vessel_mm,color="#b2cbda",opacity=.10)
    section = mesh.slice(normal=(0,0,1)).scale(1000,inplace=False)
    p.add_mesh(section,scalars="speed_mm_s",cmap="viridis",clim=(0,mesh["speed_mm_s"].max()),
               scalar_bar_args=dict(title="Speed [mm/s]",title_font_size=24,label_font_size=20,color="white"))
    p.add_mesh(lines.scale(1000,inplace=False),scalars="speed_mm_s",cmap="viridis",line_width=2,
               show_scalar_bar=False,clim=(0,mesh["speed_mm_s"].max()))
    name = "figure_02_simvascular_flow"
    finish3d(p,out/(name+".png"))
    manifest.append(dict(name=name,png=str((out/(name+".png")).relative_to(ROOT)),
        caption="Native steady CFD, z=0 speed slice and illustrative VTK streamlines. Streamline rendering does not validate local CFD accuracy."))

    control = require("phase_screening")["cases"]["B2"]
    field = ComplexAcousticField.load(ROOT / control["field"])
    z = int(np.argmin(abs(field.axes_m[2])))
    pressure = field.pressure_Pa[:,:,z]
    extent = [field.axes_m[0][0]*1000,field.axes_m[0][-1]*1000,field.axes_m[1][0]*1000,field.axes_m[1][-1]*1000]
    outline = vessel.slice(normal=(0,0,1))
    outline_mm = outline.points*1000
    segments = []
    packed = outline.lines
    offset = 0
    while offset < len(packed):
        n = packed[offset]
        segments.append(outline_mm[packed[offset+1:offset+1+n],:2])
        offset += n+1
    fig,axes = plt.subplots(1,2,figsize=(15,6),layout="constrained")
    for ax,values,title,cmap,limits in zip(axes,[abs(pressure)/1000,np.angle(pressure)],
            ["Amplitude [kPa]","Phase [rad]"],["magma","twilight"],[(0,40),(-np.pi,np.pi)],strict=True):
        image = ax.imshow(values.T,extent=extent,origin="lower",aspect="equal",cmap=cmap,vmin=limits[0],vmax=limits[1])
        ax.add_collection(LineCollection(segments,colors="#cfefff",linewidths=.9))
        ax.set(xlabel="x [mm]",ylabel="y [mm]",title=title)
        fig.colorbar(image,ax=ax,shrink=.7)
    fig.suptitle("03  |  B2 complex acoustic field at z = 0")
    manifest.append(save(fig,out,"figure_03_acoustic_amplitude_phase","Actual B2 field, fixed phase/amplitude pattern; constant scale across the ROI. Lumen outline overlaid, no pressure clipping in simulation."))

    fig,axes = plt.subplots(2,2,figsize=(15,11),subplot_kw=dict(projection="3d"),layout="constrained")
    for ax,(label,row) in zip(axes.flat,primary.items(),strict=True):
        ax.set_facecolor("#05090f")
        with h5py.File(ROOT / row["run"] / "particles/trajectories.h5") as h:
            positions,active = h["positions_m"][:]*1000,h["active"][:]
            for pid in range(positions.shape[1]):
                stop = np.flatnonzero(~active[:,pid])
                path = positions[:int(stop[0])+1 if len(stop) else len(positions),pid]
                ax.plot(*path.T,color=COLORS[label],alpha=.35,lw=.7)
        for direction in [(1,np.tan(np.pi/6),0),(1,-np.tan(np.pi/6),0)]:
            path = np.array([direction])*np.linspace(0,9,100)[:,None]
            ax.plot(*path.T,color="#7895a7",alpha=.5,lw=1)
        ax.set(xlim=(-2.5,9.5),ylim=(-6,6),zlim=(-1.4,1.4),xlabel="x [mm]",ylabel="y [mm]",zlabel="z [mm]")
        ax.set_title(f"{label}  |  A {row['N_A']} / B {row['N_B']} / U {row['N_unassigned']}")
        ax.set_box_aspect((12,12,2.8))
        ax.view_init(28,-66)
        for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
            axis.set_pane_color((.02,.04,.06,1.))
            axis._axinfo["grid"].update(color=(.3,.4,.5,.3))
    fig.suptitle("04  |  Primary-only trajectories, identical 256 initial bubbles")
    manifest.append(save(fig,out,"figure_04_primary_trajectories","All first-seed trajectories, saved at 20 ms cadence, stopped at immutable first crossing. U means no classified exit by 6 s, not numerical failure."))

    statistics = require("statistics")
    fig,axes = plt.subplots(2,2,figsize=(15,11),layout="constrained")
    for i,label in enumerate(primary):
        for ax,key in zip(axes[0],["f_A","unassigned_fraction"],strict=True):
            metric = statistics["cases"][label][key]
            if metric["mean"] is None:
                ax.text(i,.5,"Undefined",ha="center",rotation=90)
                continue
            low,high = metric["bootstrap_95CI"]
            ax.errorbar(i,metric["mean"],yerr=[[metric["mean"]-low],[high-metric["mean"]]],fmt="o",color=COLORS[label],capsize=7,ms=10)
            ax.scatter(np.full(len(metric["values"]),i)+np.linspace(-.10,.10,len(metric["values"])),metric["values"],color=COLORS[label],alpha=.5,s=15)
    for ax in axes[0]:
        ax.set_xticks(range(4),list(primary))
        ax.set_ylim(-.02,1.02)
        ax.grid(axis="y",alpha=.35)
    axes[0,0].axhline(.5,color="#84949e",ls="--",lw=1)
    axes[0,0].set_title("A fraction among assigned exits")
    axes[0,0].set_ylabel("N_A / (N_A + N_B)")
    axes[0,1].set_title("Unassigned fraction of initial population")
    axes[0,1].set_ylabel("N_unassigned / 256")
    for offset,key,color,title in [(-.14,"A_fraction_of_initial",COLORS["B2"],"Arrived at A"),
                                    (.14,"B_fraction_of_initial",COLORS["B3"],"Arrived at B")]:
        for i,label in enumerate(primary):
            metric = statistics["cases"][label][key]
            low,high = metric["bootstrap_95CI"]
            axes[1,0].errorbar(i+offset,metric["mean"],yerr=[[metric["mean"]-low],[high-metric["mean"]]],
                fmt="o",color=color,capsize=6,ms=8,label=title if i==0 else None)
    axes[1,0].set_xticks(range(4),list(primary))
    axes[1,0].set(ylim=(-.02,1.02),ylabel="Fraction of initial population",title="Actual arrivals by the observation cutoff")
    axes[1,0].legend(fontsize=11)
    axes[1,0].grid(axis="y",alpha=.35)
    clock = np.linspace(0,6,301)
    for label in primary:
        curves = []
        for row in require("statistics/paired_runs")["rows"]:
            if row["case"]!=label:
                continue
            events = read(Path(row["run"])/"BRANCH_METRICS.json")["crossings"]
            times = np.sort([e["crossing_time_s"] for e in events if e["branch"] in ["BRANCH_A","BRANCH_B"]])
            total = sum(row[key] for key in ["N_A","N_B","N_unassigned","N_wall_loss"])
            curves.append(np.searchsorted(times,clock,side="right")/total)
        axes[1,1].plot(clock,np.mean(curves,axis=0),color=COLORS[label],label=label,lw=2)
    axes[1,1].set(xlabel="Time [s]",ylabel="Cumulative assigned exits / initial population",ylim=(-.02,1.02),
                  title="Arrival-time distribution, seed mean")
    axes[1,1].grid(alpha=.35)
    axes[1,1].legend(fontsize=11)
    fig.suptitle("05  |  Ten paired seeds; seed-block bootstrap 95% CI")
    manifest.append(save(fig,out,"figure_05_branch_fractions_ci","Top small dots: individual seed outcomes; error bars: seed-block bootstrap of means, 20,000 draws. Conditional exit preference, censoring and actual A/B arrival fractions are distinct. Cumulative curves are mean observed exits divided by initial count, retaining all unassigned particles in the denominator; no imputed exits or censored-particle deletion."))

    fig,axes = plt.subplots(1,2,figsize=(13,6),layout="constrained")
    for ax,label in zip(axes,["B2","B3"],strict=True):
        branch = read(Path(primary[label]["run"]) / "BRANCH_METRICS.json")
        events = branch["bifurcation_crossings"]
        for outcome,color in [("BRANCH_A",COLORS["B2"]),("BRANCH_B",COLORS["B3"]),("UNASSIGNED",COLORS["UNASSIGNED"])]:
            points = np.array([e["position_m"] for e in events if branch["assignment"][str(e["particle_id"]) ]==outcome])
            if len(points):
                ax.scatter(points[:,1]*1000,points[:,2]*1000,s=23,c=color,alpha=.8,label=outcome.replace("BRANCH_",""))
        ax.set(xlabel="y at first x = 0 crossing [mm]",ylabel="z [mm]",title=label,xlim=(-1.5,1.5),ylim=(-1.5,1.5),aspect="equal")
        ax.grid(alpha=.3)
        ax.legend(loc="upper right",fontsize=10)
    fig.suptitle("06  |  First bifurcation crossings, colored by observed outcome")
    manifest.append(save(fig,out,"figure_06_cross_section","Interpolated first outward x=0 events; later exit assignment shown only as a label. Particles that never reached x=0 are absent; no final-coordinate branch guesses."))

    fig,axes = plt.subplots(2,2,figsize=(14,9),layout="constrained")
    for label in ["B2","B3"]:
        comparison = secondary["cases"][label]
        for mode,linestyle in [("primary","-"),("secondary","--")]:
            rows = read(Path(comparison[mode+"_run"]) / "metrics.json")["rows"]
            for ax,key,scale,title in zip(axes.flat,["mean_nearest_neighbor_m","largest_cluster_fraction","near_contact_fraction","BRANCH_A_normalized_occupancy"],
                    [1e6,1,1,1],["Mean nearest neighbor [um]","Largest cluster / active population","Near-contact particles / active population","Branch A proximal ROI / initial population"],strict=True):
                ax.plot([r["time_s"] for r in rows],[np.nan if r.get(key) is None else r[key]*scale for r in rows],
                        color=COLORS[label],ls=linestyle,label=label+(" primary" if mode=="primary" else " + secondary"),lw=1.5)
                ax.set(xlabel="Time [s]",title=title)
                ax.grid(alpha=.3)
    axes[0,0].legend(fontsize=10)
    fig.suptitle("07  |  Secondary ON/OFF on the same seed")
    manifest.append(save(fig,out,"figure_07_secondary_effect","One paired seed, descriptive secondary comparison. Active population changes on exit; geometric clusters and contacts are not fusion. Complete regime/contact/branch-specific tables in SECONDARY_COMPARISON.json."))

    history = stage("history_branch")
    fig,axes = plt.subplots(1,2,figsize=(15,7),layout="constrained")
    labels = [str(r["particle_id"]) for r in history["rows"]]
    mapping = {"BRANCH_A":0,"BRANCH_B":1,"UNASSIGNED":2,"WALL_LOSS":3}
    result = np.array([[mapping[r["production_branch"]],mapping[r["history_reference_branch"]]] for r in history["rows"]])
    from matplotlib.colors import ListedColormap
    image = axes[0].imshow(result,cmap=ListedColormap([COLORS["B2"],COLORS["B3"],COLORS["UNASSIGNED"],"#f57878"]),vmin=-.5,vmax=3.5,aspect="auto")
    axes[0].set_xticks([0,1],["Production","History reference"])
    axes[0].set_yticks(range(len(labels)),labels)
    axes[0].set_ylabel("Automatically selected particle ID")
    bar = fig.colorbar(image,ax=axes[0],ticks=[0,1,2,3],fraction=.04)
    bar.ax.set_yticklabels(["A","B","U","W"])
    comparisons = [history["branch_agreement_rate"],history["reference_dt_assignment_agreement"],history["reduced_population_control_agreement"],history["B_vs_C_history_only_branch_agreement"]]
    axes[1].barh(range(4),comparisons,color=["#59dbca","#54a8ff","#d2c2ee","#eeb770"])
    axes[1].set_yticks(range(4),["Production vs history","Reference dt refinement","Reduced-population control","Reference inertia vs history"])
    axes[1].set(xlim=(0,1.12),xlabel="Branch agreement fraction")
    for i,value in enumerate(comparisons):
        axes[1].text(value+.02,i,f"{value:.1%}",va="center")
    axes[1].grid(axis="x",alpha=.3)
    fig.suptitle("08  |  Critical-particle robustness: adiabatic two-scale reference")
    manifest.append(save(fig,out,"figure_08_history_branch","Validation-only selected-sample reference, not a fully microcycle-resolved moving history solution. A/B are observed outlet outcomes, U is unassigned at cutoff, and W is wall loss; an A/B-to-U difference is not an observed opposite-branch switch. Production comparison includes inertia, averaging and different interpolation errors; B-vs-C tests history within the same reference construction. Wall hydrodynamics remain unvalidated."))
    write(out / "FIGURE_MANIFEST.json",manifest)
    return manifest
