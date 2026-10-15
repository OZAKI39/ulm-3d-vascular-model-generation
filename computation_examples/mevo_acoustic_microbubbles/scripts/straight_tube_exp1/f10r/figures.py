"""English data-only F10R research figures; no inferred trajectories."""
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import ListedColormap
from mevo_mb.control.scheduled_pzt import Schedule
from .common import BASE, OLD, Q, CONFIG, ROOT, read, write, sha, identity

PURPLE="#70559c";TEAL="#258e95";RED="#b45a60";GREY="#697580"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.labelcolor":"#262b32","text.color":"#262b32",
                     "xtick.color":"#262b32","ytick.color":"#262b32","axes.spines.top":False,"axes.spines.right":False,
                     "pdf.fonttype":42,"savefig.facecolor":"white","axes.titlesize":11})
OUT=BASE/"figures";AN=BASE/"analysis";outputs=[];inputs=set()


def js(path):
    inputs.add(path);return read(path)


def arrays(name,old=False):
    path=(Q/"analysis" if old else AN)/(name+".npz")
    inputs.add(path);return dict(np.load(path))


def finish(fig,name,title,foot=""):
    fig.suptitle(title,fontsize=15,y=.995)
    fig.text(.5,.012,foot+"\nEXPLORATORY · HOMOGENEOUS_WATER_CONTROL_PROXY · PHYSICS_UNCALIBRATED · training_ready=false",
             ha="center",va="bottom",fontsize=7.4,color="#444950")
    fig.tight_layout(rect=(0,.060,1,.965),h_pad=2.5,w_pad=2.8)
    for ext in ["png","pdf"]:
        p=OUT/(name+"."+ext);fig.savefig(p,dpi=220);outputs.append(p)
    plt.close(fig)


def geometry(ax):
    geo=js(OLD/"CONFIG_SNAPSHOT.json")["geometry"]
    theta=np.linspace(0,2*np.pi,64);xx=np.linspace(-16.25,16.25,50)
    ax.plot_surface(np.broadcast_to(xx,(64,50)),.8*np.cos(theta)[:,None]*np.ones((1,50)),
                    .8*np.sin(theta)[:,None]*np.ones((1,50)),color=RED,alpha=.30,linewidth=0,shade=True)
    cap_r,cap_t=np.meshgrid(np.linspace(0,.8,20),theta)
    for end in [-16.25,16.25]:
        ax.plot_surface(np.full_like(cap_r,end),cap_r*np.cos(cap_t),cap_r*np.sin(cap_t),color=RED,alpha=.45,linewidth=0,shade=True)
    for x in [-10.85,10.85]:
        z=np.linspace(0,4.5,24)
        ax.plot_surface(x+.75*np.cos(theta)[:,None]*np.ones((1,24)),.75*np.sin(theta)[:,None]*np.ones((1,24)),
                        np.broadcast_to(z,(64,24)),color=RED,alpha=.28,linewidth=0,shade=True)
    rr,tt=np.meshgrid(np.linspace(0,3.5,24),theta)
    for label,element in zip(["Y−","Y+","X−","X+"],geo["pzt"]):
        n=np.array(element["normal_local"]);center=np.array(element["center_local_m"])*1e3
        u=np.cross(n,[0,0,1]);v=np.cross(n,u)
        disk=center[:,None,None]+u[:,None,None]*rr*np.cos(tt)+v[:,None,None]*rr*np.sin(tt)
        ax.plot_surface(*disk,color="#889aa5",alpha=.65,shade=True,linewidth=0)
        ax.quiver(*center,*(n*2),color="#354654",linewidth=1,arrow_length_ratio=.3)
        ax.text(*(center-n*.8),label,fontsize=10,color="#27333a")
    for x in [-16.25,16.25]:
        ax.plot([x]*5,[-4.5,4.5,4.5,-4.5,-4.5],[-4.5,-4.5,4.5,4.5,-4.5],color=GREY,lw=.6,alpha=.7)
    for y in [-4.5,4.5]:
        for z in [-4.5,4.5]:ax.plot([-16.25,16.25],[y,y],[z,z],color=GREY,lw=.6,alpha=.7)
    ax.scatter([0],[0],[0],color=TEAL,s=30,label="Initial group COM")
    ax.set(xlabel="X (mm)",ylabel="Y (mm)",zlabel="Z (mm)",xlim=(-18,18),ylim=(-6,6),zlim=(-5,6))
    ax.set_box_aspect((32.5,12,11));ax.view_init(24,-64)
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.set_alpha(0)
    ax.set_title("a  Actual four-disc layout; 7 mm diameter")
    ax.text2D(.02,.04,"Water-filled, no flow; X ends sealed\nPZT discs shown as surfaces\nPlanned A: X−/X+; B: Y−/Y+",transform=ax.transAxes,fontsize=8)


def fig1():
    fig=plt.figure(figsize=(15.5,9.5));grid=fig.add_gridspec(2,3,width_ratios=[1.45,1.05,1])
    geometry(fig.add_subplot(grid[:,0],projection="3d"))
    plans=js(BASE/"RUN_SCHEDULES.json")
    times=np.linspace(0,.002,401,endpoint=False)
    amp=[];phase=[];labels=[]
    for short,key in [("One-hot","ONE_HOT_X_SEQUENCE"),("Pair","X_PAIR"),("Four","X_SEGMENTED")]:
        s=Schedule.from_record(plans[key]);a=[];p=[]
        for t in times:
            seg=next(v for v in s.segments if v.start_s<=t<v.start_s+v.duration_s)
            a.append(np.array(seg.amplitude_Pa_equivalent)*seg.enabled)
            p.append((np.rad2deg(seg.phase_rad)+180)%360-180)
        amp.extend(np.array(a).T);phase.extend(np.array(p).T)
        labels.extend([short+" · "+ch for ch in ["Y−","Y+","X−","X+"]])
    for pos,data,title,cmap,lim,unit in [(grid[0,1],amp,"b  Registered source amplitudes","Purples",(0,1500),"Pa-equivalent"),
                                      (grid[1,1],phase,"c  Registered phase offsets","twilight",(-180,180),"Phase (degree)")]:
        ax=fig.add_subplot(pos);im=ax.imshow(data,aspect="auto",extent=[0,2,12,0],cmap=cmap,vmin=lim[0],vmax=lim[1],interpolation="nearest")
        ax.set(yticks=np.arange(12)+.5,yticklabels=labels,xlabel="Ideal scheduled time (ms)",title=title)
        ax.axvline(1,color="#444",ls="--",lw=.6);ax.axhline(4,color="white");ax.axhline(8,color="white")
        fig.colorbar(im,ax=ax,pad=.02,label=unit,fraction=.046)
    ax=fig.add_subplot(grid[0,2]);d=js(BASE/"diagnostics/SINGLE_PZT.json")
    ax.bar([r["channel"] for r in d["rows"]],[r["pressure_magnitude_Pa"] for r in d["rows"]],color=[TEAL,TEAL,PURPLE,PURPLE])
    ax.set(yscale="log",ylabel="Group-center |P| (Pa)",title="d  One disc on: actual pressure")
    ax.text(.02,.96,"Each source: 1500 Pa-equivalent\nY/X pressure efficiency = 1/202",transform=ax.transAxes,va="top",fontsize=9)
    ax.set_ylim(5,10000)
    ax=fig.add_subplot(grid[1,2]);peaks=js(BASE/"diagnostics/AXIS_PAIR_PRESSURE_MAXIMA_REVIEWED.json")
    for axis,color in [("X",PURPLE),("Y",TEAL)]:
        d=peaks[axis];rows=d["rows"]
        ax.plot([r["delta_plus_minus_deg"] for r in rows],[(r["selected_peak_m"]-d["initial_peak_m"])*1e6 for r in rows],"o-",ms=3,color=color,label=f"{axis}: {d['branch_slope_um_per_degree']:.3f} µm/°")
    ax.set(xlabel="Phase(+) − phase(−) (degree)",ylabel="Axial pressure-peak shift (µm)",title="e  Field-only continuous peak branches")
    ax.legend(fontsize=8);ax.text(.02,.04,"Y: outer maximum, not the weak center peak\nPressure maxima are unverified trap candidates",transform=ax.transAxes,fontsize=7.6)
    finish(fig,"Fig_F10R_1_segmented_control","Simulation control sequences for four PZT discs",
           "IDEAL_SCHEDULED_SEGMENTS · Pa-equivalent is not Vpp · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · No hardware execution\nAll Y results: HOMOGENEOUS_WATER_PROXY_SPECIFIC")


def fig2():
    fig,axes=plt.subplots(2,3,figsize=(14.8,9.1));a,b,c,d,e,f=axes.ravel()
    runs=js(AN/"RUN_ANALYSIS.json");colors={"OFF_DT1":GREY,"X_COMPAT_DT1":"#303840","X_PAIR_DT1":PURPLE,"X_SEGMENTED_DT1":TEAL,"ONE_HOT_SEQUENCE_DT1":RED}
    labels={"OFF_DT1":"OFF","X_COMPAT_DT1":"F10Q replay","X_PAIR_DT1":"X pair","X_SEGMENTED_DT1":"Two segments","ONE_HOT_SEQUENCE_DT1":"One-hot sequence"}
    for name,color in colors.items():
        if name not in runs:continue
        m=arrays(name);a.plot(m["time_s"]*1e3,m["centroid_m"][:,0]*1e6,label=labels[name],color=color,lw=1.6)
    a.axhspan(.125,.275,color=TEAL,alpha=.09);a.axhline(.2,color=GREY,ls="--",lw=.7)
    a.set(xlabel="Time (ms)",ylabel="COM X (µm)",title="a  Fixed X target: 0.2 ± 0.075 µm");a.legend(fontsize=8)
    names=[k for k in colors if k in runs]
    b.barh([labels[k] for k in names],[runs[k]["target_error_m"]["X"]*1e6 for k in names],color=[colors[k] for k in names])
    b.axvline(.075,color=RED,ls="--",lw=1,label="Tolerance")
    b.set(xlabel="Endpoint target error (µm)",title="b  Observed endpoint error");b.legend(fontsize=8)
    phase=js(AN/"X_PHASE_RESULTS.json");rows=phase["field"]["rows"]
    c.plot([r["delta_Xplus_deg"] for r in rows],[r["Fx_ratio"] for r in rows],"o-",color=PURPLE,label="Center primary force")
    x=np.array([-10,10]);c.plot(x,1+.06*x,"--",color=GREY,label="User 6%/° hypothesis")
    base=runs["X_COMPAT_DT1"]["endpoint_m"][0]
    c.scatter([r["shift_deg"] for r in phase["actual_2ms"]],[r["endpoint_m"][0]/base for r in phase["actual_2ms"]],marker="x",s=65,color=TEAL,label="Actual 2 ms X displacement",zorder=5)
    c.set(xlabel="Only X+ phase change (degree)",ylabel="Ratio to original X control",title="c  Phase sensitivity: 3.73%/° locally");c.legend(fontsize=7.5)
    for name in ["X_COMPAT_DT1","X_PAIR_DT1","X_SEGMENTED_DT1"]:
        if name not in runs:continue
        m=arrays(name);d.plot(m["time_s"]*1e3,(m["rg_m"]-m["rg_m"][0])*1e9,color=colors[name],label=labels[name])
    d.set(xlabel="Time (ms)",ylabel="Rg change (nm)",title="d  Translation and morphology separated");d.legend(fontsize=8)
    rob=js(AN/"X_INITIAL_ROBUSTNESS.json")
    names=[k for k in ["X_COMPAT_DT1","X_rotate_z_plus5","X_rotate_z_minus5"] if k in rob]
    e.bar(["Original","Rotate +5°","Rotate −5°"][:len(names)],[rob[k]["endpoint_m"][0]*1e6 for k in names],color=[GREY,TEAL,PURPLE][:len(names)])
    e.axhspan(.125,.275,color=TEAL,alpha=.09);e.axhline(.2,color=GREY,ls="--",lw=.8)
    e.set(ylabel="Actual endpoint X (µm)",title="e  Three prescribed initial orientations",ylim=(0,.3))
    e.text(.02,.90,"Rotated states: single-resolution observations\nNot independent experimental samples",transform=e.transAxes,fontsize=7.5,va="top")
    comps=js(AN/"COMPARISONS.json")
    names=[n for n in ["X_COMPAT_DT1","X_PAIR_DT1","X_SEGMENTED_DT1","X_PHASE_+10_DT1"] if n+"__DT05" in comps and n+"__FINE" in comps]
    vals=[max(comps[n+"__"+s]["fixed_id_m"] for s in ["DT05","FINE"])*1e9 for n in names]
    f.bar(range(len(names)),vals,color=PURPLE);f.axhline(50,color=RED,ls="--",lw=.8,label="Fixed-ID limit: 50 nm")
    f.set(xticks=range(len(names)),xticklabels=["Replay","Pair","Segments","Phase +10°"][:len(names)],ylabel="Max fixed-ID difference (nm)",yscale="log",title="f  Actual dt/cache differences")
    f.legend(fontsize=8)
    finish(fig,"Fig_F10R_2_x_navigation_robustness","X navigation, control modes and bounded robustness",
           "Target and initial states fixed before trajectories; single-resolution observations are separate from numerical admission")


def fig3():
    fig,axes=plt.subplots(2,3,figsize=(14.8,9.1));a,b,c,d,e,f=axes.ravel()
    gate=js(AN/"FINAL_GATE.json");key=gate["new_Y_selected"]
    entries=[("F10Q Y / dt", "NAV_Y_DT1",False,PURPLE),("F10Q Y / half dt","NAV_Y_DT05",False,RED),("F10Q Y / fine cache","NAV_Y_FINE",False,GREY)]
    if key:entries.append(("New low-pressure Y",key+"_DT1",True,TEAL))
    for label,name,new,color in entries:
        m=arrays(name,old=not new)
        a.plot(m["time_s"]*1e3,m["centroid_m"][:,1]*1e6,label=label,color=color)
        b.plot(m["time_s"]*1e3,(m["rg_m"]-m["rg_m"][0])*1e9,label=label,color=color)
    a.axhline(.2,color=GREY,ls="--",lw=.8);a.axhspan(.125,.275,color=TEAL,alpha=.08)
    a.set(xlabel="Time (ms)",ylabel="COM Y (µm)",title="a  Y target remains 0.2 µm");a.legend(fontsize=7.3)
    b.set(xlabel="Time (ms)",ylabel="Rg change (nm)",title="b  Historical morphology disagreement")
    if key:
        r=js(AN/"RUN_ANALYSIS.json")[key+"_DT1"]
        b.text(.03,.08,f"New Y: shape / fixed-ID checks PASS\nActual ΔRg = {(r['rg_initial_final_m'][1]-r['rg_initial_final_m'][0])*1e9:.3f} nm\nTarget NOT REACHED",transform=b.transAxes,fontsize=8)
    for suffix,style in [("DT05","-"),("FINE","--")]:
        m=arrays("compare_NAV_Y_"+suffix,old=True)
        c.plot(m["time_s"]*1e3,m["fixed_id_m"]*1e9,color=PURPLE,ls=style,label="F10Q "+suffix)
        d.plot(m["time_s"]*1e3,m["Q_relative"],color=PURPLE,ls=style)
        if key and (AN/("compare_"+key+"_DT1__"+suffix+".npz")).exists():
            m=arrays("compare_"+key+"_DT1__"+suffix)
            c.plot(m["time_s"]*1e3,m["fixed_id_m"]*1e9,color=TEAL,ls=style,label="New Y "+suffix)
            d.plot(m["time_s"]*1e3,m["Q_relative"],color=TEAL,ls=style)
    c.axhline(50,color=RED,ls=":",lw=.9);c.set(xlabel="Time (ms)",ylabel="Fixed-ID difference (nm)",title="c  No ID relabeling or shape simplification");c.legend(fontsize=7.3)
    d.axhline(.01,color=RED,ls=":",lw=.9);d.set(xlabel="Time (ms)",ylabel="Relative shape-tensor difference",title="d  Same original Q tolerance: 0.01")
    if key:
        m=arrays(key+"_DT1");coords=m["centered_positions_m"][-1]*1e6
        for i,(xx,yy,zz) in enumerate(coords):
            e.add_patch(Circle((xx,yy),1.78,facecolor=TEAL,edgecolor="#164952",alpha=.32,lw=.6))
        # Both actual Z layers project onto the same XY circles; label all fixed IDs.
        initial=m["centered_positions_m"][0]
        used=set()
        for i in range(len(coords)):
            if i in used:continue
            group=np.flatnonzero(np.linalg.norm(initial[:,:2]-initial[i,:2],axis=1)<1e-12)
            used.update(group.tolist());location=coords[group,:2].mean(0)
            e.text(*location,"/".join(map(str,group)),ha="center",va="center",fontsize=6.3,color="#173b43")
        e.set(xlim=(-9,9),ylim=(-9,9),aspect="equal",xlabel="Centered X (µm)",ylabel="Centered Y (µm)",title="e  New Y: actual 32-bubble projection")
        e.text(.03,.02,"True radii; overlapping XY projections of two Z layers",transform=e.transAxes,fontsize=7.0)
    else:e.text(.5,.5,"NEW Y: NOT_RUN",ha="center",transform=e.transAxes)
    runs=js(AN/"RUN_ANALYSIS.json");right=f.twinx()
    for name,path,color,label in [("NAV_Y_DT1",Q/"b_navigation_v1/runs/NAV_Y_DT1/trajectory.h5",PURPLE,"F10Q Y")]+(
        [(key+"_DT1",BASE/"navigation_v1/runs"/(key+"_DT1")/"trajectory.h5",TEAL,"New Y")] if key else []):
        inputs.add(path)
        with h5py.File(path) as h:time_s=h["time_s"][:];ex=h["extra_metrics"][:]
        f.plot(time_s[:len(ex)]*1e3,ex[:,0],color=color,label=label+" contacts")
        right.plot(time_s[:len(ex)]*1e3,ex[:,4],color=color,ls="--",alpha=.6)
    f.set(xlabel="Time (ms)",ylabel="Active contact pairs",title="f  Contact activity and actual pressure")
    right.set_ylabel("Maximum bubble pressure (Pa)");f.legend(fontsize=7.3)
    finish(fig,"Fig_F10R_3_y_morphology_sensitivity","Y motion and morphology: historical versus lower-pressure control",
           "All Y conclusions: HOMOGENEOUS_WATER_PROXY_SPECIFIC; no extrapolation to the real PDMS apparatus. Dashed panel f: pressure.")


def fig4():
    fig=plt.figure(figsize=(14.5,9));grid=fig.add_gridspec(2,2,width_ratios=[1.5,1])
    a=fig.add_subplot(grid[:,0]);b=fig.add_subplot(grid[0,1]);c=fig.add_subplot(grid[1,1])
    gate=js(AN/"FINAL_GATE.json");runs=js(AN/"RUN_ANALYSIS.json");cost=js(AN/"COSTS.json")
    names=[n for n in ["OFF_DT1","X_COMPAT_DT1","ONE_HOT_SEQUENCE_DT1","X_PAIR_DT1","X_SEGMENTED_DT1","X_PHASE_-6_DT1","X_PHASE_+6_DT1","X_PHASE_+10_DT1","X_rotate_z_plus5","X_rotate_z_minus5",(gate["new_Y_selected"] or "NOT_RUN")+"_DT1"] if n in gate["groups"]]
    cols=["BASIC_NUMERICS_PASS","GROUP_POSITION_NUMERICALLY_ADMITTED","GROUP_SHAPE_NUMERICALLY_ADMITTED","GROUP_TOPOLOGY_NUMERICALLY_ADMITTED","PARTICLE_ID_NUMERICALLY_ADMITTED"]
    table=[];texts=[]
    for name in names:
        g=gate["groups"][name];tt=[g[k] for k in cols]
        axis="Y" if name.startswith("Y_") else "X"
        target="PASS" if g["TARGET_REACHED"][axis] else ("FAIL" if g["GROUP_POSITION_NUMERICALLY_ADMITTED"]!="NOT_RUN" else "NOT_RUN")
        tt.append(target);texts.append(tt);table.append([{"PASS":1,"FAIL":2,"NOT_RUN":0}[s] for s in tt])
    h=gate["historical_Y"];texts.append([h[k] for k in cols]+["FAIL"]);table.append([{"PASS":1,"FAIL":2,"NOT_RUN":0}[s] for s in texts[-1]])
    labels=[n.replace("_DT1","").replace("X_COMPAT","X replay").replace("ONE_HOT_SEQUENCE","One-hot X sequence").replace("X_rotate_z_","X rotation ") for n in names]+["Historical F10Q Y"]
    a.imshow(table,cmap=ListedColormap(["#e5e7ea","#c2dedc","#e7c3c6"]),vmin=0,vmax=2,aspect="auto")
    for i,row in enumerate(texts):
        for j,value in enumerate(row):a.text(j,i,"NOT RUN" if value=="NOT_RUN" else value,ha="center",va="center",fontsize=8)
    a.set(yticks=range(len(labels)),yticklabels=labels,xticks=range(6),xticklabels=["Basic","COM","Shape","Graph","Fixed ID","Target"],title="a  Independent admission layers")
    a.tick_params(axis="x",rotation=30)
    camps=[]
    for p in sorted(BASE.glob("*/CAMPAIGN.json")):
        d=js(p);camps.append((p.parent.name,d["wall_s"]))
    b.barh([x[0].replace("_v1","") for x in camps],[x[1] for x in camps],color=PURPLE)
    b.set_ylim(-1.45,len(camps)-.4)
    b.set(xlabel="Actual process wall time (s)",title="b  Finite scientific execution cost")
    b.text(.03,.03,f"Total numerical + analysis: {cost['total_scientific_s']:.2f} s\n{cost['trajectory_attempts']} attempts; {cost['failed_attempts']} trajectory failures\nHybrid GPU field/force/walls + CPU QP",transform=b.transAxes,fontsize=8)
    c.axis("off");c.set_title("c  Scope and unmeasured quantities",loc="left")
    lines=["Three simulation control modes: implemented","Exact F10Q X replay: PASS (0 coordinate difference)",
           "Numerical thresholds: unchanged","Ideal scheduled events; measured hardware delay: NOT RUN", "Instrument commands / real voltage output: NOT RUN",
           "Pa-equivalent → Vpp calibration: NOT AVAILABLE","PDMS reflection and streaming: NOT SIMULATED",
           "Long-time moving-trap navigation: NOT RUN","Experimental validation: NOT VERIFIED","AI training admission: false"]
    c.text(.01,.97,"\n\n".join(lines),transform=c.transAxes,va="top",fontsize=9)
    finish(fig,"Fig_F10R_4_control_admission","Control admission, numerical scope and computational cost",
           "NOT RUN is distinct from FAIL; short/rotated-state observations do not inherit unperformed numerical comparisons\nAll Y results: HOMOGENEOUS_WATER_PROXY_SPECIFIC")


def main():
    OUT.mkdir(exist_ok=True)
    for f in [fig1,fig2,fig3,fig4]:f()
    write(OUT/"FIGURE_PROVENANCE.json",dict(identity=identity(),script_sha256=sha(__import__('pathlib').Path(__file__)),
          input_sha256={str(p.relative_to(ROOT)):sha(p) for p in sorted(inputs)},
          figures_sha256={p.name:sha(p) for p in outputs},data_source="Actual archived calculations only",
          no_synthetic_trajectories=True,style="White background, muted purple/teal/red, English text, no smoothing"))
    print("Saved",len(outputs),"PNG/PDF outputs")


if __name__=="__main__":main()
