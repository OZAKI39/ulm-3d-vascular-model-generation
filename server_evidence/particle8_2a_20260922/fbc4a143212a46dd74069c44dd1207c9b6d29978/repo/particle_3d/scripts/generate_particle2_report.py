#!/usr/bin/env python3
"""Eleven Particle-2 scientific figures with durable CSV/JSON plotting data."""
from pathlib import Path
import argparse
import csv
import importlib.util
import json
from math import erf,sqrt,log
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch

PACKAGE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PACKAGE/"src"))
from particle_3d.rbc_distribution import load_contract,digest
from particle_3d.rbc_orientation import shape_axis_distance
from particle_3d.audit import read_frozen
from particle_3d.particle2_audit import check_dependencies,P2_BRANCH

REPORT=PACKAGE/"reports/particle2";DATA=REPORT/"data";FIGURES=REPORT/"figures"
COLORS=["#287ca6","#dc8532","#489a65","#a365a3","#d0524c"]
NAMES=[
    "00_particle2_scope_and_literature.png","01_rbc_diameter_volume_distribution.png",
    "02_rbc_derived_geometry_distribution.png","03_rbc_joint_geometry_scatter.png",
    "04_static_flow_orientation_validation.png","05_rigid_rotation_orientation_validation.png",
    "06_simple_shear_jeffery_orbits.png","07_distribution_wide_jeffery_validation.png",
    "08_real_fem_rbc_orientation_trajectory.png","09_real_fem_orientation_diagnostics.png",
    "10_particle2_validation_timestep_comparison.png"]
for candidate in ["/mnt/c/Windows/Fonts/msyh.ttc","/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"]:
    if Path(candidate).is_file():
        font_manager.fontManager.addfont(candidate)
        plt.rcParams["font.sans-serif"]=[font_manager.FontProperties(fname=candidate).get_name(),"DejaVu Sans"]
        bold=Path("/mnt/c/Windows/Fonts/msyhbd.ttc")
        if bold.is_file():font_manager.fontManager.addfont(str(bold))
        break
plt.rcParams.update({"font.size":11,"axes.unicode_minus":False,"savefig.dpi":180})


def read_json(name):return json.loads((DATA/name).read_text())


def write_json(name,value):
    (DATA/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False,
                                    default=lambda x:x.item() if isinstance(x,np.generic) else x)+"\n")


def read_csv(name):
    with (DATA/name).open(newline="") as stream:return list(csv.DictReader(stream))


def write_csv(name,rows):
    with (DATA/name).open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator="\n")
        writer.writeheader();writer.writerows(rows)


def col(rows,name):return np.array([float(r[name]) for r in rows])


def vectors(rows,prefix,axes="xyz",suffix=""):
    return np.column_stack([col(rows,prefix+a+suffix) for a in axes])


def save(fig,index,footer,note,sources):
    fig.text(.5,.025,footer,ha="center",va="bottom",fontsize=10,
             bbox=dict(facecolor="#eff3f5",edgecolor="none",pad=7))
    fig.savefig(FIGURES/NAMES[index],bbox_inches="tight");plt.close(fig)
    write_json(f"{index:02d}_figure_sources.json",dict(figure=NAMES[index],sources=sources,
               source_sha256={s:digest(DATA/s) for s in sources},timestep_role="VALIDATION_ONLY",
               manual_visual_review="PENDING_USER_REVIEW"))
    if index:
        (REPORT/f"{index:02d}_step_notes.md").write_text(f"# Particle-2 步骤 {index:02d}\n\n{note}\n\n![审核图](figures/{NAMES[index]})\n")


def stage0():
    p0,p1=check_dependencies(PACKAGE.parent)
    c=load_contract()
    sources=json.loads((REPORT/"literature/sources.json").read_text())
    write_json("00_scope_literature.json",dict(contract=c,literature=sources,branch=P2_BRANCH,
        particle0_commit=p0["git_commit"],particle1_commit=p1["git_commit"],
        mesh_sha256=p1["mesh_sha256"],flow_sha256=p1["flow_sha256"],
        particle_count_per_dynamic_replay=1,finite_size_wall_clearance="NOT_VALIDATED_PARTICLE3"))
    fig,ax=plt.subplots(figsize=(14,9));ax.axis("off")
    boxes=[
        (.02,.69,.46,.22,"MEASURED / REPORTED","#deebf7",
         "Moss 2025 PREPRINT, Table 3: mature C57BL/6Case RBC\nD mean 6.79 µm; SD 0.93; median 6.67; IQR 1.33\nTable 4: MCV 47.9 fL; ±2.6 is between mice"),
        (.53,.69,.45,.22,"REPORTED CROSS-CHECKS","#deebf7",
         "Rivera 2013: C57BL/6J MCV 47.8 / 48.4 fL\nDe Franceschi 2005: MCV 49.3 fL; RDW 13.2%\nDifferent substrains / cohorts are not identical"),
        (.06,.36,.88,.22,"MODEL ASSUMPTIONS — USER-AUTHORIZED V0","#ffead5",
         "D ~ Normal(6.79, 0.93²) µm; V ~ Normal(47.9, 7.0892²) fL\nCV(V)=14.8% is MODEL CHOICE; D/V independent due to missing joint data\n3σ guards: D [4.00,9.58] µm; V [26.6324,69.1676] fL; reject c≥a"),
        (.13,.08,.74,.17,"DERIVED GEOMETRY","#dff1e5",
         "a=b=D/2; c=3V/(πD²); r=c/a; λ=(r²−1)/(r²+1)<0\nRigid oblate spheroid; all particle dynamics in SI; short axis p=R(q)e3")]
    for x,y,w,h,title,color,detail in boxes:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=.015",facecolor=color,edgecolor="#647889"))
        ax.text(x+w/2,y+h-.035,title,ha="center",va="top",fontsize=12,weight="bold")
        ax.text(x+w/2,y+h-.09,detail,ha="center",va="top",fontsize=10,linespacing=1.8)
    for x in [.27,.74]:ax.annotate("",xy=(.5,.60),xytext=(x,.675),arrowprops=dict(arrowstyle="->",lw=2))
    ax.annotate("",xy=(.5,.27),xytext=(.5,.34),arrowprops=dict(arrowstyle="->",lw=2))
    fig.suptitle("Particle-2: literature → V0 distribution → single rigid RBC orientation",fontsize=17)
    fig.subplots_adjust(bottom=.17,top=.92)
    save(fig,0,"P0 / P1 accepted dependencies; frozen FEM READ ONLY. 100,000 geometries are VALIDATION ONLY.\nNormality, independence, CV and guards are model choices. No production timestep or inlet orientation distribution.","",["00_scope_literature.json"])


def population_rows():return read_csv("C57BL6_RBC_GEOMETRY_VALIDATION_100000.csv")


def stage1():
    pop=population_rows();c=load_contract();fig,axes=plt.subplots(1,2,figsize=(14,7));plot=[]
    for ax,name,unit,mean,sd,guard,title in [
        (axes[0],"D_um","µm",6.79,.93,[4,9.58],"Diameter: reported mean/SD; Normal approximation"),
        (axes[1],"V_fL","fL",47.9,7.0892,[26.6324,69.1676],"Volume: MCV anchor + MODEL CV 14.8%")]:
        x=col(pop,name);density,edges=np.histogram(x,bins=55,range=guard,density=True)
        centers=(edges[:-1]+edges[1:])/2
        latent=np.exp(-.5*((centers-mean)/sd)**2)/(sd*np.sqrt(2*np.pi))
        truncated=latent/erf(3/np.sqrt(2))
        ax.bar(centers,density,width=np.diff(edges),color="#91b9cc",alpha=.8,label="Final accepted geometry population")
        ax.plot(centers,latent,"--",color="#888888",label="Latent Normal model")
        ax.plot(centers,truncated,color="#b15627",label="Guard-only truncated target (before shape rejection)")
        ax.set(xlabel=f"{name.split('_')[0]} ({unit})",ylabel=f"Probability density ({unit}^-1)",title=title)
        ax.set_ylim(0,1.6*max(density.max(),truncated.max()))
        ax.grid(alpha=.2);ax.legend(fontsize=9)
        plot.extend(dict(variable=name,bin_left=float(a),bin_right=float(b),accepted_density=float(d),
                         latent_density=float(l),guard_only_target_density=float(t)) for a,b,d,l,t in zip(edges[:-1],edges[1:],density,latent,truncated))
    axes[0].text(.03,.72,"Reported: mean 6.79; SD 0.93 µm\nmedian 6.67; IQR 1.33 µm\nMoss Table 3: PREPRINT",transform=axes[0].transAxes,fontsize=10)
    axes[1].text(.03,.72,"Model: mean 47.9 fL; CV 14.8%\nSD 7.0892 fL is derived\n±2.6 fL is NOT cell-volume SD",transform=axes[1].transAxes,fontsize=10)
    write_csv("01_distribution_histograms.csv",plot)
    # CDF smoke on latent and guard-only prefixes, deliberately before c<a.
    ledger=read_csv("candidate_ledger.csv");checks=[]
    guarded=[r for r in ledger if r["status"] in ["ACCEPT","SHAPE"]]
    for name,mean,sd in [("D_raw_um",6.79,.93),("V_raw_fL",47.9,7.0892)]:
        for label,rows in [("latent",ledger),("guard_only_before_shape",guarded)]:
            x=np.sort(col(rows,name));n=len(x)
            cdf=np.array([.5*(1+erf((v-mean)/(sd*sqrt(2)))) for v in x])
            if label!="latent":cdf=(cdf-.5*(1+erf(-3/sqrt(2))))/erf(3/sqrt(2))
            error=max(float(np.max(np.arange(1,n+1)/n-cdf)),float(np.max(cdf-np.arange(n)/n)))
            bound=sqrt(log(2/1e-6)/(2*n))
            checks.append(dict(variable=name,population=label,N=n,cdf_max_distance=error,dkw_alpha=1e-6,predefined_dkw_bound=bound,passed=error<=bound))
    write_json("01_distribution_cdf_checks.json",checks)
    fig.subplots_adjust(bottom=.22,wspace=.26,top=.87)
    save(fig,1,"N=100,000; accepted distribution includes 3σ guards AND c<a rejection.\nGuards are MODEL limits, not biological min/max. Accepted SD need not equal latent SD.",
         "左图与右图分别看直径和体积，蓝柱是最终接受样本，曲线是形状筛选前的模型参照。原始正态与仅通过范围保护的样本分别通过 CDF 检查；最终样本不强求 SD 等于潜在模型 SD。",
         ["01_distribution_histograms.csv","01_distribution_cdf_checks.json","C57BL6_RBC_GEOMETRY_VALIDATION_100000.csv"])


def stage2():
    pop=population_rows();fig,axes=plt.subplots(2,3,figsize=(15,9));plot=[]
    for ax,key,label,unit,factor in zip(axes.ravel(),["a_um","c_um","c_um","r","jeffery_lambda"],["a=b","c","Full thickness 2c","r=c/a","Jeffery λ"],["µm","µm","µm","dimensionless","dimensionless"],[1,1,2,1,1]):
        values=col(pop,key)*factor;density,edges=np.histogram(values,bins=50,density=True)
        ax.bar((edges[:-1]+edges[1:])/2,density,width=np.diff(edges),color="#67a087",alpha=.85)
        ax.set(xlabel=f"{label} ({unit})",ylabel="Density (1/µm)" if unit=="µm" else "Probability density (dimensionless)",title=f"Derived {label}");ax.grid(alpha=.2)
        plot.extend(dict(variable=label,unit=unit,bin_left=float(a),bin_right=float(b),density=float(d)) for a,b,d in zip(edges[:-1],edges[1:],density))
    axes[1,2].axis("off");axes[1,2].text(.05,.95,"DERIVED FROM V0 MODEL\n\na=b=D/2\nc=3V/(πD²)\nfull thickness=2c\nr=c/a,   0<r<1\n−1<λ<0\n\nThese are not direct\nliterature measurements.",va="top",fontsize=11,linespacing=1.25)
    write_csv("02_derived_geometry_histograms.csv",plot)
    fig.subplots_adjust(bottom=.16,hspace=.4,wspace=.3)
    save(fig,2,"All 100,000 accepted geometries: a=b>c>0; volume reconstructs V; all core axes converted to m.\nc / thickness / r / λ are MODEL-DERIVED, not measured single-cell distributions.",
         "五个直方图分别展示宽度半轴、短轴、完整厚度、厚薄比和 Jeffery 系数。薄球的 r 小、λ 更接近 −1；全部满足 a=b>c>0。这些是模型推导统计。",
         ["02_derived_geometry_histograms.csv","accepted_statistics.json"])


def stage3():
    pop=population_rows();ledger=read_csv("candidate_ledger.csv");rejected=[r for r in ledger if r["status"]=="SHAPE"]
    points=[]
    for rows,status,n in [(pop,"ACCEPT",2200),(rejected,"SHAPE_REJECT",150)]:
        for i in np.linspace(0,len(rows)-1,min(n,len(rows)),dtype=int):
            row=rows[i];d=float(row["D_um"] if status=="ACCEPT" else row["D_raw_um"]);v=float(row["V_fL"] if status=="ACCEPT" else row["V_raw_fL"])
            c=3*v/(np.pi*d*d)
            points.append(dict(candidate_id=int(row["candidate_id"]),D_um=d,V_fL=v,full_thickness_um=2*c,r=2*c/d,status=status))
    write_csv("03_joint_scatter_subsample.csv",points)
    fig,axes=plt.subplots(1,3,figsize=(16,6.5))
    for ax,key,label in zip(axes,["V_fL","full_thickness_um","r"],["V (fL)","Full thickness 2c (µm)","r=c/a (dimensionless)"]):
        for status,color,marker in [("ACCEPT","#2b8b99","."),("SHAPE_REJECT","#cf514c","x")]:
            subset=[r for r in points if r["status"]==status]
            ax.scatter(col(subset,"D_um"),col(subset,key),s=8 if status=="ACCEPT" else 18,alpha=.3 if status=="ACCEPT" else .8,color=color,marker=marker,label=status)
        ax.set(xlabel="D (µm)",ylabel=label);ax.grid(alpha=.2);ax.legend(fontsize=9)
    axes[1].plot([4,9.58],[4,9.58],"k--",lw=1,label="2c=D shape boundary")
    axes[2].axhline(1,color="k",ls="--",lw=1,label="c=a shape boundary")
    axes[1].legend(fontsize=9);axes[2].legend(fontsize=9)
    fig.suptitle("Independent latent D/V → guards → reject c≥a (retain original candidate values)",fontsize=16)
    fig.subplots_adjust(bottom=.22,top=.85,wspace=.27)
    save(fig,3,f"Plot: 2,200 accepted points + {min(150,len(rejected))} shape-rejected candidates. Statistics use all 100,000 accepted cells.\nShape conditioning can introduce D/V dependence. Red candidates are rejected whole; no clipping or shrinking.",
         "红叉是通过 D/V 范围保护但 c≥a 的原始候选，整对拒绝。图只抽取少量点便于阅读，统计用全部样本。接受后的相关结构由筛选产生，不是独立性已被实验确认。",
         ["03_joint_scatter_subsample.csv","candidate_ledger.csv","accepted_statistics.json"])


def stage4():
    rows=read_csv("04_static.csv");ids=sorted({r["rbc_id"] for r in rows},key=int);fig,axes=plt.subplots(1,3,figsize=(15,6.5))
    one=[r for r in rows if r["rbc_id"]==ids[0]]
    for a,color in zip("xyz",COLORS):axes[0].plot(col(one,"time_s"),col(one,"p_"+a),color=color,label="p_"+a)
    axes[0].set(xlabel="Time (s)",ylabel="Short-axis component (dimensionless)",title="Same fixed short axis for all five geometries")
    for i,rid in enumerate(ids):
        group=[r for r in rows if r["rbc_id"]==rid];label=f"r={float(group[0]['r']):.3f}"
        axes[1].plot(col(group,"time_s"),col(group,"rotation_matrix_error"),label=label,color=COLORS[i])
        axes[2].plot(col(group,"time_s"),col(group,"quaternion_norm_error"),label=label,color=COLORS[i])
    axes[1].set(xlabel="Time (s)",ylabel="Max rotation matrix component change",title="Check: unchanged orientation")
    axes[2].set(xlabel="Time (s)",ylabel="|norm(q)−1| (dimensionless)",title="Check: unit quaternion")
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=9)
    fig.subplots_adjust(bottom=.22,wspace=.3)
    save(fig,4,"Uniform translation; G=0 → E=W=0 → Ω=0. Five actual distribution geometries; no mean-RBC shortcut.\nAll dt values are VALIDATION_ONLY; biological inlet orientation distribution is not specified.",
         "流速恒定而梯度为零，RBC 中心移动，但方向保持不变。五种厚薄比例的旋转矩阵变化均为零，四元数保持单位长度。",
         ["04_static.csv","04_static_metrics.json","selected_geometries.json"])


def stage5():
    rows=read_csv("05_rotation.csv");ids=sorted({r["rbc_id"] for r in rows},key=int);one=[r for r in rows if r["rbc_id"]==ids[0]]
    fig,axes=plt.subplots(1,3,figsize=(15,6.5))
    for a,color in zip("xyz",COLORS):
        axes[0].plot(col(one,"time_s"),col(one,"p_"+a),color=color,label="computed p_"+a)
        axes[0].plot(col(one,"time_s"),col(one,"analytic_p_"+a),color=color,ls="--",lw=1,label="analytic "+a)
    axes[0].set(xlabel="Time (s)",ylabel="p component (dimensionless)",title="Known z-axis rigid rotation")
    for i,rid in enumerate(ids):
        group=[r for r in rows if r["rbc_id"]==rid];label=f"r={float(group[0]['r']):.3f}"
        axes[1].plot(col(group,"time_s"),col(group,"axis_vector_error"),color=COLORS[i],label=label)
        axes[2].plot(col(group,"time_s"),col(group,"quaternion_norm_error"),color=COLORS[i],label=label)
    axes[1].set(xlabel="Time (s)",ylabel="Max |p−p_exact| component",title="Check: analytic orientation")
    axes[2].set(xlabel="Time (s)",ylabel="|norm(q)−1|",title="Check: unit quaternion")
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.subplots_adjust(bottom=.22,wspace=.3)
    save(fig,5,"E=0 → Ω_RBC=Ω_fluid for every geometry. q=(w,x,y,z), body→world; world increment left-multiplies.\nDirection, multiplication order and q/−q equivalence have independent permanent tests. VALIDATION_ONLY.",
         "整体刚体旋转时没有拉伸项，五种 RBC 都跟随相同流体旋转。三个 p 分量与解析曲线重合，误差在预先按步数确定的舍入界内。",
         ["05_rotation.csv","05_rotation_metrics.json"])


def stage6():
    rows=read_csv("06_shear_dt2.csv");metrics=[read_json(f"06_shear_dt{i}_metrics.json") for i in range(3)]
    fig,axes=plt.subplots(1,2,figsize=(14,7))
    for j,m in enumerate(metrics[-1]):
        group=[r for r in rows if int(r["rbc_id"])==m["rbc_id"]]
        axes[0].plot(col(group,"time_s"),col(group,"unwrapped_angle_rad")/np.pi,color=COLORS[j],label=f"r={m['r']:.3f}")
    markers=["o","s","^"]
    for i,group in enumerate(metrics):
        axes[1].scatter([m["analytic_period_s"] for m in group],[m["measured_period_s"] for m in group],s=45,marker=markers[i],label=f"dt={group[0]['validation_dt_s']:.3g} s")
    limit=max(m["analytic_period_s"] for m in metrics[0])*1.1
    axes[1].plot([0,limit],[0,limit],"k--",label="Analytic exact period")
    axes[0].set(xlabel="Time (s)",ylabel="Unwrapped directed-axis angle / π",title="Different sampled r → different tumbling rate")
    axes[1].set(xlabel="Analytic shape-axis period (s)",ylabel="Measured shape-axis period (s)",title="Check: π/γ × (r+1/r)")
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=9)
    fig.subplots_adjust(bottom=.22,wspace=.25)
    save(fig,6,"Each half-turn p→−p repeats the ellipsoid shape. Directed-vector full period = 2 × shape-axis period.\nγ=20 s^-1 is SYNTHETIC_VALIDATION_PARAMETER; three dt values are VALIDATION_ONLY.",
         "五个代表样本来自 r 的 5%、25%、50%、75%、95% 附近。越薄的球翻转周期越长；数值周期与解析周期相符。这里 p 到 −p 已是形状重复，有向向量的整圈周期还要乘二。",
         ["06_shear_dt2.csv"]+[f"06_shear_dt{i}_metrics.json" for i in range(3)])


def stage7():
    metrics=read_json("07_distribution_wide_metrics.json");r=np.array([m["r"] for m in metrics])
    curve=np.linspace(r.min(),r.max(),400);period=np.pi/20*(curve+1/curve)
    write_csv("07_analytic_period_curve.csv",[dict(r=float(a),analytic_period_s=float(b)) for a,b in zip(curve,period)])
    fig,axes=plt.subplots(1,2,figsize=(14,7))
    axes[0].plot(curve,period,color="#555555",label="Analytic shape-axis period")
    axes[0].scatter(r,[m["measured_period_s"] for m in metrics],s=22,color="#247d91",label="64 distribution samples")
    axes[1].scatter(r,[m["relative_error"] for m in metrics],s=22,color="#c07737",label="Relative period error")
    axes[0].set(xlabel="r=c/a (dimensionless)",ylabel="Shape-axis period (s)",title="Same solver across 64 stratified random geometries")
    axes[1].set(xlabel="r=c/a (dimensionless)",ylabel="|T_measured−T_exact|/T_exact",title="Check: period accuracy across shape distribution")
    for ax in axes:ax.grid(alpha=.2);ax.legend()
    fig.subplots_adjust(bottom=.22,wspace=.3)
    save(fig,7,f"64 equal-count r strata; fixed selection seed=2026092064; geometry values taken from original population.\nMax period relative error={max(m['relative_error'] for m in metrics):.3e}; no interaction between these independent validation replays.",
         "按 r 排序分成 64 个等数量区间，每区用固定种子选一个原始样本。所有样本使用同一四元数求解器，薄和厚两端都通过周期及单位长度检查。",
         ["07_distribution_wide_metrics.json","07_distribution_wide_metrics.csv","07_analytic_period_curve.csv","orientation_sweep_selection.json"])


def real(g,d):return read_csv(f"08_real_g{g}_dt{d}.csv")


def stage8():
    rows=real(2,2)
    _,_,_,boundaries=read_frozen(PACKAGE.parent/"formal_3D_flow_solver/FEM_SimVascular")
    spec=importlib.util.spec_from_file_location("p0_plot_readonly",PACKAGE/"scripts/generate_particle0_report.py")
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    fig=plt.figure(figsize=(14,12));ax=fig.add_subplot(111,projection="3d",computed_zorder=False)
    module.vessel(ax,boundaries,.25)
    xyz=vectors(rows,"",suffix="_m")*1e6;p=vectors(rows,"p_")
    ax.plot(*xyz.T,color="#e98a2a",lw=2,label="RBC center path (all five geometries coincide)")
    indices=np.linspace(50,len(rows)-150,10,dtype=int)
    ax.quiver(*xyz[indices].T,*(p[indices]*3).T,color="#184b8a",arrow_length_ratio=.3,lw=1.4,zorder=5,label="Short axis +p (orientation, not velocity)")
    ax.quiver(*xyz[indices].T,*(-p[indices]*3).T,color="#184b8a",arrow_length_ratio=.3,lw=1.4,zorder=5)
    ax.scatter(*xyz[0],s=60,color="#d43844",zorder=6,label="Start: canonical tetra 208001")
    ax.scatter(*xyz[-1],s=90,marker="*",color="black",zorder=6,label="End: OUTLET_02")
    ax.legend(loc="upper left",bbox_to_anchor=(-.05,1.01),fontsize=10)
    fig.suptitle("Real frozen FEM: center path and short-axis orientation (median-r sample)",fontsize=15,y=.98)
    fig.subplots_adjust(bottom=.17,top=.86)
    save(fig,8,"Short-axis arrows show ±p; their display length is not an RBC radius. RBC GLYPH SIZE MAY BE VISUALLY EXAGGERATED.\nFINITE_SIZE_WALL_CLEARANCE = NOT VALIDATED. All 15 center replays exit OUTLET_02; no wall model or deformation.",
         "图展示中位 r 样本在最细验证 dt 下的短轴方向。五种几何的中心轨迹一致，但旋转不同；箭头是短轴，不是速度。中心未穿 WALL 不代表整个椭球未碰墙，有限尺寸间隙未验证。",
         ["08_real_g2_dt2.csv","08_real_g2_dt2.json","08_real_sweep.json","08_frozen_provenance.json"])


def stage9():
    fig,axes=plt.subplots(2,2,figsize=(15,10));median=real(2,2)
    for a,color in zip("xyz",COLORS):axes[0,0].plot(col(median,"time_s"),col(median,"p_"+a),color=color,label="p_"+a)
    sources=[]
    for g in range(5):
        rows=real(g,2);label=f"r={float(rows[0]['r']):.3f}";sources.append(f"08_real_g{g}_dt2.csv")
        axes[0,1].plot(col(rows,"time_s"),col(rows,"quaternion_norm_error"),color=COLORS[g],label=label,lw=.8)
        axes[1,0].plot(col(rows,"time_s"),col(rows,"omega_norm_s_inv"),color=COLORS[g],label=label,lw=.8)
    axes[1,1].plot(col(median,"time_s"),col(median,"tetra_id"),".",ms=1.2,label="Canonical tetra id")
    axes[0,0].set(xlabel="Time (s)",ylabel="Short-axis component",title="Median-r RBC: p_x, p_y, p_z")
    axes[0,1].set(xlabel="Time (s)",ylabel="|norm(q)−1|",title="All five shapes: unit quaternion")
    axes[1,0].set(xlabel="Time (s)",ylabel="|Ω_RBC| (s^-1)",title="Geometry-dependent angular velocity")
    axes[1,1].set(xlabel="Time (s)",ylabel="Tetra id (dimensionless)",title="Gradient transitions at tetra changes")
    for ax in axes.ravel():ax.grid(alpha=.2);ax.legend(fontsize=9)
    fig.subplots_adjust(bottom=.15,hspace=.38,wspace=.28)
    save(fig,9,"Finest VALIDATION_ONLY replay. All p/q/Ω finite; no gradient, vorticity or Ω smoothing.\nP0 gradients are constant within each tetra and can jump across faces; tetra id is not a physical scalar.",
         "四个面板分别看短轴分量、四元数长度误差、角速度和单元编号。所有值有限；角速度有细碎变化，与原始 tetra 梯度切换一致，没有为美观做平滑。",
         sources+["08_real_sweep.json"])


def stage10():
    fig,axes=plt.subplots(2,2,figsize=(15,10));comparisons=[];summaries=[]
    for g in range(5):
        reference=real(g,2);rt=col(reference,"time_s");rp=vectors(reference,"p_")
        shape_r=float(reference[0]["r"]);label=f"r={shape_r:.3f}"
        for d in range(3):
            rows=real(g,d);t=col(rows,"time_s");p=vectors(rows,"p_");keep=t<=rt[-1]
            interp=np.column_stack([np.interp(t[keep],rt,rp[:,j]) for j in range(3)])
            gap=shape_axis_distance(p[keep],interp)
            if d==2:gap=np.zeros(len(t))  # the reference compared with itself, by definition
            if d<2:axes[0,d].plot(t[keep],gap*180/np.pi,color=COLORS[g],label=label,lw=1)
            summary=dict(geometry_index=g,r=shape_r,dt_index=d,validation_dt_s=float(rows[0]["validation_dt_s"]),
                         max_axis_difference_rad=float(np.max(gap)),median_axis_difference_rad=float(np.median(gap)),
                         metric="acos(abs(dot(normalized p, normalized interpolated reference p)))",
                         reference="finest validation replay; not analytic truth",
                         includes_translation_path_difference=True,production_particle_timestep_frozen=False)
            summaries.append(summary)
            comparisons.extend(dict(geometry_index=g,dt_index=d,time_s=float(a),axis_difference_rad=float(b)) for a,b in zip(t[keep],gap))
        selected=[s for s in summaries if s["geometry_index"]==g]
        dt=np.array([s["validation_dt_s"] for s in selected])*1e6
        axes[1,0].plot(dt,[s["max_axis_difference_rad"]*180/np.pi for s in selected],"o-",color=COLORS[g],label=label)
        axes[1,1].plot(dt,[s["median_axis_difference_rad"]*180/np.pi for s in selected],"o-",color=COLORS[g],label=label)
    for ax,title in zip(axes[0],["Coarse vs finest validation replay","Middle vs finest validation replay"]):
        ax.set(xlabel="Time (s)",ylabel="Unoriented shape-axis difference (degrees)",title=title)
    axes[1,0].set(xlabel="Validation dt (µs)",ylabel="Maximum axis difference (degrees)",title="Measured max differences; finest is reference")
    axes[1,1].set(xlabel="Validation dt (µs)",ylabel="Median axis difference (degrees)",title="Measured median differences; finest is reference")
    for ax in axes.ravel():ax.grid(alpha=.2);ax.legend(fontsize=9)
    write_csv("10_axis_comparison.csv",comparisons);write_json("10_axis_comparison_summary.json",summaries)
    fig.suptitle("VALIDATION-ONLY TIMESTEP STUDY\nNOT PRODUCTION TIMESTEP SELECTION",fontsize=18)
    fig.subplots_adjust(bottom=.15,top=.85,hspace=.4,wspace=.29)
    save(fig,10,"Distance uses acos(|p1·p2|): p and −p are the same spheroid shape. No raw quaternion-component subtraction.\nDifferences include center-path changes and original gradient jumps. No production-convergence claim; no FEM timestep study.",
         "将每种形状的粗、中步长结果与最细验证重放在共同时间比较。使用 p 与 −p 等价的轴夹角，不比较原始 q 分量。报告最大和中位差异，包含中心路径变化，不宣称正式收敛或选定 production dt。",
         ["10_axis_comparison.csv","10_axis_comparison_summary.json","08_real_sweep.json"]+[f"08_real_g{g}_dt{d}.csv" for g in range(5) for d in range(3)])


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--stage",choices=["all"]+[str(i) for i in range(11)],default="all");args=parser.parse_args()
    FIGURES.mkdir(parents=True,exist_ok=True)
    for index in range(11) if args.stage=="all" else [int(args.stage)]:
        print(f"Generating Particle-2 figure {index:02d}",flush=True);globals()[f"stage{index}"]()


if __name__=="__main__":main()
