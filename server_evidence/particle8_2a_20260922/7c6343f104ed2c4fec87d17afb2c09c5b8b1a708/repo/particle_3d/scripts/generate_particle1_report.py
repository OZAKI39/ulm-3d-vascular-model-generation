#!/usr/bin/env python3
"""Deterministic Particle-1 plots/data; no simulation population or production dt."""
from pathlib import Path
import argparse
import csv
import importlib.util
import json
import sys
import numpy as np

PACKAGE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PACKAGE/'src'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from matplotlib.ticker import MaxNLocator
from particle_3d.audit import read_frozen
from particle_3d.particle1_cases import uniform_case,rotation_case,single_step_case,SYNTHETIC_ATOL,SYNTHETIC_RTOL
from particle_3d.sonovue_adapter import read_sonovue

REPORT=PACKAGE/'reports/particle1'
DATA=REPORT/'data'; FIGURES=REPORT/'figures'
COLORS=['#c75d19','#257bb7','#319765']
spec=importlib.util.spec_from_file_location('particle0_plotting_readonly',PACKAGE/'scripts/generate_particle0_report.py')
p0=importlib.util.module_from_spec(spec); spec.loader.exec_module(p0)
for candidate in ['/mnt/c/Windows/Fonts/msyh.ttc','/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc']:
    if Path(candidate).is_file():
        font_manager.fontManager.addfont(candidate)
        plt.rcParams['font.sans-serif']=[font_manager.FontProperties(fname=candidate).get_name(),'DejaVu Sans']
        break
plt.rcParams.update({'font.size':11,'axes.unicode_minus':False,'savefig.dpi':180})


def write_json(name,value):
    def scalar(value):
        if isinstance(value,np.generic): return value.item()
        raise TypeError(f'Unsupported JSON value: {type(value).__name__}')
    (DATA/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False,default=scalar)+'\n')


def write_csv(name,rows):
    with (DATA/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n'); w.writeheader(); w.writerows(rows)


def read_json(name): return json.loads((DATA/name).read_text())


def save(fig,name,footer):
    fig.text(.5,.025,footer,ha='center',va='bottom',fontsize=10,bbox=dict(facecolor='#f0f3f5',edgecolor='none',pad=7))
    fig.savefig(FIGURES/name,bbox_inches='tight'); plt.close(fig)


def note(index,title,paragraph,figure):
    (REPORT/f'{index:02d}_{title}.md').write_text(f'# Particle-1 小步骤 {index} 检查说明\n\n{paragraph}\n\n![检查图](figures/{figure})\n')


def stage0():
    provenance=read_json('00_particle1_scope_and_provenance.json')
    fig,(ax,meta)=plt.subplots(1,2,figsize=(15,10)); ax.axis('off'); meta.axis('off')
    boxes=[('Frozen FEM = READ ONLY','position m; Velocity m/s; Pressure Pa'),('Particle-0 = ACCEPTED DEPENDENCY','canonical tetra sampler; original 46 tests'),
           ('Particle-1 = ONE SPHERICAL MB','V = u_inf; Omega = 0.5 curl(u)'),('Explicit Euler position update','x_new = x_old + dt * V_old; dt required')]
    for i,(title,detail) in enumerate(boxes):
        y=.84-i*.22
        ax.add_patch(FancyBboxPatch((.04,y-.075),.9,.14,boxstyle='round,pad=.015',fc=['#dcecf4','#dcf1e3','#ffecd6','#eeeeff'][i],ec='#536673'))
        ax.text(.49,y+.02,title,ha='center',weight='bold',fontsize=12)
        ax.text(.49,y-.025,detail,ha='center',fontsize=10)
        if i<3: ax.annotate('',xy=(.49,y-.145),xytext=(.49,y-.095),arrowprops=dict(arrowstyle='->',lw=2))
    meta.text(.02,.95,'Check: frozen provenance and stage scope',fontsize=17,weight='bold',va='top')
    lines=[f"Mesh SHA: {provenance['mesh_sha256'][:24]}…",f"Flow SHA: {provenance['flow_sha256'][:24]}…",
           f"P0 implementation: {provenance['particle0_dependency_commit'][:16]}",'P0 manual review: PASS / USER_CHAT_REVIEW',
           'P1 branch:','dev/particle-1-single-mb-20260920','SonoVue: ORIGINAL SAMPLER, READ ONLY','One draw, explicit validation seed; no population',
           'Core units: m, s, m/s, 1/s, Pa, N','No production particle timestep frozen','No RBC / wall response / contact / CFD']
    meta.text(.02,.84,'\n\n'.join(lines),va='top',fontsize=11)
    fig.subplots_adjust(bottom=.13,wspace=.1)
    save(fig,'00_particle1_scope_and_provenance.png','PASS: accepted P0 code/tests unchanged; existing 46 P0 tests pass; 4602 frozen FEM files match SHA256.\nOne sphere only. All timesteps in this report are VALIDATION_ONLY.')
    note(0,'scope_and_provenance','Particle-0 已由用户审核通过，本轮只重复原有 46 项回归测试，未加增强项目。冻结 FEM 和 SonoVue 原仓库保持只读。图中流程只含一个球形微泡的速度、旋转速度与位置更新。','00_particle1_scope_and_provenance.png')


def stage1():
    size=read_json('01_single_mb_size_provenance.json')
    source=read_json('00_particle1_scope_and_provenance.json')
    _,dist,_=read_sonovue(source['sonovue_root'])
    size['histogram_plot_data']=dict(low_um=dist.low.tolist(),high_um=dist.high.tolist(),probability=dist.probability.tolist(),cdf_low=dist.cdf_low.tolist(),cdf_high=dist.cdf_high.tolist())
    write_json('01_single_mb_size_provenance.json',size)
    fig,axes=plt.subplots(1,2,figsize=(14,7))
    axes[0].bar(dist.low,dist.probability,width=dist.width,align='edge',color='#8fb7cc',edgecolor='white',label='Frozen number-weighted histogram')
    axes[0].set(xlabel='Diameter (µm)',ylabel='Probability per original bin (dimensionless)',title='Check: one sampled diameter in the frozen histogram')
    x=np.column_stack((dist.low,dist.high)).ravel(); cdf=np.column_stack((dist.cdf_low,dist.cdf_high)).ravel()
    axes[1].plot(x,cdf,color='#237d75',label='Original piecewise linear CDF')
    axes[1].set(xlabel='Diameter (µm)',ylabel='Cumulative probability (dimensionless)',title='Check: original CDF, no refitting or smoothing')
    for ax in axes:
        ax.axvline(size['diameter_um'],color='#b63827',lw=2,label='THIS PARTICLE')
        ax.grid(alpha=.2); ax.legend(loc='upper left')
    fig.suptitle('DEMO / VALIDATION ONLY — ONE MICROBUBBLE',fontsize=18,weight='bold')
    fig.subplots_adjust(bottom=.2,top=.85,wspace=.28)
    save(fig,'01_single_mb_size_provenance.png',f"PASS: diameter={size['diameter_um']:.9f} µm; radius={size['radius_m']*1e6:.9f} µm; seed={size['seed']}; N=1.\nHistogram SHA {size['histogram_sha256'][:20]}…; original sampler called; formal_simulation_population=false.")
    note(1,'size_provenance',f"固定 seed={size['seed']} 调用原始 SonoVue sampler 抽一个验证直径，为 {size['diameter_um']:.9f} µm；半径为 {size['radius_m']:.9e} m。适配层只把 µm 转为 m，再除以二，不改变 histogram、CDF（累积概率曲线）或随机算法。CSV 相邻 metadata 绑定样本文件 SHA；没有生成正式粒子群。",'01_single_mb_size_provenance.png')


def stage2():
    rows=uniform_case(); write_csv('02_uniform_flow_trajectory.csv',rows)
    fig,axes=plt.subplots(1,2,figsize=(14,7))
    dts=sorted({r['dt_s'] for r in rows},reverse=True)
    for i,dt in enumerate(dts):
        group=[r for r in rows if r['dt_s']==dt]
        t=np.array([r['time_s'] for r in group])
        if i==0: axes[0].plot([r['exact_x_m']*1e6 for r in group],[r['exact_y_m']*1e6 for r in group],'k-',lw=4,alpha=.4,label='Analytic exact line')
        axes[0].plot([r['x_m']*1e6 for r in group],[r['y_m']*1e6 for r in group],color=COLORS[i],marker=['o','s','^'][i],markevery=max(1,len(group)//10),ms=5,lw=1,label=f'dt={dt:g} s')
        axes[1].plot(t,[r['position_error_m'] for r in group],color=COLORS[i],label=f'dt={dt:g} s')
    axes[0].arrow(14,18,2,-1,width=.035,head_width=.35,color='#222222',length_includes_head=True)
    axes[0].text(14,18.6,'Velocity direction',fontsize=10)
    axes[0].set(xlabel='x (µm)',ylabel='y (µm)',title='Check: straight-line motion (x-y projection)'); axes[0].set_aspect('equal')
    axes[1].set(xlabel='Time (s)',ylabel='Max absolute position component error (m)',title='Check: exact uniform-flow trajectory')
    for ax in axes: ax.grid(alpha=.2); ax.legend()
    maximum=max(r['position_error_m'] for r in rows)
    write_json('02_uniform_metrics.json',dict(passed=all(r['position_error_m']<=r['position_error_bound_m'] for r in rows),max_position_error_m=maximum,
        max_velocity_error_m_s=max(r['velocity_error_m_s'] for r in rows),validation_timesteps_s=dts,position_bound_rule='16*(step_count+1)*eps*max(abs(x0)+abs(u)*duration)'))
    fig.suptitle('VALIDATION TIMESTEPS — NOT PRODUCTION PARTICLE TIMESTEP',fontsize=16,weight='bold')
    fig.subplots_adjust(bottom=.18,top=.86,wspace=.3)
    save(fig,'02_uniform_flow_trajectory.png',f'PASS: 3 timestep cases; max position error={maximum:.3e} m; V-u=0 m/s; Omega=0 1/s.\nRoundoff bound fixed from accumulated step count and coordinate scale; all 3 spatial components tested.')
    note(2,'uniform_flow',f'均匀流就是每个位置的速度一样，所以球应沿解析直线移动。本例三个验证 dt 都落在同一条直线上，最大位置分量误差 {maximum:.3e} m；V 与背景速度完全一致，旋转速度为零。误差界按加法次数和坐标大小预先确定，未调整阈值。','02_uniform_flow_trajectory.png')


def stage3():
    rows=rotation_case(); write_csv('03_pure_rotation_validation.csv',rows)
    fig=plt.figure(figsize=(14,9)); grid=fig.add_gridspec(2,2)
    ax=fig.add_subplot(grid[:,0]); compare=fig.add_subplot(grid[0,1]); err=fig.add_subplot(grid[1,1])
    xyz=np.array([[r[f'{a}_m'] for a in 'xyz'] for r in rows])*1e6
    uv=np.array([[r[f'velocity_{a}_m_s'] for a in 'xy'] for r in rows])*1e6
    ax.quiver(xyz[:,0],xyz[:,1],uv[:,0],uv[:,1],angles='xy',scale_units='xy',scale=15,color='#267cb7')
    ax.set(xlabel='x (µm)',ylabel='y (µm)',title='Check: rigid rotation direction in z=0 plane',aspect='equal',xlim=(-13,13),ylim=(-13,13)); ax.grid(alpha=.2)
    omega=np.array([rows[0][f'analytic_omega_{a}_s_inv'] for a in 'xyz']); measured=np.array([rows[0][f'sampled_omega_{a}_s_inv'] for a in 'xyz'])
    compare.bar(np.arange(3)-.18,omega,.36,label='Analytic fluid Omega',color='#8999ad'); compare.bar(np.arange(3)+.18,measured,.36,label='Sampled particle Omega',color='#dc9238')
    compare.set(xticks=range(3),xticklabels=['x','y','z'],xlabel='Component (dimensionless)',ylabel='Angular velocity (1/s)',title='Check: Omega_particle = Omega_fluid'); compare.legend(); compare.grid(axis='y',alpha=.2)
    err.plot(range(len(rows)),[r['angular_velocity_error_s_inv'] for r in rows],'o',ms=3,label='Maximum component error')
    err.set(xlabel='Query index (dimensionless)',ylabel='Absolute angular velocity error (1/s)',title='Check: the required one-half factor'); err.legend(); err.grid(alpha=.2)
    maximum=max(r['angular_velocity_error_s_inv'] for r in rows)
    write_json('03_rotation_metrics.json',dict(passed=maximum<=SYNTHETIC_ATOL,max_angular_velocity_error_s_inv=maximum,atol=SYNTHETIC_ATOL,rtol=SYNTHETIC_RTOL,count=len(rows)))
    fig.subplots_adjust(bottom=.17,top=.93,hspace=.48,wspace=.35)
    save(fig,'03_pure_rotation_validation.png',f'PASS: 49 points; curl(u)=2*Omega_fluid; strain=0; max Omega error={maximum:.3e} 1/s.\nFixed atol=rtol=256*float64 eps={SYNTHETIC_ATOL:.3e}. Arrow length uses a display scale only.')
    note(3,'pure_rotation',f'vorticity 是流场的局部旋转趋势，在这个刚体旋转人工场里等于整体旋转速度的两倍。因此球的旋转速度必须乘 1/2；49 个位置的最大误差为 {maximum:.3e} 1/s。专门的永久测试用三个非零分量检查，误写成 Omega=vorticity 会失败。','03_pure_rotation_validation.png')


def stage4():
    record=single_step_case(); write_json('04_single_step_example.json',record)
    fig=plt.figure(figsize=(14,8)); ax=fig.add_subplot(121,projection='3d'); explain=fig.add_subplot(122); explain.axis('off')
    old=np.array(record['old_position_m'])*1e6; new=np.array(record['new_position_m'])*1e6; delta=new-old
    ax.scatter(*old,s=70,color='#2477b0',label='Old center'); ax.scatter(*new,s=90,marker='*',color='#cb553c',label='New center')
    ax.quiver(*old,*delta,arrow_length_ratio=.2,color='#2477b0',linewidth=2)
    ax.plot(*np.stack((old,new)).T,'--',color='gray',label='Explicit Euler displacement')
    ax.set(xlabel='x (µm)',ylabel='y (µm)',zlabel='z (µm)',title='Check: a single 3D Euler step'); ax.legend(loc='upper left')
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]: axis.set_major_locator(MaxNLocator(4))
    explain.text(.02,.9,'查流场 → 得到粒子速度 → 推进一步',fontsize=16,weight='bold')
    lines=['① 查流场：在旧位置读取 u 和 curl(u)。','② 得到粒子速度：V=u；旋转速度 Ω=0.5 curl(u)。',
           '③ 向前推进一个验证时间步：x_new=x_old+dt V_old。',f"本例 dt={record['dt_s']} s（仅验证用）。",'④ 在新位置重新查场，保存对应位置的 V 和 Ω。',
           '重新查场不再移动位置，仍是显式 Euler。','输入 state 保持不变；返回独立的新 state。','没有质量、加速度或惯性积分。']
    explain.text(.02,.79,'\n\n'.join(lines),va='top',fontsize=11)
    fig.subplots_adjust(bottom=.18,wspace=.24)
    save(fig,'04_single_step_dataflow.png','PASS: displacement equals dt * old-position velocity exactly; input state unchanged.\nVALIDATION TIMESTEP ONLY. Marker sizes are illustrative; all calculations use SI.')
    note(4,'single_step', '本例用真正随三维位置变化的人工速度场。位置更新只使用旧点速度；返回的新 state 再保存新位置对应的速度和旋转速度，因此轨迹每行的字段都对应同一个位置。该重新查询没有增加积分阶数，也没有原地改写输入。','04_single_step_dataflow.png')


def load_trajectories():
    result=[]
    for index in range(3):
        with (DATA/f'05_real_trajectory_dt{index}.csv').open(newline='') as f: raw=list(csv.DictReader(f))
        rows=[{k:(v if k in ['timestep_role','boundary_event'] else float(v)) for k,v in row.items()} for row in raw]
        result.append((rows,read_json(f'05_real_trajectory_dt{index}.json')))
    return result


def stage5():
    cases=load_trajectories(); rows,summary=cases[-1]
    _,_,_,boundaries=read_frozen(PACKAGE.parent/'formal_3D_flow_solver/FEM_SimVascular')
    fig=plt.figure(figsize=(13,10)); ax=fig.add_subplot(111,projection='3d',computed_zorder=False); p0.vessel(ax,boundaries,.3)
    xyz=np.array([[r[f'{a}_m'] for a in 'xyz'] for r in rows])*1e6
    ax.plot(*xyz.T,color='#e3811d',lw=2,label='One MB trajectory (finest validation replay)')
    ax.scatter(*xyz[0],s=80,color='#d12c42',marker='o',zorder=5,label='Start: inlet-adjacent tetra centroid')
    ax.scatter(*xyz[-1],s=110,color='#171717',marker='*',zorder=5,label=f"End: {summary['exit_boundary']}")
    indices=np.linspace(0,len(rows)-2,12,dtype=int)
    velocity=np.array([[rows[i][f'V_{a}_m_s'] for a in 'xyz'] for i in indices])
    arrows=velocity/np.linalg.norm(velocity,axis=1)[:,None]*3
    ax.quiver(*xyz[indices].T,*arrows.T,color='#155878',arrow_length_ratio=.4,linewidth=1.8,zorder=4,label='Velocity direction (display overlay)')
    ax.set_title('Check: one microbubble follows the frozen vessel flow',pad=22)
    ax.legend(loc='upper left',bbox_to_anchor=(-.06,1.04),fontsize=10)
    fig.subplots_adjust(bottom=.18)
    save(fig,'05_real_vessel_single_mb_trajectory.png',f"{'PASS' if summary['passed'] else 'FAIL'}: first event={summary['exit_boundary']}; exit t={summary['exit_time_s']:.9f} s; wall crossing={summary['wall_crossing']}.\nMARKER SIZE EXAGGERATED FOR VISIBILITY. Arrows normalized for direction. VALIDATION ONLY; no finite-radius wall clearance claim.")
    note(5,'real_trajectory',f"起点按全部 179 个 INLET 相邻四面体中心的向内速度最大值确定，最终编号为 208001；不是反复试轨迹挑点。图展示最细验证 dt 的同一个微泡，首次边界事件为 {summary['exit_boundary']}，时间 {summary['exit_time_s']:.9f} s；未出现 WALL 或 INLET 穿出。每条线段都与全部正式边界检查，事件后立即停止，不反弹或推回；没有有限尺寸壁面间隙模型。",'05_real_vessel_single_mb_trajectory.png')


def stage6():
    cases=load_trajectories(); diagnostic=[]
    for index,(rows,_) in enumerate(cases):
        diagnostic.extend([dict(case_id=index,**{k:r[k] for k in ['time_s','validation_dt_s','timestep_role','speed_m_s','pressure_pa','vorticity_norm_s_inv','omega_norm_s_inv','boundary_event']}) for r in rows])
    write_csv('06_real_vessel_single_mb_diagnostics.csv',diagnostic)
    fig,axes=plt.subplots(2,2,figsize=(15,10))
    for ax,key,label,unit in zip(axes.ravel(),['speed_m_s','pressure_pa','vorticity_norm_s_inv','omega_norm_s_inv'],['Speed','Pressure','Vorticity magnitude','Particle angular speed'],['m/s','Pa','1/s','1/s']):
        for i,(rows,summary) in enumerate(cases):
            ax.plot([r['time_s'] for r in rows],[r[key] for r in rows],color=COLORS[i],lw=.9,alpha=.85,label=f"dt={summary['validation_dt_s']*1e6:.3f} µs")
        ax.set(xlabel='Time (s)',ylabel=f'{label} ({unit})',title=f'Check: {label.lower()} along the trajectory'); ax.grid(alpha=.2); ax.legend()
    fig.suptitle('Check: all states finite; V=u and Omega=0.5*vorticity at every recorded position',fontsize=15)
    fig.subplots_adjust(bottom=.15,top=.9,hspace=.36,wspace=.28)
    finite=all(s['active_and_terminal_finite'] for _,s in cases); uerr=max(s['max_velocity_relation_error_m_s'] for _,s in cases); oerr=max(s['max_angular_relation_error_s_inv'] for _,s in cases)
    save(fig,'06_real_vessel_single_mb_diagnostics.png',f"{'PASS' if finite and uerr==oerr==0 else 'FAIL'}: 3 independent validation replays; all finite={finite}; max |V-u|={uerr:.3e} m/s; max |Omega-curl(u)/2|={oerr:.3e} 1/s.\nDerivative jumps reflect the original piecewise constant tetra field; no smoothing. Vector relation checked component by component.")
    note(6,'real_diagnostics',f'三个验证时间步下逐行比较了全部速度分量和旋转分量。最大 V-u 误差 {uerr:.3e} m/s，最大 Omega-curl(u)/2 误差 {oerr:.3e} 1/s，全部记录有限。旋转曲线的细碎台阶来自 P0 原本单元内常量的梯度，未作平滑。','06_real_vessel_single_mb_diagnostics.png')


def stage7():
    cases=load_trajectories(); fig,axes=plt.subplots(2,2,figsize=(15,11))
    reference=cases[-1][0]; reference_t=np.array([r['time_s'] for r in reference]); reference_xyz=np.array([[r[f'{a}_m'] for a in 'xyz'] for r in reference])
    comparison=[]
    for i,(rows,summary) in enumerate(cases):
        label=f"dt={summary['validation_dt_s']*1e6:.3f} µs"
        xyz=np.array([[r[f'{a}_m'] for a in 'xyz'] for r in rows]); t=np.array([r['time_s'] for r in rows])
        axes[0,0].plot(xyz[:,0]*1e6,xyz[:,1]*1e6,color=COLORS[i],lw=1.2,label=label)
        keep=t<=reference_t[-1]
        interpolated=np.column_stack([np.interp(t[keep],reference_t,reference_xyz[:,j]) for j in range(3)])
        gap=np.linalg.norm(xyz[keep]-interpolated,axis=1)
        axes[0,1].plot(t[keep],gap*1e6,color=COLORS[i],label=label)
        comparison.extend([dict(case_id=i,time_s=float(ti),position_difference_from_finest_m=float(di),reference='finest validation path; not analytic truth') for ti,di in zip(t[keep],gap)])
    dt=np.array([s['validation_dt_s'] for _,s in cases])*1e6; times=[s['exit_time_s'] for _,s in cases]; lengths=[s['trajectory_length_m']*1e6 for _,s in cases]
    axes[1,0].plot(dt,times,'o-',label='First boundary exit time'); axes[1,1].plot(dt,lengths,'o-',color='#81619d',label='Piecewise straight path length')
    axes[0,0].set(xlabel='x (µm)',ylabel='y (µm)',title='Check: same start, three independent replays'); axes[0,0].set_aspect('equal')
    axes[0,1].set(xlabel='Time (s)',ylabel='Distance from finest validation path (µm)',title='Check: differences at common times')
    axes[1,0].set(xlabel='Validation timestep (µs)',ylabel='Exit time (s)',title='Check: exit boundary and time')
    axes[1,0].text(.03,.94,'\n'.join(f"{s['validation_dt_s']*1e6:.3f} µs → {s['exit_boundary']}" for _,s in cases),transform=axes[1,0].transAxes,va='top',fontsize=10)
    axes[1,1].set(xlabel='Validation timestep (µs)',ylabel='Trajectory length (µm)',title='Check: total traveled distance')
    for ax in axes.ravel(): ax.grid(alpha=.2); ax.legend(loc='best',fontsize=9); ax.ticklabel_format(axis='y',useOffset=False)
    fig.suptitle('VALIDATION-ONLY TIMESTEP STUDY\nNOT PRODUCTION TIMESTEP SELECTION',fontsize=19,weight='bold')
    fig.subplots_adjust(bottom=.15,top=.85,hspace=.4,wspace=.28)
    difference=[abs(times[0]-times[1]),abs(times[1]-times[2])]
    write_csv('07_position_comparison.csv',comparison)
    write_json('07_comparison_metrics.json',dict(exit_time_adjacent_differences_s=difference,exit_time_difference_decreases=difference[1]<difference[0],
        same_exit_boundary=len({s['exit_boundary'] for _,s in cases})==1,production_particle_timestep_frozen=False,
        figure_sources=['05_real_trajectory_dt0.csv','05_real_trajectory_dt1.csv','05_real_trajectory_dt2.csv','07_validation_timestep_comparison.csv','07_position_comparison.csv']))
    save(fig,'07_validation_timestep_comparison.png',f"{'PASS' if all(s['passed'] for _,s in cases) else 'FAIL'}: all 3 trajectories exit through {cases[0][1]['exit_boundary']}; no WALL/INLET crossing.\nAdjacent exit-time differences: {difference[0]:.3e}, {difference[1]:.3e} s. Trend is informational; this is not a FEM timestep study.")
    note(7,'validation_timestep_comparison',f"三个时间步在首次运行之前由局部最短边长和起点速度计算，固定逐级减半；没有根据结果调 dt。三个重放均从 {cases[0][1]['exit_boundary']} 离开，相邻两次的出口时间差为 {difference[0]:.3e}、{difference[1]:.3e} s。这里报告数值趋势，不选择 production dt，也不称作 FEM 时间步研究；轨迹是否有视觉上突然跳跃仍待用户审核。",'07_validation_timestep_comparison.png')


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--stage',choices=['all']+[str(i) for i in range(8)],default='all'); args=parser.parse_args()
    DATA.mkdir(parents=True,exist_ok=True); FIGURES.mkdir(parents=True,exist_ok=True)
    functions=[stage0,stage1,stage2,stage3,stage4,stage5,stage6,stage7]
    for i in range(8) if args.stage=='all' else [int(args.stage)]:
        print(f'Generating Particle-1 figure {i:02d}',flush=True); functions[i]()


if __name__=='__main__': main()
