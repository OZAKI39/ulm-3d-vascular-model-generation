"""Six PNG/PDF review figures and a twelve-section, data-generated review."""
import hashlib,json,re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'verification/mechanical_equilibrium_v1'
V=OLD/'adr_completion';R=ROOT/'results/streaming_erosion_equilibrium_v1';FIG=V/'figures'
COLORS=['#60ccf0','#ffb454','#b8e986','#e28bb6']


def read(p):return json.loads(p.read_text())


def write(p,text):
    with p.open('x') as f:f.write(text)


def summary(p):return read(p/'SUMMARY.json')


def history(p):
    with np.load(p/'history.npz') as z:return dict(z)


def save(fig,name):
    for ext in ['png','pdf']:fig.savefig(FIG/f'{name}.{ext}',dpi=180,facecolor='black',bbox_inches='tight',pad_inches=.2)
    plt.close(fig)


def style():
    font_manager.fontManager.addfont('/mnt/c/Windows/Fonts/arial.ttf')
    plt.rcParams.update({'font.family':'Arial','font.weight':'normal','font.size':13,'axes.titlesize':18,
        'axes.titleweight':'normal','axes.labelweight':'normal','figure.facecolor':'black','axes.facecolor':'black',
        'text.color':'white','axes.labelcolor':'white','xtick.color':'white','ytick.color':'white',
        'axes.edgecolor':'#64778b','grid.color':'#253447','axes.grid':True,'legend.frameon':False,'pdf.fonttype':42})


def plots():
    small=summary(V/'cases/small_adr')['final'];ref=summary(OLD/'cases/small_scipy')['final']
    base=R/'dt_erosion_0p01';coarse=summary(base/'equilibrium_000')['final']
    fig,axs=plt.subplots(1,2,figsize=(15,6),layout='constrained')
    for ax,size,value,title in zip(axs,['small','coarse'],[small,coarse],['Small clot: 36 particles','Baseline clot: 360 particles']):
        for ms,col in zip([2,4,8,16],COLORS):
            h=history(OLD/f'cases/{size}_legacy_{ms}ms')
            ax.plot(h['mechanical_time_s']*1000,h['maximum_epsilon_plus'],color=col,lw=1.1,label=f'{ms} ms window')
            ax.plot(ms,h['maximum_epsilon_plus'][-1],'o',color=col,ms=6)
        ax.axhline(value['maximum_epsilon_plus'],color='white',lw=1.8,label='ADR equilibrium')
        if size=='small':ax.axhline(ref['maximum_epsilon_plus'],color='#ef6f89',ls='--',lw=1.4,label='SciPy reference')
        ax.set(title=title,xlabel='Mechanical integration time (ms)',ylabel='Maximum tensile strain')
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax.legend(fontsize=10,loc='lower right')
    fig.suptitle('Fixed window convergence',fontsize=23);save(fig,'fixed_window_convergence')

    fig,ax=plt.subplots(figsize=(12,6),layout='constrained')
    for p,label,col in [(V/'cases/small_adr/equilibrium','Small clot',COLORS[0]),(base/'equilibrium_000','Baseline clot',COLORS[1])]:
        h=history(p);ax.semilogy(h['equilibrium_iteration'],h['relative_force_residual'],label=label,color=col,lw=1.6)
    ax.axhline(1e-4,color='white',ls='--',label='Acceptance: 1e-4')
    ax.set(title='Force residual history',xlabel='ADR iteration (not exposure time)',ylabel='Relative free-particle force residual')
    ax.legend();save(fig,'force_residual_history')

    fig,axs=plt.subplots(1,2,figsize=(15,6),layout='constrained')
    rows=[small,ref]+[summary(OLD/f'cases/small_legacy_{ms}ms')['final'] for ms in [2,4,8,16]]
    labels=['ADR','SciPy','2 ms','4 ms','8 ms','16 ms']
    for ax,key,scale,title,ylabel in [(axs[0],'maximum_displacement_m',1e9,'Maximum displacement','Displacement (nm)'),
                                      (axs[1],'maximum_epsilon_plus',1,'Maximum tensile strain','Tensile strain')]:
        ax.bar(labels,[r[key]*scale for r in rows],color=['white','#ef6f89',*COLORS]);ax.set(title=title,ylabel=ylabel)
    axs[1].ticklabel_format(axis='y',style='sci',scilimits=(0,0));fig.suptitle('ADR versus reference: same small clot',fontsize=23)
    save(fig,'adr_vs_reference')

    h=history(base/'equilibrium_000');t=h['equilibrium_iteration'];fig,axs=plt.subplots(1,3,figsize=(18,6.4),layout='constrained')
    axs[0].semilogy(t[1:],np.maximum(h['adr_pseudo_kinetic_measure'][1:],1e-40),color=COLORS[1]);axs[0].set(title='Numerical motion diagnostic',ylabel='Pseudo kinetic measure (not physical energy)')
    axs[1].semilogy(t,h['relative_force_residual'],color=COLORS[0]);axs[1].axhline(1e-4,color='white',ls='--');axs[1].set(title='Free-particle residual',ylabel='Relative force residual')
    axs[2].plot(t,h['maximum_epsilon_plus'],color=COLORS[2]);axs[2].set(title='Maximum tensile strain',ylabel='Tensile strain');axs[2].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    for ax in axs:
        ax.axvline(t[-1],color='white',ls=':',lw=1.4);ax.set_xlabel('ADR iteration');ax.set_xlim(0,t[-1]*1.15)
    fig.suptitle(f'Pre-relaxation: damage = 0 and exposure = 0 throughout\nEROSION CLOCK STARTS HERE: after iteration {int(t[-1])} passes every criterion',fontsize=19)
    save(fig,'pre_relaxation_driver')

    fig,ax=plt.subplots(figsize=(12,6.4),layout='constrained')
    for n,step,col,mark in zip(['0p02','0p01','0p005'],[.02,.01,.005],COLORS,['o','s','.']):
        rows=read(R/f'dt_erosion_{n}/history.json');assert all(r['equilibrium_passed'] for r in rows)
        ax.plot([r['erosion_exposure_time_s'] for r in rows],[r['mean_damage'] for r in rows],label=f'Exposure step: {step:g} s',color=col,marker=mark,ms=5,lw=1.4)
    ax.set(title='Erosion time step consistency',xlabel='Erosion exposure time (s)',ylabel='Mean particle damage')
    ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax.legend(loc='upper left')
    fig.suptitle('All mechanical states equilibrated before update',fontsize=19)
    save(fig,'erosion_time_step_consistency_equilibrated')

    z=np.load(base/'states.npz');x=z['x'][-1]*1000
    fig=plt.figure(figsize=(17,7.6),layout='constrained')
    for j,(name,values,cmap) in enumerate([('Tensile driver',z['particle_driver'][-1],'viridis'),('Particle damage',z['damage'][-1],'inferno')]):
        ax=fig.add_subplot(1,2,j+1,projection='3d')
        p=ax.scatter(x[:,0],x[:,1],x[:,2],c=values,s=37,cmap=cmap,vmin=float(values.min()),vmax=float(values.max()),depthshade=False,linewidths=0)
        ax.view_init(elev=23,azim=-63);ax.set_box_aspect(np.ptp(z['X'],axis=0))
        ax.set(xlabel='x (mm)',ylabel='y (mm)',zlabel='z (mm)',title=f'{name}\nActual min: {values.min():.3e}   max: {values.max():.3e}')
        for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
            axis.pane.fill=False;axis.pane.set_edgecolor((0,0,0,0));axis._axinfo['grid']['color']='#253447'
        cb=fig.colorbar(p,ax=ax,shrink=.65,pad=.03);cb.set_label('Strain' if j==0 else 'D (brightness does not mean D = 1)')
    fig.suptitle('Final equilibrium at 0.1 s exposure | Actual displacement scale = 1×',fontsize=21)
    save(fig,'final_driver_and_damage_maps')


def report():
    gate=read(V/'SMALL_REFERENCE_GATE.json');slow=read(V/'SLOW_STEP_COMPARISON.json');audit=read(V/'STATE_AUDIT.json')
    base=R/'dt_erosion_0p01';s=summary(base);pre=summary(base/'equilibrium_000');final=summary(base/'equilibrium_010')
    small=summary(V/'cases/small_adr');ref=summary(OLD/'cases/small_scipy');test=read(V/'TEST_SUMMARY.json')
    status=read(V/'FINAL_STATUS.json');hw=read(V/'HARDWARE.json');cg=gate['comparisons'][0]
    rows=[]
    for size in ['small','coarse']:
        for ms in [2,4,8,16]:
            a=summary(OLD/f'cases/{size}_legacy_{ms}ms')
            rows.append(f"| {36 if size=='small' else 360} | {ms} | {a['final']['maximum_epsilon_plus']:.10e} | {a['final']['relative_force_residual']:.4e} | {'PASS' if a['converged'] else 'NOT PASSED'} |")
    methods=[('ADR',small),('SciPy',ref)]+[(f'{ms} ms',summary(OLD/f'cases/small_legacy_{ms}ms')) for ms in [2,4,8,16]]
    comparison_rows=[]
    for name,ss in methods:
        f=ss['final'];comparison_rows.append(f"| {name} | {f['maximum_displacement_m']:.8e} | {f['maximum_epsilon_plus']:.8e} | {f['mean_epsilon_plus']:.8e} | {f['stored_energy_J']:.8e} | {f['surface_traction_mean_Pa']:.6f} / {f['surface_traction_max_Pa']:.6f} |")
    slow_rows=[]
    for r in slow['comparisons']:
        a=r['final'];slow_rows.append(f"| {r['case']} | {a['mean_damage']:.10e} | {a['maximum_damage']:.10e} | {a['maximum_bond_damage']:.10e} | {r['equilibrium_iterations']} | {max(r['maximum_relative_differences_at_shared_exposures'].values())*100:.7f}% |")
    energy_rows=[]
    for n in ['0p02','0p01','0p005']:
        ss=summary(R/f'dt_erosion_{n}');energy_rows.append(f"| {n} | {ss['damage_energy_increase_count']} / {ss['accepted_damage_updates']} | {ss['signed_damage_energy_drop_J']:.9e} |")
    files=['pd_clot/regularization/equilibrium_adr.py','pd_clot/regularization/erosion_equilibrium_runner.py',
        'configs/streaming_erosion_equilibrium_v1.json','configs/streaming_erosion_equilibrium_guarded_v1.json',
        'tests/test_equilibrium_adr.py','tests/test_erosion_pre_relaxation.py','scripts/verify_adr_completion.py',
        'scripts/audit_adr_completion.py','scripts/review_adr_completion.py']
    text=f'''# 1. 这一轮为什么做

之前开流后的第一段 2 ms 还在晃动，应变漂移约 8.48%。这时直接开始记损伤，会把启动过程混入侵蚀。
现在先让力学达到平衡，之后每次损伤更新前再次检查平衡。计算机迭代次数不再计入治疗曝光秒数。
**结论：本开发算例的机械平衡与慢步比较通过；损伤能量耗散问题仍未解决。**

本报告是 `adr_completion/` 新增的完成版。上级同名报告保留为当时缺论文的历史记录，不能作为最新结论。
原 BraVa CFD、微泡批次、原 cyclic fatigue、原 erosion_v1 及全部旧结果均未改动。

# 2. 我们采用什么已有方法

直接依据用户提供的 Kilic & Madenci (2010) 全文 §2.1，印刷页 196–197，公式 (11)–(14)、(18)–(19)。
已提取 11 页正文并渲染核对公式；PDF 路径和 SHA 在 `resources/PAPER_IDENTITY.json`。
[论文 DOI](https://doi.org/10.1016/j.tafmec.2010.08.001)。

原文 (14) 规定虚拟质量应覆盖刚度绝对行和。当前 NOSB 无现成切线矩阵，故对现有 free 总力作中央差分，
逐列累加行和，包括随位置变化的流体牵引。没有使用不适用于当前 NOSB 的论文键模型公式 (15)。
虚拟质量取行和界的 1.1 倍，差分扰动 h×10⁻⁵，结束时以半扰动重新检查行和界；没有修改力学 kernel。
这里采用 (14) 的坐标系相关形式，不声称具有论文基于 (15) 的框架不变性质。

首次严格按正文写出的 probe 出现负局部 Rayleigh 商，在迭代 16 后停止；完整失败数据保留在 `initial_small_probe/`。
随后核实 **Madenci & Oterkus (2014) 作者提供的 Fortran 附录**：非正商使阻尼取零；dtau=1 时 c>2 则使用1.9。
据此补上原作者保护分支，而非擅自取绝对值。来源、文件 SHA、行号和实现前声明见
`DAMPING_GUARD_ADDENDUM.md`、`resources/SPRINGER_CODE_IDENTITY.json`。
[Springer 官方附录](https://extras.springer.com/?query=978-1-4614-8464-6)。
正式标准小算例中负商保护 {small['negative_rayleigh_guard_count']} 次、上限保护 {small['damping_cap_guard_count']} 次；并不掩盖这些触发。

Peridigm 仅作成熟软件设计参照：核对了 BSD-3-Clause 许可、固定自由度及 QuasiStatic/NOX 残量设计、
correspondence 准静态示例；既未 vendor Peridigm/Trilinos，也未复制其求解器。
[源码固定版本](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/core/Peridigm.cpp)，
[2023 论文](https://doi.org/10.1007/s42102-023-00100-0)。
SciPy 直接调用现有 [least_squares](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html)，
只作为小问题独立求根算法；不是第二个材料模型，更不是实验验证。

# 3. 什么叫“平衡”

同时满足三项：自由粒子的未平衡总力足够小、数值运动量足够小、拉伸应变不再明显变化。
每 50 次迭代检查，连续 5 次通过：相对 free 总力 ≤10⁻⁴、末 50 步 driver drift ≤1%、
伪动能诊断/|储能| ≤10⁻⁶。分母下限分别为10⁻¹² N、10⁻¹² strain、10⁻²⁴ J，最大50000步。
所有收敛容差先声明再运行，首个失败后也没有放宽。

ADR 使用虚拟时间步1及虚拟质量，按论文中央差分更新位移，并以局部刚度估计阻尼。
它是 **NUMERICAL EQUILIBRATION DEVICE**，不能解释为 clot 黏度。
物理密度始终1056 kg/m³；pseudo half velocity 不当作物理速度输出，physical kinetic energy 为 null。
ADR 全过程 D、integrity、拓扑、曝光和周期数均冻结。具体公式、单位和局限见 `METHOD_AUDIT.md` 及保护补充。

# 4. 代码改了什么

**此前已有文件修改数：0。** 新增文件：

{chr(10).join('- `'+f+'`' for f in files)}

新接口支持 ADR，并转发旧 legacy 和小规模 SciPy 路径；旧 `equilibrium.py` 保留为当时缺论文的历史快照。
新 runner 每个阶段保存 x、D、力、逐迭代历史、虚拟质量及独立时钟。最后一次损伤更新后还做一次平衡，
仅用于最终状态和图，曝光仍停在0.1 s。任一未通过的阶段保存状态并停止损伤推进。

硬件：{hw['cpu_model']}，{hw['logical_cpu_count']}逻辑CPU，{hw['available_physical_cores']}可用物理核心，
RAM {hw['ram_total_bytes']/1024**3:.2f} GiB。当前内存预算使两个新 sweep 各用
{read(V/'SMALL_BATCH.json')['workers']} / {read(V/'EROSION_BATCH.json')['workers']} 个进程；每进程BLAS/OMP/MKL线程数1，
按1 GiB/worker并保留至少20%可用内存估算。逐job wall/CPU/peak RSS 已保存；无GPU重写。
旧窗口上一阶段已按冻结条件并行重跑，本次复用其不可变结果，没有重复耗费计算。

# 5. 测试

旧78项全部保留通过，前置独立参考阶段新增15项保留通过，本次新增19项通过。
**全量 {test['total_tests_passed']} passed，0失败、0跳过。**
覆盖了小问题解析解/独立参考、残余力下降、固定点、输入不变、质量界、作者保护分支、
额外迭代、静态导入场、旧路径逐元素复现、初始D和曝光为零、周期元数据独立、
初始和中途失败阻止后续损伤更新、同坐标 signed energy 重算。
日志见 `NEW_TESTS_FIRST.log`、`FULL_TESTS.log`。
另对所有 {sum(c['equilibrium_states'] for c in audit['cases'])} 个保存平衡状态重算了残余力，
对所有35次正式损伤跳变重算了能量及损伤增量，见 `STATE_AUDIT.json`。

# 6. 小算例独立验证

36粒子，12固定，72自由坐标，550 bonds；材料和mean flow均未调整。
SciPy：{ref['scipy_termination_reason']}，相对力残差 {ref['final']['relative_force_residual']:.7e}。
ADR：{small['iterations']}次迭代，相对力残差 {small['final']['relative_force_residual']:.7e}。
最大应变相对差 {cg['relative_max_driver_difference']:.7e}；RMS位移差/参考最大位移
{cg['rms_difference_over_reference_max']:.7e}，均小于预设1%。
RMS位移差 {cg['rms_displacement_difference_m']:.7e} m，最大位移差 {cg['maximum_displacement_difference_m']:.7e} m，
relative L2位移差 {cg['relative_L2_displacement_difference']:.7e}。

| 方法 | 最大位移(m) | 最大拉伸应变 | 平均拉伸应变 | 储能(J) | 牵引均值/最大值(Pa) |
|---|---|---|---|---|---|
{chr(10).join(comparison_rows)}

额外500次迭代和static resolved-field变体也通过；完整对照在 `SMALL_REFERENCE_GATE.json`。
小问题的参考值只用于同几何比较，不拿来替代360粒子问题的参考解。

# 7. Fixed-window sensitivity

全部从同一未变形状态、v=0、D=0出发，dt=1 µs、原阻尼3500 s⁻¹，不接续前一档构型。

| 粒子数 | 窗口(ms) | 最大拉伸应变 | 相对力残差 | 联合条件 |
|---|---|---|---|---|
{chr(10).join(rows)}

360粒子 ADR 无损伤平衡 driver={pre['final']['maximum_epsilon_plus']:.10e}。
2/4 ms 未通过的历史数据完整保留；8/16 ms趋向ADR平衡，图中小算例另有SciPy线。

# 8. Pre-relaxation 结果

360粒子基线从未变形状态开始，**{pre['iterations']}次ADR迭代**后通过全部门槛；这是虚拟迭代，不是{pre['iterations']}秒治疗。
最终相对残差 {pre['final']['relative_force_residual']:.7e}，L2={pre['final']['force_residual_L2_N']:.7e} N。
整段D=0、erosion exposure=0。只有连续5次联合检查通过之后，才允许第一次非零曝光损伤更新。
图 `pre_relaxation_driver` 明确标出 EROSION CLOCK STARTS HERE。

# 9. Equilibrated erosion rerun

新结果在 `results/streaming_erosion_equilibrium_v1/`，原v1结果不覆盖。
C_t=0.01 s⁻¹、m=1、epsilon_ref=10⁻⁴、threshold=0、D_break=1，mean flow及材料不变。
基线慢步0.01 s，10次更新，曝光0.1 s；更新后最终平衡不增加曝光。
最终粒子平均D={s['final']['mean_damage']:.10e}，最大粒子D={s['final']['maximum_damage']:.10e}，
最大bond D={s['final']['maximum_bond_damage']:.10e}；broken bonds={s['final']['broken_bond_count']}，
detached fraction={s['final']['detached_volume_fraction']:g}。
**没有观察到断裂或脱离**，这不是“模型失败”，也不能称为已验证断裂功能。

最终相对总力残差 {final['final']['relative_force_residual']:.7e}；
L2={final['final']['force_residual_L2_N']:.7e} N，RMS={final['final']['force_residual_RMS_N']:.7e} N，
max={final['final']['force_residual_max_N']:.7e} N。
现有ADR范围是锚定且拓扑不变的开发算例；遇到断裂导致范围改变会保存并停止，不冒充脱落后的物理动力学。

# 10. 时间步比较

三档曝光都严格结束于0.1 s，ADR迭代没有混入秒数。
比较所有共同曝光时刻的六项预声明指标：平均/最大粒子D、最大bond D、最大driver、最大位移、储能。
以0.005 s为参考，最大相对差 **{slow['maximum_relative_difference']*100:.7f}%**，低于预声明5%；
仅说明这三档、本算例的开发一致性，不构成网格客观性或一般时间收敛阶证明。

| 慢步(s，p表示小数点) | 平均粒子D | 最大粒子D | 最大bond D | ADR总迭代 | 六指标最大相对差 |
|---|---|---|---|---|---|
{chr(10).join(slow_rows)}

六组图均为黑底、Arial常规字重、英文、无箭头，PNG和PDF各一份：

- [Fixed window convergence](figures/fixed_window_convergence.png)：等更久后稳定到哪里。
- [Force residual history](figures/force_residual_history.png)：残余力如何降至门槛以下。
- [ADR versus reference](figures/adr_vs_reference.png)：同几何不同算法的位移与应变。
- [Pre-relaxation driver](figures/pre_relaxation_driver.png)：何时才允许启动曝光。
- [Erosion time step consistency](figures/erosion_time_step_consistency_equilibrated.png)：平衡后的慢步曲线。
- [Final driver and damage maps](figures/final_driver_and_damage_maps.png)：真实1×位移，各色条使用各自实际min/max。

driver图为邻接bond拉伸应变的influence权重平均；粒子D是邻接bond损伤的算术平均，
都沿用已有聚合定义，不把两者的色值混为一谈。颜色亮不代表D=1。

# 11. Energy warning 是否仍存在

**仍存在：ENERGY_DISSIPATION_NOT_ESTABLISHED。** 每次保存同一x下的U_before、U_after及
delta_U=U_before−U_after，保留负号；没有将虚拟阻尼耗散记成物理耗散，没有修改damage law。

| 曝光步长 | 储能增加次数/更新数 | 累计signed delta_U(J) |
|---|---|---|
{chr(10).join(energy_rows)}

在更新前已经通过机械平衡时，损伤跳变仍导致储能增加，因此这一现象不能仅归因于启动振荡。
它表明当前NOSB完整性更新的热力学耗散尚未建立；不因此否认已验证的力平衡，也不借平衡通过掩盖能量问题。
本轮不修复它，下一轮也没有自动启动。

# 12. 当前状态标签

```text
MECHANICAL_EQUILIBRIUM: {status['mechanical_equilibrium']}
ADR: VERIFIED FOR CURRENT DEVELOPMENT CASE, WITH DOCUMENTED AUTHOR GUARDS
SCIPY REFERENCE: PASSED FOR SMALL NO-DAMAGE CASE
EROSION UPDATES: ALL 35 GATED BY EQUILIBRIUM
SLOW-STEP COMPARISON: PASSED PREDECLARED DEVELOPMENT ACCEPTANCE
EROSION ENERGY: ENERGY_DISSIPATION_NOT_ESTABLISHED
FAILURE / DETACHMENT: NOT OBSERVED
PHYSICAL CALIBRATION: NOT PERFORMED
MESH OBJECTIVITY: NOT ASSESSED / STILL OPEN
WALL ATTACHMENT: STILL OPEN
REAL MICROBUBBLE CFD: NOT USED
```

本轮交付限于这一机械平衡问题。没有增加材料、壁面、疲劳、声流或断裂物理。
源码、配置、历史数据和图的数据来源均可追溯；保留检查见 `provenance/PRESERVATION_CHECK.json`。
'''
    write(V/'REVIEW_MECHANICAL_EQUILIBRIUM.md',text)


def main():
    log=(V/'FULL_TESTS.log').read_text();m=re.search(r'(\d+) passed',log)
    assert m and int(m.group(1))==112 and 'FAILED' not in log and 'ERROR ' not in log
    write(V/'TEST_SUMMARY.json',json.dumps(dict(old_tests_passed=78,previous_reference_tests_passed=15,new_adr_tests_passed=19,total_tests_passed=112,failures=0,skipped=0),indent=2))
    assert read(V/'STATE_AUDIT.json')['status']=='PASS'
    FIG.mkdir(exist_ok=False);style();plots();report()
    sources=list(R.glob('*/states.npz'))+list((V/'cases').glob('*/equilibrium/state.npz'))+list((OLD/'cases').glob('*/states.npz'))
    write(FIG/'IDENTITY.json',json.dumps(dict(background='black',font='Arial',font_weight='normal',arrows=False,
        displacement_scale=1.,data_smoothing=False,source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}),indent=2))
    print(V/'REVIEW_MECHANICAL_EQUILIBRIUM.md')


if __name__=='__main__':main()
