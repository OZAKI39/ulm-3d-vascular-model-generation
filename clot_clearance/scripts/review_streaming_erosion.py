"""Read-only analysis of erosion states; new PNG/PDF artifacts and review.

Uses actual saved damage and coordinates. No force arrows, bubble glyphs,
displacement amplification, interpolated simulation states or new physics.
"""
import csv,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager,colors
from matplotlib.lines import Line2D
from pd_clot.regularization.erosion import StreamingErosionDamage

ROOT=Path(__file__).resolve().parents[1]
V=ROOT/'verification/streaming_erosion_v1'
R=ROOT/'results/streaming_erosion_demo_v1'
OUT=V/'figures'
PORTAL=ROOT/'visualization/streaming_erosion_v1'
CASES=[('dt_erosion_0p02',.02),('coarse',.01),('dt_erosion_0p005',.005)]
PALETTE=['#ffb454','#60ccf0','#e28bb6']


def read(path):return json.loads(path.read_text())


def output(path,text):
    with path.open('x') as f:f.write(text)


def scientific(axis,which='y'):
    axis.ticklabel_format(axis=which,style='sci',scilimits=(0,0),useMathText=True)


def save(fig,name):
    for ext in ['png','pdf']:
        path=OUT/f'{name}.{ext}'
        if path.exists():raise FileExistsError(path)
        fig.savefig(path,dpi=170,facecolor='black',bbox_inches='tight',pad_inches=.15)
    plt.close(fig)


def configure():
    font=Path('/mnt/c/Windows/Fonts/arial.ttf')
    if font.exists():font_manager.fontManager.addfont(str(font))
    resolved=Path(font_manager.findfont(font_manager.FontProperties(family='Arial'),fallback_to_default=False))
    plt.rcParams.update({'font.family':'Arial','font.weight':'normal','font.size':13,
        'axes.titlesize':19,'axes.labelsize':14,'axes.labelweight':'normal','axes.titleweight':'normal',
        'figure.facecolor':'black','axes.facecolor':'black','savefig.facecolor':'black',
        'text.color':'white','axes.labelcolor':'white','xtick.color':'white','ytick.color':'white',
        'axes.edgecolor':'#687b90','grid.color':'#24364c','grid.alpha':.8,'axes.grid':True,
        'legend.frameon':False,'legend.fontsize':12,'pdf.fonttype':42,'ps.fonttype':42})
    return dict(family='Arial',weight='normal',sha256=hashlib.sha256(resolved.read_bytes()).hexdigest())


def rate_law(c):
    eps=np.array([0.,.5e-4,1e-4,2e-4]);t=np.linspace(0,.1,101);steps=np.arange(6)*.02
    damage=StreamingErosionDamage(4,c['streaming_erosion']);sample=[damage.D.copy()]
    for _ in range(5):damage.advance(eps,.02);sample.append(damage.D.copy())
    sample=np.array(sample);expected=steps[:,None]*.01*eps[None,:]/1e-4
    fig,ax=plt.subplots(figsize=(11,6),layout='constrained')
    for k,col in enumerate(['#bccad8',*PALETTE]):
        ax.plot(t,t*.01*eps[k]/1e-4,color=col,lw=2,label=f'Tensile strain = {eps[k]:.1e}')
        ax.plot(steps,sample[:,k],'o',mfc='black',mec=col,ms=7,mew=1.5)
    ax.set(title='Erosion rate law',xlabel='Erosion exposure time (s)',ylabel='Bond damage D',xlim=(0,.1))
    scientific(ax);ax.legend(loc='upper left');ax.text(.98,.06,'Lines: analytic solution\nOpen circles: numerical updates',transform=ax.transAxes,ha='right',va='bottom',fontsize=12)
    save(fig,'erosion_rate_law')
    with (OUT/'rate_law_samples.csv').open('x',newline='') as f:
        w=csv.writer(f);w.writerow(['erosion_exposure_time_s','epsilon_plus','numerical_damage','analytic_damage'])
        for i,time in enumerate(steps):
            for k,e in enumerate(eps):w.writerow([time,e,sample[i,k],expected[i,k]])
    return float(np.max(abs(sample-expected)))


def step_plot(histories):
    fig,axes=plt.subplots(1,2,figsize=(14,6),layout='constrained')
    for (name,step),col in zip(CASES,PALETTE):
        h=histories[name];t=np.array([r['erosion_exposure_time_s'] for r in h]);D=np.array([r['mean_damage'] for r in h])
        axes[0].plot(t,D,color=col,marker='o',ms=4,lw=1.5,label=f'Exposure step = {step:g} s')
        if step!=.005:
            times=np.arange(1,6)*.02;coarse=np.array([min(h,key=lambda r:abs(r['erosion_exposure_time_s']-tt))['mean_damage'] for tt in times])
            ref=np.array([min(histories['dt_erosion_0p005'],key=lambda r:abs(r['erosion_exposure_time_s']-tt))['mean_damage'] for tt in times])
            axes[1].plot(times,100*(coarse-ref)/np.maximum(abs(ref),1e-12),'-o',color=col,lw=1.5)
    axes[0].set(title='Mean particle damage',ylabel='Mean particle damage',xlabel='Erosion exposure time (s)')
    axes[1].set(title='Difference from 0.005 s',ylabel='Relative mean damage difference (%)',xlabel='Erosion exposure time (s)')
    axes[1].axhline(0,color='#adb9c6',lw=.6);scientific(axes[0]);axes[0].legend(loc='upper left')
    fig.suptitle('Erosion time step comparison — inconclusive settling',fontsize=21,fontweight='normal')
    save(fig,'erosion_time_step_consistency')


def settling_plot(windows):
    fig,axes=plt.subplots(1,3,figsize=(17,5.5),layout='constrained');cmap=plt.colormaps['viridis'];norm=colors.Normalize(1,len(windows))
    for window in windows:
        samples=window['samples'];t=np.array([r['window_time_s'] for r in samples])*1000;col=cmap(norm(window['macro_step']))
        axes[0].semilogy(t,[r['kinetic_energy_J'] for r in samples],color=col,lw=1.5)
        axes[1].plot(t,[r['maximum_epsilon_plus'] for r in samples],color=col,lw=1.5)
    for ax in axes[:2]:ax.set_xlabel('Mechanical time within macro window (ms)')
    axes[0].set(title='Kinetic energy',ylabel='Kinetic energy (J)');axes[1].set(title='Maximum tensile strain',ylabel='Maximum tensile strain');scientific(axes[1])
    indices=[w['macro_step'] for w in windows];drift=[100*w['summary']['driver_relative_drift'] for w in windows]
    axes[2].bar(indices,drift,color=['#f07e7e' if d>5 else '#60ccf0' for d in drift]);axes[2].axhline(5,color='white',ls='--',label='5% tolerance')
    axes[2].set(title='End-window driver drift',xlabel='Macro step index',ylabel='Relative driver drift (%)');axes[2].legend();axes[2].set_xticks([1,5,10])
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=axes[:2],orientation='horizontal',fraction=.06,pad=.16,shrink=.6)
    cb.set_label('Macro step index');cb.set_ticks([1,5,10])
    fig.suptitle('Mechanical settling — final 20% of each 2 ms window',fontsize=21)
    save(fig,'mechanical_settling')


def spatial_plot(z,c,values,name,title,label):
    x=z['x'][-1]*1000;fixed=z['fixed'];source=np.array(c['streaming']['bubble_center_m'])*1000
    assert source[0]<c['clot']['origin_m'][0]*1000
    vmin,vmax=float(values.min()),float(values.max());norm=colors.Normalize(vmin,vmax if vmax>vmin else vmin+1e-12)
    fig=plt.figure(figsize=(12,9.4),facecolor='black');ax=fig.add_axes([.015,.21,.92,.72],projection='3d',facecolor='black')
    ax.set_proj_type('ortho');ax.view_init(elev=24,azim=-68)
    ax.set(xlim=(-.97,.73),ylim=(-.49,.49),zlim=(-1.18,-.21),xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)')
    ax.set_box_aspect([1.7,.98,.97],zoom=1.13)
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
        axis.label.set_color('white');axis.line.set_color('#42556c');axis.set_pane_color((0,0,0,1));axis._axinfo['grid']['color']='#24364c'
    ax.tick_params(colors='white',labelsize=12,pad=1)
    for mask,marker in [(~fixed,'o'),(fixed,'s')]:
        ax.scatter(*x[mask].T,c=values[mask],norm=norm,cmap='inferno',s=44,marker=marker,
                   edgecolors='#d9e0e9' if marker=='s' else '#788697',linewidths=.55 if marker=='s' else .3,depthshade=False)
    ax.scatter(*source,marker='D',s=80,facecolors='none',edgecolors='white',linewidths=1.4,depthshade=False)
    fig.text(.045,.955,title,fontsize=26,weight='normal')
    fig.text(.045,.908,f'Actual min = {vmin:.5e}   |   Actual max = {vmax:.5e}',fontsize=14)
    handles=[Line2D([],[],marker='D',mfc='none',mec='white',ls='none',ms=9,label='Prescribed streaming center'),
             Line2D([],[],marker='s',mfc='none',mec='white',ls='none',ms=9,label='Fixed region'),
             Line2D([],[],color='none',label='Displacement scale: 1×')]
    fig.legend(handles=handles,loc='center',bbox_to_anchor=(.5,.153),ncol=3,fontsize=12,columnspacing=1.3)
    cax=fig.add_axes([.24,.073,.52,.022]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='inferno'),cax=cax,orientation='horizontal',format='%.2e')
    cb.set_label(label,fontsize=14);cb.ax.tick_params(labelsize=11)
    save(fig,name);return dict(minimum=vmin,maximum=vmax)


def report(c,z,histories,study,test_summary,figure_meta):
    h=histories['coarse'];end=h[-1];summary=read(R/'coarse/SUMMARY.json');audit=read(V/'comparisons/STATE_AUDIT.json')
    maxrel=study['maximum_relative_difference'];drift=summary['maximum_driver_relative_drift']
    windows=read(R/'coarse/settling_windows.json');failed=[w['macro_step'] for w in windows if w['summary']['settling_status']!='DRIVER_SETTLING_CHECK_PASSED']
    table=[]
    for name,step in CASES:
        row=histories[name][-1];s=read(R/name/'SUMMARY.json')
        table.append(f"| {step:g} | {row['erosion_exposure_time_s']:.3f} | {row['mean_damage']:.8e} | {row['maximum_damage']:.8e} | {row['broken_bond_count']} | {row['detached_total_volume_m3']:.3e} | {row['maximum_displacement_m']:.7e} | {row['maximum_epsilon_plus']:.7e} | 未通过（首步） |")
    energy=[]
    for name,step in CASES:
        row=histories[name][-1]
        energy.append(f"| {step:g} | {row['mechanical_time_s']:.3f} | {row['damage_stored_energy_drop_J']:.8e} | {row['elastic_energy_J']+row['stabilization_energy_J']:.8e} | {row['fraction_of_damage_inside_influence_region']:.6f} | {row['damage_weighted_distance_from_source_m']*1000:.6f} |")
    comparisons='\n'.join(f"| `{key}` | {100*value['maximum_common_time_relative_difference']:.6f}% | {100*value['maximum_final_relative_difference']:.6f}% |" for key,value in study['by_metric'].items())
    code_files={
        'configs/streaming_erosion_demo_coarse.json':'固定几何、材料、mean-only 流场和预先声明的侵蚀参数。',
        'pd_clot/regularization/erosion.py':'实现单独的拉伸应变率损伤，并将末段漂移检查与损伤驱动分开。',
        'pd_clot/regularization/erosion_runner.py':'复用原力学和流场，在独立机械时钟与暴露时钟下运行并保存完整状态。',
        'pd_clot/regularization/erosion_energy.py':'保留损伤前后储能变化的符号，明确报告储能增加而不冒称正耗散。',
        'tests/test_streaming_erosion_damage.py':'验证解析解、不可逆性、阈值、恒定驱动步长一致性和带符号能量记录。',
        'tests/test_streaming_erosion_runner.py':'验证 mean-only 与静态导入流场、双时钟、旧周期参数独立性和粒子安全约束。',
        'scripts/verify_streaming_erosion.py':'执行三档固定参数计算，并独立核对保存的状态和共同暴露时刻差异。',
        'scripts/review_streaming_erosion.py':'从保存的状态生成五组 PNG/PDF、浏览页和本报告。',
        'results/streaming_erosion_demo_v1/.gitignore':'保留本轮状态与 VTP 文件，不修改旧忽略规则。',
        'verification/streaming_erosion_v1/.gitignore':'保留本轮验证输入和失败诊断状态，不修改旧忽略规则。'}
    files='\n'.join(f'- `{p}`：{description}' for p,description in code_files.items())
    damage=figure_meta['final_damage_map'];driver=figure_meta['spatial_erosion_driver']
    text=f'''# 1. 这一轮做了什么

本轮在原有工程上新增了独立的 STREAMING_EROSION 模式，未修改旧 cyclic 路径。
损伤由 mean flow 下 bond 的末端拉伸应变推动，按独立暴露时间积分。
演示参数事先固定，没有为得到断裂而调整载荷、损伤系数或阈值。
三个算例均完成 0.1 s 暴露，未出现断键或脱离。
首个 0.002 s 机械窗口没有通过预设的准稳态检查，所以这些损伤轨迹保留了尚未充分稳定的末端响应。
本轮同时保留了小算例和三档演示中损伤导致储能略增的证据，未将它隐藏为正耗散。
这是一轮软件与数值开发，不是血栓侵蚀的物理验证。

# 2. 模型公式

```text
stretch_ij = |x_j - x_i| / |X_j - X_i|
epsilon_plus_ij = max(stretch_ij - 1, 0)
rate_ij = C_t_s_inv × max(epsilon_plus_ij / epsilon_ref - threshold_ratio, 0)^m_t
D_new = min(1, D_old + erosion_exposure_step_s × rate_ij)
g = 1 - D; active iff D < 1
```

`X` 是参考位置，`x` 是实际当前位置；`stretch`、`epsilon_plus` 和损伤 `D` 都无量纲。
`C_t_s_inv = 0.01 s⁻¹` 是开发用速率，`epsilon_ref = 10⁻⁴` 是参考拉伸应变，`threshold_ratio = 0` 是无量纲阈值，`m_t = 1` 是指数。
`erosion_exposure_step_s` 是慢时间步长，基线为 0.01 s；`D_break = 1`，只有 D 达到 1 才断键。
所有侵蚀参数都是 **DEVELOPMENT / UNCALIBRATED PLACEHOLDERS**，Gc coupon calibration 不等于 streaming erosion calibration。

`mechanical_time_s` 只累计真正积分的机械窗口；`erosion_exposure_time_s` 只累计上式的慢暴露步长。
每个宏步机械窗口固定 0.002 s，机械 dt = 10⁻⁶ s，最后 20% 记录动能、最大和平均拉伸应变；以末端构型的应变更新 D，**不以时间窗 max-min 作为损伤驱动**。
漂移取“各 active bond 在检查窗内最大变化 / 窗内全局峰值应变”以及标量最大、平均应变的相对变化三者中的最大值，分母下限为 10⁻¹²，要求不超过 5%。
这只是预设的驱动稳定性检查；动能全部保存，检查通过也不等于真实物理准静态得到证明。
若未通过，则记录 `MECHANICAL_STEADY_STATE_NOT_ESTABLISHED`，仍以明确标为未验证准稳态的末端应变完成诊断积分；不延长窗口，不重调参数。

只借鉴了形变影响时间依赖退化的思路；没有移植论文的细胞模型、参数、空间腐蚀或生化过程。[Bajd & Serša (2013) 原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3870808/)

# 3. 代码改了哪些文件

以下均为新增文件；旧文件修改数为 0，最终核对见 `provenance/PRESERVATION_CHECK.json`。

{files}

新计算写在 `results/streaming_erosion_demo_v1/`；`coarse` 即 0.01 s 基线，`dt_erosion_0p01` 是指向它的相对链接，不是重复计算。
新验证、日志、比较、PNG/PDF、输入 SHA256 和本报告写在 `verification/streaming_erosion_v1/`；浏览页写在 `visualization/streaming_erosion_v1/`。
每档计算保存 CONFIG、IDENTITY、SUMMARY、history JSON/CSV、states NPZ、逐步 VTP、慢时间 PVD、source_snapshot 和 SHA256；图中每个点都来自这些保存数据。

# 4. 测试结果

- Old tests：48 / 48 passed；开始开发时单独回归通过，最终全量测试仍全部通过。
- New tests：{test_summary['new_tests_passed']} / {test_summary['new_tests_passed']} passed。
- Total：**{test_summary['total_tests_passed']} passed**，最终 failures = 0；命令为 `python -m pytest tests -q -p no:cacheprovider --basetemp build_streaming_erosion_v1/full_tests`。
- 覆盖零驱动、常驱动解析解、压缩不损伤、不可逆上限、D_break 强制为 1、三档慢步长解析一致性、DeltaN 与频率元数据独立性、mean-only 小算例、静态 resolved field 的末端应变驱动和双时钟。
- 初次新测试有 23 passed、5 fixture errors：同一小算例在固定位置损伤后储能增加 1.21846804×10⁻¹⁸ J，旧正耗散账本拒绝继续；原日志和失败状态完整保存在 `tests/NEW_TESTS_FIRST.log` 与 `tests/first_attempt_energy_rejection/`。
- 修复仅增加侵蚀专用的**带符号储能账本**及两项能量测试；未改材料、侵蚀参数、流场、旧账本和能量残差容限；能量增加仍会发出 `ENERGY_DISSIPATION_NOT_ESTABLISHED`，没有被截断或消除。
- 最终 50 条 warnings 中 47 条来自 VTK/NumPy 弃用提示，3 条是明确保留的准稳态或能量诊断；三档完整状态独立审计均 PASS。

测试日志：`tests/BASELINE_OLD_TESTS.log`、`tests/FULL_TESTS.log`；数值审计：`comparisons/STATE_AUDIT.json`。
恒定驱动解析图的最大绝对误差为 {figure_meta['constant_driver_maximum_absolute_error']:.3e}；这项结论独立于机械积分。

# 5. 数值结果

mean/max damage 指**粒子相邻原始 bond 损伤的算术平均**后的粒子均值/最大值，不等同于最大 bond D。
同一等体积 360 粒子、12592 bond、60 固定粒子的粗网格用于全部三档；基线最大 bond D = {end['maximum_bond_damage']:.8e}。

| 暴露步长 (s) | 总暴露 (s) | 平均粒子 D | 最大粒子 D | 断键数 | 脱离体积 (m³) | 最大位移 (m) | 最大拉伸应变 | Settling |
|---|---|---|---|---|---|---|---|---|
{chr(10).join(table)}

**NO FAILURE / NO DETACHMENT OBSERVED**。粒子数量与总体积保持，固定粒子位置误差精确为 0；没有人为切键或粒子删除。
基线最低 J = {summary['minimum_J_or_intrinsic_stretch']:.8f}，最大能量账本相对残差 = {audit[1]['max_energy_residual']:.6e}。
基线未通过的宏步为 {failed}，最差漂移 {100*drift:.6f}% > 5%；后续窗口通过不能撤销首步失稳对已累计损伤的影响。

| 暴露步长 (s) | 实际机械累计时间 (s) | 损伤处带符号累计储能下降 (J) | 最终储能 (J) | 半径 0.525 mm 内损伤占比 | 损伤加权离源距离 (mm) |
|---|---|---|---|---|---|
{chr(10).join(energy)}

储能下降定义为固定位置更新 D 前的储能减去更新后的储能；负值就是储能增加，不能称耗散。
该量只记录实际机械路径，没有乘以慢时间/机械时间比例，也没有采用 coupon Gc；逐 bond 分摊仅按 ΔD 比例作诊断，不是断裂能标定。
基线累计储能增加为 {end['damage_energy_increase_J']:.8e} J，超过原 roundoff 阈值的增加次数为 {end['damage_energy_increase_count']}，状态为 `{summary['damage_energy_status']}`。
本轮没有碎片，因此没有治疗时间中的碎片输运结论；原 transport 只随实际机械窗口执行。

# 6. Erosion time-step comparison

**{study['status']}**。
三档总慢暴露都是 0.1 s，机械窗口都为 0.002 s；宏步数不同使机械累计时间分别为 0.010、0.020、0.040 s，它们不能解释成治疗时间。
在共同慢时刻 0.02、0.04、0.06、0.08、0.10 s，以 0.005 s 档作为数值参考，计算 `abs(case-reference) / max(abs(reference), floor)`。
损伤、应变、比例分母下限为 10⁻¹²，位移/距离为 10⁻¹⁵ m，能量为 10⁻³⁰ J；初始零状态不参与相对差比较。

| 指标 | 共同慢时刻最大相对差 | 最终时刻最大相对差 |
|---|---|---|
{comparisons}

所有指标、共同慢时刻的最大相对差为 **{100*maxrel['relative_difference']:.6f}%**，出现在 `{maxrel['metric']}`、步长 {maxrel['step_s']:g} s、暴露 {maxrel['erosion_exposure_time_s']:.2f} s。
即使曲线接近，三档首个机械窗口均未通过预声明检查，所以按预声明规则不能据此断言慢时间积分已收敛，也不宣称网格客观性。
详表见 `comparisons/time_step_differences.csv` 与 `comparisons/TIME_STEP_STUDY.json`。

# 7. 图怎么看

1. [Erosion rate law](figures/erosion_rate_law.png)：线是恒定拉伸驱动的解析解，空心点是数值更新；正常情况是零驱动不损伤，其余点落在解析线上。当前最大误差 {figure_meta['constant_driver_maximum_absolute_error']:.3e}，支持率方程的软件实现。
2. [Erosion time step consistency](figures/erosion_time_step_consistency.png)：左图比较三档平均损伤，右图放大相对 0.005 s 档的差异，横轴都是慢暴露时间。当前共同慢时刻平均损伤最大相对差为 {100*study['by_metric']['mean_damage']['maximum_common_time_relative_difference']:.6f}%，但首步未稳定，因此结论仍为 INCONCLUSIVE。
3. [Mechanical settling](figures/mechanical_settling.png)：看每个 2 ms 窗口最后 20% 的动能与最大应变，以及末段漂移是否低于 5% 虚线；这里横轴明确是窗口内机械时间。当前第 1 步漂移 {100*drift:.4f}%，后续窗口通过，所以初始瞬态影响没有被误报为已排除。
4. [Spatial erosion driver](figures/spatial_erosion_driver.png)：颜色是原始影响权重加权的相邻 bond 末端拉伸应变，白色菱形标出 clot 左侧的预设流场中心，方块是固定区；合理的图应与实际加载及锚定约束一起判断，不能仅凭热点判物理正确。当前范围 {driver['minimum']:.5e}–{driver['maximum']:.5e}，按真实位置和 1× 位移绘制，没有箭头或蓝色微泡。
5. [Final damage map](figures/final_damage_map.png)：看真实粒子损伤及色条数值，亮色只代表本图的最大值而非 D=1。当前真实 min/max = {damage['minimum']:.5e}/{damage['maximum']:.5e}，最大 bond D 也远低于 1，未发生断裂。

五图均同时提供同名 PDF；保持黑底、Arial 常规字重、英文标注，浏览页为 `visualization/streaming_erosion_v1/index.html`。
图中的“源”是预设局部声流场中心，管流牵引仍是分布载荷；它不是求解得到的微泡表面，也不是画在 clot 上方的合力中心。

# 8. 当前科学状态

```text
STREAMING_EROSION v1: SOFTWARE / NUMERICAL DEVELOPMENT
PHYSICAL CALIBRATION: NOT PERFORMED
MESH OBJECTIVITY: NOT ASSESSED IN THIS STEP
WALL-ATTACHMENT MESH CONSISTENCY: STILL OPEN
CYCLIC_FATIGUE: UNCHANGED
COMBINED EROSION + FATIGUE: NOT IMPLEMENTED
REAL MICROBUBBLE CFD: NOT USED IN THIS VERIFICATION
PARAMETERS: DEVELOPMENT / UNCALIBRATED PLACEHOLDERS
EROSION TIME-STEP STUDY: INCONCLUSIVE
MECHANICAL_STEADY_STATE_NOT_ESTABLISHED
NO FAILURE / NO DETACHMENT OBSERVED
FRAGMENT RESOLUTION: NOT_ASSESSED_NO_DETACHMENT
EROSION ENERGY: ENERGY_DISSIPATION_NOT_ESTABLISHED
```

旧科学状态保持：SOFTWARE VERIFIED；fracture-energy coupon 为 ENERGY REGULARIZED FOR SPECIFIED OPENING PROTOCOL ONLY；旧 main streaming case 为 NOT CONVERGED；旧 N=25,000 算例为 NO FAILURE / NO DETACHMENT；旧 forced-failure case 为 FRAGMENT_RESOLUTION_INADEQUATE。
本轮未运行真实微泡 CFD、材料标定、空间腐蚀反馈、组合疲劳或生化过程；能量账本闭合、平滑图像和测试通过都不能证明局部物理精度。

# 9. 下一步建议

只做一个独立的**机械松弛窗口敏感性检查**：冻结当前 mean flow、材料和侵蚀参数，比较预先指定的较长窗口与本轮 0.002 s，检验首步末端应变能否稳定；本轮结果保持不变。
'''
    output(V/'REVIEW_STREAMING_EROSION.md',text)
    return list(code_files)


def main():
    OUT.mkdir(exist_ok=False);PORTAL.mkdir(parents=True,exist_ok=False);font=configure()
    c=read(R/'coarse/CONFIG.json');histories={name:read(R/name/'history.json') for name,_ in CASES}
    study=read(V/'comparisons/TIME_STEP_STUDY.json');tests=read(V/'tests/TEST_SUMMARY.json')
    with np.load(R/'coarse/states.npz') as data:z={k:data[k] for k in data.files}
    meta=dict(font=font,source_states_sha256=hashlib.sha256((R/'coarse/states.npz').read_bytes()).hexdigest(),
              background='black',force_arrows=False,bubble_glyph=False,displacement_scale=1.,
              source_position_m=c['streaming']['bubble_center_m'],damages_are_saved_values=True)
    meta['constant_driver_maximum_absolute_error']=rate_law(c);step_plot(histories);settling_plot(read(R/'coarse/settling_windows.json'))
    meta['spatial_erosion_driver']=spatial_plot(z,c,z['particle_erosion_driver'][-1],'spatial_erosion_driver','Spatial erosion driver','Weighted bond tensile strain')
    meta['final_damage_map']=spatial_plot(z,c,z['damage'][-1],'final_damage_map','Particle damage','Particle damage D')
    code=report(c,z,histories,study,tests,meta)
    output(OUT/'FIGURE_IDENTITY.json',json.dumps(meta,indent=2))
    titles=[('erosion_rate_law','恒定驱动解析核对'),('erosion_time_step_consistency','慢暴露步长比较：INCONCLUSIVE'),
            ('mechanical_settling','机械末段漂移：首步未通过'),('spatial_erosion_driver','末端拉伸应变与源位置'),('final_damage_map','真实粒子损伤：未断键、未脱离')]
    cards='\n'.join(f'<section><h2>{title}</h2><a href="../../verification/streaming_erosion_v1/figures/{name}.pdf">PDF</a><img src="../../verification/streaming_erosion_v1/figures/{name}.png" alt="{title}"></section>' for name,title in titles)
    output(PORTAL/'index.html',f'''<!doctype html><html lang="zh"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Streaming erosion v1</title><style>body{{margin:24px auto;padding:0 20px;max-width:1350px;background:#000;color:white;font:17px Arial,sans-serif}}h1,h2{{font-weight:normal}}h2{{font-size:22px}}a{{color:#72cff5}}img{{display:block;width:100%;height:auto}}section{{margin:38px 0}}p{{line-height:1.6}}</style><h1>Streaming erosion v1</h1><p>SOFTWARE / NUMERICAL DEVELOPMENT · 未标定参数<br>总慢暴露 0.1 s；首步 MECHANICAL_STEADY_STATE_NOT_ESTABLISHED；步长比较 INCONCLUSIVE。</p><p><a href="../../verification/streaming_erosion_v1/REVIEW_STREAMING_EROSION.md">审核报告</a> · <a href="../../verification/streaming_erosion_v1/comparisons/TIME_STEP_STUDY.json">完整比较数据</a></p>{cards}</html>''')
    output(V/'provenance/ADDED_CODE_FILES.json',json.dumps(code,indent=2))
    print(json.dumps(dict(figures=str(OUT),report=str(V/'REVIEW_STREAMING_EROSION.md'),figure_identity=meta),indent=2))


if __name__=='__main__':main()
