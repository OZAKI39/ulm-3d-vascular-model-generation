"""Human review of the completed, paper-independent verification subset.

Missing ADR, pre-relaxation and erosion results are explicitly reported as
NOT RUN. No legacy result is relabeled as an equilibrated erosion result.
"""
import json,re,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT=Path(__file__).resolve().parents[1]
V=ROOT/'verification/mechanical_equilibrium_v1'
FIG=V/'figures'
PALETTE=['#ffb454','#60ccf0','#e28bb6','#b8e986']


def read(p):return json.loads(p.read_text())


def textfile(path,text):
    with path.open('x') as f:f.write(text)


def summary(name):return read(V/'cases'/name/'SUMMARY.json')


def history(name):
    with np.load(V/'cases'/name/'history.npz') as z:return {k:z[k] for k in z.files}


def style():
    font=Path('/mnt/c/Windows/Fonts/arial.ttf');font_manager.fontManager.addfont(str(font))
    plt.rcParams.update({'font.family':'Arial','font.weight':'normal','font.size':12,'axes.labelsize':13,
        'axes.titlesize':18,'axes.titleweight':'normal','axes.labelweight':'normal','figure.facecolor':'black',
        'axes.facecolor':'black','text.color':'white','axes.labelcolor':'white','xtick.color':'white','ytick.color':'white',
        'axes.edgecolor':'#687b90','grid.color':'#24364c','axes.grid':True,'legend.frameon':False,
        'pdf.fonttype':42,'ps.fonttype':42})


def save(fig,name):
    for ext in ['png','pdf']:fig.savefig(FIG/f'{name}.{ext}',dpi=170,facecolor='black',bbox_inches='tight',pad_inches=.16)
    plt.close(fig)


def scientific(ax):ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useMathText=True)


def plots():
    ref=summary('small_scipy')['final'];fig,axs=plt.subplots(1,2,figsize=(14,5.8),layout='constrained')
    for ax,size,label in zip(axs,['small','coarse'],['Small clot: 36 particles','Baseline clot: 360 particles']):
        for ms,col in zip([2,4,8,16],PALETTE):
            h=history(f'{size}_legacy_{ms}ms');ax.plot(h['mechanical_time_s']*1000,h['maximum_epsilon_plus'],color=col,lw=1.1,label=f'{ms} ms window')
            ax.plot(ms,h['maximum_epsilon_plus'][-1],'o',color=col,ms=7)
        if size=='small':ax.axhline(ref['maximum_epsilon_plus'],color='white',ls='--',lw=1.2,label='SciPy equilibrium')
        ax.set(title=label,xlabel='Mechanical integration time (ms)',ylabel='Maximum tensile strain');scientific(ax)
    axs[0].legend(loc='lower right',fontsize=10);axs[1].text(.98,.04,'ADR not run: primary paper needed\nNo coarse SciPy solve',transform=axs[1].transAxes,ha='right',va='bottom',fontsize=11)
    fig.suptitle('Fixed window convergence — no damage',fontsize=21)
    save(fig,'fixed_window_convergence')

    fig,ax=plt.subplots(figsize=(11,5.8),layout='constrained')
    for size,col in [('small','#60ccf0'),('coarse','#ffb454')]:
        h=history(f'{size}_legacy_16ms');ax.semilogy(h['mechanical_time_s']*1000,h['relative_force_residual'],color=col,lw=1.4,label=f'{size.capitalize()} clot, fixed damping')
    ax.axhline(1e-4,color='white',ls='--',lw=1.2,label='Force tolerance: 1e-4')
    ax.set(title='Force residual history — legacy method only',xlabel='Mechanical integration time (ms)',ylabel='Relative free-particle force residual',ylim=(1e-13,10))
    ax.legend(loc='upper right');save(fig,'force_residual_history')

    labels=['SciPy','2 ms','4 ms','8 ms','16 ms'];rows=[ref]+[summary(f'small_legacy_{ms}ms')['final'] for ms in [2,4,8,16]]
    fig,axs=plt.subplots(1,2,figsize=(13,5.6),layout='constrained')
    for ax,key,factor,title,ylabel in [(axs[0],'maximum_displacement_m',1e9,'Maximum displacement','Maximum displacement (nm)'),
                                     (axs[1],'maximum_epsilon_plus',1,'Maximum tensile strain','Maximum tensile strain')]:
        ax.bar(labels,[r[key]*factor for r in rows],color=['#dddddd',*PALETTE]);ax.set(title=title,ylabel=ylabel)
    scientific(axs[1]);fig.suptitle('Small clot: SciPy versus fixed windows — ADR unavailable',fontsize=20)
    save(fig,'reference_vs_fixed_windows')

    h=history('coarse_legacy_16ms');t=h['mechanical_time_s']*1000;first=summary('coarse_legacy_16ms')['first_converged_iteration']*.001
    fig,axs=plt.subplots(1,3,figsize=(17,5.5),layout='constrained')
    axs[0].semilogy(t[1:],h['physical_kinetic_energy_J'][1:],color=PALETTE[0]);axs[0].set(title='Physical-mass kinetic energy',ylabel='Kinetic energy (J)')
    axs[1].semilogy(t,h['relative_force_residual'],color=PALETTE[1]);axs[1].axhline(1e-4,color='white',ls='--');axs[1].set(title='Free-particle force residual',ylabel='Relative force residual')
    axs[2].plot(t,h['maximum_epsilon_plus'],color=PALETTE[2]);axs[2].set(title='Maximum tensile strain',ylabel='Maximum tensile strain');scientific(axs[2])
    for ax in axs:
        ax.axvline(first,color='#b8e986',ls=':',lw=1.4);ax.set_xlabel('Mechanical integration time (ms)')
    fig.suptitle(f'No-damage startup — first legacy criterion pass at {first:.2f} ms; erosion clock stays at zero',fontsize=19)
    save(fig,'no_damage_startup_diagnostics')


def report():
    compare=read(V/'COMPARISONS.json');hardware=read(V/'HARDWARE.json');test=read(V/'tests/TEST_SUMMARY.json')
    ref=summary('small_scipy');resolved=summary('small_scipy_resolved');coarse=summary('coarse_legacy_16ms')
    rows=[]
    for size in ['small','coarse']:
        for ms in [2,4,8,16]:
            s=summary(f'{size}_legacy_{ms}ms');f=s['final']
            rows.append(f"| {36 if size=='small' else 360} | {ms} | {f['maximum_epsilon_plus']:.10e} | {f['relative_force_residual']:.7e} | {f['physical_kinetic_energy_J']:.5e} | {'PASS' if s['converged'] else 'NOT PASSED'} |")
    comparisons=[]
    for r in compare['small_fixed_window_vs_scipy']:
        comparisons.append(f"| {r['window_ms']} | {r['rms_displacement_difference_m']:.5e} | {r['maximum_displacement_difference_m']:.5e} | {100*r['rms_difference_over_reference_max']:.6e}% | {100*r['relative_max_driver_difference']:.6e}% | {'PASS' if r['declared_comparison_passed'] else 'NOT PASSED'} |")
    energies=[]
    for name,label in [('small_scipy','SciPy')]+[(f'small_legacy_{ms}ms',f'{ms} ms') for ms in [2,4,8,16]]:
        r=summary(name)['final'];energies.append(f"| {label} | {r['maximum_displacement_m']:.7e} | {r['mean_epsilon_plus']:.7e} | {r['stored_energy_J']:.7e} | {r['surface_traction_mean_Pa']:.7f} | {r['surface_traction_max_Pa']:.7f} |")
    old=read(ROOT/'results/streaming_erosion_demo_v1/coarse/SUMMARY.json')
    text=f'''# 1. 这一轮为什么做

上一轮从未变形 clot 突然开流，第一步结束时还在晃动，应变漂移约 8.48%，不能据此认为损伤驱动已经平衡。
本轮先问“剩余的力是否足够小、运动是否停下、拉伸应变是否稳定”，并保持材料、mean flow、机械 dt、阻尼和侵蚀参数不变。
**目前只完成了无需 ADR 论文的验证部分；整轮任务尚未完成，状态为 NEED_PAPER_INPUT。**

# 2. 我们采用什么已有方法

- **Kilic & Madenci ADR：未实现。** 指定论文 DOI 为 [10.1016/j.tafmec.2010.08.001](https://doi.org/10.1016/j.tafmec.2010.08.001)。出版社接口只返回元数据，作者机构页面只有摘要，OpenAlex 与 Semantic Scholar 未提供开放全文；本地还检查了 24 个非生成、去重 PDF 的前两页，没有找到该文。检索证据见 `resources/LOCAL_PDF_SEARCH.json`、`resources/openalex`、`resources/semantic_scholar`、`resources/elsevier`；没有根据二手文章或记忆填写虚拟质量、阻尼或更新公式。
- **Peridigm：设计参考，没有移植求解器。** 核查 BSD-3-Clause 许可证，阅读官方 `src/core/Peridigm.cpp` 的 QuasiStatic/NOX、自由度处理和力残差，以及 ElasticCorrespondenceQSTension 算例；固定源码提交为 `f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903`。其 quasi-static 路径求解静态力平衡，显式路径使用 velocity-Verlet；damped Newton 不能冒称 Kilic–Madenci ADR。[官方源码](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/core/Peridigm.cpp)、[2023 论文 §2.6](https://doi.org/10.1007/s42102-023-00100-0)。
- **SciPy：直接调用已有 least_squares。** 使用 trust-region reflective、三点数值 Jacobian 和 linear loss，只优化 36 粒子小算例中的 72 个自由坐标；固定粒子不进入未知量。位移除以网格间距、残差除以初始外力尺度只是数值缩放，不更改力平衡根。[官方文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html)

SciPy 的 xtol、ftol、gtol 均为 10⁻¹²，max_nfev=500，配置运行前已冻结。
独立性是“不同求解算法求同一力平衡方程”，不是第二个材料模型，也不是物理验证。
未复制整个 Peridigm 或 Trilinos；保存的许可证和示例用于审计，不作为本项目求解代码。

# 3. 什么叫“平衡”

自由粒子的 `R = F_internal + F_surface` 应接近零；两项在本项目中已是 N，不能再乘一次粒子体积。
记录 L2、RMS 和最大粒子残余力，RMS 定义为 L2 / sqrt(自由粒子数)。
相对残差为 `||R_free|| / max(||F_external_free||, 1e-12 N)`，预设门槛为 **10⁻⁴**。

对 fixed damping 同时要求物理质量计算的动能 / max(|储能|,10⁻²⁴ J) ≤ **10⁻⁶**，末段驱动漂移 ≤ **1%**。
每 50 个机械步检查一次，在这 50 步的应变范围中同时检查 bond 向量和最大/平均标量漂移；分母下限为 10⁻¹²。
需要连续 **5 次**通过全部条件，不只检查应变变化。
SciPy 是代数求根，不产生物理速度、机械积分时间或 ADR pseudo time；保存为 null 而不是假造动能史，并单独重算最终力残差。

配置见 `configs/mechanical_equilibrium_verification_v1.json`，预声明见 `provenance/PREDECLARED_PLAN.json`；没有跑完后放宽容限。
ADR 如后续实现，只能标记为 NUMERICAL EQUILIBRATION DEVICE，不能解释成血栓黏度，也不能把虚拟质量写入物理密度。

# 4. 代码改了什么

以下均为新增；旧源码、旧 78 项测试和历史结果未修改。

- `pd_clot/regularization/equilibrium.py`：提供统一平衡接口、冻结拓扑力回调、原固定阻尼法和小规模 SciPy 基准；ADR 调用明确抛出 NEED_PAPER_INPUT。
- `configs/mechanical_equilibrium_verification_v1.json`：冻结联合门槛、数值缩放和 2/4/8/16 ms 窗口设置。
- `tests/test_equilibrium_reference.py`：增加 15 项前置验证，核对参考根、固定点、静态输入、时钟及旧路径逐元素一致性。
- `scripts/verify_mechanical_equilibrium.py`：按照硬件预算并行运行 10 个无损伤算例，保存逐步历史、独立日志、状态、源码快照和哈希。
- `scripts/review_mechanical_equilibrium.py`：生成实际已完成部分的四组 PNG/PDF 和本报告。
- `verification/mechanical_equilibrium_v1/.gitignore`：保留本轮状态和静态输入文件。

未修改旧 `erosion_runner.py`；由于 ADR 依据和验证未完成，还没有创建或接入正式 equilibrated erosion runner。
硬件记录：{hardware['cpu_model']}，{hardware['logical_cpu_count']} 个逻辑 CPU、{hardware['available_physical_cores']} 个可用物理核心，RAM {hardware['ram_total_bytes']/1024**3:.2f} GiB，初始可用 {hardware['ram_available_bytes']/1024**3:.2f} GiB；GPU 仅记录信息，现有 kernel 仍用 CPU。
两个 sweep 都使用 3 个进程，每个 BLAS/OMP/MKL 线程数为 1；按 1 GiB/worker 和至少 20% 内存余量预估，各 job 有独立目录和日志，实际峰值 RSS 保存于 SUMMARY。

# 5. 测试

- Old tests：**78 passed**；在新增测试前单独回归，最终全量回归仍通过。
- New tests：**{test['new_tests_passed']} passed**。
- Total：**{test['total_tests_passed']} passed**；failures = {test['failures']}；未删除、改弱、skip 或 xfail 旧测试。
- 新测试实际覆盖：零载荷残余力、两方法固定坐标、SciPy 独立根判据、legacy 与上一轮保存 x/v 逐元素相同、两方法 D=0/暴露时间=0、未通过条件的显式状态、缺论文时拒绝 ADR、禁止 SciPy 自动生产使用、位置相关表面载荷、DeltaN/频率元数据独立性、static resolved field，以及物理密度/输入不变。
- 尚未编写或执行 ADR 收敛、ADR 对照、ADR 迭代数独立性，以及正式 pre-relaxation/erosion runner 的损伤推进门控测试；这些不是已通过项目，也没有用 skipped tests 充数。

日志：`tests/BASELINE_78_TESTS.log`、`tests/NEW_TESTS_FIRST.log`、`tests/FULL_TESTS.log`。
窗口未通过平衡条件是被保留的数值观察，不是删除或掩盖的 pytest 失败。

# 6. 小算例独立验证

36 粒子、12 固定粒子、72 自由坐标、550 bonds；几何尺寸沿用上一轮 small fixture，材料和流场与原设置一致。
SciPy manufactured reference：**{ref['status']}**；相对力残差 **{ref['final']['relative_force_residual']:.8e}**，L2 **{ref['final']['force_residual_L2_N']:.8e} N**。
Termination reason：`{ref['scipy_termination_reason']}`；SciPy nfev={ref['scipy_nfev']}，计入数值 Jacobian 的真实残差评估次数={ref['actual_residual_evaluations']}，最终 Jacobian 条件数={ref['final_jacobian_condition']:.6g}。
Static resolved-field reference 也通过，相对力残差 **{resolved['final']['relative_force_residual']:.8e}**；它只用于静态力平衡，不解释成 cyclic fatigue。

| fixed window (ms) | RMS 位移差 (m) | 最大位移差 (m) | RMS差/参考最大位移 | 最大应变相对差 | 联合比较 |
|---|---|---|---|---|---|
{chr(10).join(comparisons)}

RMS 和 L2 比较只统计自由粒子；门槛为相对残差≤10⁻⁴、最大应变相对差≤1%、RMS差/参考最大位移≤1%。
8/16 ms 在这个小问题上满足预声明对照条件，但 **ADR vs SciPy 尚未进行**；小算例不能直接替代 360 粒子的独立基准。

| 方法 | 最大位移 (m) | 平均拉伸应变 | 储能 (J) | 平均牵引 (Pa) | 最大牵引 (Pa) |
|---|---|---|---|---|---|
{chr(10).join(energies)}

完整比较含 relative L2 displacement difference，见 `COMPARISONS.json`。

# 7. Fixed-window sensitivity

每个算例都从同一未变形 X、v=0、D=0 开始，dt=1 µs、damping=3500 s⁻¹；只改变窗口长度，不接续前一档最终构型。
全部窗口记录每个机械步的残余力、动能、最大/平均应变和最大位移，见各 case 的 `history.npz` 与 `history.csv.gz`。

| 粒子数 | 窗口 (ms) | 最大拉伸应变 | 相对残余力 | 最终动能 (J) | 联合门槛 |
|---|---|---|---|---|---|
{chr(10).join(rows)}

360 粒子 16 ms 最终力残差 L2={coarse['final']['force_residual_L2_N']:.8e} N，RMS={coarse['final']['force_residual_RMS_N']:.8e} N，max={coarse['final']['force_residual_max_N']:.8e} N。
本轮没有改变窗口直到“好看”；2/4 ms 的未通过结果全部保留，8/16 ms 是事先声明的两档。

图的阅读方式：

- [Fixed window convergence](figures/fixed_window_convergence.png)：看延长等待后应变稳定到哪里；左图是小算例并有 SciPy 参考线，右图是原 360 粒子基线，不能拿左图参考值去判断右图。
- [Force residual history](figures/force_residual_history.png)：看残余力是否降到 10⁻⁴ 虚线以下；横轴是原阻尼积分的机械时间，**不是 ADR iteration**，图明确标为 legacy only。
- [Reference vs fixed windows](figures/reference_vs_fixed_windows.png)：只比较同一个小算例的位移和应变，各用一个纵轴；没有画未计算的 ADR 柱。
- [No-damage startup diagnostics](figures/no_damage_startup_diagnostics.png)：看基线 clot 的动能、残余力和应变如何稳定；竖线是 legacy 条件首次连续通过，不是侵蚀时钟启动线，整段暴露时间都为 0。

四图均有同名 PDF，黑底、Arial 常规字重、英文标注；没有修改科学数据来平滑曲线。

# 8. Pre-relaxation 结果

**正式 pre_erosion_equilibration：NOT RUN。**
无损伤基线 first_converged_iteration={coarse['first_converged_iteration']}，对应 {coarse['first_converged_iteration']*1e-3:.2f} ms；这是原 fixed damping 法的诊断，不是 ADR 的虚拟时间或正式 pre-relaxation 验收。
本轮所有已跑状态 D=0、erosion_exposure_time_s=0；没有启动侵蚀时钟。
启动条件仍为：完成 ADR/小算例参考验证，并且 pre-relaxation 的力残差、运动量和驱动稳定性全部通过；随后每个损伤更新前均需通过同样门槛。

# 9. Equilibrated erosion rerun

**NOT RUN。** 因 ADR 方法来源与验证仍缺失，未写入任何声称“已平衡”的侵蚀结果。
旧 `results/streaming_erosion_demo_v1/` 完整保留；没有将旧 mean/max D 改名为本轮结果，也没有继续调用旧“未稳定仍更新损伤”的路径。
正式新 runner 的停机与保存状态分支尚待实现，本轮不宣称该分支已验证。

# 10. 时间步比较

**本轮 equilibrated erosion slow-step study：NOT RUN / NOT ASSESSED。**
未比较新的 0.02/0.01/0.005 s 侵蚀曲线，因此最大相对差异为 N/A；上一轮的 1.612081% 不属于本轮结果。
ADR residual 图、EROSION CLOCK STARTS HERE 图、equilibrated slow-step 图和最终 driver/damage maps 尚未生成；不会用占位数值或 D=0 验证图冒充正式结果。

# 11. Energy warning 是否仍存在

**ENERGY_DISSIPATION_NOT_ESTABLISHED 保持。**
本轮只做 D=0 力学验证，没有损伤跳变，故没有新的 U_before_damage、U_after_damage 或 delta_U_damage 结果。
上一轮基线累计损伤处储能增加为 {old['final']['damage_energy_increase_J']:.8e} J；本轮没有修改损伤方程、NOSB 材料、integrity 更新或带符号能量审计。
机械稳定和力残差变小不能消除这一旧警告；区分启动瞬态与损伤跳变的能量贡献，需要后续真正完成 equilibrated erosion 后才可判断。

# 12. 当前状态标签

```text
TASK: PARTIALLY COMPLETED; NEED_PAPER_INPUT
ADR: NOT IMPLEMENTED — PRIMARY METHODS UNAVAILABLE
SCIPY REFERENCE: VERIFIED FOR SMALL NO-DAMAGE CASE ONLY
LEGACY WINDOW COMPARISON: COMPLETED; 2/4 ms NOT PASSED, 8/16 ms PASSED
MECHANICAL_EQUILIBRIUM: NOT VERIFIED
PRE_EROSION_EQUILIBRATION: NOT RUN
EQUILIBRATED EROSION / SLOW-STEP STUDY: NOT RUN
PHYSICAL CALIBRATION: NOT PERFORMED
MESH OBJECTIVITY: NOT ASSESSED / STILL OPEN
WALL ATTACHMENT: STILL OPEN
REAL MICROBUBBLE CFD: NOT USED
EROSION ENERGY: ENERGY_DISSIPATION_NOT_ESTABLISHED
CYCLIC_FATIGUE / DAMAGE LAW: UNCHANGED
```

需用户提供 Kilic & Madenci (2010), DOI 10.1016/j.tafmec.2010.08.001 的 PDF 本地路径或上传论文。
这是用户提示词第 1、19 节的明确前置条件；取得原文前暂停 ADR 及其依赖步骤，不凭记忆补公式，不进入新的 clot physics 开发。
'''
    textfile(V/'REVIEW_MECHANICAL_EQUILIBRIUM.md',text)


def main():
    log=(V/'tests/FULL_TESTS.log').read_text();match=re.search(r'(\d+) passed',log)
    if not match:raise RuntimeError('Full tests not complete')
    total=int(match.group(1));assert total==93 and 'FAILED' not in log and 'ERROR ' not in log
    textfile(V/'tests/TEST_SUMMARY.json',json.dumps(dict(old_tests_passed=78,new_tests_passed=15,total_tests_passed=total,failures=0,skipped=0,ADR_tests='Not implemented pending required paper; not counted as passed'),indent=2))
    FIG.mkdir(exist_ok=False);style();plots();report()
    identity=dict(source_cases={p.parent.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (V/'cases').glob('*/states.npz')},
        figure_scope='No-damage legacy/SciPy verification only; no ADR or equilibrated erosion figures',
        background='black',font='Arial',font_weight='normal',arrows=False,science_data_altered=False)
    textfile(FIG/'IDENTITY.json',json.dumps(identity,indent=2));print(str(V/'REVIEW_MECHANICAL_EQUILIBRIUM.md'))


if __name__=='__main__':main()
