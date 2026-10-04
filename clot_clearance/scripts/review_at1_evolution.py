"""Assemble the human review without promoting failed coupled acceptance gates."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import subprocess
import xml.etree.ElementTree as ET
import numpy as np
from scripts.verify_at1_evolution import ROOT, OUT, config, write_json

RESULTS = ROOT/'results/pfpd_at1_evolution_v1'


def read(path): return json.loads(Path(path).read_text())


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(2**20),b''): h.update(chunk)
    return h.hexdigest()


def preservation():
    pre = read(OUT/'provenance/PRE_RUN_MANIFEST.json'); repo = Path(pre['repo_root'])
    after = {}; missing = []; modified = []
    for name, record in pre['files'].items():
        path = repo/name
        if not path.is_file(): missing.append(name); continue
        sha = digest(path); after[name] = dict(sha256=sha,bytes=path.stat().st_size)
        if sha != record['sha256']: modified.append(name)
    symlink_changes = []
    for name,target in pre.get('symlinks',{}).items():
        p = repo/name
        if not p.is_symlink() or os.readlink(p) != target: symlink_changes.append(name)
    check = dict(passed=not(modified or missing or symlink_changes),
                 checked_existing_regular_files=len(pre['files']), modified_files=modified,
                 missing_files=missing,symlink_changes=symlink_changes,
                 production_modified_files=sum('/pd_clot/' in n for n in modified),
                 historical_modified_files=len(modified),
                 current_HEAD=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                 checked_utc=datetime.now(timezone.utc).isoformat())
    write_json(OUT/'provenance/PRESERVATION_CHECK.json',check)
    return pre,after,check


def main():
    c = config(); unit = read(OUT/'analytic/UNIT_GATES.json'); fixed = read(OUT/'fixed_x/SUMMARY.json')
    runs = read(OUT/'staggered/ALL_RUNS.json'); extra = read(OUT/'staggered/EXTRA_ITERATIONS.json')
    xml = ET.parse(OUT/'tests/FINAL.xml').getroot(); counts = {}
    for mode in ['old','new']:
        cases = [a for a in xml.iter('testcase') if ('test_pfpd_at1_' in a.attrib['classname']) == (mode=='new')]
        counts[mode] = dict(total=len(cases),passed=sum(len(a)==0 for a in cases),
                            failed=sum(a.find('failure') is not None or a.find('error') is not None for a in cases),
                            skipped=sum(a.find('skipped') is not None for a in cases))
    counts['failed_tests'] = [a.attrib for a in xml.iter('testcase') if a.find('failure') is not None]
    write_json(OUT/'tests/COUNTS.json',counts)
    diagnostics = {}
    for name in runs:
        rows = read(RESULTS/name/'STAGGERED_HISTORY.json'); load = read(RESULTS/name/'LOAD_PATH.json')
        last = rows[-1]['load_step']; latest = [r for r in rows if r['load_step']==last]
        dx = [r['maximum_delta_s'] for r in latest if r['iteration']>0]
        gains = [b/a for a,b in zip(dx[:-1],dx[1:]) if a>1e-9 and b<1e-4]
        diagnostics[name] = dict(accepted_load_steps=len(load),requested_load_steps=len(c['schedules']['tensile_cube' if name=='tensile_cube' else 'loading_unloading']),
            last_phase_audited_iteration=latest[-1],observed_increment_gain_mean=float(np.mean(gains)) if gains else None,
            gain_note='Observed finite-iteration ratios, not a spectral-radius proof or a physical localization diagnosis',
            final_energy_balance_residual_J=load[-1]['energy_balance_residual_J'],
            final_residual_over_internal_energy=load[-1]['energy_balance_residual_J']/load[-1]['total_J'])
    write_json(OUT/'staggered/FAILURE_DIAGNOSIS.json',diagnostics)
    write_json(OUT/'energy/ACCOUNTING.json',{name:read(RESULTS/name/'LOAD_PATH.json') for name in runs})
    pre,after,preserve = preservation()
    config_sha = digest(OUT/'provenance/PREDECLARED_CONFIG.json')
    config_unchanged = config_sha == (OUT/'provenance/PREDECLARED_CONFIG.sha256').read_text().strip()
    gates = dict(A_old_169=counts['old']['passed']==169 and counts['old']['failed']==0,
        B_new_analytic=all(not('analytic' in r['classname']) for r in counts['failed_tests']),
        C_independent_minimum=unit['status']=='PASSED',D_virgin_threshold=unit['virgin_threshold'],
        E_irreversibility_of_executed_updates=all(r['minimum_accepted_phase_increment']>=0 for r in runs.values()),
        F_fixed_x_total_energy=all(r['maximum_phase_energy_increase_J']<=c['tolerances']['phase_energy_absolute_J'] for r in runs.values()),
        G_kinematics=all(r['maximum_kinematic_change']==0 for r in runs.values()),
        H_coupled_unloading_no_healing=runs['loading_unloading']['status']=='CONVERGED',
        I_staggered_with_extra_iteration_check=read(OUT/'staggered/COUPLED_GATE.json')['passed'],
        J_high_damage_stops=runs['tensile_cube']['status']=='HIGH_DAMAGE_GATE_REACHED',
        K_production_unchanged=preserve['passed'],predeclared_config_unchanged=config_unchanged)
    write_json(OUT/'ACCEPTANCE_GATES.json',gates)
    statuses = {
        'AT1 PAPER / LITERATURE REVIEW':'INCOMPLETE',
        'AUTHOR PFPD UPDATE FOUND':'NO',
        'PFPD_AT1_EVOLUTION_V1':'VERIFIED_FOR_EARLY_DAMAGE_H_EQUALS_1_REGIME' if all(gates.values()) else 'NOT_VERIFIED',
        'ANALYTIC VS NUMERIC MINIMIZATION':'PASSED' if unit['status']=='PASSED' else 'FAILED',
        'DAMAGE INITIATION THRESHOLD':'VERIFIED' if gates['D_virgin_threshold'] else 'NOT_VERIFIED',
        'IRREVERSIBILITY':'VERIFIED' if gates['H_coupled_unloading_no_healing'] else 'NOT_VERIFIED',
        'FIXED-X TOTAL ENERGY MINIMIZATION':'VERIFIED' if gates['F_fixed_x_total_energy'] else 'NOT_VERIFIED',
        'EARLY-DAMAGE KINEMATICS':'UNCHANGED' if gates['G_kinematics'] else 'NOT_UNCHANGED',
        'STAGGERED MECHANICS-PHASE':'VERIFIED_FOR_TENSILE_DEVELOPMENT_CASE' if gates['I_staggered_with_extra_iteration_check'] else 'NOT_VERIFIED',
        'HIGH-DAMAGE h(s) COUPLING':'NOT_IMPLEMENTED','FULL PFPD FRACTURE':'NOT YET VERIFIED',
        'STREAMING_EROSION COUPLING':'NOT IMPLEMENTED','PRODUCTION MODEL':'UNCHANGED' if preserve['passed'] else 'CHANGED',
        'PHYSICAL CLOT CALIBRATION':'NOT PERFORMED','REAL MICROBUBBLE CFD':'NOT USED'}
    write_json(OUT/'STATUS.json',statuses)
    code = sorted([str(p.relative_to(ROOT)) for directory in ['pd_clot','scripts','tests']
                   for p in (ROOT/directory).rglob('*.py')
                   if str(p.relative_to(Path(pre['repo_root']))) not in pre['files']])
    write_json(OUT/'provenance/ADDED_CODE_FILES.json',code)
    author = read(OUT/'resources/AUTHOR_UPDATE_CHECK.json'); paper = read(OUT/'resources/PAPER_IDENTITY.json')
    resources = read(OUT/'resources/AUTHOR_RESOURCE_SNAPSHOT.json')['repositories']
    figures = [
        ('at1_energy_landscape','为什么门槛以下不损伤：能量最低点在 s=0；圆点为约束最小值。'),
        ('at1_analytic_vs_numeric','解析曲线与独立数值点重合；虚线标出 Y/Yc=1。'),
        ('loading_unloading_irreversibility','单键卸载后相场不回退；各类能量独立绘制。图中明确注明立方体卸载尚未到达。'),
        ('fixed_x_energy_audit','弹性能下降、裂纹能增加，总能量下降；K/F 的最大变化为零。'),
        ('staggered_convergence','绿色是单调路径最后接受步；红/橙线是失败路径的交替迭代，不能把曲线增长解读为收敛。'),
        ('tensile_force_displacement','只连已接受反力点。阈值竖线为解析 1% 应变；圆圈为粗载荷表中首个超过相场容差的损伤点。'),
        ('phase_map_before_high_damage_gate','显示最后接受状态的真实位移，放大率 1；颗粒颜色只是键相场的加权后处理，不是侵蚀 D。')]
    e = fixed['energies_J']; mono = runs['tensile_cube']; lur = runs['loading_unloading']; ext = runs['loading_unloading_extra_iterations']
    repo_rows = '\n'.join(f"| [{name}]({r['url'].removesuffix('.git')}) | `{r['commit']}` | MIT |" for name,r in resources.items())
    status_rows = '\n'.join(f'| {name} | `{value}` |' for name,value in statuses.items())
    figure_rows = '\n'.join(f'- [{name}.png](figures/{name}.png) / [PDF](figures/{name}.pdf)：{desc}' for name,desc in figures)
    code_rows = '\n'.join('- `'+name+'`' for name in code)
    report = f'''# PFPD_AT1_EVOLUTION_V1 人工审核报告

**结论：NOT_VERIFIED。** 单键与固定位置 AT1 验证通过，单调拉伸按预期停在高损伤门限；完整加载—卸载—再加载及额外交替迭代检查失败。没有修改 Gc、载荷表或门限来消除失败，也没有把失败测试改成跳过。

## 1. 这一轮做什么

本轮新增独立 `PFPD_AT1_CONSISTENT_CANDIDATE`，不是生产模型替换。真实项目为 `{ROOT}`；Git 根为 `{pre['repo_root']}`；HEAD `{pre['head']}`。Git HEAD 之外已有大量未提交成果也纳入保护。开始和结束均核对适用 AGENTS.md，未发现适用文件。

使用 3×3×3 小立方体、27 个物质点，底层与顶层为夹持点，中间层 9 点保留 27 个自由度。点间距 0.125 mm，δ=0.25125 mm，夹持中心距离 0.25 mm。SVK 的 μ=1000 Pa、K=666.6666667 Pa，ν=0，以获得预期的纯拉伸仿射基准。相场从零开始；原有侵蚀 D 和历史 s_proxy 均未接入。

演示参数 `Gc_demo={c['phase']['Gc_demo_J_per_m2']:.12g} J/m²` 由预设 1% 工程应变的 Ψ0 设计，**DEMONSTRATION FRACTURE ENERGY / NOT CLOT CALIBRATION**。它不是实验拟合。载荷、Gc 和容差于候选运行前保存在 [PREDECLARED_CONFIG.json](provenance/PREDECLARED_CONFIG.json)，SHA256 `{config_sha}`；结束检查保持一致：{config_unchanged}。

新增代码及测试：

{code_rows}

输出只进入本轮 verification/results 目录。保留既有 {preserve['checked_existing_regular_files']} 个文件的 SHA256；生产修改 {preserve['production_modified_files']}，历史修改 {preserve['historical_modified_files']}，丢失 {len(preserve['missing_files'])}。见 [PRESERVATION_CHECK.json](provenance/PRESERVATION_CHECK.json)。

## 2. 为什么不用 Eq.3.25

本地 [PFPD v1](https://arxiv.org/abs/2604.09215) 的 SHA256 是 `{paper['sha256']}`。原文 γ(s)=s、g(s)=(1−s)² 的能量声明，与印刷的 Eq.3.25 关系存在上一轮已记录的歧义；旧 `UNRESOLVED_IN_V1` 结论未改写。本轮用户明确选择固定 x、h=1 下的 AT1 约束最小化，因此不把此实现叫作 Eq.3.25 的精确复现。

本轮能量密度 e=(1−s)²Y+2Yc·s，导数 q=−2(1−s)Y+2Yc，二阶导数 2Y≥0。处女键 Y≤Yc 取 0，否则试值为 1−Yc/Y；在 [s_old,1] 上取约束最小值。Eq.3.25 没有进入 solver。

正式编码前重新检查 arXiv（仅 v1）、作者仓库 80 分支、16 标签、16 release、156 PR、18 个作者仓库索引和 97 个源码 ref。当前 main `{author['main_sha']}`；Zenodo 最新可读版本 v0.5.4，record 22030918。相场关键词命中的是 2025-12-18 的 EquivalentStress 草稿 `f02c376956b03ce1083e163a026f6dc4528cab60`：其应力驱动软化和重构退化不能解决本任务 PFPD 相场律。没有发现正式 PFPD 更新，不能据此声称证明所有不可访问资源均不存在。SimZoo 无法克隆、Zenodo通用搜索返回400，限制已记录。

外部资源实际使用范围：

| 官方仓库 | 固定提交 | License |
|---|---|---|
{repo_rows}

PhaseFieldX 仅运行无 FEM 依赖的 AT1 裂纹几何函数和二次退化函数，交由 SciPy 求局部最小值。其 cw=8/3 与 PD c0 不同，先用 Gc_F=Gc_PD·cw/c0、ell=δ 换算，使局部裂纹系数相同。NewFrac 只核对 AT1、上下界及交替迭代做法，没有移植 FEM solver。

[Miehe 等 (2010)](https://doi.org/10.1016/j.cma.2010.04.011) 只取得期刊摘要和公开章节概要，尝试的高校全文链接失效。没有据摘要补写核心公式，因此完整文献综述状态仍标为 **INCOMPLETE**。本轮已实现的核心能量由本地 PFPD 原图、用户选择、直接代数和官方开源 AT1 函数支撑。详见 [AT1_OPEN_SOURCE_REVIEW.md](resources/AT1_OPEN_SOURCE_REVIEW.md) 与 [AUTHOR_UPDATE_CHECK.json](resources/AUTHOR_UPDATE_CHECK.json)。

## 3. 什么是 s、Y、Yc、g、h、sc

| 变量 | 含义和本轮取值 |
|---|---|
| s_pf_bond | 键断裂相场，初始 0；不是侵蚀 D |
| Y | 未损伤键能量密度 Ψ0，仅在规定纯拉伸域作为 Ψ0_plus |
| Yc | Gc/(2c0δ)={unit['Yc_Pa']:.10g} Pa；c0 从已有 normalization 模块读取 |
| g | (1−s)²，降低弹性能和相应内力 |
| h | 始终等于 1，保持变形测量的权重 |
| sc | 0.95，高损伤停止门限；达到或超过时不接受更新 |
| D_pf_particle | 仅作图，把键 s 用 BAQP 能量积分权重聚合到端点 |

每次相场更新前检查键 Green 应变最小特征值 ≥−1e−9。这里只判断纯拉伸适用域，没有用特征值拆分能量，也没有引入新的张压分裂。

## 4. analytic AT1 是否和数值最小化一致

16 个 Y/Yc 与 4 个 s_old 共 64 组全部通过。SciPy `minimize_scalar(method='bounded')` 独立最小化实际标量能量，然后显式比较两个端点；没有向数值目标传入解析答案。最大 |Δs|={unit['max_s_error']:.10g}，小于事先固定的 5e−8。PhaseFieldX 纯函数经归一化换算后的局部 oracle 最大误差同为 {unit['max_external_s_error']:.10g}。

原始数值解、端点候选、最终选择及逐项 KKT 保存在 [CASES.json](minimization/CASES.json)。内部最优点导数为零、不可逆下界处导数非负、更新不越界且固定 Y 总能量不增，均通过对应容差。

## 5. threshold 是否存在

处女单键在 Y/Yc≤1 时严格得到 s=0，略大于 1 时产生正相场，阈值测试 **VERIFIED**。三维算例在解析阈值的浮点重构中出现约 4.50e−14 的相场，这低于 1e−10 的相场容差，原值完整保存，没有截断。图中区分解析阈值和粗加载表首次采到的可分辨损伤点。

## 6. irreversible loading/unloading 是否通过

单键完整路径通过：加载到 Y/Yc=3、卸载到零、再加载到4，相场不下降，再加载超过峰值后继续增大。所有实际接受的多键更新也都满足 s_new≥s_old。

**完整耦合立方体路径尚未通过。** 标准路径只接受前 {lur['accepted_steps']}/19 步，到 ε=0.014 失败；额外两次收敛迭代路径只接受 {ext['accepted_steps']}/19 步，到 ε=0.012 失败。都在预设卸载之前停止。因此不能声称已验证立方体的卸载、再加载反力或残余刚度，也不绘制虚构的滞回圈；总体验收中的 IRREVERSIBILITY 标为 NOT_VERIFIED，并单列单键已通过的范围。

## 7. fixed-x energy 是否通过

固定位置 ε=0.012 的多键试验：

| 能量 | 更新前 (J) | 更新后 (J) |
|---|---:|---:|
| 弹性 | {e['elastic_before']:.12e} | {e['elastic_after']:.12e} |
| 裂纹 | {e['crack_before']:.12e} | {e['crack_after']:.12e} |
| 候选总能量 | {e['total_before']:.12e} | {e['total_after']:.12e} |

总能量下降 {e['total_before']-e['total_after']:.12e} J。全部执行过的耦合相场子步，最大正的总能量变化为 {max(r['maximum_phase_energy_increase_J'] for r in runs.values()):.12e} J，低于预设 1e−24 J 绝对容差；没有修正原始能量值。该结果只证明固定 x 的相场子问题，不证明耦合整体稳定。

## 8. K/F 是否仍保持不变

在固定 x 更新前后，K、Fbar、Ftilde、两端 shape derivatives 均 bitwise identical：最大 ΔK=0，最大 ΔF=0。所有执行的耦合相场子步也为零。原 BAQP 源文件保持原 SHA，不把 s 放进重构权重。已有归一化仅做轻量回归：同一 c0 函数得到 0.25，已有裂纹能积分与本轮加权值差 {fixed['existing_crack_energy_difference_J']:.5e} J。上一轮 Griffith 全套结果保持原样，未重跑重计算。

## 9. staggered mechanics-phase 是否收敛

机械部分复用已有 ADR `force_callback`，没有复制 ADR 算法。边界适配器只把移动夹持点传给 ADR，BAQP 真正参考坐标不变；自由点从前一步位置开始求解。没有流体力。旧接口中的 pipe 半径仅是宽松坐标保护字段，没有创建流体模型；ADR 伪时间、伪质量不解释为物理动力学。

单调拉伸的 7 个接受步均满足机械残差≤1e−9、相场增量≤1e−10、交替总能量变化≤1e−10。其最大接受机械残差为 {mono['maximum_accepted_force_residual']:.6e}。

但规定加载—卸载路径及额外交替迭代对照失败：标准路径失败步的相场增量从约 2.56e−10 逐渐放大到 4.88e−4；每次固定相场 ADR 虽收敛，下一次相场更新又破坏机械平衡。中段增量放大比约 1.93（这里只是观测到的迭代比，未经谱稳定性证明）。最终最小 Green 应变分别为 −1.115206881e−9 和 −1.773933946e−9，越过预声明适用域门限，立即保存 FAILED.json 和状态。

这表明“单调路径在容差下接受”不足以证明交替解对迭代次数可靠。目前尚不能判定原因全部属于数值迭代或软化局部化；没有通过强制对称、改动能量、降低精度、截断相场或放宽纯拉伸门限来隐藏它。额外迭代只完成共同前缀的比较，不能把该局部比较标成完整通过。见 [FAILURE_DIAGNOSIS.json](staggered/FAILURE_DIAGNOSIS.json)。

## 10. 是否触及 sc gate

单调路径最后接受 ε=0.04，最大 s={mono['maximum_accepted_s']:.12f}。ε=0.05 时提议最大 s={mono['event']['proposed_maximum_s']:.12f}>0.95，返回 **HIGH_DAMAGE_GATE_REACHED** 与原因 **HIGH_DAMAGE_KINEMATIC_COUPLING_NOT_IMPLEMENTED**。旧接受状态、当前机械平衡位置、Y、旧/试/新提议及能量均保存在结果 load_007；提议未接受，未截断，也未越界继续计算。这一预期门限通过；它与加载路径的真实 FAILED 状态分开记录。

## 11. energy accounting

Ψ0 是未退化密度；gΨ0 是退化弹性能密度；2Yc·s 是裂纹能密度。全局积分使用既有无向键系数 2·wbar·Vi·Vj，单位为 J。裂纹能没有称为 stored elastic energy。机械加载输入可以增加候选内能；只有固定 x 的相场子步要求总能量下降。

夹持反力为顶层内力的负和，外功按已接受反力与夹持位移作梯形积分。没有人为闭合。单调路径最大 |内能−外功|={mono['maximum_energy_balance_residual_J']:.12e} J；最后残差相对内能为 {diagnostics['tensile_cube']['final_residual_over_internal_energy']:.6%}。该差值包含离散加载/积分和求解误差，本轮没有把它解释成已证明的物理耗散。失败步不并入接受曲线或外功累计。全表见 [ACCOUNTING.json](energy/ACCOUNTING.json)。

## 12. 测试结果

| 范围 | 通过 | 失败 | 跳过 | 总数 |
|---|---:|---:|---:|---:|
| 原有测试 | 169 | 0 | 0 | 169 |
| 新增 AT1 | 39 | 3 | 0 | 42 |
| 合计 | 208 | 3 | 0 | 211 |

正式耦合前 34 项新增单键/固定位置测试全部通过。最终 3 项失败均明确指向完整耦合路径：同时收敛、卸载再加载、额外交替迭代不改变接受状态。这些失败断言保持有效，没有删除、xfail 或跳过。`pytest` 最终退出码为 1。见 [FINAL.log](tests/FINAL.log)、[FINAL.xml](tests/FINAL.xml)、[COUNTS.json](tests/COUNTS.json)。

已有测试包含生产、旧 BAQP 与旧验证结果回归。新增测试还覆盖只读输入、非法值、阈值、KKT、独立最小化、力的能量导数、纯拉伸域、高损伤门限、机械失败前不更新相场、迭代上限失败和秩亏保存状态。

硬件记录见 [HARDWARE.json](provenance/HARDWARE.json)：12 物理核、24 逻辑核、总内存约 8.18 GB。独立 oracle 与三个载荷路径由 ProcessPoolExecutor 运行，每个 BLAS worker 1 线程；按可用内存扣除 20% 总 RAM 和父进程预算选择 worker，上限 3。实际预算分别保存在 provenance 的 SCALAR/COUPLED_WORKER_BUDGET.json。未改写 GPU 求解器。

## 13. 图怎么看

七组图均提供 PNG/PDF，黑底、Arial 常规字重、英文、无箭头，全部直接使用保存结果。三维图位移倍率为 1；没有为了外观改科学数组。

{figure_rows}

## 14. 当前科学状态

| 标签 | 状态 |
|---|---|
{status_rows}

[ACCEPTANCE_GATES.json](ACCEPTANCE_GATES.json) 中 H（完整耦合卸载）与 I（包括额外迭代检查的交替收敛）失败，因此不授予 VERIFIED_FOR_EARLY_DAMAGE_H_EQUALS_1_REGIME。保留上一轮 PFPD v1 公式歧义，以及生产模型既有验证边界。

## 15. 唯一下一步

人工审核后，只开展当前纯拉伸小立方体的交替迭代稳定性诊断：以已保存失败状态区分相场—位移软化的局部化与数值固定点不稳定，先确定可验证的求解条件。此次到此停止，不自动进入高损伤、完整断裂、侵蚀耦合或物理标定。
'''
    (OUT/'REVIEW_PFPD_AT1_EVOLUTION_V1.md').write_text(report)
    # Managed new assets only; external clone objects, caches and test temp dirs
    # are reproducible scratch and are explicitly excluded from this inventory.
    added = {}
    for directory,dirs,files in os.walk(Path(pre['repo_root']),followlinks=False):
        dirs[:] = [d for d in dirs if d!='.git' and not Path(directory,d).is_symlink()]
        if '/pfpd_at1_evolution_v1/' in directory:
            dirs[:] = [d for d in dirs if d not in ['external','runtime_cache'] and not d.endswith('_tmp')]
        for name in files:
            p = Path(directory,name); rel = str(p.relative_to(Path(pre['repo_root'])))
            if rel in pre['files'] or p.is_symlink() or name=='POST_RUN_MANIFEST.json': continue
            added[rel] = dict(sha256=digest(p),bytes=p.stat().st_size)
    write_json(OUT/'provenance/POST_RUN_MANIFEST.json',dict(repo_root=pre['repo_root'],
        head=preserve['current_HEAD'],files=after,added_managed_files=added,
        excluded_new_scratch=['resources/external','provenance/runtime_cache','tests/*_tmp'],
        manifest_itself_excluded=True))
    print(json.dumps(dict(status=statuses,tests=counts,preservation=preserve,
                          added_managed_files=len(added),new_code_files=len(code)),indent=2))


if __name__ == '__main__': main()
