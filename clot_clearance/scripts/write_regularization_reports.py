"""Build reports and a local viewer from completed, audited result files."""
import hashlib
import html
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RUNS=ROOT/'results/streaming_regularized_demo'
VIS=ROOT/'visualization/streaming_regularized_demo'
VER=ROOT/'verification/regularization'


def read(path):return json.loads(path.read_text())


def table(headers,rows):
    return '| '+' | '.join(headers)+' |\n|'+'|'.join(['---']*len(headers))+'|\n'+''.join('| '+' | '.join(map(str,row))+' |\n' for row in rows)


def main():
    data=read(VER/'comparisons/COMPARISONS.json');rows={r['run']:r for r in data['runs']}
    assert len(rows)==18 and all(read(RUNS/n/'AUDIT.json')['status']=='PASS' for n in rows)
    qc=read(VIS/'QC.json');assert qc['status']=='PASS'
    assert '48 passed' in (VER/'FINAL_TESTS_SCOPED.log').read_text()
    base=rows['coarse'];history=read(RUNS/'coarse/history.json');final=history[-1]
    protected=read(ROOT/'provenance/regularization_20260930/PREEXISTING_SHA256.json')
    changed=[p for p,h in protected.items() if not (ROOT/p).is_file() or hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h]
    assert not changed,changed
    (VER/'FINAL_PRESERVATION.json').write_text(json.dumps(dict(status='PASS',protected_files=len(protected),changed_files=changed),indent=2)+'\n')
    mesh_table=table(['Case','h (mm)','Particles','Initial bonds','ψcrit (Pa)','Gc measured (J/m²)','Gc error (%)','First failure / detach N'],[
        [n,f"{rows[n]['spacing_mm']:.8g}",rows[n]['particle_count'],rows[n]['initial_bonds'],f"{rows[n]['calibrated_critical_energy_density_Pa']:.9g}",f"{rows[n]['Gc_measured_J_m2']:.10g}",f"{100*rows[n]['calibration_relative_error']:.7g}",'Not observed / not observed'] for n in ['coarse','medium','fine']])
    volume_table=table(['Case','Attached (mm³)','Detached (mm³)','Resolved (mm³)','Debris (mm³)','Singleton (mm³)','Rank 0/1 detached (mm³)','Damage dissipation (J)'],[
        [n,f"{rows[n]['total_volume_m3']*rows[n]['attached_volume_fraction']*1e9:.7g}",0,0,0,0,0,f"{rows[n]['damage_dissipation_estimate_J']:.9g}"] for n in ['coarse','medium','fine']])
    sensitivity_table=table(['Case','ΔN','dt (µs)','Mean D','Δ mean D (%)','Δ dissipation (%)','Max displacement (µm)','Final signed energy residual (%)'],[
        [n,rows[n]['DeltaN'],f"{rows[n]['dt_s']*1e6:g}",f"{rows[n]['mean_damage']:.9g}",f"{100*(rows[n]['mean_damage']/base['mean_damage']-1):.6g}",f"{100*(rows[n]['damage_dissipation_estimate_J']/base['damage_dissipation_estimate_J']-1):.6g}",f"{rows[n]['maximum_displacement_m']*1e6:.8g}",f"{100*rows[n]['relative_numerical_energy_residual']:.7g}"] for n in ['coarse','cycle_500','cycle_250','dt_half','medium','fine']])
    locality_table=table(['Case','Damage inside R (%)','Material initially inside R (%)','Damage-weighted distance (mm)','Fixed volume (%)'],[
        [n,f"{100*rows[n]['fraction_of_damage_inside_influence_region']:.6g}",f"{100*read(RUNS/n/'history.json')[-1]['initial_volume_in_influence_region_fraction']:.6g}",f"{rows[n]['damage_weighted_distance_from_source_m']*1e3:.7g}",f"{100*rows[n]['fixed_volume_fraction']:.6g}"] for n in ['coarse','medium','fine']])
    energy_keys=['elastic_energy_J','stabilization_energy_J','kinetic_energy_J','external_work_J','damping_dissipation_J','damage_dissipation_estimate_J','transport_work_J','numerical_energy_residual_J']
    energy_table=table(['Final coarse quantity','J'],[[k,f'{final[k]:.12g}'] for k in energy_keys])
    sweep_table=table(['Run','Gc demo (J/m²)','m_E','C_E','U (m/s)','Final mean D','Max D','First failure / detach N'],[
        [n,r['Gc_demo_J_m2'],r['m_E'],r['C_E'],r['streaming_amplitude_m_s'],f"{r['mean_damage']:.7g}",f"{r['maximum_damage']:.7g}",'Not observed'] for n,r in rows.items() if n.startswith(('sweep_','surface_'))])
    figures=sorted((VIS/'figures').glob('*.png'))
    figure_links='\n'.join(f'- [{p.stem} PNG]({p}) · [PDF]({p.with_suffix(".pdf")})' for p in figures)
    source_files=sorted([*ROOT.glob('pd_clot/regularization/*.py'),*ROOT.glob('tests/test_regularization_*.py'),*ROOT.glob('configs/streaming_regularized*.json'),
                        *[ROOT/'scripts'/n for n in ['render_regularized_damage.py','check_regularized_damage.py','verify_regularized_interfaces.py','write_regularization_reports.py']]])
    index='# Regularization source and output index\n\nAll paths are absolute. Existing source/config/result files were preserved.\n\n## New source, tests and configuration\n\n'
    index+='\n'.join(f'- [{p.name}]({p})' for p in source_files)
    index+='\n\n## Visual outputs\n\n'+figure_links
    index+=f'\n- [MP4]({VIS}/regularized_damage_transport.mp4)\n- [GIF]({VIS}/regularized_damage_transport.gif)\n- [Viewer]({VIS}/index.html)\n- [Preview]({VIS}/preview.png)\n- [QC]({VIS}/QC.json)\n'
    index+='\n## Run outputs\n\n'+''.join(f'- [{n}]({RUNS/n}): CONFIG.json, CALIBRATION.json, IDENTITY.json, AUDIT.json, SUMMARY.json, history.csv/json, states.npz, components.json, lineage.json, failure_ledger.jsonl, fracture_energy_ledger.jsonl, transport_diagnostics.json, clearance_ledger.json, particles.pvd, bonds.pvd, vtk/*.vtp, source_snapshot/.\n' for n in rows)
    (VER/'SOURCE_OUTPUT_INDEX.md').write_text(index)
    commands='''```bash
cd /home/lzy/projects/clot_clearance
export PYTHONDONTWRITEBYTECODE=1
export NUMBA_CACHE_DIR="$PWD/build_regularization/numba_cache"
export MPLCONFIGDIR="$PWD/build_fragmentation/matplotlib"
export PYTHONPATH="$PWD"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
CL_PY=/home/lzy/projects/computation_examples/_shared/python_environment/bin/python

# Coupon: use a new output directory; units are m and J/m².
"$CL_PY" -B -m pd_clot.regularization.coupon --config configs/streaming_regularized_demo_coarse.json --spacing 0.000125 --Gc-demo 0.01 --output verification/regularization/replay_coupon
# Demo using that new calibration.
"$CL_PY" -B -m pd_clot.regularization.runner --config configs/streaming_regularized_demo_coarse.json --calibration verification/regularization/replay_coupon/CALIBRATION.json --output results/streaming_regularized_demo_replay/coarse
# Medium/fine coupons are independently recalibrated by the mesh command.
"$CL_PY" -B -m pd_clot.regularization.campaign --study mesh --output-root results/streaming_regularized_demo_replay
# Cycle tests retain the baseline material and coarse calibration.
"$CL_PY" -B -m pd_clot.regularization.campaign --study cycle --output-root results/streaming_regularized_demo_replay
"$CL_PY" -B -m pd_clot.regularization.campaign --study timestep --output-root results/streaming_regularized_demo_replay
"$CL_PY" -B -m pd_clot.regularization.campaign --study sweep --output-root results/streaming_regularized_demo_replay
"$CL_PY" -B -m pd_clot.regularization.campaign --study sweep_surface --output-root results/streaming_regularized_demo_replay
# Current tests only; preserved snapshot directories are not a test suite.
"$CL_PY" -B -m pytest tests -q -p no:cacheprovider
"$CL_PY" -B scripts/verify_regularized_interfaces.py --mode drag --output verification/regularization/drag_replay
"$CL_PY" -B scripts/verify_regularized_interfaces.py --mode import --output verification/regularization/import_replay
# Presentation replay requires a fresh output directory.
"$CL_PY" -B scripts/render_regularized_damage.py --run results/streaming_regularized_demo/coarse --output visualization/streaming_regularized_demo_replay
"$CL_PY" -B scripts/check_regularized_damage.py --output visualization/streaming_regularized_demo_replay
"$CL_PY" -B -m pd_clot.regularization.figures --output visualization/streaming_regularized_demo_replay/figures
```
'''
    metadata='''```json
{
  "length_unit": "m", "velocity_unit": "m/s", "pressure_unit": "Pa", "time_unit": "s",
  "normal_convention": "outward_solid",
  "arrays": {"velocity": "velocity", "pressure": "pressure", "velocity_gradient": "velocity_gradient"},
  "validation": {"maximum_relative_divergence": 0.05, "maximum_relative_gradient_inconsistency": 0.25}
}
```
'''
    zh=f'''# 能量正则化开发报告

本轮沿用现有 `/home/lzy/projects/clot_clearance`，已完成新模块、18 组候选/敏感性计算、验证、报告和无箭头动画。**没有获得达到首次断键/脱落的主算例，也没有证明网格收敛。** 未为了得到碎裂画面提高载荷或降低断键阈值。

## 状态与验收边界

| 状态 | 结论 |
|---|---|
| SOFTWARE VERIFIED | 48 项当前测试通过；18 组主研究的数值/图/身份审计通过；旧强制碎裂重放数组逐位一致。限于已测试软件行为，不代表物理验证。|
| ENERGY REGULARIZED | 仅指本轮定义的循环张开 coupon 实测能量标定和不可逆损伤层；不是任意加载路径下已证明的客观断裂能。|
| NOT CONVERGED | 中/细档平均损伤比粗档下降 56.37% / 67.44%，损伤耗散下降 65.19% / 76.08%；没有失效/脱落事件供比较。|
| NOT_ASSESSED_NO_DETACHMENT | 新模型没有脱落物，不能判定碎片分辨率充分，也不能将 0/0 粉碎指标写成零。|
| FRAGMENT_RESOLUTION_INADEQUATE | 只对重新诊断的旧强制碎裂算例成立：P_singleton=78.899%，P_lowrank=82.569%。其载荷不同，不能作为同条件模型优劣比较。|

提示词软件验收第 5 项“粗、中档至少到首次脱落”未达到；含首次断键的 dt 窗口、断裂局域性、有限尺寸碎片占比与输运收敛亦未获得证据。按提示词要求报告 NOT CONVERGED，不为通过这些指标继续调弱材料。

## 保留与身份

工程没有 Git 仓库，使用源码快照与 SHA-256 确认身份；未发现适用 AGENTS.md。受保护的 {len(protected)} 个旧源码/配置/结果/动画文件全部未变，兼容链接保留。旧 mild、forced-failure 数据未覆盖，BraVa 旧微泡批次及其他停止任务未启动。每个新算例保存配置、标定、源码快照及其哈希。

证据：[保护核对]({VER}/FINAL_PRESERVATION.json)、[旧结果逐数组重放核对]({VER}/LEGACY_REPRODUCIBILITY.json)、[48 项测试]({VER}/FINAL_TESTS_SCOPED.log)、[新文件与全部输出绝对路径索引]({VER}/SOURCE_OUTPUT_INDEX.md)。第一次全目录 pytest 误收集了受保护的测试副本，报同名模块冲突，保留在 FINAL_TESTS.log；指定 `tests` 后 48 项通过，未删除快照。40 条提示来自 VTK/NumPy 数组 shape API 弃用，不是数值测试失败。

原冻结入口保留原样；新兼容入口 `pd_clot.regularization.legacy` 对 D_break<1 发出警告并保存 `forced_failure_verification=true`。新科学候选强制 D_break=1。旧代码中的 Q_ref 等字段仅因复制配置而保留，新损伤模式不使用它们。

## 几何、场与位置

直管 R=1.25 mm、长度 6 mm、流量 18 mL/min、黏度 0.00345312 Pa·s、流体/固体密度 1056 kg/m³。矩形壁面血栓尺寸 1.25×0.75×0.75 mm，总体积 0.703125 mm³；G=1000 Pa、K=50000 Pa，沿用原 NOSB 材料与稳定化实现。论文原件见 [PDF]({ROOT}/references/1-s2.0-S2666496826000567-main.pdf)，原核验说明见 [PAPER_READING.md]({ROOT}/references/PAPER_READING.md)；本轮不是原论文公式的完全复现。

中心现在为 **(-0.825, 0, -0.725) mm**，位于 clot 左面 x=-0.625 mm 上游 0.2 mm，与 clot 的 y/z 中点齐平。这是新计算的场参数，不是移动旧动画图标。白色空心菱形仅标记制造场的局部化中心；它不是顶部向下的点力、实际微泡表面或牵引方向。

`ManufacturedStreamingProvider` 使用 A_y=-UL exp(-r²/(2L²)) 的旋度，局部 u=U exp(-r²/(2L²))(-dz/L,0,dx/L)，叠加原 Poiseuille u/p/grad；U=0.01 m/s，L=0.35 mm，时间系数 1+0.9 sin(2π·500t)，局部压力幅度为 0。局部场解析散度为零，具有涡旋结构和空间衰减。牵引统一由 t=[-pI+μ(grad u+grad uᵀ)]n 计算，脱落后使用同一场的速度。这是人为构造的验证场，没有求解微泡边界或 Navier–Stokes 方程；1 MHz/2 µm 兼容配置字段不构成已求解的超声或泡半径效应。

## A. 断裂能与损伤律

用户授权使用演示值，选取 **Gc_demo=0.01 J/m²**，不是实验标定值。每档使用相同 0.5 mm 立方 coupon：两半刚性块位移控制循环张开、对称弱平面，仅跨平面键按同一损伤律演化，达到 D=1 才失效；没有手动切键。投影裂面面积 2.5×10⁻⁷ m²（不是两个面的面积之和）。标定通过标量求根调整 ψcrit，使实际 NOSB 储能下降累计/裂面面积达到 Gc_demo。以独立六点 Gauss 力–位移积分核验外功。低 C_E 扫参允许延续相同幅度增量最多三倍步数，直到跨平面键全部失效；粗、中、细基准都保持原提前完成的路径。

{mesh_table}

coupon 外功相对残差范围约 2.15×10⁻¹⁵–2.85×10⁻¹⁴；目标耗散约 2.5×10⁻⁹ J。配置容差为 5%，实测误差远小于该容差。W_cycle,ij=0.5 E_young a_ij² V_ij，a 是循环键伸长的半极差；V_ij 从初始加权邻域分配并满足求和等于材料体积。Wcrit,ij=ψcrit V_ij。ΔD=ΔN C_E[max(W_cycle/Wcrit−threshold,0)]^m_E；基准 C_E=0.001、m_E=1、threshold=0，无非零疲劳门槛。损伤不可逆，D_break=1。spacing/horizon/material/Gc/损伤参数不一致时拒绝复用标定。

这一能量驱动是方向应变能代理，不是 NOSB 自由能的精确逐键分解。指定 coupon 的能量一致不能推出主流场加载下的网格客观性。每个宏步实际储能降为损伤耗散；按 ΔD×驱动分配到各键的累计耗散是诊断估计。每次真正失败会写 failure_ledger.jsonl 和 fracture_energy_ledger.jsonl（失效周期、之前损伤、驱动、累计耗散、h、horizon、支撑秩）；本轮主算例没有失败，所以两个失败账本为空，这是正确结果。

## B. 网格敏感性与体积

三档 horizon/h=3.01、同一物理几何/材料/流场/Gc、25,000 代表周期；每档单独标定。所有颗粒保留，最终均连接固定基底，全部支撑秩为 3。未观测到首次失效或脱落的 N 在 JSON 为 null，不是 N=0。

{volume_table}

**固定层离散是明确的混杂因素。** 统一物理规则为 z<origin_z+0.125 mm，中心恰位于上边界者排除（1e-10 h 舍入容差）；粗/中/细固定体积分数为 1/6、1/9、1/6。因此尽管连续几何相同，边界约束积分并非相同，不能将所有差异归因于损伤律。10% 仅为描述性比较容差，不是数学收敛证明。脱落体积都为零不能证明碎片或断裂收敛。

## C. 周期跳跃与时间步长

{sensitivity_table}

ΔN=1000/500/250 不改变材料或重新标定；在相同代表 N 比较，并检查损伤单调。由于没有事件，首次断键/脱落的周期偏移无法计算，不能写成“零偏移”。dt/2 覆盖完整 25,000 周期，平均损伤差约 0.0014%，但没有包含断裂的时间窗，因此只能支持当前未断裂响应的时间步敏感性判断。

周期跳跃的实际模拟代表力学时长分别为 0.052/0.102/0.202 s，包含初始暖启动；改变 ΔN 也改变实际积分时长。25,000/500=50 s 不能代替当前已模拟的流固运动时间，外功不能凭空乘 ΔN。后续碎片输运的物理时间映射仍需要独立解决。

## D. 局域性

参考坐标中以中心半径 R=0.525 mm 定义影响区域，键用参考中点。没有冻结/删除远端键来限制损伤。

{locality_table}

粗档约 17.78% 材料所在区域内集中了 50.06% 的损伤，可描述为损伤对该区域偏聚；仍有约一半损伤在外部，不能称为全部局限于近源区域。断键比例和脱落来源比例因没有对应事件而为 null。

## E. 能量账本

{energy_table}

残差定义 R=U_elastic+U_stabilization+K−E_initial+D_damping+D_damage−W_surface−W_transport。外功沿真实代表积分路径梯形积分，阻尼/输运按实际动能改变量记账。粗档末态相对残差 -0.012664%，dt/2 为 -0.003166%；小残差是代码账本证据，不能证明局部力学或 CFD 物理精度。全时程见 energy_history.png/PDF 和 history.csv。

## F–G. 碎片、低秩、输运

诊断阈值为至少 8 颗粒、体积≥1.5625×10⁻¹¹ m³、至少 12 条内键、至少 80% 体积支撑秩≥2。它们是数值分辨要求，不是物理血栓尺寸阈值。单颗粒单独分类，不合并、不删除；其他小连通块为 UNDER_RESOLVED_DEBRIS。P_singleton>20% 或 P_lowrank>30% 会标记 FRAGMENT_RESOLUTION_INADEQUATE 并发出警告；低秩 intrinsic fallback 仅是带诊断的数值延续。新算例 P 值均未定义，所有脱落、resolved/debris/singleton/低秩体积与各类清除量均为零，没有可报告的新碎片尺寸分布。

旧结果重新诊断得到 P_singleton=78.899%、P_lowrank=82.569%，展示颗粒数、体积、等效直径三个分布。新旧加载强度、位置和损伤律均不同，不能把“新算例没有碎片”宣称为粉碎问题已改善。

可选 `transport.mode=fragment_drag` 对 resolved 连通块计算质量、体积、质心、等效球直径、迎风面积、平均同场流速和 Re，使用低 Re Stokes 阻力及 Schiller–Naumann 修正，Re≥1000 使用 Cd=0.44 的球形代理。每半步冻结系数的解析松弛，按质量分配总冲量，保持内部相对速度。参考公式：[OpenFOAM Foundation 源码](https://cpp.openfoam.org/v12/SchillerNaumann_8C_source.html)。未实现力矩、不适用于定量预测不规则可变形血栓。低分辨碎屑仍用旧粒子松弛且分开标记。主算例默认保留 relaxation，因无脱落未触发输运。

独立预分类连通块测试对照 ODE，dt=10/5 µs 的质心速度误差为 1.13×10⁻⁸/5.65×10⁻⁹ m/s，守恒检查通过；这是模块验证，没有在主模型中人为制造碎片。resolved/under-resolved 清除量按首次越线时类别记录且相加等于总量，越过清除平面不能称为临床溶栓。

## 参数敏感性

{sweep_table}

Gc、m_E、C_E 改变时按新值重新执行 coupon；速度改变时复用同一能量标定。因此 m_E/C_E 图是**保持目标 coupon Gc 的联合标定比较**，不是固定 ψcrit 的单参数偏导。额外四角计算与基准/OAT 构成 3×3 Gc–U 网格，response_surface 显示全部九个采样点；没有连续拟合、择优保留或为碎裂调参。所有 18 组均无失效或脱落，最大粒子损伤约 0.00533（m_E=0.8）。每组 hash、事件、体积、能量残差保存在 [CSV]({VER}/comparisons/study_summary.csv)。

## H. 导入器与最小下一步

`ResolvedStreamingFieldProvider` 支持 VTU/legacy VTK 的线性四面体原单元 P1 插值，以及显式声明 convex_hull 的 VTP/CSV 采样插值。u/p 当前要求节点值。体数据缺 grad 时按 P1 一致计算；平面数据必须给完整三维 grad，不能推断法向导数。拒绝缺单位、NaN、越界、不可覆盖的表面、过大散度及不一致梯度，保存观察范围与可选范围阈值诊断；只声明 outward_solid 并不构成法向几何方向的证明。当前为静态快照，不支持跨时间插值，非凸表面/孔洞不适用 convex_hull 近似。

四种格式的独立仿射测试通过；SI 最大误差 u=8.67×10⁻¹⁹ m/s、p=1.14×10⁻¹³ Pa、grad=7.22×10⁻¹⁶ s⁻¹。另一个弱均匀合成 VTU 已通过实际 PD runner 一宏步接入审计，见 [接口报告]({VER}/resolved_import_with_runner/IMPORT_VERIFICATION.json) 和 [runner 审计]({VER}/resolved_import_with_runner/runner_case/AUDIT.json)。它们不是实际微泡 CFD。

最小下一步：取得一份覆盖血栓全部受载表面及预期运动范围的**实际微泡解析体网格 VTU 快照**（节点 u、p，注明压力基准和坐标配准，最好同时提供 grad），给出以下显式元数据，在复制的新配置内设 `streaming.provider="resolved_field"`、`field_path`、`metadata_path`，先做导入诊断和一宏步对照；无需修改 PD 求解器。静态平均场的循环驱动含义仍需重新定义，不能把静态导入自动解释为超声疲劳。

{metadata}

材料、疲劳律均未实验/声致溶栓标定；制造场不是包膜微泡 CFD；尚无双向流体–PD 耦合；表面面积/法向、固定层离散、稳定化与低秩延续仍是精度限制。完成软件开发不等于完成所要求的科学验收。

## 动画、图片与复现

[本地可视化入口]({VIS}/index.html) · [MP4]({VIS}/regularized_damage_transport.mp4) · [GIF]({VIS}/regularized_damage_transport.gif) · [QC]({VIS}/QC.json)。视频 1920×1080、60 fps、751 帧、12.5167 s，黑底、Arial 常规字重、单损伤主图、无箭头/大蓝泡/Force 子图/N 顶标、固定相机与真实位移。源标记静止仅表示空间位置，不用其亮度伪造载荷相位。色标固定 0–0.0005，黄色不是 D=1；粒子显示损伤是邻域键完整度推导的值，断键阈值作用在键 D 上。

26 个原始状态每隔 30 显示帧严格对齐，位置只做线性插值、损伤保持前一宏步原值；固定点逐位不变，全部粒子保留。运动最大仅约 0.261 µm，某些相邻解码帧相同并不代表编码掉帧；60 fps 时间戳恒定，损伤宏步跳变也没有被虚构中间科学状态掩盖。四阶段图用 N=0/8000/17000/25000，C/D 明确写未到断裂/输运，局部流线均为代表峰值相位的说明叠图。

{figure_links}

复现命令（新目录存在则换名；不会覆盖未完成结果）：

{commands}
'''
    en=f'''# Energy-regularization development report

The existing clot_clearance project was extended with an isolated namespace. Eighteen study runs, verification modules, reports and a new movie are complete. **No main run reached bond failure or detachment; the model is NOT CONVERGED.** Parameters were not weakened to produce a fragmentation movie.

## Status and acceptance

SOFTWARE VERIFIED for the tested implementation: 48 tests pass, all 18 run audits pass, and every saved legacy forced-failure array is bitwise reproducible. ENERGY REGULARIZED applies to the specified numerical opening coupon and damage-energy scale only. It does not establish fracture-energy objectivity under general loading. Main mean damage changes by −56.37%/−67.44% and dissipation by −65.19%/−76.08% on medium/fine versus coarse spacing, hence NOT CONVERGED. New fragment resolution is NOT_ASSESSED_NO_DETACHMENT; 0/0 pulverization ratios are null, not zero. The independently reclassified legacy case is FRAGMENT_RESOLUTION_INADEQUATE (singleton 78.899%, low rank 82.569%). Different old/new loads prevent a matched comparison.

The requested first-detachment acceptance criterion, a time-step window containing first failure, fragmentation localization, finite-fragment dominance and transport convergence are **not satisfied**. Absence of events was retained rather than changing material parameters to pass.

## Preservation and provenance

No Git repository or applicable AGENTS.md was found. All {len(protected)} protected pre-existing source/config/result/movie files match their SHA-256 values; compatibility links and unrelated stopped BraVa jobs remain untouched. Every new run retains its configuration, coupon identity, source snapshot and hashes. [Source/output index]({VER}/SOURCE_OUTPUT_INDEX.md), [preservation]({VER}/FINAL_PRESERVATION.json), [legacy replay]({VER}/LEGACY_REPRODUCIBILITY.json), [tests]({VER}/FINAL_TESTS_SCOPED.log). The initial unrestricted pytest discovery also collected preserved test copies and failed module-name collection; the explicit `tests` suite then passed. VTK emits 40 NumPy shape deprecation warnings. No failed log or source snapshot was removed.

The frozen legacy entry is unchanged. The new `pd_clot.regularization.legacy` wrapper warns and labels D_break<1 as forced_failure_verification. New regularized runs reject any D_break other than one. Inherited Q_ref and other legacy-law fields are not used by the new damage mode.

## Configuration and consistent field

Pipe radius 1.25 mm, length 6 mm, Q=18 mL/min, μ=0.00345312 Pa·s, density 1056 kg/m³. Clot dimensions 1.25×0.75×0.75 mm, volume 0.703125 mm³; original G=1000 Pa, K=50000 Pa and NOSB stabilization retained. The [provided paper]({ROOT}/references/1-s2.0-S2666496826000567-main.pdf) and [prior reading notes]({ROOT}/references/PAPER_READING.md) remain references, not a claim of exact paper reproduction.

The new center is **(−0.825,0,−0.725) mm**, 0.2 mm upstream of the left face at x=−0.625 mm, vertically centered. It is an actual new field input. A hollow diamond marks the localization center, not a downward point load or a resolved bubble surface.

ManufacturedStreamingProvider takes curl(A), A_y=−UL exp(−r²/(2L²)); local u=U exp(−r²/(2L²))(−dz/L,0,dx/L). U=0.01 m/s, L=0.35 mm, modulation 1+0.9 sin(2π·500t), local pressure amplitude zero. The analytic divergence-free vortex is added to the existing pipe u,p,grad. Both t=[−pI+μ(grad u+grad uᵀ)]n and transport use this same field. No Navier–Stokes or coated-bubble problem is solved. Legacy 1 MHz/2 µm fields do not imply resolved ultrasound/bubble physics.

## A. Fracture-energy calibration

The user authorized a demonstration value: **Gc_demo=0.01 J/m²**, not an experimental value. Each spacing uses a separate 0.5 mm cubic, displacement-controlled two-rigid-half opening coupon. Only the symmetric weak-plane crossing bonds accumulate cyclic damage; no bond is manually cut. Projected crack area is 2.5e−7 m², counted once. A scalar root solve calibrates ψcrit against actual NOSB potential-energy loss divided by that area; independent six-point Gauss force integration measures external work. Low-C_E sweeps may extend the same opening increment to three times the initial step budget until the weak plane fails, without changing the baseline completed path.

{mesh_table}

Independent work residuals are 2.15e−15–2.85e−14, with about 2.5e−9 J dissipated. The configured Gc tolerance is 5%. W_cycle,ij=0.5 E_young a_ij² V_ij, where a is half the stretch range and original weighted-family partition volumes sum to material volume. Wcrit,ij=ψcrit V_ij. ΔD=ΔN C_E max(W_cycle/Wcrit−threshold,0)^m_E, baseline C_E=.001, m_E=1, threshold=0 (no nonzero fatigue threshold). Damage is irreversible and complete at D=1. Calibration reuse rejects inconsistent spacing, horizon, material, fatigue parameters or Gc.

This is an interpretable directional strain-energy proxy, not an exact NOSB per-bond free-energy decomposition. Coupon success is protocol-specific. At fixed coordinates the actual potential decrease is accumulated; per-bond attribution proportional to ΔD×driver remains an estimate. Failure and fracture-energy JSONL ledgers store event cycle, previous damage, driver, allocated dissipation, spacing, horizon and ranks. They are empty in these main runs because no failures occurred.

## B. Spacing sensitivity

All spacings keep physical geometry, material, source/field, Gc, total N=25,000 and horizon/h=3.01. All particles remain attached and rank three. Null event N means not observed by the exposure endpoint.

{volume_table}

The fixed slab uses z<origin_z+0.125 mm, excluding centers exactly at its upper boundary with a 1e−10 h roundoff tolerance. Fixed volume is 1/6,1/9,1/6 on coarse/medium/fine. Thus the boundary discretization is a disclosed confound even with identical continuous geometry. The 10% comparison tolerance is descriptive, not a mathematical convergence criterion. Zero detached volume does not validate fragmentation convergence or size distributions.

## C. Cycle-jump and time-step sensitivity

{sensitivity_table}

No material recalibration occurs between ΔN=1000/500/250. Damage is monotonic and histories are compared at common represented N. Failure/detachment shifts are unobservable, not zero. dt/2 spans the full exposure; it supports sensitivity of this unfractured response only. Mechanical proxy durations are 0.052/0.102/0.202 s, including warmup. Changing ΔN also changes simulated mechanical time; N/f=50 s is not the integrated transport time. External work is never multiplied by ΔN.

## D. Localization

The reference-coordinate sphere has R=0.525 mm; bonds use reference midpoints. Distant bonds are neither deleted nor frozen.

{locality_table}

At coarse spacing, 50.06% of damage lies in 17.78% of the material volume. This is preferential localization, not confinement of all damage. Failure and detached-origin localization remain undefined.

## E. Energy

{energy_table}

R=U_elastic+U_stabilization+K−E_initial+D_damping+D_damage−W_surface−W_transport. Surface work follows actual trapezoidal force/displacement integration; damping/transport use actual kinetic-energy changes. Final coarse residual is −0.012664%, dt/2 −0.003166%. Small residuals verify bookkeeping, not local physical accuracy. Energy histories are supplied in CSV and PNG/PDF.

## F–G. Fragments and transport

Resolved-fragment diagnostics require at least 8 particles, volume 1.5625e−11 m³, 12 active internal bonds and 80% volume with rank≥2. These are numerical thresholds, not physical clot-size thresholds. Singletons are a separate category; other unsupported components are debris. No particles are deleted/merged. P_singleton>.20 or P_lowrank>.30 triggers a visible inadequacy warning. Intrinsic reduced-rank continuation is diagnosed, not silently accepted as 3-D NOSB material. All new detached/clearance volumes are zero; no new fragment-size distribution exists. Legacy count/volume/equivalent-diameter distributions are shown with different-loading labels.

Optional fragment_drag computes component volume/mass/COM, equivalent sphere diameter/area, same-field velocity sampling and Re. It uses Stokes drag with Schiller–Naumann correction below Re=1000, and Cd=.44 above; see [OpenFOAM Foundation source](https://cpp.openfoam.org/v12/SchillerNaumann_8C_source.html). Frozen-coefficient COM relaxation distributes total impulse conservatively by mass; no torque is implemented. Under-resolved debris retains labeled per-particle relaxation. The main case retains relaxation and has no detached transport to exercise. An independent preclassified component ODE comparison gives 1.13e−8/5.65e−9 m/s errors at dt=10/5 µs, without manufacturing a fragment in the PD simulation. Clearance is classified once at first crossing; plane crossing is not clinical thrombolysis.

## Parameter sensitivity

{sweep_table}

Gc/m_E/C_E changes each trigger a new coupon calibration; U changes reuse the same scale. m_E/C_E curves therefore compare jointly recalibrated models at fixed target coupon Gc, not partial derivatives at fixed ψcrit. Four extra corner runs complete the 3×3 Gc–U sample grid. All outcomes are retained, no response fit or visual selection is used. No run fails/detaches; largest particle D≈.00533 at m_E=.8. [All configurations, hashes and metrics]({VER}/comparisons/study_summary.csv).

## H. Importer and limitations

ResolvedStreamingFieldProvider supports original linear-tetrahedral VTU/legacy VTK P1 interpolation and explicitly declared convex-hull VTP/CSV interpolation. Velocity and pressure are nodal. Missing volume gradients are consistently derived; planar input must supply a full 3-D gradient. Explicit units are mandatory. NaN, out-of-domain queries, missing coverage, excessive divergence and gradient inconsistency are rejected. Observed ranges and optional range thresholds are reported. Outward-solid normal metadata is a declaration, not geometric proof. Static snapshots only; no temporal interpolation or nonconvex scattered-domain reconstruction.

Independent four-format affine tests pass, maximum SI errors u=8.67e−19, p=1.14e−13, grad=7.22e−16. A separate weak uniform synthetic VTU also passed one actual PD macro step through the config-selected importer; [verification]({VER}/resolved_import_with_runner/IMPORT_VERIFICATION.json), [run audit]({VER}/resolved_import_with_runner/runner_case/AUDIT.json). These fixtures are not real microbubble CFD.

Smallest next step: obtain one actual resolved-microbubble tetrahedral VTU covering the loaded surface and anticipated motion, with nodal u,p, explicit pressure reference/coordinate registration and preferably grad; provide the metadata below. In a copied config select `streaming.provider="resolved_field"`, `field_path`, `metadata_path`, then inspect import diagnostics and one macro step. The PD solver need not change. The meaning of cyclic forcing from a static mean field still needs definition.

{metadata}

Clot material and fatigue law are not experimentally/sonothrombolysis calibrated. Manufactured streaming is not resolved microbubble flow, and no two-way fluid–PD feedback is solved. Surface quadrature/normals, fixed-slab discretization, stabilization and low-rank continuation remain limitations.

## Visual deliverables and reproduction

[Viewer]({VIS}/index.html), [MP4]({VIS}/regularized_damage_transport.mp4), [GIF]({VIS}/regularized_damage_transport.gif), [QC]({VIS}/QC.json). Black/Arial regular, single damage view, no arrows/large blue bubbles/Force panel/N header. Fixed camera, displacement scale one. 1920×1080, 60 fps, 751 frames, 12.5167 s. Diamond position is fixed and does not encode phase. Damage colorbar is explicitly 0–.0005; yellow is not complete damage one. Particle damage derives from neighboring bond integrity; D_break applies to bonds.

Every 30 display frames exactly matches one of 26 original states. Positions alone are linearly interpolated, damage holds its recorded macro value; all particles and exact anchors are preserved. Tiny real displacements (~.261 µm maximum) can produce identical adjacent decoded frames. Constant video cadence does not turn display interpolation into new solver samples. Four-stage panels use N=0/8000/17000/25000 and explicitly state that failure/transport stages were not reached; streamlines are an illustrative representative-peak-phase overlay.

{figure_links}

Use fresh output directories. Current-suite testing explicitly excludes preserved copies:

{commands}
'''
    (ROOT/'REGULARIZATION_REPORT_ZH.md').write_text(zh)
    (ROOT/'REGULARIZATION_REPORT_EN.md').write_text(en)
    cards=''.join(f'<figure><a href="figures/{p.name}"><img loading="lazy" src="figures/{p.name}" alt="{html.escape(p.stem)}"></a><figcaption>{html.escape(p.stem.replace("_"," "))} · <a href="figures/{p.stem}.pdf">PDF</a></figcaption></figure>' for p in figures)
    webpage=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Particle damage</title>
<style>body{{background:#000;color:#eee;font:18px Arial,sans-serif;margin:0 auto;padding:28px;max-width:1500px}}h1,h2{{font-weight:400}}a{{color:#9ce0ea}}video{{width:100%;max-height:82vh;background:#000}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,480px),1fr));gap:28px}}figure{{margin:0}}img{{width:100%}}p{{line-height:1.5;max-width:1000px}}figcaption{{font-size:15px;color:#c1cad5}}button,select{{font:inherit;background:#101923;color:white;border:1px solid #526278;padding:5px}}nav{{display:flex;gap:24px;flex-wrap:wrap}}</style>
<h1>Particle damage</h1><video id="movie" controls loop playsinline preload="metadata" poster="preview.png"><source src="regularized_damage_transport.mp4" type="video/mp4"></video>
<label>Playback speed <select id="rate"><option value="0.5">0.5×</option><option selected value="1">1×</option><option value="2">2×</option></select></label>
<p>Actual damage range: 0–0.0005. Yellow does not mean complete failure. The hollow diamond marks the manufactured streaming region 0.2 mm upstream of the clot. No bond failure or detachment occurred in this run.</p>
<p>Display: 60 fps, true displacement, fixed camera. Position interpolation only; damage retains recorded macro-step values. The manufactured field is a verification model. Mesh convergence and physical fragmentation have not been established.</p>
<nav><a href="regularized_damage_transport.mp4">Download MP4</a><a href="regularized_damage_transport.gif">GIF</a><a href="../../REGULARIZATION_REPORT_ZH.md">中文报告</a><a href="../../REGULARIZATION_REPORT_EN.md">English report</a><a href="QC.json">Video checks</a><a href="figures/data/study_summary.csv">Study data</a></nav>
<h2>Study figures</h2><div class="grid">{cards}</div><script>document.getElementById('rate').onchange=e=>document.getElementById('movie').playbackRate=Number(e.target.value)</script></html>'''
    (VIS/'index.html').write_text(webpage)
    print(json.dumps(dict(reports=['REGULARIZATION_REPORT_ZH.md','REGULARIZATION_REPORT_EN.md'],new_source_files=len(source_files),figures=len(figures),protected_files=len(protected),runs=len(rows)),indent=2))


if __name__=='__main__':main()
