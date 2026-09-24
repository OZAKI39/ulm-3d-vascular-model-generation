"""Package an honest hard-stop report; empty downstream tables are not results."""
from pathlib import Path
import csv
import json
import subprocess
import hashlib
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
R = ROOT/'particle_3d/reports/particle9a3_interior_inlet'
D = R/'data'


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False)+'\n')


def main():
    gate = json.loads((D/'interior_injection_section.json').read_text())
    if gate['status'] != 'NO_VALID_INTERIOR_INJECTION_SECTION':
        raise RuntimeError('This report is only for a failed section search')
    rows = list(csv.DictReader((D/'interior_section_search.csv').open()))
    first = next(r for r in rows if r['refinement'] == '1')
    best = min(rows, key=lambda r: float(r['flux_relative_error']))
    source = json.loads((ROOT/'particle_3d/contracts/SONOVUE_DIAMETER_MAX4UM_V1.json').read_text())
    bins = source['bins']
    mass = np.array([b['retained_probability'] for b in bins])
    edges = [bins[0]['low_um']]+[b['high_um'] for b in bins]
    percentiles = np.interp([.1, .5, .9], np.r_[0., np.cumsum(mass)], edges)
    mean = float(sum(b['retained_probability']*(b['low_um']+b['high_um'])/2 for b in bins))
    source_summary = dict(contract='SONOVUE_DIAMETER_MAX4UM_V1', unchanged=True,
        support_um=[edges[0], edges[-1]], density='original piecewise-uniform histogram, conditioned on D<=4 um',
        mean_um=mean, P10_P50_P90_um=percentiles.tolist(), concentration_m3=8.5e12)
    write_json(D/'source_distribution_summary.json', source_summary)
    blocked = 'NOT_EXECUTED_SECTION_GATE_FAILED'
    tables = {
        'diameter_accessible_flux.csv': ['diameter_um', 'Q_feasible_m3_s', 'accessible_flux_fraction', 'status'],
        'entering_diameter_distribution.csv': ['diameter_um', 'source_pdf_per_um', 'entering_pdf_per_um', 'status'],
        'joint2000_events.csv': ['event_id', 'birth_time_s', 'diameter_um', 'x_m', 'y_m', 'z_m', 'section_triangle', 'status'],
        'joint2000_point_outcomes.csv': ['event_id', 'outcome', 'status'],
        'smoke30_outcomes.csv': ['event_id', 'outcome', 'stationary', 'penetration', 'handoff_violation', 'status'],
        'flux_integration_convergence.csv': ['diameter_um', 'Q_default_m3_s', 'Q_finer_m3_s', 'relative_difference', 'status'],
    }
    for filename, fields in tables.items():
        with (D/filename).open('w', newline='') as f:
            csv.writer(f).writerow(fields)
    write_json(D/'joint_basin_flux.json', dict(status=blocked, J_O1=None, J_O2=None, J_O3=None, J_no_exit=None,
        interpretation='Null means not calculated; it does not mean zero flux.'))
    raw_rate = 8.5e12*gate['reference_inlet_Q_m3_s']
    write_json(D/'injection_rate.json', dict(status=blocked, concentration_m3=8.5e12,
        Q_inlet_reference_m3_s=gate['reference_inlet_Q_m3_s'], old_C_times_Q_inlet_per_s=raw_rate,
        Q_raw_section_m3_s=None, Q_effective_joint_m3_s=None, accessible_number_flux_fraction=None,
        lambda_enter_per_s=None, mean_event_interval_s=None, first_500_window_s=None))
    write_json(D/'downstream_artifacts_status.json', dict(status=blocked,
        empty_csv_files=list(tables), empty_means='Not computed because the required interior section does not exist in this search.',
        missing_selected_vtp='interior_injection_section.vtp is intentionally absent; rejected_section_witness_*.vtp are diagnostic only.',
        figures_not_generated=['02_accessible_flux_vs_diameter', '03_section_basin_and_feasible_regions',
                               '04_injection_population_summary', '05_synthetic_circular_tube', '06_method_comparison'],
        not_implemented_after_stop=['upstream finite-size connectivity certificate', 'adaptive Q_feasible(D) integration',
                                  'joint size/position sampler', 'physical finite-size schedule', 'joint basin audit',
                                  'synthetic circular Poiseuille finite-size benchmark', 'downstream permanent tests']))
    tests = json.loads((D/'test_results.json').read_text())
    protection = json.loads((D/'protected_verification.json').read_text())
    write_json(D/'gates.json', dict(status=gate['status'], human_review='PENDING', production_ready=False,
        valid_section=False, section_flux_validation='FAIL', independent_p1_integration='PASS',
        implemented_section_unit_tests=tests, protected_inputs=protection,
        finite_size_upstream_connectivity='NOT_RUN_SECTION_GATE_FAILED',
        joint_population_tests=blocked, audit2000=blocked, smoke30=blocked, formal500='NOT_AUTHORIZED_THIS_STAGE'))
    report = f'''# 一句话结论

**NO_VALID_INTERIOR_INJECTION_SECTION。2D 思想向 3D 升级的前置截面检查未通过，因此本轮尚未产生新的 entering population。** 已开发真实四面体切面与流量验证代码，在服务器搜索了 {gate['primary_candidates']} 个主候选和 {gate['refined_candidates']} 个加密候选，没有合格截面。严格按本轮“无有效截面即 STOP”的要求，未继续实现下游 joint sampler，未生成 2000 个事件、30 条轨迹或新的正式 500。人工审核为 **PENDING**。

# 为什么不再在 open inlet cap 生成完整球

人工入口 cap 相当于把水管截断后留下的计算边界；它不应被误认为真实的挡板。内部截面则是在根血管里设置一个观察面，统计已经从上游过来的有限尺寸微泡。这个想法必须先有一个位置合法、血流量也与入口相符的截面。本轮尚未过这一关，不能把后续的粒径加权公式当成已验证实现。

已读取旧 Method C 审核报告与指定 2D Git 源文件。实际 blob 与要求一致。2D 到 3D 的对应及不能直接移植的部分见 [TWO_D_TO_THREE_D_MAPPING_ZH.md](TWO_D_TO_THREE_D_MAPPING_ZH.md)。特别是旧 2D 的“connected clearance”函数只查询局部壁距，不能直接充当 3D 上游可达证明。

# 新的内部截面在哪里

**没有选中的正式截面。** `data/interior_injection_section.json` 的 `selected_section=null`，不会把失败候选冒充正式入口，也没有输出一个误导性的正式 `interior_injection_section.vtp`。

几何来源经原 FEM `source_contract.json` 追溯至同一 ROI 的 `roi_core.swc`。实际入口 cap 面积质心连接到保存的 ROI root；这一段属于原血管几何已有的人工延长管段，其来源有 `port_classification.csv` 支持。其后沿原 SWC 节点 1–15，首个分叉为节点 **15**，子节点为 16、17。没有采用找不到中心线时的 inward-normal fallback。

- 实际入口中心：({gate['inlet_center_m'][0]*1e6:.6f}, {gate['inlet_center_m'][1]*1e6:.6f}, {gate['inlet_center_m'][2]*1e6:.6f}) µm。
- 首个分叉位置：(92, 49, 111) µm。
- 沿 root 到首个分叉的弧长：**{gate['root_length_to_first_junction_m']*1e6:.6f} µm**。
- 主候选间距使用冻结求解策略中的 mesh `h10={gate['candidate_spacing_m']*1e6:.6f} µm`，每个区间中点为候选，法向取局部中心线切向；另用一半间距做一次分辨率检查。序列在读取结果之前就已确定，不随流量或球径调节。
- 第一个主候选距 inlet **{float(first['arclength_m'])*1e6:.6f} µm**，距首个分叉 **{float(first['distance_to_first_junction_m'])*1e6:.6f} µm**，但流量失败。图中的橙色面明确标为 rejected。

无限延伸的平面有时还会切到远处另一条血管。本实现只取包含 authoritative root 中心点的连通分量作为实际有限截面，并同时记录整张无限平面的分量数；没有把远处血管流量加进 root。检查截面边缘确实位于 WALL、是否触及 cap、是否在首个分叉之前，以及 daughter 轴线是否穿过同一截面。**这些是本轮前置几何检查，不等同于有限球上游连通性证书；后者未执行。** 主候选中 {gate['primary_geometry_pass_count']} 个通过这些几何检查，4 个同时触碰开放边界。所有 204 个均因流量误差被拒绝。

候选搜索是在声明的中心线候选族内进行的，不声称数学上排除了所有可能平面。没有为了逼近 cap 获得偶然数值吻合而无限缩小偏移，也没有越过首个分叉或为了容纳更大微泡移动截面。

![真实血管与被拒绝截面](figures/01_interior_injection_section.png)

图 01：左图为真实冻结血管、对应源 root 中心线和首个分叉；右图放大首个主候选。橙色截面未通过流量检查，不是选中的正式 injection section。原要求中“亮色选中截面”的图意在本轮改为明确展示失败证据。

# 这个截面的流量对不对

采用与现有微泡 `FrozenFEMField` 相同的 canonical tetra 点位、拓扑和 float64 节点速度。原文件坐标存储为 float32，读入后提升为 float64 计算，不改变几何数值。VTK 切割保留 double 坐标，P1 法向速度在切割三角形上精确积分；有负速度时先精确切分正半区。

独立重算 inlet cap 得到 **{gate['reference_inlet_Q_m3_s']:.15e} m³/s**，与保存 inlet Q 一致。采用原 `policy.json` 的 `mass_limit=1e-6`，即 **0.0001%** 相对容差；这是依本轮要求沿用的原全局流量验证阈值，不是新发明的宽松内部容差。

| 对象 | 距 inlet 的 root 弧长 (µm) | 正向 Q (m³/s) | 相对 inlet 误差 | 结果 |
| --- | ---: | ---: | ---: | --- |
| 保存并独立复算的 inlet | 0 | {gate['reference_inlet_Q_m3_s']:.12e} | 0 | reference |
| 最佳主候选（也是第一个主候选） | {float(first['arclength_m'])*1e6:.6f} | {float(first['positive_Q_m3_s']):.12e} | {float(first['flux_relative_error'])*100:.6f}% | FAIL |
| 加密搜索中最佳候选 | {float(best['arclength_m'])*1e6:.6f} | {float(best['positive_Q_m3_s']):.12e} | {float(best['flux_relative_error'])*100:.6f}% | FAIL |

加密后最小误差仍为容差的 **{gate['best_flux_error_refined']/gate['flux_relative_tolerance']:.1f} 倍**。另对首个、最佳、中部、末端及加密首个候选执行独立 NumPy 四面体边相交/多边形积分（首个与最佳重复，是 5 个审计记录、4 个不同位置），绕过 VTK 切割；与 VTK P1 积分的最大差为 inlet Q 的 **{gate['maximum_independent_method_difference']:.3e}**，远小于此处的流量不匹配。原始记录见 `data/independent_section_flux_validation.json`。

**全局 inlet/outlet 总流量平衡，不保证保存的 P1 插值速度在任意内部截面上都给出相同流量。** 这里检查的是微泡实际使用的冻结 P1 场。本轮证据说明它无法通过本任务要求的内部截面流量门槛，不能仅凭这些数据确定责任在求解、导出还是冻结表示的哪个环节。没有缩放内部速度、平滑梯度、修改 FEM 或提高容差来让检查通过。

![全部候选流量与误差](figures/07_interior_section_flux_gate.png)

图 07：上图显示候选真实 P1 流量与 inlet 的比值；下图用对数轴显示相对误差和原容差。靠近分叉的候选还会受到几何截面与开放边界接触影响；每一项拒绝原因均保存在 CSV，不能只看上图曲线选择新入口。

# 不同尺寸真正能使用多少血流

**未计算。** D=1.0/1.5/2.0/2.5/3.0 µm 的 accessible flux fraction 均为 N/A；N/A 不是 0。截面最大可达刚性球直径也未计算，不能与旧 open-cap **3.1281 µm** 宣称改善或退步。不得在失败截面上继续算一个可用的粒径分布。

# source 与真正 entering size distribution 有什么区别

冻结 source 合同保持 `SONOVUE_DIAMETER_MAX4UM_V1`，支撑区间 0.75–4 µm，保留原直方图的分箱内均匀密度，并对 D≤4 条件化；未 clipping、未改原始 SonoVue。由现有合同直接计算的 source 平均直径为 **{mean:.6f} µm**，P10/P50/P90 为 **{percentiles[0]:.6f}/{percentiles[1]:.6f}/{percentiles[2]:.6f} µm**。

计划的新 entering 分布应当是 source PDF 乘 `Q_feasible(D)` 后归一化，不能继续用 Method C 对每个尺寸反复找位置强制成功。但**本轮 entering PDF 尚未建立**，源分布统计不能当作 entering 统计。

# 实际进入率是多少

源浓度保持 **8.5×10¹² m⁻³**（modeled D≤4 µm source）。旧 raw 参考 `C×Q_inlet` 为 **{raw_rate:.12f} s⁻¹**。新的 `C×∫f(D)Q_feasible(D)dD`、降低比例、平均间隔及正式 500 的物理时间窗均 **N/A**。没有调浓度保持旧 event times。

# O1/O2/O3 的 joint accessible flux

**全部未计算，JSON 为 null。** 没有生产 joint sampler，也没有 outlet-basin 审核运行，不能把 null 当作 O1/O3 通量为零。没有使用 outlet quota。

# 2000 个入口样本结果

计划 2000，实际生成 **0**，point tracer 实际运行 **0**。未进行采样分布、worker 无关性或 joint basin 一致性验证，这些状态是 NOT_EXECUTED，不能写 PASS。

# 30 条真实轨迹

计划 30，实际运行 **0**。stationary、penetration、handoff violation、solver failure 和 outlet 统计均 N/A。零条运行不等于“30 条通过”或“零穿透验证通过”。正式 500 实际运行 **0**。

# 与 Method B/C 相比

| 历史入口 population | point O1 | point O2 | point O3 | no exit |
| --- | ---: | ---: | ---: | ---: |
| Method B accepted 500 | 9 | 469 | 20 | 2 |
| Method C audit 2000 | 1 | 1999 | 0 | 0 |
| Method C formal 500 的入口 point audit | 0 | 498 | 1 | 1 |
| P9-A.3 | N/A | N/A | N/A | N/A |

历史数字来自上轮 `PARTICLE9A2_INLET_REVIEW_ZH.md`，只作为背景，未重算或覆盖。本轮没有有效新入口 population，因此不能完成三代入口分布效果对照，也没有混入后续动力学出口统计。

# 目前说明血管几何什么问题

本轮**不能判断 P9-A.3 是否仍出现 O2 collapse，更不能把其原因归为 open-cap 或 rigid-size/root-geometry**，因为尚未执行联合通量和 point audit。已确定的阻断原因是“当前冻结 P1 场的内部截面流量不满足现有阈值”。若继续推进，需要先单独审查这一流量一致性问题；本轮禁止改变冻结物理，故在这里停止。

# 实际完成的代码、测试、数据和运行

新分支：`dev/particle9a3-interior-joint-flux-inlet`。新增 `interior_section.py`、`run_particle9a3_section_gate.py`、`render_particle9a3_section_gate.py`、本报告生成器，以及 3 个测试文件/1 个 fixture。没有修改现有科学模块。

前置阶段的 **7 项永久测试**在 WSL 和服务器均通过：首分叉前候选范围、分叉后拒绝、WALL 截面边界、触 cap 拒绝、常速解析面积通量、含正负速度的线性 P1 正通量精确积分，以及 1% 流量损失必须拒绝。它们只验证已实现的截面门槛；**不是全部 12 类 joint-model 测试已经完成，也不是 circular Poiseuille finite-size benchmark 已通过**。后续算法及测试按 STOP 条件暂未实现。

服务器目录：`{gate['execution_root']}`。使用 **6 workers**，OMP/OpenBLAS/MKL/NumExpr 均为 1。服务器截面搜索及独立积分的实测墙钟时间 **{gate['runtime_s']:.3f} s**，不包含准备、传输、写代码或画图时间；未使用 GPU。服务器结果与本地交付一致性校验另见 `data/server_delivery_verification.json`。

预先保护的 **{protection['checked_files']} 个文件 SHA 全部不变**；原 Git index tree、HEAD 不变，原工作区已有修改完整保留。没有 commit/push。`logs/git_diff_stat.txt` 保存真实 `git diff --stat`（它包含此前已有未提交修改，不包含本轮未跟踪文件）；本轮新文件统计另存 `logs/new_files_stat.txt`。

数据入口：`data/interior_section_search.csv`（612 行）、`data/interior_injection_section.json`、`data/root_topology.json`、`data/independent_section_flux_validation.json`、被拒绝截面的 VTP。所要求的下游 CSV 仅保留表头，相关 JSON 用 null 和 NOT_EXECUTED 明确表示未计算；参见 `data/downstream_artifacts_status.json`。没有创建虚假轨迹或虚假数值。

图 01、07 提供英文白底 300 dpi PNG/PDF；02–06 所需数据未产生，故不绘制假图。图 01 可由冻结 WALL/INLET、`root_topology.json` 和保存的 rejected VTP 重绘；图 07 直接从 search CSV/section JSON 重绘。图像已检查，状态文字与实际失败结果一致。复现命令见 [REPRODUCE.md](REPRODUCE.md)。

# 是否建议启动新的 500 条

**NOT_READY**。阶段状态为 **NO_VALID_INTERIOR_INJECTION_SECTION**，不是 production complete。没有合格截面就不继续运行 2000/30，更不启动新的正式 500。

# 人工审核

- [ ] interior section 位置合理
- [ ] section flux 合理
- [ ] finite-size accessible flux 合理
- [ ] 3D synthetic tube benchmark 合理
- [ ] entering size distribution 合理
- [ ] O1/O2/O3 joint flux 可接受
- [ ] 30 smoke 合理
- [ ] 同意下一轮生成正式 500
'''
    (R/'PARTICLE9A3_INTERIOR_INLET_REVIEW_ZH.md').write_text(report)
    (R/'REPRODUCE.md').write_text('''# 复现本轮截面门槛

当前状态是 `NO_VALID_INTERIOR_INJECTION_SECTION`。只有前置截面搜索/审核已实现；以下命令不会运行 2000 事件或任何微泡轨迹。

在仓库根目录（WSL）执行：

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=particle_3d/src
.venv/bin/python -m pytest particle_3d/tests/particle9a3_inlet -q
.venv/bin/python particle_3d/scripts/run_particle9a3_section_gate.py --workers 6
.venv/bin/python particle_3d/scripts/render_particle9a3_section_gate.py
```

服务器使用 `/root/particle8_2_runs/env/bin/python` 替换 `.venv/bin/python`。只在新的 P9-A.3 目录执行；搜索输出含实测 host/runtime，重新运行会更新该阶段的输出，不会覆盖旧 B/C 结果。审核报告生成器还需要保存的 `test_results.json` 和 `protected_verification.json`，它不会自行声称测试通过。

源几何追溯在 reference/source_contract.json 和 roi_core.swc。2D 参考 blob 与 SHA 在 data/two_d_source_identity.json。输入校验调用原 read_frozen，候选相对流量容差直接读取冻结 policy.mass_limit；不能用修改阈值的方式重现 PASS。
''')
    (R/'OPEN_RESULTS.html').write_text('''<!doctype html><html lang="zh"><meta charset="utf-8"><title>P9-A.3 截面审核</title>
<style>body{max-width:1100px;margin:32px auto;font:17px/1.6 system-ui;color:#203041}img{width:100%}code{background:#eee}a{color:#17679b}</style>
<h1>NO_VALID_INTERIOR_INJECTION_SECTION</h1><p>204 + 408 个候选均未通过原流量阈值。2000 入口事件、30 smoke、新正式 500 均未运行。人工审核 PENDING。</p>
<p><a href="PARTICLE9A3_INTERIOR_INLET_REVIEW_ZH.md">中文审核报告</a> · <a href="data/interior_section_search.csv">候选 CSV</a> · <a href="data/interior_injection_section.json">截面 Gate JSON</a> · <a href="REPRODUCE.md">复现</a></p>
<img src="figures/01_interior_injection_section.png" alt="血管与被拒绝截面"><img src="figures/07_interior_section_flux_gate.png" alt="截面流量失败证据"></html>''')
    (R/'logs/git_diff_stat.txt').write_text(subprocess.check_output(['git', 'diff', '--stat'], cwd=ROOT, text=True))
    paths = [ROOT/'particle_3d/src/particle_3d/interior_section.py',
             *sorted((ROOT/'particle_3d/scripts').glob('*particle9a3*.py')),
             *sorted((ROOT/'particle_3d/tests/particle9a3_inlet').glob('*.py'))]
    (R/'logs/new_files_stat.txt').write_text('\n'.join(f'{p.relative_to(ROOT)} | {len(p.read_text().splitlines())} lines' for p in paths)+
        f'\n{len(paths)} new Python files; {sum(len(p.read_text().splitlines()) for p in paths)} lines. Not staged; no existing scientific source modified.\n')
    write_json(D/'summary.json', dict(status=gate['status'], human_review='PENDING', section_search=gate,
        source_distribution=source_summary, old_raw_event_rate_per_s=raw_rate, new_event_rate_per_s=None,
        old_B_point_split=[9,469,20,2], old_C_audit2000_point_split=[1,1999,0,0], new_point_split=None,
        stationary_count=None, penetration_count=None, executed_trajectories=0,
        report='PARTICLE9A3_INTERIOR_INLET_REVIEW_ZH.md', figures=['figures/01_interior_injection_section.png','figures/07_interior_section_flux_gate.png']))
    print('Blocked-stage report and explicit NOT_EXECUTED downstream schemas written.')


if __name__ == '__main__': main()
