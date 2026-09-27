"""Generate the Chinese evidence report from frozen design and actual CFD records."""
from pathlib import Path
import argparse
import csv
import json
import re
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from network_1d0d.audit import sha256,write_json,write_csv


def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
        ['| '+' | '.join(map(str,row))+' |' for row in rows])


def run(case,source,design,report):
    case,source,design,report=map(Path,(case,source,design,report))
    read=lambda p:json.loads(p.read_text())
    d=read(design/'design_summary.json');pre=read(report/'final_preflight.json')
    audit=read(report/'independent_validation.json');assert audit['status']=='PASS'
    ex=read(case/'reports/execution.json');c=audit['comparison'];m=audit['measurements']
    h=read(report/'best_feasible_fem_handoff.json');pol=read(case/'policy.json')
    w=read(report/'wss_postprocess.json');assert w['status']=='PASS'
    b,o=d['baseline'],d['optimized'];ports=('O1','O2','O3')
    manifest=read(case/'frozen_flow/manifest.json')
    cert=ROOT/'reports/balanced_three_outlet_forward_v1/feasibility_certificate.json'
    assert sha256(cert)==d['infeasibility_certificate_sha256']
    assert sha256(ROOT/'network_1d0d/parameter_fit.py')==d['provenance']['code_sha256']['parameter_fit.py']
    testlog=(report/'logs/preflight_pytest.log').read_text()
    testline=re.findall(r'\d+ passed in [\d.]+s',testlog)[-1]
    with (design/'bound_sensitivity.csv').open() as f:bounds=list(csv.DictReader(f))
    with (report/'wss_comparison.csv').open() as f:ws=list(csv.DictReader(f))
    with (report/'final_3D_ports.csv').open() as f:port3d=list(csv.DictReader(f))
    fractions=[('0D baseline',[b['prediction']['outlet_flow_fraction'][p] for p in ports],b['metrics']),
        ('0D primary design',[o['prediction']['outlet_flow_fraction'][p] for p in ports],o['metrics']),
        ('旧 H0 3D（实际 VTU 积分）',c['old_H0_3D_fractions'],c['old_H0_3D_metrics']),
        ('新 3D（唯一 forward）',c['new_3D_fractions'],c['new_3D_metrics'])]
    ftab=table(['工况','O1 %','O2 %','O3 %','J','极差','总体 std','最大等分偏差'],
        [[label,*[f'{x*100:.8f}' for x in f],*[f'{mm[k]:.12g}' for k in ('J_balance','range','std','max_abs_deviation')]] for label,f,mm in fractions])
    write_csv(report/'all_balance_metrics.csv',[dict(case=label,**{p:v for p,v in zip(ports,f)},**mm) for label,f,mm in fractions])
    # Preserve the declared table-only form of bound sensitivity.
    btab=table(['s_max','最佳 s_O1','最佳 s_O2','O1/O2/O3 %','J','极差','active bounds','用途'],
        [[row['s_max_dimensionless'],f"{float(row['s_O1']):.10g}",f"{float(row['s_O2']):.10g}",
          '/'.join(f"{100*float(row['f_'+p]):.6f}" for p in ports),
          f"{float(row['J_balance']):.12g}",f"{float(row['range']):.12g}",
          row['active_s_O1']+'/'+row['active_s_O2'],row['role']] for row in bounds])
    htab=table(['端口','0D Q (m³/s)','real-cut P (Pa)','R_ext (Pa·s/m³)','raw cap P (Pa)','applied cap P (Pa)'],
        [[row['port'],f"{row['signed_Q_m3s']:.15g}",*[f'{row[k]:.15g}' for k in ('pressure_realcut_Pa','R_extension_Pa_s_m3','pressure_cap_raw_Pa','pressure_cap_shifted_Pa')]] for row in h['ports']])
    dtab=table(['端口','实际 Qout (m³/s)','实际 Qout/Qin %','3D−0D 百分点','cap 面积平均流体 P (Pa)'],
        [[r['port'],r['flow_m3_s'],f"{float(r['fraction'])*100:.8f}",f"{float(r['difference_from_0D'])*100:+.8f}",r['area_average_pressure_pa']] for r in port3d])
    wtab=table(['工况','区域','面片数','面积 μm²','面积均值 Pa','面积 P5/P50/P95 Pa','最小/最大 Pa'],
        [[r['case'],r['region'],r['face_count'],f"{float(r['area_um2']):.6g}",f"{float(r['area_mean_Pa']):.6g}",
        '/'.join(f"{float(r[k]):.6g}" for k in ('P5_area_Pa','P50_area_Pa','P95_area_Pa')),
        '/'.join(f"{float(r[k]):.6g}" for k in ('min_Pa','max_Pa'))] for r in ws])
    profile=(case/'run/petsc_profile.txt').read_text(errors='replace')
    gpu_lines=[line.rstrip() for line in profile.splitlines() if any(word in line for word in ('Using PETSc','GPU Mflop','GPU %F','MatMult ','KSPSolve ','PCApply ','MatLUFactorNum','PCSetUpOnBlocks'))]
    gpu_samples=[float(row['gpu'].split(',')[0]) for row in ex['GPU_samples']]
    highmem=max(float(row['gpu'].split(',')[1]) for row in ex['GPU_samples'])
    optimizer=o['optimizer'];cross=d['crosscheck']
    evidence_hashes={'frozen_balance_design.yaml':sha256(design/'frozen_balance_design.yaml'),
        'new frozen VTU':audit['frozen_flow_sha256'],'old H0 frozen VTU':c['old_native_flow_sha256'],
        'volume mesh':manifest['mesh_sha256'],'solver binary':ex['solver_sha256'],
        'PETSc library':ex['PETSc_library_sha256'],'case solver XML':manifest['solver_XML_sha256'],
        'solver.log':sha256(case/'run/solver.log'),'parameter_fit.py unchanged':sha256(ROOT/'network_1d0d/parameter_fit.py'),
        'feasibility certificate unchanged':sha256(cert)}
    hash_table=table(['文件/对象','SHA256'],[[k,f'`{v}`'] for k,v in evidence_hashes.items()])
    report_text=f'''# 有限设计范围内最佳可实现三出口分流：0D 设计与唯一 3D 前向验证

## 结论

固定完整 A 模型与 `[0.5,2.0]²` 合成设计范围内，得到 `s_O1={o['parameters']['s_O1']:.15g}`、`s_O2={o['parameters']['s_O2']:.15g}`。0D 主目标 J 比 `(1,1)` 改善 {100*d['improvement']['relative_J_reduction']:.8f}%，O1 上界激活。唯一新 CFD 已实际运行且通过独立日志、原稳态门槛和最终场检查；状态为 **{c['status']}**，3D J 相对实际旧 H0 改善 {100*c['relative_J_reduction']:.8f}%。没有 CFD 回调优化，没有第二次 CFD，没有 Particle/RBC 运行或生产晋升。

精确等分：**NO / NOT REQUIRED**。旧不可达性证书及 `BALANCED_0D_TARGET_NOT_REACHED` 原样保留。本轮结果是声明范围内经网格搜索与双初值检查一致的设计，不是生理最优、实验参数辨识或网格无关性证明。

## 1. 版本、方法与科学约束

- 实际 starting commit：`{pre['starting_commit']}`；历史基线 `f43c09b1e55cbd9702a65c600a464154fcf346f3` 未回退。
- 工作树：`{ROOT.parent.resolve()}`；分支 `sync/current-h0-dt1ms-wss-audit-20260927`。
- 完整图：7422 节点、7421 边，原始几何与端口映射哈希已封存。
- 所有优化 forward 调用原 `parameterized_hydraulics.solve_parameterized_operating_point()`；完整图解析阻力 → `G P=b` → ROI 提取。无 ROI-only surrogate、无独立端口 R/P 调参。
- `workflow_kind=DESIGN_OPTIMIZATION`，`parameter_interpretation=DESIGN_VARIABLES`，`data_identification_claim=NONE`。旧辨识保护未修改。
- `mu=0.00345312 Pa·s`，`Pd=0 Pa`，`O3=legacy_reference`，附加 O3 terminal resistance 为 `None`，ROI/J2→O1/J2→O3 几何固定。
- 0D ROI Qin `{d['roi_target_flow_m3_s']:.16g} m³/s`；production FEM Qin `{pol['Q_target_m3_s']:.16g} m³/s` 原样保留，不把网络流量的末位浮点差异写回入口。
- bounds 来源 `SYNTHETIC_DESIGN_ENVELOPE`：两项 `[0.5,2.0]`，仅为当前数值设计范围，**不是生理/测量范围或置信区间**。
- 旧证书状态 `{read(cert)['status']}`，哈希核验及对应 pytest 通过。必要约束 `Q_O1/Q_O3 <= {d['ratio_bound_Q_O1_over_Q_O3']:.17g}` 排除了 1/3 等分。更大必要区域投影 `{d['theoretical_necessary_region_projection']}` 不代表有限参数可实现解。

## 2. 纯 0D 优化与基线

目标 `J=sum((f_i-1/3)^2)`；守恒时等于 `3*population_variance(f)`。优化器 `scipy.optimize.least_squares`，log 坐标、原生有限 bounds、三点 Jacobian，`ftol=xtol=gtol=1e-12`，没有参数 clipping。NumPy {d['provenance']['numpy_version']}；SciPy {d['provenance']['scipy_version']}。

{ftab}

主初值 `(1,1)`；优化器 `{optimizer['message'].replace(chr(96),'')}`，nfev={optimizer['nfev']}、njev={optimizer['njev']}、optimality={optimizer['optimality']:.8g}。实际主优化 forward 评价 {o['forward_evaluations']} 次，primary 全流程含网格/第二初值/重复共 {d['primary_forward_evaluations']} 次。`nfev` 不含全部有限差分调用，因此不同于 forward 总次数。

`BOUND_ACTIVE_OPTIMUM`：`s_O1: UPPER`，`s_O2: NONE`。在有限范围内落在上界是允许的结果，结果依赖设计范围。J 改善 {d['improvement']['delta_J']:.15g}，极差降低 {d['improvement']['range_reduction']:.15g}，std 降低 {d['improvement']['std_reduction']:.15g}。

全网质量审计：最大相对残差 {o['prediction']['mass_audit']['max_relative_residual']:.8g}，ROI 相对残差 {o['prediction']['mass_audit']['roi_relative_residual']:.8g}；无端口倒流，所有有效半径、边阻力有限且为正。

25×25 对数网格共 {cross['grid_points']} 点，网格最好 J={cross['grid_best']['J_balance']:.15g}；从其启动第二次 bounded 局部优化。两个局部结果参数绝对差 `{cross['parameter_absolute_difference']}`，J 差 {cross['objective_absolute_difference']:.8g}，交叉检查 **{cross['status']}**。重复主优化的 {d['deterministic_reproducibility']['compared_forward_evaluations']} 次评价载荷逐项完全一致。此检查不等同于全局最优数学证明。

### 纯 0D 上界敏感性（仅表格）

下界始终 0.5。1.25、1.5、3 的结果仅为 0D 诊断；唯一 CFD 始终使用预先声明的上界 2，未根据结果更换 envelope。

{btab}

冻结文件：`{(design/'frozen_balance_design.yaml').resolve()}`。状态 `FROZEN_DESIGN`，SHA256 `{evidence_hashes['frozen_balance_design.yaml']}`。设计在 FEM handoff 与唯一 CFD 启动前已封印；后续没有改写。

## 3. 单向 FEM handoff

从冻结设计取得 real-cut P/Q；只读取当前延伸段几何，用显式 `mu=0.00345312 Pa·s` 和现有渐变半径解析阻力积分。沿每段使用原 H0 的 40 个截面中点取样及端部半格常半径封闭规则。

```text
R_ext = (8*mu/pi) integral(ds/r(s)^4)
P_cap_raw[i] = P_realcut[i] - R_ext[i]*Q[i]
P_cap_applied[i] = P_cap_raw[i] - min(P_cap_raw)
```

{htab}

共同减去 `{h['gauge_subtracted_constant_Pa']:.15g} Pa`（即加上 `{h['gauge_added_constant_Pa']:.15g} Pa`）；平移前后两两压差最大误差 `{h['max_pressure_difference_preservation_error_Pa']:.8g} Pa`。raw cap 负值是参考零点下的表压，不裁剪。对每个出口使用相同常数平移，保持驱动压差。

纯解析 H0 handoff 回归已运行：R_ext、raw cap P、applied cap P 的最大绝对差均为 **0**；没有为此运行 CFD。详见 `handoff_H0_regression.json`。

## 4. 唯一 CFD：配置、GPU 与运行证据

本地算例：`{case.resolve()}`。

服务器：`root@50.115.148.16:4159`；算例 `/workspace/flow_best_feasible_balance_20260928/{case.name}`；private supervisor `best_balance_gpu`，autorestart=false。`reports/one_scientific_case_dispatch.json` 排他登记启动，完成记录为 `reports/dispatch_completion.json`。

**NEW SCIENTIFIC CFD RUN COUNT = 1**。其他 envelope 的 CFD、校准重试 CFD、粒子与 RBC 调用：均 0；**CFD-based retuning = NO**。

预检通过后才启动：H0 XML 仅许可 3 个出口 traction `Value` 字段变化；实际 O1 从 585.2864332598318 到 {pre['applied_cap_pressure_pa']['O1']:.15g} Pa，O2 从 2932.0152710782013 到 {pre['applied_cap_pressure_pa']['O2']:.15g} Pa，O3 两者均 0。其余 XML 字段、体网格/外表面/壁面、入口、policy、PETSc options 均逐文件/逐字段核验不变。

求解器为原 `svmultiphysics`，P1/P1+VMS，刚性不可压缩牛顿流体，`rho=1056 kg/m³`，无滑移壁面；网格 {manifest['node_count']} 节点、{manifest['tetra_count']} 四面体。RTX 4090，MPI=1、OMP=1；沿用原 PETSc CUDA 矩阵和原预条件配置。没有为本任务改变并行数、后端或性能配置。

实际 dt 为 `{pol['dt_s']:.17g} s`，取自固定原算例 policy（分支名中的 dt1ms **不代表本算例步长**）。每 {pol['save_interval_steps']} 步采样，连续 {pol['steady_last_intervals']} 个间隔满足速度变化 ≤{pol['velocity_change_limit']:g}、流量变化 ≤{pol['flow_change_limit']:g}、全局质量残差 ≤{pol['mass_limit']:g}，线性/非线性原门槛不变。从零场开始，无重启或缩放旧场。

实际 exit code `{ex['exit_code']}`，wall time `{ex['wall_time_s']:.6f} s`（{ex['wall_time_s']/60:.3f} min），停止器合格 step `{ex['stop']['step']}`，最终原生输出 step `{ex['final_step']}`，最终物理时间 `{ex['final_step']*pol['dt_s']:.12g} s`。两者可能因写 STOP_SIM 的安全退出多一个时间步而不同。

GPU 使用证据：solver.log 的矩阵类型 `seqaijcusparse`；5 秒采样最大 GPU 利用率 {max(gpu_samples):g}%，显存峰值 {highmem:g} MiB。实际 PETSc profile 摘录：

```text
{chr(10).join(gpu_lines)}
```

这是 CPU/GPU 混合实现：FEM 装配及 ILU 因子构建仍有 CPU 工作；profile 中 MatMult、PCApply 的 GPU 浮点运算占比为 100%，MatLUFactorNum 和 PCSetUpOnBlocks 为 0%。后端保持原 GMRES(100)、右预条件 ASM overlap=2、子域 ILU(2)、CUDA matrix/vector 配置。不能把显存占用或峰值利用率解读为整个 CFD 都在 GPU 上执行。完整执行 command、host、软件/库哈希及采样保存在算例 `reports/` 和 `run/petsc_profile.txt`。

## 5. 独立验收与实际 3D 分流

本地重新解析原生 solver.log，与远端历史逐项一致，线性/非线性门槛通过。重新读取所有保存的稳态快照，独立重算间隔，合格 step `{audit['independently_qualifying_steps']}`。最终网格坐标/拓扑与输入一致，checkpoint step/time 验证通过。

- 实际 Qin `{m['Q_in_m3_s']:.16g} m³/s`；入口相对目标误差 `{m['epsilon_Q']:.8g}`。
- 全局质量残差 `{m['epsilon_mass']:.8g}`，补偿求和残差 `{m['epsilon_mass_compensated']:.8g}`；各出口积分流量均为正，无净倒流，未据此声称每个面片都不存在局部回流。
- 最大速度 `{m['velocity_max_m_s']:.12g} m/s`，压力范围 `{m['pressure_range_pa']}` Pa；全部有限。
- 壁面最大速度 `{m['wall_velocity_max_m_s']:.8g} m/s`，无滑移通过。
- 最终场相对停止合格场的速度 L2 变化 `{audit['final_to_stop_relative_velocity_change']:.8g}`、归一化端口流量变化 `{audit['final_to_stop_normalized_flow_change']:.8g}`，均通过原门槛。

{dtab}

分流采用带符号 Qout/实际 Qin，没有 clipping 或强行归一化。表中 cap 面积平均流体压力由解场积分，区别于施加的法向牵引压力；两者不必逐点相同。3D–0D 差异如实保留，既不是调参信号，也不作为必须为零的验收门槛。

实际旧 H0 3D 分流由其正式 frozen VTU 重新积分，未硬编码四舍五入数字。比较状态 **{c['status']}**：仅在新 J 和极差均低于旧值时才标记改善。J 相对降低 {100*c['relative_J_reduction']:.8f}%，极差减少 {c['range_reduction']:.12g}，std 减少 {c['std_reduction']:.12g}（正值为改善，负值为恶化）。旧/新完整指标见第 2 节和 `all_balance_metrics.csv`。

新 frozen flow：`{audit['frozen_flow']}`。其 `manifest.json` 包含设计 envelope、0D 参数/指标、cap pressures、实际 3D 指标、0D–3D 差异、旧 H0 比较、输入/输出哈希及最终 step；`particle_production_promoted=false`。

## 6. 已有 WSS 链路复用

复用并归档已验证的 `flow_solver_support/wss.py` 与 `wss_case.py` 原字节源码，哈希见 `evidence/wss_core/provenance.json`，旧解析张量及回归证据同时保留。没有重跑圆管 CFD。本轮对旧/新实际 FEM 场重新运行该链路：从壁面相邻四面体的 P1 完整速度梯度张量，计算单位外法向的切向黏性牵引模长：

```text
tau = mu * (grad(u) + grad(u).T) @ n
WSS = norm(tau - dot(tau,n)*n)
```

只处理壁面 tag 1，排除入口/出口截面；原坐标 m、速度 m/s、黏度 Pa·s 得到 Pa。统计使用原始面片量，面积加权均值与经验累积面积分位数。J1/J2 区域按面片中心落入半径 5 μm 的球选择，球心分别 `(92,49,111)` 和 `(130.04,82.04,87.18)` μm。全壁面统计含人工延伸段，不能称其为仅真实血管统计。

{wtab}

`figures/wss_raw_comparison.png/pdf` 使用相同正交视角、相同网格、原始面片值与线性色标 0–55 Pa，无平滑。VTP 另外保留原链路的面积加权节点显示值，但本图和 CSV 不使用它。0D 网格目标与真实 3D 分流比较见 `figures/balance_design_and_forward.png/pdf`；bound sensitivity 仅表格。

## 7. 测试、执行与复核

启动前实际执行命令（从仓库根目录）：

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=ulm_3D_vascular \\
/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv/bin/python -B -m pytest \\
ulm_3D_vascular/tests/network_parameterized \\
ulm_3D_vascular/tests/network_1d0d \\
ulm_3D_vascular/tests/network_h0 -q -p no:cacheprovider
```

结果 **{testline}**。覆盖旧完整测试与新增的不可达性保留、辨识 guard 不变、原生 bounds、无限半径拒绝、设计改善、方差关系、trace 确定性、双初值一致、bound-active、冻结 roundtrip/篡改/覆盖拒绝、零 CFD 依赖、H0 handoff 回归、XML 受控改变及唯一启动保护。早期一次新增测试发现 NumPy scalar YAML 序列化问题，已在正式冻结前修复；失败日志保留，不冒充首轮全部通过。

完整输入和下述只读后处理入口随提交保存；独立验收的 frozen 输出禁止覆盖，复核时应在独立副本中执行。**不要重启 scientific CFD**；唯一 dispatch 文件存在时启动器会拒绝重复执行。

```text
scripts/optimize_balanced_0d_design.py --config configs/parameterized_hydraulics/best_feasible_balance_design.yaml --output <new-audit-dir> --freeze
scripts/prepare_best_balance_fem.py --help
scripts/validate_best_balance_fem.py --case <case-copy> --source-case <old-H0> --flow-root <FEM-root> --report <report-copy>
scripts/postprocess_best_balance.py --case <case> --source-case <old-H0> --flow-root <FEM-root> --design <design-dir> --report <report-dir>
scripts/report_best_balance_forward.py --case <case> --source-case <old-H0> --design <design-dir> --report <report-dir>
```

数值证据：`design_trace.csv`、`coarse_grid.csv`、`bound_sensitivity.csv`、`all_balance_metrics.csv`、`final_3D_ports.csv`、`wss_comparison.csv`；独立间隔检查见 `independent_validation.json`。原生中间快照留在本地和服务器 case 的 `run/1-procs/`；Git 保存正式最终场、checkpoint、完整求解日志和验收摘要，不重复提交所有中间大文件。

## 8. 关键哈希

{hash_table}

## 9. 尚未解决的科学限制

1. 这是给定完整 A、O3 参考闭合和合成 bounds 下的设计，不是生理真值或实验 calibration；O1 上界激活表明结果依赖 envelope。
2. 25×25 网格和双初值交叉检查增加数值可信度，不构成严格全局最优证明；已有不可达性证书的必要区域投影也不是已实现的有限设计。
3. 0D 圆截面/充分发展阻力近似、延伸段阻力补偿和 3D 分叉几何不同，因此允许分流偏差；本轮没有耦合迭代或 CFD 反馈校正。
4. 保持原单张网格、P1/P1+VMS 和相同时间策略，只完成受控边界改变的 forward 比较；没有新的网格独立性、内部截面质量误差研究、时间步或模型不确定性量化。全局质量守恒不能证明局部梯度/WSS 准确。
5. WSS 面片差异受 P1 梯度、曲面离散和原网格分辨率影响；本轮不据图像平滑程度诊断物理正确性，也不为颜色均匀调压或平滑。
6. 本轮未与实验流量/压力/WSS 比较，未运行 CORE500、microbubble 或 RBC，未把新流场晋升为轨迹生产输入。
'''
    (report/'BEST_FEASIBLE_BALANCE_FORWARD_REPORT_ZH.md').write_text(report_text)
    registry=read(report/'scientific_case_registry.json')
    registry.update(status='COMPLETED_AND_INDEPENDENTLY_VALIDATED',solver_launches=1,new_scientific_CFD_run_count=1,
        comparison_status=c['status'],frozen_flow=audit['frozen_flow'])
    write_json(report/'scientific_case_registry.json',registry)
    write_json(report/'delivery_status.json',dict(status='PASS',design_status=d['status'],CFD_status=audit['status'],
        balance_comparison=c['status'],new_scientific_CFD_run_count=1,CFD_based_retuning=False,
        CORE500_calls=0,microbubble_calls=0,RBC_calls=0,particle_production_promoted=False,
        preflight_tests=testline,protected_artifact_hashes_verified=True))
    print(report/'BEST_FEASIBLE_BALANCE_FORWARD_REPORT_ZH.md')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('case','source-case','design','report'):parser.add_argument('--'+key,type=Path,required=True)
    a=parser.parse_args();run(a.case,a.source_case,a.design,a.report)
