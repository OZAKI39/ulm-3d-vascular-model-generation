"""Save small, pure-0D validation evidence and its Chinese review report."""
from pathlib import Path
import argparse
import json
import sys
from dataclasses import asdict
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from network_1d0d.audit import write_json,write_csv,sha256
from network_1d0d.parameter_config import read_configuration
from network_1d0d.parameterized_hydraulics import load_geometry_cache,ParameterizedHydraulicSpec as Spec,solve_parameterized_operating_point as solve
from network_1d0d.parameter_fit import Targets,Observation,fit_parameters,UnidentifiableParameterizationError
from network_1d0d.parameter_freeze import load_frozen_artifact


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--report-directory',type=Path,required=True)
    args=p.parse_args();out=args.report_directory
    loaded=read_configuration(ROOT/'configs/parameterized_hydraulics/legacy_h0_regression.yaml')
    cache=load_geometry_cache(loaded['geometry_path'],loaded['ports_path'])
    radius=[]
    for factor in [.95,1.,1.05]:
        r=solve(cache,Spec(s_O2=factor))
        radius.append(dict(s_O2_dimensionless=factor,
            **{p+'_fraction_dimensionless':float(f) for p,f in zip(('O1','O2','O3'),r.signed_outlet_fractions)},
            mass_residual_dimensionless=r.mass_audit['max_relative_residual'],
            minimum_radius_m=float(r.effective_radius_m.min()),minimum_resistance_pa_s_m3=float(r.edge_resistance_pa_s_m3.min())))
    assert radius[0]['O2_fraction_dimensionless']<radius[1]['O2_fraction_dimensionless']<radius[2]['O2_fraction_dimensionless']
    write_csv(out/'radius_sensitivity.csv',radius)
    R=1e17;r=solve(cache,Spec(o3_mode='terminal_resistance',terminal_resistance_O3_pa_s_m3=R))
    closure=dict(terminal_resistance_O3_pa_s_m3=R,Pd_pa=0.,P_O3_pa=float(r.port_pressure_pa[3]),
                 Q_O3_m3_s=float(r.port_flow_m3_s[3]),R_times_Q_pa=float(R*r.port_flow_m3_s[3]),
                 closure_error_pa=float(r.port_pressure_pa[3]-R*r.port_flow_m3_s[3]),
                 mass_residual_dimensionless=r.mass_audit['max_relative_residual'],provenance=r.closure_provenance)
    write_json(out/'O3_closure_validation.json',closure)
    target=Targets('SYNTHETIC_TARGET',tuple(Observation('flow_fraction',p,float(f),.005) for p,f in zip(('O1','O2','O3'),r.signed_outlet_fractions)),'pure_0D_underdetermined_example')
    try:
        fit_parameters(cache,r.spec,target,free_parameters=('s_O1','s_O2','terminal_resistance_O3_pa_s_m3'))
    except UnidentifiableParameterizationError as e:
        rejection=dict(status='REJECTED_AS_REQUIRED',audit=e.audit,forward_evaluations_before_rejection=0)
    else:raise AssertionError('Underdetermined fit unexpectedly accepted')
    write_json(out/'identifiability_rejection.json',rejection)
    legacy=json.loads((out/'legacy_regression/fit_summary.json').read_text())
    fit=json.loads((out/'synthetic_radius_fit/fit_summary.json').read_text())
    truth=json.loads((out/'synthetic_radius_fit/synthetic_truth.json').read_text())
    frozen=load_frozen_artifact(out/'synthetic_radius_fit/frozen_parameterized_boundary.yaml')
    for k in ('s_O1','s_O2'):assert abs(fit['final_parameters'][k]-truth['parameters'][k])<1e-7
    summary=dict(legacy_regression=legacy['regression'],legacy_prediction=legacy['prediction'],
        synthetic_recovery=dict(true_parameters=truth['parameters'],fitted_parameters=fit['final_parameters'],
            initial_residual=fit['initial_residual'],final_residual=fit['final_residual'],
            initial_norm=fit['initial_objective_norm'],final_norm=fit['final_objective_norm'],nfev=fit['nfev'],
            trace_evaluations=fit['total_forward_evaluations'],status=fit['status']),
        radius_sensitivity=radius,O3_closure=closure,underdetermined_guard=rejection,
        CFD_calls=0,particle_simulation_calls=0,particle_based_fitting=False,
        frozen_artifact=str((out/'synthetic_radius_fit/frozen_parameterized_boundary.yaml').resolve()),
        frozen_content_sha256=frozen['artifact_sha256'],evidence_script_sha256=sha256(__file__))
    write_json(out/'VALIDATION_SUMMARY.json',summary)
    table='\n'.join(f"| {row['s_O2_dimensionless']:.2f} | {row['O2_fraction_dimensionless']:.10f} | {row['mass_residual_dimensionless']:.3e} |" for row in radius)
    old=legacy['prediction'];ports='\n'.join(f"| {p} | {old['real_cut_pressure_pa'][p]:.11f} | {old['port_flow_m3_s'][p]:.12e} |" for p in ('INLET','O1','O2','O3'))
    tests=(out/'logs/pytest_all.log').read_text().strip().splitlines()[-1]
    report=f'''# 参数化完整 0D 边界模型：实现与数值验收

目标仓库 `OZAKI39/ulm-3d-vascular-model-generation`，分支 `sync/current-h0-dt1ms-wss-audit-20260927`。开发起点 `52fcb419dadbcca7566f263768e7889c2d33f1db`；工作树为 `/home/lzy/projects/temp_storage/github_sync_20260927`。本次只开发纯 0D，不更改当前 H0/FEM/Particle/RBC 工作流。

## 交付代码

新增 `network_1d0d/parameterized_hydraulics.py`（SI 几何缓存、掩码、全网前向、O3 closure）、`parameter_fit.py`（加权残差、对数参数、秩检查、trace）、`parameter_config.py`（显式单位 YAML、命令工作流）、`parameter_freeze.py`（来源校验、不可覆盖及内容封印）。

新增脚本 `scripts/fit_a_parameterized_0d.py`、`freeze_a_parameterized_0d.py` 和本报告生成器 `report_parameterized_0d_validation.py`；新增 `configs/parameterized_hydraulics/` 两份 YAML 和带来源 SHA 的 H0 回归后备快照；新增 `tests/network_parameterized/`；新增方法文档 `docs/PARAMETERIZED_0D_BOUNDARY_MODEL.md` 及本报告目录。

唯一修改的既有文件为 `utils/cfd_preprocess/one_d_flow.py`：只把标量 `edge_resistance()` 委托到权威解析实现，并保留正阻力检查。旧 `solve_global_flow()` 未改。原实现对一般 taper 在浮点容差内等价，但 r1/r0=1+1e-8 时相对差约 −5.53e-10；统一公式避免近等半径消差及未来漂移。详见 `legacy_resistance_before_unification.json`。

未修改 `idealized_h0.py`、`hydraulic_resistance.py`、`network_solver.py`、`boundary_conditions.py`、`extension_transfer.py`、现有 H0 reports，亦未修改 `formal_3D_flow_solver/`、`ulm_flow_mean_2p0_mmps/`、`ulm_particle_formal_p9a5/`（含 CORE500 与 RBC）、`vascular_printing/`、`server_evidence/`。

## 数学定义、缓存和每次 evaluation

优化变量是 x1=log(s_O1)、x2=log(s_O2)，可选 x3=log(R_O3/R_ref)；默认 μ=0.00345312 Pa·s 固定。O1/O2 是保存的下游几何节点半径尺度；O3 R 是 `MODEL_CLOSURE_PARAMETER`，不是虚构的几何树。P、Q 均为求解状态。

预存 {len(cache.node_ids)} 节点、{len(cache.edge_nodes)} 边的长度、原始半径/阻力、ROI/下游掩码、源/终端/端口和物理符号。O1/O2 分别缩放 {int(cache.downstream_node_masks['O1'].sum())}/{int(cache.downstream_node_masks['O2'].sum())} 个下游节点；端口 real-cut 半径保持不变。

每次 evaluation：复制有效半径 → 缩放相应节点 → 用 `linear_radius_resistance()` 重算全部几何边的解析 taper R → 如有需要增加一个 O3 虚拟阻力边和 distal 参考节点 → 原 `solve_network()` 组装并联立求解整个 G P=b → 正向匹配 ROI Q → 加回共同 Pd → 提取有符号 P/Q、质量和反流审计 → 构造残差向量。没有使用旧 H0 solver 的返回值替代新路径，也没有对第一条连接切口的边粗暴乘 s^-4。

所要求旧文件已逐一读取，其复用关系和限制逐项见方法文档“文件复用与隔离”。特别是 `load_h0()` 会检查历史绝对路径（含 FEM 文件）清单，因此新优化器以几何-only loader 加载保存的 H0 SI 图；原 loader 和原 H0 solver 仍独立通过旧测试。

## 已运行数值证据

**H0 回归**：端口压力最大差 {legacy['regression']['maximum_pressure_error_pa']:.3e} Pa；分流最大差 {legacy['regression']['maximum_fraction_error_dimensionless']:.3e}；入口 Q 差 {legacy['regression']['inlet_flow_error_m3_s']:.3e} m³/s。

| 端口 | real-cut 压力 Pa | 正方向流量 m³/s |
|---|---:|---:|
{ports}

O1/O2/O3 分流依次 {old['outlet_flow_fraction']['O1']:.10f} / {old['outlet_flow_fraction']['O2']:.10f} / {old['outlet_flow_fraction']['O3']:.10f}。最大相对质量残差 {old['mass_audit']['max_relative_residual']:.3e}；全网 {old['mass_audit']['global_relative_residual']:.3e}；ROI {old['mass_audit']['roi_relative_residual']:.3e}。这些是 H0 regression，不是实验真值。

质量残差同时与保存金标准核对：全网相对残差差值 {old['mass_audit']['global_relative_residual']-legacy['regression']['legacy_mass_audit']['global_relative_residual']:.3e}，来自等价求和顺序的浮点舍入；未人为强制归零。其余主要质量指标在记录精度内一致，均远低于 1e-9 门槛。

**合成恢复**：真值 (0.92, 1.07)，初值 (1, 1)，仅拟合三出口比例，sigma_f=0.005。恢复 ({fit['final_parameters']['s_O1']:.15f}, {fit['final_parameters']['s_O2']:.15f})。标准化残差范数 {fit['initial_objective_norm']:.12g} → {fit['final_objective_norm']:.6e}，SciPy nfev={fit['nfev']}；全部 trace 前向评估 {fit['total_forward_evaluations']} 次（含数值 Jacobian）。状态 {fit['status']}；初末数据 Jacobian 秩均为 2。真值、每次 trace、完整节点/边状态和收敛审计均已保存。

初始 residual = `{fit['initial_residual']}`；最终 residual = `{fit['final_residual']}`。

**半径敏感性**：

| O2 radius scale | O2 signed fraction | 最大相对质量残差 |
|---:|---:|---:|
{table}

方向严格递增，所有半径和阻力 finite positive。O1/O2 掩码与旧 `external_radius_variant()` 一致，切口及 ROI 节点未动；第一 taper edge 的阻力不等于原值统一乘 s^-4。

**O3 closure**：R={R:.1e} Pa·s/m³，Pd=0 Pa；P_O3={closure['P_O3_pa']:.12g} Pa，Q_O3={closure['Q_O3_m3_s']:.12e} m³/s；R·Q={closure['R_times_Q_pa']:.12g} Pa；闭合误差 {closure['closure_error_pa']:.3e} Pa。O3 不再固定，虚拟 distal 节点才固定；最大相对质量残差 {closure['mass_residual_dimensionless']:.3e}。该算例仅为纯 0D 一致性检查。

**可辨识性**：3 个自由参数、2 个独立 flow constraints、0 pressure、0 prior 的配置在任何 fit forward 前拒绝；保存 `identifiability_rejection.json`。测试还覆盖 Q 与 fraction 不重复计数、固定参考压力不算观测、计数足够但数值 Jacobian 秩不足仍拒绝、加压力可辨识三参数、prior 仅称正则化。legacy 模式下不存在的 O3 R 明确为 inactive，不虚称 anchor。

## 测试、冻结与审计

实际命令（工作目录为仓库根目录）：

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=ulm_3D_vascular \\
 /home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv/bin/python -B -m pytest \\
 ulm_3D_vascular/tests/network_parameterized ulm_3D_vascular/tests/network_1d0d \\
 ulm_3D_vascular/tests/network_h0 -q -p no:cacheprovider
```

结果：**{tests}**。包含原 H0 独立求解/leave-one-terminal-out、解析阻力、端口方向反转、SI/质量、合成恢复、trace 字节级确定性、冻结 roundtrip/篡改拒绝/禁止覆盖、非法目标来源、依赖静态检查及实际拟合时禁止进程启动/CFD 场读取的运行时检查。旧 FEM 相关测试仅检查已有文件及 XML，不创建/启动算例。没有运行生产 CFD 或粒子测试。

最终冻结文件：`{summary['frozen_artifact']}`。

状态 FROZEN，canonical 内容 SHA256 `{frozen['artifact_sha256']}`。冻结前只在固定最终参数上做一次纯 0D 校验；只读文件、exclusive-create、防覆盖、来源及内容哈希校验均已测试。`development_runs/` 保留第一次审阅前的输出及未改动冻结件，仅 inactive/anchor 审计标签随后修正，正式文件使用根部 `synthetic_radius_fit/`。

**任何 CFD 调用：NO。任何 particle-based parameter fitting：NO。粒子模拟调用次数：0。**既有 3D 报告只用于理解历史，任何 3D Q/P/WSS 未进入辨识目标或参数更新。

optional FEM pressure handoff 本版未实现：现有扩展段 helper 使用默认 μ、写死 H0 provenance，直接接入新模型会增加科学语义风险。本轮不要求该可选项，也不创建 FEM case。

## 结论边界

已验证的是完整稳态 Newtonian 电阻网络实现、自身合成恢复、回归与审计机制。没有真实实验校准、全局唯一性证明或 CFD 验证。无 compliance、inertance、RCR 或 Womersley impedance；只支持共同 distal gauge。测量协方差、半径/μ/闭合参数相关性、真实边界与 SWC 几何误差仍是模型限制。后续 3D 只能读取冻结参数作独立 forward prediction，禁止据其差异 retune。
'''
    (out/'IMPLEMENTATION_REPORT_ZH.md').write_text(report)
    print(json.dumps(summary['synthetic_recovery'],indent=2))

if __name__=='__main__':main()
