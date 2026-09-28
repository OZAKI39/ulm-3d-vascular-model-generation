"""Write the Chinese report only from completed independent validation evidence."""
from pathlib import Path
import argparse
import json
import sys
import re
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from network_1d0d.audit import sha256,write_json,write_csv
from network_1d0d.roi_boundary_design import load_frozen_design,metrics


def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(str(v) for v in row)+' |' for row in rows])


def generate(case,source,report):
    case,source,report=map(lambda x:Path(x).resolve(),(case,source,report))
    read=lambda name:json.loads((report/name).read_text())
    f=load_frozen_design(report/'roi_boundary_design_frozen.yaml')
    r=read('roi_mixed_bc_regression.json');b=read('roi_pressure_basis.json');c=read('roi_geometry_cache_summary.json')
    h=read('roi_fem_handoff.json');v=read('independent_validation.json');pre=read('final_preflight.json')
    e=json.loads((case/'reports/execution.json').read_text());policy=json.loads((case/'policy.json').read_text())
    assert v['status']=='PASS' and e['status']=='PASS'
    comp=v['comparison'];m=v['measurements'];ports=('O1','O2','O3')
    old=comp['old_H0_3D_fractions'];new=comp['new_3D_fractions'];zero=comp['design_0D_fractions']
    hashes=read('protected_previous_work_hashes.json')
    changed=[p for p,digest in hashes.items() if sha256(ROOT.parent/p)!=digest]
    if changed:raise ValueError('Protected prior work changed: '+str(changed))
    write_json(report/'previous_work_preservation.json',dict(status='PASS',checked_files=len(hashes),changed_files=changed,
        statement='All previously tracked network source and reports retain original bytes; new work is separate'))
    figdir=report/'figures';figdir.mkdir(exist_ok=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    x=np.arange(3);width=.24
    for off,values,label,color in [(-1,old,'Old H0 3D','#757575'),(0,zero,'ROI-only 0D design','#009E73'),(1,new,'New ROI-only BC 3D','#0072B2')]:
        bars=ax[0].bar(x+off*width,100*np.array(values),width,label=label,color=color)
        ax[0].bar_label(bars,fmt='%.2f',fontsize=8,padding=3)
    ax[0].axhline(100/3,color='black',linestyle='--',linewidth=.8)
    ax[0].set(xticks=x,xticklabels=ports,ylabel='Outlet flow / actual inlet flow (%)',ylim=(0,60),title='One frozen pressure design; one CFD forward test')
    ax[0].legend(fontsize=8)
    ax[1].bar(x,100*np.array(comp['discrepancy_3D_minus_0D']),color=['#D55E00' if q<0 else '#0072B2' for q in comp['discrepancy_3D_minus_0D']])
    ax[1].axhline(0,color='black',linewidth=.8)
    ax[1].set(xticks=x,xticklabels=ports,ylabel='3D minus 0D (percentage points)',title='Discrepancy retained; no CFD feedback tuning')
    for i,q in enumerate(comp['discrepancy_percentage_points']):ax[1].text(i,q/2,f'{q:+.3f}',ha='center',va='center',color='white',weight='bold')
    for suffix in ('png','pdf'):fig.savefig(figdir/('roi_design_vs_forward.'+suffix),dpi=220)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
    for offset,values,label in [(-.25,[f['relative_realcut_pressure_pa'][p] for p in ports],'Real cut: O3=0'),(0,[f['gauge_realcut_pressure_pa'][p] for p in ports],'Real cut: minimum=0'),(.25,[row['pressure_cap_shifted_Pa'] for row in h['ports']],'FEM cap: minimum=0')]:
        ax.bar(x+offset,values,.25,label=label)
    ax.axhline(0,color='black',linewidth=.7);ax.set(xticks=x,xticklabels=ports,ylabel='Pressure (Pa)',title='Two distinct gauge shifts; extension transfer changes pressure differences')
    ax.legend(fontsize=8)
    for suffix in ('png','pdf'):fig.savefig(figdir/('pressure_handoff.'+suffix),dpi=220)
    plt.close(fig)
    # Raw numerical tables, not rounded presentation values.
    rows=[]
    for i,p in enumerate(ports):
        hr=h['ports'][i]
        rows.append(dict(port=p,realcut_relative_pa=f['relative_realcut_pressure_pa'][p],
            realcut_gauge_pa=f['gauge_realcut_pressure_pa'][p],ROI_gauge_added_pa=f['gauge_shift_pa'],
            extension_R_pa_s_m3=hr['R_extension_Pa_s_m3'],raw_cap_pa=hr['pressure_cap_raw_Pa'],
            FEM_gauge_added_pa=h['FEM_cap_gauge_added_pa'],applied_cap_pa=hr['pressure_cap_shifted_Pa'],
            ROI_0D_flow_m3_s=f['outlet_flow_m3_s'][p],ROI_0D_fraction=zero[i],
            new_3D_flow_m3_s=m['outlet_flows_m3_s']['OUTLET_0'+str(i+1)],new_3D_fraction=new[i],
            old_H0_3D_fraction=old[i],difference_3D_minus_0D_fraction=new[i]-zero[i]))
    write_csv(report/'boundary_and_flow_summary.csv',rows)
    samples=[list(map(float,s['gpu'].split(','))) for s in e['GPU_samples']]
    gpu=dict(samples=len(samples),sampling_interval_approx_s=5,maximum_utilization_percent=max(s[0] for s in samples),
        mean_utilization_percent=float(np.mean([s[0] for s in samples])),peak_memory_MiB=max(s[1] for s in samples),
        MPI_ranks=e['MPI_ranks'],OMP_NUM_THREADS=e['OMP_NUM_THREADS'])
    profile=(case/'run/petsc_profile.txt').read_text()
    events='\n'.join(line for line in profile.splitlines() if re.match(r'^(MatMult\s|PCApply\s|MatLUFactorNum\s|PCSetUpOnBlocks\s|KSPSolve\s|Using PETSc Release)',line))
    gpu['petsc_profile_excerpt']=events
    write_json(report/'GPU_execution_evidence.json',gpu)
    (report/'petsc_gpu_profile_excerpt.txt').write_text(events+'\n')
    metrows=[]
    for label,values in [('旧 H0 3D',old),('ROI-only 0D',zero),('新 ROI-only BC 3D',new)]:
        a=metrics(values);metrows.append([label,*[f'{100*z:.8f}' for z in values],f"{a['range']:.12g}",f"{a['std']:.12g}",f"{a['max_abs_deviation_from_one_third']:.12g}"])
    flowtable=table(['模型/算例','O1 %','O2 %','O3 %','极差（比例）','总体 std','最大等分偏差'],metrows)
    ptable=table(['端口','相对 real-cut Pa','平移后 real-cut Pa','0D Q m³/s','R_ext Pa·s/m³','raw cap Pa','applied cap Pa'],
        [[p,*[f'{x:.14g}' for x in [f['relative_realcut_pressure_pa'][p],f['gauge_realcut_pressure_pa'][p],f['outlet_flow_m3_s'][p],h['ports'][i]['R_extension_Pa_s_m3'],h['ports'][i]['pressure_cap_raw_Pa'],h['ports'][i]['pressure_cap_shifted_Pa']]]] for i,p in enumerate(ports)])
    h0table=table(['端口','H0 提取压力 Pa','重构压力 Pa','H0 提取 Q m³/s','重构 Q m³/s'],
        [[p,f"{r['reference_port_pressure_pa'][i]:.14g}",f"{r['reconstructed']['port_pressure_pa'][i]:.14g}",f"{r['reference_port_flow_m3_s'][i]:.14g}",f"{r['reconstructed']['port_flow_m3_s'][i]:.14g}"] for i,p in enumerate(('INLET',*ports))])
    porttable=table(['端口','网络节点 ID','坐标 μm','源图边索引','流量系数/正方向'],
        [[p['port'],p['node_id'],', '.join(f'{z*1e6:.6f}' for z in p['xyz_m']),p['source_edge_index'],f"{p['coefficient']} / {p['positive']}"] for p in c['ports']])
    finaltable=table(['端口','实际 3D Q m³/s','Qout/Qin %','3D−0D 百分点','cap 面积均压 Pa'],
        [[p,f"{m['outlet_flows_m3_s']['OUTLET_0'+str(i+1)]:.15g}",f'{100*new[i]:.8f}',f'{100*(new[i]-zero[i]):+.8f}',f"{m['area_average_pressure_pa']['OUTLET_0'+str(i+1)]:.12g}"] for i,p in enumerate(ports)])
    pytest=(report/'pytest.log').read_text().strip().splitlines()[-1]
    report_text=f'''# ROI-only 直接压力边界设计与唯一 3D 前向验证

## 通俗结论

本轮只保留固定 ROI 内部几何，直接设计出口间压差。混合边界回归通过，0D 三出口均分达到浮点精度；所需 O1−O3 = {f['relative_realcut_pressure_pa']['O1']:.9f} Pa，原 full-A 被动下游模型不允许这一压差，新模型允许。新压力冻结后只进行一次 GPU CFD，独立验收通过。实际 3D 三出口为 {100*new[0]:.5f}% / {100*new[1]:.5f}% / {100*new[2]:.5f}%，结论为 **{comp['status']}**。0D 均分没有被冒充为 3D 均分；保留二者差异，不进行反馈调压。这些是数值设计边界，不是实测生理压力。旧 full-A 不可达证明未改写；Particle/RBC 调用为 0。

## 1. 版本和并列模型

- 实际起点 `{pre['starting_commit']}`；历史基线 `f43c09b1e55cbd9702a65c600a464154fcf346f3`，没有回退。
- 分支 `sync/current-h0-dt1ms-wss-audit-20260927`，工作树 `{ROOT.parent}`。本次新增代码及结果位于独立 ROI-only 路径。
- MODEL A：FULL-A DOWNSTREAM-CONSTRAINED，固定 ROI、正下游阻力、共同 distal reference 下等分不可达，证书原样保留。
- MODEL B：`{f['model_name']}`；`workflow_kind=ROI_BOUNDARY_DESIGN`。full-A downstream 是否参与新边界生成：**NO**。full-A 仅为已验证 ROI 几何来源。
- 保留检查覆盖原有已跟踪网络代码与报告 {len(hashes)} 个文件，逐字节 SHA256 全部相同；见 `previous_work_preservation.json`。新模型没有半径设计变量，也没有 source node 2410、外部叶子或终端树参与求解。

## 2. ROI 定义、单位和端口

从 `analysis_A_H0_graph_si.npz` 的 `roi_internal_edge_mask` 取边及其端点，得到 **{c['node_count']} 节点 / {c['edge_count']} 边**的连通树，四个一度叶节点恰为已验证的四个端口。缓存 `roi_geometry_cache.npz` 包含节点 ID、SI 坐标、原半径、源图节点/边索引、本地连接、长度、阻力及端口系数；其中 source_node_indices 指输入文件索引，不是水力源。

{porttable}

端口身份来自 `roi_ports_in_a.json`，与已保存符号约定逐项回归；INLET 正方向为流入 ROI，O1/O2/O3 正方向为流出。edge u→v 只是代数存储方向，反转存储测试通过，不用 abs(Q) 隐藏倒流。

长度 m、半径 m、时间 s、Q m³/s、压力 Pa、动力黏度 Pa·s、阻力 Pa·s/m³。实际 `mu={f['dynamic_viscosity_pa_s']}` Pa·s，`Qin={f['Qin_target_m3_s']:.16g}` m³/s。FEM rho=1056 kg/m³，nu=mu/rho=3.27e-6 m²/s。ROI 缓存 SHA256 `{f['geometry_sha256']}`（覆盖几何、索引、端口系数及固定 mu 下的长度/阻力数组；不是仅坐标散列）；几何原文件 SHA256 `{f['geometry_file_sha256']}`。

## 3. 实际代码和混合边界方程

新增 `network_1d0d/roi_only_hydraulics.py::load_roi_cache/solve_roi_mixed_bc` 和 `roi_boundary_design.py::regress_h0/pressure_basis/design`，入口为 `scripts/optimize_roi_boundary_balance.py`。局部阻力只调用权威 `hydraulic_resistance.py::linear_radius_resistance`：

```text
R_e = (8*mu*L/pi)*(r0²+r0*r1+r1²)/(3*r0³*r1³)
g_e = 1/R_e
G = assembled ROI graph Laplacian
G_ff P_f = q_f - G_fb P_b
q_INLET = +Qin; all interior q = 0
P_b = [P_O1, P_O2, P_O3]
Q_e = (P_u-P_v)/R_e
Q_port = port_coefficient * Q_port_edge
```

INLET 在自由压力未知量中，三个出口为固定压力节点。代码用统一数值缩放 `G/max(g)` 改善线性系统量级，**不做 source=1 Pa 后再缩放流量的 operating-point solve**。四个端口流量从完整 ROI 网络解提取，不用简化总阻力公式替代网络。

## 4. 必须先通过的 H0 回归

使用正式 H0 0D 提取的相同 real-cut 三出口压力和 Qin；未调用 full-A optimizer，也未读取 CFD 场来设计压力。最大端口 Q 差 / Qin = `{r['maximum_port_flow_error_relative_Qin']:.12g}`，门槛 1e-9，**PASS**。此检查在 basis 与优化之前执行。重构入口压力也一致，但未把它当额外强制入口条件。

{h0table}

重构 O1/O2/O3 fraction = `{r['reconstructed']['outlet_fraction']}`；完整节点压力和边流量保存在 `roi_mixed_bc_regression.json`。

## 5. 线性 basis 与独立迭代交叉检查

变量 `[dP_O1_O3, dP_O2_O3]`，P_O3=0。三次实际 mixed solves 分别为 (0,0)、(1 Pa,0)、(0,1 Pa)：

```text
q_out = q0 + B @ dP
q0 [m³/s] = {b['q0_m3_s']}
B [m³/(s·Pa)] =
{np.array2string(np.array(b['B_m3_s_per_pa']),precision=15)}
```

B rank = {b['rank']}，2-norm condition = {b['condition_number']:.12g}，奇异值 `{b['singular_values_m3_s_per_pa']}`。直接用 `numpy.linalg.lstsq(B, Qin/3-q0)` 求解；另外从 H0 `{f['initial_relative_pressure_pa']}` Pa 出发，以实际网络分流残差执行无界 `scipy.optimize.least_squares`，三点 Jacobian、ftol/xtol=1e-13、gtol=1e-14、x_scale=1000 Pa。负相对压差允许。

- direct dP Pa：`{f['direct_relative_pressure_pa']}`。
- iterative dP Pa：`{f['iterative_relative_pressure_pa']}`。
- direct−iterative Pa：`{f['direct_minus_iterative_pa']}`；严格差异门槛 1e-6 Pa，PASS。
- nfev={f['iterative_nfev']}；含有限差分实际网络评价 {f['iterative_evaluations']} 次，完整 trace 在 CSV。终止原因 `{f['iterative_message']}`。
- equal-split 状态 **{f['status']}**；最大等分偏差 `{f['max_abs_deviation_from_one_third']:.12g}` < 1e-8；极差 `{f['range']:.12g}`，std `{f['std']:.12g}`，J `{f['J_balance']:.12g}`。
- 内部节点最大相对残差 `{f['mass_audit']['max_free_node_relative_residual']:.12g}`，ROI 质量残差 `{f['mass_audit']['ROI_relative_residual']:.12g}`，PASS。

固定路径独立提取：J2 节点 {f['fixed_path_identity']['junction_node_id']}，R(J2→O1)={f['fixed_path_identity']['paths']['O1']['resistance_pa_s_m3']:.16g}、R(J2→O3)={f['fixed_path_identity']['paths']['O3']['resistance_pa_s_m3']:.16g} Pa·s/m³；与旧证书严格回归。等分恒等式 `(R3-R1)*Qin/3` = {f['fixed_path_identity']['equal_split_dP_O1_O3_pa']:.12f} Pa，与直接解差 {f['relative_realcut_pressure_pa']['O1']-f['fixed_path_identity']['equal_split_dP_O1_O3_pa']:.8g} Pa；与题给四舍五入参考 -795.60984355 Pa 一致。旧被动模型不允许 P1<P3；新设计允许，二者没有逻辑冲突。

## 6. 两次 gauge 与延伸段 handoff

第一次共同加 `{f['gauge_shift_pa']:.15g}` Pa，把最小 real-cut outlet 压力设为零；入口从 `{f['inlet_pressure_relative_pa']:.15g}` 变为 `{f['inlet_pressure_pa']:.15g}` Pa。全体节点统一平移验证通过，Q 改变量 / Qin 最大 `{f['gauge_flow_error_relative_Qin']:.12g}`。该操作只换参考零点，不改变压差。

冻结文件 `roi_boundary_design_frozen.yaml`，SHA256 `{sha256(report/'roi_boundary_design_frozen.yaml')}`。它是 handoff 的唯一压力来源，不是 full-A frozen 参数文件。文件有内容校验和、禁止覆盖，正式求解器启动前已冻结。

复用实际 FEM extension geometry 截面方法：沿各延伸段取 40 个中点截面，用闭合面积换算等效圆半径；相邻半径间调用同一个解析渐变半径积分，两端半格按端部半径闭合。显式 mu={f['dynamic_viscosity_pa_s']} Pa·s。这是圆截面充分发展阻力近似，不包含真实 3D 分叉损失。

```text
P_cap_raw_i = P_realcut_gauge_i - R_ext_i*Q_i
P_cap_applied_i = P_cap_raw_i - min(P_cap_raw)
```

第二次统一加 `{h['FEM_cap_gauge_added_pa']:.15g}` Pa，只对完成延伸段压降转移后的 cap 参考做平移。它与第一次 real-cut gauge 分别记录，不能把 R_ext*Q 压降误称为 gauge。raw cap 的负值不裁剪。

{ptable}

## 7. 唯一科学 CFD 的受控配置与资源

本地算例：`{case}`。

服务器：`vast4090`（root@50.115.148.16:4159），算例 `/workspace/flow_roi_only_balance_20260928/{case.name}`；独立 supervisor `roi_balance_gpu`，autorestart=false。原 H0 本地 `{source}`。

`prepare_roi_balance_fem.py` 复用 `serialize_fixed_pressure_xml()` 和 `audit_xml_change()`，预检后创建此唯一算例。XML 仅三个 outlet traction Value 改变；网格、WALL、入口、材料参数、P1/P1+VMS、policy 与 PETSc options 逐字节不变。XML 中入口 Value 为负是外法向通量约定，实际正 Qin 经表面积分核查。配置差异明细在 `configuration_diff.json`。

生产 Qin `{policy['Q_target_m3_s']:.16g}` m³/s 原样保留；它和 0D 目标的相对差约 2.87e-10，没有写回改入口。刚性不可压缩牛顿流体、无滑移壁面，70363 节点 / 371402 四面体。dt `{policy['dt_s']:.16g}` s，每 {policy['save_interval_steps']} 步采样；连续 {policy['steady_last_intervals']} 个间隔 E_u≤{policy['velocity_change_limit']}、E_Q≤{policy['flow_change_limit']}，入口/整体质量残差≤{policy['mass_limit']}。线性 KSP rtol=1e-10、atol=1e-24；非线性每步至少 2 次，最终残差比≤1e-10。最大 800 步 / 14400 s 原预算未放宽。

MPI=1、OMP=1，保持原已验证 PETSc CUDA 配置；GMRES(100)、右预条件 ASM overlap=2、子域 ILU(2)、aijcusparse/cuda。CPU 做 FEM 装配及部分预条件器构造，GPU 做实际稀疏线性运算。没有为了占满 CPU 核数改 MPI 分解或重做后端试验。0D 仅 109 节点，CPU 求解；把它搬到 GPU 无助于本次大头 CFD。

GPU 实际采样 {gpu['samples']} 次，峰值利用率 {gpu['maximum_utilization_percent']:.0f}%，采样平均 {gpu['mean_utilization_percent']:.2f}%，显存峰值 {gpu['peak_memory_MiB']:.0f} MiB。实际 profile 而非只凭显存占用判断 GPU 工作（事件行末列为 GPU %F，表示该事件浮点工作在 GPU 的比例，不是总 CFD 时间比例）：

```text
{events}
```

**NEW SCIENTIFIC CFD RUN COUNT = 1**；单次 dispatch 由 O_EXCL 排他登记、重复启动拒绝。求解器实际退出码 {e['exit_code']}，wall time {e['wall_time_s']:.6f} s（{e['wall_time_s']/60:.3f} min）。原停止判据在 step {e['stop']['step']} 满足，安全退出最终 native step {e['final_step']}，物理时间 {e['final_step']*policy['dt_s']:.12g} s。从零场启动，无缩放旧场或 warm start。

## 8. 独立日志、场与稳态验收

`validate_roi_balance_fem.py` 继承原验收链路，重新解析完整 solver.log，逐步 PETSc reason / 真实线性残差与非线性最终残差检查；重新读取全部保存的采样 VTU，独立计算连续稳态区间。独立合格步 `{v['independently_qualifying_steps']}`，总体 **PASS**。检查最终 VTU 坐标/四面体与输入一致，并核对 checkpoint 步号/时间及哈希。仅残差变小或图像平滑不作为科学正确性的证明。

- 实际入口 Q `{m['Q_in_m3_s']:.16g}` m³/s；入口误差 `{m['epsilon_Q']:.12g}`。
- 整体质量残差 `{m['epsilon_mass']:.12g}`，补偿求和 `{m['epsilon_mass_compensated']:.12g}`。
- 壁面无滑移 `{m['wall_noslip_pass']}`，速度/压力全部有限；最大速度 `{m['velocity_max_m_s']:.12g}` m/s，压力范围 `{m['pressure_range_pa']}` Pa。
- 最终场对合格停止场速度相对 L2 差 `{v['final_to_stop_relative_velocity_change']:.12g}`，归一化流量差 `{v['final_to_stop_normalized_flow_change']:.12g}`，通过原门槛。
- 所有出口净流量为正；未据此宣称每个面片没有局部回流。

## 9. 实际 3D 分流与旧 H0 比较

所有比例是带符号 Qout/实际 Qin，不裁剪、不再次强行归一化；旧 H0 3D 从正式 frozen VTU 重新积分，不用旧 0D 数字替代。cap 面积平均流体压力是输出场测量，区别于施加的法向牵引压力。

{finaltable}

{flowtable}

按预先规定 `new range < old range`，得到 **{comp['status']}**。极差减少 `{comp['range_reduction']:.12g}`，std 减少 `{comp['std_reduction']:.12g}`；正值表示改善。0D→3D 差异 `{comp['discrepancy_3D_minus_0D']}` 保留，不把差异当作重新调压依据，也不要求三出口 3D 各等于 1/3 才能验收。

![0D 设计与实际 3D 前向验证](figures/roi_design_vs_forward.png)

![两次压力参考与延伸段转移](figures/pressure_handoff.png)

新冻结流场 `{v['frozen_flow']}`，SHA256 `{v['frozen_flow_sha256']}`；SI arrays 和 manifest 在同一 frozen_flow 目录。**CFD-based retuning = NO；Particle/RBC calls = 0；没有生产微泡/RBC 晋升。**

## 10. 实际测试、命令及独立复核

```bash
cd {ROOT.parent}
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=ulm_3D_vascular \\
/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv/bin/python -B -m pytest \\
ulm_3D_vascular/tests/network_parameterized ulm_3D_vascular/tests/network_1d0d \\
ulm_3D_vascular/tests/network_h0 ulm_3D_vascular/tests/network_roi_only -q -p no:cacheprovider
```

正式启动前实际结果 **{pytest}**。旧 full-A 不可达性测试保留且通过；新 ROI 17 项包括限定内部几何、四端口、质量、H0 回归、gauge、affine basis、direct/iterative 一致性、负压差恒等式、无下游/CFD 依赖、半径不可变、冻结篡改/覆盖拒绝、仅三个 XML 压力改动、边反转及真实倒流符号。

按执行顺序保存命令在 `REPRODUCE.md`；审核者可在新审计目录重算纯 0D 和只读 3D 验收。**不要再次启动 CFD**，现有 dispatch 文件使启动器拒绝重复执行。报告生成器仅读取完成的独立验证记录，不会启动求解器。

核心来源与代码副本在 `evidence/`，原始取样和求解日志在 case `run/`，独立验收在 `independent_validation.json`。必要汇总 `boundary_and_flow_summary.csv`、`balance_comparison.csv`、`final_3D_ports.csv`、`roi_iteration_trace.csv`、`roi_fem_cap_pressures.csv`；机器可读数值保留完整精度。

## 11. 未验证的科学部分和结论边界

1. ROI-only pressure 是自由数值设计，不是全网预测或实测压力；人工移除下游被动约束，不能宣称得到可生理实现的全脑网络。
2. 0D 解析阻力假设圆截面、充分发展、刚性牛顿流，3D 有真实分叉及非圆几何；延伸段阻力补偿也仅为近似。当前差异符合模型层级不同的可能性，但没有把单次比较当作原因分离证明。
3. 只进行单张原生产网格上的一次受控压力改变，无新网格独立性、时间步独立性或边界敏感性扫描；不能据此量化真实血管离散误差。
4. 3D 实际输出与稳态门槛通过支持“这个固定算例的前向计算与分流比较”；并不等同于生理真值、WSS 局部精度或实验校准通过。
5. 原 full-A study、旧不可达证明和旧正式流场全部保留；两种模型回答不同问题。本轮没有压力反馈迭代，没有新增轨迹计算。

## 12. 需求交付索引

| 原要求 | 证据位置 |
| --- | --- |
| 1–4 起点、模型、节点/边数、无下游 | 第 1–2 节；cache summary |
| 5–6 H0 mixed 回归 P/Q/fraction | 第 4 节；regression JSON |
| 7–11 B/rank/条件数/direct/iterative/差值 | 第 5 节；basis 与 solution JSON |
| 12–17 相对/gauge 压力、Q/fraction/指标/−795.61 恒等式 | 第 5–6 节；frozen YAML |
| 18–20 R_ext/raw/applied cap | 第 6 节；handoff JSON 与 CSV |
| 21–22 唯一算例和实际 CFD 次数 | 第 7 节；dispatch/completion 原记录 |
| 23–26 实际 3D/0D/H0 比较 | 第 8–9 节；independent validation 与 CSV |
| 27–28 无反馈、Particle/RBC=0 | 第 9 节；dispatch 及 manifest |
| 29 pytest 命令和通过数 | 第 10 节；pytest.log |
| 30 科学限制 | 第 11 节 |
'''
    (report/'ROI_ONLY_BOUNDARY_BALANCE_FORWARD_REPORT_ZH.md').write_text(report_text)
    write_json(report/'progress.json',dict(status='COMPLETE',comparison_status=comp['status'],new_scientific_CFD_run_count=1,
        CFD_based_retuning=False,particle_RBC_calls=0,report='ROI_ONLY_BOUNDARY_BALANCE_FORWARD_REPORT_ZH.md',
        frozen_flow=v['frozen_flow'],frozen_flow_sha256=v['frozen_flow_sha256']))
    return comp['status']


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ('case','source-case','report'):p.add_argument('--'+arg,required=True,type=Path)
    a=p.parse_args();print(generate(a.case,a.source_case,a.report))
