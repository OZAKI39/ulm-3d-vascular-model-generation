"""Package completed diagnostics without changing any frozen scientific gate."""
from pathlib import Path
import json,hashlib,datetime
from zoneinfo import ZoneInfo
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;OUT=HERE/'error_diagnosis'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def put(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def main():
 d=read(ROOT/'V2_FIELD_ERROR_DECOMPOSITION.json');p=read(ROOT/'P1_FIELD_PARTICLE_SENSITIVITY.json');proposal=read(ROOT/'PROPOSED_V3_MESH_POLICY.json')
 frozen=read(ROOT/'FIELD_VALIDATION.json');identity=read(ROOT/'GPU_MANUFACTURED_RUN_IDENTITY.json')
 assert frozen['status']=='REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC'
 assert d['input_sha256']['FIELD_VALIDATION.json']==sha(ROOT/'FIELD_VALIDATION.json')
 assert d['input_sha256']['GPU_MANUFACTURED_RUN_IDENTITY.json']==sha(ROOT/'GPU_MANUFACTURED_RUN_IDENTITY.json')
 for row in p['comparison']:
  a=next(c for c in p['cases'] if c['case']==row['case'] and c['field']=='EXACT_ANALYTIC_FIELD')
  b=next(c for c in p['cases'] if c['case']==row['case'] and c['field']=='SAME_MESH_P1_ANALYTIC_FIELD')
  with np.load(ROOT/a['data_path']) as f:aa={k:f[k] for k in f.files}
  with np.load(ROOT/b['data_path']) as f:bb={k:f[k] for k in f.files}
  for c in [a,b]:assert sha(ROOT/c['data_path'])==c['data_sha256']
  t=np.linspace(0,.2,6401)
  def held(q):
   index=np.searchsorted(q['time_s'],t,side='left');index=np.minimum(index,len(q['time_s'])-1)
   return np.linalg.norm(q['step_velocity_m_s'][index],axis=1)
  sa,sb=held(aa),held(bb)
  row['held_integration_speed_relative_L2']=float(np.sqrt(np.trapezoid((sb-sa)**2,t)/np.trapezoid(sa**2,t)))
  row['initial_RHS_exact_m_s']=a['initial_RHS_velocity_m_s'];row['initial_RHS_P1_m_s']=b['initial_RHS_velocity_m_s']
  row['exact_local_transit_time_s']=a['transit_time_s'];row['P1_local_transit_time_s']=b['transit_time_s']
  row['accepted_velocity_check_max_m_s']=max(float(np.max(abs(np.diff(q['position_m'],axis=0)/np.diff(q['time_s'])[:,None]-q['step_velocity_m_s'][1:]))) for q in [aa,bb])
 p['additional_postprocessing']=dict(held_velocity_definition='Original accepted step velocity on its own preceding physical-time interval; sampled by timestamp on the same common grid, not by step index. Primary preregistered RHS-speed metrics unchanged.',
  original_RHS_interpolation_detail='Linear interpolation of three RHS velocity components, then speed norm; all measured transverse velocities are roundoff-scale.',
  script_sha256=sha(Path(__file__)),no_new_particle_runs=True)
 p['validation_checks']=dict(exact_geometry_tests='3 passed',decomposition_regressions='2 passed',production_source_files_unchanged=40,
  all_eight_runs_finite_no_penetration=True,all_eight_reached_observation_end=all(c['status']=='OBSERVATION_END' for c in p['cases']),
  same_contact_categories=all(c['contact_categories_equal'] for c in p['comparison']))
 put(ROOT/'P1_FIELD_PARTICLE_SENSITIVITY.json',p)
 rows=[r for r in d['groups'] if r['role']!='WALL_STRESS_DIAGNOSTIC']
 labels={'P1_SPACE_LIMIT_SUPPORTED':'A：P1_SPACE_LIMIT_SUPPORTED','DISCRETE_SOLUTION_LIMIT_SUPPORTED':'B：DISCRETE_SOLUTION_LIMIT_SUPPORTED','MIXED_ERROR':'C：MIXED_ERROR','PASS':'D：PASS'}
 table='\n'.join(f"| {r['wall_distance_m']*1e6:.6f} | {r['velocity']['total_relative_L2']*100:.6f}% | {r['velocity']['space_relative_L2']*100:.6f}% | {r['velocity']['discrete_relative_L2']*100:.6f}% | {r['velocity']['cross_term_normalized']:.9g} | {r['velocity']['alignment_cosine']:.6f} | {labels[r['classification']]} |" for r in rows)
 absolute='\n'.join(f"| {r['wall_distance_m']*1e6:.6f} | {r['velocity']['exact_L2']:.6g} | {r['velocity']['exact_RMS']:.6g} | {r['velocity']['total_absolute']['RMS']:.6g} | {r['velocity']['total_absolute']['RMS_mm_s']:.6g} | {r['velocity']['total_absolute']['median']*1e3:.6g} | {r['velocity']['total_absolute']['P95']*1e3:.6g} | {r['velocity']['total_absolute']['max']*1e3:.6g} | {r['velocity']['total_absolute']['RMS_over_Umean']*100:.6g}% |" for r in rows)
 components='\n'.join(f"| {r['wall_distance_m']*1e6:.6f} | {r['axial_absolute_m_s']['total']['RMS']:.6g} | {r['transverse_absolute']['CFD']['RMS']:.6g} | {r['transverse_absolute']['CFD']['RMS_um_s']:.6g} | {r['transverse_absolute']['CFD']['RMS_over_Umean']:.6g} | {r['transverse_absolute']['P1']['RMS']:.1g} |" for r in rows)
 derivatives='\n'.join(f"| {r['wall_distance_m']*1e6:g} | "+' | '.join('/'.join(f"{r['derivatives'][k][part]*100:.3f}" for part in ['space_relative_L2','discrete_relative_L2','total_relative_L2']) for k in ['G','s','curl'])+' |' for r in rows if r['role']=='PARTICLE_RELEVANT')
 trajectory='\n'.join(f"| {v['case']} | {v['initial_RHS_velocity_relative_error']*100:.6f}% | {v['speed_history_relative_L2']*100:.6f}% | {v['held_integration_speed_relative_L2']*100:.6f}% | {v['axial_displacement_relative_error']*100:.6f}% | {v['transit_relative_error']*100:.6f}% | {v['maximum_radial_separation_m']*1e6:g} | {v['maximum_gap_difference_m']*1e6:g} | {v['classification'].removeprefix('P1_FIELD_TRAJECTORY_IMPACT_')} |" for v in p['comparison'])
 transit='\n'.join(f"| {v['case']} | {np.linalg.norm(v['initial_RHS_exact_m_s'])*1e3:.9f} | {np.linalg.norm(v['initial_RHS_P1_m_s'])*1e3:.9f} | {v['exact_local_transit_time_s']*1e3:.6f} | {v['P1_local_transit_time_s']*1e3:.6f} |" for v in p['comparison'])
 walltable='\n'.join(f"| {v['wall_distance_m']*1e6:g} | {v['space_relative_L2']*100:.6f}% | {v['h_wall_estimate_m']*1e6:.6f} |" for v in proposal['wall_estimates'])
 r377=rows[0];v377=r377['velocity'];tangent=proposal['tangential_design_estimate']
 pro=d['radial_profile']['fields']['P1']['stable_thresholds']
 radial='\n'.join(f"| {float(k)*100:g}% | {v['first_ten_consecutive_samples_start_m']*1e6:.3f}–{v['first_ten_consecutive_samples_end_m']*1e6:.3f} | {v['first_below_and_stays_to_center_m']*1e6:.3f} |" for k,v in pro.items())
 locality=[]
 for group in d['discrete_bulk_localization']:
  rr=group['by_axial_plane'];corr=group['correlation_with_discrete_error_magnitude']
  locality.append(f"- {group['wall_distance_m']*1e6:.3f} µm：z/L=0.4/0.5/0.6 的 CFD−P1 RMS 为 "+'/'.join(f"{v['RMS']*1e3:.6f}" for v in rr)+f" mm/s；轴向/横向 RMS 为 {group['axial_absolute_m_s']['RMS']*1e3:.6f}/{group['transverse_absolute_m_s']['RMS']*1e3:.6f} mm/s；所列相关系数最大绝对值 {max(abs(v) for v in corr.values() if v is not None):.3f}。")
 receipt=read(ROOT/'monitoring/v2_p1_particle_sensitivity/EXECUTION.json')
 report=f'''# Frozen V2 field error decomposition 与 particle sensitivity

2026-10-04（Europe/Zurich）。本轮完成状态：**MIXED_FIELD_ERROR + P1_FIELD_TRAJECTORY_IMPACT_LARGE（DIAGNOSTIC_ONLY）**。旧正式判定仍为 `REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC`；本轮没有 CFD、mesh、solver、BC 或 gate 修改。分析使用冻结第71步及960原探针，粒子仅使用解析场与同网格解析 P1 场。

## 五个问题的结论

1. **0.377 µm 主要支持 P1 表示限制，且存在误差抵消。** space=9.758956%，CFD−P1=4.688421%，total=8.519477%；误差夹角余弦−0.487793。不能把百分数相减或把 total 解释为两范数的相加。
2. **bulk 并非同一原因。** 29.839/59.678 µm 的 space 仅0.380179%/0.272039%，失败主要来自 CFD−P1；149.195 µm 的 space=9.455478%，支持表示限制。298.390 µm 原本 PASS，保留 PASS，虽然较大的两项误差正在抵消。
3. **近壁 u 与导数的归一化尺度不同。** u 随距壁趋零，切向弦面带来的绝对插值亏损会被很小的 u_exact 放大；解析 G/s/curl 在近壁仍有有限尺度，且门槛为8%而非5%。测量中粒子带导数总误差约4.4–5.0%，因此通过。没有以 WSS 替代 RHS 输入。
4. **纯 P1 表示会传播到生产动力学。** 预声明0.2 s、单方位四场景中，A局部通过时间差5.148892%触发 LARGE；B/C为MODERATE，D为SMALL。没有接触类别改变、径向偏离或穿壁。LARGE 是用户定义的诊断类别，并不意味着出现失稳或已证明全部BraVa轨迹异常。
5. **若以后只允许一次V3，条件性建议为 wall altitude≤16 µm、相关中央core bulk target=120–130 µm；wall altitude单独变小不够。** 切向项还提示周向弦长≤约7.075 µm（532段时7.048 µm）。全部写入 `PROPOSAL_ONLY`，未经验证，不生成网格。29.8/59.7 µm 的离散解差异不能据此保证消失。

## 冻结身份、定义与恒等式

[原 FIELD_VALIDATION](FIELD_VALIDATION.json)、[GPU 身份](GPU_MANUFACTURED_RUN_IDENTITY.json)、final VTU/native checkpoint、mesh、960 probes和40个production粒子源码均经SHA核对。CFD probe u直接读取冻结 `FIELD_SAMPLES.npz`；从原节点重构的最大差仅 {d['CFD_reconstruction_max_roundoff_m_s']:.3g} m/s，为浮点舍入，没有替换原验收数值。

每个原mesh节点直接计算 u_exact，再用原tetra P1形函数构造 I_h u_exact及逐cell常量G。没有使用CFD节点值、壁面强制归零、gradient recovery或平滑；与上一轮为了数值rim清零的插值诊断之差仅浮点量级，不改变既有结论。

`e_total=u_CFD-u_exact`，`e_space=I_h u_exact-u_exact`，`e_discrete=u_CFD-I_h u_exact`。逐点恒等式最大残差 {max(r['velocity']['point_identity_max_abs'] for r in d['groups']):.3g} m/s；平方范数恒等式最大绝对残差 {max(r['velocity']['energy_identity_abs_residual'] for r in d['groups']):.3g} (m/s)²，均在预声明float64 roundoff范围。G/s/curl也验证同样恒等式。没有将范数相减。

cross列为 `2<e_space,e_discrete>/||u_exact||²`，是无量纲平方量，**不是百分数**。cos为两误差向量在整组探针上的夹角余弦；范数小于64 eps×参考范数时标N/A。以下全部cos<0，表示部分抵消。分类先保留原CFD PASS，再对失败组应用A/B/C，解决提示词规则的重叠；0.1 µm只在JSON保留诊断，不参与gate。

| 壁距 µm | total | space | CFD−P1 | cross | cos | 预声明分类 |
|---:|---:|---:|---:|---:|---:|---|
{table}

例如0.377 µm：`0.085194774² = 0.097589559² + 0.046884214² − 0.004463702`。A类支持的是本轮定义的**解析节点插值 oracle**达不到5%，不是对所有可能P1系数的最佳逼近误差下界的证明。

## 绝对误差与分量

L2=sqrt(96个探针向量模长平方和)，RMS=L2/sqrt(96)。所有行各96点。相对Umean采用RMS/0.1 m/s，仅诊断，不替换冻结gate。JSON另含space/discrete的全部绝对统计、m/s与mm/s数值。

| 壁距 µm | exact L2 m/s | exact RMS m/s | total RMS m/s | total RMS mm/s | median mm/s | P95 mm/s | max mm/s | RMS/Umean |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{absolute}

0.377 µm 的exact RMS只有 **0.252610 mm/s**，total RMS为 **0.021521 mm/s（21.521 µm/s）**，即Umean的0.021521%。这解释相对误差为何敏感，但不是撤销5%门槛的理由。

解析与P1场的横向速度均为0，不做相对除零；下面只列绝对值与Umean归一化值。JSON包含轴向total/space/discrete，以及横向median/P95/max和µm/s换算。

| 壁距 µm | CFD axial error RMS m/s | CFD transverse RMS m/s | CFD transverse RMS µm/s | transverse/Umean | P1 transverse RMS m/s |
|---:|---:|---:|---:|---:|---:|
{components}

## 导数精简对照

每格依次为 **space / CFD−P1 / total，相对L2百分数**；G使用完整Frobenius，s使用原解析向内法向，curl直接来自同一G，未做WSS audit。

| 粒子壁距 µm | G：space/discrete/total % | s：space/discrete/total % | curl：space/discrete/total % |
|---:|---:|---:|---:|
{derivatives}

对二次Poiseuille，P1轴向误差有精确表达：`I_h u_z-u_z = -(2 Umean/R²) Σ w_i ||node_xy_i-x_xy||²`。将每条偏移分成当前位置的径向和周向，可得到逐点可相加的两个亏损；再次说明二者**范数**不能直接相加作等式。在0.377 µm带，径向项norm为 **{r377['quadratic_space_identity']['radial_relative_L2']*100:.6f}%**，切向项norm为 **{r377['quadratic_space_identity']['tangential_relative_L2']*100:.6f}%**。切向项较大，是不能只缩小wall altitude的直接证据；不是凭网格外观推断。

## 中央径向分布与bulk定位

新增剖面与原960 gate点分开保存：z=L/2，16条预定方位径向线，{d['radial_profile']['point_count']}点，0.1–25 µm对数密集、25 µm至中心线性采样。全部严格包含，未移动或外推。空间误差到shell区域下降、进入较粗core后又上升，因此“第一次低于阈值”不代表之后一直满足。

| P1阈值 | 首个连续10采样点满足的范围 µm | 此后直到中心全部满足的最早壁距 µm |
|---:|---:|---:|
{radial}

这些是离散radial采样上的描述，既不是更改960点gate，也不是连续全域误差上界。

对29.8/59.7 µm仅做已有96点的CFD−P1分布、轴向/横向分解，以及端距、所属cell质心至25/50/75 µm shell距离、cell边轴向/径向朝向的Pearson相关性；预声明|r|≥0.5才标显著关联。结果：

{chr(10).join(locality)}

未发现对上述单一特征的明显相关定位，但中央面z/L=0.5比0.4/0.6误差较低。原点仅覆盖中央0.4–0.6L，不能据此排除端部影响、推断特定BC或求解器项造成误差；本轮没有继续solver研究。

## 生产粒子敏感性：只含exact与P1 exact

A：a=0.375 µm,gap=2 nm；B：a=0.964365028029 µm,gap=2 nm；C：同median,gap=2 µm；D：同median,gap=10 µm。两个arm均用相同exact cylinder gap/normal、radius、mu、near-wall阻力、planar affine block、handoff/contact、原 `v1_trial` bytecode与 `refine_interval`。验证适配器只替换解析几何查询及其精确连续路径证书，没有投影位置、改变lower-gap或手写另一套积分器。生产源码没有修改。

开始前固定：theta为原probe第一个方位，z=L/2，dt_max=0.125 ms，观察0.2 s，局部通过距离10 µm。使用服务器现有CPU生产路径，4个独立进程；exact sampler为CPU实现，按本轮授权不新建CUDA kernel。监督任务耗时 **{receipt['elapsed_s']:.3f} s**、峰值进程树RSS **{receipt['peak_tree_RSS_GiB']:.3f} GiB**，退出码0。四scenario×两field，共8条，每条1600个接受步、0拒绝、无穿壁，均实际通过10 µm，无censoring。

主比较按共同31.25 µs物理时间网格，线性插值即时RHS速度分量和位置，再取速度模长；不是按step index。另从原接受步重建分段常量的实际积分速度，作为附加核对，两种speed口径结论一致。0.125 ms只减少本次时间离散污染，没有证明其已收敛。

| Case | 初始RHS差 | RHS speed历史L2 | 原积分speed历史L2 | 轴向位移差 | 局部transit差 | 最大径向分离 µm | 最大gap差 µm | 诊断 |
|---|---:|---:|---:|---:|---:|---:|---:|---|
{trajectory}

| Case | exact初始speed mm/s | P1初始speed mm/s | exact局部transit ms | P1局部transit ms |
|---|---:|---:|---:|---:|
{transit}

A/B始终是 `CONTINUUM_HANDOFF_CONTACT`，C/D始终 `FAR_FIELD`，两arm一致；2 nm handoff不是分子级硬接触。径向位置和gap在本次轨迹保持一致。A的速度约降低4.897%，固定距离时间随速度倒数变化，因此时间差约5.149%，超过5%参考线；这不是突发失稳。B/C速度存在随cell变化的差异，原始P1 piecewise gradient未平滑。

总体按用户规则为 **{p['status']}**。只覆盖一个预定方位、0.2 s及10 µm局部通过时间，不能泛化为全管transit，也不能否定正式5% field gate。没有证据支持本轮写 `CURRENT_FIELD_GATE_MAY_BE_OVERCONSERVATIVE_FOR_TRAJECTORY`。

## 条件性V3建议，未实施

使用用户指定的 `h_new=h_old sqrt(0.04/E_space)`，h_old采用实测wall altitude=24.998332 µm：

| 粒子壁距 µm | P1 space | h_wall一阶估计 µm |
|---:|---:|---:|
{walltable}

取最严格16.0044 µm并向下取 **16 µm**。这仅为h²设计估计，**不是已验证收敛率**；对固定近壁probe单独改法向高度不保证h²缩放。实测切向项7.905%已经大于4%，因此16 µm不能单独当作合格方案。

切向条件估计保守地不计径向项会改善，要求 `E_radial + E_tangent*(chord_new/chord_old)² ≤4%`。得到弦长≤**{tangent['allowable_chord_m']*1e6:.6f} µm**，周向272→**{tangent['proposed_intervals_multiple_of_four']}**段时弦长**{tangent['proposed_chord_m']*1e6:.6f} µm**、理论弦高**{tangent['proposed_sagitta_m']*1e9:.6f} nm**。圆柱R/L不变，但若未来授权实施，cap rim必须共形；本轮没有改CAD、PLC、cap或mesh。

仅149.195 µm这个representation-supported失败bulk带参与估计：200 µm名义target×sqrt(4/9.455478)=130.082314 µm，建议相关中央core的名义target **120–130 µm**。不因为29.8/59.7 µm total约5–5.6%而全局缩小bulk；它们的P1 oracle已经很好，新的局部尺寸也不能保证CFD−P1差异自动消失。298.390 µm保持原PASS，不将其作为新增失败驱动加密。

[PROPOSED_V3_MESH_POLICY.json](PROPOSED_V3_MESH_POLICY.json) 状态 **PROPOSAL_ONLY**。没有启动mesher、CFD、正式oracle、dt sweep或BraVa。下一步须由用户决定是否采纳/修改提案，不能自动执行。

## 交付与验证

两张新图：[误差分解与径向剖面 PNG](gmsh_reference/figures/v2_error_decomposition.png)；[四sentinel诊断 PNG](gmsh_reference/figures/p1_particle_sensitivity.png)。均300 dpi PNG；[HTML入口](OPEN_RESULTS.html)。

机器结果：[V2_FIELD_ERROR_DECOMPOSITION.json](V2_FIELD_ERROR_DECOMPOSITION.json)、[P1_FIELD_PARTICLE_SENSITIVITY.json](P1_FIELD_PARTICLE_SENSITIVITY.json)。逐点分解、radial采样、8条接受态及共同时间比较保存在 `gmsh_reference/error_diagnosis/`，SHA写入结果；预声明 [particle protocol](gmsh_reference/error_diagnosis/PARTICLE_DIAGNOSTIC_PROTOCOL.json) 保留。

5项有针对性的测试通过（解析gap/contact、连续路径与handoff事件、向量抵消、二次函数P1梯度）。服务器任务已退出并移除TASK；原FIELD_VALIDATION、GPU identity、steady VTU/checkpoint和40个生产源码身份保持。
'''
 (ROOT/'V2_P1_REPRESENTATION_DIAGNOSIS_ZH.md').write_text(report)
 path=ROOT/'MANUFACTURED_TRAJECTORY_FIELD_VALIDATION_ZH.md';old=path.read_text();marker='<!-- V2_ERROR_DIAGNOSIS -->'
 if marker in old:old=old.split(marker)[0].rstrip()+'\n'
 path.write_text(old+f'''\n{marker}
## 后续冻结场诊断（2026-10-04，fd6c068c）

本轮没有重跑CFD或改原gate。误差分解结论为 **MIXED_FIELD_ERROR**：0.377/149.195 µm支持P1表示限制，29.839/59.678 µm支持CFD−P1离散解差异；各项存在负交叉项，不能相减百分数。0.377 µm的P1误差以切向项为主，不能仅凭缩小wall altitude保证通过。

新授权的exact-vs-P1四场景诊断已完成，**P1_FIELD_TRAJECTORY_IMPACT_LARGE（DIAGNOSTIC_ONLY）**，由A的10 µm局部通过时间误差5.148892%触发；四组无径向分离、gap差或contact类别改变。它未使用CFD场，不是正式oracle/CFD trajectory certification；正式dt收敛仍未运行。

详见[误差分解与敏感性报告](V2_P1_REPRESENTATION_DIAGNOSIS_ZH.md)。[V3条件性策略](PROPOSED_V3_MESH_POLICY.json)仅为PROPOSAL_ONLY，未生成网格。原 `REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC` 继续有效。
''')
 review=f'''# Production workflow 当前结论

2026-10-04。当前任务完成并停止：**MIXED_FIELD_ERROR + P1_FIELD_TRAJECTORY_IMPACT_LARGE（DIAGNOSTIC_ONLY）**。原正式field状态继续是 **REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC**。

| 本轮问题 | 结论 |
|---|---|
| 0.377 µm 的8.52%误差来自哪里？ | P1 space=9.759%、CFD−P1=4.688%，cos=−0.488，存在抵消；P1空间项中切向几何项比径向项大。 |
| 29.8–149.2 µm bulk是否同样受P1限制？ | 29.8/59.7 µm主要是离散解差异；149.2 µm支持P1表示限制，不能统一归因。 |
| 为什么G/s/curl通过而部分u失败？ | 近壁u_exact很小且受切向插值亏损影响，导数有有限尺度、采用8%门槛；原粒子带导数实测约4.4–5.0%。 |
| P1表示是否影响生产粒子轨迹？ | A局部transit差5.149%触发LARGE，B/C为MODERATE，D为SMALL；无径向分离或contact类别变化。仅0.2 s单方位诊断。 |
| 若以后只做一次V3，应改多大？ | 条件估计：wall altitude≤16 µm，周向弦长≤7.075 µm（532段时7.048 µm），相关中央core名义bulk target=120–130 µm；全部PROPOSAL_ONLY。 |

[完整诊断报告](V2_P1_REPRESENTATION_DIAGNOSIS_ZH.md) · [机器分解](V2_FIELD_ERROR_DECOMPOSITION.json) · [粒子敏感性](P1_FIELD_PARTICLE_SENSITIVITY.json) · [V3提案](PROPOSED_V3_MESH_POLICY.json) · [两张新图与HTML入口](OPEN_RESULTS.html)。

已完成：冻结输入SHA审计、960点三项向量分解/交叉项/绝对量/分量分析、G/s/curl简表、3840点radial profile、两个discrete-limited bulk带的零新CFD定位、4×2条production动力学诊断、5项针对性测试。粒子服务器计算约12.03 s，4CPU进程、峰值RSS1.269 GiB；无新增CFD或meshing。

原GPU CFD已在88.91 min内达到连续五区间稳态并保存第71步，CUDA/RA身份保持。原field的失败没有被覆写；[流场报告](MANUFACTURED_TRAJECTORY_FIELD_VALIDATION_ZH.md)继续有效。正式CFD-vs-oracle轨迹和三个dt的canonical convergence仍未运行，不能写READY，也不能宣称dt=0.5 ms已获支持。

服务器job已退出、TASK已移除，无自动下一任务。V3建议只是本地设计模型，不是已验证h²收敛率或PASS保证。等待用户新任务指令，不自动改gate、不生成V3、不运行CFD、不进入BraVa。

正式BraVa `cases/balanced_pressure_final` 与旧1565条population保持原身份；141观察截止、2数值失败不改称成功，新固定灌注协议未运行新正式population。未来新BraVa边界设计仍按“中心线阻力估算→出口压力”。
'''
 (ROOT/'PRODUCTION_WORKFLOW_REVIEW_ZH.md').write_text(review)
 state=dict(status='MIXED_FIELD_ERROR',VALIDATION_ONLY=True,current_authorization='fd6c068c frozen V2 decomposition + exact/P1 particle diagnostic ONLY',reason='Mixed representation and discrete-solution errors; prescribed diagnostic completed, proposal only.',formal_field_status=frozen['status'],diagnostic_particle_status=p['status'],diagnostic_scenarios=4,diagnostic_particle_runs=8,new_CFD=0,CPU_CFD_reruns=0,new_meshes=0,new_formal_oracle=0,formal_CFD_trajectory_runs=0,new_timestep_cases=0,new_BraVa=0,new_population=0,active_jobs=[],automatic_next_stage=False,automatic_BraVa=False,requires_new_user_instruction_to_continue=True,server_TASK_removed=True,completed_local=datetime.datetime.now(ZoneInfo('Europe/Zurich')).isoformat(),diagnostic_result_sha256=sha(ROOT/'V2_FIELD_ERROR_DECOMPOSITION.json'),diagnostic_particle_sha256=sha(ROOT/'P1_FIELD_PARTICLE_SENSITIVITY.json'),proposal_status='PROPOSAL_ONLY',tests={'diagnostic_geometry_and_decomposition':'5 passed in 0.83s'},production_physics_changed=False,production_dt_s=.0005,diagnostic_dt_max_s=.000125,task_target_wall_time_s=1800,task_hard_bound_s=3600,diagnostic_execution_wall_time_s=receipt['elapsed_s'],diagnostic_execution_exit_code=0,FIELD_VALIDATION_sha256=sha(ROOT/'FIELD_VALIDATION.json'),GPU_execution_identity_sha256=sha(ROOT/'GPU_MANUFACTURED_RUN_IDENTITY.json'),frozen_CFD=dict(case=identity['case'],status=identity['status'],final_step=71,wall_time_s=identity['wall_time_s'],rerun=False))
 put(ROOT/'CURRENT_STATE.json',state)
 path=ROOT.parent/'ACTIVE_DEVELOPMENT_SCOPE.md';text=path.read_text();start=text.index('当前停止点');end=text.index('默认只读取以下入口')
 text=text[:start]+'''当前停止点（2026-10-04，用户fd6c068c）：**MIXED_FIELD_ERROR + P1_FIELD_TRAJECTORY_IMPACT_LARGE（DIAGNOSTIC_ONLY）**。本轮只分析冻结GPU V2第71步，零新CFD/mesh。960点分解显示0.377/149.195 µm支持P1表示限制，29.839/59.678 µm支持离散解差异。4×2条exact-vs-P1解析场诊断完成（相同exact cylinder/生产physics，dt_max=0.125ms，0.2s）；A局部10µm通过时间误差5.149%触发LARGE，没有contact类别变化或径向分离。正式field仍为REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC；正式CFD-vs-oracle与dt convergence未运行。

权威新入口：`production_workflow/V2_P1_REPRESENTATION_DIAGNOSIS_ZH.md`、`V2_FIELD_ERROR_DECOMPOSITION.json`、`P1_FIELD_PARTICLE_SENSITIVITY.json`。`PROPOSED_V3_MESH_POLICY.json`仅为PROPOSAL_ONLY：wall altitude≤16µm、周向约532段/7.048µm弦长、相关central-core bulk target120–130µm，未实施、未验证收敛率。等待用户指令，不自动生成V3/CFD/BraVa或修改gate。服务器诊断job已退出，TASK已移除。

历史冻结：V2 mesh永久冻结，GPU production 8MPI/RTX4090/RA已在88.91min达到第70步五区间稳态，第71步正常保存；原field数值与身份不变。旧CPU1800s/step15停止仅作历史证据，不能用step10作新认证。生产粒子源码与正式灌注协议未改动。

'''+text[end:];path.write_text(text)
 path=ROOT/'OPEN_RESULTS.html';text=path.read_text();head=text[:text.index('<h1>')];old=text[text.index('<h1>'):text.rfind('</html>')]
 if '<!-- ERROR_DIAGNOSIS_CURRENT -->' in text:
  old=text.split('<details id="previous_gpu_result"><summary>此前GPU稳态与field验收</summary>',1)[1].rsplit('</details></html>',1)[0]
 page='''<!-- ERROR_DIAGNOSIS_CURRENT --><h1>冻结V2：误差分解与粒子敏感性诊断完成</h1><p><strong>MIXED_FIELD_ERROR + P1_FIELD_TRAJECTORY_IMPACT_LARGE（DIAGNOSTIC_ONLY）</strong>。最靠壁及149.2 µm支持P1表示限制，29.8/59.7 µm支持CFD−P1差异。A局部transit差5.149%，四组无径向分离或contact类别改变。原field FAIL保留；本轮零新CFD、零mesh。</p><p><a href="V2_P1_REPRESENTATION_DIAGNOSIS_ZH.md">完整诊断报告</a> · <a href="PRODUCTION_WORKFLOW_REVIEW_ZH.md">五个问题</a> · <a href="V2_FIELD_ERROR_DECOMPOSITION.json">分解JSON</a> · <a href="P1_FIELD_PARTICLE_SENSITIVITY.json">粒子JSON</a> · <a href="PROPOSED_V3_MESH_POLICY.json">仅提案的V3策略</a></p><figure><a href="gmsh_reference/figures/v2_error_decomposition.png"><img src="gmsh_reference/figures/v2_error_decomposition.png" alt="Frozen V2 error decomposition and radial profile"></a><figcaption>原gate不变，误差向量有负交叉项。</figcaption></figure><figure><a href="gmsh_reference/figures/p1_particle_sensitivity.png"><img src="gmsh_reference/figures/p1_particle_sensitivity.png" alt="Four exact versus P1 production particle diagnostics"></a><figcaption>同一解析几何/生产动力学，0.125 ms最大步长、0.2 s、10 µm局部transit。没有使用CFD场或做dt收敛证明。</figcaption></figure><p>条件性V3估计：wall altitude≤16 µm，周向弦长约7.05 µm，相关中央core bulk target120–130 µm；PROPOSAL_ONLY，未实施。任务已退出，无自动后续任务。</p>'''
 path.write_text(head+page+'<details id="previous_gpu_result"><summary>此前GPU稳态与field验收</summary>'+old+'</details></html>')
 print('Reports/state/HTML finalized; original FIELD_VALIDATION and GPU identity unchanged.')
if __name__=='__main__':main()
