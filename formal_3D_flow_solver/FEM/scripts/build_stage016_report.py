#!/usr/bin/env python3
"""Finalize the failed study only after tests and history/reference audits."""
import json,sys,xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import sha256,write_json,timestamp
R=ROOT/'reports/stage01_6';O=ROOT/'outputs/stage01_6'
def read(p):return json.loads(p.read_text())
selection=read(R/'candidate_selection.json');contract=read(R/'planar_port_contract_v2.json');policy=read(R/'acceptance_policy.json');lock=read(R/'freeze_lock.json')
assert selection['status']=='FAIL' and selection['selected_candidate'] is None
history=read(R/'history_preservation.json');references=read(R/'reference_integrity.json')
assert history['status']==references['status']=='PASS'
xml=ET.parse(R/'pytest_results.xml').getroot();suites=list(xml.iter('testsuite'))
counts={key:sum(int(s.get(key,0)) for s in suites) for key in ('tests','failures','errors','skipped')};counts['passed']=counts['tests']-counts['failures']-counts['errors']-counts['skipped']
assert counts['failures']==counts['errors']==0
rows=selection['records'];qcs=[r['surface_qc'] for r in rows];geometry=[p['geometry'] for qc in qcs for p in qc['ports'].values()]
maximums={'rim_displacement_m':max(p['rim']['maximum_rim_displacement_m'] for p in geometry),
 'projected_area_relative_error':max(max(p['relative_projected_area_error'],p['independent_projected_area_relative_error']) for p in geometry),
 'vector_area_relative_error':max(p['relative_vector_area_error'] for p in geometry),
 'projected_centroid_displacement_m':max(p['projected_centroid_displacement_m'] for p in geometry),
 'independent_centroid_difference_m':max(p['independent_projected_centroid_error_m'] for p in geometry),
 'interior_plane_deviation_m':max(p['new_interior_max_plane_deviation_m'] for p in geometry)}
base=read(ROOT/'outputs/stage01/medium/qc/geometry_qc.json');bq=base['quality'];proxy=policy['cost']['baseline']
def table(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,row))+' |' for row in rows])
caprows=[]
for row in rows:
 s=row['surface_qc'];caprows.append([row['candidate'],f"{policy['candidates'][row['candidate']]*1e6:.2f}",s['total_cap_triangles'],sum(p['interior_vertex_count'] for p in s['ports'].values()),f"{s['combined_cap_quality']['P5']:.6f}",f"{s['combined_cap_quality']['median']:.6f}",'质量通过；密度 FAIL'])
densityrows=[[name.upper(),p['original_cap_triangle_count'],policy['density']['per_port_max'][name]]+[s['ports'][name]['new_triangle_count'] for s in qcs] for name,p in contract['ports'].items()]
legacyrows=[[n.upper(),f"{p['legacy_scalar_area_m2']:.15e}",f"{p['formal_projected_area_m2']:.15e}",f"{p['relative_difference_legacy_scalar_vs_projected']:.3e}",f"{p['legacy_max_nonplanarity_m']:.3e}"] for n,p in contract['ports'].items()]
scalarrows=[[n.upper()]+[f"{s['ports'][n]['geometry']['candidate_scalar_vs_legacy_relative_difference']:.3e}" for s in qcs] for n in contract['ports']]
figures=read(R/'visualization_manifest.json');figlist='\n'.join(f'- [{name}]({name}) — {meta["role"]}' for name,meta in figures['images'].items())
report=f'''# Stage 1.6 — 平面端口契约与稀疏端盖重网格

**STAGE 1.6 STATUS: FAIL。三个候选的几何和 cap 质量通过，但全部超过冻结的端盖密度预算。没有 winner，没有生成候选四面体网格。**

## Stage 1.5 为什么失败？

Stage 1.5 的 A/B/C 保持 wall、rim、plane 和标签不变，端盖三角形质量明显改善，却未通过旧的 3D scalar triangle area 相对误差 ≤1e-12 门槛。其误差约 6.3e-10～3.0e-9；独立 longdouble 复核表明这不是 double 求和误差。原始 float32 VTP 的 rim 存在约 6e-12～1.5e-11 m 的非共面性，改变内部剖分就可能改变 3D 标量面积。

[Stage 1.5 分析](../stage01_5/area_invariant_analysis.json)和 **Stage 1.5 FAIL** 永久保留。本阶段没有把旧阈值改成 1e-8。

## Stage 1.6 改变了什么定义？

人工 CFD 端口现在由冻结的 3D rim 顶点/边、原始 x0、原始 outward normal、固定投影基和投影多边形定义。正式面积是 **A_projected_polygon**；正式中心是该平面多边形的面积中心。由有序 3D rim 直接积分得到的向量面积独立于 cap 内部剖分。

历史 scalar area 和 scalar-area centroid 继续保留为 provenance/diagnostic，每个候选仍记录实际 scalar triangle area。原始 rim 不移动、不投影替换；只有新增内部点放到原平面。WALL 是解剖表面，CAP 是人工 CFD 边界，修订人工端口面积的定义没有修改血管解剖。

[修订依据](PORT_CONTRACT_RATIONALE.md)、[v2 契约](planar_port_contract_v2.json)、[冻结策略](acceptance_policy.json)、[冻结哈希](freeze_lock.json)。Stage 0 source contract SHA256：`{contract['source_contract_sha256']}`。

{table(['端口','legacy scalar area / m²','formal projected area / m²','legacy 相对 planar 差异','rim 非共面性 / m'],legacyrows)}

候选 scalar area 相对 legacy scalar area 的**带符号变化**如下。全部公开，但不再作为端口身份 gate；两个正式面积不变量仍执行 ≤1e-12。

{table(['端口','sparse_A','sparse_B','sparse_C'],scalarrows)}

![端口定义](port_contract_explanation.png)

左图仅法向高度放大 10,000 倍，已明确标注；INLET 的真实最大偏离仅约 8.597 pm。

## 有没有改变血管几何？

没有移动 wall 或 rim。三个候选的 67,071 个 WALL 三角形坐标、连接及标签完全相同，maximum wall displacement = **0 m**，maximum rim displacement = **0 m**；rim vertex ID set、rim edge set、plane origin、reference normal、基向量和 entity mapping 完全一致。标签仍为 1=WALL、2=OUTLET_03、3=OUTLET_01、4=INLET、5=OUTLET_02，恰好 1 inlet、3 outlets。

每个曲面均为 1 个连通分量、0 条 boundary edges、0 条 nonmanifold edges，且朝向一致；每个 cap 是一个封闭 polygon loop，无孔洞、重叠或越界。覆盖检查同时使用边界身份、Euler disk 拓扑、点包含、边交叉和三角形相交面积，不以总面积相等替代拓扑检查。

两种独立面积计算为 longdouble rim shoelace 与实际投影三角形带符号面积求和。所有候选最大相对面积误差 **{maximums['projected_area_relative_error']:.3e}**；向量面积相对误差 **{maximums['vector_area_relative_error']:.3e}**。多边形中心位移 **{maximums['projected_centroid_displacement_m']:.3e} m**；三角形加权中心的独立交叉误差最大 **{maximums['independent_centroid_difference_m']:.3e} m**。新增内部点偏离原平面的最大值 **{maximums['interior_plane_deviation_m']:.3e} m**，小于 1e-15 m。原 rim 的历史非共面性原样保留。

![rim 对照](rim_overlay.png)

![边界标签](boundary_tags_selected.png)

图中 sparse_C 是失败候选的几何诊断示例，**不是 selected**。

## 为什么这次端盖更稀疏？

Stage 1.5 的三组 cap 分别有 2,539、1,787、1,397 个三角形，本次降到 1,143、975、887。运行前固定公式为：

`h(d) = h_rim + (h_center - h_rim) × min(d / (0.5 × sqrt(A_projected/π)), 1)`。

距离 d 是平面内到原 rim 线段的精确最小距离；h_rim 是原 3D rim 边长中位数。等效半径只用于尺寸过渡，没有把真实多边形拟合成圆。通过 [Gmsh 的尺寸回调接口](https://gmsh.info/doc/texinfo/#gmsh_002fmodel_002fmesh_002fsetSizeCallback)提供该确定性尺寸，二维算法固定为 6，每条原 rim 边约束为两个端点，不拆边。三个 center targets 为 0.35、0.45、0.55 µm。

策略在 **{lock['timestamp']}** 冻结，早于所有远端候选；SHA256 为 `{lock['policy_sha256']}`。结果出来后未改公式、中心尺寸、过渡距离或任何预算。这一冻结方案虽然比 Stage 1.5 稀疏，但仍保留较多近 rim 小单元，**没有达到用户要求的稀疏预算**。

## 三个 sparse candidate

{table(['候选','中心尺寸 / µm','cap triangles','新增 Steiner points','cap P5','cap median','结论'],caprows)}

{table(['端口','Stage 1','允许最大值','sparse_A','sparse_B','sparse_C'],densityrows)}

总数上限是 800，而逐端口上限之和为 764；两套预算都必须满足。三候选四个端口分别超预算，总数也均超 800。因此全部拒绝，不进入 volume meshing。

![剖分对照](cap_triangulation_before_after.png)

左右使用相同投影基、尺度和视角。右侧是预定最大中心尺寸的 sparse_C，明确标注 REJECTED。

## 新端盖是否仍然质量足够？

足够。按实际 3D 面积计算 `q_tri = 4√3 × Area / Σ(edge²)`，三个候选都没有 q_tri<0.1 的三角形，P5≥0.45、median≥0.70，所有面积为正且有限。原端盖 P5≈0.082、median≈0.091。二维质量改善不能覆盖密度预算失败。

![密度与质量](cap_density_quality_tradeoff.png)

## 生成体网格后改善了吗？

**没有执行，不能证明改善。** 表面 gate 在体网格之前失败，因此没有任何 Stage 1.6 tetra、tetra QC、candidate P2 proxy 或 selected mesh。体网格路径仍锁定 Stage 1 medium 的 Algorithm3D=1、bulk=5e-7 m 和相同优化选项，但未对失败候选启动。没有 HXT、Netgen、局部 volume refinement 或扩大搜索。

{table(['指标','Stage 1 medium','Stage 1.6 selected'],[['tetra',147569,'未生成'],['minimum minSICN',bq['gmsh_min_sicn']['minimum'],'N/A'],['P1',bq['gmsh_min_sicn']['P1'],'N/A'],['P5',bq['gmsh_min_sicn']['P5'],'N/A'],['median',bq['gmsh_min_sicn']['median'],'N/A'],['q<0.1',153,'N/A'],['cap-adjacent q<0.1',129,'N/A']])}

![体质量](tetra_quality_before_after.png)

![低质量单元位置分类](low_quality_count_by_boundary.png)

分类继续使用 tetra center 到 exterior triangle center 的最近邻，与 Stage 1 一致，是位置诊断，不声称精确拓扑邻接。基线分类为 WALL=24、INLET=51、OUTLET_01=26、OUTLET_02=16、OUTLET_03=36。

## 新最差单元在哪里？

新体网格不存在，因此没有新最差单元。原基线 worst cell 为 **137690**，minSICN=0.012171307762867093，最近边界为 INLET。图中保留基线实际 worst 20 和最差单元形状；位置标记放大与相机放大均已标注，没有改变 tetra 相对几何。

![最差单元](worst_elements_before_after.png)

![内部剖面](tetrahedral_cutaway_selected.png)

内部剖面取自真实 Stage 1 medium tetra；右侧明确注明 selected 不存在。没有用旧网格冒充新结果。

## 后续 FEM 会不会明显更贵？

本次无法给出新体网格成本结论。基线只按拓扑统计，没有创建 FEM 空间：

{table(['计数','Stage 1 medium'],list(proxy.items()))}

`P2 scalar = N_vertex + N_edge`，`P2 velocity = 3 × (N_vertex + N_edge)`，`P1 pressure = N_vertex`。冻结预算：tetra≤200,000，P2 velocity proxy≤1.35×801,405={1.35*801405:.2f}。这只是拓扑自由度成本代理，未估算具体求解时间或内存。

![成本比较](mesh_cost_comparison.png)

所有新体网格相关值都写为 N/A/未执行，不能从较好的 cap 质量推断未来 FEM 更便宜。

## 为什么选择 winner？

**没有选择 winner。** 第一层完整 surface gate 已淘汰全部候选。冻结规则要求先满足 geometry、cap quality、cap density、tetra quality、tetra budget 和 P2 budget，再依次按最少 P2 DOF、最少 tetra、最少 cap-adjacent low、最少 total low、最高 P1 排序。排序的后续指标没有被计算或虚构。

[选择记录](candidate_selection.json)、[质量对比](quality_comparison.json)、[成本对比](cost_comparison.json)。未生成 `outputs/stage01_6/selected/`；仅 winner 才允许执行的 1-rank/2-rank DOLFINx round-trip 均未执行。

## 自动测试

完整 pytest：**{counts['passed']} passed、{counts['failures']+counts['errors']} failed、{counts['skipped']} skipped**，见 [pytest_results.xml](pytest_results.xml)。Stage 1.6 新增二十个永久测试模块，共 75 passed、2 skipped。新增跳过项仅为无 selected 时禁止执行的两项 round-trip；历史跳过项另有三项。测试通过不等于候选符合实验门槛。

恶意回归覆盖：rim 坐标移动、拆边、改法向、改投影面积、反转向量面积、孔洞、重叠、cap 超密、质量合格但 tetra>200,000、质量优秀但 P2>1.35×基线，以及两个质量合格候选必须选择较低 DOF。体积审计在冻结基线上重算了原 minSICN 和最近边界计数；新候选 volume 路径未获得运行验证。

第一次完整测试被 SIGTERM（退出 143）中断，日志保留、未计为通过；随后完整重跑产生上述正式 XML。没有修改历史测试或降低断言。

历史完整性：{history['file_count']} 个 Stage 0/1/1.5/2 reports、outputs、logs、inputs 文件的集合、大小、SHA256 全部相同；所有预存 `src/fem3d/*.py`（含 Stage 2 solver core）哈希相同。开始前已有的 Stage 1 报告格式修改保持原字节，未被本阶段提交。两个只读参考工程均 **0 modified / 0 deleted / 0 added**（完整 {references['entries']} 项清单）。证据：[history_preservation.json](history_preservation.json)、[reference_integrity.json](reference_integrity.json)。

## 人工应该审核什么？

WSL 生成并逐图检查了下列十张图。由于没有 winner，四张依赖新体网格的对比图保留真实基线并显示“未生成”；两个表面图只展示明确标识的失败候选。文件名中的 selected 不代表存在 selected artifact。

{figlist}

## 还没有证明什么？

没有证明稀疏 cap 在冻结密度预算内可行，没有证明新 tetra 质量改善、P2 成本满足预算或 DOLFINx selected round-trip 通过。没有求真实 vascular flow，没有创建 velocity/pressure/lambda 空间，没有 Stokes、Navier–Stokes、WSS 或流线计算；没有 mesh convergence，也没有证明 Stage 3 solution 收敛。

远端仅运行三次 CPU cap meshing/QC（GPU used=false）；WSL 是唯一源码来源，沿用现有 probe/sync/run/fetch，所有新证据进入 stage01_6。环境、命令、源代码/配置/输入 SHA、stdout/stderr、运行时间和 RSS 随远端回传日志及候选 metadata 保留。

## 最终状态

**STAGE 1.6 STATUS: FAIL。** 失败原因为 sparse_A/B/C 全部违反逐端口及总 cap 密度预算。平面端口契约的几何验证和二维质量通过不足以把该阶段升级为 CONDITIONAL PASS。

未添加 sparse_D/E，未调参或放宽预算。Stage 1.5 仍为 FAIL，Stage 2 core 保持冻结。停止于 Stage 1.6，没有进入 Stage 3。
'''
(R/'REPORT.md').write_text(report)
write_json(R/'final_status.json',{'timestamp':timestamp(),'status':'FAIL','reason':'All three candidates exceed per-port and total cap density budgets','tests':counts,'geometry_maximums':maximums,'selected_candidate':None,'history_status':history['status'],'reference_status':references['status'],'volume_generated':False,'fem_solved':False,'stage3_started':False,'policy_sha256':sha256(R/'acceptance_policy.json')})
lines=['Stage 1.6 completed.','','planar port contract:','    status = geometry PASS (overall stage FAIL)',f"    max rim displacement = {maximums['rim_displacement_m']} m",f"    max projected area error = {maximums['projected_area_relative_error']:.17e}",f"    max vector area error = {maximums['vector_area_relative_error']:.17e}",'','Stage 1 medium:','    cap triangles = 191','    tetra = 147569','    q<0.1 = 153','    cap-adjacent q<0.1 = 129','    P2 velocity DOF proxy = 801405']
for row in rows:
 lines += ['',row['candidate']+':',f"    cap triangles = {row['surface_qc']['total_cap_triangles']}",'    tetra = NOT RUN (surface density gate FAIL)','    P2 velocity DOF proxy = NOT RUN','    q<0.1 = NOT RUN (tetra); cap q_tri<0.1 = 0']
lines += ['','selected:','    candidate = NONE','    cap triangles = N/A','    tetra = N/A','    P2 velocity DOF proxy = N/A','    cap-adjacent q<0.1 = N/A','    total q<0.1 = N/A','    P1 = N/A','    P5 = N/A','    median = N/A','','tests:',f"    passed = {counts['passed']}",f"    failed = {counts['failures']+counts['errors']}",f"    skipped = {counts['skipped']}",'','report:','    reports/stage01_6/REPORT.md','','human review required:']
lines += ['    '+n for n in figures['images']]
lines += ['','STAGE 1.6 STATUS:','    FAIL','','STOP — Stage 3 not started.']
summary='\n'.join(lines)+'\n';(R/'terminal_summary.txt').write_text(summary);print(summary)
