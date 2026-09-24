"""Write Chinese reports from audited artifacts; contains no hydraulic solver."""
from pathlib import Path
import argparse
import json
import xml.etree.ElementTree as ET


def write_reports(out):
    data = out/'data'
    read = lambda name: json.loads((data/name).read_text())
    s, p, mapping, topology = map(read, ['final_summary.json', 'a_network_provenance.json', 'roi_ports_in_a.json', 'external_network_topology.json'])
    a = p['original_audit']; ports = s['roi_ports']
    port_table = '| Port | 原始 parent→child | t | 真实位置 / μm | 原半径 / μm |\n|---|---|---:|---|---:|\n'
    extension_table = '| Port | 到当前 FEM cap 轴向长 / μm | 当前 cap 等效半径 / μm | 20 截面全长度阻力估计 / Pa·s·m⁻³ |\n|---|---:|---:|---:|\n'
    for q in ports:
        port_table += f"| {q['name']} | {q['original_edge'][0]}→{q['original_edge'][1]} | {q['fraction_parent_to_child']:.12f} | {', '.join(f'{x:.6f}' for x in q['real_cut_xyz_um'])} | {q['radius_um']:.9f} |\n"
        extension_table += f"| {q['name']} | {q['fem_cap']['real_cut_to_cap_axial_um']:.6f} | {q['fem_cap']['equivalent_radius_um']:.6f} | {q['extension']['current_length_piecewise_linear_station_radius_R_Pa_s_m3']:.9e} |\n"
    tests = out/'logs/pytest.xml'
    test_info = '尚未读取测试记录'
    if tests.exists():
        suites = list(ET.parse(tests).getroot().iter('testsuite'))
        count = sum(int(x.get('tests', 0)) for x in suites); failures = sum(int(x.get('failures', 0))+int(x.get('errors', 0)) for x in suites)
        test_info = f'{count} 项测试，{failures} 项失败；见 logs/pytest.xml 和 logs/pytest.txt。'
        s['numerical_tests_status'] = 'PASS' if failures == 0 else 'FAIL'
        s['tests'] = dict(count=count, failures=failures, junit='logs/pytest.xml')
    (data/'final_summary.json').write_text(json.dumps(s, ensure_ascii=False, indent=2)+'\n')
    block = '结构根 2410 已精确找到，但现有记录只把它作为假定入口；没有找到将它标注为水力源点的数据集证据。外围端点的水力身份也未确认。按本轮源点 STOP 条件，未执行 A 的 operating-point 求解，也未创建或运行新 3D case。'
    principle = '这不是生理 ground truth；即使后续补齐边界定义，结果也只是 steady Newtonian、geometry-informed 的边界模型。'
    (out/'A_NETWORK_PROVENANCE_ZH.md').write_text(f'''# A-Network v1：来源审计

原始 A 到当前 FEM 的文件链已恢复，几何映射 PASS。水力源点识别未通过，阶段为 `{s['status']}`。

## 原始文件和两个“A”的区别

原始文件：`{s['original_A_path']}`。

SHA256：`{s['original_A_sha256']}`。

完整发布 SWC 有 **{a['nodes']} 节点、{a['edges']} 边、{a['components']} 连通分量、{a['cycle_rank']} 独立环**。
历史界面用于 ROI 的 analysis A 只取第 42 号分量：**7419 节点、7418 边、1 分量、0 环**，中心线总长 12013.420916 μm。其余 42 个分量保留为 reference，并不是错误或被修复掉的血管；全部组件列在 [a_components.csv](data/a_components.csv)。原文件有 43 个结构根，不能把 node 1 当作当前 A 的入口。

## 脚本、配置及实际加载路径

`s1-1_swc_roi_generate_mouse.py` → `utils.swc_roi_yaml_config` → `utils.rodent_vasculature.pipeline` → `catalog`/`swc_io`/`swc_analysis` → `utils.sampling.pipeline`/`roi_extraction`。

历史实际配置为 `{p['saved_config']}`：input_dir 指向原始数据目录、cohort=raw-analysis、sample_id=fMOST_0_5_6_0_0_6_0001_02_01。原始路径是从这个配置与运行的 preprocess_manifest 共同解析，并与 FEM source_contract 交叉核验，未凭文件名猜测。

**当前默认 YAML 已改为 `{p['current_default_input']}`，不能用重新运行默认脚本代替历史来源。** 当前源代码哈希和历史保存配置分开记录；当前辅助代码也不被假称为历史运行时逐字节快照。已读取两个入口的本地 import 闭包，52 个源文件及函数位置见 [source_helper_inventory.json](data/source_helper_inventory.json)。

原 xyz 是 voxel，实际换算是 **x×1、y×1、z×2 μm**；半径原值已为 μm，不再乘 z spacing。原图像为 192³，采样盒约定使用像素中心范围 `[0,191]×[0,191]×[0,382] μm`；体素外缘相差半个体素。原 SWC 物理 bbox 为 `{a['bbox_um']}`。历史导出的 analysis `.swc` 仍保存 voxel xyz；本轮以原始坐标转换并核对 NPZ `points_um`，没有重复乘 2。

`swc_analysis.select_analysis_swc` 按总中心线长度最大选分量 42。其原 node ID、parent、坐标和半径与 normalized NPZ 全部精确相同。图展示分支可能派生 1 μm 重采样曲线，但 ROI 生成重读原 normalized SWC，不消费派生曲线；本 ROI 输入未平滑/重采样。

## 裁剪、选择和重编号

seed=42，farthest_point 锚点，min_distance=45 μm，最多 80 候选；radius_plus_structure、robust scaler、kmeans k=5，coverage_balanced 选择。当前代表 anchor=3274，cluster=0，selection_rank=1。

中心 `{p['roi_center_um']}` μm；bbox `{p['roi_bbox_um']}` μm；尺寸 80×80×120 μm。先对原边作精确盒裁剪，再只保留含锚点的连通片。历史记录盒内有 12 个局部分量；最终保留 106 个原节点、3 个虚拟切点、108 段边（总长 189.264567 μm）。

原节点是 **3225–3307、4462–4484**；analysis edge IDs 为 **814–897、1975、2051–2073**。完整原边 `(u,v)`、raw edge index、裁剪起止 t 和局部 ID 在 [roi_edge_provenance.csv](data/roi_edge_provenance.csv)。

ROI NPZ 保留 original IDs；canonical `roi_core.swc` 重新编号为 1–109，SWC node 1 是 local node 106 的入口虚拟切点，绝不是原 A node 1。所有对应关系见 [roi_node_provenance.csv](data/roi_node_provenance.csv)。

{port_table}

O3 是 node 4484（原边终点 t=1），不是盒裁剪点。其 `TRUE_TERMINAL` 只表示原图 degree=1。三个中间切点已在新网络副本中插入虚拟节点 -10001、-10002、-10003；原 SWC 未改动。完整 SI 图副本为 [a_graph_with_exact_cuts_si.npz](data/a_graph_with_exact_cuts_si.npz)。

## 全部 graph crossing

对完整原始文件逐边扫描：盒面交点共 **35**，其中 analysis 分量 **26**。只有 **3** 个属于实际保存 ROI，另有 O3 原端点，合计 **4** 个 FEM 端口。其他盒面交点属于没有选入当前 ROI 的血管片段；它们不构成当前 ROI 的遗漏端口。也检查了所有保留原节点的 incident edges，没有额外连接。

若把“ROI”另定义为盒内所有血管，35 个 cut crossings 就必须全部处理，不能继续套用当前四端口 FEM。当前模型采用历史保存的连通 ROI。逐条记录见 [all_A_roi_bbox_crossings.csv](data/all_A_roi_bbox_crossings.csv)。

## STL 与人工延长段

`s2_swc_stl_model_generate.py` → saved ROI NPZ → canonical SWC → vascular H5 → Ultraliser。原 ROI SWC 半径未变；**H5 feed radius=原半径×0.91**，H5 第四列写 diameter=2×feed radius，坐标 float32 量化。Ultraliser：6 voxel/μm、polylines-with-spheres、DMC、solid、adaptive optimization 5、smooth 5、Laplacian 10；输出 `{p['lineage']['reconstruction_run']}/geometry/lumen_surface_um.stl`。

随后局部端口切开，VMTK boundarynormal + thinplatespline 人工延长，局部 collar 重网格并封 distal cap；核心来源和 cap 身份由最终 surface contract 固定。真实 cut 与人工 FEM cap 分开记录：

{extension_table}

这里阻力使用历史 20 个截面等效圆半径，补齐两端半 bin 后按分段线性半径精确积分；它是几何估计。旧代码 `utils/cfd_surface_prepare/vmtk_qc.py:973` 对 t=0.025…0.975 做 trapz，覆盖 95% 长度，因此旧表与本轮全长度估计相差约 5.3%。旧代码及结果未修改。

## 当前 FEM 的实际网格链

最终 VMTK tagged surface → FEM Stage 0 contract → Stage 1.6 planar port/Stage 1.7 cap 内部三角化优化 → `inputs/fem_reference/exterior_surface.npz`（SHA 31746b2c…）→ `sv.modeling.PolyData` → `sv.meshing.TetGen` 同时 surface+volume meshing → `outputs/sv1/SV_MESH` → frozen_reference → mean-2p0-mmps 原样复制网格。

SV 生成尺寸为 3.917947590640251e-7 m，当前 **70363 nodes、371402 tetra**，网格 SHA `{p['actual_sv_mesh_manifest']['sha256']}`。SV surface remesh 相对旧表面有可测变化：壁面双向采样距离 max=0.066919 μm、P95=0.013988 μm，包围体积差 0.676911%；不能说从最初 STL 到当前 FEM 完全没有几何变化。当前四个 cap 面积和重心已直接从 production VTP 重算。

旧 20 截面资料早于 SV surface remesh，不能被视为当前 3D 延长段的精确阻力。后续若开放 BC 转移，应在当前冻结网格上重新取截面；本轮因源点条件未通过，没有输出可用压力 BC。

来源文件、配置、网格和生产流场哈希见 [protected_input_hashes.json](data/protected_input_hashes.json)，导入链及 SV 几何差异完整值见 [a_network_provenance.json](data/a_network_provenance.json)。本轮未修改上述任何源文件。

## 数据集与根/终端证据

数据集论文描述的是切块形态重建与半径/中心线标注，体素尺度为 1×1×2 μm，并讨论上游标注阶段的平滑与采样。[原论文](https://www.frontiersin.org/journals/neuroinformatics/articles/10.3389/fninf.2026.1809341/full)、[发布说明](https://zenodo.org/records/19184697)。本地记录未找到当前 sample 的动静脉源汇/水力端点标签。

{block} 2410 的物理位置为 (60.95,124.99,108.96) μm，距图像采样中心边界最近 60.95 μm，处于内部。单凭 parent=-1、半径较大或过去的 assumed inlet 都不能在本轮将其升级成已确认水力源点。

共有 231 个 degree-1 节点（analysis A 为 123 个），保守保留 UNKNOWN_TERMINAL。另列几何证据：122 个端点半径球触及图像中心范围边界，31 个在一体素距离内，31 个为深内部候选。它们是可复核的几何筛选，不是已验证的 SOURCE_ROOT/OUTER_BOUNDARY_TERMINAL 分类。节点 4484 距边界 26.99 μm，也不能自动设 p=0。
''')
    (out/'ROI_NETWORK_BOUNDARY_MODEL_ZH.md').write_text(f'''# ROI 外部网络模型：目前只完成拓扑审计

{block}

已用原始半径构建全部 9785 条边的 SI 阻力：

`R_e = (8 μ L / π) (r0²+r0·r1+r1²) / (3 r0³ r1³)`，μ=0.00345312 Pa·s。

这是 r(s) 线性变化的精确积分，r0=r1 时连续退化为 `8μL/(πr⁴)`，不使用平均半径。半径正、长度正、没有零长边；没有删环或强制树化。原 SWC 一节点单 parent 的格式无法表达任意多父连接，观测 0 环不证明生理网络没有吻合环。

删除 ROI 内部 108 段真实边后，保留三个精确切点及 O3。包含四个端口的外部组件互不连接，外部 degree 分别为 **1、1、1、0**，见 [external_network_topology.json](data/external_network_topology.json)。盒内未选中血管仍属于外部网络，人工 extension 从未加入真实 A 图。

这说明拓扑上没有三个 outlet 之间的外部连接。要把 O1/O2 化为有限 R，仍需要确定各自远端如何接到参考压力。O3 外部没有边：若将其视为零压 reservoir，就是固定压力约束/R=0 的理想化；若视为未解析的盲端，就是 no-flow，不能输出有限 R3。两种含义不可混淆。

通用 sparse solver 与 Schur 工具已实现并通过解析测试。完整图只建立 scipy CSR matrix，没有 dense N×N。数值解、条件数、节点守恒残差目前均未计算，见 network_solver_audit.json。未把拓扑无连接画成已验证的 conductance heatmap。

合法的 exterior reduction 应保留 ports，内部 sparse 消元：`Y=G_BB−G_BI G_II⁻¹ G_IB`。源点若是非零固定压力，ROI 端口关系一般是 **q_B=Y p_B+f_source**；只有把源点也保留为边界变量，或正确处理齐次参考后，才能写齐次映射。通用测试已覆盖 affine forcing、对称性、被动性、full/reduced 等价和奇异 Z 不可逆情形。

本案例 Y、Z、R1/R2/R3、coupling magnitude 均为 null。CSV 只有 schema，无伪造数字。{principle}
''')
    (out/'ROI_FIXED_PRESSURE_BC_ZH.md').write_text(f'''# 固定压力 BC：BLOCKED_DO_NOT_USE

{block}

**本轮没有可用于 FEM 的 O1/O2/O3 压力建议值。** JSON 中 raw/shifted pressure、pressure differences、source solution SHA 都是 null；mapping SHA 已保存。旧 preprocessing 的压力不能冒充本轮解，它使用未经本轮确认的 structural-root/all-leaves 假设和不同 operating point。

若后续源点和终端条件明确，先设源点单位压差，计算有符号真实 ROI inlet cut 流量 `Q_unit`，确认 `Q_unit>0`，再用 `λ=1.551359160440232e-14/Q_unit` 缩放整网压力和流量。目标是 ROI inlet cut，不是 A source total flow；不允许取 abs 掩盖反向流。

全网所得 p 位于真实 cut。人工延长段应按 `p_cap=p_real−R_ext Q_outward` 转到 FEM cap；入口的 Q_outward<0，符号自然反转。这是冻结 operating point 的估计，3D 求解后实际分流变动会影响修正精度。应先按当前网格复核 extension profile，再作 gauge shift。

可选择所有 cap 压力减去同一个最小值；每一对压力差必须保持不变。0 Pa 仅为参考压力。gauge invariance、extension 符号和精确变半径阻力均有永久测试。

{extension_table}

表内数值仅为阻力估计，不是压力建议值。source/terminal gate 未通过，因此没有复制 mean-2p0-mmps，也没有 config diff、新稳态结果或三组分流比较图。
''')
    (out/'NEXT_MODEL_EXTENSIONS_ZH.md').write_text('''# 后续模型扩展（本轮未实现）

先解决源点和终端身份，再讨论 Pries–Secomb apparent viscosity、Fåhræus–Lindqvist 效应，以及 hematocrit/phase separation。它们分别改变有效黏度、管径相关阻力和分叉血细胞分配，会使半径与流量反馈耦合，需要新的模型参数、适用范围和独立测试。

不能把这些模型和 outlet boundary 环境同时改变后，将分流差异全部归因于 ROI 裁剪。本轮一直固定 μ=0.00345312 Pa·s、ρ=1056 kg/m³，没有实现任何血细胞比容输运、非牛顿黏度或微泡模型。

±5%/±10% 是声明的测量敏感性幅度，不是从本 sample 实测推得的误差分布。没有找到可以用于逐节点/逐分支误差模型的重复半径测量或校准不确定度，因此标记 NO_EMPIRICAL_RADIUS_UNCERTAINTY_AVAILABLE。
''')
    global_rows = ''.join(f"| {x['radius_factor']:.2f} | {x['exact_R_factor']:.6f} |\n" for x in s['radius_sensitivity']['global_scaling'])
    (out/'A_NETWORK_1D0D_REVIEW_ZH.md').write_text(f'''# 一句话结论

**原始 A 图和 ROI 精确来源已恢复，当前阶段是 `{s['status']}`。** {block} 因而还不能回答 O2=85.2% 是否主要由 equal-pressure 截断造成，也不能提供网络分流。

# 完整 A 是怎么变成水管网络的

把每段原 SWC 边视为一根有阻力的水管。半径恒定时 `R=8μL/(πr⁴)`：半径稍变，阻力就会明显变。本轮逐边保留线性半径变化并用稳定解析积分；所有计算内部采用 m、Pa、m³/s。

完整发布文件有 9828 节点/9785 边/43 分量/0 环；历史 analysis A 有 7419/7418/1/0。没有孤立节点、重复 ID、重复坐标、零长边或非正半径。度数分布 1:231、2:9454、3:141、4:2；共 143 个分支节点。最短边 {a['min_edge_length_um']:.6f} μm，最小半径 1 μm。原始单 parent 格式的 0 环不等同生理无环。

# ROI 在完整 A 的哪里

![完整 A 与 ROI](figures/01_full_A_roi_provenance.png)

锚点 3274，中心 (130.04,82.04,87.18) μm，80×80×120 μm 盒内只保留锚点连通片。图中星号是**结构根**，没有标成已确认的水力源点。

{port_table}

盒内其他血管共有更多 crossings，但保存的 ROI 只有 3 个 cut 和 O3 原端点，没有第 5 个漏报端口。三个人工切点已精确插入网络副本。

# 完整 A 的边界条件是什么

提出的 outer terminal p=0 只是 gauge/reference。原文件的 43 个结构根、选中分量的唯一结构根 2410 都可以确定；**谁是水力源点尚未确定**。旧脚本显式用 assumed inlet 和 all structural leaves=0，这些旧假设不满足本轮更严格的源点/终端证据要求。

231 个原图端点暂留 UNKNOWN_TERMINAL，保留图像边界距离、球半径是否接边、是否在一体素以内等筛选证据。不是把所有端点叫 physiological outlet，也不是宣称没有真实 outer terminal；当前缺少确认它们身份的数据。

# 为什么要把 ROI inlet 调到 2.0 mm/s

为使网络与原 3D 比较同一个 operating point，需匹配真实 ROI inlet cut 的 Q={s['roi_target_Q']:.15e} m³/s，即约 15.513592 pL/s。单位源压解后按 `λ=Q_target/Q_ROI_unit` 缩放；它不是要求整个 A 的入口流量等于 ROI 流量。该逻辑已通过解析测试，尚未用于 A 的未确定边界。

# 完整 A 给出的 ROI pressure 和 flow

| Port | 压力 Pa | 流量 m³/s | 分流比例 |
|---|---|---|---|
| INLET | 未计算 | 未计算（目标见上） | 未计算 |
| O1 | 未计算 | 未计算 | 未计算 |
| O2 | 未计算 | 未计算 | 未计算 |
| O3 | 未计算 | 未计算 | 未计算 |

未计算不是 0；A source pressure/total inflow 同样为 null。

# 当前 zero-pressure ROI 和完整 A 有多大区别

当前 3D 的已知比例仍是 **4.257917% / 85.205026% / 10.537057%**。没有合法 A baseline，不能计算百分点差；Figure 04 未生成。压力图/流量图 Figure 02、03 也不生成；清单见 data/uncomputed_outputs.json。

# O2 为什么多或者为什么不再多

目前无法判定。不能把原比例归因于 geometry、topology 或 zero-pressure BC 中的任何单一因素，也没有预设 O2 必须下降。

# 三个 outlet 能不能用独立 resistance

移除真实 ROI 内部段后，端口外部组件彼此分离；没有外部 outlet-to-outlet 拓扑路径。O1/O2 的远端压力条件仍未确定；O3 外部 degree=0，**没有数据支持一个有限 R3**。Y/Z 数值和 normalized coupling 尚不可用，Figure 05 未生成；原理和拓扑证据见 [ROI_NETWORK_BOUNDARY_MODEL_ZH.md](ROI_NETWORK_BOUNDARY_MODEL_ZH.md)。

# fixed-pressure BC 是什么

本轮 O1/O2/O3 均为 **未定义，禁止使用**。精确 mapping SHA 和 extension 估计已就绪，压力差、gauge shift 结果和 solution SHA 都为空。见 [ROI_FIXED_PRESSURE_BC_ZH.md](ROI_FIXED_PRESSURE_BC_ZH.md)。

# 如果跑了新 3D FEM

未创建、未运行，原因是前置源点/终端条件失败；3D validation split=null，1D–3D 差异不可计算，Figure 06 不生成。没有修改原 production case、几何、网格或 Particle。

# A 自己也是截断数据意味着什么

A 来自一个图像子块，原图端点可能是图像截断或重建终止，并非实际微循环的生理终末。原始数据的 43 个组件有历史来源解释，不是本轮求解失败后偷偷删掉的分量；其中 42 个不参与原 ROI 生成。O3=4484 位于图像内部，尤其不能自动接到公共 p=0 reservoir。

# radius uncertainty

**NO_EMPIRICAL_RADIUS_UNCERTAINTY_AVAILABLE**。已计算几何阻力敏感性，并保存确定性的分支正负 5% 原节点半径图案，没有随机抽一次。分叉共享节点使用最小 incident branch ID，保证每个原节点只有一个半径。source 未识别，所以没有声称完成 full-A 的压力/分流敏感性求解。

| 全局半径倍率 | 精确阻力倍率 |
|---:|---:|
{global_rows}

如果边界定义固定、所有半径统一缩放，并重新匹配 ROI inlet Q，线性模型的压力尺度会乘上述因子而分流比例不变；这是解析性质，不是本案例已计算结果。分支 ±5% 的阻力变化范围约 0.822702–1.227738；水力分流敏感性仍被 source gate 阻断。

# terminal uncertainty

MODEL_A 的可信 outer 集合、MODEL_B 的严格 outer 集合都无法认证。若将所有 uncertain 端点设 no-flow，将没有有效 sink；无法匹配正 ROI inlet target。两个模型均明确 invalid/not solved，没有捏造终端敏感性分流。几何候选数量见 data/terminal_sensitivity.json，不能与可信终端数量混用。

# 这是不是生理 ground truth

{principle} 另有圆管近似、原始半径/截断拓扑、原 SWC-to-STL radius feed 缩放、SV surface remesh、人工延长段估计等限制。未来即便通过边界 gate，1D 与 3D 的曲率、分叉、非圆截面与入口发展也可能导致真实差异。

# 测试、复现与耗时

{test_info}

本地 WSL、1 worker、4 项线程环境变量均为 1；源图/映射/阻力审计 {s['runtime_seconds']:.3f} s，peak RSS {s['peak_memory_MiB']:.1f} MiB（这是脚本耗时，不含人工式来源调查、测试和绘图）。未使用 Vast/GPU。绘图耗时另见 logs/render_metrics.json。

新增科学计算包 `network_1d0d/`，执行/绘图/报告脚本分别在 `scripts/`，永久测试在 `tests/network_1d0d/`。复现命令见 README.md。

Git 分支 `dev/a-network-1d0d-boundary-v1`，起始 HEAD `fd21a850a16d0ba17ef3d521123badd53864a3a2`。进入时已有 1133 个 tracked 改动，本轮没有覆盖；完整前后 git status/diff --stat 在 logs/，本轮新增文件统计单列。

最后检查时检测到已有 `utils/sampling/sampling_io.py` 和 `utils/swc_roi_yaml_config.py` 又发生变化（本轮未编辑它们），当前 tracked diff 增为 1134 个文件。已保存 source helper 审计时/当前哈希，见 logs/concurrent_worktree_changes.json；没有覆盖这些期间变化。40 个受保护的原始/配置/几何/生产输入文件仍全部哈希一致。来源追踪的当前 helper inventory 是审计时快照，不冒称最终工作树逐字节快照。

# 下一步

1. 补充当前 sample 的源点/外围端点标注，或明确声明以结构根 2410 作为理想化源点并选择终端闭合策略；后者属于新的模型假设，不是恢复出来的实测事实。
2. 条件明确后运行单位压差、真实 ROI cut flow scaling、终端/半径敏感性及 Schur audit；复核当前网格的 extension 截面。
3. 全部 gate 通过后，才创建新固定压力 3D case，并如实比较分流。
''')
    (out/'README.md').write_text('''# A-Network Steady 1D/0D Boundary Model v1

本轮结论与结果：[中文主报告](A_NETWORK_1D0D_REVIEW_ZH.md)。状态 A_NETWORK_SOURCE_AMBIGUOUS；已完成几何来源/拓扑/精确端口/阻力审计与通用库解析测试，未生成 A 水力解或 3D BC。

在项目根目录使用已具备 numpy/scipy/pyvista/matplotlib/PyYAML/pytest 的 Python。此次复用只读的 `/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python`，没有修改 Particle 代码或环境。

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python scripts/run_a_network_1d0d_boundary_v1.py
python -m pytest -q tests/network_1d0d --junitxml=reports/a_network_1d0d_boundary_v1/logs/pytest.xml
python scripts/render_a_network_1d0d_boundary_v1.py
python scripts/report_a_network_1d0d_boundary_v1.py
```

审计入口针对当前被冻结来源，默认输出仅写本报告目录；可通过 --output 改为新目录。它不运行旧 s1/s2，不写已有生产 case，也不根据旧 assumed-inlet 资料自动越过 source gate。

所有 null 表示未计算，空 hydraulic CSV 只定义 schema。Figure 01 有 PNG/PDF；Figure 02–06 因前置条件失败不生成。data/uncomputed_outputs.json 逐项说明。日后源点确认后仍需添加经过审阅的 case boundary 声明与完整求解驱动，不能仅改 summary 状态。
''')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1]/'reports/a_network_1d0d_boundary_v1')
    write_reports(parser.parse_args().output)
