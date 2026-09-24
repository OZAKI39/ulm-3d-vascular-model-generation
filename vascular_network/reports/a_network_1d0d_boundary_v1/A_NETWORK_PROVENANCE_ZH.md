# A-Network v1：来源审计

原始 A 到当前 FEM 的文件链已恢复，几何映射 PASS。水力源点识别未通过，阶段为 `A_NETWORK_SOURCE_AMBIGUOUS`。

## 原始文件和两个“A”的区别

原始文件：`/home/lzy/projects/ulm_3D_vascular/vessel_model/T - A high-resolution dataset of mouse brain vasculature/raw_data/analysis_data/analysis_data/swc/fMOST_0_5_6_0_0_6_0001_02_01.swc`。

SHA256：`4d244a5ed73abb820e1c8014a27e58389f132bd1f47eb18e91a55a380dcd972a`。

完整发布 SWC 有 **9828 节点、9785 边、43 连通分量、0 独立环**。
历史界面用于 ROI 的 analysis A 只取第 42 号分量：**7419 节点、7418 边、1 分量、0 环**，中心线总长 12013.420916 μm。其余 42 个分量保留为 reference，并不是错误或被修复掉的血管；全部组件列在 [a_components.csv](data/a_components.csv)。原文件有 43 个结构根，不能把 node 1 当作当前 A 的入口。

## 脚本、配置及实际加载路径

`s1-1_swc_roi_generate_mouse.py` → `utils.swc_roi_yaml_config` → `utils.rodent_vasculature.pipeline` → `catalog`/`swc_io`/`swc_analysis` → `utils.sampling.pipeline`/`roi_extraction`。

历史实际配置为 `/home/lzy/projects/ulm_3D_vascular/outputs/rodent_vasculature/all_run_20260825_133152/source_swc_roi_generate.yaml`：input_dir 指向原始数据目录、cohort=raw-analysis、sample_id=fMOST_0_5_6_0_0_6_0001_02_01。原始路径是从这个配置与运行的 preprocess_manifest 共同解析，并与 FEM source_contract 交叉核验，未凭文件名猜测。

**当前默认 YAML 已改为 `vessel_model/T - A high-resolution dataset of mouse brain vasculature_vmd_branch_aic`，不能用重新运行默认脚本代替历史来源。** 当前源代码哈希和历史保存配置分开记录；当前辅助代码也不被假称为历史运行时逐字节快照。已读取两个入口的本地 import 闭包，52 个源文件及函数位置见 [source_helper_inventory.json](data/source_helper_inventory.json)。

原 xyz 是 voxel，实际换算是 **x×1、y×1、z×2 μm**；半径原值已为 μm，不再乘 z spacing。原图像为 192³，采样盒约定使用像素中心范围 `[0,191]×[0,191]×[0,382] μm`；体素外缘相差半个体素。原 SWC 物理 bbox 为 `[[1.0, 0.58, 0.98], [190.23, 190.4, 381.1]]`。历史导出的 analysis `.swc` 仍保存 voxel xyz；本轮以原始坐标转换并核对 NPZ `points_um`，没有重复乘 2。

`swc_analysis.select_analysis_swc` 按总中心线长度最大选分量 42。其原 node ID、parent、坐标和半径与 normalized NPZ 全部精确相同。图展示分支可能派生 1 μm 重采样曲线，但 ROI 生成重读原 normalized SWC，不消费派生曲线；本 ROI 输入未平滑/重采样。

## 裁剪、选择和重编号

seed=42，farthest_point 锚点，min_distance=45 μm，最多 80 候选；radius_plus_structure、robust scaler、kmeans k=5，coverage_balanced 选择。当前代表 anchor=3274，cluster=0，selection_rank=1。

中心 `[130.04, 82.04, 87.18]` μm；bbox `[[90.03999999999999, 42.040000000000006, 27.180000000000007], [170.04, 122.04, 147.18]]` μm；尺寸 80×80×120 μm。先对原边作精确盒裁剪，再只保留含锚点的连通片。历史记录盒内有 12 个局部分量；最终保留 106 个原节点、3 个虚拟切点、108 段边（总长 189.264567 μm）。

原节点是 **3225–3307、4462–4484**；analysis edge IDs 为 **814–897、1975、2051–2073**。完整原边 `(u,v)`、raw edge index、裁剪起止 t 和局部 ID 在 [roi_edge_provenance.csv](data/roi_edge_provenance.csv)。

ROI NPZ 保留 original IDs；canonical `roi_core.swc` 重新编号为 1–109，SWC node 1 是 local node 106 的入口虚拟切点，绝不是原 A node 1。所有对应关系见 [roi_node_provenance.csv](data/roi_node_provenance.csv)。

| Port | 原始 parent→child | t | 真实位置 / μm | 原半径 / μm |
|---|---|---:|---|---:|
| INLET | 3224→3225 | 0.034013605442 | 103.266599, 46.999184, 147.180000 | 1.587390476 |
| O1 | 3307→3308 | 0.371794871795 | 113.417436, 122.040000, 100.410256 | 1.173538462 |
| O2 | 3238→4386 | 0.890909090909 | 90.040000, 48.091273, 110.910909 | 1.152856364 |
| O3 | 4483→4484 | 1.000000000000 | 164.010000, 96.360000, 82.780000 | 1.255000000 |


O3 是 node 4484（原边终点 t=1），不是盒裁剪点。其 `TRUE_TERMINAL` 只表示原图 degree=1。三个中间切点已在新网络副本中插入虚拟节点 -10001、-10002、-10003；原 SWC 未改动。完整 SI 图副本为 [a_graph_with_exact_cuts_si.npz](data/a_graph_with_exact_cuts_si.npz)。

## 全部 graph crossing

对完整原始文件逐边扫描：盒面交点共 **35**，其中 analysis 分量 **26**。只有 **3** 个属于实际保存 ROI，另有 O3 原端点，合计 **4** 个 FEM 端口。其他盒面交点属于没有选入当前 ROI 的血管片段；它们不构成当前 ROI 的遗漏端口。也检查了所有保留原节点的 incident edges，没有额外连接。

若把“ROI”另定义为盒内所有血管，35 个 cut crossings 就必须全部处理，不能继续套用当前四端口 FEM。当前模型采用历史保存的连通 ROI。逐条记录见 [all_A_roi_bbox_crossings.csv](data/all_A_roi_bbox_crossings.csv)。

## STL 与人工延长段

`s2_swc_stl_model_generate.py` → saved ROI NPZ → canonical SWC → vascular H5 → Ultraliser。原 ROI SWC 半径未变；**H5 feed radius=原半径×0.91**，H5 第四列写 diameter=2×feed radius，坐标 float32 量化。Ultraliser：6 voxel/μm、polylines-with-spheres、DMC、solid、adaptive optimization 5、smooth 5、Laplacian 10；输出 `/home/lzy/projects/ulm_3D_vascular/outputs/model_generate/ultraliser_anchor003274_20260825_133350/geometry/lumen_surface_um.stl`。

随后局部端口切开，VMTK boundarynormal + thinplatespline 人工延长，局部 collar 重网格并封 distal cap；核心来源和 cap 身份由最终 surface contract 固定。真实 cut 与人工 FEM cap 分开记录：

| Port | 到当前 FEM cap 轴向长 / μm | 当前 cap 等效半径 / μm | 20 截面全长度阻力估计 / Pa·s·m⁻³ |
|---|---:|---:|---:|
| INLET | 15.708205 | 1.571326 | 2.242282455e+16 |
| O1 | 11.797444 | 1.177505 | 5.281376695e+16 |
| O2 | 11.696299 | 1.161572 | 5.470624510e+16 |
| O3 | 12.692440 | 1.267740 | 4.249596780e+16 |


这里阻力使用历史 20 个截面等效圆半径，补齐两端半 bin 后按分段线性半径精确积分；它是几何估计。旧代码 `utils/cfd_surface_prepare/vmtk_qc.py:973` 对 t=0.025…0.975 做 trapz，覆盖 95% 长度，因此旧表与本轮全长度估计相差约 5.3%。旧代码及结果未修改。

## 当前 FEM 的实际网格链

最终 VMTK tagged surface → FEM Stage 0 contract → Stage 1.6 planar port/Stage 1.7 cap 内部三角化优化 → `inputs/fem_reference/exterior_surface.npz`（SHA 31746b2c…）→ `sv.modeling.PolyData` → `sv.meshing.TetGen` 同时 surface+volume meshing → `outputs/sv1/SV_MESH` → frozen_reference → mean-2p0-mmps 原样复制网格。

SV 生成尺寸为 3.917947590640251e-7 m，当前 **70363 nodes、371402 tetra**，网格 SHA `1a7a69dcba52f0475b2280775921e4dedcdfefa67b7965096d86cb9dba38bdf9`。SV surface remesh 相对旧表面有可测变化：壁面双向采样距离 max=0.066919 μm、P95=0.013988 μm，包围体积差 0.676911%；不能说从最初 STL 到当前 FEM 完全没有几何变化。当前四个 cap 面积和重心已直接从 production VTP 重算。

旧 20 截面资料早于 SV surface remesh，不能被视为当前 3D 延长段的精确阻力。后续若开放 BC 转移，应在当前冻结网格上重新取截面；本轮因源点条件未通过，没有输出可用压力 BC。

来源文件、配置、网格和生产流场哈希见 [protected_input_hashes.json](data/protected_input_hashes.json)，导入链及 SV 几何差异完整值见 [a_network_provenance.json](data/a_network_provenance.json)。本轮未修改上述任何源文件。

## 数据集与根/终端证据

数据集论文描述的是切块形态重建与半径/中心线标注，体素尺度为 1×1×2 μm，并讨论上游标注阶段的平滑与采样。[原论文](https://www.frontiersin.org/journals/neuroinformatics/articles/10.3389/fninf.2026.1809341/full)、[发布说明](https://zenodo.org/records/19184697)。本地记录未找到当前 sample 的动静脉源汇/水力端点标签。

结构根 2410 已精确找到，但现有记录只把它作为假定入口；没有找到将它标注为水力源点的数据集证据。外围端点的水力身份也未确认。按本轮源点 STOP 条件，未执行 A 的 operating-point 求解，也未创建或运行新 3D case。 2410 的物理位置为 (60.95,124.99,108.96) μm，距图像采样中心边界最近 60.95 μm，处于内部。单凭 parent=-1、半径较大或过去的 assumed inlet 都不能在本轮将其升级成已确认水力源点。

共有 231 个 degree-1 节点（analysis A 为 123 个），保守保留 UNKNOWN_TERMINAL。另列几何证据：122 个端点半径球触及图像中心范围边界，31 个在一体素距离内，31 个为深内部候选。它们是可复核的几何筛选，不是已验证的 SOURCE_ROOT/OUTER_BOUNDARY_TERMINAL 分类。节点 4484 距边界 26.99 μm，也不能自动设 p=0。
