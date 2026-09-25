# Step 1 血管几何合同包

本包固定实际旧 CFD 冻结网格对应的完整连续血管表面：`geometry/cfd_surface_axis_aligned_inlet_m.stl`。它含人工入口/出口延长段和四个封帽，已整体旋转到入口外法向约 +Z 的 CFD 坐标系。它是后续 HemoCell 几何测试应评审的候选；本轮没有运行 HemoCell。

STEP1_GEOMETRY_CONTRACT = PASS

来源绝对路径：`/home/lzy/projects/ulm_3D_vascular/outputs/cfd_flow/axis_aligned_inlet_geometry_anchor003274_20260829_111451/geometry/cfd_surface_axis_aligned_inlet_m.stl`。原文件与本包副本 SHA256 均为 `840da5e1c43ec31ac70ec781cbf75b32940bc73538b5dba83c58af1ddf5e36fb`。当前 CFD 配置虽仍写 s4 source_surface_run，真正被重用的 Base mesh 输入与containment记录指向本文件，并通过25项历史引用哈希核对。原始lumen和未旋转s4 STL都只是上游阶段，见 geometry_lineage.md。

单位是 **m**，不是凭数值猜测。s4 在 `vmtk_qc.py:tag_and_export_final_surface` 对 `meter_mesh.vertices` 乘1e-6；后续保存的刚体变换 scale=1，m/um矩阵单位明确。本包文件不能再乘1e-6。原始fMOST实验标定是否准确属于独立未验证项。

上游为fMOST标注SWC → 单连通分量/ROI003274 → Ultraliser lumen → s3边界包 → s4局部切口/TPS延长/重网格 →已有open表面继续封帽 → 42.264604756361436°整体刚体旋转。历史轴对齐导出脚本未在当前源码中找到；保存变换和产物/消费哈希均可查，不伪造执行记录。

| port / source identity | center in final CFD frame (m) | outward normal | area (m²) | origin/entity |
|---|---|---|---|---|
| inlet / cut_000 | [0.00010367370563330126, 5.7554233041039264e-05, 0.00015880639888138026] | [-7.007266199383469e-07, -6.027718016794421e-07, 0.9999999999995729] | 7.819753007034155e-12 | CUT_PORT / 4 |
| outlet_01 / cut_001 | [0.000114478874948153, 0.00015394596532480393, 0.00016617610422877752] | [0.018388404769920114, 0.8046101061009614, 0.5935186970350784] | 4.416140919927327e-12 | CUT_PORT / 3 |
| outlet_02 / cut_002 | [8.0800058383901e-05, 7.966215707512321e-05, 0.00011273836139137068] | [-0.9010100696103183, -0.2746666844492976, -0.3357663874101153] | 4.3150707517979834e-12 | CUT_PORT / 5 |
| outlet_03 / terminal_000 | [0.00017810767760060634, 0.0001383893240340422, 0.00012967037042376975] | [0.9913380590746697, 0.12599709605873874, -0.03706189977364442] | 5.110621093985814e-12 | TRUE_TERMINAL / 2 |


中心使用最终 cap 三角形的面积加权重心。boundary_contract.json 还保留原始s3切面中心、计划extension端点、s4非加权三角中心、s4 distal loop重心，彼此不可混用。normal为几何外法向，入口物理假定流向相反；浮点STL末位偏差不被强制抹成精确[0,0,1]。

1个入口cut_000（global edge814）；出口cut_001（edge897）、cut_002（edge1975）、terminal_000（global node4484/edge2073）。前三个是ROI截边，最后一个是源SWC结构叶。INLET/OUTLET来自ASSUMED角色，生理流向未确认。

真实来源的重建管腔与人工部分：CORE共45722面，其中45512个FAR_CORE面在s4局部操作中保留，210个核心侧collar面已经重网格；人工延长壁共21349面；CAP共191面。所有部分之后一起整体旋转。`metadata/region_contract.json`给出每端口extension归属；CORE/EXTENSION近接口标签含最近区域归属算法，不能称整段core都是原始实验表面未改。

已有拓扑QC及本轮轻量读取支持单连通、闭合、0边界边、0非流形边。历史逐面角度/aspect与配置阈值有差异，cap planarity硬阈值验收及最终人工评审未确认，自交未在本轮复算。详见 qc/existing_project_qc_summary.json 与 qc/independent_geometry_check.json。

旧Seeder实际入口使用矩形canoND数值平面，面积不能替代本表中的物理cap面积。旧1D和后期CFD流量目标也不同；均未映射到HemoCell。

HemoCell现有getFlagMatrixFromSTL按refDirN缩放到LU，并默认在X两端开口；当前几何是+Z入口和三出口，直接调用不具备已确认的适用条件。状态为 **NOT_READY_FOR_HEMOCELL_VOXELIZATION_TEST**，实际体素化未执行。

核查入口：[geometry_contract.json](metadata/geometry_contract.json)、[boundary_contract.json](metadata/boundary_contract.json)、[unit_contract.json](metadata/unit_contract.json)、[source_paths.json](metadata/source_paths.json)、[provenance.json](provenance/provenance.json)。详细中文报告另存于本包provenance/reports与上一级run目录。来源metadata原字节快照位于provenance/source_records，可离线核对；原STL保留原名，SHA256SUMS可核验包内全部文件。

未验证项：

- 历史axis-alignment STL导出器的当前源码文件/精确实现位置未找到；变换、产物、消费哈希存在。
- 原始fMOST采集spacing及半径物理标定、解剖轴/手性/世界坐标配准未独立证明。
- 端口生理流向、源SWC结构terminal是否对应真实解剖终端、样本实测流量未证明。
- 最终封帽表面的人工评审完成状态，以及planarity/angle/aspect等阈值差异如何用于未来验收尚未确认。
- HemoCell选轴、refDirN与物理dx合同、四端口体素化接口及坐标映射未实施；现有helper默认打开X两端不适用。
- 历史原始图像/SWC生成者与采集链完整身份未保存到可独立验证的本地证据。
- 本轮未复验自交、几何修复、体素化、流体求解；历史初始preflight失败与后续通过阶段分别保留。

NEXT_RECOMMENDED_ACTION = review Step 1 geometry contract manually
