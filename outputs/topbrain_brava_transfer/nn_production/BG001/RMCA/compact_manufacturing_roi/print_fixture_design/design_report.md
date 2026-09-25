# BALANCED 血管芯与浇注盒设计审核

## 1. 这次做了什么

本轮围绕已经验收的血管芯设计了一个顶部开放的浇注盒。入口和出口的位置来自原有中心线及其编号对应表，延长段沿端点几何方向起步，再通向侧壁。原血管与延长段构成牺牲芯，盒体单独保存为另一个部件。图片用于人工检查位置、弯曲和装配关系，本轮没有切片或提交打印。当前仅完成 3/4 个端口延长，其余端口因间距约束被保留为待调整。

## 2. 输入是什么

- 已验收血管 STL：`/home/lzy/projects/ulm_3D_vascular/outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi/BG001_RMCA_BALANCED_print_candidate.stl`。
- 紧凑血管中心线：`/home/lzy/projects/ulm_3D_vascular/outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi/candidates/BALANCED/compensated.swc`。
- 生成最终表面的拟合中心线：`/home/lzy/projects/ulm_3D_vascular/outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi/candidates/BALANCED/VascularMD/compensated_vmd_smooth.swc`。
- 来源清单：`/home/lzy/projects/ulm_3D_vascular/outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi/compact_manifest.json`。
- 固定打印坐标变换：`/home/lzy/projects/ulm_3D_vascular/outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi/print_transform.json`。

SWC 用来确认真实端点、父子关系、原分支编号和半径；STL 用来计算盒体尺寸、碰撞与最终实体。本轮不重新识别 MeVO，也不重新运行配准、NN、半径补偿、VascularMD 建模或朝向搜索。端口使用拟合中心线中与原封口记录一致的实际半径，早期补偿值和原始值一并保存在端点表中。

![处理前：原始血管芯](QC/01_source_vascular_core.png)

## 3. 盒子多大

- 内部：**78.29 × 81.85 × 56.39 mm**（X × Y × Z）。
- 外部：**84.29 × 87.85 × 59.39 mm**。
- X、Y 两侧分别留 12、12 mm；原血管下方留 10 mm，上方留 12 mm。
- 侧壁厚 3 mm，底板厚 3 mm，顶部无盖。

“顶部开放”指盒内空腔可从上方进入；盒壁和底板作为一个有厚度的实体，其三角网格仍应闭合，不能有破面。沿用原打印坐标，因此盒底可能位于 Z=0 以下；这不是重新选择打印姿态，后续切片时再分别定位两个部件。

## 4. 有几个入口和出口

来源拓扑确认 **1 个入口、3 个出口**，未根据 STL 外观猜测身份。本次仅有 **3/4 个端口**的延长通过检查；端点没有被删除。

I 表示入口，O 表示出口。端点箭头是中心线向模型外侧的几何方向，不是实测血流方向。

## 5. 每个端口在哪里

| 端口 | 原 SWC 编号 / 分支 | 盒壁 | 路线 | 延长长度 mm | 半径 mm | 最小弯曲半径 mm |
|---|---|---|---|---:|---:|---:|
| I1（入口） | 2039 / 137_138 | -Y | 平滑曲线 | 26.39 | 0.6026 | 6.20 |
| O1（出口） | 2081 / 138_139 | +Y | 平滑曲线 | 41.84 | 0.5851 | 6.19 |
| O2（出口） | 2123 / 140_141 | -X | 平滑曲线 | 32.55 | 0.6263 | 6.19 |
| O3（出口） | 2127 / 140_142 | 未找到可接受布局 | 待调整 | — | 0.6311 | — |

盒壁名称属于当前打印坐标，与左右脑解剖方向无关。每段长度从原端点算到最外端，不包含向原封口内侧的布尔连接重叠。

端口在盒外继续延伸 10 mm。穿墙孔的余量为**径向 0.2 mm**，即孔直径比端口直径大 0.4 mm；孔不是软管接头。
同壁孔中心间距至少 8 mm，且不少于两半径加两倍径向余量；孔中心距顶部、底部和侧角至少 8 mm。

![处理后：血管芯与通过检查的端口](QC/06_core_with_ports.png)

## 6. 有没有风险

- 穿墙孔留有径向 0.20 mm 装配余量，浇注前需要确认定位和防漏方式。
- 当前固定打印坐标保持不变；两个部件尚未进行新一轮切片或打印误差验证。
- 仅用中心线几何方向布置端口，不据此推断生理血流。
- 入口/出口中最小盒壁边距为 8.021 mm，接近 8 mm 下限，人工审核时需留意实际打印公差。
- ENDPOINT_CLEARANCE_CONFLICT: O3 起始截面距分支 140_141 仅 1.8647 mm，小于要求 3 mm；改换盒壁不能解决。
- PORT_LAYOUT_NEEDS_MANUAL_REVIEW: An endpoint has no feasible route
- INCOMPLETE_ENDPOINT_EXTENSION
- 最大端点方向转角：90.23°。该值是起始切线与盒壁法线的总方向差，并非接口处的尖角。
- 路线至非连接区血管的最小表面间距：4.094 mm。
- 任意两个端口的最小表面间距：38.580 mm。
- 同壁端口最小中心间距：不适用 mm。

![起点间距冲突](QC/11_endpoint_clearance_conflict.png)

橙色圆盘是端口必须包含的起始截面，红线表示它到旁边血管的最短距离。由于这段距离在延长起点就已经不足，改变目标盒壁、加大盒子或修改后续曲线都无法消除该冲突。本轮保持原血管、端点半径及 3 mm 约束，不擅自接受这个例外。

已输出的部分芯体包含 1 个连通实体，闭合检查为 True，退化面为 0 个；原血管体积减损的数值检验结果为 2.92e-07 mm³，低于 0.001 mm³ 数值容差。盒体已开 3 个孔，期望完整设计应为 4 个。

**设计参数的含义。**

端点方向通过约 5 mm 的中心线估计，按 0.25 mm 重采样后作直线回归，真实窗口和原采样点数逐口记录。只有方向差不超过 45° 且间距满足要求时才考虑沿切线直连，否则使用 CadQuery 曲线扫掠。弯曲半径至少 6 mm，并以 1.02 倍安全系数筛选；曲率按不超过 0.25 mm 的弧长间隔查询 OCC。血管及端口之间要求 3 mm 表面间距，另留两倍 0.02 mm 网格离散余量。

连接处允许向原封口内侧重叠 0.6 mm。碰撞检查仅豁免同一父分支端点附近 5 mm 的连接区；其他分支及更远父分支仍参与检查。表面归属依据已有 SWC，以 0.2 mm 间隔关联；一个面仅在三个顶点都属于此连接区时才豁免。表面最短距离由 FCL 计算，不以最近网格顶点距离代替。

候选侧壁按距离、方向变化、碰撞和边缘条件评分；权重依次为 0.35、0.35、0.2、0.1。距离除以盒内对角线，角度除以 180°，不满足条件的候选被拒绝。优先入口和出口相对，若曲线总转向超过偏好值 120° 则尝试其他布局；出口最多占 2 个相邻面。所有候选（含被拒绝者）的分数、坐标和原因保存在 CSV。

为避免把目标孔强制放在端点同一高度，程序还尝试沿端点方向顺势转弯。试用的弯曲半径为 [6.2, 8.0, 10.0, 12.0] mm，各用 17 个圆弧导点交给 CadQuery 拟合，转弯后直达盒壁。另一组候选在壁面横向偏移 [0, -10, 10] mm、竖向偏移 [0, 8, -8] mm；两点样条的端切线幅度取端点间距的 [0.7, 1.0, 1.5, 2.0] 倍。每面保留至多 8 个不同目标位置进入组合比较，曲线在穿墙前至少直行 3 mm。

几何离散角度精度为 0.1 rad，来源坐标核对容差为 0.0001 mm，面积不大于 1e-10 mm² 的三角形按退化面处理。仅当合并结果出现数值退化面时，才使用 Manifold 官方简化功能，容差为 1e-05 mm；不修改输入 STL，并重新检查体积保留和连通性。完整数值和图像显示参数保存在 [本次配置](resolved_design_config.yaml)。表内空的距离或弯曲半径表示不适用，不代表零。

VMTK 检查：现有 VMTK 可调用，但已验收 STL 没有开放边界，不能直接作 flow extension；采用 CadQuery 局部延长，避免切开或重建原主体。

## 7. 当前能不能进入下一步

**NEEDS_ADJUSTMENT**

当前仍有约束未满足，需调整本轮结构设计后再审核；失败原因保存在 failure_summary.md。

这不是打印批准。两件套没有增加可拆壁或卡扣，刚性芯能否穿入多个固定侧孔、穿墙余量如何密封，以及实际打印误差，仍需人工确认。

本次 STL/STEP 文件名含 PARTIAL_REVIEW_ONLY，仅供检查已完成部分，不能作为全部端口完成的设计。没有导出完整最终芯体或完整最终盒体，3MF 暂不生成。

审核入口：[装配图](QC/08_assembly_preview.png)；[顶视图](QC/09_top_view.png)；[交互装配](assembly/assembly_preview.vtm)；[端口表](tables/port_layout_summary.csv)；[完整检查](geometry_qc.json)；[源文件哈希](protected_source_hashes.json)。

VTM 在 ParaView 中打开后可分别控制 vascular_core、inlet_ports、outlet_ports、casting_box 四层。3MF 若成功导出则只含两个物理部件。

本轮验证：新增合成与几何测试 **25 项通过**，相关回归测试 **35 项通过**，合计 **60 项通过、0 项失败**。测试过程出现 72 条 VTK/NumPy 弃用提示，不影响几何结果。STEP 已重新导入验证为 1 个有效实体；VTM 四个图层均可读取。完整记录见 [测试结果](validation_results.json)。

保护复核：**1157 个来源文件的 SHA256 全部保持不变**，包括已验收 STL、配套 SWC、原清单、语义精化来源、VascularMD 上游以及 s1-2、s1-3。
