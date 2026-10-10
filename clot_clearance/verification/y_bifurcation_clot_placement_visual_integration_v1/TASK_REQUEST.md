# Codex 开发任务：Y_BIFURCATION_CLOT_PLACEMENT_AND_VISUAL_INTEGRATION_V1

你已熟悉服务器与代码仓库。本轮必须**真正读取文件、实施软件开发、构建血栓几何、生成 LAMMPS 数据与完整三维可视化**，不是仅给研究计划。本轮结果必须由服务器的**实际 Y 形血管网格和已保存的稳态 CFD 数据**产生，绝不能用看起来相似的合成 Y 管替代实际资产。

## 0. 开始前理解现在的项目阶段

项目已从高保真 B1/B2 phase-field 和完整 PD2 fracture 开发转向**便宜、可配置、可真实断键的 LAMMPS 血栓碎裂筛选器**。已完成：

- 使用 LAMMPS `bpm/peri` LPS、三层可断裂壁面锚点、当前位置 streaming provider；
- Liu et al. 2024 同源血栓材料的未损伤小应变映射：`K=33333.333333... Pa`、`G=940 Pa`；
- 工程密度 `rho=1000 kg/m³`，`dx=50 µm`，`delta=3.01*dx`；bulk 临界伸长 `s00=0.15`，interface `s00=0.140624859375`；interface 是 **Liu 同文献猪皮剥离模型的 proxy**，绝非实际血管内壁标定；
- 官方 LAMMPS `stable_30Sep2026` 和独立修改版 `contact_distance_factor=1.0dx`。该接触模型已在 Gaussian/Opposing 少量真实工况中得到有限范围的探索支持，但并非生物接触距离；
- 旧的力–时间地图和全部历史 evidence 只读保存。老 `1.35dx` 与新 `1.0dx` 不能混成同一模型的数据；
- 大块碎片统计已经区分真正游离块、附着大主体、200 µm 文献关注尺度和 300 µm 的工程展示门槛，但不具有临床精度。

**本轮是“真实血管场景的几何、输入资产、可视化接入”轮次，不是新的 fracture 数值轮次。**

### 本轮唯一主要成果

> 把一个贴合真实母管圆柱内壁的三维附壁血栓放入**已有实际 CFD 求解所用的 Y 形分叉几何**，以相同坐标/单位生成合法的 LAMMPS 血栓 bulk / wall / interface 几何、数据文件、载荷/坐标契约；将真实未放血栓的**冻结背景 CFD**叠加在同一幅三维场景里，并提供可审查的实体渲染与原始粒子/连接渲染。

直接回答：“此血栓放置是否几何可行、哪里合适、什么形状比较像附壁血栓、怎样接流场及以后接微泡 streaming？”

不得将“显示在同一个画面中”冒充“已完成含血栓 CFD / FSI / 壁面剪应力计算”。

## 1. 严格限制与历史安全

禁止：

- 运行或恢复 B1/B2、phase-field fracture、PD2 fracture、PD2/FEM 比较；
- 新做任何血栓断裂/fragmentation 动力学计算、force-duration map、Gaussian/Opposing/Vortex 参数扫、200/300 µm 统计扩展；
- 重跑 Y 形 CFD、重算微泡、修改 PZT / acoustics、加真实流场双向反馈；
- 新建高保真材料、本构软化、VES、改变 contact 修正、重新标定 Liu 的 `s00`；
- 合成 4k/10k/20k 性能压力测试、GPU/MPI scaling、网格无关性 campaign；
- 以纯解析 Y 形几何取代已有 SimVascular/TetGen 的实际三维几何；
- 误将 BraVa BG001 RMCA 四端口案例的 **18 mL/min** 或三出口设置用于本 Y 形案例。

保留所有未提交工作。工作开始先保存 Git 分支/HEAD/status、相关已有目录及数据哈希。严禁 `reset --hard`、`clean`、`stash`、危险 checkout/覆盖、删除既有结果。新增代码优先写入血栓项目；平行微泡项目必须**只读访问**，不得修改其 config、mesh、steady CFD 或 import verification。无需 commit/push/PR。

## 2. 两个项目位置与真实资源定位

### 平行项目：Y 形血管 + 真实冻结背景 CFD

根目录：
`/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles/`

以下是用户给出的实际入口，必须在服务器中读取并逐级解析它们引用的真实文件：

- `/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles/configs/stage4/bifurcation.yaml`
- `/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles/results/stage4/simvascular/mesh/STATE.json`
- `/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles/results/stage4/simvascular/steady/STATE.json`
- `/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles/results/stage4/simvascular/import/STATE.json`

### 当前血栓项目

先在服务器定位现有 `clot_clearance/`、`pd_clot/lammps_fragmentation/`、verification/results 和独立 LAMMPS `contact_distance_factor=1.0dx` 变体。历史工作树根目录可能为：
`/home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930`

这个路径只是定位线索，先实际确认，不得以不存在的相对路径假装工作完成。

不得修改平行项目代码/原始数据；新适配器应使用一个独立 source-config JSON/YAML 记录外部只读路径、各数据文件 hash、单位、坐标变换与来源。

### 必须真正核实的身份

Y 形血管理想化但是真实已求解的那套源资产：

- 母管：沿 +x，x ∈ [-14,0] mm，`R_mother=1.3 mm`；
- 子管：长度 12 mm、半径 1.0 mm，偏转 `±30°`，终点约 `(10.392,±6,0) mm`；
- 分叉在原点、母管/分支中心线在 xy 平面、完整三维圆形截面；
- 隐式平滑并集过渡长度参数 `0.20 mm`，端口为平面截面；几何生成过程曾强制表面左右镜像，但四面体网格不必单元级镜像；
- CFD：水 `rho_f=1000 kg/m³`、`mu_f=0.001 Pa·s`、固定无滑移壁面、不可压缩三维 Navier–Stokes、恒定流量入口；
- 入口截面平均速度 `5 mm/s`，实际总流量约 `1.58124 mL/min`；两个出口零参考压力自然牵引，不能人为强制各 50%；
- TetGen/SimVascular：目标边长约 0.25 mm、约 26,712 节点与 146,230 四面体；
- 稳态速度最大约 10.074 mm/s；入口截面平均压力约 0.70494 Pa；实际两支流量约 0.79056/0.79068 mL/min；
- 原求解结果为冻结**无血栓**背景 CFD，steady final step 121，时间约 1.21 s；微泡 P1 四面体插值检查过 1000 个采样点。

以上均为**待从实际 state/mesh/solution 核实的用户交接数据**。逐项将 `USER_HANDOFF / ACTUAL_FILE_VERIFIED / MISSING_OR_DIFFERENT` 标记在 ASSET_AUDIT。不能只抄本任务数字就写“验证通过”。如实际资产与交接描述不符，优先报告证据和路径，而不是为了匹配数字篡改文件。

如果找不到真正的 Y 形网格或 CFD 数据，禁止用解析/合成场假冒，明确标记相应资产缺失；仅做能够由真实现存资产支撑的部分并停止需要该资产的功能。

## 3. 冻结第一版血栓位置（只放一个，不做位置扫）

**目标：分叉前、母管正 z 方向的内壁，附壁非闭塞型血栓。**

中心：

```
clot_axial_center_m = [-0.006, 0.0, 0.0]  # 中心线投影，仅定义母管局部坐标
nominal_wall_footpoint_m = [-0.006, 0.0, +0.0013]
nominal_flow_direction = [+1, 0, 0]
nominal_wall_to_lumen_direction = [0, 0, -1]
```

- 母管已求解的**真实内表面**是最终参考，不可仅凭 `z=1.3 mm` 判定贴合；从实际网格提取母管内壁、用最近面/局部圆柱坐标定位该壁面脚点；
- 如果真实内壁表面相对理想圆柱的偏差超过 `dx` 或几何定位歧义，报告并修正局部投影规则，不能静默挪血栓位置；
- 选择 `x=-6 mm` 是为了远离 x=-14 mm 入口和 x=0 分叉复杂连接，不要在本轮把血栓再放入两子管或分叉尖端；
- 不要从旧半球“平移+弯曲”旧 bonds；必须从新实体域重新生成规则 cubic-lattice material points 和 bonds。

## 4. 生成贴合圆柱内壁的附壁隆起，而不是悬空平底半球

母管名义半径 `Rv=1.3 mm`。定义：

```
r = sqrt(y*y + z*z)
theta = atan2(y,z)             # +z 壁面对应 theta=0
s2 = ((x-x0)/a)^2 + ((Rv*theta)/b)^2
h(x,theta) = H*(1-s2)^(3/2)  if s2<1, else 0
```

名义几何初值（本轮一次冻结）：

```
x0 = -6.0 mm
a  = 0.5625 mm   # 轴向 footprint 半长度
b  = 0.5625 mm   # 贴壁周向 footprint 半宽
H  = 0.5625 mm   # 朝管腔中心的最高隆起
```

血栓实体：

`Rv - h(x,theta) <= r <= Rv`，同时必须在实际 CFD 母管管腔内部。用于贴合时以**实际母管内壁表面/径向交点**校正墙面位置，保留上述 profile 为局部高度定义，不重新发明整个 Y 形血管。

重要：`a=b=H=0.5625 mm` 仅表示与原来半球类似的**名义尺寸量级**，新形状体积不必等于原半球体积，必须测量并报告新体积；不许说“新旧血栓体积一致”。

要求：

- 墙面接触足迹有实面积且持续连通；中心向管腔内部突出，边缘高度平滑趋零；
- 避免墙面穿透、虚假的切平面悬空、外露面错误穿出血管、血栓接触对侧壁和遮住三端口；
- 全面覆盖曲面附着区域，不生成薄到完全不可能被 `dx` 表达的悬空单粒子“胡须”；允许合理的 voxel stair-step 并如实描述；
- 坚持 3D 体积血栓，不能用 STL 表面的点代替实体材料点；
- 输出血栓体积、bbox、中心、表面面积估计、近壁 footprint、最大伸入管腔高度、截面占据比例（明确其定义），并与原半球体积单独对照；
- 首版只用一个形态 `CURVED_MURAL_DOME_V1`。后续可支持轴向拉长和不对称形态，但本轮不批量新建多个形态。

## 5. 新的 LAMMPS cubic-lattice body + 曲壁锚点

复用已证明工作的 `geometry_adapter`、`LAMMPS writer`、材料/bonds/API，添加 `Y_CYLINDRICAL_MURAL_CLOT` 作为 geometry variant；不是重写 LAMMPS 或 PD。

固定：

```
dx       = 50e-6 m
delta    = 3.01*dx
K        = 33333.333333333336 Pa
G        = 940 Pa
rho_clot = 1000 kg/m^3           # engineering assumption
s00_bulk = 0.15
s00_interface = 0.140624859375
alpha    = 0
contact_model = MODIFIED_1P0DX   # 既有独立 LAMMPS 变体，非官方默认
```

要求：

1. 使用同一个全局等间距 cubic grid，生成曲面附壁血栓的体积节点；每个体积分量真实可数，坐标以 **m** 写入，名义粒子体积 `dx³`、质量 `rho_clot*dx³`；
2. 在**真实管壁之外的固体侧**生成三层固定锚点，尽可能沿曲壁局部外法向分布；锚点位置仍落在与血栓一致的 global cubic lattice（按 mask 选择），不能为贴合曲面随意移动格点，使 LAMMPS 不再满足 cubic spacing；
3. `BULK` bonds 只连接血栓体积点；`INTERFACE` bonds 只跨血栓贴壁近场与壁面锚点；`WALL` bonds 只连接锚点。使用参考点间距/实际接触区域/δ 创建；禁止跨越空腔或血栓薄缘产生非物理穿空连接；
4. 锚点不计入 clot 质量、体积、fragment count；固定 wall layer 不是 deformable vessel wall 或 FSI；
5. 必须区分 `clot`, `wall_anchor`, `bulk_bond`, `interface_bond`, `wall_bond`, `exposed_surface`, `attachment_surface` 等显式 tag；
6. 使用现有 `1.0dx` 接触变体，不改变其源码/参数、不重校准 s00；
7. 曲壁界面刚度/断裂强度和原来的平面三层锚点 geometry 不同，历史 Liu 的宏观 interface calibration 不能自动视为已在此新几何重新验证：结果字段必须设 `CURVED_INTERFACE_STRENGTH_VALIDATION = NOT_ESTABLISHED`。本轮只是**几何装配**；
8. 实施一个 LAMMPS `read_data` + **`run 0` 读取/初始化检查**，验证新几何 LPS/bond/interface/tag/anchor group 被正确加载，且无缺失键、重复坐标、NaN 或致命错误。`run 0` 不是动力学科学案例，不输出“碎裂成功”结论。

若 `run 0` 发现当前 LAMMPS 对象标签或固定壁层不可兼容，只修必要文件配置/adapter；不要因此启动 fracture 算例。

## 6. 实际 Y 型网格与新血栓的空间交集检查

必须读取实际 Y 形管壁三角网格、CFD 体四面体网格（若格式需要，允许使用 `meshio / pyvista / scipy`），并进行少量**真正必要**的几何检查：

- 原始入口/出口数为 3，几何 bbox/体积、分叉方向符合资源身份；
- 血栓体素/实体内部样点位于实际 CFD **原管腔**内部；锚点在管腔外/固体侧，clot-wall 无显著穿透；
- clot 在管腔内占据有限体积，能够构造局部 exposed surface 与壁面 attachment patch；
- clot 坐标与真实 mesh 的 m/mm 单位一致，且足迹实际贴近管壁（记最大/中位距离、相交面积或 gap）；
- 含血栓与不含血栓的 domain 身份明确区分；
- LAMMPS 材料节点/锚点 group、三类 bonds、单节点体积/质量、分组和总连接情况有限且合理；
- 实体血栓不与分叉管腔边界/出口错误相交。

检查仅服务于**装配可信**，不展开 voxel resolution convergence、CFD flow convergence、contact physics validation 或 mesh independence campaign。

## 7. 用现有 CFD 实际数据叠加显示；绝不冒充含血栓流场

正式读取 `steady/STATE.json` 指向的**真实求解速度和压力**，使用原有 P1 四面体插值/导入机制或等价只读解析，不静默改写数据。重要：

- **该 CFD 是在没有血栓的原 Y 型管腔上算出的**，属于 `UNOBSTRUCTED_STEADY_BACKGROUND_FLOW`；
- 现有数据允许用于三维流线/速度模、背景流向、血栓空间定位及坐标系统一；
- 新血栓内部原 CFD 数据虽然存在，但**该区域现在物理上是固体**，渲染时必须遮罩或截断旧流线，不能显示水流穿过实体血栓；
- 即使在渲染中遮掉流线，也不能声称流线已经绕开血栓；标题或图注必须醒目标出 `Original unobstructed background CFD; not recomputed with clot`；
- 血栓表面准确 pressure traction / wall shear stress **尚不可从此不含血栓的旧 CFD 认证**；不在本轮直接计算旧背景流的“真实血栓表面剪切力”并施加到 LAMMPS；
- 之前的 1000 点插值验证属于原平行项目历史证据，本轮不重复 1000 点压力测试；只做一组能识别单位、索引及插值有效性的轻量只读 spot checks；
- 保留 CFD 全域速度模/压力的原始单位和 source range；禁止手工画假的微泡群/流场涡旋/声流来提高视觉效果。

将背景 CFD 的入口、两出口、速度场和新血栓**使用同一 xyz 坐标系**可视化。对流线交叉新血栓的情形可以 mask/clip 为 unavailable，而不是声称真实绕流；记录剪切/压力缺口供以后 `CLOT_PRESENT_CFD` 任务处理。

## 8. 为未来背景流和 bubble streaming 冻结空间加载契约

新增一个与 LAMMPS 接口兼容的、**尚不施加真实载荷**的 geometry-aware surface mapping 数据文件。

对曲面血栓 exposed surface 保存：

- 全局 `position_xyz_m`、外露节点 ID、`normal_xyz`、`A_i_m2`（面积权重及其计算方式）、`wall_distance_m`；
- 每个节点相对于当前真实血管局部坐标的 `axial_tangent_xyz`、`circumferential_tangent_xyz`、`wall_to_lumen_normal_xyz`；
- `local_flow_direction`, `geometry_id`, `CFD mesh/field ID`, source frame, unit convention；
- 分辨初始外露面、贴壁面及断裂后新增表面（后者仍不支持动态更新）。

未来 provider schema：

```
load_source = BACKGROUND_Y_FLOW | ANALYTIC_STREAMING | BUBBLE_STREAMING | CFD_STREAMING
surface_node_id
position_xyz_m
force_xyz_N or traction_xyz_Pa and nodal_area_m2
time_s or time series
source_geometry_id
source_flow_id
coordinate_frame
source_units
validity_domain
```

特别要求：`BACKGROUND_Y_FLOW` 目前状态必须是 `REQUIRES_CLOT_PRESENT_CFD_FOR_TRACTION`，不得把旧无血栓场默认为可用于新血栓表面的准确 traction；
`BUBBLE_STREAMING`/`CFD_STREAMING` 为以后接入的 provider，不在本轮计算。沿真实曲壁切向投影，不能继续直接沿用平底半球原点的径向方向作为 streaming tangent。

面积 `A_i` 不能简单把所有体积粒子都算成外露表面，也不能同时在 attachment 和 exposed 表面重复计数。该映射在本轮只用作**接口与可视化**，不用于真实力学结果。

## 9. 正式可视化：医学观感 + 可审查来源

构建**同一个真实场景**的两种视图和相关数据：

### A. 几何/展示视图

- 真实 Y 管腔表面（来自网格），半透明；三条管道真实比例，母管方向及左右支路清楚；
- 贴在 x≈−6 mm、+z 内壁的光滑隆起型血栓，使用接近 clot 的深红色/暗红色；
- 局部剖切或透明视窗，使血栓不被蓝色血管壁挡住；
- 标签明确 `idealized Y-vessel`、`curved mural clot (synthetic demo morphology)`、长度/坐标/尺寸；
- 渲染表面须由**当前实际 voxel/mesh**得到，允许只在**视觉层**轻度平滑/抗锯齿，不得改变 LAMMPS 节点/bond。不能为了好看发明另一个血栓轮廓；
- 给出一个局部放大视图，直观看到贴壁、边缘和进入管腔的高度。

### B. 真实 CFD 背景叠加视图

- 对真实 steady velocity 作流线/箭头、速度模着色、可选压力切面；
- 显示原入口到两支路的完整流动方向；
- 复用上面完全相同的 Y 几何与血栓位置；
- 显著注明 `UNOBSTRUCTED BACKGROUND CFD — NOT CLOTTED FLOW FIELD`；
- 旧流线不能穿过新血栓仍被宣称绕流，实际渲染时可 clip/mask；
- 不把未知的 near-clot wall shear stress 画成已经模拟完成的物理结果。

### C. 求解审查视图

- 显示 cubic LAMMPS clot material points、三层 curved wall anchors、bulk/interface/wall bonds、exposed/attachment groups；
- 血栓体积、锚点位置、接触边界与实际三角网格对照；
- 图注用真实粒子数/bond 数，不捏造数量。

至少生成以下 4 张图：

1. `y_bifurcation_clot_overview.png` 和 PDF —— 完整 Y 管 + 血栓位置、三个端口与尺寸；
2. `curved_mural_clot_wall_fit.png` 和 PDF —— 真实管壁与新血栓剖面/局部放大，含最大壁面 gap/穿透检查；
3. `original_background_flow_with_clot_overlay.png` 和 PDF —— 真实 CFD 速度流线 + 血栓 overlay，明确“不是含血栓 CFD”；
4. `lammps_curved_interface_audit.png` 和 PDF —— 真实 LAMMPS points/bonds/groups/wall anchors 审查图。

此外：

- 可输出用于 ParaView 的 `vessel_surface.vtp / clot_solid.vtp / clot_points.vtu / anchors.vtu / background_flow.vtu`（根据格式合理选择），加载即可交互旋转、剖切；
- 生成一个短的镜头环绕 MP4/GIF（仅静态几何可视化，**不是断裂动力学**），以方便用户直接看到 3D 相对位置；
- 图中 mm 轴刻度、比例尺、legend、相机角度必须统一并便于阅读；血管与血栓的相对比例真实，不得为了可见性额外放大血栓或位移；
- 所有可视化由实际已保存几何/CFD 数组离线读取，不在绘图函数里触发任何 CFD 或 LAMMPS 动力学求解；
- 至少保存 figure data/源数组、渲染脚本及可复现日志；检查图像非空、标注未被裁切、内容身份与来源一致。

## 10. 最少必要的测试、绝不延伸

本轮必要测试只围绕这些事项：

1. 外部真实网格/CFD 的路径、格式、hash、坐标单位、域范围能被正确读取；
2. 新血栓全在真实原管腔内；贴真实母管曲壁且 anchors 在管腔外；no duplicate nodes、有限量、合理体积/bbox；
3. `BULK/INTERFACE/WALL` 三类 bonds 类型正确，不连穿流体空隙；分组各自非空；
4. LAMMPS 使用现有变体 `run 0` 可加载此新数据文件；**不时间推进**；
5. 用少数实际节点 spot check 背景 P1 速度/压力数组读回与原工程量纲一致，场外不擅自外推；
6. 四张主要图 + ParaView 文件/短环绕片的来源、数据完整性、图形 QA 通过。

禁止任何数值压力测试、PD fracture tests、Y CFD 再求解、虚构流场或重新校准材料。用户明确要求最少有效验证。

## 11. 新增模块、建议目录和可复用性

优先沿现有 `pd_clot/lammps_fragmentation/` 扩展，名称按项目实际 layout 做最小变动：

```
y_bifurcation_sources.py     # 只读读取平行项目 mesh/steady/STATE/config
y_mural_clot_geometry.py    # 曲壁附壁血栓 SDF/体素/墙面投影
y_wall_anchors.py           # 三层曲壁锚点 & INTERFACE bonds
y_surface_coordinates.py    # 外露面与局部切向/法向/面积
y_background_flow_view.py   # 旧原始背景 CFD 只读绘图/遮罩
y_scene_visualization.py    # 原位血管/血栓/求解审查场景
```

只添加最小模块，优先复用现有 `geometry_adapter.py`、`lammps_writer.py`、runner 与可视化工具。不要另外复制出独立的第二套 LAMMPS solver。

请求参数文件应允许将来更改位置、footprint、H 与素材文件；本轮正式案例只冻结并执行上述一个参数组合。

## 12. 资源与性能

此轮主要是几何和三维可视化，无需新 CFD / fracture；先复用当前可用 Python 环境。

允许安装需要的 Python 包及系统工具（例如 `pyvista`, `vtk`, `meshio`, `trimesh`, `scipy`, `imageio`, `ffmpeg`），尽量隔离安装、不破坏平行项目环境或现有 LAMMPS。可利用 CPU 或现有渲染后端，不必做 GPU benchmark。RAM 不足时分块处理实际 tetra/point arrays，别生成无关超大 mesh。

官方参考链接：

- PyVista：https://docs.pyvista.org/
- PyVista UnstructuredGrid：https://docs.pyvista.org/api/core/_autosummary/pyvista.UnstructuredGrid.html
- Trimesh：https://trimesh.org/
- Meshio：https://github.com/nschloe/meshio
- SimVascular：https://simvascular.github.io/
- LAMMPS BPM/PERI：https://docs.lammps.org/bond_bpm_peri.html
- LAMMPS BPM/PERI HowTo：https://docs.lammps.org/Howto_peri.html
- LAMMPS `read_data`：https://docs.lammps.org/read_data.html
- ParaView：https://www.paraview.org/
- FFmpeg：https://ffmpeg.org/

## 13. 输出目录与必交证据

建议：

```
verification/y_bifurcation_clot_placement_visual_integration_v1/
results/y_bifurcation_clot_placement_visual_integration_v1/
```

输出至少：

- `TASK_REQUEST.md` / `POLICY_SNAPSHOT.md` / `STARTING_GIT_STATE.md`；
- `SOURCE_ASSET_REGISTER.json`：两个项目路径、CFD/mesh/config 实际文件、hash/mtime/大小和读取身份；
- `ACTUAL_Y_ASSET_AUDIT.md/json`：真实 mesh、入口出口、steady 字段、单位、速度/压力范围与 handoff 对照；
- `PLACEMENT_CONFIGURATION.json`：x=-6 mm、+z wall footpoint、曲面 footprint/height、坐标法向和几何身份；
- `CURVED_MURAL_CLOT_GEOMETRY.npz`：真实 material points、mask、node volumes、surface/attachment mapping；
- `CURVED_WALL_ANCHORS.npz`、`BOND_GROUPS.npz`；
- `Y_CLOT_LAMMPS.data`、`in.y_clot_geometry_check.lammps`、`LAMMPS_RUN0.log`；
- `GEOMETRY_AND_WALL_FIT_CHECK.json`、`LUMEN_INTERSECTION_AND_CLEARANCE.json`、`FLOW_COORDINATE_CHECK.json`；
- `SURFACE_LOAD_INTERFACE_SCHEMA.json`、`SURFACE_NODES_AND_NORMALS.npz`；
- 4 张主要 PNG/PDF、ParaView 数据资产、短 MP4/GIF、`FIGURE_DATA/`、`VISUALIZATION_QA.json`；
- `PERFORMANCE_REAL_TASK.json`（只记录几何/渲染真正耗时和占用）；
- `FINAL_STATUS.json`、`REVIEW_Y_BIFURCATION_CLOT_PLACEMENT_AND_VISUAL_INTEGRATION_V1.md`、`evidence.zip`（轻量证据包，大 CFD 原文件可只保存 hash/path，不重复打包）。

所有正式图至少一份使用真实 mesh、真实血栓点云和真实 steady CFD 数据；没有的字段要按情况注明缺失，不允许补假值。

## 14. FINAL_STATUS 明确分层

至少给出：

```
REAL_Y_VESSEL_MESH = VERIFIED / MISSING / CONFLICT
UNOBSTRUCTED_STEADY_CFD_IMPORTED = VERIFIED / MISSING / CONFLICT
Y_CFD_SOURCE_UNCHANGED = PASS / FAIL
MOTHER_CURVED_WALL_PLACEMENT = PASS / FAIL
MURAL_CLOT_SHAPE = READY / FAIL
LAMMPS_CUBIC_POINTS = READY / FAIL
CURVED_WALL_INTERFACE_BONDS = READY / FAIL
LAMMPS_DATA_RUN0 = PASS / FAIL
GEOMETRY_WALL_CLEARANCE = PASS / FAIL
BACKGROUND_CFD_OVERLAY = PASS / FAIL
SURFACE_GEOMETRY_LOAD_SCHEMA = READY / FAIL
CLOT_PRESENT_CFD_SOLVED = NOT_RUN
CLOT_SURFACE_TRACTION_VALIDATED = NOT_ESTABLISHED
CURVED_INTERFACE_STRENGTH_VALIDATED = NOT_ESTABLISHED
CLOT_FRAGMENTATION_DYNAMICS = NOT_RUN
BUBBLE_STREAMING_CONNECTED = NOT_YET
PHASE_FIELD_B1_B2 = TERMINATED
PD2_FRACTURE = STOPPED
BIOLOGICAL_QUANTITATIVE_VALIDATION = NOT_ESTABLISHED
PATIENT_SPECIFIC_CLOT = NOT_ESTABLISHED
PRIMARY_USE = EXPLORATORY_GEOMETRY_AND_BACKGROUND_FLOW_INTEGRATION
```

所有源代码、JSON、最终报告必须统一使用 `CLOT_PRESENT_CFD_SOLVED` 这个状态键。

## 15. 成功与停止规则

这一轮成功意味着：

1. 真正读到了另一项目的 Y 形血管和冻结背景 CFD，记录来源；
2. **只建立一个**分叉前母管 +z 内壁血栓；附着表面贴合实际曲壁而不是平底穿壁/悬空；
3. 重建真实 3D 的 cubic material point cloud、三层固定 curved anchors、bulk/interface/wall bonds；
4. LAMMPS `run 0` 成功读取；无新动力学或材料重标定；
5. 三维交互/图片真正显示 Y 管腔、贴壁血栓、原稳态速度/流线；清楚标示“原始无血栓背景流，未重算”；
6. 同时生成给科学审查的原始节点/bonds 视图和可观看的医学实体视图；
7. 已定义未来真实 streaming/CFD 的表面节点法向、切向、面积和坐标输入契约；
8. 提交图、资产、校验日志、review；随即停止。

即使成功也**不得**把它描述为“流体已经绕过血栓”“已经算出了血栓表面真实 shear”“微泡已经在该血栓旁聚集”或“血栓被血流碎裂”。

## 16. 最终中文 REVIEW 必须回答

1. 是否找到了原始 Y 血管/三维体 CFD 的**真实文件**？有哪些路径与 hash？与 handoff 是否吻合？
2. 当前血栓放在真实母管的哪里？贴壁误差是多少？是否侵入壁面/超出 lumen？
3. 与旧平面半球比较，新血栓的体积、footprint、厚度和形状如何变化？
4. 新血栓有多少 LAMMPS nodes、bulk/interface/wall bonds 和三层 wall anchors？`run 0` 是否成功？
5. 从真实 CFD 叠加图中能看出哪些位置/流向信息？哪些不能推断？
6. 是否还保留了原稳态 CFD 的真实进口流量与近各半分流？若没有重算，就只报告原流场，而不是新血栓后的分流。
7. 当前展示图与审查图是否来自同一份真实坐标与 geometry 数据？视觉平滑是否只用于渲染？
8. 未来微泡/CFD provider 需要输出哪些 surface-node 力/牵引数据才能驱动现有 LAMMPS？
9. 现有 `1.0dx` 接触、Liu bulk/猪皮 proxy interface 在新曲壁上的使用范围和未验证事项是什么？
10. 下一阶段若要真正得到血栓在背景血流中的 traction，需要修改哪些 mesh/BC，为什么当前现有 CFD 不够？

本轮完成后等待人工审查，不自动继续 CFD、bubble 或 fracture case。