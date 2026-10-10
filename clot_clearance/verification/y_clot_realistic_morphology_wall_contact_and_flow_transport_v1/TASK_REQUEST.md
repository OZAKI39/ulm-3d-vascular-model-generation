# CODEX TASK — Y_CLOT_REALISTIC_MORPHOLOGY_WALL_CONTACT_AND_FLOW_TRANSPORT_V1

你已经知道如何访问和使用服务器。本轮请**直接执行代码开发、有限数量的真实计算、三维可视化和结果审查**，而不是只给实施计划。你有权安装、拉取所需软件包与源代码；优先隔离环境，不覆盖既有安装。

## 0. 当前阶段、目标、为什么做

此前已完成：

- 从平行项目真实文件读取 Y 形分叉血管三维网格、WALL/端口表面与**无血栓状态下**的稳态 CFD；
- 将曲面贴壁血栓置于母管 `x≈−6 mm`、`+z` 内壁，已完成真实壁面贴合、均匀立方 PD 粒子、三层 wall anchors、typed bulk/interface/wall bonds 与 LAMMPS `run 0`；
- 后续 V2 血栓增大到 `a/b/H = 0.68/0.62/0.65 mm`，形成 2035 个 clot nodes 与 1522 个 anchors；
- 同一真实 Y 几何上完成 C0 零载、C1 Gaussian 100 Pa、C2 Gaussian 200 Pa、C3 Opposing 100 Pa 的人工 streaming/shear 碎裂实验；
- 已证明真实断键、附着质量下降、少量小碎屑分离，但没有 ≥200 µm 的游离大块；
- 当前动画三个实质性问题：**血栓形状仍偏光滑/规则，碎屑穿过完整血管壁，碎屑不会随已有背景血流运动；粒子方块的外观不自然。**

本轮任务：

> **在同一已求解 Y 形血管中，构建一个更自然、可复现、贴合原生壁面的 V3 不规则血栓；用完整 Y 形血管壁对 LAMMPS 动力学中的真实可动节点施加非穿透排斥；将短时断裂和较长时间的背景流场碎片输运分开计算；在忠实保留断键拓扑的前提下升级医学风格的三维展示。**

这是一轮 **low-cost qualitative visualization + mechanics integration**，不是高保真流固耦合。最高状态：

`Y_CLOT_REALISTIC_MORPHOLOGY_WALL_CONTACT_AND_FLOW_TRANSPORT_V1 = ACTIVE`

不要因为渲染好看而制造不存在的裂纹、游离大块或支路选择。先修**物理穿壁和碎片输运**，再修**视觉外观**。

## 1. 冻结的项目政策与禁止事项

继续保持：

```
PRIMARY_FRAGMENTATION_SOLVER = LAMMPS_BPM_PERI_LPS
CONTACT_MODEL = MODIFIED_1P0DX / CONDITIONAL_EXPLORATORY_DEFAULT
FRACTURE_SCOPE = EXPLORATORY_ONLY
B1_B2 = TERMINATED_BY_PROJECT_POLICY_CHANGE
PD2_FRACTURE = STOPPED
PHASE_FIELD_FRACTURE = ARCHIVED_OPTIONAL_REFERENCE_ONLY
CLOT_PRESENT_CFD = NOT_SOLVED
BUBBLE_STREAMING = NOT_YET_CONNECTED
```

**严禁**：恢复 B1/B2、重启高阶 PD2、增加新 fracture law 或软化模型、重校准 Liu 材料卡、扫描材料/接触强度、重跑 37 个旧 map 点、进行 GPU/MPI 压力测试、无关网格收敛研究、重新求原 Y CFD、未经许可覆盖平行项目的任何文件。仅运行本轮实现是否真实有效所需的少量计算。

开始先记录当前 branch/HEAD/dirty files 与结果目录。当前可能含大量未提交工作；不得使用 `git reset --hard`、`git clean`、`git stash`、可能覆盖未提交内容的 checkout 或移动目录。新模块与结果追加写入；旧历史保持原字节与 SHA。

## 2. 真实资源：绝不能重画假的 Y 管

优先沿用上一轮已核验的只读资产、副本与配置；缺少时从实际平行项目读取：

```
/workspace/mevo_acoustic_microbubbles/configs/stage4/bifurcation.yaml
/workspace/mevo_acoustic_microbubbles/results/stage4/simvascular/mesh/20261002T141441884393Z/SV_MESH/mesh-complete.mesh.vtu
/workspace/mevo_acoustic_microbubbles/results/stage4/simvascular/steady/20261002T141810826990Z/run/8-procs/result_121.vtu
```

网格同目录 `mesh-surfaces/` 包含原生 `WALL`、入口和两个出口。服务器平行项目不能访问时，优先使用上两轮已下载/校验的本地只读镜像，不得凭几何尺寸制造替身。

源数据身份：理想化、对称 Y 型、母管 R=1.3 mm，子管 R=1.0 mm，三端口；母管沿 x 轴，第一血栓位于 x≈−6 mm、+z 内壁。稳态无血栓 CFD 为水的不可压缩层流，入口平均 5 mm/s，原入口约 1.58124 mL/min，速度最大约 10.074 mm/s，**以原生实际文件读取值为准：26,712 节点、146,230 四面体。** 背景流文件是时间步 121、物理计算时间 1.21 s。严格区分这个算例与 BraVa 真实四端口模型。

保存 source register、路径、时间戳、hash、坐标单位及 VTK 数组名。可以复用历史核验，无需再做一次全文件大规模 hash campaign；只检查本轮实际加载的资源标识。

**不得把原背景 CFD 当成含血栓 CFD；它不能产生当前血栓面的真实剪切牵引。** 本轮对血栓碎裂仍用人工 Gaussian shear，对脱落碎片仅用原背景速度做单向输运演示。

## 3. V3 血栓：比 V2 更复杂，但力学形状仍要可解析

不要只把 V2 的随机噪声放大；也不要改变 Y 血管、血栓轴向位置和壁面身份。

从现有曲壁附着生成器演进出可冻结、可配置的：

`CURVED_MURAL_MULTILOBE_CLOT_V3`

几何目标：

- 保持同一原生 WALL 上的连续附着足迹，主要方向沿 +x 主流；
- 在大致 `a≈0.68 mm, b≈0.62 mm, H≈0.65 mm` 尺度上构造可辨认的 **2–3 个宽尺度非对称隆起**、较厚上游肩部、渐薄下游尾部、非镜像侧缘；
- 具有自然不规则边缘，但无尖刺、薄丝、悬空岛、非物理洞或自交；
- **真正进入力学几何的最小特色尺度**不得小于几倍 `dx=50 µm`，重要隆起/凹陷优先达到 `~0.2–0.3 mm` 量级；`δ=3.01dx` 邻域连续性要保留；
- 不规则性用确定性多尺度基函数、平滑椭圆扰动与固定 seed 定义，所有控制点、幅值、尺度在正式动力学前冻结；
- 表面微观纹理/颜色/roughness 可以单独用于渲染，但不冒充真实微米级力学结构；
- 为隔离形态因素，V3 连续实体体积尽量在 V2 `0.285677 mm³` 的 ±10% 以内；如自然形态做不到，写明原因并报告真实体积差异，不以谎报等体积通过；
- 严格保持新 clot volume 全部位于原生流体管腔内、与壁面无显著间隙；固定锚点位于真实 WALL 外侧；
- 重新用全局固定 50 µm cubic lattice 生成全部 material points / three-layer anchors / typed bonds，不弯曲旧点云并继承旧 bonds；
- 删除跨真实自由面/错误穿腔的非法键，报告新积分体积与连续参考体积的偏差，并明确这种差异仍受 dx 影响。

**只生成一个冻结 V3，不做形态大参数 sweep**。需要图示 V2/V3 同相机、同实际比例比较；若新增的 lobes 在正常尺寸下肉眼根本看不出，正式求解前允许一次确定性形态修正，之后严禁看碎裂结果再调形。

## 4. LAMMPS 材料、边界与载荷保持原身份

继续复用已有可工作的代码和同一 LAMMPS `1.0dx` contact variant；不创建第二套 solver：

```
K = 33333.333333333336 Pa
G = 940 Pa
rho = 1000 kg/m^3       # engineering assumption
s00_bulk = 0.15
s00_interface = 0.140624859375
alpha = 0
dx = 50 um
delta = 3.01 dx
contact_distance_factor = 1.0
numerical_drag = 500 / s   # not real fluid drag
```

Liu 2024 whole-blood LPS initial modulus mapping 与 porcine-skin peeling interface proxy 身份必须保留，不得改称 patient-specific 血管壁粘附。血栓节点均可动，wall anchors 固定，interface 可断。

本轮动力学只需要：

- `C0`: 新 V3 的 zero-load 短时稳定性检查（要实际积分，以排除壁面碰撞/几何初始预载导致自裂）；
- `C1`: 新 V3 **Gaussian shear 200 Pa**，沿当前空间位置计算载荷，使用上一轮人工脉冲定义及约 `1.1605 ms` 的绝对脉宽、`dt=0.5 µs`、总观察至约 `2.33 ms`。载荷大小是 engineering demonstration，非当前 CFD 推出的 traction。

若 `C1` 产生真实游离片，立即进入输运。若完全未产生游离片，不准增加载荷来强制制造；保留 `NO_DETACHED_FRAGMENTS`，可以利用上一轮 V2 的**已保存真实脱落状态**独立验证后续输运并单独标为 `V2_REUSED_TRANSPORT_DEMO`，不能与 V3 的碎裂视觉拼成同一个场景，也不得冒充 V3 结果。

本轮最多 **2 次正式 LAMMPS 动力学运行**（C0/C1）；基础装配/run0 和必需的最小软件检查不计科学工况。不跑 100/200 Pa 多模式矩阵，不做 37 工况地图。

## 5. 核心物理任务 A：让完整 Y 形血管壁真正阻挡碎屑

本轮最重要的修复是 **native Y WALL triangle mesh -> time-stepping contact force**，而不是画面裁剪。

必须使用真实 `WALL` 三角面与原生 tetrahedral lumen membership；入口与两个出口是**开放端口**，不是封闭碰撞帽。原生 WALL 表面本身有端口开口，不能把未封闭侧壁直接假设为闭合 signed-volume。按真实 volume mesh 判断点是否在管腔内，并按最近原生 WALL 三角形获取距离及朝**管腔内部**的法向；对近端口的粒子用单独的 crossing/exit 逻辑，不施加假端口墙力。

接触方案优先级：

1. 使用 LAMMPS 已支持的 `fix external` 外部逐节点力回调，在每个必要时间步从当前 clot/mobile node 位置计算真实壁面最近点、距离与一致的内向法向，施加**仅向内推的短程排斥**；wall anchors 排除。可使用 Python LAMMPS library + `mpi4py` + VTK cell locator / 预计算窄带场。确认其在实际 MPI 实例中有效，不要以离线 Python 后处理冒充实时力。若 Python callback 明显拖慢，用已有本地 C++扩展/预计算插值加速，但不重写 LAMMPS 的整个 pair/bond 模型。
2. `fix wall/region` 只在确实能严格表达同一个原生 Y WALL 时使用；**不能**为了省事以几根理想圆柱取代真实曲面后声称原生壁碰撞已解决。

接触参数（有效材料点半径、近壁阈值、罚刚度、可选仅接近时阻尼）均标注为工程数值选择，必须与现有显式 dt 稳定性匹配；首次接入不得让已贴壁的完整 clot 在零载下被人为弹起/断裂。不能改变 `bpm/peri` 的既有 `1.0dx` **粒子间**接触因子来混淆两种接触。新壁面接触应独立命名，例如 `NATIVE_Y_WALL_REPULSION_V1`。

必要的实际验证：

- C0 中初始、全过程无异常 wall-force-induced 断键、离壁弹跳或显著穿壁；
- C1 采用实时壁力，保存每个时刻的 native WALL clearance、penetration depth、接触节点数、法向总反力、求解中的出界节点 ID、持续时间以及是否合法通过端口；
- 基于**真实动力学坐标**进行壁面近穿越/扫掠检查，不能用渲染时隐藏、删除或夹回粒子掩盖穿模；
- 不能把位于 `INLET/OUTLET_A/OUTLET_B` 之外的合法离域解释为穿壁；
- 出现穿壁要报告真实最大值、节点质量与时刻并尝试一次有明确数值理由的接触参数修正。仍未达工程几何容差则 `WALL_CONTACT=NOT_ESTABLISHED`，不得宣称成功，也不能在动画中暗中裁掉穿壁粒子。

对单个材料中心，可采用明确的 `r_effective` 数值容差；其表示 PD 体素的碰撞代理，并非测得的血栓粒子半径。报告真实 `distance_to_wall / dx` 分布、允许的近壁重叠阈值和事件证据。不要求实现昂贵的全固体表面接触 FSI。

## 6. 核心物理任务 B：独立的长时间 ONE-WAY 碎片输运

**不要**用 LAMMPS 的 `0.5 µs` 显式时间步硬算 `0.2–5 s` 的血流输运。分成两个物理角色：

```
LAMMPS: bond break / interface failure / detached components / COM and velocities
          t ~ 0–2.33 ms
                ↓ documented single handoff checkpoint
LIGHTWEIGHT TRANSPORT: background CFD velocity + simplified drag + native Y wall confinement
          t_handoff → up to 5 s or port crossing
```

在固定冻结的 handoff 时刻，必须从 **真实 surviving-bond connected components** 导出：

`component_id, member_atom_ids, detached_from_wall, V, mass, COM_xyz_m, COM_velocity_mps, equivalent_diameter, actual_handoff_time, topology_parent_id, geometry_id, contact_model_id`。

- 只输运真正与所有 wall anchors 断开连接的 component。仍附着的主体不得被移动到流线上。
- 对 `≥3` 节点的 resolved detached components，用质心/等效水动力直径做简化**刚体平移**与拖曳；速度从 LAMMPS 实际 COM 初速度继承，不得瞬间强行替换成流速。
- 对 `1–2` 节点碎屑，保留总质量和来源并单独表示为 `NUMERICALLY_UNRESOLVED_DEBRIS` 的流动示踪/质量包。它们不是可信的 50 µm 或 60 µm 生物微粒，不得直接拿其网格尺寸推导真实粒径风险。
- 可使用简单 Stokes relaxation 模型 `m dv/dt = 3π μ_f d_h (u(x)-v)` 并明确 Re 适用范围；如果 Re 超出线性阻力适用范围，可采用有来源的常用有限-Re 阻力修正并标记，**不要**构建完整两相 CFD。水：`rho_f = 1000 kg/m³, mu_f=0.001 Pa s`。允许对各碎块不同 size 使用相同已声明 law。
- 每个运动步查询原生稳态 CFD 的 P1 tetra velocity interpolation。来源必须是实际 `result_121.vtu`；不得造解析 Poiseuille 流线取代已算数据。
- 使用原生 Y WALL 的几何约束处理粒子/组件 near-wall contact、停留、回弹或有限滑移（工程约束，非测得黏附），禁止穿壁；允许经入口/出口开放平面出域，停止跟踪并报告事件；将到达距离子管分叉 10 mm 的既有分类截面和 12 mm 物理出口分别报告，不混为一谈。
- 初始 fragment 可能贴近 wall、旧 CFD 的近壁梯度未被独立认证；允许 `NEAR_WALL_UNCERTAIN`、`STALLED`、`NOT_REACHED_BRANCH`。**不能为了让动画看见分支选择就提高背景血流速度或把碎块投影到中心流线。**
- 优先对阻力速度松弛采用稳定半解析/隐式更新，空间积分使用合理的自适应步长，避免因小碎块惯性时间很短而采用微秒级大循环；transport 总时长可到 5 s，实际到出口、停止或越域即结束。
- 不引入随机“血流冲击”或自创涡旋来获得更好看的轨迹。

所有输出明确：

`TRANSPORT_MODEL=ONE_WAY_UNOBSTRUCTED_BACKGROUND_CFD_DEMONSTRATION`

**不是**含血栓 CFD、微泡诱导 streaming、流固双向反馈、真实粒径预测或临床捕获率。

## 7. 真实壁面约束与输运的接口要统一

建议建立复用的 `native_y_wall_query`，供短时 LAMMPS wall contact 和长时 fragment transport 共同查询：

- input: `xyz_m`, optional previous `xyz_m`, radius/clearance proxy;
- output: nearest native WALL point, inward unit normal, signed clearance w.r.t. lumen, exit-port status, triangle/asset ID, validity;
- ports: independently classify `INLET`, `OUTLET_A`, `OUTLET_B`, plus user-defined bubble classification planes;
- use tetra-membership or robust geometry query to decide lumen status; close-wall surface is not itself watertight.

同一位置、同一壁面、同一 SI 单位约定。重用已有 Y source/asset register，不在 renderer 里重写一套不同的几何判据。

## 8. 医学风格的高质量可视化：必须基于真实组件

现在每个 material point 一个红色 cube 看起来非常假。新增两个可以切换/分别导出的展示层：

### A. Scientific solver audit view（必保留）

- 全部真实 cubic material points、surviving bonds、断键事件、wall anchors、当前载荷箭头；
- component / detached / attached 身份；
- 透明原生 Y WALL；
- 真坐标、真时刻、不缩放位移；
- 用于审查任何新增美化是否改写了拓扑。

### B. Medical presentation view（重点升级）

- 使用原生 Y 血管 mesh，浅色半透明血管、局部 cutaway/clip 提升可见性；
- V3 血栓为暗红/深红、有温和湿润高光、非均匀粗糙度/色阶与 2–3 个明确隆起的**连续实体**，而非红色方块；
- 对动态碎裂，**按每个真实 surviving-bond component 单独**从其粒子体素提取/重建封闭外壳，可应用有限的 per-component 表面平滑、抗锯齿、法向修正；任何两个互不连接的 component 不得跨组件合并，也不能填补真实裂缝；
- 小于模型可解析尺度的 1–2 node debris 使用细点或小云团，可调透明度，但不得隐藏其统计质量或假装其形态有物理分辨率；
- 静态/动态展示材质微纹理是视觉层，不允许增加 solver 中不存在的 3D 几何大块；
- 输运阶段的 resolved fragments 可把同一时刻的组件几何作为**刚体平移**展示，不得无依据地添加变形、旋转或另造碎块。未解析小碎屑作为明确标记的 tracer/cloud；
- 视角至少有：全 Y 透明外景、血栓附着局部近景、碎屑在分叉附近的输运视角；统一真实尺度，附上流向、mm 标尺与 `BACKGROUND CFD / ONE-WAY DEMO` 标签；
- 若为了医学美感使用 Blender 或更好的离线渲染器，可以单独隔离安装，但原始 PyVista/ParaView 求解审查视图必须保留，且最终视频必须能追溯到原生时间和 component IDs。

绝不能用 renderer 遮蔽真实穿壁材料点来让其“看似留在血管内”。Presentation 可对**血管壁做截面裁切**，但所有血栓节点/碎片统计必须完整；需要同时提供“完整未裁切物理审查视图”。

## 9. 可视化交付必须形成连贯的故事

至少生成：

1. `V2_vs_V3_morphology_and_wall_fit.png/pdf`：相机一致、实际尺寸、形态差异、原生壁面 fit、连续体积/节点体积；
2. `native_y_wall_contact_and_penetration_audit.png/pdf`：新血栓、真实 Y WALL、合法端口事件、最小 clearance 与穿壁计数随时间；
3. `fracture_and_detachment_in_y.png/pdf`：一个新 V3 实际 case 的初始、中期、末态；保留/游离质量、断键和组件身份；
4. `one_way_fragment_transport_routes.png/pdf`：原 CFD 速度、真实 handoff 分量、轨迹、支路分类与停滞/出界状态，不假装 CFD 已按血栓更新；
5. `solver_vs_medical_view_same_frames.png/pdf`：相同组件、相同时刻的 raw solver 对照和 medical presentation；
6. `y_clot_realistic_fracture_and_transport.mp4` + GIF preview：视频先演示 LAMMPS 碎裂，再用明显章节标识切换到 one-way seconds-scale flow transport；不得把两种时间步混为一个“连续高保真 FSI”。

另输出可交互 HTML 或 ParaView VTM（任选已验证可靠的一种），附数据。图中文本使用英文科研标签，REVIEW 用中文。不得伪造没有发生的 ≥200/300 µm 大碎块或“碎块顺利进入 A/B”的故事。

## 10. 必须保留的尺寸/风险代理

继续保留项目定义：

- Detached ≥200 µm（按 reference volume equivalent diameter，至少 34 个 50 µm 节点）：文献碎屑统计的**参考尺度**，不是普适危险阈值；
- Detached ≥300 µm（至少 114 节点）：当前 dx 下的**工程可辨大块尺度**，不是患者危险阈值；
- 1–2 nodes = `NUMERICALLY_UNRESOLVED_DEBRIS`；
- 附着主体永远不计入 free large fragment；
- 保存首次 detached 时间、exact-membership lifespan、大小、质量、COM 初速度、后续 transport 身份；
- 断键与碎屑增多不能独立解读为治疗更好，特别是大游离块可能增加远端风险。

若 V3 仍只有小碎屑，请坦率给出；视频只展示真实小碎屑/云团，不用人为拼成大块。

## 11. 工程结构与最小必要检查

复用现有 `pd_clot/lammps_fragmentation/`，建议新增模块：

```
y_clot_morphology_v3.py
native_y_wall_query.py
native_y_wall_contact.py
fragment_handoff.py
background_flow_sampler.py
one_way_fragment_transport.py
medical_surface_reconstruction.py
medical_renderer.py
```

如果需隔离 LAMMPS library build 或 `fix external` adapter，新增在原隔离 `tools/`，记录固定版本和 SHA；不得覆盖现有接触变体二进制/旧原生安装。避免把所有 wall query 写成每一时间步串行遍历整张 Y WALL 的 O(Nwall×Nnode) Python 循环，优先 locator/BVH/窄带预处理/批量查询。

只做与本轮成功直接相关的 checks：

- V3 几何在真实 Y volume 内、表面不自交、一个 clot bulk 主体、合法 wall/interface bonds、初始 wall contact 无非物理预载；
- C0 实际 zero-load stability；
- C1 真实有 wall contact 参与时间积分，壁面/端口分类可靠，穿壁不靠可视化裁剪解决；
- 手交 detached component 数据完全来自原生 surviving-bond topology，mass/COM/velocity/ID 有回溯；
- 5 s 内（或较早结束）的 one-way solver 对有效域、壁面和开放端口行为正确；
- 关键渲染帧没有跨 component 平滑粘接或虚构碎块，原视图与医学视图可用相同时间索引核对；
- 记录 wall/transport/visualization 实际时间和内存即可，不做独立 stress/benchmark campaign。

不要额外运行：新 fracture calibration、pressure-only、full CFD、microbubble acoustic field、mesh/horizon convergence、GPU benchmark、多形态/多载荷 sweep、新 PD2 或 phase-field。

## 12. Stage gates：不得把未完成模型包装成完整演示

按以下顺序推进：

**Gate A — V3 Native geometry**

- 在真实 Y 上贴合并冻结一个有明显自然不规则度的 V3；初始颗粒与接触 topology 合法；保存原生点/三角面数据。失败则不跑新动力学。

**Gate B — Real native Y wall contact**

- 用 C0/C1 证明壁面排斥在实际 LAMMPS 时间积分中工作，独立核查越壁及开放端口事件。
- 若越壁仍严重：保留失败轨迹和诊断，允许一次有依据的数值接触修正；不能通过裁剪“通过”。无法完成则 `WALL_CONTACT=FAIL_OR_NOT_ESTABLISHED`，视频只作有限技术演示且不得宣称正确管壁运动。

**Gate C — Fragment handoff + transport**

- 真实游离碎屑则从 V3 C1 输出 handoff；如果零游离则按第 4 节保留状态，并可独立使用旧 V2 真实手交轨迹测试输运，但不能拼接。
- 必须明确哪些碎块流向 A/B、哪些停滞/仍在入口区域、哪些越出有效 CFD 域或未达到分类截面。无人到分叉不算失败，不得为观感做人工运输。

**Gate D — Medical presentation without false geometry**

- 输出 raw / medical 同帧对照、全 Y 和局部展示、断裂/输运分章节 MP4。质量必须比上一轮的红色 cubic glyph 更自然；如果可视化无法不跨组件修复，则保留 raw evidence 并明确限制。

每个 gate 单独 PASS/PARTIAL/FAIL，失败不能静默切换到美化-only 路线。

## 13. 官方技术入口和允许安装的软件

务必优先核查当前固定 LAMMPS 版本对以下接口的支持，勿假设其他版本的 API 自动一致：

- [LAMMPS fix external callback](https://docs.lammps.org/fix_external.html)
- [LAMMPS Python module / set_fix_external_callback](https://docs.lammps.org/Python_module.html)
- [LAMMPS Python installation and MPI](https://docs.lammps.org/Python_install.html)
- [LAMMPS wall/region（仅分析适用性，不能假冒真实三角面）](https://docs.lammps.org/fix_wall_region.html)
- [LAMMPS bpm/peri bonds](https://docs.lammps.org/bond_bpm_peri.html)
- [LAMMPS bpm/peri contact](https://docs.lammps.org/pair_bpm_peri.html)
- [PyVista](https://docs.pyvista.org/)
- [PyVista select interior / closed-surface caveat](https://docs.pyvista.org/api/core/_autosummary/pyvista.datasetfilters.select_interior_points)
- [VTK](https://vtk.org/doc/nightly/html/classvtkStaticCellLocator.html)
- [SciPy integrate](https://docs.scipy.org/doc/scipy/reference/integrate.html)
- [Trimesh](https://trimesh.org/)
- [OVITO Python](https://www.ovito.org/manual/python/)
- [FFmpeg](https://ffmpeg.org/)
- [Blender（如确实必要）](https://www.blender.org/)

允许安装/拉取任意必要包，记录版本，隔离安装，不改变平行项目环境。

## 14. 交付文件和目录

建议：

```
verification/y_clot_realistic_morphology_wall_contact_and_flow_transport_v1/
results/y_clot_realistic_morphology_wall_contact_and_flow_transport_v1/
```

至少输出：

```
TASK_REQUEST.md
STARTING_GIT_STATE.md
SOURCE_ASSET_REGISTER.json
POLICY_SNAPSHOT.md
V3_MORPHOLOGY_DEFINITION.md
V3_GEOMETRY_CONFIG.json
V3_WALL_FIT_AND_LUMEN_CHECK.json
V3_MATERIAL_POINTS_AND_BONDS.npz
V3_LAMMPS.data
WALL_CONTACT_MODEL.md
WALL_CONTACT_DIAGNOSTICS.csv/json
OPEN_PORT_EVENTS.csv
BACKGROUND_CFD_LIMITATION.md
FLOW_SAMPLER_AND_DOMAIN_VALIDITY.json
FRAGMENT_HANDOFF.csv/json
FRAGMENT_TRANSPORT_TRAJECTORIES.csv/npz
FRAGMENT_TRANSPORT_EVENTS.csv/json
FRAGMENT_MASS_AND_RISK_PROXIES.json
PERFORMANCE_REAL_CASES.json
VISUALIZATION_PROVENANCE_AND_QA.json
FINAL_STATUS.json
REVIEW_Y_CLOT_REALISTIC_MORPHOLOGY_WALL_CONTACT_AND_FLOW_TRANSPORT_V1.md
figures/ + figure_data/
animations/ + movie_frame_manifest/
evidence.zip (lightweight; big native arrays preserved in-place and hashed)
```

每个科学求解必须有 frozen request/config、原生 LAMMPS input/log/dump/events、exit code、runtime、source geometry/flow/contact identities。不可覆盖前轮资产。

## 15. FINAL_STATUS 枚举

至少：

```
V3_MORPHOLOGY = PASS / PARTIAL / FAIL
NATIVE_Y_WALL_FIT = PASS / FAIL
LAMMPS_V3_RUN0 = PASS / FAIL
ZERO_LOAD_DYNAMICS = PASS / FAIL
NATIVE_Y_WALL_CONTACT = PASS / PARTIAL / FAIL
NATIVE_Y_WALL_PENETRATION = WITHIN_DECLARED_TOLERANCE / EXCEEDS_TOLERANCE / NOT_MEASURED
OPEN_PORT_CLASSIFICATION = PASS / PARTIAL / FAIL
FRACTURE_C1 = EXECUTED / NOT_EXECUTED
DETACHED_COMPONENTS = OBSERVED / NOT_OBSERVED
FRAGMENT_HANDOFF = PASS / NOT_APPLICABLE / FAIL
ONE_WAY_BACKGROUND_FLOW_TRANSPORT = PASS / PARTIAL / NOT_RUN
BRANCH_CLASSIFICATION = OBSERVED_A / OBSERVED_B / NOT_REACHED / UNKNOWN
REALISTIC_MEDICAL_RENDER = PASS / PARTIAL / FAIL
SOLVER_VIEW_AND_MOVIE_TRACEABILITY = PASS / FAIL
CLOT_PRESENT_CFD_SOLVED = NOT_RUN
BUBBLE_STREAMING_CONNECTED = NOT_RUN
BIOLOGICAL_QUANTITATIVE_VALIDATION = NOT_ESTABLISHED
PATIENT_SPECIFIC = NOT_ESTABLISHED
FRAGMENT_SIZE_CLINICAL_PREDICTION = NOT_ESTABLISHED
```

## 16. 最终中文 REVIEW 必答

1. V3 是如何在实际 Y WALL 上得到比 V2 更自然的多隆起/非对称形态的？连续与离散体积各多少，能否肉眼区别？
2. 完整原生 Y 壁碰撞如何在**求解时**实际施加，是否覆盖远离附着面的下游管壁？是否仍有穿壁？列最大深度、发生位置、总质量和开放端口分类。
3. 真实 LAMMPS 断裂结果如何？附着质量、游离质量、断键率、碎屑类别、是否 ≥200/300 µm？
4. 选哪个真实 time/state 做 handoff？保留了多少质量、多少 resolved fragments 和多少 unresolved debris？
5. 旧空血管 CFD 在何时、怎样被 P1 查询用于**单向**输运？哪些位置超出可信有效域？
6. 有多少碎片/质量到达 A/B 分类截面或真实出口？若没有，是因速度、近壁滞留、时间不足还是未知？不得假造分流。
7. 碎片是否会再次穿壁？壁面接触/停滞是否在 transport ODE 内实现，而非最后帧修图？
8. 医学展示是否仍完全忠实于求解节点、每个 connected component 与原始大小？为什么看起来比旧方块真实？
9. LAMMPS 两例、wall-contact 查询、transport、后处理和视频分别花多少时间？没有做哪些无关验证？
10. 以后要获得真正血栓表面 traction、真实微泡 streaming、动态 wet surfaces、定量碎片输运各还缺什么物理模块？

## 17. 停止条件

本轮应在以下任务范围完成后停止：**一个冻结 V3 + 一个零载动态 + 一个 Gaussian 200 Pa fracture 动态 + 一个忠实的 native-Y-wall contact + 一条真实 fragment handoff + 一个 one-way background-CFD transport demo + 对应医学风格视频/审查图。**

禁止擅自继续扩展多模式、提高载荷、补做 CFD、加入新流体/疲劳/侵蚀/软化材料、优化 GPU、重开 PD2/B1/B2。完整交付后等待人工审查。