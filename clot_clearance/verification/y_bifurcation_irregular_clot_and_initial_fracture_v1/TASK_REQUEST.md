# Codex 开发任务：Y_BIFURCATION_IRREGULAR_CLOT_AND_INITIAL_FRACTURE_V1

**执行要求：你已熟悉服务器与现有仓库。本轮要真正修改代码、生成新形态、运行少量真实 LAMMPS 显式破坏计算、交付可视化、结果和审查报告；不要只给方案。** 用户已明确放弃 B1/B2 高保真相场 fracture 和完整 PD2 fracture 开发。本轮继续低成本、可视、可调整材料和载荷的 LAMMPS 断键/碎块路线。

## 0. 当前项目阶段、最终目标和本轮必要性

最终项目仍是：**微泡群集形态/运动 → streaming、shearing 表面载荷 → 血栓附着/破碎/游离碎块 → 低成本筛选、surrogate、后续控制**。当前微泡声流尚未落地，因此本轮使用清楚标记的**人工解析剪切载荷**；不是血流引起的真实 traction，更不是实验校准的治疗效应。

已完成并经过人工审查的阶段：

- LAMMPS `bpm/peri` 的 LPS 破裂、真实删除 bonds、游离连通分量和碎块动力学；
- Liu 2024 同一研究体系的**未损伤小应变弹性映射**和 bulk/interface 失效强度工程映射；
- 三层固定 wall anchors 与 breakable interface；
- 当前位置解析 streaming provider；
- 官方原生断后接触 `1.35dx` 会显著放大粉碎级联，独立接触变体 `contact_distance_factor=1.0dx` 在少量 Gaussian / Opposing 场景中表现较合理；
- 200 µm 文献关注尺度及 300 µm 当前离散网格工程可辨认尺度的游离碎块/寿命统计；
- **真实 Y 管腔**的三维数据读入、母管 `x≈−6 mm`、`+z` 内壁的曲面血栓几何、cubic LAMMPS 点云、曲面三层锚点与 `run 0`；真实无血栓稳态 CFD 的空间叠加可视化。

刚完成的 `REVIEW_Y_BIFURCATION_CLOT_PLACEMENT_AND_VISUAL_INTEGRATION_V1.md` 中原始曲壁 demo 是 `CURVED_MURAL_DOME_V1`，`a=b=H=0.5625 mm`，连续体积 **0.192860 mm³**、LAMMPS 离散体积 **0.170500 mm³**、clot points **1364**、anchors **1125**、bulk bonds **57346**、interface bonds **16186**、wall bonds **37691**；初始 LAMMPS `run 0` 通过，**没有破坏动力学**。这些是已完成基线，必须保留，不要改写历史。

**本轮唯一核心目标**：

1. 将这个母管附壁血栓做得**稍大、更像自然形成的附壁不规则隆起**，仍贴合已求解真实 Y 网格内壁；
2. 用**新的真实 cubic point cloud 和曲面界面**在 Y 场景中运行 **0 Pa 基线 + 三个明确的人工 shear 初步破坏工况**；
3. 输出漂亮但诚实的“完整 Y 管上下文 + 血栓真实断裂运动”视频，以及原始粒子/断键审查视图；
4. 定量区别 **附着主体、暂时游离块、可辨认 ≥200/≥300 µm 游离块、1–2 节点未解析碎屑**。

本轮不是研究哪种碎裂理论最正确，也不是新的 CFD、微泡、临床碎裂预测或性能压测。

## 1. 最高约束：历史、版本、边界与停止规则

开始时先读取最新 `REVIEW_Y_BIFURCATION_CLOT_PLACEMENT_AND_VISUAL_INTEGRATION_V1.md`、上一轮 `FINAL_STATUS.json`、`GEOMETRY_SUMMARY.json`、`SOURCE_ASSET_REGISTER.json`、`SURFACE_NODES_AND_NORMALS.npz`、原 LAMMPS deck、最新版 Liu card、修改接触变体及 200/300 µm 碎块统计模块。**优先复用已保存快照和实际函数，不能从旧图或提示词猜代码接口。**

- 血栓工作树位置线索：`/home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930/clot_clearance/`。实际工作前核对路径。
- 平行项目真实 Y 管资源：服务器 `/workspace/mevo_acoustic_microbubbles/`；已有本地镜像可能在 `/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles/`。**只读**，不得修改其 CFD/mesh/config。
- 原始实际体网格身份：`mesh-complete.mesh.vtu`，**26712 节点 / 146230 tetra**；稳态 `result_121.vtu`、实际 WALL/端口表面与原 `bifurcation.yaml`。先读取上轮的原始路径与 SHA256 登记，避免再次搬运大文件。若源文件缺失，报告真实 blocker，不造合成 Y 管或假流场。
- 不得用 BraVa 四端口旧算例替代当前 Y 三端口模型；**18 mL/min 不是该模型入口流量**。当前无血栓背景 Y 流量约 **1.58124 mL/min**。
- Git 开始记录 branch/HEAD/status/diff；多轮可能仍有未提交改动。禁止 `git reset --hard`、`git clean`、`git stash`、覆盖式 checkout、删除任何旧 verification/results 或覆盖现有工作。新任务输出单独目录。不要自发 commit/push/PR。
- B1/B2、FEM phase-field、完整 PD2 fracture 均保持已停止状态；本轮不调用。
- **不重算含血栓 CFD，不重跑原 Y CFD，不跑声场/微泡/PZT，也不宣称流固双向耦合。**
- 不调 Liu 材料参数、不重做 s00 强度校准、不写新本构、不改 1.0dx 接触 C++、不做 VES/softening、不重开过往 37-case map。
- 不进行 mesh/horizon convergence、GPU benchmark、MPI scaling、4k/10k/20k 压测、几何或 shear 大范围扫参、pressure-only case、无关 unit test campaign。

## 2. 新血栓：只冻结一版稍大且低频不规则的 V2

任务名称：`CURVED_MURAL_IRREGULAR_V2`。

**位置保持**母管分叉前 x=-6.0 mm、y=0 的 +z 内壁，中心线投影 `(-6,0,0) mm`，附着中心在真实 WALL 表面法向拟合得到，不直接拿理想圆柱 `z=1.3 mm` 当最终内壁。母管平均主流 `+x`。新 clot 必须向管腔内部隆起，而不是穿入血管壁。

本轮冻结建议名义尺寸（只是合成应用几何，不是患者解剖测量）：

```yaml
geometry_id: CURVED_MURAL_IRREGULAR_V2
centerline_projection_mm: [-6.0, 0.0, 0.0]
wall_side: POSITIVE_Z_MOTHER_WALL
axial_half_length_a_mm: 0.68
circumferential_half_width_b_mm: 0.62
max_lumen_protrusion_H_mm: 0.65
cubic_spacing_dx_um: 50
horizon_over_dx: 3.01
irregularity_amplitude_fraction_of_local_height: 0.08-0.12  # 上限，不是随机噪声
geometry_seed: 20261008 # 仅在算法需要确定性 seed 时；不得在看碎裂结果后修改
```

几何构造必须在真实 WALL 三角形局部 patch 上进行：沿局部轴向 `x` 与圆周弧长参数化 footprint；旧隆起 `h=H(1-q²)^(3/2)` 可保留为 envelope，但允许低频、光滑的**周向轻微不对称、上游较厚/下游渐细尾部、约 2–3 个缓和的表面凸起**。不规则 footprint 和顶部轮廓都要可见，不只是全局缩放椭圆。所有不规则项应在边界随 envelope 归零，维持连续贴壁；不造尖锐毛刺、空穴、单体素细须或尺度小于约 `2*delta` 的高频表面噪声。

具体解析函数由 Codex 采用可追溯、可复现的**单一参数化**确定，先存 `MORPHOLOGY_DEFINITION.md / CONFIG.json` 并冻结，再建网格和运行力学。不要为了“碎得好看”回调几何参数。不得把当前形态称为真实 CT/MRI 分割或真实患者血栓。

新形态实际体积无需强制等于旧半球，也不要为体积目标而改变物理参数。输出和旧 V1 的 side-by-side：

- `a/b/H` 与实际轴向、周向 footprint 跨度；
- 连续闭合实体体积、LAMMPS 节点积分体积、二者差率；
- `clot points / exposed nodes / anchors / bonds`；
- 真实截面 occupancy（几何截面，不是血流阻力或临床狭窄率）；
- 与实际 WALL 的 attachment gap/穿透；与入口/分叉/对侧壁 clearance。

旧 V1 离散体积曾比连续体积小约 **11.59%**。新几何要如实记录差异，理想控制在约 **15%** 内，但不为迎合数字升级 mesh refinement 研究；若几何质量或相交问题无法排除，应停止动力学而非带着非法模型继续。

## 3. 重新建立 cubic LAMMPS 点云与曲壁三层 anchors

不得简单拉伸或弯曲旧 LAMMPS nodes 后复用旧 bond list。必须从新连续实体生成**全局等距 50 µm cubic lattice**，选择在实际闭合实体内的材料点，`vfrac=dx³`、`m=rho*dx³`；`delta=3.01dx`。复用已存在 `y_mural_clot_geometry.py`、`y_wall_anchors.py`、`geometry_adapter.py`、typed bond writer 等模块。

- 使用实际 Y `WALL` 表面重建三层固定 wall proxy anchors，保留 global equal-spaced grid，wall anchor 位于 lumen 外，不许把粒子拖动到真实壁面破坏 cubic spacing；
- BULK / INTERFACE / WALL 三种 bond types 保持；interface 键跨真实附着 patch，禁止误连到裸露 clot 面、穿过流体空隙或穿出血管壁；wall bonds 只在 anchors 内；
- wall proxies 用于局部固定锚定，不等于真正可变形血管壁，也不代表整个 Y 管道拥有碎块碰撞边界；
- 新 `surface nodes` 的位置、实际法向、顺血流切向、周向切向、面积权重必须全部重算；面积权重应来自同一新暴露连续表面一次性分配，不混入 attachment patch；
- 以当前质量、体积、连接分组为准写 LAMMPS `read_data` deck，并执行一次 `run 0` 检查几何/内力与读取合法性。无须额外 tiny-coupon 或解析 patch campaign。

如果发现材料节点/锚点在真实 Y lumen 意义上越界，或存在严重薄须/孤立点，只允许在**执行动态计算前**基于固定通用几何规则修复、记录前后差异。不要事后为结果修改几何。

## 4. 材料、接触与破坏模型完全复用，禁止重新标定

保留唯一 baseline 材料和接触身份：

```yaml
solver: LAMMPS
bond_model: BPM_PERI_LPS
K_Pa: 33333.333333333336
G_Pa: 940.0
rho_kg_m3: 1000.0                    # ENGINEERING_DEMO_ASSUMPTION
bulk_s00: 0.150000000000
interface_s00: 0.140624859375
alpha_failure: 0.0
contact: MODIFIED_1P0DX            # 工程候选，非真实生物接触距离
numerical_damping_per_second: 500.0
units: SI
```

复用已构建的独立 `1.0dx` 变体和 SI 力密度→节点 N 转换，不修改旧官方 `1.35dx` 可执行文件。开始时检查**修改版 SHA、材料/接触真实运行参数**；若不符则停止并记录，不能默默用原生版本代替。

来源身份：Liu 2024 whole-blood clot **小应变弹性映射** + 同文献 **porcine-skin peeling interface proxy**。后者**不是实测人体血管内壁粘附**；几何从平面变曲面以后，之前 calibrated failure-strength 不等于在新几何中已被重新验证。`EXPLORATORY_ONLY` 是全程硬约束。

## 5. 本轮正式 4 个动力学工况——只做这 4 个

**Case C0 — ZERO_LOAD_DYNAMIC**

- 相同 V2 geometry、相同 LPS/anchors/contact，所有外部流场力为零；
- 运行与主工况相同的观察窗口，判断新曲壁装配有没有自行漂移、异常起始断键、负体积或严重接触问题；
- 这是与新几何直接相关的必要 sanity gate，不是合成压力测试。

C0 若发现明显非零外力下才可能解释的预应力假象、无载大量断键、lost atoms/NaNs，应停止 C1–C3；只排查实际装配/初始化错误，不“调强材料直到不动”。

**Case C1 — LOCAL_GAUSSIAN_SHEAR_A100**

- 参考平均切向 traction budget **100 Pa**；
- 载荷作用于新 clot 初始暴露表面；高斯中心选靠近上游坡面/顶点的、可追溯的单一空间点，避免将高斯中心随 particle ID 飞走；
- 力的方向为实际外露面局部切平面内、以 `+x` 轴向为主的方向，必须使用新 V2 normals/tangent，不可复用旧半球“相对原点径向”的单位向量。

**Case C2 — LOCAL_GAUSSIAN_SHEAR_A200**

- 与 C1 相同几何、pattern、中心、宽度、时间/方向，只把 mean traction budget 改为 **200 Pa**；
- 直接比较是否增加局部断键、附着主体减少、真实 detached mass 和 >=200/300 µm 大块，而不是默认认为更大的力会生成更大的块。

**Case C3 — OPPOSING_SHEAR_A100**

- 与 C1 相同 mean traction budget **100 Pa**、相同几何/时间窗口；
- 在同一新 clot 表面定义两个相反切向 lobes，实际方向投影到各处新外露面切平面；
- 与 C1 比较不同 spatial topology 对裂开/脱落位置、碎块和大块风险的影响。

**共同加载约定**：

- 使用**当前位置空间场** provider，`F_i(t)=A_i*tau(x_i(t),t)`，离场自然衰减；严格记录参照场坐标、峰值、实际载荷预算和 SI 单位；
- C1/C2/C3 统一用历史已用过的**绝对时间脉冲**：约 `1.1605 ms`（25% 平滑上升、50% 平台、25% 平滑下降），`dt=0.5 µs`，观察结束 `2.33 ms`，约 `4660 steps`。这是为了避免同轮混入脉宽变量；旧的 `2T_cross` 只属于旧尺度参考，新几何可以另报告其 Tcross，不要误称新几何恰好为 2Tcross；
- 三种 loaded cases 的名义 budget 定义均为 `Σ_i|F_i(reference)| / Σ_i A_i`。同 budget 不等于同局部峰值，也不意味着相同反力/净力/力矩。输出 actual peak traction 与力积分；
- 两种 100 Pa pattern 在**相同**新 clot material/geometry/area/time/contact 下运行。

人工 100/200 Pa 是 `ENGINEERING_ANALYTIC_STREAMING_LOAD`，绝不是现有 CFD 或微泡模型算得的 surface traction。**禁止把原无血栓 CFD 的压力/速度梯度直接作为本轮血栓应力。**

本轮**至多 4 次正式动力学求解**。一个案例真实地没有断键/没有游离大块也是合格的物理输出；不得因为不碎再加 300/500/1000 Pa，也不新增 pressure-only/vortex/材料扫描。加载失败时最多做针对此失败的求解器/输入软件修复及同 case 重新执行，并保留失败记录；不得修改材料/几何/载荷以追求破裂画面。

## 6. 让结果回答“破坏了什么”，而不只是“断键多少”

所有分量来自真实 surviving-bond graph；原始节点和断键事件不隐藏、不删、不得通过动画着色制造断开。

至少每个保存时刻输出：

- bulk broken fraction、interface broken fraction、first-break location/time；
- 附着到 wall anchor 的 clot 组件及其体积/质量比例；
- 真正不再通过 interface 连接 wall 的 **detached clot-only components**：数量、质量、体积、COM、速度、距原附着区移动距离；
- 末态大主体 fraction、小碎屑（1–2 node unresolved debris）分数、真实 contact 数、动能；
- V1 同样尺度的 `FI_free` 和旧 `>=3 nodes` 工程分类，以保持历史兼容，但不得将其命名为“有明显体积的大碎块”；
- 更直观的文字类目，例如 `NO_FAILURE / ATTACHED_CORE_WITH_LOCAL_DEBRIS_SHEDDING / WHOLE_CLOT_DETACHMENT / DETACHED_RESOLVED_FRAGMENTS / PULVERIZATION / HYBRID / UNCERTAIN`，需保存分类规则和原始连续指标。

**重点：游离碎块尺寸与寿命，严格复用上一轮已写好的 lineage 模块。**

- > =200 µm：文献**报告参考尺度**，在 50 µm cubic nodes 下体积等效球直径满足条件约需 **>=34 nodes**；
  >
- > =300 µm：项目**工程大块可辨识尺度**，约需 **>=114 nodes**（113 nodes 还未达到）；
  >
- 记录每个时刻的同时大块数量/占初始 clot mass 比例、峰值/终值、每块首次出现/末次出现、**完全相同成员集合**连续存在的 episode 以及分裂前后 parent/child lineage；
- 大块尺寸由 `d_eq=(6V/π)^(1/3)` 求体积等效球直径，不冒充最长轴实际直径；如某一块是细长结构，还需保存简单 PCA 三轴/短轴薄片标志；
- 附着主体永远不算“detached large fragment”。碎屑 1–2 个节点为 `NUMERICALLY_UNRESOLVED_DEBRIS`，不能称为真实 50 µm 生物粒子；
- 即使只在中间出现 >=200 µm，而终态为零，也必须报告“短暂出现”，不能只展示末态；
- 大块更多并不意味着治疗更安全：≥200/300 µm 的游离质量、持续时间属于潜在远端栓塞**风险代理**，不是患者栓塞概率。

如实际 fragmented components 会飞出局部锚点范围，可继续算自由运动，但**当前完整 Y 管腔不是固体碰撞边界，且没有流体阻力/真实下游输运**。不能把动画上游离碎块驶入 A/B 支路的画面称为求解出的分支选择或捕获率。

## 7. 初始力学完整性与极少量必要检查

只做：

1. **真实几何**：新 clot 连通、位于原 Y lumen 内、footprint 贴实际 WALL、surface/attachment 不重叠、wall anchors 在 lumen 外；无非法孤点和穿腔 bonds；图与数值保持相同真实坐标；
2. **真实材料点**：总质量=ρΣV，点体积/连续体积/暴露面积差异有据，bulk/interface/wall bond counts 合法；
3. **LAMMPS 装配**：`read_data + run 0` 一次通过、预应力异常则定位，运行器能读取实际接触变体；
4. **C0 零载时间推进**：有限状态、未无故整体破坏；
5. **C1/C2/C3**：载荷节点/N/Pa/面积、当前空间位置、脉宽一致；保存的失效拓扑与 LAMMPS 原生 broken events 一致、没有 NaN/lost-atom/early abort；如果真实出现接触，则简单记录 `rmin/dx`，非要求再做 contact torture test；
6. **图和动画**：每个展示物理量可回溯到实际保存坐标、连接与接触身份；不夸大微粒尺寸/速度。

不做通用刚度有限差分大测试、材料扫、mesh convergence、压力/流场重算、性能 stress test。仅记录真实四例 wall time、steps/s、峰值内存和渲染成本。

## 8. 可视化务必更贴近实际——同时保持科学诚实

本轮用户最重视**看上去是一个血管里的血栓，且能看见真实分裂过程**，而不是仅看彩色散点。

为同一实际坐标生成两类并行视图：

**A. 医学/物理呈现视图**

- 真实 Y 形三端口血管半透明或局部剖切，标明 INLET/A/B，管长和 clot 位置保持实际尺度；
- 展示新 V2 的完整、不规则附壁隆起以及同一位置的原 V1 并排静态对比图（不叠加成两个同时存在的血栓）；
- 动态时给**真实 clot component**按附着/游离身份涂色，使用 voxel / node-based 与真实 surviving graph 一致的 rendering shell；每个 disconnected component **单独重建**，不能用平滑把断口桥接起来；
- 1–2 node unresolved debris 使用弱化 glyph/点显示，不渲染成临床真实微粒；不删点、不人为扩大碎块；
- 全血管视图和局部近景可在一个同步画面中呈现；同物理时刻、相同时间轴、固定相机/同尺度对比三个 loaded cases。不得为“视觉效果好”移动 particles、夸大位移或插入假断裂帧。

**B. 数值审查视图**

- 真实 cubic material nodes、三层 anchors、BULK/INTERFACE/WALL bonds、断开的 bond events、保存的 component IDs；
- 对照实际 WALL surface 和 exposed surface loads；
- 能看到断键过程和数值结构，颜色/统计不应和医学展示矛盾。

原 Y 稳态流场仅作可选背景流线 overlay，始终加显眼提示：

`UNOBSTRUCTED BACKGROUND CFD — NOT CLOTTED FLOW / NOT ACTUAL CLOT TRACTION`。

任何 detached particles 穿过示意的完整 Y vessel wall，都必须解释**完整 Y 管只是可视化背景，没有求解碎块—管壁碰撞**。不要裁掉这些粒子以制造“沿支路运输”的假象。

### 至少输出以下图与动画

1. `y_irregular_clot_v1_vs_v2.png/.pdf`：真实管壁上 V1 和 V2 的实体大小、体积、footprint、厚度、不规则性，带局部剖面和尺度；
2. `y_irregular_clot_material_and_interface.png/.pdf`：新 LAMMPS 点、三层 anchors、三类 bonds、暴露表面、载荷方向；
3. `y_clot_fracture_response_comparison.png/.pdf`：C0/C1/C2/C3 bulk/interface break、attached mass、debris、fragment count 随物理时间变化，图内附真实预算/峰值；
4. `y_clot_detached_fragment_sizes.png/.pdf`：真正 detached >=200 / >=300 µm 的时变数量、峰值质量比例、出现时间/寿命、末态是否存在；未出现则显示 `NOT OBSERVED`，不能凭空制造；
5. `y_clot_initial_fracture_synchronized.mp4` + 简短 GIF：Gaussian100 / Gaussian200 / Opposing100 同时刻对照，至少一个全 Y 视角 + 局部放大真实断裂；
6. 可复用的 `ParaView`/`PyVista` 资产，包括真实血管 mesh、新参考 clot surface、每时刻真实点/bonds/fragment tags，另保留至少一段对主要 case 的 solver-view 小片段或关键帧。

只有一个视频也可以，前提是它能看清全血管上下文及局部裂纹/碎块，且三例可对比；不为了多一个视频花大量 GPU 渲染时间。保存 200–300 左右真实 solver 帧已足够，避免每 timestep dump 和巨型输出。

## 9. 软件复用与可安装依赖

继续扩展现有 `pd_clot/lammps_fragmentation/`，复用先前 `y_bifurcation_sources`、`y_mural_clot_geometry`、`y_wall_anchors`、`y_surface_coordinates`、`lammps_writer`、`current_position_provider`、`fragment_lifetimes`、`visualization` 等实际模块（若名字有差异，以当前仓库为准）。**不要复制第二套 solver 或自行写一个粒子破裂器**。

建议增加/调整最少模块：

- `y_irregular_morphology.py`：确定性几何与 curvature-conforming footprint，清楚记录工程形态；
- `y_fragmentation_cases.py`：使用同一新的材料点/表面积数据生成 4 个参数冻结的 case；
- `y_fragment_visualization.py`：从现有 solver trajectory+lineage 组装 Y context + local fracture shell；
- 必要时最小扩展 `y_surface_coordinates.py` 和已存在的破坏统计结果字段。

允许安装/拉取任何真正必要的包，但必须优先复用已有环境、在 isolated project env 安装，不污染平行微泡项目或现有 LAMMPS 正式版本。官方资源：

- LAMMPS `bond_style bpm/peri`: https://docs.lammps.org/bond_bpm_peri.html
- LAMMPS `pair_style bpm/peri`: https://docs.lammps.org/pair_bpm_peri.html
- LAMMPS `read_data`: https://docs.lammps.org/read_data.html
- LAMMPS `fix addforce`: https://docs.lammps.org/fix_addforce.html
- LAMMPS 官方仓库: https://github.com/lammps/lammps
- PyVista: https://docs.pyvista.org/
- VTK: https://vtk.org/
- Trimesh: https://trimesh.org/
- Meshio: https://github.com/nschloe/meshio
- ParaView: https://www.paraview.org/
- OVITO: https://www.ovito.org/manual/python/
- FFmpeg: https://ffmpeg.org/

使用当前固定 `stable_30Sep2026` + 单独的 `1.0dx` 接触变体，不为这轮下载/编译另一个 LAMMPS；只有原二进制真实缺失才允许恢复它的原本 build。

## 10. 交付结构

建议：

```text
verification/y_bifurcation_irregular_clot_and_initial_fracture_v1/
results/y_bifurcation_irregular_clot_and_initial_fracture_v1/
```

至少交付：

```text
TASK_REQUEST.md
STARTING_GIT_STATE.md
POLICY_SNAPSHOT.md
EXTERNAL_SOURCE_ASSET_REFERENCES.json      # 复用上轮 hash/path，不搬运无关 CFD 大文件
GEOMETRY_V2_CONFIG.json
MORPHOLOGY_DEFINITION.md
V1_V2_GEOMETRY_COMPARISON.json
GEOMETRY_WALL_CLEARANCE_AND_VOLUME.json
Y_CLOT_V2_GEOMETRY.npz
Y_CLOT_V2_SURFACE_NODES_AND_NORMALS.npz
Y_CLOT_V2_BONDS_AND_ANCHORS.npz
Y_CLOT_V2_LAMMPS.data
LAMMPS_RUN0.log
MATERIAL_CONTACT_IDENTITY.json
CASE_CONFIG_FROZEN.json
cases/C0_ZERO_LOAD_DYNAMIC/{request.json,in.lammps,log.lammps,result.json,...}
cases/C1_GAUSSIAN_A100/{request.json,in.lammps,log.lammps,result.json,...}
cases/C2_GAUSSIAN_A200/{request.json,in.lammps,log.lammps,result.json,...}
cases/C3_OPPOSING_A100/{request.json,in.lammps,log.lammps,result.json,...}
FRAGMENT_BREAK_AND_ATTACHMENT_TIMELINES.csv
DETACHED_FRAGMENT_SIZE_AND_LIFETIME.csv
FRAGMENT_LINEAGE_SUMMARY.json
CASE_RESPONSE_COMPARISON.csv/json
PERFORMANCE_REAL_CASES.json
FIGURE_DATA/
VISUALIZATION_QA.json
figures/*.png + *.pdf
animations/y_clot_initial_fracture_synchronized.mp4
animations/y_clot_initial_fracture_preview.gif
FINAL_STATUS.json
REVIEW_Y_BIFURCATION_IRREGULAR_CLOT_AND_INITIAL_FRACTURE_V1.md
evidence.zip
```

轻量 evidence.zip 必须包含关键 config、source hashes、deck、log、原始断键/部分轨迹、统计、图和源码快照。大型原始 CFD 留在已有只读来源，保留路径与 SHA；大 trajectory 留原位并加入哈希清单。图形/视频脚本只读 saved arrays，不允许渲染时偷偷再求解。

## 11. FINAL_STATUS 应至少包括

```text
REAL_Y_MESH_REUSED = PASS/FAIL
OLD_Y_CFD_UNMODIFIED = PASS/FAIL
NEW_IRREGULAR_MURAL_GEOMETRY = PASS/FAIL
V1_V2_GEOMETRY_COMPARE = PASS/FAIL
TRUE_CURVED_WALL_ATTACHMENT = PASS/FAIL
GEOMETRY_CLEARANCE = PASS/FAIL
LAMMPS_CUBIC_V2 = PASS/FAIL
CURVED_INTERFACE_BONDS = PASS/FAIL
LAMMPS_RUN0_V2 = PASS/FAIL
CONTACT_MODEL = MODIFIED_1P0DX_EXPLORATORY_CANDIDATE
MATERIAL_CARD = LIU2024_LPS_SMALL_STRAIN_PLUS_PORCINE_SKIN_INTERFACE_PROXY
ZERO_LOAD_DYNAMIC = PASS/FAIL
GAUSSIAN_100 = EXECUTED/PARTIAL/FAIL
GAUSSIAN_200 = EXECUTED/PARTIAL/FAIL
OPPOSING_100 = EXECUTED/PARTIAL/FAIL
TRUE_BOND_BREAK = OBSERVED/NOT_OBSERVED
TRUE_DETACHED_COMPONENT = OBSERVED/NOT_OBSERVED
DETACHED_GE_200UM = OBSERVED/NOT_OBSERVED
DETACHED_GE_300UM = OBSERVED/NOT_OBSERVED
PERSISTENT_LARGE_DETACHED_FRAGMENT = OBSERVED/NOT_OBSERVED
UNRESOLVED_DEBRIS = TRACKED/FAIL
Y_CONTEXT_FRACTURE_VISUALIZATION = PASS/FAIL
CLOT_PRESENT_CFD_SOLVED = NOT_RUN
CLOT_SURFACE_TRACTION_VALIDATED = NOT_ESTABLISHED
FULL_VESSEL_FRAGMENT_COLLISION = NOT_MODELED
FLUID_FRAGMENT_TRANSPORT = NOT_MODELED
BUBBLE_STREAMING_CONNECTED = NOT_YET
CURVED_INTERFACE_STRENGTH_VALIDATED = NOT_ESTABLISHED
BIOLOGICAL_QUANTITATIVE_VALIDATION = NOT_ESTABLISHED
PATIENT_SPECIFIC = NOT_ESTABLISHED
B1_B2 = TERMINATED
PD2_FRACTURE = STOPPED
PRIMARY_USE = EXPLORATORY_Y_VESSEL_MURAL_CLOT_FRAGMENTATION
```

`TRUE_BOND_BREAK=NOT_OBSERVED` 是一种允许的科学结果，不是请求失败；不能因用户想看碎裂而调高载荷。`LAMMPS_RUN0`、`ZERO_LOAD_DYNAMIC`、几何合法性则属于硬软件 gate。

## 12. 用中文写最终 REVIEW，直接回答下列问题

1. 实际来源 Y 三维网格和无血栓 CFD 的 hash 是否仍未变化？是否避免替换成假管道？
2. 新血栓实际 a/b/H、体积、轴向/周向尺寸和外观形态与 V1 相比变化多少？不规则项具体是什么？
3. 真实附着壁面最大 gap、越界/管壁交叉/入口分叉 clearance 如何？
4. 新 V2 LAMMPS 有多少 clot nodes、三层 anchors、BULK/INTERFACE/WALL bonds？连续与积分体积是否差得较大？
5. `run 0` 与零载动态是否通过？若异常原因何在？
6. Gaussian100/Gaussian200/Opposing100 真实发生了什么？断键、附着质量减少、游离质量与碎屑差多少？不能保证三例都碎。
7. 是否出现 ≥200/≥300 µm 的真正游离块？最大同时数量、最大尺寸、峰值质量、持续时间、最后有没有留下？
8. 三种工况谁更多地释放碎屑、谁更多保留附着主体？同预算比较有哪些空间局部峰值和形状上的限制？
9. MP4 是否能清楚看到真实 Y 管上下文、形态稍大更不规则的 clot、并且只有真实 component 断开才可视化为碎块？是否完整保留未解析颗粒？
10. 真实计算墙钟时间、后处理时间、内存分别是多少？是否保持低成本？
11. 当前无血栓 CFD、猪皮 proxy interface、无完整血管碰撞/无输运以及没有真实 bubble traction 分别限制了哪些解释？
12. 在用户下一阶段若要获得真实流体载荷时，缺哪一步（含血栓稳态 CFD 或 bubble streaming field）？**只报告后续接口缺口，不自动执行**。

## 13. 本轮完成与停止条件

只有**一个冻结 V2 几何**、**一份 LAMMPS assembly**、**四个声明的 C0–C3 动态 case**；做完并导出图/动画/证据之后立即停止等待人工审查。不进行额外 amplitude bracket、shape search、stability paper、PD2/FEM review、GPU tuning、CFD 重算、bubble 集群计算、fragment transport。若 C0 或几何 gate 真正失败，保存证据并停止 loaded cases，不通过修改材料和载荷回避错误。

**本轮对最终控制路线的价值是：首次在实际 Y 形血管上下文中，有一个较自然的不规则附壁血栓，可以低成本计算人工 shear 下的真实断键、脱粘与碎屑运动；它仍只是探索性固体模型，不是已验证的“真实背景血流/微泡让血栓破裂”。**