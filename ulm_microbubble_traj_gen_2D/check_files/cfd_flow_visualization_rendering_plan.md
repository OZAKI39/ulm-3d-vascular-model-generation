# CFD 风格流场与流线可视化渲染计划

## 1. 文档目的

本文档汇总当前项目的 CFD 风格流场可视化目标、工具选型、数据转换、流线生成、渲染图层、交互交付和验收标准。

截至 2026-07-11，本文档定义的二维/2.5D 首版已经接入运行链：使用 PyVista/VTK 代替原 Plotly 自动流场面板，但不改变现有 PhiFlow/SciPy 求解、Taichi 粒子推进及其正确注释。本文档同时保留旧实现的问题分析，作为迁移依据。

---

## 2. 需求理解

### 2.1 目标视觉效果

期望效果接近专业 CFD 后处理软件生成的流线图，而不是在二维热图上叠加少量普通折线。目标画面应包含：

- 连续、平滑、根根分明的流线；
- 流线具有稳定半径、抗锯齿、光照和空间层次感；
- 沿流线按照真实速度大小连续着色；
- 所有正式流线必须从 root 入口连续积分到一个有效出口，中途不得断裂；
- 多条 root 种子流线共同覆盖主干、分叉、所有终端和细小支路；
- 速度、压力等标量场保持独立颜色条；
- 固体或管腔外区域严格透明；
- 管壁或血管几何作为灰色、半透明的上下文存在；
- 可旋转、缩放、平移，并能查询局部数值；
- 初始场和最终收敛场使用一致的可视化规则，便于比较。

### 2.2 参考图与当前数据的物理差异

参考图属于真正的三维 CFD 体积流场：流体存在于三维物体周围，流线可分布在不同深度并形成空间涡旋。

当前项目保存的是二维 X-Z 平面管腔流场：

```text
velocity.shape = (nx, nz, 2)
velocity[..., 0] = vx
velocity[..., 1] = vz
fluid domain = lumen_mask
```

因此必须区分两种目标。

#### 目标 A：忠实的二维 / 2.5D CFD 风格

- 保留当前二维物理结果；
- 将二维流线路径渲染成具有半径和光照的 Stream Tube；
- 使用正交视图或轻微倾斜相机产生空间层次；
- 流线路径仍严格位于同一个物理平面；
- 不虚构不存在的 Y 向速度和三维涡旋。

这是当前阶段的推荐目标。

### 2.3 Root-to-outlet 连续性是强制条件

正式 CFD 渲染中的每一条流线都必须满足：

```text
root inlet -> continuous streamline -> one valid outlet
```

不允许把以下路径作为正式渲染流线：

- 从血管中部或分叉处开始的路径；
- 尚未到达出口便因显示密度、路径占位或去重而截断的路径；
- 因最大步数、最大长度、低速阈值、NaN、插值失败或碰壁提前停止的路径；
- 由多段不连续折线拼接而成、但中间没有真实积分连续性的路径；
- 起点不属于 root 入口或终点不属于已标记出口的路径。

一条数学流线在分叉处不会自行分裂。因此“从 root 覆盖所有出口”必须由多条独立流线共同完成，每条流线对应一条连续的 root-to-one-outlet 路径，而不是让一条线在分叉处分成多条。

#### 目标 B：真正的三维 CFD 风格

若需要真实的不同深度流线、空间涡旋和三维速度结构，则必须：

1. 构建三维血管流体网格；
2. 求解三维速度 `(vx, vy, vz)` 和压力；
3. 在三维流体体积内积分流线。

把二维结果简单拉伸或复制到三维只能产生视觉厚度，不能形成真实三维流动。

---

## 3. 旧实现基线与当前实现

### 3.1 已移除的旧调用链

```text
generate_microbubble_trajectories.py
    -> utils.runner.run_generation()
    -> 求解初始场与最终收敛场
    -> vis_utils.plotly_flow.render_interactive_flow_fields()
    -> Matplotlib axes.streamplot() 计算路径
    -> Plotly Scattergl 显示流线
```

### 3.2 当前 PyVista/VTK 调用链

```text
generate_microbubble_trajectories.py
    -> utils.runner.run_generation()
    -> vis_utils.pyvista_flow.render_cfd_flow_fields()
    -> vis_utils.vtk_flow_grid.build_vtk_stage_grid()
    -> lumen threshold + vtkImageDataLIC2D
    -> vis_utils.vtk_streamlines.trace_root_to_outlets()
    -> outlet backward 诊断定位 root 射击点
    -> root forward RK4 重新积分并做连续性/覆盖验收
    -> PyVista 1×2 标量面、Stream Tube、Glyph
    -> VTI / VTP / CSV / PNG / 离线 HTML
```

### 3.3 旧流线计算方式

当前流线使用完整分辨率速度场，没有进行流线积分网格降采样。核心参数为：

```text
density = 4.0
minlength = 0.03
integration_direction = both
broken_streamlines = True
```

Matplotlib 首先在整个矩形域上建立 `120 x 120` 的 StreamMask，再自动选择候选种子。每个 StreamMask 网格最多允许一条流线通过，后续流线靠近已有流线时可能提前终止。

### 3.4 旧实现真实结果统计

当前典型结果：

- 原始网格约 `1555 x 1337`；
- 管腔约占整个矩形域 `5.30%`；
- 血管段数量为 69；
- 初始场约生成 43 条流线；
- 最终场约生成 44 条流线；
- 最终流线触达约 67/69 个 `vessel_id`。

当前问题不是速度为零，也不是 Plotly 丢弃路径，而是：

1. 全矩形均匀播种与稀疏细血管几何不匹配；
2. StreamMask 的全局占位规则会阻止共享主干的后续流线；
3. `broken_streamlines=True` 会在已有路径附近提前截断；
4. `minlength=0.03` 会删除约 76–89 μm 以下的短路径；
5. Plotly 只放大预先计算好的静态路径，缩放后不会重新播种；
6. 白色 `1.25 px` 细线在全域双子图中视觉存在感较弱。

因此，后续方案不能只提高 Matplotlib `density`，必须更换流线积分和播种体系。

---

## 4. 工具选型

### 4.1 推荐工具栈

| 工具 | 主要职责 | 推荐程度 |
|---|---|---:|
| ParaView | CFD 正确性检查、交互调参、LIC、出版级渲染 | 首选 |
| PyVista / VTK | 集成进 Python、流线积分、Stream Tube、数据导出 | 首选 |
| PyVista / VTK + trame | 完整网页应用、图层控制、动态探针 | 条件推荐 |
| Plotly | 旧自动流场面板 | 已从成功运行链移除 |
| Tecplot 360 | 商业 CFD 后处理与出版级出图 | 有许可证时备选 |
| VisIt | 超大规模 HPC 数据和并行积分 | 大规模场景备选 |

### 4.2 ParaView 的职责

ParaView 用于先建立高质量视觉基准并验证物理正确性：

- `Threshold`：根据 `lumen_mask` 移除固体单元；
- `Surface LIC`：以连续纹理表达整个流体域的局部方向；
- `Stream Tracer With Custom Source`：使用任意点集作为流线种子；
- `Tube`：将一维流线转换成具有半径的管状几何；
- `Glyph`：沿流线放置方向箭头或圆锥；
- 标量着色：按速度、压力或积分时间着色；
- 正交相机、光照、抗锯齿、透明表面和高分辨率截图。

ParaView 阶段的主要目的不是最终自动交付，而是确定：

- 合理的种子数量；
- Stream Tube 半径；
- LIC 强度与透明度；
- 速度色阶；
- 管壁透明度；
- 相机角度；
- 全局与局部视图的视觉密度。

### 4.3 PyVista / VTK 的职责

PyVista / VTK 用于将经过 ParaView 验证的流程自动化：

- 从 NumPy 构建规则网格或结构化网格；
- blank / 移除固体单元；
- 根据自定义种子独立积分流线；
- 使用 RK2、RK4 或 RK45；
- 保存 SeedId、积分时间、速度和终止原因；
- 生成 Stream Tube、方向 Glyph 和 LIC；
- 输出 `.vti`、`.vtp`、`.vtkhdf` 或 Web 场景。

推荐优先使用自定义种子的 `vtkStreamTracer` / `streamlines_from_source()`。`streamlines_evenly_spaced_2D()` 可用于美学填充，但不能单独保证所有细分支均被覆盖。

### 4.4 trame 的职责

若最终交付允许运行 Python 服务，trame 用于：

- 在浏览器中显示 VTK 场景；
- 客户端或服务器端渲染；
- 图层开关；
- 点击和悬停探针；
- 动态相机与多视图；
- 根据缩放级别切换流线密度；
- 使用 GPU 处理更复杂的 Stream Tube 和 LIC。

trame 不是双击即可打开的单文件 HTML，需要应用或服务部署。

### 4.5 Plotly 的移除范围

`generate_microbubble_trajectories.py` 成功路径中的旧 Plotly 流场模块和对应测试已经移除。当前自动可视化不再使用 Matplotlib `streamplot()` 或 Plotly `Scattergl` 计算/显示流线。项目中独立的 Matplotlib 粒子动画、孔洞与窄腔诊断工具不属于这条旧流场链，因此保持不变。

---

## 5. VTK 数据模型与导出

### 5.1 二维场坐标映射

PyVista 的二维等间距流线滤镜要求数据位于 VTK 的 X-Y 平面，因此将当前 X-Z 平面映射为：

```text
VTK X = physical X
VTK Y = physical Z
VTK Z = 0
```

速度映射为：

```text
VTK velocity = (vx, vz, 0)
```

该映射只改变 VTK 坐标命名，不改变物理 X/Z 数值和比例。

### 5.2 推荐保存的数组

```text
velocity                 point vector, [μm/s]
speed                    point scalar, [μm/s]
pressure                 point scalar, projection potential [a.u.]
lumen_mask               point/cell scalar, bool or uint8
wall_mask                point/cell scalar
vessel_id                point/cell scalar
distance_to_wall_um      point scalar
wall_shear_stress_pa     point scalar
wall_shear_display_mask  cell scalar, closed-wall display band
inlet_label              point/cell scalar
outlet_label             point/cell scalar
```

初始场没有独立求解的初始压力。初始压力面板只能显示明确标注的零参考场，不能误用最终压力。`wall_shear.py` 同样只在最终接受速度场上计算一次二维面内壁面剪切应力代理，因此 initial VTI 的同名数组全部为 NaN，不能制作或解释为 initial WSS。

### 5.3 推荐文件格式

| 文件 | 内容 |
|---|---|
| `initial_flow_field.vti` | 初始速度、速度大小、掩膜和相关属性 |
| `final_flow_field.vti` | 最终速度、压力、速度大小和相关属性 |
| `initial_streamlines.vtp` | 初始场通过连续性检查的正式 root-to-outlet 流线 |
| `final_streamlines.vtp` | 最终场通过连续性检查的正式 root-to-outlet 流线 |
| `initial_streamlines_diagnostic.vtp` | 初始场不完整或诊断路径，不参与正式渲染 |
| `final_streamlines_diagnostic.vtp` | 最终场不完整或诊断路径，不参与正式渲染 |
| `initial_root_to_outlet_continuity.csv` | 初始场每条路径的入口、出口、完整性、终止原因和分支序列 |
| `final_root_to_outlet_continuity.csv` | 最终场每条路径的入口、出口、完整性、终止原因和分支序列 |
| `final_wall_shear_stress_cfd.html` | 最终二维面内壁面剪切应力代理的离线交互场景 |
| `final_wall_shear_stress_cfd.png` | 同一 WSS 场景的静态预览 |
| `flow_visualization.pvsm` | ParaView 管线与视觉参数 |

规则均匀网格优先使用 `vtkImageData` / `.vti`。若后续提取真实管腔单元或使用非规则网格，可改用 `vtkUnstructuredGrid` 或 `.vtkhdf`。

### 5.4 固体处理

标量渲染和流线积分应使用不同但一致的掩膜策略：

- 标量场：`~lumen_mask` 为无效或透明；
- LIC：固体区 alpha 为 0；
- 流线积分：最好实际 blank 或删除固体单元；
- 低速无效区：不生成 LIC 噪声，也不作为流线种子；
- `wall_mask` 是近壁流体带，不是完整固体域，不能直接从流体中删除。

实际提取管腔流体网格比只在完整矩形外部填 NaN 更稳健，因为积分器离开真实流体单元时可自然终止。

### 5.5 壁面剪切应力代理的显示约束

`wall_shear.py` 计算的是 `abs(mu * t·[(grad u + grad u^T)n])`，单位为 Pa。它是二维面内 WSS proxy，不是乘过面积的剪切力 N，也不是完整三维无滑移壁面重建。正式 WSS 场景只显示 `lumen_mask & wall_mask` 的近壁流体带，并排除 `inlet_label`、`outlet_label` 对应的开放端面；管腔内部、固体和开放面只作透明或浅灰上下文。色谱使用非负顺序色谱，显示上限取有效壁面值 P99.5 以抑制栅格尖角峰值，但 VTI 与 hover 始终保留原始 Pa 数值。

---

## 6. Root-to-outlet 流线种子方案

正式流线只允许从 root 入口播种。血管中部、分叉和出口种子只能用于诊断、反向追踪或寻找入口射击位置，不能直接作为最终显示路径的起点。

### 6.1 Root 入口种子面

首先根据求解器已有的 root vessel、`inlet_label`、入口法向和开放边界面构造入口种子面。

入口种子必须：

1. 位于 root 入口边界向管腔内偏移少量距离的位置，避免第一步立即离开计算域；
2. 严格位于 `lumen_mask` 内；
3. 具有有限、非零速度；
4. 速度方向指向管腔内部；
5. 保存入口面坐标、局部速度、SeedId 和初始流量权重。

种子分布可使用：

- 等物理距离播种；
- 按入口直径自适应播种；
- 等流量播种，使相邻种子近似代表相同流量；
- 中心种子加近壁分层种子，以同时覆盖中心高速区和近壁低速区。

入口种子数量不能由全矩形 `density` 决定，而应由 root 截面宽度、目标流线间距、出口数量和出口覆盖结果共同决定。

### 6.2 正向 root-to-outlet 积分

每个入口种子使用独立、仅向前的 VTK 流线积分：

```text
integration_direction = forward
```

积分过程中：

- 不使用全局 StreamMask 阻止后续路径共享主干；
- 允许多条流线在 root 主干中平行或接近；
- 允许不同出口路径共享部分上游主干；
- 在分叉处完全依据连续速度场自然选择子分支；
- 不因靠近已有流线而提前终止；
- 只有到达已标记出口时，该路径才被视为成功的正式流线。

### 6.3 出口覆盖与自适应补种

第一次正向积分完成后，按 `outlet_label` 对成功路径分组，并统计：

```text
outlet_id -> number of complete root-to-outlet streamlines
```

若某个出口没有完整流线到达，则依次执行：

1. 增加 root 入口种子密度；
2. 在能进入该出口吸引域的入口区间进行局部细分；
3. 从缺失出口向上游做反向诊断积分；
4. 根据反向路径与 root 入口的交点确定新的入口种子；
5. 从该 root 种子重新执行 forward 积分；
6. 只有重新得到连续的 root-to-outlet 正向路径后，才加入正式渲染集合。

反向积分路径可以帮助寻找 root 入口射击位置，但最终保存和显示的路径必须按 root 到出口方向排序，并通过连续性验证。

### 6.4 拓扑保证流线

拓扑保证的目标从“每个分支有一条局部流线”升级为：

```text
每个有效出口至少有一条从 root 连续到达的流线
所有 69 个 vessel_id 均被至少一条完整 root-to-outlet 路径触达
```

在树状网络中，只要每个出口都有一条完整路径，所有通向这些出口的父级血管应同时被覆盖。仍需用 `vessel_id` 序列进行显式验证，防止栅格映射或分叉插值造成遗漏。

当前实现优先使用求解器 `BoundaryFluxFields.outlet_ids` 建立 `outlet_label -> terminal vessel_id` 权威映射，不再依赖出口附近栅格标签的多数投票。每条候选路径只按非 junction 弧长统计血管 ID，并要求所有具有可观测中心流路的祖先段按 root→outlet 顺序完整出现；兄弟分支、漏段、乱序或回跳均判为诊断路径。只有中心流路完全被 `junction_core_mask` 吞没的极短祖先段，才允许由已经命中该出口的权威父子链补记覆盖。

### 6.5 宽度自适应填充种子

美学填充流线同样必须从 root 入口开始，不能从中间血管直接播种。

入口种子密度根据下游需求自适应：

- 对通向细小低流量出口的入口吸引域，至少保留一个种子；
- 对宽主干和高流量出口，增加平行 root-to-outlet 流线；
- 对近壁区域，可增加低速层种子；
- 对已经具有足够路径密度的出口，不再继续增加种子；
- 每条新增路径仍必须完整到达某个出口，否则只进入诊断集合。

### 6.6 完整路径去重

去重只能作用于已经通过 root-to-outlet 连续性检查的完整路径：

1. 先保留每个出口至少一条拓扑保证路径；
2. 再比较相同出口路径之间的空间距离和共享长度；
3. 删除几乎重合且不增加覆盖率的路径；
4. 不允许截取部分路径作为替代；
5. 不允许删除某个出口的最后一条完整路径；
6. 不允许删除某个 `vessel_id` 的最后一条完整覆盖路径。

### 6.7 诊断路径与正式路径分离

中部、分叉或出口播种产生的路径必须标记为 diagnostic，只用于：

- 检查局部速度方向；
- 定位插值或分叉问题；
- 反向寻找入口种子；
- 分析某个出口为何无法从 root 到达。

这些路径不能进入正式 Stream Tube 渲染，也不能计入 root-to-outlet 覆盖率。

---

## 7. 流线积分策略

### 7.1 推荐积分器

当前正式路径使用以网格间距定义的固定步长 RK4，使初始场与最终场采用完全一致、可重复审计的积分规则。RK45 可保留为后续精度对照，但不能改变开放面命中和 root-to-outlet 验收语义。

### 7.2 初始参数原则

- 步长以 `grid_spacing_um` 的一定比例定义；
- 最小、最大步长均使用物理单位；
- 最大积分长度必须大于最长 root-to-outlet 血管路径，并使用物理长度而不是归一化轴长度；
- 正式 root 种子统一使用 forward；
- 中部或出口诊断种子可使用 backward 或 both，但其原始结果不进入正式渲染；
- 到达已标记出口时正常终止并记为成功；
- 离开流体网格但未命中出口时记为失败；
- 速度低于物理阈值时终止并记为失败；
- 达到最大步数或最大长度时终止并记为失败；
- 任何失败路径不得作为正式 root-to-outlet 流线显示；
- 保存每条流线的终止原因。

### 7.3 必须保存的路径属性

```text
SeedId
path_id
root_inlet_id
destination_outlet_id
source_vessel_id
IntegrationTime
ReasonForTermination
semantic_termination
is_complete_root_to_outlet
velocity vector
speed
pressure
sampled vessel_id
```

其中 `ReasonForTermination` 保存 VTK 原始终止原因，`semantic_termination` 保存项目判定结果。VTK 可能把流线穿过开放出口记为离开计算域；只有路径末端命中 `outlet_label`、穿过对应开放边界面且方向向外时，项目才将其语义分类为 `OUTLET_REACHED`。

这些属性用于检查入口来源、目标出口、分支覆盖、路径连续性、流线穿墙、异常终止和颜色映射。

### 7.4 路径连续性检查

每条正式路径必须通过以下逐项检查：

1. 第一个点属于指定 root 入口内侧区域；
2. 最后一个点属于某个 `outlet_label` 区域；
3. 路径点顺序为 root 到 outlet；
4. 相邻路径点物理距离不超过积分步长容差；
5. 所有路径点和插值线段均位于流体域内；
6. 采样的 `vessel_id` 序列在血管图上构成连续父子路径；
7. `semantic_termination == OUTLET_REACHED`；VTK 原始终止原因可以是离开计算域，但末端必须命中已标记开放出口；
8. 路径中不存在 NaN、Inf 或人为拼接跳点；
9. 路径在转换成 Tube 前仍保持单一连续 PolyLine cell。

任一检查失败，该路径只能进入诊断输出，不能进入正式 CFD 渲染。

---

## 8. 渲染图层

### 8.1 图层一：管壁或血管几何

- 使用浅灰色或半透明材质；
- 避免完全遮挡内部流线；
- 2D 方案可显示管腔边界和轻微厚度；
- 2.5D 方案可叠加现有 SWC 圆柱血管几何；
- 若使用三维几何，必须明确速度场仍是二维平面场，除非已完成三维求解。

### 8.2 图层二：标量背景

支持切换：

- 速度大小；
- 压力投影势；
- 涡量；
- 壁面切应力；
- 其他诊断量。

要求：

- 固体透明；
- 保留颜色条；
- 使用真实物理坐标；
- 初始和最终阶段可以共享或独立色阶，具体由比较目标决定。

### 8.3 图层三：LIC

LIC 用于让每个有效流体像素都显示局部方向。

推荐样式：

- 灰白或黑白纹理；
- 作为半透明层叠加在标量颜色上；
- 固体区域透明；
- 极低速度区域减弱或关闭；
- 纹理长度与网格分辨率和局部血管宽度匹配。

LIC 只能表达局部轴向，不能区分 `v` 和 `-v`，所以不能取代方向箭头或代表性流线。

### 8.4 图层四：Stream Tube

将一维流线转换成具有半径的管状几何：

- 只渲染通过 root-to-outlet 连续性检查的正式路径；
- 每个 Tube 必须对应一个从 root 到单一出口的连续 PolyLine；
- 不允许在分叉处用不连续线段视觉拼接成“完整”路径；
- 初始半径建议从 `0.2–0.5 x grid_spacing_um` 范围调试；
- 沿路径按真实速度大小连续着色；
- 使用平滑法线；
- 开启抗锯齿；
- 使用适度环境光、漫反射和高光；
- 避免半径过大导致分叉处严重遮挡；
- 颜色条单位保持 `[μm/s]`。

若全局视图仍拥挤，可采用线宽分级或只把拓扑保证流线转换成 Tube，其余填充流线使用较细普通线。

### 8.5 图层五：方向标记

可选择：

- 沿流线按固定物理距离放置 Cone Glyph；
- 使用小型箭头；
- 让少量粒子沿路径运动；
- 在入口和分叉处重点显示方向标记。

方向标记数量必须独立降采样，不能在每个路径点都放置 Glyph。

---

## 9. 相机、光照与颜色

### 9.1 二维忠实视图

- 使用正交相机；
- 视线垂直于流场平面；
- X/Z 比例固定为 1:1；
- 适合定量比较和压力/速度读数。

### 9.2 2.5D 展示视图

- 相机略微倾斜；
- Stream Tube 保留轻微空间厚度；
- 管壁使用半透明；
- 不改变路径的实际平面坐标；
- 适合展示 CFD 风格和空间层次。

### 9.3 颜色策略

- 速度流线可使用 Turbo、Jet 或其他高对比连续色谱；
- 压力背景优先使用以零为中心的冷暖色谱；
- 速度和压力颜色条必须分离；
- 避免用光照改变标量颜色到难以读取的程度；
- 初始/最终比较图若强调绝对差异，应共享同一色阶；
- 若强调各自内部结构，可使用独立稳健百分位色阶，但必须明确标注。

---

## 10. 多尺度与 LOD

完整物理范围约为 3 mm，而典型血管直径约为 15 μm。在双子图全局视图中，细血管只有几个屏幕像素，任何工具都无法同时显示其内部全部流线细节。

建议设计三档 Level of Detail：

### 全局视图

- LIC；
- 每个出口至少一条完整 root-to-outlet 代表性流线；
- 少量方向箭头；
- 管壁轮廓；
- 避免大量 Stream Tube 遮挡。

### 中等缩放

- 显示宽度自适应填充流线；
- 增加方向 Glyph；
- 显示分支 ID 或入口/出口标记。

### 局部放大

- 更密集的流线；
- 速度矢量；
- 精确速度、压力和 `vessel_id`；
- 分叉局部诊断；
- 流线终止原因。

trame 应用可根据相机范围动态切换。若必须使用离线 HTML，则预计算三套路径几何，并根据缩放事件切换可见性。

---

## 11. 交互与交付形式

### 11.1 ParaView 分析产物

- `.vti/.vtkhdf` 场数据；
- `.vtp` 流线；
- `.pvsm` 状态文件；
- 高分辨率 PNG；
- 旋转或时间序列视频。

### 11.2 Python 自动产物

- PyVista/VTK 预计算场景；
- 流线统计 CSV/JSON；
- root-to-outlet 连续性和出口覆盖报告；
- `.vtp` 路径文件；
- 自动截图和视频。

### 11.3 完整网页应用

使用 PyVista/VTK + trame：

- 图层开关；
- 初始/最终切换；
- 速度/压力切换；
- 点击或悬停查询；
- 相机同步；
- 动态 LOD；
- 服务端或客户端渲染。

当前首版已经提供数值 hover 探针入口：

```text
python -m ulm_microbubble_traj_gen.vis_utils.trame_flow_viewer \
    --result-dir <results/时间戳目录> \
    --stage final \
    --open-browser
```

该页面直接读取同目录 VTI/VTP，鼠标悬停时查询原始 cell-centered VTI，并显示物理 X/Z、`vx`、`vz`、`|v|`、压力语义、压力值和 `vessel_id`。它不会从 Tube 的插值颜色反推数值。

最终 WSS 场景使用同一个查看器入口：

```text
python -m ulm_microbubble_traj_gen.vis_utils.trame_flow_viewer \
    --result-dir <results/时间戳目录> \
    --stage final \
    --view wall-shear \
    --open-browser
```

WSS hover 额外显示原始 `wall_shear_stress_pa` 和 `wall_shear_semantics`。鼠标落在管腔内部或入口/出口开放面时不会伪报壁面 WSS，而是明确提示该位置不属于 closed-wall display band。

### 11.4 单文件离线 HTML

若必须双击打开：

1. VTK/PyVista 离线计算 LIC、正式流线、Tube、Glyph 和标量几何；
2. 使用 PyVista `export_html()` 输出 trame/VTK.js OfflineLocalView；
3. HTML 保留旋转、缩放、平移和颜色条；
4. 不在浏览器中实时重新积分 CFD 流线；
5. OfflineLocalView 自身不执行数值 hover；精确 hover 使用同目录提供的 `trame_flow_viewer`，也可在 ParaView 中读取 VTI/VTP。

离线 HTML 适合分享和几何交互，但无法完整替代 ParaView/trame 服务端应用的动态探针与 LOD。

---

## 12. 分阶段实施计划

### 阶段 0：确认物理范围

- 确认首版目标是二维/2.5D，而不是真正三维 CFD；
- 确认最终交付优先级：桌面分析、离线 HTML 或 Web 应用；
- 确认初始/最终颜色范围是否共享；
- 确认压力仅为 projection potential `[a.u.]`。

### 阶段 1：VTK 导出与 ParaView 原型

1. 导出初始和最终 `.vti`；
2. 在 ParaView 中 Threshold 管腔；
3. 验证速度和压力标量；
4. 配置 Surface LIC；
5. 只从 root 入口设置少量自定义种子，验证流线可连续到达出口且不穿墙；
6. 添加 Tube 和 Glyph；
7. 确定视觉参数；
8. 保存 `.pvsm`。

这一阶段用于快速验证方向和视觉目标，避免过早编写大量自动化代码。

### 阶段 2：Root-to-outlet 拓扑保证流线

1. 根据 root vessel、`inlet_label` 和入口法向构造入口种子面；
2. 所有正式种子只放在 root 入口内侧；
3. 使用 VTK forward 独立积分，不使用全局路径占位互斥；
4. 根据 `outlet_label` 记录每条路径到达的出口；
5. 对缺失出口提高入口种子密度并进行局部细分；
6. 必要时从缺失出口反向诊断，以确定新的 root 入口射击种子；
7. 检查每条正式路径起于 root、终于出口且中间连续；
8. 检查每个出口至少有一条完整路径；
9. 检查 69/69 分支被完整 root-to-outlet 路径覆盖；
10. 检查所有路径及线段均位于管腔内；
11. 将不完整路径分离到诊断输出；
12. 输出正式 `.vtp`、诊断 `.vtp` 和连续性报告。

### 阶段 3：美学填充与 CFD 风格渲染

1. 在 root 入口按出口覆盖、流量和下游血管宽度增加种子；
2. 仅对完整 root-to-outlet 路径实施积分后去重；
3. 添加 Stream Tube；
4. 沿路径按速度着色；
5. 添加 Glyph；
6. 添加 LIC；
7. 完成全局、中等、局部三档 LOD；
8. 确定相机、光照、背景和颜色条。

### 阶段 4：Python 自动化

- 将 ParaView 验证过的参数迁移到 PyVista/VTK；
- 接入 `run_generation()` 成功路径；
- 自动生成 VTK 场、流线、截图和统计；
- 对初始和最终场使用一致流程；
- 增加单元测试与真实场回归测试。

### 阶段 5：最终交付

按实际需求选择：

- ParaView 高质量分析视图；
- PyVista 本地交互视图；
- trame Web 应用；
- PyVista/VTK 预计算场景的离线 HTML。

---

## 13. 验收标准

### 13.1 物理与几何正确性

- [ ] 每条正式流线的第一个点位于 root 入口内侧；
- [ ] 每条正式流线的最后一个点位于有效出口；
- [ ] 每条正式流线从 root 到出口全程连续；
- [ ] 正式流线只使用 forward 物理方向显示；
- [ ] 中部、分叉或出口诊断种子生成的原始路径不进入正式渲染；
- [ ] 所有流线点均位于 `lumen_mask` 内；
- [ ] 流线线段不穿过固体或管壁；
- [ ] X/Z 坐标和间距与求解网格一致；
- [ ] X/Z 显示比例始终为 1:1；
- [ ] 初始流线来自 `initial_velocity_xz_um_s`；
- [ ] 最终流线来自收敛后的 `velocity_xz_um_s`；
- [ ] 最终压力来自求解器 pressure；
- [ ] 初始压力只显示明确标注的零参考场。

### 13.2 分支与空间覆盖

- [ ] root 入口具有足够的有效种子；
- [ ] 每个有效出口至少有一条完整 root-to-outlet 流线；
- [ ] 69/69 `vessel_id` 均被完整 root-to-outlet 路径覆盖；
- [ ] 每个分叉的父段和子段均位于至少一条完整路径中；
- [ ] 通往细血管出口的路径至少保留一条；
- [ ] 通往宽血管/高流量出口的路径具有合理数量的平行流线；
- [ ] 不完整路径数量、失败原因和对应出口均有报告；
- [ ] 流线覆盖率达到预设物理距离阈值。

### 13.3 渲染效果

- [ ] 流线根根分明且无明显锯齿；
- [ ] Stream Tube 半径不会严重遮挡分叉；
- [ ] 流线颜色与真实速度一致；
- [ ] LIC 覆盖全部有效非零速度管腔区域；
- [ ] LIC、标量颜色与 Stream Tube 之间对比清晰；
- [ ] 方向箭头数量合理；
- [ ] 管壁上下文清晰但不过度遮挡；
- [ ] 颜色条单位和含义正确。

### 13.4 交互与性能

- [ ] 平移、缩放和旋转流畅；
- [ ] 放大后可查询 `vx`、`vz`、`|v|`、pressure、`vessel_id`；
- [ ] 初始/最终视图可切换或并列比较；
- [ ] 全局视图不会因过多 Tube/Glyph 卡顿；
- [ ] LOD 切换不会导致明显跳变；
- [ ] 输出文件大小和加载时间符合交付要求。

---

## 14. 风险与权衡

### 14.1 流线越多不等于越清晰

过密 Stream Tube 会产生严重遮挡。全域方向应主要由 LIC 表达，流线负责展示拓扑路径和方向。

### 14.2 LIC 不能表达正反方向

LIC 只能显示局部轴向，因此必须配合箭头、Cone Glyph、移动粒子或带方向标记的代表性流线。

### 14.3 二维结果不能产生真实三维涡旋

2.5D 只是一种渲染表现。所有出平面运动必须来自真实三维求解，不能由可视化层推断或虚构。

### 14.4 单文件 HTML 与完整 CFD 交互存在冲突

单文件 HTML 便于分享，但不适合实时 VTK 积分、动态 LIC 和复杂探针。最高质量 Web 体验通常需要 trame 服务。

### 14.5 全局统一流线间距不适合直径跨度大的血管

宽血管和细血管必须在 root 入口使用出口感知、宽度自适应的播种策略。不能通过在中间分支直接增加正式种子来伪造分支覆盖。

### 14.6 并非所有 root 种子都能到达出口

部分入口种子可能因数值插值、近壁低速、分叉速度误差或网格离散而提前终止。必须将这些路径标记为失败，并通过入口局部细分或出口反向诊断寻找新的 root 种子，不能将失败路径截断后继续显示。

### 14.7 一条流线不能在分叉处分裂

root 到多个出口的覆盖必须由多条独立、连续的流线共同完成。任何在分叉处复制路径或把一条 PolyLine 人为拆成多条子路径的做法，都不能被解释为真实流线。

---

## 15. 推荐结论

当前项目的推荐路线为：

```text
ParaView 建立并验证 CFD 视觉基准
    -> PyVista / VTK 自动生成掩膜流体网格、拓扑流线、Stream Tube 和 LIC
    -> 按交付需求选择 PyVista 离线 HTML、trame Web 应用或 ParaView
```

首版应实现忠实的二维/2.5D CFD 风格，不虚构三维速度。真正解决当前“流线少且不够直观”问题的关键不是单纯增加 Matplotlib `density`，而是：

1. 使用 LIC 覆盖全部有效流体区域；
2. 只在 root 入口进行正式播种，并按出口覆盖自适应细分入口种子；
3. 使用 VTK forward 独立积分每条 root-to-outlet 流线；
4. 对缺失出口使用反向诊断定位新的 root 入口种子；
5. 拒绝所有未到达出口或不连续的正式路径；
6. 使用 Stream Tube、速度着色、Glyph、光照和抗锯齿完成 CFD 风格渲染；
7. 使用多尺度 LOD 同时保证全局清晰与局部细节。

---

## 16. 官方参考资料

- ParaView Surface LIC：<https://docs.paraview.org/en/latest/UsersGuide/displayingData.html>
- ParaView Stream Tracer 与 Custom Source：<https://docs.paraview.org/en/latest/UsersGuide/filteringData.html>
- PyVista `streamlines_from_source`：<https://docs.pyvista.org/api/core/_autosummary/pyvista.datasetfilters.streamlines_from_source>
- PyVista `streamlines_evenly_spaced_2D`：<https://docs.pyvista.org/api/core/_autosummary/pyvista.datasetfilters.streamlines_evenly_spaced_2d>
- PyVista 2D Streamlines 示例：<https://docs.pyvista.org/examples/01-filter/streamlines_2d>
- VTK `vtkStreamTracer`：<https://vtk.org/doc/nightly/html/classvtkStreamTracer.html>
- VTK `vtkImageDataLIC2D`：<https://vtk.org/doc/nightly/html/classvtkImageDataLIC2D.html>
- PyVista trame：<https://docs.pyvista.org/user-guide/jupyter/trame.html>
- PyVista HTML 导出：<https://docs.pyvista.org/api/plotting/_autosummary/pyvista.plotter.export_html>
- VisIt Integral Curve System：<https://visit-sphinx-github-user-manual.readthedocs.io/en/v3.5.0/using_visit/Operators/OperatorTypes/ICS/index.html>
