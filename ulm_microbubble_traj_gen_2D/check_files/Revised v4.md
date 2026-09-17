> From: https://chatgpt.com/c/6a43a3c2-4534-83ee-9625-4c3b2ed4b811

# chatgpt response

下面给出一版完整、固定的 lumen mask 改进方案。它的目标不是单纯让 mask 看起来更粗、更连通，而是在不破坏 DCCO 原始血管半径和拓扑结构的前提下，降低 `narrow_cells`，保持 `holes=0`，并让生成的二维 lumen mask 能够可靠服务于后续二维流场求解和微泡轨迹生成。

## 总体原则

当前 lumen mask 的生成应坚持一个基本原则：**mask 是 DCCO 血管图在二维流场网格上的受控离散化，而不是后期图像修补出来的形状**。DCCO/CCO 本身允许在二维或三维 vascular territory 中生成血管网络，且每段血管具有半径、近端点、远端点等明确几何属性。  因此，mask 生成的首要任务不是“把线画粗”，而是忠实地把这些有半径的二维血管段转化为可解析的 lumen 区域。

由于后续二维流场求解依赖不可压条件、边界通量和壁面无穿透，mask 中的细颈、孔洞、锯齿边界、错误端帽和分叉处裂缝都会直接影响速度场。DCCO 原始血流模型采用低 Reynolds 数、黏性主导的 Poiseuille 型血流假设，血管半径对阻力和流量具有核心意义。 因此，不能通过随意扩大半径、全局 dilation 或大尺度 closing 来降低 `narrow_cells`。这些操作虽然会让 mask 更粗，但会改变血管半径的物理含义。

改进后的固定思路是：**先保证参与二维 PDE 求解的血管在当前网格 spacing 下可解析，再在连续几何层面构造无孔、无细颈的 lumen polygon，最后通过受控 rasterization 转成 mask。**

---

## 1. 先做可解析性筛查，而不是事后修补

在生成 mask 之前，必须先检查每一段血管的物理直径在当前网格上是否足够宽。设网格间距为 $\Delta x$，第 $i$ 段血管物理半径为 $R_i$，则它的理论像素直径为：

$$
D_{i,\mathrm{px}}=\frac{2R_i}{\Delta x}.
$$

正式进入二维流场求解的血管段应满足：

$$
D_{i,\mathrm{px}}\geq 10.
$$

这里用 10 而不是 8，是因为 rasterization 后边界离散、斜向血管和分叉拼接可能损失 1–2 个有效像素。后续 narrow 判定仍可用 8 个像素作为底线，但生成阶段应留出余量。

因此，`min_mask_radius` 不应作为正式流场 mask 中的主要补救手段。它可以作为几何连通性诊断或可视化保护，但不应让一根物理直径只有 3–4 个像素的血管在 mask 中被强行扩成 8–10 个像素，然后仍然在物理参数中使用原始小半径。这样会造成几何域和血流参数不一致。

固定规则应是：**若某个 ROI 内存在参与流场求解但 $D_{i,\mathrm{px}}<10$ 的关键血管段，则该 ROI 不合格，应提高空间分辨率或重新选择 ROI；不通过扩大该段 mask 半径来伪装修复。**

---

## 2. 连续几何层面构造 lumen，而不是依赖像素形态学

每段血管仍应先在连续几何空间中生成 corridor polygon。也就是说，取该血管在 X-Z 平面上的中心线段，以其物理半径 $R_i$ 生成带宽通道。这个通道的宽度应来自 DCCO 原始半径，而不是后期 mask 半径。

内部连接点和外部边界点要分开处理。对于内部 junction，即父血管和子血管相接的位置，端部可以使用圆端帽，因为这里的目标是保证分叉核心连续、无裂缝。对于入口和出口端点，应使用平端帽，因为它们在流场求解中是开放截面，需要形成清晰的入口或出口边界，而不是被圆帽包住。

因此，端帽规则固定为：

```text
内部 junction 端点：round cap
入口 / 出口端点：flat cap
```

这一步很重要。入口和出口如果继续使用圆端帽，会让边界变成封闭圆头，后续再人为从圆头中找一条出口线，容易制造边界通量误差和局部 divergence。平端帽能让入口/出口边界更明确。

---

## 3. 分叉区域采用重叠式连续连接

当前已有 junction disk 和 junction convex hull，这是合理的，但要进一步明确其功能。分叉区域的目标不是简单“填满一点”，而是避免父血管与子血管在连续几何中只是刚好相接，导致 raster 后出现单像素细颈。

对每个 junction，设连接到该节点的最大血管半径为 $R_{\max}$。定义分叉核心半径为：

$$
R_j=R_{\max}+2\Delta x.
$$

在每条连接血管方向上，从 junction 向外保留一段重叠连接长度：

$$
L_j=\max(2R_j,6\Delta x).
$$

这意味着每条进入 junction 的血管都不是“刚好碰到节点”，而是伸入一个共同的局部腔体。最终分叉区域由三部分构成：

```text
junction disk
+ incident vessel overlap corridors
+ smooth junction hull
```

其中 smooth junction hull 可以继续使用凸包思想，但不应过度外扩。它只用于平滑分叉核心，避免小孔、尖角和狭颈。这样既能保持 `holes=0`，又能减少 junction 处的 `narrow_cells`。

---

## 4. 几何并集后直接消除内部孔洞

所有 corridor、junction disk 和 junction hull 合并后，应形成一个整体 lumen polygon。由于 DCCO 生成的是树状血管结构，原文也强调血管树由无交叉、无物理互穿的段组成。 在这种设定下，lumen polygon 内部出现完全封闭的小孔通常不是生理结构，而是几何 union 或 rasterization 的伪影。

因此，在连续几何层面应进行两步修正：

第一，检查 polygon 有效性，修复自交、裂缝或无效环。

第二，移除 lumen polygon 内部的 interior holes。也就是说，连续几何最终输出前应保证：

```text
polygon.is_valid = true
polygon interior holes = 0
```

这样，后续 raster mask 中的 `holes=0` 不再主要依赖二值图像填孔，而是在连续几何源头上得到保证。

---

## 5. Rasterization 使用超采样面积判据，而不是直接依赖 all_touched

当前启用 `all_touched` 有助于保持细血管连通，但它会让斜向边界变厚，并可能在某些局部产生不均匀边界。正式方案应改为超采样 rasterization。

固定做法是：先在目标网格的 8 倍分辨率上 rasterize 连续 lumen polygon，然后对每个目标 cell 统计其高分辨率子像素中有多少比例属于 lumen。若该比例大于或等于 0.5，则将目标 cell 标记为 lumen：

$$
\mathrm{occupancy}(c)\geq 0.5 \Rightarrow c\in \mathrm{lumen}.
$$

这种方法比 `all_touched=True` 更稳定。它既不会像 cell-center 判据那样容易断裂，也不会像 all_touched 那样系统性加粗边界。由于 junction 已经在连续几何层面被加强，连通性不应再主要依赖 all_touched。

因此，正式 rasterization 固定为：

```text
8× supersampling
occupancy threshold = 0.5
no global all_touched at target resolution
```

如果 raster 后仍然出现孔洞或细颈，应回到连续几何层修 junction 或提高分辨率，而不是在像素层做大尺度 dilation/closing。

---

## 6. Raster 后只允许保守清理

完成 rasterization 后，只进行最小必要的二值清理。允许的操作包括：

```text
保留最大连通分量
填补内部孔洞
删除极小孤立碎片
```

不允许把全局 dilation、全局 closing 或大尺度 erosion 作为正式修复手段。原因很简单：这些操作会移动血管壁，改变局部直径，而血管半径在 DCCO hemodynamics 中具有物理意义。

填孔可以保留，因为树状 lumen 内部的小孔通常是几何或 rasterization 伪影。清理后必须重新计算：

```text
components
holes
narrow_cells
narrow_segments
```

只有 `components=1` 且 `holes=0`，mask 才能进入流场求解。

---

## 7. Segment 属性分配要避免误判 junction

你当前做法是把每个 lumen cell 归属给最近血管段，然后 junction 新增 cells 从最近已归属 cell 继承属性。这个做法能保证每个 cell 都有物理属性，但对 `narrow_cells` 诊断可能造成误判。尤其在分叉核心区域，一个较宽的 junction 腔体可能被归属给某条细子支，从而继承了很小的 graph diameter，被错误判为 narrow。

改进后仍然保留 `segment_id`，不新增复杂实体，但需要改变 junction cells 的诊断方式。对于普通血管段区域，局部有效直径由原始 DCCO 半径和 mask 实际宽度共同决定。对于 junction 区域，不使用单一子支半径判断 narrow，而使用 SDF/medial-axis 得到的实际局部宽度判断。

换言之，junction cells 可以继续继承最近血管段的速度、黏度和流量属性，用于速度初始化；但在 `narrow_cells` 诊断中，junction 区域应依据实际 mask 宽度，而不是某一条 incident segment 的原始半径。

---

## 8. narrow_cells 的最终定义

`narrow_cells` 不应表示“靠近壁面”的 cell。任何血管都有壁面，靠近壁面的 cell 很多；如果用 SDF 小于某阈值来定义 narrow，会把正常边界层也算进去，指标就失去意义。

正式定义应为：

> **narrow cell 是 lumen 中属于局部有效血管直径小于 8 个网格点的 cell。**

对普通血管段，先计算 graph 直径像素数：

$$
D_{\mathrm{graph,px}}=\frac{2R_i}{\Delta x}.
$$

同时，根据 mask 的 medial-axis 或 SDF 估计实际局部直径：

$$
D_{\mathrm{mask,px}}=2\,\mathrm{EDT}_{\mathrm{centerline}}.
$$

局部有效直径取两者较小值：

$$
D_{\mathrm{eff,px}}=\min(D_{\mathrm{graph,px}},D_{\mathrm{mask,px}}).
$$

当：

$$
D_{\mathrm{eff,px}}<8
$$

该区域计为 narrow。

这样可以同时捕捉两种问题：原始 DCCO 半径在当前网格下太小，以及 rasterization 后出现了非物理窄颈。

最终输出不应只有 `narrow_cells`，而应包括：

```text
narrow_cells
narrow_fraction
narrow_segments.csv
min_resolved_diameter_px
p5_resolved_diameter_px
median_resolved_diameter_px
narrow_in_inlet_or_outlet
narrow_in_junction_core
```

其中 `narrow_in_inlet_or_outlet` 和 `narrow_in_junction_core` 比全局 narrow fraction 更重要。少量末端小分支 narrow 可以标记为低置信区域；入口、出口和主分叉 narrow 则必须修复。

---

## 9. 固定质量门控

生成 mask 后必须通过以下门控，才能进入二维流场求解：

```text
components = 1
holes = 0
global narrow_fraction < 1%
narrow_in_inlet_or_outlet = 0
narrow_in_main_junction_core = 0
p5_resolved_diameter_px >= 8
```

其中 `global narrow_fraction < 1%` 是硬性合格线；`<0.5%` 可以作为优化目标。你当前最新结果约为 0.70%，如果 narrow 主要出现在末端小分支而不是入口、出口或主分叉，则可以认为 mask 已经可用于纯流动轨迹生成；如果 narrow 集中在关键流动路径，则仍需修复。

这里的逻辑是：**全局 narrow_cells 不能完全为 0 不是问题；关键区域不能 narrow 才是重点。**

---

## 10. 自动修正闭环

完整流程不应是“生成一次 mask 然后直接求解”，而应是一个质量闭环：

```text
生成连续 lumen polygon
→ 超采样 rasterization
→ 保守二值清理
→ 计算 SDF 和局部有效直径
→ 输出 holes / narrow diagnostics
→ 若 holes>0，回到连续几何层修复 interior rings 和 junction
→ 若 junction narrow，增加 junction overlap 或提高分辨率
→ 若 terminal / inlet / outlet narrow，修正端帽和边界截面
→ 若真实血管物理直径不足，提升分辨率或缩小 ROI
→ 重新 rasterize
```

注意，这个闭环的修正优先级是固定的：

1. **几何层修复优先于像素层修复。**
2. **提高分辨率优先于强行扩大血管半径。**
3. **修 junction overlap 优先于全局 dilation。**
4. **入口/出口使用 flat cap，不能用圆帽掩盖边界问题。**

---

## 11. 最终完整方案表述

改进后的 lumen mask 生成流程应表述为：

本方法首先从 DCCO 血管图中读取每段血管的中心线、物理半径、父子连接关系、入口和出口信息。由于后续二维流场求解需要在网格上解析血管腔结构，所有参与求解的血管段先进行分辨率筛查，要求其物理直径在当前网格上至少覆盖 10 个像素。未达到该标准的区域不通过扩大 mask 半径来伪装修复，而应通过提高空间分辨率或重新选择可解析 ROI 处理。

随后，在连续几何空间中将每段血管中心线缓冲为带半径的 corridor polygon。内部 junction 使用圆端帽、junction disk、重叠 corridor 和平滑 hull 共同构造，以保证分叉核心连续、无细颈、无小孔。入口和出口端点采用 flat cap，以形成明确的开放边界截面。所有 corridor 和 junction 几何合并后，进行几何有效性检查，并移除 interior holes，使连续 lumen polygon 本身满足无孔洞要求。

连续几何完成后，采用 8 倍超采样 rasterization，并以面积占比 0.5 作为 lumen 判据生成二维 mask。raster 后仅进行保守清理，包括保留最大连通分量、填补内部小孔和删除极小碎片，不使用全局 dilation 或 closing 改变血管半径。随后根据 SDF、medial axis 和 DCCO 原始半径共同计算局部有效血管直径，并以 8 个像素作为 narrow 判据生成 `narrow_cells` 和 `narrow_segments.csv`。

最终 mask 必须满足单连通、无孔洞、关键入口/出口和主分叉无 narrow、全局 narrow fraction 小于 1% 的质量门控，才能进入二维流场求解。这样生成的 lumen mask 既保留了 DCCO 血管半径和分叉拓扑的生物物理基础，又避免了栅格化带来的孔洞、细颈和局部不可解析区域。

这套方案的核心思想是：

> **holes=0 主要靠连续几何层保证；narrow_cells 降低主要靠解析度、junction 重叠和超采样 rasterization 保证；物理半径不通过后期形态学操作随意改变。**