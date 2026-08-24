
当前 v3 已经成功解决 explicit junction 的：

self-intersection
internal faces
internal caps
non-manifold

问题。

不要推翻现有 v3 hybrid reconstruction。

本轮 v4 只处理两个剩余问题：

1. JUNCTION_PORT_REGION_CONFLICT
2. junction area metric / possible junction bulging

==================================================
TASK A — CFD Core ROI / Context Domain Separation
==================================================

当前 representative ROI 必须保持完全不变。

定义：

CORE_ROI = 原 Sampling Project 选出的 representative ROI

CFD_DOMAIN = 用于 CFD geometry / flow 的计算区域

要求：

CORE_ROI ⊂ CFD_DOMAIN

Sampling 的：

feature
cluster
representative identity

全部基于 CORE_ROI，

不得因为 CFD 要求重新计算。

---

对于每个原 CUT_PORT：

计算其到最近 junction 的沿 centerline geodesic distance。

调用现有：

collar
overlap
spatial-separation
port-conflict

检查。

如果原 CUT_PORT 满足全部条件：

保留。

如果不满足：

禁止缩短 collar。

禁止改变 source radius。

禁止移动 source SWC。

而是：

沿对应 global SWC branch 向 CORE_ROI 外继续追踪真实 centerline。

寻找新的 CFD cut position。

新位置必须同时满足：

1. junction collar 可以完整存在；
2. implicit overlap 可以完整存在；
3. port extension 区域不与 junction 冲突；
4. spatial separation PASS；
5. 不位于另一个 bifurcation core。

使用已有约束动态计算 required distance。

不要硬编码：

7.37 μm。

---

如果向外延伸过程中遇到下一个 bifurcation：

将该 bifurcation 及其相关 branch 纳入 CFD context，

然后继续寻找下一组安全 CFD ports。

---

保存：

original_cut_port
new_cfd_port
added_centerline_length
added_branch_ids
reason_for_extension

原 CUT_PORT 不删除，

标记：

CORE_ROI_BOUNDARY

新 port 标记：

CFD_BOUNDARY_PORT

---

CFD 求解未来作用于：

CFD_DOMAIN

但 MB / ULM statistics 默认只评价：

CORE_ROI

==================================================
TASK B — 修正 Junction Area QC
===============================

当前：

A/source A = 11.856

不能直接作为 geometry failure。

原因：

junction 附近一个截面可能同时穿过多个 branch lumen。

因此重新实现：

BRANCH_LOCAL_CROSS_SECTION_QC

---

对于最终 surface 的每个 point / triangle，

计算其距离每一个 source branch centerline 的：

radius-normalized distance：

d_norm_i = d_i / r_i

定义：

owner_branch =
argmin(d_norm_i)

同时计算：

ownership_margin =
second_min(d_norm) - min(d_norm)

---

ownership_margin 太小时：

标记：

AMBIGUOUS_JUNCTION_SURFACE

这些 surface 不进入：

branch radius fidelity。

---

对 branch i 的截面：

只使用：

owner_branch == i

且：

ownership confidence PASS

的 surface，

计算：

A_branch_local(s)

以及：

r_eq_branch_local(s)
====================

sqrt(A_branch_local / π)

---

输出：

old_total_plane_area
new_branch_local_area
source_area
ownership_confidence

重新评价 Junction 49。

==================================================
TASK C — Junction Core 不使用 branch-radius fidelity 硬约束
============================================================

明确区分：

BRANCH_FIDELITY_REGION

和：

JUNCTION_CORE_REGION。

---

BRANCH_FIDELITY_REGION：

具有 source SWC radius ground truth。

继续严格评价：

r_CFD vs r_SWC。

---

JUNCTION_CORE_REGION：

SWC 只提供：

topology
centerline
junction node radius

但不提供真实 3D bifurcation wall。

因此：

不要要求：

A_junction / πr² ≈ 1

作为硬 PASS 条件。

这里只评价：

watertight
boundary edges
nonmanifold
self-intersection
internal faces
internal caps
degenerate triangles
surface smoothness
no extreme artificial throat

area ratio 作为 descriptive metric 保存，
而不是单独作为 FAIL 条件。

==================================================
TASK D — 只有新 QC 仍确认真实 bulge 才改 geometry
==================================================

如果新的：

branch-local area QC

仍然显示 junction 外的 branch：

A_local / A_source

严重异常，

才启用：

CONTROLLED_LOCAL_IMPLICIT

不要立即默认启用。

---

CONTROLLED_LOCAL_IMPLICIT：

不再让每条 branch polyball 在 junction 后方无限产生 spherical support。

定义：

junction core field：

phi_J(x)
========

||x - x_J|| - r_J

对于 incident branch i：

定义 outward tangent：

t_i

branch implicit field 只在：

(x - x_J) · t_i >= -epsilon

的区域有效。

其余位置：

phi_i = +inf

最终：

phi =
min(
phi_J,
phi_1_clipped,
phi_2_clipped,
...
)

即：

single controlled junction core
+
direction-clipped branch fields。

---

不要使用：

sphere mesh + Boolean tubes。

所有内容仍在同一个 local implicit field 中完成。

因此继续要求：

self intersections = 0
internal faces = 0
internal caps = 0

==================================================
TASK E — v4 Acceptance
=======================

Junction 13：

原：

distance to CUT_PORT ≈ 2.162 μm

不允许通过缩短 collar 解决。

要求：

自动生成 expanded CFD domain，

所有 final CFD ports：

JUNCTION_PORT_REGION_CONFLICT = 0。

---

Junction 49：

先输出：

old total-plane area ratio
new branch-local area ratio

明确判断：

原 11.856 是否主要来自 multi-branch cross-section contamination。

---

如果 branch-local geometry 实际正常：

保留当前 v3 hybrid geometry。

不要为了降低 junction total area ratio 修改 geometry。

---

如果 branch-local geometry 仍异常：

才运行 controlled local implicit，

比较：

v3
vs
v4 controlled implicit。

==================================================
必须保留 v3 已达到的指标
========================

self intersections = 0

internal faces = 0

internal caps = 0

boundary edges = 0

nonmanifold edges = 0

surface components = 1

collar radius error 不明显恶化

radius P95 不明显恶化。

==================================================
最终报告必须回答
================

1. J13 是否通过 CFD context extension 解决？
2. 新 CFD boundary 比原 CORE boundary 向外增加多少 μm？
3. 是否增加了真实 global SWC branches？
4. CORE ROI 是否保持完全不变？
5. J49 原来的 11.856 area ratio 中，
   有多少是 multi-branch section contamination？
6. branch-local radius / area error 是多少？
7. J49 是否真的存在需要进一步修改 geometry 的 artificial bulge？
8. 如果不需要：
   明确说明为什么 v3 hybrid 可以接受。
9. 如果需要：
   controlled local implicit 是否改善？
10. 当前完整 representative ROI 是否达到 CFD-ready PASS？
