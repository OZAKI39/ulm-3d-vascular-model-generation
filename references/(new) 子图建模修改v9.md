
当前 v8 已确认：

- saw-tooth 最早出现在 S0 raw FlyingEdges；
- remeshing 不是根因；
- cross-branch PolyBall ownership switching 与缺陷高度相关；
- 但 defect→branch-switch coverage 只有 26.808%；
- branch-level smooth union 无法完全消除 silhouette saw-tooth。

请基于当前版本代码继续执行 v9。

本轮重点不是继续增加 k/r，
而是检查并解决：

【same-branch segment-to-segment PolyBall switching】

因为当前代码：

1. geometry.centerline_smoothing 默认 false；
2. branch 坐标通过 np.interp 形成 piecewise-linear centerline；
3. unified PolyBall field 对每一个 straight segment 分别求值后取 hard minimum；
4. v8 smooth union 只发生在 branch fields 之间，
   不会消除同一 branch 内 segment switching。

==================================================
PHASE 1 — Segment-level ownership diagnosis
============================================

先不要修改正式 geometry。

扩展 PolyBall evaluator，使其能够返回：

winner_segment_id
winner_branch_id
segment_parametric_t

注意：

不要破坏现有 evaluate API；
新增 diagnostic method 或 optional return_owner。

---

对 S0 raw FlyingEdges surface 每个 vertex：

计算：

winner_segment_id
winner_branch_id
second_segment_id
second_branch_id
ownership_margin

---

将 surface adjacency switching 分为：

TYPE_A:
CROSS_BRANCH_SWITCH

winner_branch 不同

TYPE_B:
SAME_BRANCH_ADJACENT_SEGMENT_SWITCH

branch 相同，
segment 在 source polyline 中相邻

TYPE_C:
SAME_BRANCH_NONADJACENT_SEGMENT_SWITCH

branch 相同但 segment 不相邻

---

分别计算：

switch edge count

gradient-jump:
mean
P95
P99
max

surface dihedral:
mean
P95
P99
max

---

与现有 saw-tooth defect overlay。

输出：

defect_near_cross_branch_fraction

defect_near_same_branch_adjacent_segment_fraction

defect_near_same_branch_nonadjacent_segment_fraction

---

如果 same-branch segment switching
解释了主要剩余 defect：

结论：

ROOT_CAUSE_CONFIRMED =
PIECEWISE_LINEAR_CENTERLINE_SEGMENT_CREASE

==================================================
PHASE 2 — Source centerline tangent audit
==========================================

对每条 resampled branch 的相邻 segments：

t_i =
normalized(p_{i+1}-p_i)

计算：

theta_i =
acos(t_i · t_{i+1})

输出：

branch_id
source node neighborhood
segment ids
theta_deg
local radius
distance_to_junction

---

将：

theta_deg

与：

segment-switch gradient jump
surface dihedral
saw-tooth defect

做相关性分析。

==================================================
PHASE 3 — 不使用当前 moving-average smoothing
==============================================

不要简单打开：

geometry.centerline_smoothing = true

现有 _light_smooth 只允许保留作为 legacy control。

正式新增：

CFD_DERIVED_SPLINE_CENTERLINE

不得修改：

raw_points_um
raw_radius_um
source SWC
source topology

==================================================
PHASE 4 — C1 spline centerline
===============================

每条 branch 构建：

continuous parametric curve c(s)

推荐优先实现：

centripetal Catmull-Rom
或等价的 cubic Hermite interpolation。

要求：

1. 所有原始 SWC branch points 保留为插值约束；
2. branch 两端点精确不变；
3. bifurcation node 精确不变；
4. CUT_PORT / CFD port 相关 endpoint 精确不变；
5. 不改变 topology；
6. 不允许 curve self-intersection。

---

如果使用 SciPy CubicSpline：

必须增加 overshoot / deviation QC。

不要无条件接受 natural cubic spline。

==================================================
PHASE 5 — Radius field
=======================

radius 继续使用 PCHIP：

r(s)

不得使用可能 overshoot 的 unrestricted cubic radius spline。

==================================================
PHASE 6 — Adaptive spline discretization
=========================================

不要只按固定：

Δs = 0.35 r

采样。

新增三个约束：

max_spacing

Δs <= alpha r

max_tangent_angle

angle(t_i,t_{i+1}) <= theta_max

max_sagitta

distance(
spline midpoint,
straight-segment chord
)
<= eta * r

---

测试：

theta_max:
0.5°
1.0°
2.0°

eta:
0.01
0.02

不要预设 winner。

==================================================
PHASE 7 — Centerline fidelity QC
=================================

由于 spline 是 CFD-derived geometry，
必须证明没有篡改 source anatomy。

计算：

source point → spline distance

spline → source polyline distance

Hausdorff distance

P95 distance

branch length change

junction position error

endpoint position error

---

junction/endpoints：

必须 exactly 0。

==================================================
PHASE 8 — Build v9 hard PolyBall
=================================

使用：

dense C1 spline-derived polyline
+
PCHIP radius

建立 unified PolyBall。

暂时仍使用：

hard min。

先回答：

仅消除 piecewise-linear tangent jumps
是否已经明显减少 saw-tooth。

输出：

V7 original hard-min
vs
V9 smooth-centerline hard-min

同相机：

flat
wireframe
silhouette。

==================================================
PHASE 9 — Recompute segment-level defect
=========================================

重点比较：

same-branch segment-switch gradient P99

surface dihedral P99

defect coverage

silhouette roughness

如果明显下降：

证明中心线折线化是重要根因。

==================================================
PHASE 10 — Junction branch smooth union
========================================

只有完成 smooth-centerline 后，

再对真实 cross-branch switching 做 local smooth union。

不要对整个 junction 球形区域无条件 blend。

新增 competition-aware support。

计算：

phi_sorted_1
phi_sorted_2

competition_margin =
phi_2 - phi_1

只有：

competition_margin
<
competition_threshold * local_radius

时，

smooth union active。

同时必须处于：

true junction neighborhood。

==================================================
PHASE 11 — k sensitivity
=========================

在 smooth-centerline 基础上重新测试：

k/r =
0.10
0.20
0.30
0.40
0.50

因为此前 k sensitivity 是在
piecewise-linear branch fields 上完成的，
不能直接用于新模型。

---

评价：

silhouette defect
gradient jump
radius fidelity
junction volume
hydraulic resistance

==================================================
PHASE 12 — 最终可选 constrained fairing
========================================

如果：

smooth spline
+
competition-aware branch union

后仍存在少量肉眼可见 saw-tooth，

不要继续无限修改 implicit equations。

新增：

junction_local_fairing

只作用于：

detected defect band
+
1–2 ring neighborhood。

---

outer band boundary vertices：

FIXED

branch fidelity region：

FIXED

port：

FIXED

---

优化目标：

minimize local curvature / normal variation

同时限制：

vertex displacement

cross-sectional radius

junction volume

hydraulic resistance。

---

非常重要：

fairing 后：

禁止重新 Newton-project 到 hard-min field。

否则 crease 会回来。

==================================================
PHASE 13 — Fairing safety
==========================

要求：

self-intersection = 0
internal face = 0
internal cap = 0
boundary edge = 0
nonmanifold = 0
component = 1

并比较：

radius P95
branch hydraulic resistance
junction volume

==================================================
PHASE 14 — v9 最终比较
=======================

至少比较：

A.
v7 piecewise-linear hard PolyBall

B.
v8 piecewise-linear smooth branch union

C.
v9 smooth-centerline hard PolyBall

D.
v9 smooth-centerline + smooth branch union

E.
如需要：
v9 + constrained local fairing

---

最终表：

method
cross-branch gradient P99
same-branch segment gradient P99
surface normal P99
silhouette roughness
visible saw-tooth count
radius P95
hydraulic error
volume change
triangle count
runtime
topology status

==================================================
PHASE 15 — Acceptance
======================

不要继续使用：

“silhouette improvement > 5%”

单一规则。

正式 ACCEPT 必须同时：

1. current visible saw-tooth 在同相机 flat/wireframe 中消失；
2. segment-switch defect 明显消失；
3. cross-branch switch 得到平滑；
4. radius P95 < 1%；
5. hydraulic error 在既定容差内；
6. topology QC 全部 PASS；
7. source centerline deviation 在允许范围内；
8. junction/endpoints 精确保持。

==================================================
代码组织
========

不要把所有逻辑继续塞进 unified_polyball.py。

建议新增：

utils/cfd_lumen/
    smooth_centerline.py
    segment_ownership_qc.py
    junction_fairing.py
    v9_pipeline.py
    v9_qc.py
    v9_visualization.py

unified_polyball.py 只扩展必要接口。

==================================================
最终必须回答
============

1. 剩余 defect 中多少由 same-branch segment switching 解释？
2. 最大 centerline tangent kink 在哪里、多少度？
3. C1 spline 后 tangent kink 降到多少？
4. 仅 smooth centerline 是否已经消除大部分 saw-tooth？
5. cross-branch smooth union 还有多少额外收益？
6. 是否最终需要 mesh fairing？
7. 最终视觉锯齿是否真正消失？
8. source geometry 改变量是多少？
9. radius/hydraulic fidelity 是否保持？
10. 最终推荐 backend 是什么？
