
当前 v5 已经证明：

self-intersection = 0
internal faces = 0
internal caps = 0
boundary edges = 0
nonmanifold edges = 0
degenerate triangles = 0
surface components = 1

radius P95 ≈ 0.327%
collar radius error ≈ 0.37%

因此：

v5 已经满足 topology validity 和 radius fidelity。

但是人工检查最终 STL 后仍发现：

1. hybrid junction 的 implicit→explicit collar 处存在明显锯齿状 seam；
2. CORE branch→CFD extension 附近仍存在环形 shoulder / step。

这些缺陷虽然没有造成 open crack 或 non-manifold，
但说明 surface fairness / first-order continuity 仍然不足。

本轮定义为：

v6 Continuous-Field Transition Refinement

============================================================
一、本轮最重要的原则
====================

不要再修改：

source SWC
source radius
CORE ROI
CFD context extension topology
junction topology
CUT_PORT mapping

不要重新启用：

legacy junction sphere

不要用：

更多 transition rings

作为主要解决方案。

不要通过：

aggressive smoothing

掩盖问题。

正式目标：

【从算法上移除产生 visible seam 的离散 surface stitching。】

============================================================
二、首先修改验收逻辑
====================

当前 P95/P99 是较大区域统计，
容易掩盖只发生在 interface ring 上的局部缺陷。

新增：

INTERFACE-SPECIFIC QC

分别明确记录：

A.
CORE→PORT_EXTENSION_INTERFACE

B.
IMPLICIT→EXPLICIT_JUNCTION_INTERFACE

对这些 interface edge 单独计算：

dihedral angle mean
dihedral angle P95
dihedral angle P99
dihedral angle max

不要再只使用 whole-port / whole-junction P95/P99。

同时计算：

interface_edge_count

interface_length

============================================================
三、Port：诊断 C1 continuity
============================

目前已经确认：

position continuity PASS
radius continuity PASS
area continuity PASS

因此继续检查：

CENTERLINE TANGENT CONTINUITY

在原 CORE CUT_PORT 两侧：

source side:
拟合最后若干真实 centerline points

extension side:
拟合最开始若干 CFD extension points

不要只使用最后一条 edge。

推荐使用：

weighted PCA / local linear regression

估计：

t_source

t_extension

计算：

theta_tangent =
acos(t_source · t_extension)

保存：

port_tangent_jump_deg

============================================================
四、Port radius derivative continuity
=====================================

用 CUT_PORT 前：

最后 3～5 个 resampled source radius

估计：

(dr/ds)_source

extension 初始为 constant radius 时：

(dr/ds)_extension = 0

计算：

radius_slope_jump =
|(dr/ds)_source - (dr/ds)_extension|

如果明显非零，

则 port shoulder 很可能是：

radius C1 discontinuity。

============================================================
五、Port v6 几何修改
====================

不要修改 CORE ROI 内 source geometry。

只允许修改：

CFD_EXTENSION

区域。

---

5.1 centerline

extension direction 不再直接使用单个 endpoint tangent。

使用：

source branch 最后若干点

进行局部 line / polynomial fit，

得到稳定：

t_fit

使用该 t_fit 作为 extension 主方向。

---

如果：

theta_tangent_jump

仍然较大，

建立一个短：

TANGENT CONDITIONING ZONE

只位于 CFD-derived extension 内。

长度：

L_tangent_blend

建议初始：

1D～2D

配置化。

使用 cubic Hermite curve：

起点：
CUT_PORT

起始 tangent：
source fitted tangent

终点：
straight extension axis

终点 tangent：
straight extension tangent

要求：

C1 continuity。

============================================================
六、Port radius conditioning
============================

如果 source radius slope：

m0 = dr/ds

明显不为零，

不要在 CUT_PORT 后瞬间变成：

constant r。

在 extension 内使用：

RADIUS CONDITIONING ZONE。

要求：

r(0) = r_cut

r'(0) = m0

r'(L_blend) = 0

之后：

constant radius。

可以使用 cubic Hermite radius profile。

不要修改 CUT_PORT 以前的真实 SWC radius。

保存：

radius_profile_before
radius_profile_after

============================================================
七、Port 验收
=============

重新输出：

port_centerline_tangent_profile.png
port_radius_profile.png
port_interface_dihedral.png

重点比较：

v5 interface max dihedral
v6 interface max dihedral

以及：

area error
radius fidelity

必须保持。

============================================================
八、Junction：停止 boundary-loop stitching
==========================================

当前：

implicit junction surface
→ extract loop
→ intermediate rings
→ explicit branch surface

会在 surface level 建立 zipper seam。

v6 正式默认：

NO SURFACE LOOP STITCHING

保留旧实现：

transition_backend = loop_stitch_legacy

仅用于 regression。

正式：

transition_backend = continuous_implicit_field

============================================================
九、新方法：一个连续 local scalar field
=======================================

对于每一个 incident branch，

已有：

junction core polyball field：

phi_J(x)

同时构建一个：

analytic / polyball branch tube field：

phi_B(x)

该 field 必须对应当前：

centerline + radius

的 branch tube。

============================================================
十、定义 Transition Coordinate
==============================

沿当前 incident branch centerline 定义：

s

其中：

s = 0

位于 junction core 一侧，

s = L_transition

位于 explicit branch 一侧。

构建：

w(s)

满足：

w(0) = 0
w(L) = 1

w'(0) = 0
w'(L) = 0

使用 quintic smoothstep：

w(t) =
6t^5 - 15t^4 + 10t^3

t = s/L

============================================================
十一、Continuous field blend
============================

在 transition collar 中：

phi_transition(x)
=================

[1-w(s(x))] phi_J(x)
+
w(s(x)) phi_B(x)

注意：

s(x)

应通过最近 branch centerline projection 得到，

不是简单使用 Cartesian coordinate。

============================================================
十二、Transition 最外端必须变成纯 branch field
==============================================

在：

s >= L_transition

必须：

phi_transition = phi_B

因此 transition 最外端的 zero-set：

必须与 branch tube 的 zero-set 一致。

即：

不再存在：

implicit loop
vs
explicit loop

两个独立 boundary loops。

============================================================
十三、Junction core
===================

在 core：

w = 0

保持当前 v5 已经验证正确的：

local implicit junction。

不要重新调整 core topology。

============================================================
十四、生成方式
==============

为每个：

junction + transition collars

建立一个 local grid。

一次性计算：

continuous scalar field

然后：

ONE marching-cubes extraction

得到：

junction core
+
all transition collars

组成的单一 surface patch。

不要：

每条 collar 分别生成 surface 后再拼。

============================================================
十五、Outer overlap
===================

local implicit patch 在：

w = 1

以后继续向 explicit branch 方向延长：

0.5D～1D

作为：

PURE BRANCH FIELD OVERLAP

此时：

phi_local = phi_B

所以 local surface 与 explicit tube
理论上描述的是同一个管道。

============================================================
十六、最终 Merge
================

优先方案：

在：

PURE BRANCH FIELD REGION

进行裁切。

因为：

local implicit surface
和
explicit tube

在此处来自同一：

centerline + radius field。

---

不要再次在 transition 区做 stitching。

如果仍需要 surface merge：

只能发生在：

w = 1 的纯 branch region。

============================================================
十七、优先尝试同截面 clip + weld
================================

在 pure branch region 选：

merge plane。

分别裁切：

local patch
explicit tube

得到：

两个近乎一致的 circular loops。

---

计算：

loop center error
loop radius error
loop Hausdorff distance

如果均低于 tolerance：

将两个 loop 投影至一个共同 target loop，

然后：

vertex weld

而不是：

增加 transition triangles。

这样 seam 位于：

纯圆管区域，

不会出现在 junction transition。

============================================================
十八、如果 loop 不足够一致
==========================

则使用：

small local manifold Boolean

但只发生在：

pure branch field overlap。

禁止 Boolean 进入：

junction core
transition collar。

============================================================
十九、重要：提高 local implicit resolution 只能用于 convergence
===============================================================

不要默认认为：

提高 marching cubes resolution

就能解决 seam。

因为当前主要问题是：

surface stitching architecture。

但是 continuous-field 方法完成后，

对：

cells_across_min_diameter

测试：

16
20
24

观察：

normal continuity
triangle count
runtime

============================================================
二十、Local mesh quality
========================

marching cubes 后只对：

continuous implicit patch

进行：

degenerate cleanup
small sliver cleanup

不要全局 smoothing。

============================================================
二十一、新增 Transition-specific QC
===================================

对每一条 junction collar，

将 surface 分成：

CORE
TRANSITION
PURE_BRANCH
EXPLICIT

分别计算：

normal jump P95
P99
max

mean roughness

triangle aspect ratio P95
max

============================================================
二十二、Silhouette QC
=====================

当前图片中的锯齿已经改变 silhouette。

新增：

多视角 silhouette smoothness test。

至少从：

3 个正交 view
+
3 个斜视角

提取 surface silhouette。

计算：

silhouette curvature variation

以及：

large corner count。

不要只依赖 shading image。

============================================================
二十三、v5 seam provenance 可视化
=================================

生成：

v5_interface_edges.png

明确将：

loop-stitch interface edges

用高亮线显示。

确认人工看到的锯齿
是否确实与：

stitch interface

位置重合。

如果不重合，

必须重新诊断，
不能盲目修改 transition backend。

============================================================
二十四、v6 同相机 comparison
============================

对当前截图位置输出：

v5
v6

完全相同：

camera
lighting
flat shading
smooth shading

四套图。

必须包含：

wireframe overlay。

============================================================
二十五、Acceptance 规则必须升级
===============================

不要再仅凭：

watertight
P95
P99

宣布 smooth PASS。

新的 smooth PASS 必须满足：

topological QC 全部 0

并且：

interface-specific seam QC PASS。

============================================================
二十六、暂定不是硬编码的比较指标
================================

v6 应至少做到：

junction interface P99
<
v5 junction interface P99

junction interface max
显著降低

port interface tangent jump
显著降低

port interface normal max
显著降低

silhouette roughness
明显下降

同时：

radius P95
不得明显恶化

collar radius error
不得明显恶化。

============================================================
二十七、严禁事项
================

不要：

继续增加 loop-stitch intermediate rings

用全局 smoothing 掩盖 seam

修改 source SWC

修改 source radius

重新抽 ROI

重新使用 Mask 参与正式 surface reconstruction

恢复 junction sphere

只看 smooth shading 宣布问题解决

============================================================
二十八、必须保留 v5 所有优点
============================

v6 必须继续：

self intersection = 0

internal faces = 0

internal caps = 0

boundary edges = 0

nonmanifold edges = 0

degenerate triangles = 0

components = 1

============================================================
二十九、Synthetic Tests
=======================

必须加入：

straight branch
→ straight extension
带 source taper

用于检查 radius slope discontinuity。

curved branch
→ CFD extension

用于检查 tangent discontinuity。

Y junction
continuous implicit transition

acute Y junction

unequal-radius Y junction

curved daughter branches

============================================================
三十、输出
==========

保存：

outputs/model_generate/<run_id>/v6/

至少：

diagnostics/
    port_tangent_jump.csv
    port_radius_slope.csv
    junction_interface_dihedral.csv
    silhouette_qc.csv

figures/
    v5_interface_edges.png
    port_tangent_profile.png
    port_radius_profile.png
    junction_field_transition.png
    junction_interface_dihedral.png
    v5_v6_flat_shading.png
    v5_v6_smooth_shading.png
    v5_v6_wireframe.png
    silhouette_comparison.png

============================================================
三十一、最终必须汇报
====================

PORT：

当前台阶对应位置的：

centerline tangent jump
是多少？

radius slope jump
是多少？

v6 修改后分别是多少？

interface max dihedral：
v5 vs v6

area/radius fidelity 是否保持？

---

JUNCTION：

当前图中的锯齿是否和 loop-stitch interface 重合？

如果重合：

确认：
LOOP_STITCH_ZIPPER_ARTIFACT

continuous implicit-field transition 后：

interface P95
P99
max

分别是多少？

silhouette roughness 改善多少？

---

TOPOLOGY：

所有 defect count 是否仍然为 0？

---

PERFORMANCE：

runtime
triangles

v5 vs v6。

---

最终回答：

当前模型是否同时满足：

1. topology validity
2. radius fidelity
3. port C1 continuity
4. junction surface continuity
5. CFD-ready meshing suitability

不要仅因为测试通过就输出 PASS。

必须结合：

interface-local metrics
+
silhouette
+
同相机实际可视化

共同判断。
