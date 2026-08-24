
当前 v4 已经解决 CFD lumen 的拓扑缺陷：

self-intersection = 0
internal faces = 0
internal caps = 0
boundary edges = 0
nonmanifold edges = 0
degenerate triangles = 0
surface components = 1

radius P95 ≈ 0.328%
collar radius error ≈ 0.236%

因此：

【本轮禁止重新修改 source SWC、ROI、radius、global context extension 或 junction topology。】

当前任务改为：

【在保持 v4 已有全部拓扑和半径精度指标的前提下，消除 port 与 hybrid junction 处明显的表面拼接痕迹，提高 CFD lumen 的局部法向连续性和表面光顺性。】

本轮称为：

v5 Surface Continuity Refinement

============================================================
一、当前需要解决的两个问题
==========================

问题 A：

CUT_PORT extension 接口仍存在明显环形视觉接缝。

已有诊断证明：

endpoint position error = 0
endpoint radius error = 0
port area difference < 0.2%
positive overlap exists

因此：

禁止通过修改 radius 解决。

目标是：

从几何构造上消除 separate branch + cylinder seam。

---

问题 B：

local implicit junction 与 explicit branch 的 collar transition
仍存在明显锯齿状 / 突兀表面。

当前：

watertight = true
manifold = true

但：

normal jump P99 ≈ 28.94°

图片中可以观察到明显 surface seam。

因此：

当前模型只是拓扑连续，
还没有达到理想的 surface-normal continuity。

============================================================
二、本轮核心原则
================

必须保持：

source SWC unchanged

CORE_ROI unchanged

CFD_DOMAIN topology unchanged

source radius unchanged

CUT_PORT position logic unchanged

junction local implicit topology unchanged

不得为了光滑：

删除 branch
缩小 radius
修改 SWC centerline
扩大 junction sphere
重新使用 legacy sphere Boolean

本轮只修改：

surface construction
surface transition
surface tessellation
surface normals

============================================================
三、第一项修改：Port 不再使用独立 cylinder Boolean
==================================================

当前逻辑大致为：

explicit branch tube
        +
separate extension cylinder
        ↓
Boolean union
        ↓
visible circular seam

正式 v5 改为：

【centerline extension before tube generation】

============================================================
四、Port centerline extension
=============================

对于最终 CFD boundary port：

已知：

x_cut
r_cut
n_out
L_ext

不要再创建：

Cylinder(radius=r_cut)

作为独立 solid。

而是在对应 branch centerline 后追加真实几何点：

p_0 = x_cut

p_1 =
x_cut + Δs n_out

p_2 =
x_cut + 2Δs n_out

...

p_N =
x_cut + L_ext n_out

这些点属于：

CFD_DERIVED_EXTENSION

不是 source SWC。

必须通过 metadata 明确区分。

---

对应 radius：

r(p_k) = r_cut

保持 constant。

也就是说：

source branch
→ CUT_PORT
→ straight constant-radius extension

成为：

【一条连续 centerline】

然后整个对象一次性通过：

vtkTubeFilter

生成。

============================================================
五、Port tube generation
========================

原：

branch tube
+
extension cylinder

改为：

extended branch centerline
        ↓
one vtkTubeFilter
        ↓
one continuous tube

只在：

最终 CFD port

位置 cap。

不要在原 CORE ROI CUT_PORT 位置 cap。

============================================================
六、Port extension 必须保持 source traceability
===============================================

新增 centerline points 保存：

point_type =
CFD_EXTENSION

source_cut_port_id

source_global_edge_id

distance_from_core_boundary

original_core_cut_position

最终 boundary：

boundary_type =
CFD_BOUNDARY_PORT

原 CUT_PORT 保留：

boundary_type =
CORE_ROI_BOUNDARY

但它不再是 CFD surface boundary。

============================================================
七、Port continuity QC
======================

虽然 seam 理论上应该消失，
仍然必须定量检查。

沿 branch：

CORE cut - 1D
CORE cut
CORE cut + 1D
CORE cut + 2D

连续计算：

A(s)

r_eq(s)

surface normal variation

要求：

radius continuity 不比 v4 恶化。

同时新增：

port_normal_jump_p95
port_normal_jump_p99

============================================================
八、Port 可视化
===============

输出：

port_v4_vs_v5.png

左：

v4 separate cylinder Boolean

右：

v5 continuous-centerline tube

必须用：

相同 camera
相同 shading
相同 tube sides

进行对照。

============================================================
九、第二项修改：Junction seam 不能继续仅靠 Boolean overlap
==========================================================

当前：

local implicit junction patch
+
explicit branch tube
+
Manifold Boolean

虽然 topology 已经正确，
但 collar 位置仍有明显 saw-tooth seam。

不要重新改 implicit core。

问题主要位于：

【implicit → explicit transition band】

============================================================
十、重新定义 Junction 三个区域
==============================

对于每个 incident branch：

定义：

REGION 1:
JUNCTION_CORE

完全 local implicit。

---

REGION 2:
TRANSITION_COLLAR

负责：

implicit surface
→
explicit branch surface

平滑过渡。

---

REGION 3:
EXPLICIT_BRANCH

保持现有高精度 vtkTube surface。

即：

junction
    ↓
[ IMPLICIT CORE ]
    ↓
[ TRANSITION COLLAR ]
    ↓
[ EXPLICIT TUBE ]

============================================================
十一、不要在复杂 junction core 做 stitching
===========================================

所有 transition 都必须位于：

single-tube-like region。

现有 collar 条件继续保留：

远离 junction
空间分离
没有下一个 bifurcation
没有 CUT_PORT conflict

============================================================
十二、优先方法：Boundary-loop stitching
=======================================

正式 v5 不再优先依赖：

implicit patch 与 explicit tube overlap 后 Boolean union

来生成最终 seam。

推荐实现：

【clip → boundary loops → resample → phase align → stitch】

============================================================
十三、在 collar 处切开两套 surface
==================================

对 branch i：

定义 collar plane：

Π_i

plane normal：

local branch tangent t_i

---

将：

local implicit junction surface

在 collar plane 外侧裁掉。

保留：

junction-side implicit surface。

得到一个 boundary loop：

L_i^implicit

---

将：

explicit branch tube

在 collar plane 内侧裁掉。

保留：

branch-side explicit surface。

得到：

L_i^explicit

============================================================
十四、提取 Boundary Loops
=========================

对：

L_i^implicit

和：

L_i^explicit

分别检查：

必须：

single closed loop

不得：

多个 loop
open loop
self-intersection

否则：

COLLAR_LOOP_EXTRACTION_FAILED

============================================================
十五、统一 loop vertex count
============================

两条 loop 沿周长参数化：

u ∈ [0,1)

重新采样为：

N_loop

个点。

建议：

N_loop = tube_sides

默认：

32

如果 convergence 使用：

48

则同步使用 48。

============================================================
十六、Loop phase alignment
==========================

两个圆周 loop 即使形状相似，
起始 vertex 可能不同。

必须寻找 circular shift：

k*

使：

Σ_j
||p_j^implicit -
p_{j+k}^explicit||²

最小。

同时检查是否需要：

reverse ordering

以保持一致 surface orientation。

============================================================
十七、Transition strip
======================

不要只用一排 triangles 直接连接两个 loops。

建议构建：

N_transition

个中间 rings。

例如：

N_transition = 4～8

配置化。

---

令：

α_m ∈ [0,1]

中间 ring：

p_m,j =
(1-w_m) p_implicit,j
+
w_m p_explicit,j

但不要只做线性 interpolation。

使用 smoothstep：

w(t) =
3t² - 2t³

或者 quintic smoothstep：

w(t) =
6t^5 - 15t^4 + 10t^3

使 transition 两端的一阶变化更平滑。

============================================================
十八、更好的 transition direction
=================================

如果条件允许，

不要仅插值 xyz。

沿：

branch-local frame

分解成：

center
radial vector
radius
angle

分别插值。

目标是避免：

loop center shift

和：

cross-section twisting

产生锯齿。

============================================================
十九、Transition strip triangulation
====================================

相邻 rings：

ring_m
ring_m+1

采用一致的 quad strip。

每个 quad 再按稳定规则 triangulate。

避免：

随机 diagonal orientation

造成 zipper-like triangles。

所有 strip triangle：

orientation 必须一致。

============================================================
二十、不要在 transition strip 内产生 cap
========================================

最终 hybrid surface：

implicit junction patch
+
transition strip
+
explicit branch surface

之间是：

surface stitching

而不是 closed-solid + cap Boolean。

因此：

transition 内不应存在任何 internal cap。

============================================================
二十一、如果 loop stitching 实现失败
====================================

保留当前：

Manifold Boolean overlap

作为：

fallback

但不是正式默认。

配置：

hybrid_transition:
    backend: loop_stitch

    fallback_backend: manifold_boolean

============================================================
二十二、Junction Surface Normal QC
==================================

新增：

local dihedral angle

对共享 edge 的两个 triangles：

n_1
n_2

计算：

θ =
acos(n_1 · n_2)

---

分别统计：

whole surface

junction core

transition collar

explicit branch

输出：

normal_jump_mean
normal_jump_p95
normal_jump_p99
normal_jump_max

============================================================
二十三、当前 28.94° P99 不能再只记录
=====================================

当前：

normal jump P99 ≈ 28.94°

是明显需要优化的 surface-quality signal。

v5 必须比较：

before vs after。

不要一开始人为规定最终必须低于多少度。

先通过：

synthetic tube
synthetic Y
真实 J49

得到正常分布。

随后确定正式 QC threshold。

============================================================
二十四、Surface roughness QC
============================

增加：

mean curvature

或离散 Laplacian magnitude

只用于：

surface roughness comparison。

计算：

transition collar

相对于：

nearby explicit tube

的 roughness ratio。

目的：

识别视觉上的锯齿 seam。

============================================================
二十五、Local remeshing
=======================

stitch 后只允许：

TRANSITION COLLAR

局部 remesh。

不要 remesh 整个 ROI。

推荐评估：

pyvista / vtk
或
trimesh / pymeshlab
或
CGAL-compatible Python package

可使用成熟 package。

---

目标：

triangle size 在 transition 区域不要剧烈变化。

不得改变：

collar cross-sectional area
source topology

============================================================
二十六、Local smoothing
=======================

如果 stitching 后仍存在局部 faceting，

允许：

transition-only
constrained smoothing

例如：

Taubin smoothing
或
Windowed Sinc smoothing

但是：

explicit-region vertices fixed

implicit-core deep vertices fixed

只允许中间 transition rings 移动。

============================================================
二十七、Smoothing 必须满足 volume/radius constraint
===================================================

每次 smoothing 后重新计算：

collar area

equivalent radius

要求：

radius fidelity 不明显恶化。

保存：

before_smoothing
after_smoothing

============================================================
二十八、Junction Core 不要过度 smoothing
========================================

v4 已证明：

J49 core topology 正确。

不要为了让图片圆滑，
大范围 smooth junction core。

尤其不要导致：

bifurcation geometry
人工变圆或收缩。

============================================================
二十九、Visual normals 与 Geometry normals 分开
===============================================

最终 CFD geometry：

保存真实 mesh。

另外为 visualization：

重新计算 point normals

使用：

consistent normals
auto-orient normals

允许：

smooth shading。

但不要让 visualization normal 操作改变：

vertex coordinates。

输出：

lumen_surface_cfd.vtp

和：

lumen_surface_visualization.vtp

============================================================
三十、STL 特别说明
==================

STL 不保存 smooth vertex normals。

因此：

不同 STL viewer

可能显示不同 faceting。

正式 geometry QC：

不要以 STL smooth shading 作为唯一依据。

主要依据：

VTP geometry
triangle coordinates
dihedral angle
cross-section
CFD meshing result

============================================================
三十一、v5 Comparison
=====================

必须生成：

v4 vs v5

同一个 J49：

对比：

self intersection
internal faces
internal caps
boundary edges
non-manifold
components

radius P95

collar radius error

normal jump P95
normal jump P99

transition roughness

triangle count

runtime

============================================================
三十二、必须保持 v4 所有 topology PASS
======================================

v5 不允许破坏：

self intersections = 0
internal faces = 0
internal caps = 0
boundary edges = 0
nonmanifold edges = 0
degenerate triangles = 0
components = 1

如果 surface continuity 改善，
但上述任何一项重新失败：

v5 = FAIL

============================================================
三十三、Port v5 验收
====================

重点：

原来视觉上的 circular seam

应明显消失。

同时：

area error 不明显恶化。

输出：

port_transition_area_profile.png

port_normal_profile.png

port_v4_v5_comparison.png

============================================================
三十四、Junction v5 验收
========================

重点检查：

当前图片中类似锯齿的 seam。

必须输出：

junction_transition_wireframe.png

junction_transition_normals.png

junction_v4_v5_same_camera.png

============================================================
三十五、Synthetic Tests
=======================

新增：

curved tube → straight extension

重点测试 port seam。

implicit Y junction
→ explicit straight branches

测试 collar stitching。

acute Y junction

different radius daughters

curved daughter branch transition

============================================================
三十六、正式优先级
==================

优先完成：

P0:
Port continuous-centerline extension

原因：

最简单，
且可以从构造层面彻底删除独立 cylinder seam。

---

P1:
Junction loop-stitch transition

---

P2:
transition local remesh / constrained smoothing

只有 P1 后仍有明显 surface roughness 才启用。

============================================================
三十七、禁止事项
================

不要：

重新改 J13 context extension

重新抽 ROI

修改 source radius

修改 SWC topology

恢复 junction sphere

全局 smoothing

全局 remesh

因为视觉接缝就使用 global implicit

通过放宽 topology QC 换取光滑表面

============================================================
三十八、Definition of Done
==========================

只有同时满足：

1. Port extension 与 branch 使用连续 centerline tube；
2. 原 CORE CUT_PORT 不再作为 branch/cylinder Boolean seam；
3. final CFD port 仍为 flat planar cap；
4. port radius fidelity 保持；
5. local implicit junction 保留；
6. junction collar 默认采用 loop stitching；
7. transition 不含 internal cap；
8. self intersection = 0；
9. internal faces = 0；
10. internal caps = 0；
11. boundary edges = 0；
12. nonmanifold edges = 0；
13. components = 1；
14. radius P95 不明显差于 v4；
15. collar radius error 不明显差于 v4；
16. normal jump P99 相比 v4 明显降低；
17. 当前截图中的 saw-tooth collar seam 明显改善；
18. current port circular seam 明显改善；
19. synthetic tests PASS；
20. 所有结果保存到：
    ulm_3D_vascular/outputs

才可以称为：

v5 CFD-ready smooth-surface PASS

============================================================
最终需要向我汇报
================

请明确报告：

Port：

v4 area error
v5 area error

v4 normal jump
v5 normal jump

separate cylinder 是否已经彻底取消？

---

Junction：

v4 normal jump P95/P99
v5 normal jump P95/P99

transition triangle count

transition roughness before/after

collar radius error before/after

---

Topology：

self-intersection
internal face
internal cap
boundary edge
nonmanifold

是否仍全部为 0？

---

Geometry：

radius P95 before/after

surface volume before/after

triangle count before/after

runtime before/after

---

最后回答：

当前结果是否已经同时满足：

topological validity
radius fidelity
surface continuity

三个 CFD-ready 条件？
