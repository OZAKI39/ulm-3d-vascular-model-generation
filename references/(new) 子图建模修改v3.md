
当前 CFD lumen 根因诊断已经完成，请基于诊断证据修改现有模型生成方法。

本轮任务不是重新设计整个 SWC → CFD lumen pipeline，而是：

【保留目前表现良好的 explicit branch / port reconstruction，只替换存在明确缺陷的 bifurcation/junction reconstruction。】

最终采用：

Explicit branch tubes
        +
Explicit port extensions
        +
Local implicit junction patches
        ↓
Hybrid vascular lumen
        ↓
CFD geometry QC

============================================================
一、已经确认的诊断结果
======================

必须以以下实际结果作为修改依据，不要重新猜测根因。

---

1. CUT_PORT / extension

---

三个 CUT_PORT：

Port 0:
branch radius = 1.583824 μm
extension radius = 1.582292 μm
area difference = 0.1936%
overlap volume = 12.485381 μm³

Port 1:
branch radius = 1.170489 μm
extension radius = 1.169769 μm
area difference = 0.1230%
overlap volume = 4.837990 μm³

Port 2:
branch radius = 1.148660 μm
extension radius = 1.149154 μm
area difference = 0.0861%
overlap volume = 4.727826 μm³

并且：

endpoint position error = 0
endpoint radius error = 0
CUT_PORT radius interpolation error = 0

branch sides = 32
extension sides = 32

extension 均存在正体积 overlap。

VTK radius scalar synthetic test：

target = 2 μm
measured = 1.9935765 μm

通过。

因此：

【禁止重写 CUT_PORT reconstruction。】

目前 Port 的视觉阶梯主要属于：

PORT_NORMAL_OR_SHADING_ARTIFACT

不是实际 lumen radius discontinuity。

---

2. Junction

---

当前 explicit model：

watertight = true
surface components = 1
boundary edges = 0
non-manifold edges = 0
flipped faces = 0

但是：

self-intersection pairs = 74
ray-confirmed internal faces = 67
confirmed internal-cap faces = 50

全部严重问题集中在 junction 49。

两个 junction solid 与 adjacent tubes 的 intersection volume 均约为 junction solid volume 的 50%。

因此：

INSUFFICIENT_JUNCTION_OVERLAP

不是主因。

真正 root causes 已确认：

BOOLEAN_SELF_INTERSECTION
BOOLEAN_INTERNAL_SURFACE
INTERNAL_CAP_SURVIVED_BOOLEAN
BOOLEAN_SLIVER_TRIANGLES

同时 junction 49 存在严重 artificial bulging：

Branch 1:
A(1D) = 42.7814 μm²
A(2D) = 4.1748 μm²

Branch 4:
A(1D) = 33.6433 μm²
A(2D) = 3.9185 μm²

因此：

【当前 sphere + tube + Boolean junction construction 不再作为正式默认方案。】

---

3. Explicit vs implicit

---

Explicit:
self intersections = 74
internal faces = 67
radius P95 error = 0.3220%
port area step = 0.1758%
triangles = 30,126
runtime = 0.792 s
degenerate triangles = 0

Implicit:
self intersections = 0
internal faces = 0
radius P95 error = 1.0036%
port area step = 0.0365%
triangles = 92,948
runtime = 55.223 s
degenerate triangles = 6

因此不能：

完全保留 explicit junction；

也不能：

把整个 ROI 全部切换为 global implicit。

正式修改目标：

【Explicit outside junction + Local implicit at junction】

============================================================
二、总体方法
============

新的正式 reconstruction pipeline：

SWC ROI
   ↓
branch extraction
   ↓
explicit variable-radius branch tubes
   ↓
explicit CUT_PORT extensions
   ↓
detect bifurcation neighborhoods
   ↓
REMOVE sphere-based junction construction
   ↓
construct local implicit junction field
   ↓
extract local junction surface
   ↓
merge with explicit branches in simple collar regions
   ↓
hybrid watertight lumen
   ↓
existing CFD geometry QC

即：

Ω =
Ω_explicit_branches
∪
Ω_local_implicit_junctions
∪
Ω_port_extensions

============================================================
三、最重要的修改原则
====================

1. 保留当前 branch explicit tube implementation。
2. 保留当前 CUT_PORT extension implementation。
3. 不再使用 junction sphere 作为正式 junction primitive。
4. 不再在复杂 bifurcation 中直接做：

parent tube
+
daughter tube
+
sphere
→ Boolean

5. implicit reconstruction 只作用于 bifurcation 附近的小局部区域。
6. 不允许局部 implicit field 修改 SWC topology。
7. 非邻接 branch collision QC 必须继续在 reconstruction 前执行。
8. 不修改 source SWC。
9. 不修改 source radius。
10. 不为了消除 junction 问题 aggressive smoothing centerline。

============================================================
四、代码组织
============

继续使用：

ulm_3D_vascular/model_generate.py

作为 workflow / CLI。

不要堆代码进入 model_generate.py。

在：

ulm_3D_vascular/utils/cfd_lumen/

现有模块基础上修改。

建议新增：

local_implicit_junction.py

如果已有：

junction reconstruction module

则合理拆分/重构，不要建立重复逻辑。

建议最终结构：

geometry_preprocess.py
collision_qc.py
branch_tube.py
port_geometry.py

local_implicit_junction.py
hybrid_merge.py

surface_qc.py
geometry_diagnostics.py
export.py
visualization.py

============================================================
五、删除正式流程中的 Junction Sphere
====================================

找到当前：

junction sphere

生成代码。

保留原 backend 用于：

legacy_explicit

诊断对照。

但正式默认：

junction_backend = local_implicit

配置：

junction:
    backend: local_implicit

同时允许：

junction:
    backend: legacy_sphere

仅用于 regression test。

不要删除旧方法，
但它不能再作为正式默认结果。

============================================================
六、定义 Junction Neighborhood
==============================

对于每一个真实 bifurcation node：

J

获取：

incident branches:

B_1 ... B_m

以及 junction coordinate：

x_J

和 SWC junction radius：

r_J

对每一条 incident branch，

沿 centerline 从 junction 向 branch 内部寻找一个：

COLLAR point

记为：

c_i

============================================================
七、Collar 的意义
=================

Junction 附近：

使用 implicit reconstruction。

远离 junction：

继续使用高精度 explicit tube。

因此每一条 branch 都需要定义一个过渡位置：

junction
   ↓
implicit region
   ↓
COLLAR
   ↓
explicit branch region

示意：

explicit branch
====================|
                    | collar
                 ___|___
                /       
               / implicit
              / junction  
                region

============================================================
八、Collar Distance
===================

第一版定义：

L_collar,i =
λ × D_i

其中：

D_i = 2 r_i

建议默认：

junction_collar_diameters = 2.0

即：

L_collar ≈ 2D

但是必须配置化。

测试：

1.5D
2D
3D

---

collar 位置必须满足：

1. 已经离开 junction core；
2. branch cross-section 基本稳定；
3. 不包含下一个 bifurcation；
4. 不越过 CUT_PORT；
5. 相邻 incident branches 在 collar 位置已经空间分离。

如果 2D 位置不满足上述条件，

允许自动向外扩展：

直到：

min(branch separation)

local radii sum + tolerance

但不得越过下一个 topology node。

============================================================
九、Local Implicit Field：不要使用简单 Junction Sphere
======================================================

局部 implicit geometry 应直接从：

SWC centerline + radius

生成。

正式实现采用：

【variable-radius polyball / swept-sphere field】

而不是：

junction sphere。

============================================================
十、Polyball Field
==================

对于 centerline sample：

c(s)

及：

r(s)

定义局部 tube field：

φ(x)
=====

min_s
[
||x - c(s)|| - r(s)
]

其中：

φ < 0:
inside lumen

φ = 0:
vessel wall

φ > 0:
outside lumen

---

注意：

不要只使用单个 SWC junction node 做 sphere。

必须使用：

parent centerline segment
+
all daughter centerline segments

共同构造 local field。

============================================================
十一、Local implicit 输入范围
=============================

每个 junction 只使用：

junction
→
各 branch collar

之间的 centerline。

另外为了保证 merge overlap，

每条 branch 再向 collar 外侧延长：

L_overlap

建议：

junction_overlap_diameters = 0.5

即：

local implicit patch 比 collar 多覆盖约 0.5D 的 explicit branch。

这样 hybrid merge 发生在：

简单、近似单管的区域，

而不是复杂 bifurcation core。

============================================================
十二、为什么这样做
==================

当前 Boolean failure 出现在：

复杂的三管/多管 junction core。

新的方法将：

复杂 junction

交给 implicit reconstruction。

Boolean/merge 只发生在：

single tube-like collar region。

这样可以避免：

sphere
+
multiple capped tubes

同时 Boolean 所产生的：

internal caps
self intersections
sliver triangles。

============================================================
十三、Local Grid
================

禁止重新建立整个：

80 × 80 × 120 μm

ROI 的高分辨率 implicit grid。

每一个 junction 单独建立：

local bounding box。

bbox 只覆盖：

junction center
+
centerline to collars
+
最大 local radius
+
padding

建议：

padding =
1.5 × max local radius

配置化。

============================================================
十四、Local Grid Resolution
===========================

全局 implicit 的 radius P95 error 是：

1.0036%

主要原因之一是 grid discretization。

因为现在只做局部 junction，

可以采用更高 resolution。

定义：

h =
D_min / N_D

其中：

D_min

是 junction incident branches 中的最小 diameter。

建议测试：

N_D = 12
16
20
24

默认第一版：

cells_across_min_diameter = 16

不要硬编码。

============================================================
十五、性能优化
==============

Local implicit field 必须高效。

使用：

float32

不要默认 float64。

centerline samples 使用：

scipy.spatial.cKDTree

或高效 segment spatial query。

field 分块计算。

禁止：

每个 voxel
×
全部 centerline samples

的 Python nested loop。

优先：

NumPy vectorization
Numba
cKDTree candidate search

============================================================
十六、不要直接使用全局 smooth-min
=================================

第一版：

标准 union：

φ(x) =
min_i φ_i(x)

作为 baseline。

如果 junction 表面明显存在数学 sharp crease，

再实现：

smooth-min

作为可选模式。

例如：

smooth_union = false

默认先 false。

不要一开始引入额外 smoothing 导致：

artificial bulging。

============================================================
十七、Marching Cubes
====================

使用：

skimage.measure.marching_cubes

在：

φ = 0

提取 local junction surface。

输出：

junction_<id></id>_implicit_raw.vtp

并保存局部 field metadata：

grid spacing
grid dimensions
bbox
runtime

============================================================
十八、局部 Implicit Patch 不能带 PORT
=====================================

CUT_PORT reconstruction 完全使用现有 explicit 方法。

Local implicit patch 只处理：

TRUE BIFURCATION REGION

不得因为 junction patch 扩展到 CUT_PORT。

如果：

junction collar
与
CUT_PORT extension region

发生 overlap，

说明 ROI 太短或参数不合理。

应报告：

JUNCTION_PORT_REGION_CONFLICT

而不是静默合并。

============================================================
十九、Hybrid Merge
==================

这是整个修改最关键的步骤之一。

不要让 implicit patch 和 explicit branch：

恰好只接触。

必须在 collar 外：

存在明确 overlap：

L_overlap > 0

---

每个 branch：

implicit patch
==============

           || overlap ||
                    ================= explicit branch

Boolean union 只发生在这些：

simple tubular overlap regions。

============================================================
二十、显式 branch 的内部 cap 问题
=================================

当前诊断已经发现：

INTERNAL_CAP_SURVIVED_BOOLEAN。

因此必须专门处理。

不要继续在真实 junction node 上 cap explicit branch。

新的 explicit branch primitive：

应至少延伸到：

collar - overlap

或：

collar region inside implicit solid。

它的 closed end cap 必须：

【深埋在 local implicit lumen 内部】

而不是位于 junction 表面附近。

---

Boolean 后必须再次检查：

internal cap faces = 0

如果：

> 0

则：

FAIL

============================================================
二十一、推荐两种 Hybrid Merge 实现并比较
========================================

优先实现：

Method A:
manifold Boolean union

因为现在 Boolean 只处理简单 single-tube overlap。

如果仍出现 internal cap：

再实现：

Method B:
collar loop stitching

不要直接一开始做复杂 stitching。

============================================================
二十二、Method A：Local Boolean Merge
=====================================

将：

explicit branch solids
+
local implicit junction solid

进行 manifold union。

但必须：

按 junction 单独 merge。

例如：

junction 49 patch
∪
its incident branch tubes

完成后，

再进入整体 vascular surface。

不要：

全部 ROI 30 个 primitives
一次性大 Boolean。

============================================================
二十三、Boolean 顺序
====================

每个 junction：

implicit junction patch

作为 base。

逐个 incident branch 或 balanced tree union。

记录：

before
after

每一步：

triangles
volume
self-intersection
internal-face count

便于定位问题。

============================================================
二十四、Local Remeshing
=======================

Implicit patch marching cubes 会产生较多 triangles。

不要全局 remesh。

只允许：

local junction patch

进行：

clean
remove degenerate faces
optional isotropic remesh

目标：

减少 sliver triangles。

============================================================
二十五、禁止 aggressive mesh smoothing
======================================

如果要 smooth：

只允许 local junction。

且必须设置：

boundary collar vertices fixed

避免：

explicit branch radius
被拉动。

任何 smoothing 前后必须重新测：

collar radius fidelity。

============================================================
二十六、Port 阶梯只修改显示，不修改 CFD geometry
================================================

当前 Port 数据已经证明几何连续。

因此不要修改：

CUT_PORT radius
extension radius
overlap
extension length

如果希望验收图片视觉更平滑，

可以：

recompute point normals

并对 visualization-only VTP：

smooth_shading=True

同时保留：

CFD STL

原几何不变。

---

可以增加：

tube_sides = 48

作为 optional geometry convergence test。

但不要因为截图有台阶就默认提高。

当前 32 sides 已经产生：

<0.2%

port area mismatch。

因此：

32 remains valid baseline。

============================================================
二十七、新增 Hybrid-specific QC
===============================

现有 QC 全部保留。

另外必须增加：

---

1. Internal face

---

要求：

internal_face_count = 0

---

2. Internal cap

---

要求：

internal_cap_face_count = 0

---

3. Self intersection

---

要求：

self_intersection_pairs = 0

---

4. Junction degenerate triangles

---

要求：

degenerate_triangle_count = 0

或经过明确清理后为 0。

---

5. Boundary edge

---

0

---

6. Non-manifold edge

---

0

---

7. Surface component

---

1

============================================================
二十八、Collar Radius Fidelity
==============================

对于每个 junction incident branch，

在：

collar - 0.5D
collar
collar + 0.5D

分别测量：

A(s)

定义：

r_eq =
sqrt(A/π)

比较 source：

r_SWC

输出：

collar_radius_error.csv

---

Hybrid reconstruction 必须证明：

implicit → explicit

过渡处不存在新的：

radius step。

============================================================
二十九、Junction Area Profile
=============================

当前 explicit junction 出现：

42.78 μm²
vs
4.17 μm²

这种异常膨大。

新版本必须沿每条 incident branch：

从 3D away
→
junction
→
3D away

连续采样：

A(s)

输出：

junction_area_profile.csv

以及：

junction_area_profile.png

不要只取：

1D
2D

两个点。

---

图中标记：

junction center
collar
implicit region
explicit region

这样可以直观看出：

是否出现人工 abrupt bulge。

============================================================
三十、Junction Volume Metric
============================

保存：

local junction volume

对比：

legacy sphere explicit
global implicit
new local implicit hybrid

不要设置生理硬阈值。

只用于 A/B comparison。

============================================================
三十一、Local Grid Convergence
==============================

选择当前最有问题的：

junction 49

进行：

N_D =
12
16
20
24

四档测试。

比较：

runtime
triangles
self-intersection
internal faces
radius fidelity
collar radius error
local volume
area profile

目标：

找到：

geometry quality
vs
runtime

的平衡。

============================================================
三十二、正式 Comparison
=======================

必须在同一个真实 ROI 上比较：

A. legacy explicit sphere
B. global implicit
C. new hybrid local implicit

保存：

reconstruction_comparison.csv

至少包括：

method

runtime

triangles

watertight

boundary_edges

nonmanifold_edges

self_intersections

internal_faces

internal_caps

degenerate_triangles

radius_P95_error

port_area_step

junction_max_area_ratio

surface_volume

============================================================
三十三、预期目标
================

新的 hybrid method 目标不是：

在所有指标上绝对优于 explicit 和 global implicit。

而是：

保留 explicit 的：

高速
低 radius error
低 triangle count
高 port fidelity

同时获得 implicit junction 的：

0 self-intersection
0 internal surface
0 internal cap

即希望达到：

runtime
远低于 55 s

同时：

self intersections = 0
internal faces = 0
internal caps = 0

============================================================
三十四、输出
============

继续保存到：

ulm_3D_vascular/outputs/model_generate/<run_id>/

新增：

hybrid/
    junctions/
        junction_49/
            centerline.vtp
            local_field_metadata.json
            implicit_raw.vtp
            implicit_clean.vtp
            merged.vtp

    figures/
        junction_49_field.png
        junction_49_implicit.png
        junction_49_hybrid.png
        junction_49_area_profile.png
        collar_transition.png
        legacy_vs_implicit_vs_hybrid.png

    qc/
        hybrid_surface_qc.json
        collar_radius_error.csv
        junction_area_profile.csv
        reconstruction_comparison.csv

============================================================
三十五、Synthetic tests
=======================

新增：

Test 1:
symmetric Y junction

Test 2:
asymmetric Y junction
different daughter radii

Test 3:
acute-angle bifurcation

Test 4:
large radius-ratio bifurcation

Test 5:
three-way junction

Test 6:
hybrid collar transition

每个 test 要求：

watertight
connected
no boundary edges
no nonmanifold
no self-intersection
no internal faces
no internal caps

============================================================
三十六、真实 ROI smoke test
===========================

必须优先测试当前：

junction 49

因为它已经被诊断为：

74 self-intersections
67 internal faces
50 internal-cap faces

这个 junction 是最好的 regression case。

修改成功以后要求：

junction 49:

self_intersections = 0
internal_faces = 0
internal_caps = 0
boundary_edges = 0
nonmanifold_edges = 0

并输出：

before/after

对比。

============================================================
三十七、第一阶段不要删除 global implicit
========================================

global implicit backend 保留。

用途：

reference reconstruction
fallback
A/B comparison

正式配置改为：

reconstruction:
    branch_backend: explicit
    port_backend: explicit
    junction_backend: local_implicit

不要删除：

global_implicit

和：

legacy_explicit

============================================================
三十八、推荐 config
===================

junction:
  backend: local_implicit

  collar_diameters: 2.0

  overlap_diameters: 0.5

  bbox_padding_radius: 1.5

  implicit:
    cells_across_min_diameter: 16
    dtype: float32

    smooth_union: false

  remesh:
    enabled: false

hybrid_merge:
  backend: manifold

surface_qc:
  require_watertight: true
  require_single_component: true
  require_zero_boundary_edges: true
  require_zero_nonmanifold_edges: true

  require_zero_self_intersections: true
  require_zero_internal_faces: true
  require_zero_internal_caps: true

不要把参数散落在代码里。

============================================================
三十九、Definition of Done
==========================

本轮只有同时满足下面条件才算完成：

1. 正式 pipeline 不再使用 junction sphere；
2. branch explicit reconstruction 保留；
3. CUT_PORT explicit reconstruction 保留；
4. junction 使用 local implicit；
5. global implicit 仍保留作为 comparison/fallback；
6. junction 49 self-intersection 从 74 降到 0；
7. internal faces 从 67 降到 0；
8. internal caps 从 50 降到 0；
9. boundary edges = 0；
10. nonmanifold edges = 0；
11. surface components = 1；
12. no source topology modification；
13. collar radius fidelity 被测量；
14. junction area profile 被保存；
15. port area continuity 没有比原 explicit 明显恶化；
16. runtime 明显低于 global implicit 55.223 s；
17. comparison table 完整输出；
18. synthetic junction tests 通过；
19. 当前真实问题 ROI 完成 before/after validation；
20. 所有日志和可视化结果保存在：
    ulm_3D_vascular/outputs

============================================================
四十、最终向我汇报
==================

不要只说“hybrid 已完成”。

需要明确报告：

1. junction 49：

before:
self-intersections = 74
internal faces = 67
internal caps = 50

after:
分别是多少？

2. runtime：

legacy explicit =
global implicit =
new hybrid =

3. triangle count：

三种方法分别是多少？

4. radius fidelity：

legacy explicit P95 =
global implicit P95 =
hybrid P95 =

5. port area step：

是否保持在原有约 0.2% 以内？

6. junction area profile：

原来的巨大 bulge 是否消失？

7. collar transition：

是否产生新的 step？

8. 是否仍存在：

degenerate triangles
self intersection
internal surface
non-manifold
boundary edge

9. 最终推荐：

是否可以把：

explicit branch
+
local implicit junction
+
explicit port

设为正式默认 backend？

必须根据实际运行指标回答。

============================================================
核心目标
========

不要继续优化：

“tube + sphere + Boolean”

如何看起来更漂亮。

当前诊断已经证明：

问题集中在该 junction construction。

新的正式方法应该是：

【保留 explicit 方法在简单管状区域的效率和管径精度，
只把复杂 bifurcation 区域交给 local implicit reconstruction，
并在远离 junction 的简单 collar 区域与 explicit tube 合并。】

即：

Explicit where explicit works well.
Implicit only where implicit is necessary.
