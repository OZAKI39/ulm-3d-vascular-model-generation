
当前任务不是继续增加 CFD 功能，也不是立即重写 lumen reconstruction。

当前任务是：

【针对现有 SWC ROI → CFD lumen 生成结果中出现的几何异常，完成系统性的 root-cause diagnosis，并用实际中间数据、定量指标和可视化结果证明异常究竟发生在哪一个处理阶段。】

目前人工检查最终 STL 时发现两类明显异常：

问题 A：
CUT_PORT / extension 与原血管连接处出现明显环形“台阶”或 shoulder。

问题 B：
bifurcation / junction 区域出现明显异常尖角、三角面畸变和疑似裂缝/接缝。

不要根据截图直接猜原因。

必须通过代码检查、几何测量、中间 mesh 导出和自动 QC，
逐层确定真正的 root cause。

============================================================
一、项目环境
============

项目：

ulm_3D_vascular

Python interpreter：

d:\anaconda3\envs\pmp\python.exe

所有命令必须使用：

d:\anaconda3\envs\pmp\python.exe

例如：

d:\anaconda3\envs\pmp\python.exe -m pip ...
d:\anaconda3\envs\pmp\python.exe ulm_3D_vascular\model_generate.py ...

允许安装必要的高效包。

优先复用当前已经安装的：

numpy
scipy
vtk
pyvista
trimesh
manifold3d
matplotlib
pandas

============================================================
二、重要原则
============

1. 本轮首先执行 DIAGNOSIS，不要先改变 source SWC。
2. 不要为了让模型看起来平滑而：

   - 缩小 radius；
   - 自动连接 branch；
   - 删除 branch；
   - aggressive smoothing。
3. 不要在不知道 root cause 的情况下直接把：
   tube_sides = 32
   改成：
   tube_sides = 128
   并宣称问题解决。
4. 必须区分：

   - visual faceting；
   - geometric step；
   - Boolean seam；
   - open crack；
   - non-manifold edge；
   - internal overlapping surface；
   - normal/shading artifact。
5. 所有结论必须由实际运行产生的：
   数字 + mesh + 图片
   支撑。

============================================================
三、代码组织
============

不要把诊断代码全部塞进：

ulm_3D_vascular\model_generate.py

model_generate.py 只增加必要 CLI 调用。

在现有：

ulm_3D_vascular\utils\cfd_lumen\

下增加合理的 diagnostic 模块。

建议：

geometry_diagnostics.py
port_diagnostics.py
junction_diagnostics.py
mesh_defects.py
diagnostic_visualization.py

如果现有文件已经承担相同功能，则扩展现有模块，
不要创建重复功能。

============================================================
四、输出目录
============

所有诊断结果保存到：

ulm_3D_vascular\outputs\model_generate\<run_id>\diagnostics\

每个 ROI：

diagnostics/
    <roi_id>/
        summary.json

        port_diagnostics.csv
        junction_diagnostics.csv
        mesh_defects.json

        primitives/
            branch_tubes/
            junction_solids/
            port_extensions/
            pre_boolean_combined.vtp
            post_boolean_surface.vtp

        figures/
            01_port_step_overview.png
            02_port_cross_sections.png
            03_port_radius_profile.png
            04_extension_overlap.png
            05_junction_primitives.png
            06_junction_post_boolean.png
            07_boundary_edges.png
            08_nonmanifold_edges.png
            09_surface_components.png
            10_normals_overlay.png
            11_triangle_quality.png
            12_defect_summary.png

        logs/
            geometry_diagnostics.log

不要覆盖以前的 model generation 结果。

============================================================
五、第一步：找到产生异常的真实 ROI
==================================

先检查当前 model_generate outputs。

识别用户截图所对应的 ROI，或者从当前 selected ROIs 中定位：

1. 至少一个存在明显 port step 的 ROI；
2. 至少一个存在明显 bifurcation/junction artifact 的 ROI。

不要硬编码 ROI ID。

允许 CLI：

--diagnose-roi ROI_ID

以及：

--diagnose-all

第一轮只跑有问题的 1–2 个 ROI：

workers = 1

保证日志容易追踪。

============================================================
六、问题 A：Port / Extension 阶梯的系统诊断
===========================================

需要回答：

【阶梯究竟是 tube polygon faceting，还是 branch 与 extension 的真实半径/位置不连续？】

对每一个 CUT_PORT 分别检查以下内容。

---

A1. 检查 CUT_PORT 原始数据
--------------------------

输出：

cut_port_id

exact_cut_position
cut_radius

source_global_edge

source_edge_start_position
source_edge_end_position

source_edge_start_radius
source_edge_end_radius

检查：

cut_radius

是否确实按照 source edge 两端 radius 线性插值得到。

重新独立计算：

t =
projection parameter of cut point on source edge

r_recomputed =
(1-t) r0 + t r1

计算：

radius_interpolation_error =
|r_saved - r_recomputed|

结果写入：

port_diagnostics.csv

---

A2. 检查 resampling 是否保留精确 CUT_PORT endpoint
--------------------------------------------------

这是重点排查项。

对于 CUT_PORT 对应 branch，检查：

branch endpoint position：

x_branch_end

是否满足：

||x_branch_end - x_cut|| < tolerance

并检查：

r_branch_end

是否满足：

|r_branch_end - r_cut| < tolerance

输出：

endpoint_position_error_um
endpoint_radius_error_um
endpoint_radius_relative_error

特别检查现有 resampling 是否使用：

np.arange(...)

导致最后一个 endpoint 被遗漏。

如果 endpoint 不在 resampled samples 中，
明确记录为：

ROOT_CAUSE_CANDIDATE:
RESAMPLING_DROPPED_CUTPORT_ENDPOINT

不要静默修复。

---

A3. 检查 VTK tube 实际 endpoint radius
--------------------------------------

不能只检查输入 scalar。

必须从真正生成的 branch tube 几何中，
在靠近 CUT_PORT 处测量实际截面积。

建立垂直于 local tangent 的 plane：

P_branch

在：

x_cut - ε n_out

位置与 branch tube 相交。

计算：

A_branch

对应 equivalent radius：

r_branch_mesh =
sqrt(A_branch / π)

---

A4. 检查 extension 实际 radius
------------------------------

在 extension 内部：

x_cut + ε n_out

建立截面。

计算：

A_ext

以及：

r_ext_mesh =
sqrt(A_ext / π)

比较：

r_source_cut
r_branch_mesh
r_ext_config
r_ext_mesh

输出：

source_radius_um
branch_mesh_radius_um
extension_radius_input_um
extension_mesh_radius_um

---

A5. 定量计算“台阶”
--------------------

定义：

E_radius =
|r_branch_mesh - r_ext_mesh|
/
[(r_branch_mesh + r_ext_mesh)/2]

以及：

E_area =
|A_branch - A_ext|
/
[(A_branch + A_ext)/2]

所有 port 都计算。

如果视觉上有台阶，而：

E_radius ≈ 0
E_area ≈ 0

则说明更可能是：

polygon faceting
normal/shading
mesh seam

而不是实际 lumen 半径跳变。

---

A6. 检查 tube_sides 是否一致
----------------------------

记录：

branch_tube_sides
extension_cylinder_sides

如果两者不同：

例如：

branch = 32
extension = 64

则可能产生可见环形接缝。

必须记录：

cross_section_vertex_count_branch
cross_section_vertex_count_extension

但不要直接认定这是 root cause，
因为 Boolean union 后理论上仍可能形成正确几何。

---

A7. 检查 radius 语义是否用错
----------------------------

重点检查当前代码使用：

vtkTubeFilter

时：

SetRadius
SetVaryRadiusToVaryRadiusByAbsoluteScalar
scalar array

之间的真实语义。

必须检查实际代码：

输入的 scalar 到底表示：

radius

还是：

radius multiplier / diameter / scale factor。

用一个 synthetic straight vessel：

target radius = 2.0 μm

生成 tube，

然后实际测量截面：

r_mesh

要求确认：

r_mesh ≈ 2.0 μm

如果出现：

1 μm
4 μm
或其他固定比例，

则定位为：

ROOT_CAUSE:
VTK_RADIUS_SCALAR_SEMANTICS

---

A8. 检查 extension 与 branch 是否真正重叠
-----------------------------------------

根据配置：

port_overlap_diameters

重新计算理论 overlap：

L_overlap_expected

然后根据实际 primitive bounding / centerline geometry 检查：

L_overlap_actual

检查：

extension start point
branch endpoint
branch tube volume

是否存在真实正体积 overlap。

输出：

expected_overlap_um
actual_axial_overlap_um

并用 intersection volume / manifold Boolean test 判断：

V_intersection(branch, extension)

必须报告：

intersection_volume

如果：

intersection_volume ≈ 0

说明两个 solids 只是接触，而不是重叠。

标记：

ROOT_CAUSE_CANDIDATE:
INSUFFICIENT_PORT_OVERLAP

---

A9. 导出 Boolean 前 primitive
-----------------------------

保存：

branch_tube_<id></id>.vtp
port_extension_<id></id>.vtp

并生成：

04_extension_overlap.png

图片必须显示：

branch tube 半透明
extension 半透明
CUT_PORT
branch endpoint
extension start
extension cap
local tangent

让人工能够看见：

两个 solids 是否真正 overlap。

============================================================
七、问题 B：Bifurcation / Junction 裂缝与异常面的系统诊断
=========================================================

需要回答：

【看到的白色缝隙到底是真裂缝、non-manifold seam、内部重叠面、Boolean artifact，还是单纯 normal/shading artifact？】

---

B1. 保存 Boolean 前所有 junction primitives
-------------------------------------------

对于每个 bifurcation：

保存：

parent tube
daughter tubes
junction sphere / junction solid

分别不同 mesh 文件。

生成：

05_junction_primitives.png

要求：

不同 primitive 使用不同显示颜色。

同时显示：

source centerline
junction node
source radius
branch IDs

---

B2. 测量 primitive overlap
--------------------------

对于每一个：

junction sphere ↔ branch tube

计算：

intersection volume

至少分别计算：

V_parent
V_daughter1
V_daughter2
...

输出：

intersection_volume

以及：

intersection_volume /
junction_solid_volume

如果某一 branch：

intersection volume ≈ 0

说明 junction solid 与 branch 几乎只是接触。

标记：

ROOT_CAUSE_CANDIDATE:
INSUFFICIENT_JUNCTION_OVERLAP

---

B3. 检查 junction sphere radius
-------------------------------

输出：

junction SWC radius
junction primitive radius

以及：

各 adjacent branch 在 junction 附近的：

branch radius

计算：

r_junction /
max(r_adjacent)

和：

r_junction /
min(r_adjacent)

不要修改值。

只判断是否出现：

junction sphere 显著小于 adjacent vessel

从而难以形成稳定 overlap。

---

B4. 检查 adjacent tube 是否被 cap
---------------------------------

这是重点。

检查当前 vtkTubeFilter / primitive construction：

branch tube 在 bifurcation 端是否已经存在 flat cap。

如果：

tube 本身 cap=True

然后再与：

junction sphere

做 Boolean，

必须检查该内部 cap 是否被 Boolean 完全去除。

输出：

branch_internal_cap_present_before_boolean = true/false

并保存 cap 面的位置。

如果最终 union 中仍然残留内部 cap：

标记：

ROOT_CAUSE:
INTERNAL_CAP_SURVIVED_BOOLEAN

---

B5. Boolean 前后 triangle count
-------------------------------

对每个 junction neighborhood 记录：

pre_boolean vertices
pre_boolean triangles

post_boolean vertices
post_boolean triangles

以及：

number of solids unioned

Boolean backend
Boolean runtime

Boolean warnings/exceptions

---

B6. Boundary edge detection
---------------------------

最终 mesh 上分别提取：

boundary edges

不要与 feature edges 混在一起。

输出：

boundary_edge_count

如果：

boundary_edge_count > 0

生成：

07_boundary_edges.png

只把 boundary edges 高亮。

并将每个 boundary loop 分组：

loop_id
edge_count
loop_length
centroid

如果用户截图中的白缝与 boundary edges 重合：

结论：

TRUE_OPEN_CRACK

这不是显示问题。

---

B7. Non-manifold edge detection
-------------------------------

独立提取：

non-manifold edges

输出：

non_manifold_edge_count

生成：

08_nonmanifold_edges.png

如果截图位置对应 non-manifold edges：

结论：

NON_MANIFOLD_JUNCTION

---

B8. Connected components
------------------------

最终 surface：

计算：

surface connected component count

并保存每个 component：

area
volume
triangle count

如果 >1：

FAIL

生成：

09_surface_components.png

============================================================
八、检查内部重叠面 / 残留面
===========================

即使：

watertight = true

也可能有 internal surfaces。

必须额外检查。

使用：

trimesh / manifold / ray casting

检查最终 mesh 是否存在：

duplicate coplanar faces
internal faces
self intersections
near-coincident opposing faces

输出：

duplicate_face_count
suspected_internal_face_count
self_intersection_count

如果库无法可靠给出完整 self-intersection count，
至少：

1. 在 junction bounding box 内做 triangle spatial acceleration；
2. 查找不共享 vertex 的 triangle-triangle intersections；
3. 输出 suspect pair count。

不要做全局 O(N²)。

使用：

AABB / R-tree / VTK spatial locator。

============================================================
九、Normal / Shading Artifact 检查
==================================

如果：

boundary_edge_count = 0
non_manifold_edge_count = 0
self_intersection = 0

但视觉仍有“白缝”，

则检查 normals。

输出：

is_winding_consistent
face_normal_consistency
number_of_flipped_faces

生成：

10_normals_overlay.png

显示：

surface
face normals

并重新计算一次 normals：

但不要覆盖正式 mesh。

只生成：

recomputed_normals_preview

如果重新计算 normals 后白缝消失：

结论：

DISPLAY_OR_NORMAL_ARTIFACT

而不是 geometry crack。

============================================================
十、Triangle Quality Diagnosis
==============================

分叉区域的大而尖锐三角面可能导致视觉异常和后续 CFD mesh 问题。

计算每个 triangle：

area
aspect ratio
minimum angle
maximum angle

至少输出：

triangle_area_min
triangle_area_p5
triangle_area_median

aspect_ratio_p95
aspect_ratio_max

min_angle_p5

特别对 junction neighborhood 单独计算。

生成：

11_triangle_quality.png

将：

bad aspect-ratio triangles

高亮。

如果异常位置主要由：

extremely thin sliver triangles

组成，

标记：

ROOT_CAUSE_CANDIDATE:
BOOLEAN_SLIVER_TRIANGLES

============================================================
十一、局部分叉截面诊断
======================

对每个真实 junction，

沿 parent branch：

junction - 2D
junction - 1D

以及各 daughter：

junction + 1D
junction + 2D

取垂直截面。

计算：

cross-section area

equivalent radius

输出：

A_parent_pre
A_parent_near
A_daughter_near
...

同时在 junction 本身附近计算：

A_max_local
A_min_local

目的：

检测：

artificial narrowing
artificial bulge
almost-closed throat

不要先设置生理 threshold，
只报告真实值和相对比值。

============================================================
十二、检查“裂缝”是否来自 viewer/rendering
===========================================

对同一个最终 STL/VTP：

至少使用两种不同 rendering 方式生成图：

1. PyVista smooth_shading=False
2. PyVista smooth_shading=True

并分别输出。

另外输出：

wireframe overlay。

如果：

smooth shading 下有白缝

但：

wireframe
boundary edge
non-manifold

都显示 geometry 连续，

则说明主要是 shading artifact。

============================================================
十三、Synthetic Control Experiments
===================================

为了区分代码问题和真实复杂 ROI 问题，
必须运行最小 synthetic tests。

---

Control 1
---------

Straight tube + straight extension

要求：

无 visible shoulder

测得：

E_radius ≈ 0

---

Control 2
---------

Tapered tube + constant-radius extension

CUT_PORT radius 明确定义为 taper endpoint radius。

确认接口是否连续。

---

Control 3
---------

Simple symmetric Y bifurcation

使用和真实模型相同：

tube backend
junction sphere
Boolean backend

检查是否也出现：

crack / sliver triangle / seam。

如果 synthetic Y 也失败：

说明 reconstruction implementation 有系统问题。

如果 synthetic Y 正常但真实 ROI 失败：

说明问题与：

junction geometry
branch angles
radius ratios
local point spacing

有关。

============================================================
十四、对 explicit backend 和 implicit backend 做一个诊断对照
============================================================

本轮不要全面切换 backend。

只选择同一个有问题的真实 ROI。

生成：

A:
current explicit/manifold geometry

B:
implicit fallback geometry

要求两者使用：

完全相同 SWC centerline
完全相同 radius
完全相同 CUT_PORT

比较：

watertight
non-manifold
boundary edges
surface components
radius fidelity
port step metric
junction visual quality
runtime
triangle count

保存：

backend_comparison.csv

以及：

12_explicit_vs_implicit.png

这一实验的目的不是立即宣布哪个 backend 更好，
而是判断：

异常是否主要来自 explicit tube + sphere + Boolean reconstruction。

如果：

implicit 无裂缝
explicit 有裂缝

则 root cause 很可能位于：

explicit primitive / Boolean junction construction

而不是 source SWC。

============================================================
十五、Root Cause Classification
===============================

最终每个异常必须归入以下一个或多个类别：

PORT_RESAMPLING_ENDPOINT_MISMATCH

PORT_RADIUS_MISMATCH

PORT_TUBE_SCALAR_SEMANTICS

PORT_INSUFFICIENT_OVERLAP

PORT_POLYGON_RESOLUTION_MISMATCH

PORT_NORMAL_OR_SHADING_ARTIFACT

JUNCTION_INSUFFICIENT_OVERLAP

JUNCTION_RADIUS_INCOMPATIBILITY

INTERNAL_CAP_SURVIVED_BOOLEAN

BOOLEAN_OPEN_CRACK

BOOLEAN_NONMANIFOLD

BOOLEAN_INTERNAL_SURFACE

BOOLEAN_SELF_INTERSECTION

BOOLEAN_SLIVER_TRIANGLES

NORMAL_OR_SHADING_ARTIFACT

SOURCE_GEOMETRY_COLLISION

SOURCE_CENTERLINE_GEOMETRY_ISSUE

UNKNOWN

不要只输出：

"mesh bad"

必须给出：

evidence
metric
location
likely root cause
confidence

============================================================
十六、summary.json
==================

每个 ROI 输出类似：

{
  "roi_id": "...",

  "port_issue": {
    "detected": true,
    "affected_ports": [...],
    "root_causes": [...],
    "evidence": {
      "max_endpoint_position_error_um": ...,
      "max_endpoint_radius_error": ...,
      "max_step_area_error": ...,
      "min_overlap_volume": ...
    }
  },

  "junction_issue": {
    "detected": true,
    "affected_junctions": [...],
    "root_causes": [...],
    "evidence": {
      "boundary_edges": ...,
      "nonmanifold_edges": ...,
      "self_intersections": ...,
      "min_intersection_volume": ...
    }
  },

  "surface": {
    "watertight": ...,
    "connected_components": ...,
    "nonmanifold_edges": ...,
    "boundary_edges": ...
  },

  "recommended_next_action": [...]
}

============================================================
十七、最终诊断报告
==================

运行后生成：

diagnostic_report.md

要求使用非常明确的结构。

---

1. Port step

---

回答：

视觉上的阶梯是否真实存在？

YES / NO

如果 YES：

真正的几何差值是多少？

r_branch =
r_extension =

area difference =

最可能原因：

...

证据：

...

---

2. Junction crack

---

回答：

截图中的白缝是否是真正 open crack？

YES / NO / INCONCLUSIVE

boundary edges =

non-manifold edges =

self intersections =

internal faces =

最可能原因：

...

---

3. Source SWC

---

是否有证据说明问题来自：

source centerline / radius？

YES / NO

不要因为 reconstruction 有问题就修改 SWC。

---

4. Reconstruction backend

---

问题主要位于：

resampling

tube construction

extension construction

junction primitive

Boolean union

normal/rendering

中的哪一阶段？

---

5. Explicit vs implicit

---

同一 ROI 上：

哪一种几何质量更高？

不要只根据肉眼判断。

提供指标。

---

6. Recommended fixes

---

只在诊断完成以后提出最小修复方案。

按优先级：

P0
P1
P2

不要在 root cause 未明确前大规模重构。

============================================================
十八、本轮禁止事项
==================

在最终诊断报告完成之前：

不要：

- 自动修改 SWC；
- 自动修改 radius；
- 自动把 sphere 放大 1.5 倍；
- 自动把 tube_sides 提到 128；
- 自动 smoothing junction；
- 自动切换所有 ROI 到 implicit；
- 删除有问题的 representative ROI；
- 继续 Navier-Stokes CFD；
- 继续 microbubble simulation。

当前任务是：

【找出原因。】

不是：

【让截图暂时看起来更漂亮。】

============================================================
十九、功能验收标准
==================

本轮任务只有满足以下条件才算完成：

1. 至少一个 port-step ROI 被完整诊断；
2. 至少一个 junction-artifact ROI 被完整诊断；
3. CUT_PORT endpoint position continuity 被定量检查；
4. CUT_PORT radius continuity 被定量检查；
5. branch mesh 和 extension mesh 的实际截面 radius 被测量；
6. extension overlap volume 被实际计算；
7. VTK tube radius scalar 语义通过 synthetic test 验证；
8. junction primitive overlap 被实际测量；
9. boundary edges 单独检测；
10. non-manifold edges 单独检测；
11. surface connected components 被检测；
12. normals/winding 被检测；
13. internal/self-intersection 至少进行有效检测；
14. triangle quality 被检查；
15. Boolean 前 primitives 被保存；
16. Boolean 后 mesh 被保存；
17. diagnostic overlay 图片完整输出；
18. synthetic straight tube test 通过；
19. synthetic Y-bifurcation test 完成；
20. 一个真实问题 ROI 完成 explicit vs implicit A/B comparison；
21. diagnostic_report.md 明确给出 root cause；
22. 所有结果保存到：
    ulm_3D_vascular\outputs
23. 不破坏现有 sampling 和 model_generate 功能。

============================================================
二十、最终向我汇报
==================

不要只回复：

“诊断完成”。

请明确告诉我：

A. Port 阶梯：

- 是不是真实 geometry discontinuity？
- endpoint position error 是多少？
- radius error 是多少？
- branch actual radius 是多少？
- extension actual radius 是多少？
- overlap volume 是多少？
- root cause 是什么？

B. Junction：

- 是不是真裂缝？
- boundary edge count？
- non-manifold edge count？
- self-intersection？
- internal face？
- Boolean overlap 是否充分？
- root cause 是什么？

C. Backend：

- explicit 是否是主要问题来源？
- implicit 是否消除了该问题？
- 两者 radius fidelity 和 runtime 分别是多少？

D. 下一步：

只给出基于证据的最小修复方案，
不要凭经验直接重写全部 pipeline。

最终目标是：

【用定量几何证据明确判断当前 port step 和 bifurcation defect 分别在哪一个处理阶段产生，以及是否影响 CFD-ready 条件。】
