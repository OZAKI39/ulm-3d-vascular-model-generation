
当前 v6 已经实现了 topology-valid hybrid reconstruction，
但实际可视化仍在 local-implicit / explicit-branch 最终 merge ring
出现明显 triangular saw-tooth artifacts。

这些 artifacts 已经改变 surface silhouette，
因此不是单纯 shading 问题。

本轮不要继续修改 hybrid transition。

新增 experimental backend：

unified_polyball

目标：

使用整个 CFD_DOMAIN 的 SWC centerline + radius
生成一个单一连续 implicit vascular surface，
彻底取消：

explicit branch surface
local implicit surface
surface stitching
surface Boolean merge

之间的任何 surface interface。

---

1. 优先调查并尝试安装 VMTK

---

Python:

d:\anaconda3\envs\pmp

优先通过 conda-forge / compatible package 安装 VMTK。

如果当前 Windows/Python 环境无法可靠安装，
不要破坏环境。

保留 custom cKDTree/polyball backend fallback。

---

2. 使用 vtkvmtkPolyBallLine / PolyBallModeller

---

将完整 CFD_DOMAIN centerline 转成 vtkPolyData lines。

point/cell 中保存：

Radius

要求 radius 沿 centerline line 连续插值。

生成一个：

single polyball implicit function

覆盖全部：

branches
junctions
context extensions
port extensions。

---

3. 禁止 surface-level hybrid

---

unified_polyball backend 不允许：

vtkTubeFilter branch surface
junction surface stitching
branch/junction Boolean
collar-loop stitching

最终 vascular wall 只能来自：

ONE implicit field
→ ONE iso-surface extraction。

---

4. Port

---

port extension 继续通过 centerline extension 产生。

但 centerline 应额外超过最终 CFD boundary plane：

port_tail_length = [configurable, e.g. 1D]

先让 unified polyball 形成连续 tube。

之后在：

exact CFD port plane

统一 clip surface。

再生成 planar cap。

因此 port 本身不能成为 implicit termination sphere。

---

5. Non-adjacent collision QC

---

继续使用现有 collision QC。

如果 topology-disconnected branches 的 geometric envelopes overlap：

unified_polyball = FAIL

不得因为 implicit union 自动连接。

---

6. Field extraction

---

比较：

vtkFlyingEdges3D

与当前 marching_cubes。

优先选择：

更快且几何稳定的方法。

测试：

cells_across_min_diameter:
12
16
20
24

---

7. Zero-level vertex projection

---

iso-surface extraction 后，

对每一个 wall vertex：

使用 exact polyball function：

phi(x)

以及：

gradient phi(x)

执行 1–3 次 Newton projection：

x_new =
x -
phi(x)/||grad phi(x)||² * grad phi(x)

直到：

|phi| < tolerance

或达到 max iteration。

记录：

pre_projection_radius_error
post_projection_radius_error。

---

8. Remeshing

---

对投影后的完整 unified wall surface：

使用 VMTK surface remeshing
或功能等价的高质量 isotropic remesher。

目标：

消除 marching-cubes triangle-size variation，
不是改变 vascular geometry。

禁止自由 smoothing。

---

9. Remesh 后再次 project

---

remeshing 后的 vertices：

再次投影到：

phi = 0。

因此：

geometry source remains SWC + radius。

---

10. Port faces

---

完成 wall surface 后：

clip at final CFD planes

并创建：

flat port caps。

分别赋予：

WALL
PORT_0
PORT_1
...

boundary IDs。

---

11. v6 vs v7 comparison

---

同一个 ROI anchor_003274，

输出完全相同 camera：

v6 flat
v7 flat

v6 smooth
v7 smooth

wireframe comparison

silhouette comparison。

---

12. 重点 QC

---

v7 必须继续：

self intersection = 0
internal face = 0
internal cap = 0
boundary edge = 0
nonmanifold edge = 0
degenerate triangle = 0
components = 1

同时比较：

radius P95
collar radius error
port area error

triangle count
runtime

normal jump P95/P99/max

triangle aspect ratio
silhouette roughness。

---

13. 新增 seam check

---

因为 v7 不存在 hybrid interface，

最终 geometry 中：

hybrid_interface_edge_count

必须：

0

并确认当前图片中三角锯齿所在的位置：

不存在任何 reconstruction interface。

---

14. 评价原则

---

不要因为 v7 runtime 比 v6 慢就直接否决。

当前代表 ROI 数量有限，
正确、连续的 CFD geometry 优先级高于单个 ROI 快十几秒。

同时也不要因为 visually smooth 就接受。

必须保持：

radius fidelity
topology
mesh validity。

---

15. 正式选择

---

只有真实 A/B 结果完成后才决定：

A. KEEP_V6_HYBRID

或

B. ADOPT_V7_UNIFIED_POLYBALL

如果 v7：

消除 visible saw-tooth seam，

且：

radius P95 <= [根据真实结果比较，不预设]

拓扑 QC 全部通过，

runtime 在可接受量级，

则推荐：

ADOPT_V7_UNIFIED_POLYBALL。
