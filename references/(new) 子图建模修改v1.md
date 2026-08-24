
你现在需要在现有项目中正式实现：

SWC ROI → 3D CFD lumen → CFD-ready watertight vascular surface

请不要只写方案或伪代码，而是先完整检查现有 repository，然后直接实现、运行、调试和验收。

============================================================
一、项目环境与硬性要求
======================

项目根目录中已有：

ulm_3D_vascular

Python 解释器固定使用：

d:\anaconda3\envs\pmp\python.exe

所有 Python 命令、pip 安装、测试和正式运行必须使用该解释器，例如：

d:\anaconda3\envs\pmp\python.exe -m pip ...
d:\anaconda3\envs\pmp\python.exe ...

可以安装任意合理、高效、成熟的 Python 功能包。

优先考虑：

numpy
scipy
pyvista
vtk
trimesh
manifold3d
shapely
meshio
matplotlib
pandas

如确有需要，可以增加：

scikit-image
numba
pymeshfix
gmsh

但不要为了堆依赖而安装无关包。

---

主入口必须实现于：

ulm_3D_vascular/model_generate.py

但是：

【严禁把所有代码堆进 model_generate.py】

model_generate.py 只负责：

1. CLI 参数解析；
2. 配置加载；
3. 调用各功能模块；
4. batch 调度；
5. logging；
6. 输出 run summary；
7. end-to-end workflow orchestration。

真正的几何计算、QC、surface reconstruction、port construction、
export 和 visualization 必须合理拆分至：

ulm_3D_vascular/utils/

下新建的功能目录。

建议新建：

ulm_3D_vascular/utils/cfd_lumen/

建议结构如下，可根据现有工程风格合理调整：

cfd_lumen/
    __init__.py
    config.py
    types.py
    roi_io.py
    geometry_preprocess.py
    collision_qc.py
    lumen_builder.py
    port_geometry.py
    surface_qc.py
    export.py
    visualization.py

不要为了“模块化”再拆成几十个只有一两个函数的小文件。

原则是：

一个文件解决一个明确的功能域。

============================================================
二、当前研究数据语义——必须严格遵守
====================================

当前项目已经采用 SWC-centric workflow。

正式 CFD geometry 的唯一几何/拓扑来源是：

人工修订后的 SWC ROI。

不要使用 segmentation Mask 重建 lumen。

不要：

SWC → Mask → surface

而应该：

SWC centerline + radius
        ↓
CFD lumen surface

当前 ROI 已经来自真实 analysis_swc，不是生成模型产生的。

ROI 中已经存在：

- local nodes
- local edges
- global node mapping
- global edge mapping
- exact CUT_PORT position
- CUT_PORT radius
- ROI bounding box
- branch topology

请先检查现有代码和 outputs 的真实数据结构，
复用已有 ROI extraction 和 branch graph，
不要重新实现另一套 ROI parser。

---

非常重要：

SWC parent → child 只是结构记录方向，

不是已经确认的 physiological flow direction。

因此 geometry generation 阶段只建立：

GEOMETRIC_OUTWARD_DIRECTION

而不能创建：

PHYSIOLOGICAL_FLOW_DIRECTION

也不要在本阶段自动把 CUT_PORT 分类成真实 inlet/outlet。

port 类型第一版可保存为：

CUT_PORT_UNASSIGNED

后续 1D flow / CFD boundary module 再确定：

inlet / outlet。

============================================================
三、本阶段总目标
================

对一个 representative ROI：

G_ROI

构建：

Ω_f

即三维 vascular lumen fluid domain。

最终 surface 应明确划分：

Γ_wall
Γ_port_0
Γ_port_1
...
Γ_port_N

整体流程：

Representative ROI
        ↓
ROI geometry validation
        ↓
Branch extraction
        ↓
Arc-length resampling
        ↓
Variable-radius tube construction
        ↓
Bifurcation/junction reconstruction
        ↓
CUT_PORT extension construction
        ↓
Robust solid union
        ↓
Watertight lumen
        ↓
Port patch identification
        ↓
Geometry QC
        ↓
CFD-ready VTP/STL
        ↓
optional volume-mesh verification

本阶段不求：

pressure
flow
velocity
microbubble trajectory

只负责构建可靠的 CFD geometry。

============================================================
四、Phase 0：先检查现有 Repository
==================================

编码前必须先阅读已有代码，确认：

1. representative ROI 从哪里保存；
2. selected ROI manifest 文件格式；
3. ROI nodes/edges 存储格式；
4. radius 数据单位；
5. physical coordinates 当前是否已经是 μm；
6. branch-level graph 是否已经保存；
7. CUT_PORT 数据字段；
8. CUT_PORT 是否保存：
   - global edge ID
   - intersection position
   - radius
   - boundary face
9. 当前 visualization 是否已有 PyVista/VTK helper；
10. 当前 logging/config/output run-id 机制。

优先复用现有：

- graph
- ROI
- config
- visualization
- logging

基础设施。

不得重写现有 Sampling Project。

不得改变 ROI sampling 结果。

============================================================
五、单位处理——这是硬性要求
============================

首先确认：

ROI 坐标当前采用：

μm

还是其他单位。

确认 SWC radius 当前实际单位。

不要根据猜测再次乘 voxel spacing。

必须从现有 pipeline 和真实数据检查确认。

程序需要明确维护：

geometry_unit = "um"

用于所有 reconstruction / visualization。

同时为 CFD 导出一套 SI 几何：

1 μm = 1e-6 m

因此建议输出：

lumen_surface_um.vtp

用于：

可视化
几何 QC
论文图片

以及：

lumen_surface_m.stl

用于：

CFD solver

STL 本身没有单位，
因此必须在同目录保存：

units.json

例如：

{
    "source_coordinate_unit": "um",
    "cfd_export_unit": "m",
    "scale_um_to_m": 1e-6
}

严禁产生一个不知道单位的 STL。

============================================================
六、Step 1：读取并标准化 ROI geometry
=====================================

将 ROI 表示为：

G_ROI = (V, E)

进一步整理为 branch collection：

B_i = {c_i(s), r_i(s)}

其中：

c_i(s) = [x(s), y(s), z(s)]

r_i(s) = local vessel radius

s ∈ [0, L_i]

必须保存：

branch_id
source_global_edges
source_global_nodes

不能因为 CFD reconstruction 丢失 source identity。

============================================================
七、Step 2：CFD geometry pre-QC
===============================

在创建任何三维 surface 前必须完成 geometry QC。

---

7.1 Radius validation
---------------------

要求所有实际参与 CFD reconstruction 的 radius：

r > 0
finite

遇到：

NaN
Inf
zero
negative radius

不要静默替换。

默认：

ROI geometry generation = FAIL

并记录：

roi_id
branch_id
node_id
radius
failure reason

因为：

R ∝ 1/r^4

radius 异常会严重污染之后 CFD。

---

7.2 Edge length validation
--------------------------

检查：

||p_{i+1} - p_i|| > ε

重复或近零长度 edge 应：

记录
拒绝/清理派生几何

但不得修改 reference SWC。

---

7.3 Topology preservation
-------------------------

记录输入 ROI 的：

node count
edge count
branch count
bifurcation count
connected component count

geometry generation 后不能修改这些结构关系。

禁止自动：

增加血管 branch
连接两个非相邻 branch
删除真实 branch

============================================================
八、Step 3：非邻接血管碰撞检测
==============================

这是非常重要的 CFD geometry QC。

两个 topology 上不相邻的 branches：

B_i
B_j

可能在空间上距离很近。

若：

d_min(B_i, B_j)
<
r_i + r_j

则三维管道扩张后会互相穿透，

导致本来不连接的 SWC 被错误连接。

---

实现方式要求高效：

第一阶段使用：

scipy.spatial.cKDTree

对 resampled centerline 做快速 candidate pair search。

只对空间邻近 pair 做精确 segment-to-segment distance calculation。

不要：

每个 segment 与全部 segment O(N²) 暴力比较。

---

排除：

真正 topology adjacent branches

以及共享真实 bifurcation node 的 branch pair。

对于非邻接 branch，计算：

clearance = d_min - (r_i + r_j)

定义：

clearance < -collision_tolerance

为：

HARD_COLLISION

定义接近 0 的情况：

NEAR_CONTACT

配置参数全部进入 config：

collision_tolerance_um
near_contact_tolerance_um

不要硬编码在算法内部。

---

默认策略：

HARD_COLLISION
→ 当前 ROI CFD geometry generation FAIL

不要自动缩小 radius。

不要自动改变 topology。

输出：

collision_report.csv

并生成碰撞位置 3D QC 图。

============================================================
九、Step 4：Branch arc-length resampling
========================================

对每条 branch：

p_0 ... p_n

建立：

s_0 = 0

s_k =
Σ ||p_{j+1} - p_j||

然后沿 s 重采样。

目标：

获得更加均匀的中心线点间距，
便于 variable-radius tube reconstruction。

---

默认不要做 aggressive smoothing。

原因：

当前 SWC 来自人工修订数据，
原数据论文已经包含 skeleton smoothing / refinement。

因此：

centerline_smoothing = false

作为正式默认值。

---

resampling spacing 不使用一个完全固定的绝对值。

支持：

Δs = α r_local

作为目标 spacing。

建议默认配置：

resample_radius_fraction = 0.35

并允许：

min_resample_spacing_um
max_resample_spacing_um

约束极端情况。

具体值写进 config，不得散落在代码中。

---

centerline interpolation：

优先：

linear / shape-preserving interpolation

避免 spline overshoot 改变血管路径。

radius interpolation：

优先使用：

scipy.interpolate.PchipInterpolator

因为 radius 必须保持正值并减少 overshoot。

重采样后重新检查：

r(s) > 0。

============================================================
十、Step 5：默认高速 lumen reconstruction backend
=================================================

正式默认 backend：

MANIFOLD_EXPLICIT

核心思路：

每个 branch
→ variable-radius tube

每个 bifurcation
→ junction solid

每个 CUT_PORT
→ straight extension

然后：

robust manifold Boolean union

最终形成：

一个 closed vascular lumen。

---

10.1 Variable-radius tube
-------------------------

优先使用 VTK C++ backend：

vtkTubeFilter

不要用 Python for-loop 手工逐 triangle 建管道。

给 centerline points 设置：

radius scalar

使用：

VaryRadiusByAbsoluteScalar

构建 variable-radius tube。

每个 tube 输出必须：

triangulated
capped
manifold

tube 圆周离散数量参数：

tube_sides

进入 config。

推荐默认：

tube_sides = 32

允许用户修改。

---

10.2 Junction reconstruction
----------------------------

不能简单：

tube A
+
tube B
+
tube C

只 concatenate triangles。

这会产生：

holes
non-manifold junction
overlap

对于真实 SWC bifurcation node：

v_j

在该 node 建立 junction solid。

默认采用：

sphere centered at junction

radius：

直接使用该 SWC junction node 的真实 radius

不要人为放大。

即：

r_junction = r_swc(node)

该 sphere 仅作为：

相邻真实 branches 的几何融合体。

---

只允许：

SWC topology 中真实相邻 branch

在该 junction 中融合。

不能利用空间靠近自动创建新的连接。

============================================================
十一、Step 6：CUT_PORT 几何方向
===============================

对于每个 CUT_PORT：

已知：

x_k
r_k
global_edge_id

计算局部中心线 tangent：

t_k

只表示：

GEOMETRIC_TANGENT

不是 flow direction。

---

确定 outward tangent：

n_out

方法：

从候选方向：

+t_k
-t_k

分别沿少量 ε 前进，

判断哪一个方向从 ROI bounding box 内部走向外部。

选择那个方向：

n_out

保存：

outward_tx
outward_ty
outward_tz

这个方向用于：

port extension

而不是表示血流方向。

============================================================
十二、Step 7：CUT_PORT straight extension
=========================================

直接在原 ROI box face 截断 CFD geometry 不够理想。

对每个 CUT_PORT 建立直管 extension。

定义：

D_k = 2 r_k

extension length：

L_ext = λ D_k

建议默认：

port_extension_diameters = 5.0

即：

L_ext = 5D

但必须写入 config，
并在后续允许 3D / 5D / 10D sensitivity study。

---

extension radius：

r_ext = r_k

第一版保持 constant radius。

---

为了保证 extension 与原 tube 有体积重叠，
extension 不应只从 CUT_PORT 平面向外开始。

向 ROI 内部额外延伸：

L_overlap = η D

默认：

port_overlap_diameters = 0.5

因此 cylinder 实际范围：

x_cut - L_overlap n_out

到

x_cut + L_ext n_out

最外端：

x_cap =
x_cut + L_ext n_out

作为最终 CFD boundary cap center。

---

extension cylinder 必须具有：

flat outer cap

这样后续获得真正平面的 CFD boundary。

============================================================
十三、Step 8：Robust Boolean Union
==================================

不要使用普通 fragile triangle boolean 作为核心。

安装并优先使用：

manifold3d

或者：

trimesh + manifold3d backend

将以下 solids：

all branch tubes
all junction spheres
all port extension cylinders

进行 Boolean union。

---

为了高速和稳定性：

不要从第一个 mesh 开始依次 union 100 次形成不断膨胀的 mesh。

实现 balanced union tree：

例如：

level 1:
A∪B
C∪D
E∪F

level 2:
AB∪CD

...

减少 Boolean accumulated complexity。

---

Boolean 结束后：

clean duplicate vertices
triangulate
fix normals orientation

但：

禁止用自动 repair 改变 topology。

============================================================
十四、Implicit backend 作为 fallback
====================================

不要第一版就默认使用高分辨率全 ROI voxel field，
因为可能占用大量内存和计算时间。

但是请保留一个干净的 backend interface：

backend:
    manifold
    implicit

默认：

manifold

当：

Boolean exception

或者：

surface QC FAIL

时，

如果：

allow_implicit_fallback = true

再调用 implicit backend。

---

implicit backend 可以采用：

resampled centerline + radius

建立：

φ(x) =
min_i ( ||x - c_i|| - r_i )

并通过：

scipy cKDTree

查找近邻 centerline samples，

分块计算 field，

禁止一次性创建不受控制的巨大 float64 体数据。

要求：

float32
chunked evaluation

最后使用：

skimage.measure.marching_cubes

提取：

φ = 0

surface。

---

implicit backend 同样必须服从：

non-adjacent collision QC

不能因为 field overlap 自动改变 SWC topology。

============================================================
十五、Step 9：Port patch identification
=======================================

Boolean union 后必须知道：

哪些 triangles 是 wall

哪些 triangles 是 port。

对每个 port k 已知：

cap center:

c_k

cap normal:

n_k

radius:

r_k

对于 surface triangle centroid：

x_f

计算：

plane distance：

d_plane =
|(x_f - c_k) · n_k|

以及 radial distance：

d_radial =
||(x_f-c_k)
-----------

[(x_f-c_k)·n_k]n_k||

满足：

d_plane < plane_tolerance

且：

d_radial <= r_k + radial_tolerance

的 faces 标记为：

PORT_k

其余：

WALL

---

同时检查：

每个 port patch：

1. 必须存在；
2. 必须是一个 connected surface component；
3. surface normal 应与 n_out 基本一致；
4. cap area 应接近：

A_expected = π r_k²

计算：

area_error =
|A_mesh - πr²| / (πr²)

写入 QC。

============================================================
十六、Step 10：Boundary patch metadata
======================================

每个 port 保存：

port_id
roi_id
cut_port_id
global_edge_id

original_x_um
original_y_um
original_z_um

cap_x_um
cap_y_um
cap_z_um

radius_um
diameter_um

outward_tx
outward_ty
outward_tz

extension_length_um

patch_area_um2

boundary_role

其中第一版：

boundary_role =
CUT_PORT_UNASSIGNED

不要在 geometry generation 里自动写：

INLET
OUTLET

除非现有数据中已经有一个后续显式传入的 simulation BC file。

预留：

P_1D
Q_1D

字段，

但本阶段默认：

null。

============================================================
十七、Step 11：Surface QC
=========================

CFD surface 必须经过严格自动验收。

至少实现：

---

17.1 Connectedness
------------------

最终 lumen 必须：

surface connected components = 1

如果 >1：

FAIL

---

17.2 Watertight
---------------

cap 完成后：

boundary edges = 0

使用：

trimesh.is_watertight

以及 VTK/PyVista feature-edge check

交叉确认。

---

17.3 Manifold
-------------

检查：

non-manifold edges = 0

---

17.4 Positive volume
--------------------

enclosed volume > 0

---

17.5 Winding / normals
----------------------

surface normals orientation consistent

---

17.6 Port number
----------------

expected port count

必须等于：

detected port patch count

---

17.7 Source-topology preservation
---------------------------------

输出：

source branch count
source bifurcation count
source CUT_PORT count

geometry builder 不允许生成额外 vessel connectivity。

============================================================
十八、Step 12：Radius fidelity QC
=================================

这是整个模块最重要的验收之一。

因为：

R ∝ r^-4

几何 radius 的小误差会显著影响 CFD。

---

在每条 branch 中选择若干：

远离：

junction
port

的 centerline sample。

建议默认跳过距 junction / port：

2D

以内的区域。

---

在 sample point：

c(s)

使用 tangent：

t(s)

建立垂直截面：

normal = t(s)

与 reconstructed lumen surface 求交。

得到 cross-section polygon area：

A_CFD(s)

定义 equivalent radius：

r_eq(s) =
sqrt(A_CFD(s) / π)

与 source SWC radius：

r_SWC(s)

比较：

ε_r(s) =
[r_eq(s) - r_SWC(s)] / r_SWC(s)

输出：

median absolute relative error
mean absolute relative error
P95 absolute relative error
max absolute relative error

以及每个 branch 的结果。

---

不要一开始拍脑袋设置一个“论文通过阈值”。

先通过真实 ROI 测试得到误差分布。

config 中允许用户设置：

max_radius_p95_error

若为 null：

只报告，不作为 FAIL。

之后再根据 convergence experiment 确定正式 threshold。

============================================================
十九、Step 13：几何离散 convergence
===================================

至少实现一个可配置的 geometry resolution test。

对于同一个简单 ROI，

测试例如：

tube_sides:
16
24
32
48

比较：

surface area
enclosed volume
port area
radius fidelity

默认正式运行：

32

但必须输出 convergence 结果供确认。

不要假定 32 一定最优。

============================================================
二十、Step 14：CFD SI export
============================

内部 reconstruction 使用：

μm

输出 CFD surface 时统一：

x_m = x_um × 1e-6

生成：

lumen_surface_m.stl

同时：

wall_m.stl

以及：

ports/
    port_000_m.stl
    port_001_m.stl
    ...

推荐同时输出：

lumen_surface_m.vtp

并在 VTP face data 中保存：

patch_id
patch_type
port_id

---

因为 STL 无法可靠保存 boundary metadata，

VTP 是 master surface format，

STL 是 solver compatibility format。

============================================================
二十一、Optional：Volume mesh verification
==========================================

本阶段主要目标是 surface reconstruction。

但是建议增加一个：

--build-volume-mesh

可选功能，

用于确认 lumen 确实可以 tetrahedralize。

优先考虑：

gmsh

或成熟 tetrahedral mesher。

该功能默认：

false

当启用时：

watertight surface
    ↓
volume mesh
    ↓
mesh QC

至少报告：

number of nodes
number of cells
minimum cell quality
negative-volume cell count
boundary patch recovery status

但：

不要在本阶段运行 Navier-Stokes CFD。

============================================================
二十二、性能要求
================

代码要尽可能高速。

要求：

1. NumPy vectorization；
2. SciPy cKDTree；
3. VTK C++ filters；
4. manifold3d C++ Boolean；
5. 避免 Python triangle-level loops；
6. ROI batch 可并行；
7. 大数组使用 float32，在无需 float64 时不要浪费内存。

---

Windows 下：

VTK objects 不建议多线程共享。

批量 ROI 使用：

concurrent.futures.ProcessPoolExecutor

做 ROI-level parallel。

model_generate.py 必须使用：

if __name__ == "__main__":

避免 Windows multiprocessing spawn 问题。

---

workers 参数：

--workers N

默认：

min(
    max(cpu_count - 1, 1),
    number_of_rois
)

但如果 manifold/VTK 在多进程中表现不稳定，
优先保证正确性并允许：

--workers 1

============================================================
二十三、model_generate.py CLI
=============================

至少支持：

--config PATH

--sampling-run PATH

--roi-id ROI_ID

--selected-only

--all-selected

--backend manifold|implicit

--workers N

--headless

--build-volume-mesh

--run-id

--output-dir

---

如果用户没有指定 roi-id，

默认：

处理当前 sampling run 中的 selected representative ROIs。

同时支持先：

--roi-id

单独跑一个最简单 ROI。

============================================================
二十四、建议配置文件
====================

建立：

cfd_lumen_config.yaml

示意：

geometry:
  source_unit: um
  export_unit: m

  centerline_smoothing: false

  resample_radius_fraction: 0.35
  min_resample_spacing_um: 0.2
  max_resample_spacing_um: 1.0

  tube_sides: 32

junction:
  enabled: true
  radius_scale: 1.0

ports:
  extension_diameters: 5.0
  overlap_diameters: 0.5

collision_qc:
  enabled: true
  hard_collision_tolerance_um: 0.0
  near_contact_tolerance_um: 0.2

boolean:
  backend: manifold
  allow_implicit_fallback: true

implicit_fallback:
  dtype: float32
  cells_across_min_diameter: 8
  chunk_size: 1000000
  k_nearest: 16

surface_qc:
  check_watertight: true
  check_manifold: true
  check_connected: true
  check_port_area: true
  radius_fidelity_samples_per_branch: 10
  max_radius_p95_error: null

volume_mesh:
  enabled: false

所有参数应根据实际工程合理调整。

不要照抄任何明显不适合当前数据的参数。

============================================================
二十五、Outputs 目录要求
========================

所有日志、可视化和功能验收结果必须保存至：

ulm_3D_vascular/outputs/

建议：

ulm_3D_vascular/
└── outputs/
    └── model_generate/
        └── <run_id>/
            │
            ├── logs/
            │   └── model_generate.log
            │
            ├── config/
            │   └── resolved_config.yaml
            │
            ├── manifests/
            │   ├── geometry_summary.csv
            │   ├── failed_rois.csv
            │   └── run_manifest.json
            │
            ├── figures/
            │   ├── all_roi_geometry_summary.png
            │   └── ...
            │
            ├── report/
            │   └── run_summary.json
            │
            └── rois/
                └── <roi_id>/
                    │
                    ├── source/
                    │   ├── roi_metadata.json
                    │   └── branch_mapping.csv
                    │
                    ├── centerline/
                    │   ├── raw_centerline.vtp
                    │   └── resampled_centerline.vtp
                    │
                    ├── geometry/
                    │   ├── lumen_surface_um.vtp
                    │   ├── lumen_surface_m.vtp
                    │   ├── lumen_surface_m.stl
                    │   ├── wall_m.stl
                    │   └── ports/
                    │       ├── port_000_m.stl
                    │       └── ...
                    │
                    ├── boundary/
                    │   ├── ports.csv
                    │   └── ports.json
                    │
                    ├── qc/
                    │   ├── pre_geometry_qc.json
                    │   ├── collision_report.csv
                    │   ├── surface_qc.json
                    │   ├── radius_fidelity.csv
                    │   └── units.json
                    │
                    ├── figures/
                    │   ├── 01_source_centerline.png
                    │   ├── 02_resampled_centerline.png
                    │   ├── 03_lumen_overlay.png
                    │   ├── 04_ports_and_normals.png
                    │   ├── 05_surface_wireframe.png
                    │   ├── 06_radius_fidelity.png
                    │   └── 07_cross_sections.png
                    │
                    └── mesh/
                        └── [optional]

每次运行建立独立 run_id。

禁止覆盖历史运行。

============================================================
二十六、Logging 要求
====================

不要大量使用 print()。

统一 logging。

console + file 同时输出。

至少记录：

Input sampling run
ROI IDs
Backend

Per ROI:

node count
edge count
branch count
bifurcation count
CUT_PORT count

radius min/median/max
centerline total length

resampled point count

collision candidate count
hard collision count
near-contact count

tube primitive count
junction primitive count
extension count

Boolean runtime
surface triangle count
surface vertex count

watertight status
manifold status
surface component count

enclosed volume
surface area

port count
port areas

radius fidelity metrics

export paths

runtime:
    load
    resample
    collision QC
    primitive construction
    Boolean
    patch detection
    surface QC
    visualization
    export
    total

failure reason

============================================================
二十七、必须生成的验收可视化
============================

所有图必须基于真实运行结果，不得画示意图冒充结果。

---

Figure 1：SWC source geometry
-----------------------------

显示：

raw centerline
radius variation
CUT_PORT

---

Figure 2：Raw vs Resampled Centerline
-------------------------------------

叠加：

raw centerline
resampled centerline

用于确认 resampling 没有改变路径。

---

Figure 3：Centerline + CFD lumen overlay
----------------------------------------

显示：

semi-transparent lumen
source centerline
bifurcations
CUT_PORT extensions

这是最重要的 geometry sanity check。

---

Figure 4：Ports and normals
---------------------------

显示：

每个 flat CFD cap
port ID
outward normal arrow
original CUT_PORT position

用于确认：

port plane
extension
direction

均正确。

---

Figure 5：Surface wireframe
---------------------------

观察：

junction
triangle quality
port geometry
surface discontinuity

---

Figure 6：Radius fidelity
-------------------------

横轴：

source radius

纵轴：

reconstructed equivalent radius

同时显示：

y = x

以及误差统计。

---

Figure 7：Cross-section examples
--------------------------------

选若干真实 branch locations，

显示：

CFD surface cross-section
centerline point
target circle radius

确认 lumen 没有系统性偏粗/偏细。

============================================================
二十八、Unit tests
==================

不要只依赖真实 ROI。

建立 synthetic geometry tests。

至少实现：

---

Test 1：Straight constant-radius vessel
---------------------------------------

输入：

straight centerline
constant radius

检查：

watertight
one component
two ports
radius fidelity
port area ≈ πr²

---

Test 2：Tapered vessel
----------------------

检查：

radius 沿程变化是否被保留。

---

Test 3：Curved vessel
---------------------

检查：

tube orientation
surface continuity
无错误 twisting。

---

Test 4：Y bifurcation
---------------------

检查：

three branches
one junction
one connected lumen
watertight
no non-manifold edges。

---

Test 5：Non-adjacent collision
------------------------------

构造两条 topology 不连接但几何相交的血管。

要求：

collision QC

在 Boolean 前检测并拒绝。

---

Test 6：Multiple CUT_PORTs
--------------------------

检查：

extension 数量
port patch 数量
port IDs
outward normals。

---

Test 7：Unit conversion
-----------------------

检查：

μm → m

精确为：

1e-6

---

Test 8：Determinism
-------------------

相同 input + config

输出：

surface metrics
patch IDs
QC metrics

保持一致。

============================================================
二十九、真实 ROI smoke test
===========================

完成 synthetic tests 后，

必须在当前真实 representative ROI 上跑一次 end-to-end。

优先选择：

branch 数较少、
CUT_PORT 明确、
不存在 geometry collision

的 representative ROI。

当前研究记录中已有约 5 条聚合 branch 的真实代表案例；
若当前 outputs 中仍能定位到该 ROI，
优先作为第一批 smoke test。

但不要硬编码 ROI ID。

程序应根据实际 manifest 查询。

---

第一次真实运行：

只处理 1 个 ROI：

workers = 1

完成：

load
↓
QC
↓
resample
↓
tube construction
↓
junction
↓
port extension
↓
Boolean
↓
patch labeling
↓
surface QC
↓
visualization
↓
export

确认全部 PASS 后，

再运行全部 selected representative ROIs。

============================================================
三十、Acceptance Criteria / Definition of Done
==============================================

只有满足以下条件才算本任务完成：

1. model_generate.py 可以独立执行；
2. 核心算法没有堆积在 model_generate.py；
3. utils/cfd_lumen/ 模块划分合理；
4. 能从现有 representative ROI 输出直接读取真实数据；
5. 不依赖 segmentation Mask；
6. 不修改原 SWC topology；
7. branch centerline + radius 成功重建为 3D lumen；
8. bifurcation 处形成连续 lumen；
9. non-adjacent collision 能提前检测；
10. CUT_PORT 成功生成 straight extensions；
11. CFD boundary caps 基本垂直于局部 centerline；
12. 每个 port 有唯一 patch ID；
13. VTP 可以保存 wall/port labels；
14. STL 输出为 SI meter coordinates；
15. 输出明确 units metadata；
16. 最终 surface：
    connected = true
    watertight = true
    manifold = true
17. port count 与 source CUT_PORT count 一致；
18. radius fidelity 被定量测量；
19. 所有 synthetic unit tests 通过；
20. 至少一个真实 representative ROI 完成 end-to-end PASS；
21. 必要日志全部保存至：
    ulm_3D_vascular/outputs
22. 必要可视化全部保存至：
    ulm_3D_vascular/outputs
23. 失败 ROI 不能静默跳过；
    必须保存 failure reason；
24. 固定 input/config 结果可复现；
25. 不破坏已有：
    SWC preprocessing
    ROI sampling
    visualization
    clustering
    representative selection

功能。

============================================================
三十一、实现过程中的原则
========================

1. 不要为了让 geometry PASS 而偷偷修改 source SWC。
2. 不要自动缩小 radius 来解决 collision。
3. 不要自动连接空间上靠近但 topology 不相邻的 branches。
4. 不要重新引入 segmentation Mask。
5. 不要把 CUT_PORT outward tangent 叫做 blood-flow direction。
6. 不要在 geometry generation 阶段猜 inlet/outlet。
7. 不要为了表面好看进行 aggressive smoothing。
8. 不要只检查 STL“看起来正常”；
   必须做 quantitative QC。
9. 不要输出无法追溯 source branch/global edge 的 geometry。
10. 所有重要阈值必须配置化。

============================================================
三十二、本任务暂不实现的功能
============================

本轮不要加入：

1D blood flow solver
physiological pressure
inlet/outlet physiological classification
CFD Navier-Stokes
microbubble transport
PROTEUS acoustic simulation

但 ports metadata 中为未来预留：

P_1D
Q_1D
boundary_role

字段。

============================================================
三十三、最终需要向我汇报
========================

代码执行完成后，请不要只说“完成”。

必须给出：

1. 实际新增/修改的文件列表；
2. 每个文件的职责；
3. 安装了哪些新 Python package；
4. 使用的实际 Python interpreter；
5. synthetic tests 结果；
6. 真实 ROI smoke-test：

   - ROI ID
   - branches
   - ports
   - triangles
   - surface volume
   - watertight
   - manifold
   - radius fidelity
   - runtime
7. 输出目录；
8. 关键验收图片路径；
9. 是否存在失败项；
10. 当前还不能保证什么。

如果运行中发现现有 ROI 文件结构和以上假设不同：

先检查现有代码真实数据流，
然后适配现有数据结构。

不要新造一套平行的数据格式。

最终目标是：

【直接使用当前 Sampling Project 的 representative ROI，
稳定生成可追溯、拓扑不变、管径误差可量化、
具有明确 CFD port patches 的三维 watertight vascular lumen。】
