
当前任务是：

【判断当前 fMOST mouse brain vascular dataset 中实际提供的 Mask，
是否具有足够的独立几何信息，可以作为 SWC → CFD STL reconstruction
的辅助参考，尤其用于改善 bifurcation / junction 区域。】

本轮不是立即修改正式 STL reconstruction pipeline。

必须首先完成：

Mask provenance identification
        ↓
SWC–Mask relationship analysis
        ↓
local junction quantitative comparison
        ↓
A/B/C reconstruction experiment
        ↓
evidence-based recommendation

只有有足够证据以后，才允许建议是否把 Mask 正式加入
CFD lumen reconstruction。

============================================================
一、项目环境
============

项目：

ulm_3D_vascular

Python interpreter：

d:\anaconda3\envs\pmp\python.exe

所有 Python 命令必须使用该解释器。

例如：

d:\anaconda3\envs\pmp\python.exe ...

可以安装必要的高效 Python package。

优先复用已有：

numpy
scipy
pandas
vtk
pyvista
trimesh
scikit-image
matplotlib

不要为了该实验修改当前正式 SWC sampling pipeline。

============================================================
二、论文中已经明确的事实
========================

请先在 repository 中定位并阅读论文：

A high-resolution dataset of mouse brain vasculature
for deep learning-based reconstruction

特别检查：

Materials and Methods
Human-in-the-loop data annotation
Reliability analysis

必须把以下论文事实作为分析前提：

1. preliminary semi-automatic annotation 阶段：

先获得：

vascular skeleton
+
vessel radius

然后：

binary vessel masks

是根据这些 annotation 通过计算生成的。

也就是说至少存在一种：

SWC / skeleton + radius
→ computational Mask

而不是逐 voxel 人工勾画的 Mask。

---

2. 后续存在：

3D U-Net
→ predicted vessel mask
→ automatic skeleton reconstruction
→ manual topology correction

因此论文实际上涉及：

至少两种来源不同的 Mask：

A:
semi-automatic annotation-derived Mask

B:
U-Net predicted Mask

不能把二者混为一谈。

---

3. topology correction 后：

skeleton endpoints
和
branch points

被固定，

skeleton 经过 smoothing / resampling，

然后根据 local vessel radius 使用 Gaussian kernel
形成 vascular morphology annotation map。

因此：

semi-automatic Mask / morphology annotation

本身很可能高度依赖：

corrected skeleton + radius。

---

4. 论文的人工参考 skeleton / topology
   比 fully automatic segmentation-derived topology 更可靠。

因此即使最终决定使用 Mask：

Mask 不允许：

改变 SWC topology
增加 branch
删除 branch
连接 disconnected branches。

============================================================
三、本任务首先要回答的问题
==========================

不要先假定：

Mask 有帮助。

也不要先假定：

Mask 没有帮助。

需要依次回答：

Q1.
当前项目中实际保存的 Mask 是哪一种？

Q2.
这个 Mask 是：

SWC-derived
U-Net predicted
还是无法确定？

Q3.
如果是 SWC-derived Mask，
它相对于当前 SWC + radius 是否真的包含额外几何信息？

Q4.
如果是 U-Net predicted Mask，
它与最终人工修订 SWC 在：

centerline
radius
junction volume
topology

方面有多大差异？

Q5.
Mask 是否能改善当前 STL 中仍存在的：

junction surface transition
normal discontinuity
局部视觉突兀

同时又不降低：

radius fidelity
topological validity
surface validity？

============================================================
四、Phase 1：检查当前 Dataset 文件来源
======================================

先扫描当前 dataset directory。

不要修改数据。

列出每一个 sample 对应的：

raw image
Mask
SWC

完整路径。

对 Mask 文件记录：

filename
parent folder
data split
shape
dtype
value range
binary / probability / label
voxel spacing
corresponding SWC
corresponding image

输出：

mask_inventory.csv

============================================================
五、从文件结构判断 Mask provenance
==================================

检查：

文件名
folder name
metadata
README
dataset manifest
download structure
annotation script
generation script

搜索 repository 中涉及：

mask
segmentation
semi-automatic
annotation
ground truth
prediction
U-Net
probability map
Gaussian

的代码和说明。

如果存在生成 Mask 的脚本，

确定其输入究竟是：

SWC + radius

还是：

neural-network output。

============================================================
六、建立 Mask provenance 分类
=============================

每个 Mask 必须被分类为：

SWC_DERIVED_ANNOTATION_MASK

UNET_PREDICTED_MASK

UNKNOWN_MASK_PROVENANCE

不要根据视觉自行判断。

每项结论必须给：

evidence_source
file/path
code/document evidence
confidence

输出：

mask_provenance_report.csv

============================================================
七、检查当前实际使用的 sample
=============================

重点定位当前用于：

representative ROI
model_generate
CFD lumen

的真实 sample。

确定它所对应的 Mask provenance。

最终必须明确回答：

【当前 ROI003274 或对应实际 ROI 的 Mask，
究竟属于哪一种 Mask。】

不要假设。

============================================================
八、Phase 2：SWC 与 Mask 的信息依赖性测试
=========================================

如果当前 Mask 被判断为：

SWC_DERIVED_ANNOTATION_MASK

则必须验证：

当前 Mask 是否基本可以由：

current SWC centerline + current radius

重新生成。

---

建立一个：

SWC_RECONSTRUCTED_MASK

使用当前：

SWC centerline
+
radius

进行 voxelization。

注意：

该步骤只用于 provenance / information analysis。

不是正式 STL pipeline。

---

使用与论文尽可能一致的逻辑：

centerline resampling
local-radius-aware Gaussian / tubular expansion

若无法完全复现论文算法，

明确说明 approximated implementation。

============================================================
九、比较 Original Mask 和 SWC-reconstructed Mask
================================================

统一到：

相同 shape
相同 xyz/zyx orientation
相同 physical spacing：

1 × 1 × 2 µm

然后计算：

Dice

IoU

precision

recall

volume difference

surface distance

Hausdorff distance

95th-percentile Hausdorff distance

============================================================
十、解释结果
============

如果：

Original Mask

和：

SWC-reconstructed Mask

高度一致，

则说明：

Mask largely contains information already encoded by
SWC + radius

因此它不能被称为：

independent lumen ground truth。

---

如果差异显著，

则进一步定位差异区域：

branch
junction
terminal
large vessel
capillary

分别统计。

不能只给一个 global Dice。

============================================================
十一、Phase 3：重点分析 Junction
================================

当前 STL 主要残余问题集中在：

junction surface smoothness / transition。

因此重点分析现有真实问题 junction：

优先：

J49

同时选择至少：

一个正常 junction

作为 control。

============================================================
十二、提取 Local Junction Mask
==============================

围绕 junction：

x_J

建立局部 physical bounding box。

范围必须覆盖：

junction core
+
incident branches
+
现有 collar region

不要截取整个原始 volume。

保存：

junction_<id></id>_mask.tif
junction_<id></id>_swc.vtp
junction_<id></id>_current_surface.vtp

============================================================
十三、Mask 不能直接当作 topology
================================

即使 binary Mask 在多个 branch 间连通，

仍然必须使用：

current corrected SWC

定义：

incident branches
junction topology。

Mask 只能用于：

volumetric-shape comparison。

============================================================
十四、Phase 4：建立三个 reconstruction 对照
===========================================

必须在完全相同 junction 上建立：

A. SWC_ONLY

当前 v4/v5 正式方法：

corrected SWC
+
radius
+
local implicit / explicit hybrid

---

B. MASK_ONLY

仅用于实验比较：

local Mask
→ signed distance / marching cubes
→ local surface

注意：

B 不是正式 topology 方法。

只用于回答：

“Mask 中的 volumetric junction shape 长什么样？”

---

C. SWC_MASK_ASSISTED

使用：

SWC topology + radius

作为硬约束，

Mask 只作为：

local surface-shape soft reference。

============================================================
十五、MASK_ONLY reconstruction
==============================

不要直接对 binary Mask 用最低质量 Marching Cubes 后肉眼比较。

先处理 physical spacing：

dx = 1 µm
dy = 1 µm
dz = 2 µm

建立物理空间 signed distance field：

phi_M

再从：

phi_M = 0

提取 surface。

---

允许非常有限的：

signed-distance smoothing

但必须保存：

raw_mask_surface

和：

smoothed_mask_surface

不能用过度 smoothing 美化结果。

============================================================
十六、SWC_MASK_ASSISTED 第一版不要太复杂
========================================

本轮只做一个 conservative experiment。

不要重写整个 lumen algorithm。

以当前：

valid SWC-only junction surface

为初始表面。

只允许：

JUNCTION_CORE
+
TRANSITION_COLLAR

中的 surface vertices 受 Mask 影响。

远端 explicit branch 固定。

============================================================
十七、Mask influence region
===========================

定义：

junction core:
Mask influence strongest

transition collar:
Mask influence gradually decreases

explicit branch:
Mask influence = 0

构造空间权重：

w(x)

满足：

junction center:
w → 1

collar outer boundary:
w → 0

============================================================
十八、Surface-to-Mask distance
==============================

计算当前 surface vertex：

x_i

到 Mask boundary 的 signed distance：

phi_M(x_i)

Mask fitting energy：

E_mask =
mean(phi_M(x_i)^2)

============================================================
十九、SWC radius constraint
===========================

不能为了贴 Mask：

破坏 current SWC radius。

继续计算：

r_eq(s)

与：

r_SWC(s)

定义：

E_radius =
mean[
(r_eq - r_SWC)^2
]

============================================================
二十、Surface smoothness constraint
===================================

继续使用：

normal jump
dihedral angle
或 discrete curvature

定义：

E_smooth。

最终只进行局部 constrained refinement：

E =
lambda_mask E_mask
+
lambda_radius E_radius
+
lambda_smooth E_smooth

---

第一版不要自动寻找 lambda。

测试少量明确配置：

Mask influence:
weak
medium
strong

所有结果分别保存。

============================================================
二十一、不得改变 topology
=========================

C 方法优化过程中硬性要求：

source branch count unchanged

source junction connectivity unchanged

no new lumen connection

no branch disappearance

no self-intersection

no internal face

no internal cap

============================================================
二十二、A/B/C 必须统一评价
==========================

三个方法必须使用完全相同的 evaluation pipeline。

至少比较：

self-intersections

internal faces

internal caps

boundary edges

non-manifold edges

surface components

degenerate triangles

radius P95 error

collar max radius error

normal jump P95

normal jump P99

transition roughness

surface volume

triangle count

runtime

============================================================
二十三、与 Dataset Mask 的一致性指标
====================================

必须额外计算：

STL voxelized mask
vs
dataset Mask

的：

Dice
IoU

但不要只看 Dice。

同时计算：

surface-to-mask-boundary:

mean distance
median distance
P95 distance
max distance

============================================================
二十四、重点分析 Junction-local 指标
====================================

对 J49 单独输出：

junction Dice

junction surface distance

junction normal jump

junction local volume

branch-local radius fidelity

collar radius fidelity

============================================================
二十五、必须生成 Overlay Visualization
======================================

生成：

A:
SWC-only surface
+
Mask boundary

B:
Mask-only surface
+
SWC centerline

C:
SWC+Mask surface
+
Mask boundary
+
SWC centerline

使用相同 camera。

============================================================
二十六、生成 Cross-section Comparison
=====================================

选择：

parent branch
junction core
daughter branch 1
daughter branch 2

至少四个真实位置。

同一平面显示：

dataset Mask boundary
SWC target radius
SWC-only STL
Mask-only STL
SWC+Mask STL

这是本任务非常重要的人工验收图。

============================================================
二十七、特别关注 voxel anisotropy
=================================

因为原始数据 spacing：

1 × 1 × 2 µm

所以必须报告：

Mask surface 在 z 方向是否表现出明显：

staircase
radius quantization
surface roughness。

不要因为 Mask 是实验数据就自动认为：

Mask-derived surface 更平滑或更准确。

============================================================
二十八、特别关注最细血管
========================

论文数据包含约：

2 µm radius

量级血管。

这意味着：

diameter ≈ 4 µm

在 z 方向只有约：

2 voxels

跨越直径。

因此单独分析：

small-vessel junctions

是否存在严重 voxelization error。

============================================================
二十九、Mask-only 不能自动成为 Winner
=====================================

即使：

Mask Dice = 1

也不能说明：

Mask-only STL

最好。

因为如果 dataset Mask 本身由：

SWC + radius

生成，

那么这是 circular validation。

最终报告必须明确讨论：

independence of evidence。

============================================================
三十、判定规则
==============

最后按以下逻辑判断。

---

CASE 1
------

如果 Mask provenance = SWC_DERIVED

且：

SWC reconstruction 可以高度复现该 Mask，

并且：

SWC+Mask

没有明显改善：

normal jump
roughness

那么结论：

MASK_NOT_NEEDED_FOR_FORMAL_STL_RECONSTRUCTION

Mask 仅保留：

QC / visualization reference。

---

CASE 2
------

如果 Mask provenance = SWC_DERIVED

但：

Mask-assisted junction refinement

能够显著改善：

normal jump
junction smoothness

同时不恶化：

radius fidelity
topology QC

那么结论：

MASK_USEFUL_AS_DERIVED_GEOMETRIC_REGULARIZER

注意：

不能称为：

independent ground truth。

---

CASE 3
------

如果 Mask provenance = UNET_PREDICTED

并且：

它与 SWC 在 junction volume 上提供有意义的额外信息，

同时 SWC+Mask 能改善 surface quality，

则结论：

MASK_USEFUL_AS_AUXILIARY_VOLUMETRIC_PRIOR

仍然：

topology must remain SWC-defined。

---

CASE 4
------

如果 Mask 导致：

radius error ↑
voxel staircase ↑
false vessel merging
surface roughness ↑

则：

REJECT_MASK_FOR_STL_GEOMETRY

============================================================
三十一、不要立即修改正式 pipeline
=================================

本任务结束时：

current SWC-only CFD lumen pipeline

必须保持可用。

新的 Mask experiment 默认：

experimental_only = true

除非最终 quantitative evidence 明确支持，

不要将其设为正式 default。

============================================================
三十二、代码组织
================

不要把代码塞进 model_generate.py。

在现有：

ulm_3D_vascular/utils/cfd_lumen/

下新增或复用：

mask_reference.py
mask_provenance.py
mask_surface.py
mask_comparison.py
mask_assisted_refinement.py

如果已有功能相同模块则扩展，
不要重复造轮子。

============================================================
三十三、Outputs
===============

全部结果保存到：

ulm_3D_vascular/outputs/model_generate/<run_id>/mask_evaluation/

建议：

mask_evaluation/
    provenance/
        mask_inventory.csv
        mask_provenance_report.csv

    reconstruction/
        swc_reconstructed_mask.tif

    junctions/
        junction_49/
            dataset_mask.tif
            swc_only.vtp
            mask_only.vtp
            swc_mask_assisted.vtp

    qc/
        global_mask_comparison.csv
        junction_mask_comparison.csv
        abc_comparison.csv

    figures/
        mask_provenance_summary.png
        swc_vs_mask_overlay.png
        junction_49_abc.png
        junction_49_cross_sections.png
        radius_vs_mask.png
        mask_anisotropy.png

    report/
        mask_stl_reference_assessment.md
        summary.json

============================================================
三十四、最终报告必须首先回答
============================

不要只说：

“Mask 有帮助”

或：

“Mask 没帮助”。

必须逐条回答：

当前实际使用的 Mask 来源是什么？

SWC-derived？
U-Net predicted？
Unknown？

证据是什么？

当前 Mask 与 corrected SWC 是否具有独立性？

SWC + radius 能否重新复现该 Mask？

Dice / IoU 是多少？

Mask 在普通 branch 中是否提供额外几何信息？

Mask 在 bifurcation 中是否提供额外 shape information？

Mask-only surface 是否存在明显 voxel anisotropy？

J49：

SWC-only
Mask-only
SWC+Mask

哪一个：

normal continuity 最好？

三者：

radius fidelity

分别是多少？

Mask-assisted 是否会改变 topology？

必须为 NO。

是否值得把 Mask 正式加入 STL reconstruction？

最终结论只能从以下选一个：

A.
DO_NOT_USE_MASK_IN_FORMAL_STL_PIPELINE

B.
USE_MASK_ONLY_FOR_QC_AND_VALIDATION

C.
USE_MASK_AS_LOCAL_JUNCTION_REGULARIZER

D.
USE_MASK_AS_AUXILIARY_VOLUMETRIC_PRIOR

============================================================
三十五、最重要的科学解释
========================

报告必须区分：

【Mask 与 SWC 一致】

和：

【Mask 为 SWC 提供独立验证】

不是一回事。

如果 Mask 本来就是：

SWC + radius

生成的，

那么：

高 Dice

只说明 reconstruction 与 annotation-generation method 一致。

不能声称：

STL 与真实 vascular wall 达到相同精度。

============================================================
三十六、本任务 Definition of Done
=================================

只有完成以下内容才算结束：

1. 当前真实 Mask provenance 被确认；
2. 论文描述与当前实际文件来源被对应起来；
3. SWC → Mask 重建测试完成；
4. global Mask similarity 被计算；
5. J49 local Mask 被提取；
6. A/B/C reconstruction 完成；
7. topology QC 完成；
8. radius fidelity 完成；
9. normal / roughness comparison 完成；
10. Mask anisotropy 被评估；
11. cross-section comparison 图完成；
12. 所有实际数据和图片写入 outputs；
13. 正式 SWC-only pipeline 未被破坏；
14. 最终明确给出：

Mask 在该数据集中究竟应该扮演什么角色。

============================================================
最终目标
========

本轮不是为了“把 Mask 用起来”。

而是要严格回答：

【这个特定论文数据集中提供的 Mask，
相对于已经人工修订的 SWC + radius，
究竟有没有能够改善 CFD STL geometry 的额外信息价值？】

必须用真实数据和定量实验回答，
不能凭经验判断。
