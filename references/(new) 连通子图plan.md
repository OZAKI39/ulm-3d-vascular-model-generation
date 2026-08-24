# Codex Implementation Plan：fMOST 三维脑血管 ROI 代表性采样项目

## 1. 项目目标

当前 `ulm_3D_vascular` 已经具备：

* **完整 global vessel model 的可视化界面**；
* **ROI 区域的可视化功能**；
* 能够观察从完整血管模型中截取的局部 vascular structures。

本阶段 Sampling Project 的目标不是生成新血管，也不是重新组装血管，而是：

$$
\boxed{ G_{\mathrm{full}} \rightarrow \text{大量真实 connected ROIs} \rightarrow \text{多维解剖/拓扑表征} \rightarrow \text{聚类与代表性采样} \rightarrow \text{simulation-ready ROI library} }
$$

其中：

* $G_{\mathrm{full}}$：完整 fMOST 血管模型；
* 每个 ROI 必须直接来源于真实 global vessel model；
* ROI 内血管保持其**真实几何和真实拓扑连接**；
* 最终选择一批在血管尺度和网络复杂度方面具有代表性的真实 ROI；
* 供后续：

  * blood-flow simulation；
  * microbubble transport；
  * ULM simulation；
    使用。

---

# 2. 本阶段明确不包含的内容

Codex 在本轮实现中**不要加入**：

* generative vascular model；
* synthetic vessel generation；
* segment/subtree artificial assembly；
* physiological inlet/outlet 自动推断；
* 1D blood-flow solver；
* CFD；
* microbubble trajectory simulation；
* flow-aware second-stage resampling。

这些功能只能预留接口，不能提前进入 Sampling Core。

当前阶段严格限定为：

$$
\boxed{ \text{Full real vascular network} \rightarrow \text{Representative real ROI library} }
$$

---

# 3. Codex 开始编码前必须完成 Repository Inspection

在修改代码前，Codex 应首先阅读现有工程，确认并记录：

1. global vascular model 当前的数据结构；
2. SWC / node / edge 数据的读取位置；
3. radius 信息存储形式；
4. 当前 graph 是否已经存在：

   * global node ID；
   * global edge ID；
   * parent ID；
   * branch information；
5. 当前 ROI 如何定义：

   * bounding box；
   * spatial clipping；
   * graph subgraph；
   * 或仅为可视化框；
6. ROI 当前如何从 global model 中获取数据；
7. 当前 visualization interface 如何传递：

   * ROI center；
   * ROI size；
   * selected vessel；
   * node IDs；
8. 当前是否已有：

   * graph utils；
   * logging；
   * config；
   * visualization utils；
   * test structure。

### 强制要求

**优先复用现有功能。**

不得为了 Sampling Project：

* 重写 global vessel visualization；
* 重写 ROI visualization；
* 把 sampling algorithm 塞入 UI callback；
* 破坏现有交互逻辑。

正确架构应是：

```text
Sampling Core
    ↓
Sampling Result
    ↓
Existing Visualization Interface
```

而不是：

```text
UI
 ↓
Sampling Algorithm
 ↓
Graph Logic
```

---

# 4. 建议新增的代码目录

在：

```text
ulm_3D_vascular/utils
```

下新增：

```text
ulm_3D_vascular/
│
├── utils/
│   │
│   ├── sampling/
│   │   ├── __init__.py
│   │   ├── sampling_config.py
│   │   ├── sampling_types.py
│   │   │
│   │   ├── roi_extraction.py
│   │   ├── roi_boundary.py
│   │   │
│   │   ├── radius_features.py
│   │   ├── structural_features.py
│   │   ├── roi_features.py
│   │   │
│   │   ├── feature_scaling.py
│   │   ├── clustering.py
│   │   ├── representative_selection.py
│   │   ├── overlap_control.py
│   │   │
│   │   ├── sampling_validation.py
│   │   ├── sampling_visualization.py
│   │   └── sampling_io.py
│   │
│   └── [existing utils ...]
│
├── outputs/
│   └── sampling/
│
└── [existing UI / application files]
```

如果项目已有统一：

```text
config/
tests/
scripts/
```

目录，则遵循现有工程规范，不重复创建功能相同目录。

---

# 5. Sampling 的核心数据对象

建议定义：

```python
ROIRecord
```

至少包含：

```text
roi_id
source_model_id
source_mouse_id

anchor_id
anchor_position

bbox_min
bbox_max
bbox_center
bbox_size

global_node_ids
global_edge_ids

local_node_ids
local_edge_ids

node_count
edge_count
branch_count
bifurcation_count

true_terminal_ids
cut_port_ids

radius_features
structural_features
feature_vector

cluster_id
distance_to_cluster_center

is_representative
selection_rank
```

---

## 5.1 Global traceability 是硬约束

任何 ROI 中的结构必须能够追踪回原始 global model：

$$
\boxed{ \text{ROI node} \leftrightarrow \text{Global node ID} }
$$

以及：

$$
\boxed{ \text{ROI edge} \leftrightarrow \text{Global edge ID} }
$$

不能在 clipping 后完全重新编号而丢失 global mapping。

这是后续：

* flow boundary inheritance；
* ROI 返回 global view；
* microbubble simulation；
* debug；

的基础。

---

# 6. Candidate ROI Generation

## 6.1 Anchor 不允许简单使用全部 SWC points

如果每个 SWC point 都作为 ROI center，会产生大量高度重叠样本：

```text
ROI 1
ROI 2
ROI 3
...
```

彼此仅相差几个 micrometers。

这会造成：

* 数据量虚高；
* 聚类被重复结构主导；
* computation cost 增大；
* selected ROI 空间高度冗余。

---

## 6.2 至少支持三种 anchor strategy

配置支持：

```text
random
farthest_point
poisson_disk
```

推荐默认：

```text
farthest_point
```

或 minimum-distance sampling。

示例配置：

```yaml
sampling:
  seed: 42

  anchor_mode: farthest_point

  min_anchor_distance_um: 100

  roi_size_um:
    - 250
    - 250
    - 250
```

所有参数必须通过 config 控制，不允许写死。

---

# 7. ROI Extraction 必须输出真实 Connected Subgraph

对于一个 spatial ROI：

$$
\Omega_i
$$

先提取：

$$
G_i=G_{\mathrm{full}}\cap\Omega_i.
$$

如果出现多个 disconnected components：

```text
ROI
│
├── Component A
├── Component B
└── Component C
```

默认保留：

$$
\boxed{ \text{包含 anchor 的 connected component} }
$$

作为 simulation ROI。

同时保存：

```text
raw_component_count
raw_total_vessel_length
retained_component_length
```

用于后续质量分析。

---

# 8. TRUE_TERMINAL 与 CUT_PORT 必须严格区分

这是 Sampling Project 的关键功能。

ROI 中：

$$
degree_{\mathrm{ROI}}(v)=1
$$

不能自动等于 anatomical terminal。

---

## 8.1 TRUE_TERMINAL

若：

$$
degree_{\mathrm{global}}(v)=1
$$

且该点在完整 global graph 中就是终点，则定义：

```text
TRUE_TERMINAL
```

---

## 8.2 CUT_PORT

若：

$$
degree_{\mathrm{ROI}}(v)=1
$$

但：

$$
degree_{\mathrm{global}}(v)>1
$$

或者该 branch 被 ROI boundary 截断，则定义：

```text
CUT_PORT
```

并保存：

```text
cut_port_id
global_edge_id
intersection_position
radius_at_cut
boundary_face
```

如果当前 clipping 实现不能精确得到 intersection point，则应增加几何求交模块，而不是直接使用最近节点替代。

---

# 9. ROI Representative Sampling 不是只依据 Radius

正式 Sampling Descriptor 定义为：

$$
\boxed{ \mathbf z_i=\left[ \mathbf z_i^{radius}, \mathbf z_i^{size}, \mathbf z_i^{topology} \right] }
$$

其中：

* **radius distribution**：血管尺度组成；
* **network size**：局部网络规模；
* **topology complexity**：局部连接复杂程度。

Radius 是核心，但不是唯一变量。

---

# 10. Radius Distribution Features

## 10.1 禁止直接对 SWC point radius 做普通 histogram

原因是不同 branch 的 SWC sampling interval 可能不同。

例如：

```text
Vessel A:
每 1 μm 一个 SWC point

Vessel B:
每 5 μm 一个 SWC point
```

如果直接统计 points：

$$
P(r)
$$

Vessel A 会被人为赋予 5 倍权重。

---

# 11. 使用 Arc-length Weighted Radius Distribution

对于相邻中心线点：

$$
\mathbf p_j,\;\mathbf p_{j+1}
$$

计算：

$$
\Delta s_j=\left\|\mathbf p_{j+1}-\mathbf p_j\right\|.
$$

局部 radius：

$$
\bar r_j=\frac{r_j+r_{j+1}}{2}.
$$

定义加权经验 CDF：

$$
\boxed{ F_i(r)=\frac{ \sum_j \Delta s_j \mathbf 1(\bar r_j\leq r) }{ \sum_j\Delta s_j } }
$$

由此提取：

$$
\boxed{ r_{10}, r_{25}, r_{50}, r_{75}, r_{90} }
$$

作为默认 radius feature。

同时保存但不一定全部用于 clustering：

$$
r_{05}, r_{95}, r_{\min}, r_{\max}, \mu_r, \sigma_r.
$$

---

# 12. Structural Features

## 12.1 Branch Count

定义：

$$
\boxed{ N_{\mathrm{branch}} }
$$

用于描述 ROI 的网络规模。

它直接影响：

* MB potential pathways；
* bifurcation encounter rate；
* computational complexity。

---

## 12.2 Bifurcation Count

定义：

$$
\boxed{ N_{\mathrm{bif}} }
$$

分叉是未来 MB transport 中最关键的局部结构之一，因此必须纳入默认 descriptor。

---

# 13. Total Vessel Length

定义：

$$
\boxed{L_{\mathrm{total}}=\sum_{e\in E}L_e}
$$

它描述 ROI 中实际 vascular material 的总量。

如果 ROI size 固定：

$$
L_{\mathrm{total}}
$$

可直接用于 clustering。

---

# 14. Vessel Length Density

如果未来支持不同 ROI size，则应使用：

$$
\boxed{ \rho_L=\frac{ L_{\mathrm{total}} }{ V_{\mathrm{ROI}} } }
$$

替代单纯 total length。

---

# 15. Branch / Bifurcation Density

同样，如果 ROI size 可变：

$$
\rho_{\mathrm{branch}}=\frac{ N_{\mathrm{branch}} }{ V_{\mathrm{ROI}} },
$$

$$
\rho_{\mathrm{bif}}=\frac{ N_{\mathrm{bif}} }{ V_{\mathrm{ROI}} }.
$$

因此代码必须知道 ROI size 是否固定。

---

# 16. Cycle Rank

对于 connected ROI：

$$
\boxed{ \beta_1=|E|-|V|+1 }
$$

用于表示局部独立环路数量。

解释：

$$
\beta_1=0
$$

表示 tree-like structure；

$$
\beta_1>0
$$

表示存在 loop / cross-connection。

因为 fMOST 数据可能包含 microvascular loops，所以建议默认计算。

但：

> 如果先导实验发现几乎所有 ROI 的 $\beta_1$ 都相同，则应允许从 clustering features 中关闭该变量。

---

# 17. 第一版默认 Feature Vector

如果 ROI size 固定：

$$
\boxed{ \mathbf z_i=\left[ r_{10}, r_{25}, r_{50}, r_{75}, r_{90}, N_{\mathrm{branch}}, N_{\mathrm{bif}}, L_{\mathrm{total}}, \beta_1 \right] }
$$

共 9 个核心特征。

这应当作为：

```text
radius_plus_structure
```

模式。

---

# 18. 如果 ROI Size 可变

则使用：

$$
\boxed{ \mathbf z_i=\left[ r_{10}, r_{25}, r_{50}, r_{75}, r_{90}, \rho_{\mathrm{branch}}, \rho_{\mathrm{bif}}, \rho_L, \beta_1 \right] }
$$

从而避免“大 ROI 天然拥有更多 branches”造成的 clustering bias。

---

# 19. Feature Mode 必须至少支持三种

## Mode 1 — `radius_only`

```text
r10
r25
r50
r75
r90
```

即：

$$
\mathbf z_i=\left[ r_{10}, r_{25}, r_{50}, r_{75}, r_{90} \right]
$$

### 用途

**Baseline / ablation。**

它回答：

> 仅根据 vessel caliber composition 是否已经足以得到代表性 ROI？

---

# 20. Mode 2 — `radius_plus_structure`

这是正式方案默认值：

$$
\boxed{ \left[ r_{10}, r_{25}, r_{50}, r_{75}, r_{90}, N_{\mathrm{branch}}, N_{\mathrm{bif}}, L_{\mathrm{total}}, \beta_1 \right] }
$$

或可变 ROI size 下对应 density version。

### 用途

同时描述：

* vessel caliber；
* network complexity；
* vascular amount；
* topology loops。

---

# 21. Mode 3 — `extended_morphology`

后续研究可以加入：

* mean tortuosity；
* tortuosity quantiles；
* curvature；
* bifurcation angle；
* taper；
* orientation entropy；
* nearest-vessel distance；
* vessel volume fraction。

当前版本只实现接口或基础计算，不建议全部默认参与 clustering。

原因是避免：

$$
20\sim30
$$

维特征空间导致结果难解释。

---

# 22. 默认配置建议

```yaml
features:

  mode: radius_plus_structure

  radius_quantiles:
    - 0.10
    - 0.25
    - 0.50
    - 0.75
    - 0.90

  include_branch_count: true
  include_bifurcation_count: true
  include_total_vessel_length: true
  include_cycle_rank: true

  include_tortuosity: false
  include_curvature: false
  include_bifurcation_angle: false

  scaler: robust
```

---

# 23. CUT_PORT 不作为默认生物学聚类特征

需要计算：

```text
N_cut_ports
```

但它默认不进入：

$$
\mathbf z_i
$$

因为 CUT_PORT 很大程度由 ROI clipping 位置决定，不完全反映真实 vascular biology。

它应该属于：

$$
\boxed{ \text{ROI quality / boundary metadata} }
$$

可用于 quality filtering。

例如：

```yaml
roi_quality:
  max_cut_ports: null
```

或者后期进行 sensitivity analysis。

---

# 24. Feature Scaling 是必须的

不同 feature 数值范围差异明显：

例如：

$$
r_{50}=4\,\mu\mathrm{m}
$$

而：

$$
N_{\mathrm{branch}}=100.
$$

如果直接计算 Euclidean distance，branch count 会压倒 radius。

因此必须 scale：

$$
\boxed{ \tilde z_k=\frac{ z_k-\operatorname{median}(z_k) }{ IQR(z_k) } }
$$

推荐使用：

```text
RobustScaler
```

作为默认。

Scaler 的：

```text
median
IQR
feature order
```

必须保存。

---

# 25. Feature Group Weighting 预留接口

未来可定义：

$$
d_{ij}^2=w_r \left\| \tilde{\mathbf z}_i^r- \tilde{\mathbf z}_j^r \right\|^2 + w_s \left\| \tilde{\mathbf z}_i^s- \tilde{\mathbf z}_j^s \right\|^2.
$$

其中：

$$
\mathbf z^s=\left[ N_{\mathrm{branch}}, N_{\mathrm{bif}}, L_{\mathrm{total}}, \beta_1 \right].
$$

第一版：

$$
w_r=w_s=1.
$$

不要主观规定 radius 一定比 topology 更重要。

配置可以预留：

```yaml
feature_weights:
  radius: 1.0
  structure: 1.0
```

---

# 26. Clustering Strategy

最终需要选择的是**真实 ROI**。

因此首选：

```text
K-Medoids
```

如果项目不希望增加 dependency，则使用：

```text
KMeans
↓
centroid
↓
nearest real ROI
```

即：

$$
i_k^*=\arg\min_{i\in C_k} \left\| \tilde{\mathbf z}_i-\boldsymbol\mu_k \right\|.
$$

最终 representative 必须满足：

$$
\boxed{ ROI^* \in \text{actual candidate ROIs} }
$$

不允许使用 artificial centroid 作为 vascular phantom。

---

# 27. Cluster Number 不硬编码

配置：

```yaml
clustering:
  method: kmeans

  n_clusters: 20

  exploratory_k:
    - 5
    - 10
    - 20
    - 30
    - 50
```

对于每个 $K$ 输出：

* inertia；
* silhouette score；
* cluster size。

第一版**不要求自动决定最佳 K**。

人工根据：

* statistical metric；
* visualization；
* simulation library size；

共同决定。

---

# 28. Representative Selection

每个 cluster 默认选择：

$$
1
$$

个最接近 medoid/centroid 的真实 ROI。

同时支持：

```yaml
representatives_per_cluster: 3
```

如果每 cluster 选择多个样本，则应尽可能保证：

$$
d_{\mathrm{source-space}}\ge d_{\min}
$$

避免多个 representative 实际来自几乎同一个 global 区域。

---

# 29. Sampling 必须支持两种模式

## 29.1 Distribution-preserving

目的：

> 模拟真实 fMOST 网络中的总体结构比例。

如果 cluster $k$ 有：

$$
N_k
$$

个 candidate，则：

$$
\boxed{ n_k=N_{\mathrm{target}} \frac{N_k} {\sum_jN_j} }
$$

大 cluster 多选，小 cluster 少选。

配置：

```yaml
selection_mode: distribution_preserving
```

---

# 30. Coverage-balanced

目的：

> 建立多样化 simulation benchmark。

每个 cluster 选取近似相同数量：

$$
n_k\approx constant.
$$

配置：

```yaml
selection_mode: coverage_balanced
```

这样 rare vascular structures 会得到更多权重。

---

## 两种模式不得混淆

必须在输出中明确标注：

```text
population-representative
```

或：

```text
benchmark-balanced
```

因为两者回答的是不同科学问题。

---

# 31. Spatial Overlap Control

对于 ROI $i,j$：

$$
\Omega_i,\Omega_j
$$

定义 overlap：

$$
\boxed{ O_{ij}=\frac{ |\Omega_i\cap\Omega_j| }{ \min( |\Omega_i|, |\Omega_j| ) } }
$$

或者 Jaccard：

$$
J_{ij}=\frac{ |\Omega_i\cap\Omega_j| }{ |\Omega_i\cup\Omega_j| }.
$$

配置：

```yaml
selection:
  max_selected_overlap: 0.25
```

如果当前最佳 representative 与已选 ROI 重叠过高，则选择 cluster 中第二近的真实 ROI。

---

# 32. Outputs 目录要求

用户明确要求：

> **所有功能验证所需日志和可视化结果必须保存至 `ulm_3D_vascular/outputs`。**

每次 Sampling Run 创建唯一：

```text
run_id
```

例如：

```text
20260821_121500_radius_plus_structure_k20
```

目录：

```text
ulm_3D_vascular/
└── outputs/
    └── sampling/
        └── 20260821_121500_radius_plus_structure_k20/
            │
            ├── logs/
            │   └── sampling.log
            │
            ├── config/
            │   └── sampling_config.yaml
            │
            ├── manifests/
            │   ├── candidate_rois.csv
            │   ├── selected_rois.csv
            │   └── cluster_summary.csv
            │
            ├── features/
            │   ├── roi_features.csv
            │   └── scaler.json
            │
            ├── clustering/
            │   ├── cluster_assignments.csv
            │   └── cluster_centers.csv
            │
            ├── figures/
            │   ├── global_candidate_overview.png
            │   ├── global_selected_overview.png
            │   ├── radius_distribution_global_vs_selected.png
            │   ├── branch_count_global_vs_selected.png
            │   ├── bifurcation_count_global_vs_selected.png
            │   ├── vessel_length_global_vs_selected.png
            │   ├── cluster_size_distribution.png
            │   ├── pca_feature_space.png
            │   └── silhouette_scan.png
            │
            ├── roi_previews/
            │   ├── cluster_000_rep_00.png
            │   ├── cluster_001_rep_00.png
            │   └── ...
            │
            └── report/
                └── sampling_summary.json
```

---

# 33. Logging 要求

统一使用 `logging`。

不要依赖大量：

```python
print()
```

日志至少包含：

```text
Input model
Global node count
Global edge count

ROI size
Anchor mode
Minimum anchor distance
Random seed

Candidate anchor count
Candidate ROI count
Valid ROI count
Rejected ROI count

Reject reasons:
    empty
    disconnected
    too few branches
    invalid radius
    invalid graph
    excessive overlap

Feature mode
Feature names
Feature dimension

Clustering method
Cluster number
Cluster sizes
Silhouette score

Selection mode
Selected ROI count
Spatial-overlap rejection count

Validation metrics

Output directory

Runtime:
    ROI extraction
    feature extraction
    clustering
    selection
    visualization
    total
```

---

# 34. 必须输出的 Feature Table

`roi_features.csv` 至少包括：

```text
roi_id
source_model_id
anchor_x
anchor_y
anchor_z

r10
r25
r50
r75
r90

branch_count
bifurcation_count
total_vessel_length
vessel_length_density

cycle_rank

true_terminal_count
cut_port_count

cluster_id
distance_to_center

selected
```

即使某个 feature 当前不参与 clustering，也尽量保留原始统计值。

---

# 35. 功能验收可视化

## 35.1 Global Candidate Overview

在完整 global vessel model 上显示：

* candidate ROI centers；
* candidate ROI boxes。

用于检查：

> ROI 是否过度集中在某些区域。

---

# 36. Global Selected Overview

在 global vessel model 中显示：

* selected ROI boxes；
* cluster ID。

这是 Sampling 最重要的功能验收图之一。

需要直接复用现有 visualization interface。

---

# 37. Radius Distribution Validation

比较：

$$
P_{\mathrm{candidate}}(r)
$$

与：

$$
P_{\mathrm{selected}}(r).
$$

至少生成：

1. arc-length weighted histogram；
2. weighted CDF。

不能只比较 mean radius。

---

# 38. Structural Distribution Validation

分别比较 candidate 与 selected：

### Branch count

$$
P(N_{\mathrm{branch}})
$$

### Bifurcation count

$$
P(N_{\mathrm{bif}})
$$

### Total vessel length

$$
P(L_{\mathrm{total}})
$$

### Cycle rank

$$
P(\beta_1)
$$

这样才能验证：

> 多维 sampling 是否真正保持了 network complexity。

---

# 39. Feature-space Visualization

对 scaled feature vector：

$$
\tilde{\mathbf z}_i
$$

执行 PCA。

输出：

```text
candidate points
selected representatives
cluster labels
```

目的不是用 PCA clustering，而是：

> 检查 selected ROIs 是否覆盖完整 candidate feature space。

UMAP 可以作为可选功能，不作为必须依赖。

---

# 40. Sampling Validation 独立模块

必须建立：

```text
sampling_validation.py
```

不能把 validation 逻辑写入：

```text
clustering.py
```

---

# 41. 单变量 Validation Metrics

对于连续变量，例如 radius：

### Wasserstein distance

$$
W_1(P_C,P_S)
$$

### KS statistic

$$
\boxed{ D_{\mathrm{KS}}=\sup_x \left| F_C(x)-F_S(x) \right| }
$$

### Quantile error

$$
\boxed{ E_Q=\frac{1}{K} \sum_q \frac{ |x_q^C-x_q^S| }{ |x_q^C|+\epsilon } }
$$

应分别计算：

* radius；
* branch count；
* bifurcation count；
* vessel length。

---

# 42. Cluster Coverage

定义：

$$
\boxed{ C_{\mathrm{cluster}}=\frac{ N_{\mathrm{represented\ clusters}} }{ N_{\mathrm{clusters}} } }
$$

对于：

```text
coverage_balanced
```

应接近：

$$
1.
$$

---

# 43. Spatial Redundancy Metrics

输出：

$$
O_{\mathrm{mean}}, O_{\mathrm{median}}, O_{\max}.
$$

还建议输出：

```text
nearest-selected-ROI distance
```

distribution。

---

# 44. Radius-only 与 Radius-plus-structure 必须可直接比较

这是非常重要的实验设计。

Codex 应支持用完全相同：

```text
candidate pool
seed
ROI size
K
```

分别运行：

```text
radius_only
```

和：

```text
radius_plus_structure
```

然后输出 comparison。

目的在于后续回答一个真正有价值的问题：

> **只用血管 radius distribution 是否已经足够构建代表性 phantom library？**

如果答案是 yes：

$$
\text{radius-only}
$$

未来可以作为更简单方案。

如果加入 topology 后显著改善：

* branch coverage；
* bifurcation coverage；
* MB simulation scenario diversity；

则可以证明：

$$
\text{radius + structure}
$$

的必要性。

---

# 45. Sampling Summary JSON

每次 run 输出：

```json
{
  "run_id": "...",
  "seed": 42,

  "roi_size_um": [250, 250, 250],

  "anchor_mode": "farthest_point",

  "candidate_count": 0,
  "valid_candidate_count": 0,

  "feature_mode": "radius_plus_structure",

  "features": [
    "r10",
    "r25",
    "r50",
    "r75",
    "r90",
    "branch_count",
    "bifurcation_count",
    "total_vessel_length",
    "cycle_rank"
  ],

  "cluster_method": "kmeans",
  "n_clusters": 20,

  "selection_mode": "coverage_balanced",

  "selected_count": 0,

  "validation": {
    "radius_wasserstein": null,
    "radius_ks": null,
    "branch_count_ks": null,
    "bifurcation_count_ks": null,
    "cluster_coverage": null,
    "mean_spatial_overlap": null
  }
}
```

所有数值必须由程序真实计算。

---

# 46. 与现有 UI 的集成

UI 只负责：

$$
\boxed{\text{display}}
$$

不负责 clustering。

---

## 46.1 UI 增加 Sampling Layer

允许加载：

```text
candidate_rois.csv
selected_rois.csv
```

显示模式：

```text
Show all candidate ROIs
Show selected ROIs
Show cluster X
Show ROI X
```

---

# 47. ROI Information Panel

选中一个 ROI 时显示：

```text
ROI ID
Cluster ID

Node count
Branch count
Bifurcation count

Radius:
    P10
    P25
    P50
    P75
    P90

Total vessel length

Cycle rank

True terminals
Cut ports
```

---

# 48. Cluster Inspection

选择：

```text
cluster_id = k
```

后允许查看：

* representative ROI；
* cluster size；
* cluster 内其他 candidate；
* representative distance to center。

这对人工验收 clustering 是否合理非常重要。

---

# 49. Unit Tests

至少实现以下 synthetic tests。

## Test 1 — Straight vessel

验证：

* clipping；
* radius；
* cut ports。

---

## Test 2 — Y bifurcation

```text
      ───
─────<
      ───
```

验证：

* branch count；
* bifurcation count；
* terminal type。

---

# 50. Test 3 — Boundary cutting

确保：

$$
\boxed{ \text{TRUE\_TERMINAL} \neq \text{CUT\_PORT} }
$$

---

# 51. Test 4 — Arc-length weighting

构造两个完全相同的真实 vessel geometry：

### Vessel A

```text
sampling interval = 1 μm
```

### Vessel B

```text
sampling interval = 5 μm
```

要求：

$$
F_A(r)\approx F_B(r).
$$

---

# 52. Test 5 — Branch Count

构造：

```text
one straight tube
```

与：

```text
Y tree
```

确保：

$$
N_{\mathrm{branch}}
$$

正确。

---

# 53. Test 6 — Cycle Rank

构造：

### Tree

$$
\beta_1=0
$$

### One loop

$$
\beta_1=1.
$$

确保正确。

---

# 54. Test 7 — Deterministic Sampling

相同：

```text
input
config
seed
```

必须产生相同：

```text
anchor IDs
ROI IDs
feature vectors
cluster assignments
selected ROI IDs
```

---

# 55. Test 8 — Representative Must Be Real ROI

必须满足：

```python
representative_roi_id in candidate_roi_ids
```

---

# 56. 性能要求

fMOST global graph 可能非常大，因此：

* 不允许每个 ROI 都遍历全部 graph；
* 建立 KD-tree 或 spatial index；
* ROI candidate query 使用 bounding-box / spatial search；
* radius feature 尽量 vectorize；
* graph extraction 避免重复复制大型数据；
* candidate feature calculation 支持 batch。

第一版优先：

$$
\boxed{\text{Correctness + Reproducibility}}
$$

其次才是 multiprocessing。

---

# 57. Runtime Profiling

日志中必须分别记录：

```text
model loading
anchor generation
ROI extraction
feature extraction
scaling
clustering
representative selection
validation
visualization
total runtime
```

---

# 58. 推荐 Codex 实施顺序

## Phase 0 — Repository Inspection

只阅读，不修改。

输出现有：

* graph；
* ROI；
* UI；

data flow。

---

## Phase 1 — Data Types + Config

实现：

```text
sampling_types.py
sampling_config.py
```

---

## Phase 2 — ROI Extraction

实现：

```text
roi_extraction.py
roi_boundary.py
```

完成：

$$
G_{\mathrm{full}} \rightarrow G_i.
$$

---

## Phase 3 — Radius Features

实现：

```text
radius_features.py
```

重点验证：

$$
\text{arc-length weighted radius distribution}.
$$

---

## Phase 4 — Structural Features

实现：

```text
structural_features.py
```

包括：

$$
N_{\mathrm{branch}}
$$

$$
N_{\mathrm{bif}}
$$

$$
L_{\mathrm{total}}
$$

$$
\rho_L
$$

$$
\beta_1.
$$

---

# 59. Phase 5 — Unified ROI Features

实现：

```text
roi_features.py
```

支持：

```text
radius_only
radius_plus_structure
extended_morphology
```

---

# 60. Phase 6 — Scaling + Clustering

实现：

```text
feature_scaling.py
clustering.py
```

保存 scaler。

---

# 61. Phase 7 — Representative Selection

实现：

```text
representative_selection.py
overlap_control.py
```

支持：

```text
distribution_preserving
coverage_balanced
```

---

# 62. Phase 8 — Validation

实现：

```text
sampling_validation.py
```

输出：

* KS；
* Wasserstein；
* quantile error；
* cluster coverage；
* spatial overlap。

---

# 63. Phase 9 — Outputs & Visualization

实现：

```text
sampling_io.py
sampling_visualization.py
```

所有结果必须写到：

```text
ulm_3D_vascular/outputs/sampling/<run_id>/
```

---

# 64. Phase 10 — UI Integration

只增加薄接口：

```text
sampling output
↓
UI rendering
```

禁止在 UI callback 重新计算 clustering。

---

# 65. Phase 11 — End-to-End Test

真实 global model 上完成：

```text
Load global vessel
        ↓
Generate anchors
        ↓
Extract candidate ROIs
        ↓
Compute radius features
        ↓
Compute structural features
        ↓
Scale features
        ↓
Cluster
        ↓
Select representative ROIs
        ↓
Validate
        ↓
Visualize
        ↓
Export
```

一次运行完成。

---

# 66. Definition of Done

本阶段只有同时满足以下标准才算完成：

| 功能                   | 验收要求                                    |
| ---------------------- | ------------------------------------------- |
| Candidate generation   | 可批量产生真实 ROI                          |
| Connected subgraph     | ROI topology 合法                           |
| Global traceability    | 保留 global node/edge IDs                   |
| Boundary semantics     | TRUE_TERMINAL / CUT_PORT 正确               |
| Radius features        | arc-length weighted                         |
| Branch count           | 正确                                        |
| Bifurcation count      | 正确                                        |
| Vessel length          | 正确                                        |
| Cycle rank             | 正确                                        |
| Feature modes          | 至少支持两种正式模式                        |
| Scaling                | RobustScaler 可复现                         |
| Clustering             | 可配置                                      |
| Representative         | 必须是真实 ROI                              |
| Overlap control        | 可配置                                      |
| Sampling modes         | distribution-preserving + coverage-balanced |
| Validation             | 多维 distribution validation                |
| Logs                   | 完整保存于 outputs                          |
| Figures                | 完整保存于 outputs                          |
| UI                     | 能展示 candidate / selected / cluster       |
| Unit tests             | 核心逻辑通过                                |
| Reproducibility        | 固定 seed 结果一致                          |
| Existing functionality | 不破坏现有 UI                               |

---

# 67. 为后续 Flow / MB Project 保留接口

Sampling 项目的最终 ROI 对象应能扩展为：

```python
SimulationROI
```

未来增加：

```text
boundary_conditions
flow_direction
pressure
flow_rate
velocity_field
mb_trajectories
```

但当前均保持：

```text
None
```

或未计算状态。

模块边界保持：

$$
\boxed{ \text{Sampling}:\quad G_{\mathrm{full}} \rightarrow \text{Representative Real ROI Library} }
$$

随后：

$$
\boxed{ \text{Flow}:\quad G_i\rightarrow(P,Q,v) }
$$

最后：

$$
\boxed{ \text{MB}:\quad(P,Q,v)\rightarrow\mathbf x_{\mathrm{MB}}(t) }
$$

---

# 68. 本阶段最重要的科学比较

实现完成后，应至少运行：

$$
\boxed{ \text{Radius-only Sampling} }
$$

与：

$$
\boxed{ \text{Radius + Structure Sampling} }
$$

两套完全相同 candidate pool 上的实验。

比较：

* radius distribution；
* branch-count distribution；
* bifurcation-count distribution；
* vessel-length distribution；
* cycle-rank distribution；
* feature-space coverage。

这样才能客观回答：

> **“对于后续 MB simulation，单纯使用 radius distribution 是否已经足够？”**

而不是预先假定必须使用复杂多变量方案。

---

# 69. 给 Codex 的最终任务定义

> **在不重写现有 global vessel model 与 ROI visualization interface 的前提下，为 `ulm_3D_vascular` 实现一个模块化、可重复、可量化验收的真实血管 ROI Sampling Pipeline。系统需要从完整 fMOST vascular graph 中批量提取保留 global node/edge ID 映射的 connected ROIs，并正确区分 `TRUE_TERMINAL` 与 ROI 裁切产生的 `CUT_PORT`。每个 ROI 应计算弧长加权的 vessel-radius distribution，并进一步计算 branch count、bifurcation count、total vessel length / vessel length density 以及 cycle rank 等网络结构特征。Sampling 模块必须至少支持 `radius_only` 与默认的 `radius_plus_structure` 两种 feature mode，并在经过 Robust Scaling 的多维特征空间中完成 clustering 和 representative real-ROI selection，同时支持 `distribution_preserving` 与 `coverage_balanced` 两种选择策略以及 spatial-overlap control。所有候选 ROI、特征、聚类结果、selected ROI、配置、验证指标、运行日志和功能验收可视化必须统一保存于 `ulm_3D_vascular/outputs/sampling/<run_id>/`。核心算法按照功能拆分到 `ulm_3D_vascular/utils/sampling/`，UI 仅负责调用结果和展示，不承载 sampling/clustering 核心逻辑。最终必须通过真实 global vessel model 的 end-to-end test，并提供 Radius-only 与 Radius-plus-structure 的可直接比较结果，为后续确定最简且充分的 MB-simulation phantom sampling strategy 提供定量依据。**
