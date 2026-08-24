# 基于高分辨率 fMOST 脑血管网络的多尺度代表性连通子图采样与微泡运动仿真模型构建方案

# 1. 摘要（Executive Summary）

高分辨率 fMOST 脑血管数据能够提供微米尺度的三维血管中心线、局部管径及拓扑连接，为构建解剖真实的微泡运动仿真模型提供数据基础。然而，完整脑微血管网络规模巨大，难以直接用于大批量血流动力学、微泡输运及超声仿真；简单随机截取局部区域又容易造成**血管尺度、拓扑复杂度和血流条件的系统性采样偏倚**，并将 ROI 裁切形成的人工边界误认为真实血管终点。因此，本研究拟建立一套面向微泡运动仿真的**多尺度、代表性、血流可计算的真实血管连通子图采样框架**。首先将 fMOST/SWC 数据转换为带有三维中心线、连续管径及拓扑属性的 Vascular Centerline Graph；随后采用空间窗口与拓扑邻域相结合的方式产生候选连通子图，并显式区分真实终末支与裁切边界端口。在此基础上，根据血管尺度、密度、分叉复杂度、环路特征及几何阻力等变量进行分层采样，并进一步利用标准化一维血流求解获得流量分配、速度异质性和通过时间等输运指标，实现由“形态代表性”向“形态—血流联合代表性”的二阶段筛选。最终形成具有明确边界条件和统计覆盖范围的多尺度真实脑血管 phantom library，为后续三维 CFD、microbubble transport 和 ULM 仿真提供标准化输入。

---

# 2. 立项依据与研究背景

## 2.1 痛点阐述

### 2.1.1 当前主要瓶颈已由“缺少血管几何”转变为“如何从完整真实网络中构建合适的仿真样本”

fMOST 等高分辨率三维成像已经能够获得完整或大范围脑微血管结构。现有公开数据中，血管可以被表示为**三维 skeleton、局部 radius 和 connectivity**；相关标注工具能够对分叉点、终点及连接错误进行人工修正，并以 SWC 等形式保存高精度中心线结构。 在进一步图化后，bifurcation 与 endpoint 可以作为节点，相邻节点之间的连续 skeleton segment 可以作为 graph edge。

因此，对于后续 microbubble simulation，已经没有必要通过生成模型重新创造单条血管。真正需要解决的是：

$$
\boxed{ G_{\mathrm{full}} \longrightarrow \{G_1,G_2,\ldots,G_N\} }
$$

其中 $G_{\mathrm{full}}$ 为完整真实脑血管网络，$G_i$ 为适合血流和微泡运动仿真的**真实连通 vascular subgraph**。

关键问题不再是“是否能够获得大量局部血管”，而是：

> **如何保证所采样的局部血管在尺度、拓扑、几何和血流输运特征上能够代表原始完整网络，并且在被截取以后仍具有明确、可解释的仿真边界。**

---

## 2.2 简单随机采样会导致严重的结构偏倚

脑微血管网络在空间上并非均匀分布。不同位置可能具有显著不同的：

* vessel radius；
* vascular density；
* bifurcation density；
* tortuosity；
* orientation；
* cycle/loop density；
* large-vessel/capillary composition。

如果仅在三维空间均匀随机放置固定 ROI，则高密度毛细血管区域可能贡献大量结构高度相似的样本，而低频但重要的：

* 较粗血管；
* 高不对称分叉；
* 高 tortuosity 区域；
* 多环路区域；
* 大管径—微血管过渡区；

可能被系统性低估。

因此：

$$
N_{\mathrm{sample}}\uparrow
$$

并不必然意味着：

$$
\text{population coverage}\uparrow.
$$

**采样数量与采样代表性必须被明确区分。**

---

## 2.3 ROI 裁切会人为改变血管网络拓扑

假设完整血管为：

```text
                ───────────
              /
─────────────<
              \
                ─────────────────
```

如果 ROI 边界经过两个 daughter branches：

```text
         ROI
     ┌────────────┐
─────│──────<─────│
     │       \    │
     └────────────┘
```

ROI 中会出现新的 endpoint。

然而这些 endpoint 并不是生理意义上的 terminal vessel，而是：

$$
\boxed{\text{cut boundary port}}
$$

因此必须严格区分：

$$
\text{TRUE\_TERMINAL} \neq \text{CUT\_PORT}.
$$

否则：

1. topology statistics 会发生偏差；
2. terminal density 会被高估；
3. downstream subtree size 会被错误计算；
4. 血流求解中会引入错误 outlet；
5. MB trajectory 会在人工边界处被错误终止。

因此，**边界语义保存是血管 ROI 采样区别于普通图子图采样的关键问题之一。**

---

## 2.4 仅保证形态代表性仍不能保证 microbubble transport 代表性

Microbubble motion 最终由血流决定，而不仅由几何决定。

对于一段血管，其低 Reynolds 数条件下的几何阻力近似满足：

$$
R = \frac{8\mu}{\pi} \int_0^L \frac{1}{r(s)^4}\,ds.
$$

因此半径的微小变化即可导致显著不同的流动阻力。

两个局部网络即使具有相近的：

* branch number；
* vascular density；
* total length；

仍可能因为 radius distribution 和 network connectivity 不同而表现出完全不同的：

* flow split；
* pressure drop；
* residence time；
* velocity heterogeneity；
* MB network coverage。

因此，**面向 microbubble simulation 的代表性采样不能只在 morphological feature space 中完成，而应进一步扩展到 hemodynamic/transport feature space。**

---

# 2.5 核心科学问题

基于上述分析，本研究拟回答三个问题。

### 科学问题一：如何定义能够代表完整脑微血管网络的局部采样单位？

需要确定：

* 固定物理 FOV；
* 固定 graph-hop；
* 固定 branch number；
* 固定 vascular length；

哪一种尺度更适合作为 microbubble simulation unit，或者是否需要采用多尺度联合策略。

---

### 科学问题二：有限数量的真实局部子图如何覆盖完整网络的联合形态分布？

需要解决：

$$
P_{\mathrm{sample}}(\mathbf z) \approx P_{\mathrm{full}}(\mathbf z)
$$

其中：

$$
\mathbf z=\left[ r, L, \rho_v, N_b, N_{\mathrm{bif}}, \tau, \kappa, \beta_1, \dots \right]
$$

包含血管尺度、长度、密度、分叉、tortuosity、curvature 及 cycle 等多维特征。

---

### 科学问题三：如何使所采样的局部结构同时具有代表性的血流与 MB 输运条件？

核心假说为：

> **通过“形态分层采样 + 标准化一维血流求解 + 输运特征再采样”的二阶段设计，可以使用远少于完整网络的数据量，构建在形态和血流输运空间上均具有较高覆盖度的真实 vascular phantom library。**

即：

$$
\boxed{ \text{Morphology-aware Sampling} \rightarrow \text{Flow Characterization} \rightarrow \text{Transport-aware Sampling} }
$$

---

# 3. 研究目标与内容

## 3.1 总体研究目标

建立一套**从完整 fMOST 脑微血管网络自动提取多尺度、统计代表性、边界语义明确且适合血流—微泡运动仿真的真实连通 vascular subgraph 的标准化方法**。

最终形成：

$$
\boxed{ \text{Full fMOST Graph} \rightarrow \text{Candidate Subgraphs} \rightarrow \text{Morphological Stratification} \rightarrow \text{Representative Sampling} \rightarrow \text{1D Hemodynamics} \rightarrow \text{Transport-aware Resampling} \rightarrow \text{Simulation-ready Phantom Library} }
$$

---

# 3.2 具体研究内容

## 研究内容一：建立适用于局部采样的完整 Vascular Centerline Graph

首先将原始 SWC point graph 转换为 branch-level vascular graph：

$$
G=(V,E).
$$

其中：

* (V)：junction、真实 endpoint 以及必要的结构节点；
* (E)：两个相邻关键节点之间的完整 vascular branch。

每条 edge $e_i$ 保存：

$$
e_i=\left\{ \mathbf c_i(s), r_i(s), L_i, \tau_i, \kappa_i, \mathbf t_i \right\}.
$$

分别表示：

* 三维中心线；
* 沿程半径；
* branch length；
* tortuosity；
* curvature；
* local orientation。

此外计算：

* node degree；
* local vascular density；
* bifurcation density；
* graph cycle information；
* spatial coordinate；
* source animal；
* source brain region `[若可获得]`。

**该完整 graph 作为所有采样结果的唯一母体，不再对局部结构进行人工重组或生成。**

---

## 研究内容二：建立空间—拓扑双尺度连通子图采样方法

本研究不采用单一 ROI 定义，而建立两类互补的采样单元。

### A. Physical-FOV sampling

在三维空间定义实际物理尺寸：

$$
L_x\times L_y\times L_z.
$$

建议首先测试：

$$
L_x,L_y,L_z\in\{100,250,500,1000\}\,\mu\mathrm{m}
$$

量级，最终参数依据实际 fMOST 网络尺度和 ultrasound simulation FOV 确定。

该模式更接近未来：

* CFD computational domain；
* ultrasound FOV；
* ULM imaging volume。

---

### B. Topology-conditioned sampling

随机或分层选取 anchor branch，从其邻域进行：

$$
k=1,2,3,\ldots
$$

hop expansion，或限制：

$$
N_{\mathrm{branch}}\in\{10,25,50,100,200,\ldots\}.
$$

该模式适合研究：

* bifurcation transport；
* topology complexity；
* MB path selection；
* network filling。

因此：

$$
\boxed{ \text{Physical ROI} + \text{Topological ROI} }
$$

将被共同保存，而不要求两类 ROI 完全一致。

---

## 研究内容三：建立形态—拓扑代表性分层采样体系

对于所有候选 subgraph $G_i$，构建 feature vector：

$$
\mathbf z_i=\left[ \mathbf z_{\mathrm{scale}}, \mathbf z_{\mathrm{geometry}}, \mathbf z_{\mathrm{topology}}, \mathbf z_{\mathrm{spatial}} \right].
$$

主要变量如下：

| 类别     | 代表变量                                         |
| -------- | ------------------------------------------------ |
| 血管尺度 | mean/median/min/max radius、radius quantiles     |
| 网络规模 | branch number、total length、vascular volume     |
| 局部密度 | vessel length density、volume fraction           |
| 分叉     | bifurcation number、density、angle distribution  |
| 几何     | tortuosity、curvature、taper                     |
| 拓扑     | degree、terminal number、cycle rank、path length |
| 空间组织 | orientation entropy、nearest-vessel distance     |
| 阻力代理 | $\int ds/r^4$、minimum-radius bottleneck       |

在此基础上采用：

1. **分位数分层**；
2. **多变量聚类**；
3. **rare-structure enrichment**；
4. **空间去重**；

构建形态代表性样本。

目标不是令每种结构数量完全相等，而是：

$$
P_{\mathrm{sample}}(\mathbf z)
$$

能够覆盖完整真实分布，同时避免大量重复的高频局部结构淹没低频结构。

---

## 研究内容四：建立面向 MB transport 的二阶段血流感知采样

第一阶段完成 morphology-aware sampling 后，对候选网络建立标准化 1D hemodynamic model。

进一步提取：

* mean flow；
* flow CV；
* flow-split asymmetry；
* velocity distribution；
* pressure drop；
* network resistance；
* residence-time proxy；
* accessible vascular volume；
* low-flow branch fraction。

然后在新的 feature space：

$$
\mathbf h_i=\left[ Q, v, \Delta P, R_{\mathrm{net}}, T_{\mathrm{res}}, H_Q,\ldots \right]
$$

中进行第二阶段 coverage analysis。

最终得到：

$$
\boxed{ P_{\mathrm{sample}}(\mathbf z,\mathbf h) \approx P_{\mathrm{full}}(\mathbf z,\mathbf h) }
$$

从而建立真正面向 microbubble transport 的 simulation phantom library。

---

# 4. 研究方法与技术路线

## 4.1 完整血管数据标准化

### 4.1.1 输入

主要数据：

* fMOST volume；
* SWC vascular skeleton；
* local radius；
* source mouse ID。

需要进一步确认：

| 参数                     | 状态     |
| ------------------------ | -------- |
| 动物数量                 | [待补充] |
| voxel size               | [待确认] |
| SWC 坐标单位             | [待确认] |
| radius 单位              | [待确认] |
| brain atlas registration | [待确认] |
| artery/vein label        | [待确认] |
| flow-direction label     | [待确认] |

需要注意：

**采样阶段本身不要求将 SWC parent ID 解释为真实血流方向。**

在缺乏生理 flow label 时：

$$
G
$$

首先按照 **undirected vascular graph** 处理。

真实流向将在后续边界条件和血流求解中定义。

---

# 4.2 Graph cleaning 与 branch segmentation

采用 Python + NetworkX/igraph 结合 VTK/PyVista 实现。

主要步骤：

1. 去除 duplicate nodes；
2. 修复数值误差造成的极短 segment；
3. 检查 disconnected components；
4. 定位：

   * degree = 1 endpoint；
   * degree > 2 junction；
5. 将两个相邻关键点之间所有 SWC points 合并为一条 branch；
6. 保存原始 SWC point ID，实现可追溯性。

对于 branch：

$$
e_i={p_1,p_2,\ldots,p_n},
$$

建立弧长坐标：

$$
s_j= \sum_{k<j} |\mathbf p_{k+1}-\mathbf p_k|.
$$

由此得到：

$$
\mathbf c_i(s),\qquad r_i(s).
$$

---

# 4.3 候选 ROI 的产生

## 4.3.1 Anchor selection

不建议对所有节点逐一产生高度重叠的 ROI。

可采用：

### Spatial Poisson-disk sampling

在完整网络中保持最小 anchor 间距：

$$
d_{\min}.
$$

或采用：

### Farthest-point sampling

迭代选择距现有 anchors 最远的位置。

目的是降低：

$$
\text{ROI overlap}.
$$

---

## 4.3.2 Spatial ROI extraction

对于 anchor $a_i$，建立：

$$
\Omega_i= [x_i-L_x/2,x_i+L_x/2] \times\cdots
$$

然后计算：

$$
G_i= G_{\mathrm{full}} \cap \Omega_i.
$$

只保留包含 anchor 的主要 connected component，或根据研究目的保留多个 component 并分别标记。

---

# 4.4 建立 Core–Buffer 双区域结构

为了减少 ROI 边缘截断效应，建议每个样本定义：

```text
┌─────────────────────┐
│       BUFFER        │
│    ┌───────────┐    │
│    │   CORE    │    │
│    │ analysis  │    │
│    └───────────┘    │
│                     │
└─────────────────────┘
```

其中：

* **CORE**：用于 morphology、flow 与 MB 指标评价；
* **BUFFER**：用于提供边界上下文和降低截断效应。

例如：

$$
L_{\mathrm{buffer}}=(1.2\sim1.5)L_{\mathrm{core}}
$$

作为先导实验范围，正式比例由 sensitivity analysis 决定。

这样可以避免大量指标直接在 artificial boundary 附近计算。

---

# 4.5 Cut-port detection 与边界语义保存

若原 graph edge 穿过 ROI boundary：

$$
e=(p_a,p_b)
$$

且：

$$
p_a\in\Omega,\quad p_b\notin\Omega,
$$

则计算 edge 与 ROI boundary 的交点：

$$
p_c.
$$

原 branch 在 $p_c$ 处被截断，同时创建：

$$
\boxed{\text{CUT\_PORT}}
$$

节点。

因此每个 ROI 中 node 至少分为：

```text
BIFURCATION
TRUE_TERMINAL
CUT_PORT
```

如果采用 rooted topology analysis，可进一步加入：

```text
ANCHOR
```

标签。

这一设计使：

$$
N_{\mathrm{true\ terminal}}
$$

与：

$$
N_{\mathrm{cut\ port}}
$$

能够分别统计。

---

# 4.6 Feature Engineering

每个 candidate ROI 建立统一 feature vector。

## 4.6.1 Scale features

$$
N_{\mathrm{branch}}, N_{\mathrm{node}}, L_{\mathrm{total}}, V_{\mathrm{vessel}}.
$$

## 4.6.2 Radius features

计算：

$$
r_{10}, r_{25}, r_{50}, r_{75}, r_{90}, r_{\min}, r_{\max}.
$$

不建议仅保存平均 radius，因为 microvascular radius distribution 通常存在长尾。

---

## 4.6.3 Geometry features

包括：

$$
\tau= \frac{L_{\mathrm{path}}} {L_{\mathrm{Euclidean}}},
$$

以及：

* curvature；
* taper；
* orientation entropy；
* branch-angle distribution。

---

## 4.6.4 Topology features

包括：

* node degree distribution；
* bifurcation density；
* true-terminal density；
* cycle rank。

cycle rank 可计算为：

$$
\beta_1=|E|-|V|+C,
$$

其中 $C$ 为 connected component 数。

---

## 4.6.5 Hydraulic geometry proxy

在尚未求流之前，定义：

$$
R_i^*=\int_0^{L_i} \frac{ds}{r_i(s)^4}.
$$

其中常数项暂时忽略。

进一步计算：

* network $R^*$ distribution；
* maximum $R^*$ branch；
* minimum-radius bottleneck。

该指标用于识别：

> **形态上相似但流体阻力可能显著不同的局部网络。**

---

# 4.7 第一阶段：Morphology-aware Sampling

首先从完整 candidate set：

$$
\mathcal C=\{G_1,\ldots,G_M\}
$$

中选择有限 subset：

$$
\mathcal S\subset\mathcal C.
$$

不建议简单 random sampling。

推荐采用：

### Step 1：标准化 feature space

$$
\tilde{\mathbf z}=\frac{\mathbf z-\mu}{\sigma}
$$

或 robust scaling。

### Step 2：降维用于检查而非决定

使用：

* PCA；
* UMAP；

观察 candidate distribution。

### Step 3：聚类/分层

可比较：

* K-means；
* Gaussian mixture；
* HDBSCAN；
* quantile stratification。

### Step 4：每个 stratum 内选择代表性样本

可以采用 **medoid** 而非 centroid，因为最终必须选择真实 vascular ROI。

---

# 4.8 代表性优化目标

可进一步将 sample selection 表述为一个 coverage optimization 问题。

定义候选 feature similarity：

$$
s_{ij}=\exp \left( -\frac{ \|\mathbf z_i-\mathbf z_j\|^2 }{2\sigma^2} \right).
$$

选择 $K$ 个样本最大化：

$$
\max_{\mathcal S:\,|\mathcal S|=K} \left[ \sum_{i=1}^{M}\max_{k\in\mathcal S}s_{ik} -\lambda\sum_{\substack{i,k\in\mathcal S\\i<k}}O_{ik} \right],
$$

其中：

$$
O_{ik}
$$

表示两个 ROI 的空间重叠率。

第一项鼓励：

> **覆盖完整 candidate feature space。**

第二项抑制：

> **大量高度重叠 ROI 被重复选中。**

该问题可采用 greedy facility-location algorithm 求解。

---

# 4.9 数据独立性控制

即使同一只鼠可产生数万 ROI，这些样本也不能视为独立生物样本。

因此需要保存：

```text
mouse_id
source_volume
source_region
anchor_position
global_node_ids
```

并采用两种不同评价方式。

### 对 phantom library

允许同一只鼠提供多个 ROI，但控制空间 overlap。

### 对任何机器学习模型

必须：

$$
\boxed{\text{split by animal}}
$$

而不是：

$$
\text{split by ROI}.
$$

如果动物数量有限，则采用：

* leave-one-animal-out；
* group K-fold。

---

# 4.10 第二阶段：统一一维血流求解

Morphology sampling 完成后，对候选 ROI 建立 1D network model。

每条 branch：

$$
R_i=\frac{8\mu}{\pi} \int_0^{L_i} \frac{1}{r_i(s)^4} ds.
$$

边满足：

$$
Q_i=\frac{P_a-P_b}{R_i}.
$$

节点满足：

$$
\sum_i Q_i=0.
$$

需要设置：

* blood viscosity $\mu$：[待确定模型]；
* inlet/outlet boundary conditions；
* 是否考虑 Fåhraeus–Lindqvist effect；
* 是否考虑 vessel compliance。

第一阶段建议使用 Newtonian / effective-viscosity baseline，随后增加 microvascular rheology 模型进行 sensitivity analysis。

---

# 4.11 一个建议重点研究的技术：从完整网络向 ROI 传递边界条件

如果能够在完整：

$$
G_{\mathrm{full}}
$$

上建立标准化 1D blood-flow solution，则可得到每个 point：

$$
P(x),Q(x).
$$

当截取：

$$
G_i
$$

时，每个 `CUT_PORT` 的边界条件可以直接继承：

$$
P_{\mathrm{cut}}=P_{\mathrm{full}}(x_{\mathrm{cut}})
$$

或：

$$
Q_{\mathrm{cut}}=Q_{\mathrm{full}}(x_{\mathrm{cut}}).
$$

这样：

$$
\boxed{ \text{Full-network hemodynamics} \rightarrow \text{ROI boundary conditions} }
$$

可以显著减少人为设定：

$$
P_{\mathrm{out}}=0
$$

所引入的边界偏差。

**该设计尤其适合本研究，因为所有 ROI 本来就来自同一个完整真实网络。**

如果完整网络血流边界暂时无法可靠确定，则建立：

1. standard-pressure mode；
2. standardized-flow mode；
3. resistance-boundary mode；

进行 sensitivity analysis。

---

# 4.12 Flow/Transport Feature Extraction

1D solver 完成后，对每个 ROI 计算：

### Flow magnitude

$$
Q_{\mathrm{mean}}, Q_{\max}, Q_{\min}.
$$

### Velocity

$$
\bar v_i= \frac{Q_i}{\pi r_i^2}.
$$

### Flow heterogeneity

$$
CV_Q=\frac{\sigma_Q}{\mu_Q}.
$$

### Bifurcation flow asymmetry

例如：

$$
A_Q=\frac{|Q_1-Q_2|} {Q_1+Q_2}.
$$

### Residence-time proxy

对于 path $p$：

$$
T_p=\sum_{e_i\in p} \frac{L_i}{\bar v_i}.
$$

### Flow-split entropy

$$
H_Q=-\sum_i p_i\log p_i,
$$

其中：

$$
p_i= \frac{Q_i}{\sum_jQ_j}.
$$

这些指标比单纯 geometry 更直接反映 MB transport difficulty。

---

# 4.13 第二阶段：Transport-aware Resampling

将：

$$
\mathbf z_{\mathrm{morph}}
$$

扩展为：

$$
\mathbf z_{\mathrm{joint}}=\left[ \mathbf z_{\mathrm{morph}}, \mathbf z_{\mathrm{flow}} \right].
$$

例如：

$$
\left[ r_{50}, \rho_v, N_{\mathrm{bif}}, \tau, \beta_1, R_{\mathrm{net}}, CV_Q, A_Q, T_{\mathrm{res}} \right].
$$

然后重新检查第一阶段 phantom library 是否覆盖完整 candidate set。

若出现：

```text
高阻力 / 低流速 ROI   → 覆盖不足
复杂环路 ROI          → 覆盖不足
高流量主干区域        → 覆盖不足
```

则进行 targeted enrichment。

最终得到：

$$
\boxed{ \text{anatomically representative} + \text{hemodynamically representative} }
$$

的 phantom library。

---

# 4.14 与 microbubble simulation 的接口

最终每个 sampled vascular phantom 输出：

```text
phantom_id
source_mouse
source_global_graph
3D centerline graph
radius profile
true terminals
cut ports
bounding box
core region
buffer region
morphology features
flow features
pressure BC
flow BC
```

随后进入：

$$
\text{vascular phantom} \rightarrow \text{blood flow} \rightarrow \text{MB transport}.
$$

对于 MB：

$$
\frac{d\mathbf x_{\mathrm{MB}}}{dt}=\mathbf v(\mathbf x_{\mathrm{MB}}).
$$

第一阶段可采用 branch-wise Poiseuille profile：

$$
v(r)=v_{\max} \left( 1-\frac{r^2}{R^2} \right),
$$

后续在代表性 ROI 中使用 3D CFD velocity field 对其进行校准。

---

# 4.15 采样质量评价

最终需要同时验证三个层面的代表性。

## A. Marginal distribution

比较：

$$
P_{\mathrm{full}}(x)
$$

与：

$$
P_{\mathrm{sample}}(x)
$$

对于：

* radius；
* density；
* tortuosity；
* bifurcation number；
* flow velocity。

使用：

* KS distance；
* Wasserstein distance。

---

## B. Joint distribution

对于：

$$
(r,\rho_v,N_{\mathrm{bif}},R_{\mathrm{net}},CV_Q)
$$

采用：

* MMD；
* energy distance；
* sliced Wasserstein distance。

因为：

> **边际分布一致并不保证联合分布一致。**

---

## C. Rare-structure coverage

单独统计：

* top 5% tortuous ROI；
* top 5% dense ROI；
* lowest 5% radius ROI；
* highest cycle-density ROI；
* highest resistance ROI。

评价其是否在 phantom library 中得到合理覆盖。

---

# 4.16 Sampling Size Convergence

需要回答一个实际问题：

> **到底需要多少 ROI 才够？**

令样本数：

$$
N= 50,100,200,500,1000,\ldots
$$

逐步计算：

$$
D(P_N,P_{\mathrm{full}}).
$$

若：

$$
D_N=D(P_N,P_{\mathrm{full}})\rightarrow0,
$$

则认为 sampling coverage 基本收敛。

最终选择：

$$
N^*
$$

作为在**计算成本与统计代表性之间的平衡点**。

这样 phantom 数量不是人为决定，而由收敛实验确定。

---

# 4.17 建议的软件实现

| 任务                   | 推荐工具                            |
| ---------------------- | ----------------------------------- |
| SWC parsing            | Python 自定义程序                   |
| Graph operation        | NetworkX / igraph                   |
| Spatial indexing       | scipy KDTree / PyVista              |
| Geometry processing    | VTK / PyVista                       |
| ROI clipping           | VTK                                 |
| Feature computation    | NumPy / SciPy                       |
| Clustering             | scikit-learn                        |
| Graph visualization    | PyVista + NetworkX                  |
| 1D hemodynamics        | 自定义 sparse linear solver         |
| Statistical comparison | SciPy / POT                         |
| Surface reconstruction | VTK / VMTK                          |
| Representative CFD     | Musubi / OpenFOAM / 其他已选 solver |

所有 ROI 处理步骤建议建立**deterministic pipeline + configuration file**，保证不同动物、不同尺度采样完全可复现。

---

# 4.18 技术路线的完整逻辑串联

```text
完整 fMOST / SWC 脑微血管数据
                ↓
      Graph cleaning & validation
                ↓
       Branch-level vascular graph
                ↓
────────────────────────────────
第一阶段：候选局部结构生成
────────────────────────────────
                ↓
     Spatial anchor generation
                +
      Topological anchor generation
                ↓
       Connected ROI extraction
                ↓
    CORE + BUFFER region definition
                ↓
TRUE_TERMINAL / CUT_PORT labeling
                ↓
────────────────────────────────
第二阶段：形态代表性采样
────────────────────────────────
                ↓
Morphology / topology / geometry features
                ↓
Candidate feature-space characterization
                ↓
   Stratification + rare enrichment
                ↓
 Facility-location / medoid selection
                ↓
   Morphology-representative library
                ↓
────────────────────────────────
第三阶段：血流与输运表征
────────────────────────────────
                ↓
       Standardized 1D flow
                ↓
 P / Q / velocity / resistance
                ↓
 flow split / residence / heterogeneity
                ↓
────────────────────────────────
第四阶段：输运感知再采样
────────────────────────────────
                ↓
Morphology + Hemodynamics joint space
                ↓
     Coverage / bias evaluation
                ↓
      Targeted enrichment
                ↓
────────────────────────────────
Simulation-ready vascular phantom library
────────────────────────────────
                ↓
          3D CFD calibration
                ↓
       Microbubble transport
                ↓
           RF / IQ simulation
                ↓
                ULM
```

---

## 本模块的核心逻辑凝练

**本研究中的“采样”并不是从完整血管中随意截取若干局部结构，而是将其定义为一个具有明确统计目标和物理目标的模型降维过程：**

$$
\boxed{ \text{完整真实网络} \rightarrow \text{有限数量真实子网络} }
$$

同时要求尽可能保持：

$$
\boxed{ \text{Morphology} + \text{Topology} + \text{Hemodynamics} + \text{Transport} }
$$

四个层面的代表性。

由此，后续 microbubble simulation 所研究的差异可以更可靠地归因于**血管结构和血流条件本身**，而不是由不受控制的 ROI 选择偏差造成。

因此，该模块的最终目标并不是“获得更多局部样本”，而是：

> **用尽可能有限且可控数量的真实 fMOST 连通血管子图，构建能够代表完整脑微血管形态—血流空间的标准化 anatomical phantom library，并为后续 microbubble transport 和 ULM 仿真提供具有明确边界条件的物理输入**

---

这个方案的核心动机是： 

**我们已经有了一张非常大、非常真实的小鼠脑血管“地图”，真正的问题不是再去人工造血管，而是怎样从这张大地图里挑出一批有代表性的局部血管，用来做后续血流和微泡运动仿真。**

**如果只是随便截取局部区域，可能反复抽到相似的小血管，却漏掉粗血管、复杂分叉、高弯曲或多环路等重要结构；同时，截取局部区域还会人为产生一些“断口”，这些断口并不是真实血管终点，却会直接影响血流边界和微泡运动。**

**因此，本方案把“采样”本身作为一个需要认真设计的问题：先把完整 fMOST 血管网络转换成带有连接关系、三维中心线和管径信息的血管图，再从中提取不同尺度的连通子网络，记录哪些边界是真实终点、哪些只是裁切产生的出口，并根据管径、密度、分叉复杂度、曲折程度和环路等特征，挑选能够代表整体血管多样性的样本；进一步还要计算这些局部网络的血流阻力、速度分布和流量分配，确保选出的样本不仅“长得有代表性”，而且“流得也有代表性”。**

**最终得到的不是一堆随机截取的血管，而是一套来源真实、尺度可控、结构多样、边界明确、可以直接用于微泡运动和 ULM 仿真的标准化血管模型库。**

---

You do **not necessarily need Stage 2**.

The reason I introduced Stage 2 was not that radius-based clustering is wrong, but that **radius distribution alone does not uniquely determine the flow environment experienced by microbubbles**. For your project, the cleaner strategy is actually:

> **First sample ROIs based on anatomical/geometric features such as radius distribution; then test whether those samples already span the relevant flow/MB-transport conditions. Only introduce Stage 2 if they do not.**

That distinction is important.

## Why radius distribution may already be enough

Your ultimate objective is not to reproduce every possible hemodynamic state of the brain. It is to obtain a set of **realistic fMOST-derived vascular ROIs** for MB-motion simulation.

If you cluster ROIs according to, for example,

$$
\mathbf z_i^{(r)}=\left[ r_{5},r_{25},r_{50},r_{75},r_{95}, \mu_r,\sigma_r \right],
$$

and perhaps also ROI scale,

$$
L_{\mathrm{total}},\quad N_{\mathrm{branch}},
$$

then you are already controlling one of the strongest determinants of microvascular flow because, under a Poiseuille approximation,

$$
R_e=\frac{8\mu L_e}{\pi r_e^4}.
$$

The $r^{-4}$ dependence means radius is indeed an extremely important variable.

So your intuition is reasonable:

$$
\boxed{ \text{radius-stratified real ROI sampling} }
$$

may be sufficient for building the anatomical phantom library.

---

# But radius distribution is not equivalent to hydraulic behavior

The key issue is that a radius histogram tells you **what vessel sizes exist**, but not **how those vessels are connected**.

Consider two ROIs containing exactly two vessels. Both vessels have identical radius $r$ and length $L$, so their radius distributions are identical.

For one vessel,

$$
R= \frac{8\mu L}{\pi r^4}.
$$

### ROI A: vessels connected in series

$$
\text{inlet} \rightarrow R \rightarrow R \rightarrow \text{outlet}
$$

Then

$$
R_{\mathrm{eq}}^{A}=R+R=2R.
$$

### ROI B: vessels connected in parallel

$$
\text{inlet}\; \begin{matrix} \longrightarrow R\longrightarrow\\ \longrightarrow R\longrightarrow \end{matrix} \;\text{outlet}
$$

Then

$$
\frac{1}{R_{\mathrm{eq}}^{B}}=\frac1R+\frac1R,
$$

so

$$
R_{\mathrm{eq}}^{B}=\frac{R}{2}.
$$

The two ROIs have exactly the same:

* radius histogram;
* mean radius;
* radius variance;
* even total vessel length;

but

$$
\boxed{R_{\mathrm{eq}}^{A}=4R_{\mathrm{eq}}^{B}}.
$$

So clustering by radius alone cannot distinguish them.

---

# A more biologically relevant example is the location of a narrow vessel

Suppose two ROIs have the same radius distribution:

```text
large vessels: 5
medium vessels: 10
small vessels: 3
```

In ROI A, one very narrow vessel sits on a major common pathway:

```text
                 ┌──
inlet ── narrow ─┼──
                 └──
```

Almost all downstream blood must pass through that bottleneck.

In ROI B, exactly the same narrow vessel is only a minor side branch:

```text
             ┌── narrow
inlet ───────┤
             └──────── large downstream network
```

The radius distributions are identical.

But their effects on:

$$
Q_e, \qquad v_e, \qquad \Delta P,
$$

and therefore MB trajectories, can be very different.

This is the main justification for considering topology in addition to radius.

---

# More importantly, you may not need a separate “Stage 2” to solve this

This is where I would simplify the proposal based on your current project direction.

Rather than writing:

$$
\text{Stage 1: morphology sampling} \rightarrow \text{Stage 2: flow-aware resampling},
$$

you could define the initial ROI descriptor slightly more comprehensively:

$$
\boxed{ \mathbf z_i=\left[ \text{radius distribution}, \text{branch number}, \text{total length}, \text{bifurcation density}, \text{cycle structure} \right] }
$$

For example:

$$
\mathbf z_i=\left[ r_{10}, r_{25}, r_{50}, r_{75}, r_{90}, L_{\mathrm{total}}, N_{\mathrm{branch}}, N_{\mathrm{bif}}, \beta_1 \right].
$$

Then select representative ROIs directly in this feature space.

That may already be perfectly adequate.

---

# I would therefore change the role of flow analysis

Instead of:

> “Flow analysis is required to perform the second sampling stage.”

I would make it:

> **“Flow analysis is used to validate whether anatomy-based sampling also provides adequate coverage of MB-relevant hemodynamic conditions.”**

That is much cleaner.

Your pipeline then becomes:

$$
\boxed{ \text{Full fMOST network} \rightarrow \text{candidate ROIs} \rightarrow \text{radius/topology-based representative sampling} \rightarrow \text{selected ROI library} }
$$

and only afterwards:

$$
\boxed{ \text{selected ROIs} \rightarrow \text{1-D flow} \rightarrow \text{MB transport characterization} }
$$

There is no second resampling unless the validation shows that one is actually needed.

---

# You can decide this empirically

This can be formulated as a very clean experiment.

Suppose radius/topology-based sampling produces $N$ ROIs:

$$
\mathcal S_M=\{G_1,\ldots,G_N\}.
$$

Run the inexpensive 1-D flow solver on them and calculate, for example,

$$
\mathbf h_i=\left[ R_{\mathrm{eq}}, CV_v, A_Q, T_{\mathrm{transit}} \right].
$$

Then compare the distribution of these variables against a larger reference set of candidate ROIs.

For example:

$$
D_{\mathrm{flow}}=W_1 \left( P_{\mathrm{sample}}(\mathbf h), P_{\mathrm{candidate}}(\mathbf h) \right),
$$

where $W_1$ is a Wasserstein distance.

If

$$
D_{\mathrm{flow}}
$$

is already small, then:

$$
\boxed{\text{Stage 2 is unnecessary.}}
$$

You have demonstrated that morphology-based sampling implicitly covers the relevant hemodynamic space.

If instead you discover that one region of transport space is absent—for example:

```text
               Selected ROIs
Flow regime
high velocity      ✓
medium velocity    ✓
low velocity       ✗
high heterogeneity ✗
```

then you selectively add several ROIs from those missing regimes.

That is better described as:

$$
\boxed{\text{targeted enrichment}}
$$

rather than a full second sampling stage.

---

# For your current project, I would actually recommend this simpler design

### Step 1 — Generate a large candidate ROI pool

$$
G_{\mathrm{full}} \rightarrow \{G_i\}_{i=1}^{M}
$$

while retaining real connectivity and cut-port information.

### Step 2 — Describe anatomical variability

At minimum use:

$$
\mathbf z_i=\left[ \text{radius distribution}, L_{\mathrm{total}}, N_{\mathrm{branch}}, N_{\mathrm{bif}} \right].
$$

If necessary add:

$$
\tau,\quad \beta_1,\quad \rho_v.
$$

### Step 3 — Cluster / select representative real ROIs

For example using k-medoids:

$$
\mathcal C \rightarrow \mathcal S.
$$

No generated vessels and no reassembly are required.

### Step 4 — Perform flow simulation on the selected ROIs

$$
G_i \rightarrow R_e \rightarrow Q_e,P_e,v_e.
$$

### Step 5 — Check whether flow diversity is adequately covered

Calculate:

$$
R_{\mathrm{eq}}, \quad CV_v, \quad A_Q, \quad T_{\mathrm{transit}}.
$$

### Step 6 — Only if necessary, enrich the sample set

$$
\mathcal S \rightarrow \mathcal S\cup\mathcal S_{\mathrm{rare}}.
$$

Then proceed to actual MB simulation.

---

## There is also a very useful middle ground

If you want the sampling itself to have some “flow awareness” without actually solving blood flow for every candidate ROI, you can add a cheap **hydraulic geometry descriptor**:

$$
\boxed{ R_e^*=\int_0^{L_e} \frac{ds}{r_e(s)^4} }
$$

and summarize each ROI using quantities such as

$$
\operatorname{median}(R_e^*), \qquad P_{90}(R_e^*), \qquad \max(R_e^*).
$$

This still uses only the measured geometry.

It does **not** require knowing:

* inlet;
* outlet;
* flow direction;
* physiological pressure.

But it detects some differences that a radius histogram alone misses because it incorporates both

$$
r(s)
$$

and

$$
L.
$$

For your current project, I think this may be a better compromise than building an elaborate Stage-2 sampling framework.

---

# So I would revise the original proposal

I would no longer present

> **“two-stage flow-aware sampling”**

as a necessary central method.

Instead I would write:

> **The primary phantom library will be constructed by representative sampling of real connected vascular ROIs according to vessel-caliber distributions and complementary topological/geometric descriptors. Hemodynamic simulations will subsequently be used as an independent validation layer to determine whether the anatomically selected ROIs also span the range of flow conditions relevant to microbubble transport. If systematic gaps in the hemodynamic space are identified, targeted enrichment of the phantom library will be performed.**

This is methodologically stronger because **you are not assuming in advance that a complicated Stage 2 is necessary**. You are turning that into a testable question:

$$
\boxed{ \text{Does anatomy-aware sampling already imply adequate transport coverage?} }
$$

If the answer is yes, the project becomes simpler. If the answer is no, you then have quantitative evidence justifying additional flow-aware enrichment.
