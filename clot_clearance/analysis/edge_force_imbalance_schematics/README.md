# Trace the Force Imbalance Near Edges — actual-data schematics

本组图对应用户提供的三个要点，使用已有 L075 均匀拉伸诊断的完整数值数据。这里的“实际数据”是项目已经保存的数值计算结果，不是患者测量。本次仅生成图、数据摘录和来源记录，没有求解新的力学状态、演化裂纹或改变边界条件。

## 查看与使用

| 图 | 内容 | 导出 |
|---|---|---|
| 三点汇总，16:9 | 三维残余内力、邻点依赖与层分布、尚未确定的根因 | [PNG](output/figures/edge_force_imbalance_three_points.png) · [PDF](output/figures/edge_force_imbalance_three_points.pdf) · [SVG](output/figures/edge_force_imbalance_three_points.svg) |
| 1. Uniform stretch still leaves free-node forces | 全部 2,197 个粒子的三维位置；自由节点按有符号轴向内力着色 | [PNG](output/figures/01_affine_free_node_imbalance.png) · [PDF](output/figures/01_affine_free_node_imbalance.pdf) · [SVG](output/figures/01_affine_free_node_imbalance.svg) |
| 2. Neighbor dependence and inward force distribution | 四个独立分区：真实三维粒子、逐层合力、邻点依赖路径、力分解；V3 调整 (a) 标注间距 | [PNG](output/figures/02_neighbor_dependence_and_layers.png) · [PDF](output/figures/02_neighbor_dependence_and_layers.pdf) · [SVG](output/figures/02_neighbor_dependence_and_layers.svg) |
| 3. Boundary influence or a bulk operator defect? | 27 个几何完整邻域粒子、未通过的检验与根因状态 | [PNG](output/figures/03_boundary_or_bulk_unresolved.png) · [PDF](output/figures/03_boundary_or_bulk_unresolved.pdf) · [SVG](output/figures/03_boundary_or_bulk_unresolved.svg) |

[交互三维图](output/figures/edge_force_imbalance_interactive.html)自带 Plotly，可离线打开，切换内力、完整邻域和真实依赖路径，查看粒子编号及保存的力矢量。PNG 为 600 dpi；PDF/SVG 为矢量导出；带 `_preview` 的 PNG 仅用于快速预览。图中文字均为英文。

图 2 的 V2 版使用 14 × 9.5 英寸画布（PNG：8,400 × 5,700 像素），将粒子编号放在三维粒子图外，并用引线连接；曲线图例移至绘图区上方。四个分区各有独立标题和说明空间。此次仅调整图 2 的布局，原始科学数值、符号及数据来源保持一致。

V3 将 (a) 图的三个粒子标注向主图靠近，并收紧竖向间距以缩短引线；保留文字与坐标刻度的空隙。(b)–(d) 的布局和数据保持原样。

## 数据与显示约定

- L075：13 × 13 × 13 = 2,197 个粒子，参考边长 750 µm，间距 $h=62.5\,\mu\mathrm m$，作用半径 $\delta=251.25\,\mu\mathrm m=4.02h$。
- 轴向均匀拉伸 $\epsilon=0.003$，真实总伸长 2.25 µm；相场关闭，材料完好。该状态是施加的仿射运动学状态，并非已经求得的无外力平衡状态。
- 全部三维视图采用 $\mathbf x_{\mathrm{display}}=\mathbf X+50(\mathbf x-\mathbf X)$，仅放大几何位移。颜色、CSV、百分比及力的数值不随之放大。
- 图 1 灰色方块表示上下加载面共 338 个约束粒子；其余 1,859 个自由粒子使用有符号 $f_z$ 着色，蓝负、红正。静态与交互图颜色方向一致。
- 图 1 黑箭头取自自由节点中轴线上的原始三维力矢量，绘制比例为 8 µm/nN，没有将各箭头归一化为相同长度。
- 图 2 的三维线段按真实粒子坐标绘制；下方三个节点的说明图为示意排布，间隔不按比例，但两条边的参考长度标注来自数据。

## 三个要点的证据与解释

### 1. 均匀拉伸以后，靠近加载面的自由粒子仍有非零内力

自由粒子的最大内力为 **18.10829241 nN**。加载面上的内力可以由约束反力平衡，因此不能仅凭固定面粒子的内力非零就宣称失衡。在自由节点没有附加外力的情形，非零内力才表明该仿射状态不是静力平衡。

图中保存的量是内部节点力 $\mathbf f_i$，不是已经加入约束反力后的总合力，也不是为强制维持仿射状态而构造的 $-\mathbf f_i$。最大值采用 $\|\mathbf f_i\|$，着色采用 $f_{i,z}$。

### 2. 邻点系数依赖连接表面和内部多层

层曲线使用每个参考 $z$ 平面上保存的有符号内力之和：

$$
S_\ell=\frac{\sum_{i:\,z_i=-L/2+\ell h} f_{i,z}}{R_{\mathrm{cont}}},
\qquad R_{\mathrm{cont}}=P_{zz}L^2=3.3902026874998353\times10^{-6}\,\mathrm N.
$$

这里的 $R_{\mathrm{cont}}$ 是存档中的连续体参考反力，用于归一化；不是实测总反力。图保留符号、过零和反向尾部，没有取绝对值、平滑或人为构造单调衰减。

| 距底面层数 $\ell$ | $S_\ell$ | 位置 |
|---:|---:|---|
| 0 | +0.5552293285 | 受约束的底部加载面 |
| 1 | +0.5972562715 | 第一自由层 |
| 2 | +0.0393386406 | 自由层 |
| 3 | −0.1112467634 | 自由层 |
| 4 | −0.0629636823 | 自由层 |
| 5 | −0.0153442362 | 距底面 $5h=312.5\,\mu\mathrm m>\delta$ |
| 6 | −4.7340492876 × 10⁻¹⁵ | 中面对称抵消，不能据此认定不存在边界依赖 |

曲线覆盖从底面到中面的全部七层，连线只用于引导阅读。三维内力视图保留上下两半域所有粒子。

用以下真实粒子展示一条系数依赖路径；编号是数据中的零起始 ID：

| 角色 | ID | 参考坐标 $(X,Y,Z)$，µm | 邻点数 | 自身邻域 |
|---|---:|---|---:|---|
| 表面 $k$ | 1248 | (62.5, −62.5, −375) | 152 | 截断 |
| 邻点 $j$ | 1249 | (62.5, −62.5, −312.5) | 197 | 截断 |
| 内部 $i$ | 1253 | (62.5, −62.5, −62.5) | 256 | 完整 |

$|\mathbf X_k-\mathbf X_j|=h$，$|\mathbf X_j-\mathbf X_i|=4h<\delta$，而 $|\mathbf X_k-\mathbf X_i|=5h>\delta$。因此内部点无需直接连接表面点，也可能通过邻点系数读取更远处的几何信息。箭头表示系数依赖方向；没有开展逐点扰动干预或时间传播计算。

ID 1253 的保存结果满足四项分解：

$$
f_{i,z}=f^{\mathrm{direct,self}}_{i,z}
+f^{\mathrm{direct,neighbor}}_{i,z}
+f^{Z,\mathrm{self}}_{i,z}
+f^{Z,\mathrm{neighbor}}_{i,z}
=-3.619819926167023\times10^{-10}\,\mathrm N.
$$

每一项的有符号占比为 $100f^{(c)}_{i,z}/f_{i,z}$。邻点 direct 项占 **4.14365%**，邻点 Z 项占 **95.85635%**，两个 self 项接近浮点噪声。这里 Z 是现有 BAQP 力公式中的张量项；**95.9% 汇总了所有邻点的 Z 贡献，不能归因于图中画出的单一路径。** 分解来源是现有 `self_neighbor_terms`，没有重新定义力公式。

这些证据支持“非零力具有延伸到内部多层的空间分布，力公式存在邻点系数依赖”；它们不直接证明存在随时间逐层推进的传播前沿。

### 3. 边界影响还是体内算子缺陷，目前仍无法排他判断

27 个粒子的自身几何邻域完整，但满足存档中“两跳封闭”“三跳封闭”条件的粒子均为 **0**。自身完整邻域不等于整个力计算依赖范围与边界隔离。

| 存档检验 | 实际值 | 存档判据 |
|---|---:|---:|
| 几何 core 最大节点力 / $R_{\mathrm{cont}}$ | 1.0677296492 × 10⁻⁴ | 1 × 10⁻¹⁰ |
| 几何 core / collar 最大节点力之比 | 0.0199898469 | 1 × 10⁻⁴ |

原始 `bulk_patch_test=FAIL` 保持不变。已有能量一阶变分检查和固定版本参考实现均复现这些力，当前审计状态仍为 **`ROOT_CAUSE_PRIMARY=INCONCLUSIVE`**、`BOUNDARY_CONTROL_RUN=NO`。本组图没有把结果标为“已修复”或“已确认是物理边界效应”。代码推导中的 $3\delta$ 是可能的系数依赖距离上界，不是已经验证的边界层厚度或可直接施加的边界条件。

## 来源、导出表与复现

主要来源：

- [完整粒子位置和内力](../../verification/pfpd_baqp_finite_domain_boundary_consistency_v1/L075/AFFINE_STATE.npz)
- [L075 原始汇总](../../verification/pfpd_baqp_finite_domain_boundary_consistency_v1/L075/SUMMARY.json)
- [原始层分布](../../results/pfpd_baqp_finite_domain_boundary_consistency_v1/tables/L075_layers.csv)
- [原始四项分解](../../results/pfpd_baqp_finite_domain_boundary_consistency_v1/tables/L075_coefficient_channels.csv)
- [总体科学状态](../../verification/pfpd_baqp_finite_domain_boundary_consistency_v1/SUMMARY.json)
- [已有审计报告及边界理论说明](../../verification/pfpd_baqp_finite_domain_boundary_consistency_v1/REVIEW_BAQP_FINITE_DOMAIN_BOUNDARY_CONSISTENCY_V1.md)
- [四项分解的原始实现](../../pd_clot/pfpd_candidate/boundary_force_diagnostics.py)

本轮摘录表保留 SI 单位与完整精度：[粒子表](output/tables/particles_and_internal_force.csv)、[逐层表](output/tables/lower_half_signed_layers.csv)、[依赖路径](output/tables/actual_dependency_path.csv)、[选定粒子四项分解](output/tables/selected_core_force_channels.csv)、[证据汇总](output/tables/evidence_summary.csv)。

[来源与校验清单](output/manifests/SOURCE_AND_VALIDATION.json)记录全部 12 个直接输入的 SHA256、绘图数值、缩放约定和生成前后哈希一致性。[交付检查](output/manifests/ARTIFACT_CHECK.json)记录导出格式、分辨率、表行数和交互文件检查。交互图已检查 Plotly 数据结构与离线打包，当前环境没有浏览器，未执行浏览器操作测试。

在项目 `clot_clearance` 根目录运行：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/clot_edge_force_mpl_20261005 \
/home/lzy/projects/computation_examples/_shared/python_environment/bin/python -B \
analysis/edge_force_imbalance_schematics/make_schematics.py
```

只重绘图 2 时，在上述命令末尾添加 `--neighbor-only`；其他图保持原文件。导出后运行同目录的 `verify_artifacts.py`，更新交付检查记录。

依赖 NumPy、Matplotlib、Pillow、Plotly。脚本读取存档、验证邻接距离和原始汇总、导出图表，不调用生产力学求解器。生产代码、科学结果、旧测试与旧 manifest 均未修改；此前暂停的整套论文图任务未恢复。

## English captions

**Combined.** Archived L075 affine-state evidence for finite-domain force imbalance. Free-particle internal forces remain nonzero under uniform 0.3% extension. Actual neighbor links and signed layer sums demonstrate spatial coefficient dependence extending inward; they do not establish a dynamic propagation front. The geometric core fails the recorded force criteria, but no particle has a closed two-hop or three-hop coefficient neighborhood, leaving the physical root cause unresolved.

**Figure 1.** Three-dimensional internal axial forces on all free particles in the archived L075 affine state. Gray squares identify constrained loading-plane particles, whose internal forces may be balanced by reactions. The maximum free-particle force is 18.108 nN. Geometry displacements are magnified 50 times; force colors retain their actual values, and centerline force arrows use a fixed scale of 8 µm/nN.

**Figure 2.** An actual coefficient-dependency path links a surface particle to a particle with a geometrically full family through a shared neighbor. Signed layer forces remain nonzero at five lattice spacings from the bottom face, beyond one horizon. For the selected full-family particle, the sum of all neighbor Z terms accounts for 95.856% of its signed axial force; this percentage is not an attribution to the illustrated path. The evidence is static and algebraic.

**Figure 3.** A geometrically complete family does not isolate the full force stencil from finite-domain coefficients. The 27 full-family particles fail the recorded normalized-force and core-to-collar criteria, while both closed two-hop and closed three-hop core counts are zero. Existing energy-derivative and reference-implementation checks reproduce the forces. Boundary influence versus a bulk operator defect remains unresolved.
