# BraVa P1 加密路线最终结论（已退休）

当前路线 RETIRED_NO_PRODUCTION_WINNER；不恢复18.75/17.5/15 µm full-domain CFD。历史已完成11个合格研究CFD、12个46初态singleton批次；没有production winner，没有dt收敛证据。

## 当前网格结果

|候选|tetra|四壁点数|门槛|说明|
|---|---|---|---|---|
|M1_coplanar_medium|804079|0|PASS||
|M1_fixedplc_medium|762096|0|PASS||
|M1_medium|734025|0|PASS||
|M1_coplanar_fine|1064977|0|PASS||
|M1_coplanar_fine_h020|1414875|0|PASS||
|M1_fine|796909|0|PASS||
|M1_fixedplc_fine|895647|0|PASS||
|M2_medium|1851317|见 raw QA / 尚无有效网格|FAIL|Surface and tetrahedral volume disagree|
|M2_medium_frame180|1445265|0|PASS||
|M2_fine|2457761|见 raw QA / 尚无有效网格|FAIL|Surface and tetrahedral volume disagree|
|M2_fine_frame180|1748487|0|PASS||
|M3_medium|1881909|见 raw QA / 尚无有效网格|FAIL|Surface and tetrahedral volume disagree|
|M3_medium_frame180|1476965|0|PASS||
|M3_fine|—|见 raw QA / 尚无有效网格|FAIL|{'effective_options': {'boundary_layer_inside': 'True', 'sphere_refinement_on': 'True', 'no_merge': 'True', 'global_edge_size': '0.0002', 'local_edge_size_on': 'False', 'no_bisect': 'True', 'sphere_refinement': "{'edge_…|
|M3_fine_frame180|1949538|见 raw QA / 尚无有效网格|FAIL|Overlapping/folded mesh: adjacent tetra on same side of internal face|
|M3_fine_h075_frame180|1832190|0|PASS||

`frame180` 只表示 mesher 工作坐标作 diag(-1,-1,1) 刚体变换，所有输出已逆变换回正式坐标。旧版本在原姿态生成的边界层内部有同侧重叠；不能因全体 Jacobian 为正而接受。合法的三、四层网格分别由真实 native TetGen/VMTK 生成，没有删单元、补流速或把无效体域当 CFD。

## CFD 与 replay

|候选|CFD 执行|独立后处理|轨迹批次|
|---|---|---|---|
|M1_coplanar_medium|PASS|PASS|COMPLETE|
|M1_fixedplc_medium|PASS|PASS|COMPLETE|
|M1_medium|PASS|PASS|COMPLETE|
|M1_coplanar_fine|PASS|PASS|COMPLETE|
|M1_coplanar_fine_h020|PASS|PASS|COMPLETE|
|M1_fine|PASS|PASS|COMPLETE|
|M1_fixedplc_fine|PASS|PASS|COMPLETE|
|M2_medium_frame180|PASS|PASS|COMPLETE|
|M2_fine_frame180|PASS|PASS|COMPLETE|
|M3_medium_frame180|PASS|PASS|COMPLETE|
|M3_fine_h075_frame180|PASS|PASS|COMPLETE|

每个 batch 的 source / accepted states / terminal / manifest 都单独保存。失败、supported stationary、观察截止、出口终止分别计数；没有把异常删除来改善外观。


## 25/20 µm最终对照

零速tetra两档均为0，但u/G/WSS相对L2仍为1.6801%/16.5887%/22.1126%，全部轨迹指标通过0/46。25 µm出流41、截止5；20 µm出流44、截止2，两档数值失败均0。消除零速结构不等于输运收敛。25 µm全部407106个已存状态位置的壁面查询，两档最大gap差8.51e-18 m、法向差0、接触数量变化0；这只覆盖已查询位置。

# 两个已定位异常的 P1 表达机制

这里解释当前冻结流场中的两个具体结构，不把所有 slow→fast 事件都归为同一种原因。数值由原始四节点重算，见 `00_baseline/P1_REGRESSION_DERIVATION.json` 和 `REGRESSION_TARGETS.json`。

P1 四面体速度是四个节点速度的形函数加权和：`u(x) = Σ Ni(x) ui`。

**138584：四个速度自由度都被壁面 no-slip 固定为零。** 因而整个正体积单元的速度严格为零，不只是四个顶点为零。这与单元形状质量是否很好是两个问题；该单元 mean-ratio 质量约 0.717。压力收敛、全局入口出口守恒和 CPU/GPU 一致都不能恢复缺失的速度自由度。减小粒子 dt 也不能让这个冻结场在单元内部出现非零速度。原失败轨迹的空间归属和近零速度邻域需要结合此前失败审计理解，不能解释为生物黏附或正常捕获。

**388769：三个壁点为零，只剩一个自由速度节点。** 设该节点速度为 `uf`，其对面 WALL 的内法向为 `n`，对面高度为 `h`。此时：

```
u(x) = Nf(x) uf
G = uf ⊗ ∇Nf
∇Nf = n / h
div(u) = (uf · n) / h
eta = |uf · n| / |uf|
```

因此，在这个特定的三壁点单元中，eta 接近 1 直接意味着自由节点速度几乎沿壁面法向，并非只是一项抽象比值。原始节点重算得到 `eta ≈ 0.994887`、`div(u) ≈ −827.429 s⁻¹`，自由节点速度约 0.15768 m/s，几乎全部指向壁面；单元中心的法向速度约 −0.03922 m/s。该单元 mean-ratio 质量约 0.871，不需要瘦劣单元也会发生这种近壁离散异常。

这解释了为什么要同时改进近壁分辨率和检查四壁点拓扑。P1/P1 稳定化求解的全局/弱形式约束不能直接当作每个四面体逐点无散度证明。边界层加密能否减少这种误差，仍以同一物理区域、固定采样点和微泡输运收敛的实际结果判断。

对微泡而言，现有过阻尼模型根据当前局部流场、阻力和接触状态计算速度；壁间隙增大后，泡心可能进入更快的流层，出现 slow→fast。局部离散法向速度和壁面三角片的接触切换都会影响离壁时机。部分事件在 eta<0.5 的窗口也存在，不能用单一筛查阈值判定物理真伪。需要空间和时间步收敛；给速度限幅、补 epsilon 或平滑轨迹会改变要验证的问题。

M1 的局部尺寸约束仍由 TetGen 执行全域重划分，不能将其视为只替换两个异常单元的隔离实验。ID440 在 M1 中档出流伴随路径改变，不能单凭出流就宣称原位置的连续物理解已验证。四壁点与 exact-zero 结构是否消失，要由对应网格和完整 frozen-flow 审计独立确认。

进一步的[几何对应审计（历史证据已归并；原路径 `01_mesh_capability/REGRESSION_GEOMETRY_CORRESPONDENCE_ZH.md`）确认：M1/M3 的球形加密还改变了局部壁面近似，原 138584 中心及历史 ID440/1505 失败终点落在新域外。原零速体积内 4096 个固定均匀探针中，只有约 11%–12% 留在这些新域内；M2 全部保留，且上述固定点壁距与 M0 相同。因此 M1/M3 的零速结构消失混合了几何和体网格变化，不能独立证明同一位置的输运场得到修复。M2 是本轮更受控的拓扑/近壁分辨率对照；其 CFD 与轨迹必须完成后才能评价改善程度。


## 补充：零速单元不是封闭的独立几何腔室

从原始边界三角形重新枚举可确认：138584 有 2 个 WALL 面和 2 个内部共享面，6 条边中 5 条属于 WALL 三角网格、1 条为内部边（零基节点 18690–16114）。它在几何上仍与邻近流体单元连通；但四个顶点速度都被无滑移条件约束为零，使 P1 速度在整个单元及两个内部面上恒为零。正体积、几何连通与速度自由度的有效表达必须分别检查。388769 有 1 个 WALL 面、3 个非 WALL 面。原始编号及拓扑清单见 [BOUNDARY_DOF_TOPOLOGY.json](00_baseline/BOUNDARY_DOF_TOPOLOGY.json)。这项核对没有运行或验证高阶 CFD 求解器。

## M2 中档完成后的限定解释

M2 中档已保留原几何位置并消除严格零速单元；原 138584 中心现在有约 23.44 mm/s 流速。原 388769 中心的 eta 约 0.536，abs(div) 约 1361 s⁻¹。两处新所属单元都没有 WALL 顶点，不能直接套用上文只适用于三壁点单元的单自由节点公式。审计 layer=3 仅表示 W0/W1/W2 之外，并不等于远离壁面；两位置实际壁距约 47.40 / 26.21 µm。近壁物理距离、拓扑邻接层、速度自由度与局部散度应分别评估。详见 M2 中档同位置审查（历史引用 `M2_MEDIUM_FIELD_REVIEW_ZH.md`；结论见本报告和最终摘要，旧路径已清理）。

M2 中细档场与轨迹现已完成：固定点 ALL 的 u/G 相对 L2 差异分别约 0.4761%/26.5136%，高 eta 区 G 差异约 55.7082%。两档虽均有 46 条出流，但仅 1/46 条通过全部原输运指标，6 条在两档均出流的轨迹改变出口。详见 M2 中细档场审查（历史引用 `M2_MEDIUM_FINE_FIELD_REVIEW_ZH.md`；结论见本报告和最终摘要，旧路径已清理）和轨迹审查（历史引用 `M2_MEDIUM_FINE_TRAJECTORY_REVIEW_ZH.md`；结论见本报告和最终摘要，旧路径已清理）。零速结构消除与数值失败减少，不等于局部场或粒子输运已收敛。

## 历史失败的进入位置与两种保护检查

[旧正式 1565 条结果的原因复核](../population_v1/analysis/failure_speed_review_20260929/CAUSE_ANALYSIS_ZH.md)进一步确认，ID440/1505 原先都处于相邻单元 578387，渐近逼近其与 138584 共享的内部零速面。578387 的非零梯度特征值约 −443.471477 s⁻¹，对应冻结离散场约 2.255 ms 的衰减时间尺度；这不是已经验证的生理滞留时间。

两条历史末点距离该共享面约 4.1–4.8×10⁻¹⁶ m；相对 138584 的最小重心坐标约 −4.37×10⁻¹² / −5.13×10⁻¹²，在当时约 5.38×10⁻¹² 的包含容差内。多个单元都容差包含时，定位规则取较小 canonical ID，记录因而变为 138584。没有证据支持宏观穿面跳跃或穿壁；即使保持相邻单元求值，末点流速也只有约 2.2–2.5×10⁻¹³ m/s。不能仅调整定位容差来解决上游零速结构。

历史终止保护为 `ROUNDOFF_SCALE_STAGNATION_64_ACCEPTED_SUBSTEPS`；归入零速单元后分别记录 82、71 个严格零速样本，检查调用时序使样本数不必恰好为 64。另做的 CPU **32 状态**接触支持复核是另一项检查，两者不能混称。代码还另有 `EXACT_STATIONARY_COORDINATES_32_NOMINAL_STEPS` 分支，不能据此改写这两条历史终止记录。复核接触秩为 0，末态壁间隙约 13.630/8.026 µm，远大于 2 nm 交接尺度，因此保留数值失败分类，不能改称接触支持静止或吸附。该结论仅针对旧正式批次，不把当前单泡重放与旧耦合轨迹混为一批。

## 固定原边界对照的作用

固定原始 WALL 和四个端口三角形的中档 CFD 也消除了严格零速单元，原 138584 中心仍在流体内且速度约 17.4271 mm/s；原 388769 中心 div 仍约 −835.198 s⁻¹。这排除了将原零速位置移到域外这一解释，但不是只改变一个单元的孤立因果实验。最新固定边界中细档结果及完整 gate 以动态场与出生输入审查（历史引用 `FIXED_PLC_FIELD_REVIEW_ZH.md`；结论见本报告和最终摘要，旧路径已清理）和[总报告](FINAL_VALIDATION_REPORT_ZH.md)为准。

## 身份及保留范围

[机器摘要](FINAL_P1_SUMMARY.json)。昂贵的原始mesh/CFD/replay和全部原始初态继续保存。旧1565仍为458/527/437出流、141截止、2数值失败。旧路线源码从清理前Git提交恢复，未跟踪的大型日志删除前只提取必要片段；不宣称Git保存了所有历史数据。
