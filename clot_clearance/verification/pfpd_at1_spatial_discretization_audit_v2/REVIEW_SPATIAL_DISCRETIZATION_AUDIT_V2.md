# PFPD_AT1 空间离散审计 V2

**结论：INCOMPLETE。几何与完整材料公式验证通过；corrected C3 的失稳已确认，C5/C7 在满足驻定条件之前触发纯拉伸域退出，不能报告空间收敛。**

仓库 `/home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930`；HEAD `f0e891686b9b107510ab3292aa27ce83638b3339`。本报告只审计数值开发试件。

## 1. 为什么重新定义 coarse geometry

此前的 27 点坐标外包络边长是 0.25 mm，但每个点都分配完整的立方体体积。因此积分中的材料总量对应 0.375 mm 立方体。新家族采用用户明确批准的物理边界节点定义。旧结果永久保存在原目录，其标签仍为 `LEGACY_3x3x3_FULL_VOXEL`。

## 2. legacy 0.375 和 corrected 0.25 有什么不同

两者的 27 点坐标相同；区别在积分权重和物理解释。legacy 总体积为 0.052734375 mm³，新值为 0.015625 mm³，比例 3.375。新 C3 由正式 BAQP 代码重新算出 1% 能量 1.5781640625 pJ、反力 1.2688125 microN；与 V1 隔离 proposal 的相对误差低于 1e-10，且未调整 Gc。详见 [重算对照](comparisons/CORRECTED_C3_PROPOSAL_CHECK.json) 和 [图 10](figures/legacy_vs_corrected_C3.png)。

## 3. 为什么 fixed delta 很重要

所有网格的物理 horizon 固定为 0.25125 mm。它规定非局部相互作用的距离，改变它就会改变模型。此次只改变 h，delta/h 依次为 2.01、4.02、6.03。delta 大于试件半边长，所以所有点的相互作用邻域均被边界截断；这不是自动判错的理由，也不能据此声称已验证体材料标定。

## 4. conservative nodal volumes 是什么

一维端点权重为 h/2，内部为 h，三维取乘积。内部、面、棱、角点分别是 h³、h³/2、h³/4、h³/8。上下表面使用二维同类权重，夹持点位于精确 z = ±0.125 mm。实现只在新脚本内建立 Cloud/Family，未改生产几何函数。

## 5. 三档 mesh 信息

| 网格 | 点 | h (mm) | delta/h | bonds | 平均 family | 独立求解解析阈值 |
|---|---:|---:|---:|---:|---:|---:|
| C3 | 27 | 0.125 | 2.01 | 185 | 13.703704 | 0.009999999999999856 |
| C5 | 125 | 0.0625 | 4.02 | 5392 | 86.272000 | 0.01000000000000003 |
| C7 | 343 | 0.0416666667 | 6.03 | 46059 | 268.565598 | 0.01000000000000029 |

各档 shape tensor 均满秩 3；完整的 family 最小/中位/最大值、表面/内部统计、K 条件数、点坐标及边界距离保存在 [geometry](geometry/LEVELS.json)。未使用秩正则化。图 1–3 及补充图 `family_size_vs_boundary_distance` 展示固定几何和采样变化。

## 6. volume/surface conservation 是否通过

通过。每档 V = 1.5625e-11 m³，A_top = A_bottom = 6.25e-8 m²，机器精度相对容差 2e-15。应力始终用反力除以真实面积 A，未使用节点数或外包体积推算面积。见 [图 2](figures/volume_and_surface_conservation.png)。

## 7. intact mechanics 是否通过

刚体平移、转动、单轴拉伸/压缩、剪切、各向同性伸长的 Fbar/Ftilde、能量密度、总力和总力矩检查全部通过。总能量等于解析能量密度乘固定体积。固定 Peridynamics.jl commit `b6e86563845a83b5cb06c805af19a817bde6b106`、Julia 1.10.12、MIT，三档均对照了仿射及非仿射状态的总能量、力、Fbar/Ftilde，误差均小于 1e-9。Manifest SHA256 `bb60d4987b8f59b7cc45eec6261f2471f92274f678871a0ae6763335bc06c808`，上游未升级，作者实现未复制到新生产代码。

**仿射完整材料公式正确，不等于自由节点已平衡。** 1% 时 C3/C5/C7 的仿射反力分别为 1.2688125、0.7919107060、0.5646225130 microN；对应名义应力 20.301、12.6705712960、9.0339602080 Pa。C5/C7 差异以 C7 为分母约 40.26%，以两者最大值为分母约 28.70%，均超过 2%。因此仿射反力采样检查 FAILED，不能用这组非驻定反力冒充 onset 反力。

仿射能量导数还包括自由节点虚功，见 [虚功分解](mechanics/AFFINE_VIRTUAL_WORK.json)。仅 C3 自由节点残余力接近零。这个区别解释了为何总能量一致而仅对 top 节点求和的反力并不一致。详见 [图 5](figures/reaction_stress_at_onset.png)。

## 8. analytic threshold 是否 mesh-independent

通过。保持 mu = 1000 Pa、bulk = 666.6666666666666 Pa、Gc = 1.2688439062499997e-5 J/m²、原 normalization.c0 = 0.25、h_kinematic = 1、sc = 0.95；Yc = Gc/(2 c0 delta) = 0.1010025 Pa，三档完全一致。每档独立用正式 BAQP 的仿射 driver 求 Y(epsilon) = Yc，阈值均为 0.01，差异小于 1e-13。此阈值是均匀仿射状态的局部 AT1 开启阈值，不是细网格已求出的驻定失稳阈值。

## 9. corrected C3 stability result

两条独立路径 A（粗步长 5e-4）与 B（2.5e-4）再二分，均得 `[0.0099921875, 0.0100000000]`，宽度 7.8125e-6，端点差为零。所有接受点沿用原 least_squares 成功、力残差 1e-9、KKT 1e-8、不可逆性、s < sc、最小 Green 应变 ≥ -1e-9、精确固定边界等门槛。原 `branch_stability_tracking.py` 求解函数未经修改，未新写分支求解算法。

三个差分尺度下重新生成 `CORRECTED_C3_DENSE_REFERENCE`，检查原始非对称误差、解析 Hss、Hxx 条件数和 Schur 消元残差。C3 阈值态 Hxx 条件数约 17.0325，最大原始对称误差 3.08e-09，最大 Hss 相对误差 9.39e-11。见 [稠密参考](matrix_free/C3_validation/DENSE_REFERENCE.json)。

临界锥要求强下界相场扰动为零，弱下界相场只能非负开启，并满足一阶能量方向导数近零。机械方向由 Hxx u = -Hxs p 消元。候选包含均匀、driver 权重、最低耦合模态的正/负部分、C3 物理中点投影、确定性随机、z 层和边界/内部开启。选中合法负见证曲率为 `[-0.025275867728484006, -0.025275867723120352, -0.025275867761147455]`，三档 Hv 也为负，一阶中性、相场可行和原纯拉伸可行性均通过。只有正的随机探针不会被宣称为稳定证明。

跨网格曲率统一用长度 L×0.01 = 2.5e-6 m、能量 Gc×L² = 7.930274414062498e-13 J 归一化，不比较旧 mesh-dependent 坐标的数值大小。C3 稠密与矩阵自由验证包含 10 个确定性随机方向、最低耦合模态和两类开启方向，每类均有三个尺度；最大相对误差 1.89609e-08 < 预声明 5e-5。残差 Jacobian 的转置对照也通过。开发时发现的分组归一化边界问题已修正并保留错误日志；未放宽物理域或验收门槛。

## 10. C5 stability result

C5 有 225 个自由位移自由度和 5392 个相场自由度。A、B 路径在起始 0.009 即触发 `PURE_TENSION_DOMAIN_VIOLATION`。最后保留的合法试点还具有约 0.071159 的归一化力残差，远高于 1e-9，不能称驻定或稳定。域外试探状态仅作为失败诊断保存，未被接受为物理分支。不存在合法临界区间、onset 反力、负临界锥见证或临界模态。

## 11. C7 stability result

C7 有 735 个自由位移自由度和 46059 个相场自由度，共 46794 DOF。单个稠密 Hessian 预计 17,517,427,488 bytes，超出当时可用 RAM 的 40%，因此禁止稠密分配；使用 `LinearOperator` 和 SciPy dogbox/LSMR，仍调用原驻定求解函数与残差定义。C7 的 A/B 路径也在 0.009 触发纯拉伸域退出。原始力残差约 0.0474434，不能接受。

额外从零状态检查 0.00025 小载荷，C5/C7 仍触发同一域门槛；零应变本身的平衡检查通过，但没有据此越过域退出点继续拼接分支。

| 尝试 | 应变 | 最后合法点力残差 | 被拒试点最小 Green 应变 | 原因 |
|---|---:|---:|---:|---|
| C5 A_000 | 0.009 | 0.071158991 | -7.7875746e-07 | PURE_TENSION_DOMAIN_VIOLATION |
| C5 B_000 | 0.009 | 0.071158991 | -7.7875746e-07 | PURE_TENSION_DOMAIN_VIOLATION |
| C5 low_start_0.00025 | 0.00025 | 0.0019509629 | -3.4633026e-05 | PURE_TENSION_DOMAIN_VIOLATION |
| C7 A_000 | 0.009 | 0.047443362 | -4.6756022e-08 | PURE_TENSION_DOMAIN_VIOLATION |
| C7 B_000 | 0.009 | 0.047443362 | -4.6756022e-08 | PURE_TENSION_DOMAIN_VIOLATION |
| C7 low_start_0.00025 | 0.00025 | 0.0013007526 | -2.9179253e-06 | PURE_TENSION_DOMAIN_VIOLATION |

所有被拒试点的完整 z 和 metrics 保存在对应 `branch_tracking/*/STATIONARY.json` 的 exit 字段。实际细网格驻定门槛失败后，按逻辑停止依赖它的细网格 Hessian 谱、相场候选 A–H、负见证、模态与定位计算；没有在非驻定状态上补算这些量来形成失稳结论。

## 12. critical intervals 是否趋于一致

INCONCLUSIVE。只有 C3 的合法区间。两条细网格加载路径均已尝试，但都未越过首个驻定门槛，因此不存在可以比较宽度、重叠或中点差的 C5/C7 区间。预声明 2e-5 中点容差及 1e-5 区间宽度要求保持不变。图 4 的 C5/C7 显示空缺，legacy 单列灰色且不进入新家族连线。

## 13. reaction stress 是否趋于一致

Onset 一致性 INCONCLUSIVE。C3 首个不稳定端点 R = 1.2688124999999e-06 N，nominal stress = 20.3009999999983 Pa。C5/C7 没有合格端点。上述仿射反力采样 FAILED 与此结论分开记录。总能量同时报告 elastic、crack、total，以及按固定 V 或 A 归一化的值；未用粒子数作为能量分母。

## 14. active-set behavior 怎样变化

C3 是 185 strong lower → 185 weak lower，onset fractions 为 strong 0、weak 1、interior 0、near_sc 0。图 6 只把接受的 C3 路径解释为 onset 行为。细网格失败起始点的活动集及 bond IDs/中点仍保存为诊断，但其 onset 行为未知。受用户指定枚举限制，`ACTIVE_SET_BEHAVIOR = DIFFUSE_ALL_BONDS` 的范围明确限定为 C3，不能理解为三档均已证明如此。

## 15. instability 是否仍然 diffuse

需要区分“所有 bonds 同时处于弱下界”和“负曲率见证的空间幅度均匀”。C3 的活动集切换是全域的，选中的合法负见证却集中在一侧夹持附近。没有 C5/C7 见证，不能回答这种空间分布是否在细化后保留。

## 16. critical mode 是否局域化

C3 bond 模态按能量配对权重聚合到粒子、归一化，再插值到共同 21³ 后处理网格。PR = 0.135939148，top10 fraction = 0.447041851；27 点的 top 10% 用 ceil(2.7) = 3 点。PR 接近 1 表示均匀，较小表示集中。这里称为数值变分见证，不能命名为裂纹或裂纹宽度。

A/B 选出的见证分别偏向上下两个对称夹持面；已报告原坐标相关性和 z 镜像后的相关性，而不是把镜像方向当作不同物理裂纹。跨网格相关性、中心/展宽张量收敛均为缺失。见 [共同网格诊断](localization/COMMON_GRID_COMPARISON.json)、图 8–9。

## 17. boundary association

C3 schedule A 幅度中心 z ≈ -0.10423661 mm，平均到 clamp 的距离约 0.020763386 mm。按 z 层、上下夹持距离、侧面距离分组的幅度和弱/内部 bond 比例已保存。此时可以说该见证与夹持位置有关，不能自动称为边界伪影；跨分辨率 PERSISTS/WEAKENS/STRENGTHENS 均无证据，所以总标签 INCONCLUSIVE。

## 18. legacy result 与新 family 的关系

legacy 与 corrected C3 的临界区间数值恰好相同，并不使其成为同一物理离散。体积、能量和反力差异仍存在。新家族只包含 corrected C3/C5/C7，不含 legacy 回归、外推或连接线。本轮所有历史文件 19556 个和符号链接 321 个均保持原哈希/目标；production modifications = 0，historical modifications = 0，historical tests modified = 0。

## 19. tests

原科学基线 255 passed / 3 known failed 全部保留；V1 的 6 项新增测试也保留，因此运行前实际基线是 261 passed / 3 known failed。V2 新增 30 项全部通过。最终全套为 291 passed / 3 known failed。原三个失败的 ID/message/traceback 与 V1 基线逐字比较（运行前及最终）：True / True。没有 xfail/skip 包装，也没有删除测试。测试通过证明实现/审计证据自洽，不会把失败的细网格科学门槛改写成通过。

初次 pytest 未限定目录而误收集了归档快照和外部 PhaseFieldX 测试，产生 collection errors；完整日志保留。随后按既有项目范围显式运行 `pytest tests`，没有在项目 tests 内排除任何测试。

新增文件仅为隔离脚本、独立测试、V2 verification/results。CPU、核数、RAM、磁盘与 Python/NumPy/SciPy/Numba/Julia 信息、每档内存预估和测得的 Python peak RSS 已存 provenance/各档目录。C7 重任务串行，BLAS 线程数 1，未分配完整细网格 Hessian。提供十组必需 PNG/PDF 和一组 family 距边界补充图，黑底、Arial regular、英文标注、无箭头。

## 20. 当前 scientific status

```json
{
  "SPATIAL_DISCRETIZATION_AUDIT_V2": "INCOMPLETE",
  "GEOMETRY_DEFINITION": "VERIFIED_0P25MM_BOUNDARY_NODE",
  "VOLUME_CONSERVATION": "VERIFIED",
  "SURFACE_CONSERVATION": "VERIFIED",
  "PHYSICAL_HORIZON": "FIXED_0P25125MM",
  "REFINEMENT_LEVELS": [
    "C3",
    "C5",
    "C7"
  ],
  "INTACT_MECHANICS_ACROSS_LEVELS": "VERIFIED",
  "ANALYTIC_THRESHOLD_CONSISTENCY": "PASSED",
  "STABILITY_INTERVAL_CONSISTENCY": "INCONCLUSIVE",
  "REACTION_AT_ONSET_CONSISTENCY": "INCONCLUSIVE",
  "CRITICAL_CONE_INSTABILITY": "INCONCLUSIVE",
  "ACTIVE_SET_BEHAVIOR": "DIFFUSE_ALL_BONDS",
  "ACTIVE_SET_BEHAVIOR_SCOPE": "C3 ACCEPTED ONSET ONLY; C5/C7 ONSET UNRESOLVED",
  "CRITICAL_MODE_LOCALIZATION": "INCONCLUSIVE",
  "BOUNDARY_ASSOCIATION": "INCONCLUSIVE",
  "CORRECTED_C3_VS_LEGACY": "DIFFERENT_PHYSICAL_DISCRETIZATION",
  "CURRENT INSTABILITY CLASSIFICATION": "INCONCLUSIVE",
  "MESH OBJECTIVITY": "NOT YET ESTABLISHED",
  "PFPD_AT1_EVOLUTION_V1": "REMAINS_NOT_VERIFIED",
  "PRODUCTION MODEL": "UNCHANGED",
  "AFFINE_REACTION_SAMPLING_CONSISTENCY": "FAILED",
  "INTACT_MECHANICS_SCOPE": "AFFINE F/ENERGY/BALANCE AND PINNED JULIA VERIFIED; FREE-NODE AFFINE EQUILIBRIUM FAILS C5/C7",
  "BLOCKER": "C5/C7 PURE_TENSION_DOMAIN_VIOLATION during stationary solve; no refinement instability bracket or fine cone/mode claim is admissible."
}
```

总状态 INCOMPLETE 的原因是 C5/C7 未通过分支驻定域门槛，而不是几何决策缺失或内存被耗尽。完整材料公式验证、空间一致性和真实血栓物理验证是三个不同问题。本轮没有建立 mesh objectivity，也没有将这个 0.25 mm 数值试件解释为真实血栓。

## 21. 唯一下一步

人工审核后，单独审计**固定 delta、仅物理表面夹持与当前“所有 bond 的最小 Green 应变必须非负”定义之间的相容性**：解释为何细网格仿射自由节点有残余力、平衡试探何处先离开纯拉伸域，并决定该开发试件的适用域。当前证据不能直接证明域内不存在任何平衡，也不能授权放宽门槛或更换边界条件。本轮在此停止；没有进入后失稳分支、arc-length、h(s)<1、真实断裂、CFD 或模型修复。
