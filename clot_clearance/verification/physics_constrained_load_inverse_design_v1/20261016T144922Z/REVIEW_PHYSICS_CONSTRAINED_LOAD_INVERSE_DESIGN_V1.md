# PHYSICS_CONSTRAINED_CLOT_LOAD_INVERSE_DESIGN_V1

本轮实现并运行了局部流体牵引 → 原始节点映射 → Optuna 提交 → 原生 LAMMPS 求解 → 实际拓扑/能量读回 → BoTorch 下一候选的最小闭环。完成四个有效小试样；保留一次 run 0 初始化接口拒绝。**没有新软化、新断键或真实脱离分量；没有运行新完整 V3。** 这不是已找到块状碎裂载荷的结果。

## 实际观测

| Trial | 驱动速度 mm/s | 开始关断 ms | 实际峰值剪切 Pa | 冻结代理目标 | 新软化 | 历史最大自由分量 | 原生墙钟 s |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.5434 | 0.6636 | 0.06173 | 3.9888e-07 | 0 | 0 | 96.93 |
| 2 | 1.7282 | 0.6931 | 0.06912 | 4.68681e-07 | 0 | 0 | 95.90 |
| 3 | 1.0390 | 0.6556 | 0.04155 | 2.65497e-07 | 0 | 0 | 100.93 |
| 4 | 1.9772 | 0.7218 | 0.07908 | 5.66319e-07 | 0 | 0 | 105.35 |

代理目标最高的 trial 4 仍只有弹性伸长前兆。所有 ≥14、≥34、≥114 节点分量、脱离质量、1–2 节点碎屑质量和不可逆损伤耗能均为零；不存在可报告的自由块体内部活键比例或二次裂解。不能将这些零值转化为声流普遍不能碎裂血栓的结论。

## 十五项科学问题

1. **是否有流体牵引模板？** 有：直接解一维非稳态 Stokes 扩散方程的 Fourier 模态 ODE，形成可查询时间场和真实 V3 面片 NPZ。下壁静止、上壁平移，侧向无限，初始流体静止，零压力梯度。物理来源为明确的流体边界值问题；微泡驱动尚未求解。
2. **与 Gaussian 的差别？** 没有人工指定应力幅度/空间 Gaussian 宽度或任意混合方向；参数为上壁速度和关断时刻，剪切由 μ∂u/∂y 得到，启动由黏性扩散决定。上壁运动是理想化平均流驱动，不能称作已实现的微泡 Streaming。V3 上投影为独立切平面，尚未满足全曲面流动连续性。
3. **压力、剪切、速度和时间依据？** μ=0.001 Pa·s、ρ=1000 kg/m³，H=25 μm，速度范围0.5–4 mm/s，Re_H≤0.1。稳态剪切0.02–0.16 Pa；p增量严格0 Pa。驱动0时刻开始，0.1 ms半余弦上升，0.4–0.8 ms开始0.1 ms关断，随后继续计算0.3 ms黏性尾流。H²/ν=0.625 ms，最慢模态时间0.063326 ms，5倍为0.316629 ms。不是Cho实验建立时程，也没有第一步施加稳态峰值。独立有限差分与模态解的峰值归一化应力差0.000247753。
4. **真实V3映射？** 使用含栓CFD的3611个真实暴露面片、原始网格三角形节点ID及1927原始材料ID。它与上一轮4918个几何法向面片的离散来源不同，不进行身份混同。复用历史稀疏支持；共面/病态支持仅补入邻近原始材料点。原V3库、材料和历史数据未改变。
5. **力、矩、功率？** 每面片 Σw=1、ΣwX=c，力为W(A t)，面片速度为Wᵀv。相对误差依次为1.21e-15、2.81e-15、8.94e-14。同参考点、参考面积/法向；标记FROZEN_SURFACE_ONE_WAY。权重约50.31%为负，范围[-3.821,3.939]，列绝对和最大9.929。守恒不等于局部载荷精度认证，变形后真实空间力矩/真实壁面功率尚未保证。
6. **Optuna提交与恢复？** 实际SQLite ask/tell，Python独立子进程运行原生LAMMPS。完成两条有效记录后进程退出，再加载DB及pickle采样器状态继续；每个阶段从磁盘重新加载并核对状态和值。失败初始化未作为零分样本。支持中断后重新挂接RUNNING试验；已有SUMMARY直接读回。未完成的单个原生轨迹重启从原参考态执行，不宣称本轮实现了半途机械checkpoint续算。
7. **BoTorch样本够吗？** 不够。三条有效随机初始观测后实际拟合SingleTaskGP，运行带Re及能量约束的qLogExpectedImprovement，提出并完成第四条真实原生计算；四条观测后另存NEXT_CANDIDATE.json但不运行。只验证可运行反馈接口，不声称GP预测可靠、已全局优化或高质量贝叶斯优化。
8. **真实小试样？** 四组，原始1792材料点+336固定锚点、原始预裂纹与endpoint-star、G/K/s0/sc及库完全复用。加载改为平顶面Stokes剪切力，解除旧位移夹持，材料点全部积分；底部原三层锚点固定。这是局部有限试样机制筛选，不能预测V3碎块尺寸。没有为结果新增预切或调整材料。
9. **裂纹集中/凝聚候选？** 未观察到新增损伤，因此未得到凝聚裂纹证据。冻结上半试样作为诊断区域，剩余z=0韧带为候选边界；代理为边界平均损伤减内部平均损伤，加0.01倍边界平均clip(κ/s0,0,1)，取全历史最大值，预裂键排除。只剩伸长前兆时明确如此标注。没有事后改变目标。
10. **独立多节点块体？** 没有。每个新断键事件都会基于原始活键做锚点可达性分析；本轮无新断键，末态原生gather_bonds与参考存活端点逐一核对通过。原有预裂键不计作新增断裂。
11. **有效载荷需求合理吗？** 已试速度/剪切在规定局部低Re范围内，但没有识别出能产生块体的有效需求或断裂阈值。不无限线性放大，不反推一个未经运行的临界剪切。表中是实际试验需求，不是已达成碎裂的需求。
12. **MEVO能直接实现吗？** 尚不能评价。实际读取历史F1 k-Wave复基底及当前F10S的复压力、APECSS有效叶单元、32粒子位置/有限差分速度和一/二次Bjerknes力。当前F10S为HOMOGENEOUS_WATER_CONTROL_PROXY、STRICT_ZERO背景流，覆盖3 ms。当前径向数据库外包范围R=1.779–1.781 μm、f=2.499–2.501 MHz、P=0.1–40000 Pa，还须落入有效叶。缺少机械剪切需求→声压的逆映射，不能凭粒子平移把牵引标为已实现。历史k-Wave与当前F10S装置不同。接口分别返回SUPPORTED_BY_AVAILABLE_MODEL（仅有效径向插值）、CONDITIONALLY_PLAUSIBLE（局部流体/时间范围）、OUTSIDE_VALIDATED_RANGE、NOT_YET_EVALUABLE；总体为最后一项。
13. **主要不确定性？** 微泡/群集至近壁平均速度梯度的转换、曲面流动及边缘误差、局部有符号节点映射的载荷精度、材料Gc/界面参数、尺寸收敛。当前全体弱弹性结果无法区分更强但仍物理允许的局部流动是否有效。
14. **下一轮具体机制？** 优先以现有APECSS谐波和群集状态为输入，建立单一局部近壁Stokes平均流/受迫流动问题，包含实际壁面与血栓几何，再生成空间非均匀剪切。先验证源项、时间建立及表面牵引，不将中心速度代替流速；若需要惯性射流，应另立Navier–Stokes模型，不扩大本轮低Re标记范围。
15. **成本节约多少？** 本轮四个有效小试样实测墙钟合计399.111 s、CPU合计399.002 s、最大RSS 0.613 GiB，共17339步。上一轮单个2.33 ms完整V3实测1797.796 s；几何、物理时长和检查开销不同，不报告等工作量加速倍数，也不把未执行的V3成本当作实测节约。本轮新增完整V3次数0。

## 状态与复现

```json
{
  "task": "physics_constrained_load_inverse_design_v1",
  "run_id": "20261016T144922Z",
  "TRACTION_MODEL_STATUS": "PASS_LOCAL_UNSTEADY_STOKES_APPROXIMATION",
  "FORCE_MAPPING_STATUS": "PASS_REFERENCE_FORCE_MOMENT_POWER_SIGNED_WEIGHTS",
  "OPTIMIZER_STATUS": "REAL_CLOSED_LOOP_INITIAL_DESIGN_BOTORCH_EXPLORATORY",
  "REAL_NATIVE_EVALUATION_STATUS": "FOUR_COMPLETED_ONE_RETAINED_INIT_FAILURE",
  "CHUNK_FORMATION_STATUS": "NOT_OBSERVED",
  "MEVO_FEASIBILITY_STATUS": "NOT_YET_EVALUABLE",
  "FULL_V3_RUN_STATUS": "NOT_RUN_EVIDENCE_INSUFFICIENT",
  "SCIENTIFIC_CONCLUSION": "Within the bounded low-Re plane-flow surrogate, four real native evaluations show no new damage or detached components. A persisted physical-load optimization loop works, but neither chunk-producing demand nor MEVO streaming realization has been identified.",
  "material_Gc": "UNVALIDATED",
  "interface_parameters": "UNCALIBRATED",
  "macro_chunk_size_convergence": "UNVALIDATED",
  "curved_V3_flow": "UNVALIDATED_LOCAL_PATCH_APPROXIMATION",
  "best_proxy_trial": 4,
  "new_full_V3_cases": 0
}
```

见 `REPRODUCE.md`、`OPTIMIZATION_STUDY.db`、`OPTIMIZATION_TRIALS.csv`、`NUMERICAL_RESULTS.json`、`TRACTION_MODEL_SCOPE.md`、`MEVO_FEASIBILITY.json` 和图表索引。设计冻结后只有初始化回调名称适配修正，记录于PORT_ADAPTER_FIX.json；物理域、目标、材料和原生库不变。完整V3载荷提供器已离线验证原始ID、锚点零力及背景力保留，但没有宣称做过本轮V3动力学验收。

## 文献与适用边界

Oh 2021本地全文已读：第4.1节以指定点对伸长作为裂纹路径代理，第4.2节以有限差分核对伴随灵敏度。本轮借鉴预定义区域/伸长前兆，未移植二维形状伴随到三维不可逆原生系统。[论文](https://doi.org/10.1016/j.engfracmech.2021.107837)

Cho 2015仅摘要/出版记录核实：PIV测得群集诱导射流与对称涡结构；这些测量不能使本项目粒子速度自动成为流体速度。全文建立时程和全部实验边界未核实，不采用其中峰值作首步加载。[摘要](https://pubmed.ncbi.nlm.nih.gov/26025507/)

Zhang 2022仅出版社摘要核实：FEM/PIV比较入射角对流速、表面剪切及溶栓的影响。未将含药物/微泡因素的溶栓率映射成本材料断裂阈值。[出版社](https://aapm.onlinelibrary.wiley.com/doi/10.1002/mp.15874)

Karmakar 2025本地全文已读：第1节及界面建模讨论本体凝聚破裂与界面脱粘竞争。保持既有LPS变分软化，不移植其NOSB/Ogden及界面标定参数。[论文](https://doi.org/10.1002/cnm.70118)

实现参考官方[Optuna持久化](https://optuna.readthedocs.io/en/stable/tutorial/20_recipes/001_rdb.html)与[BoTorchSampler](https://optuna-integration.readthedocs.io/en/stable/reference/generated/optuna_integration.BoTorchSampler.html)，实际使用Optuna 5.0.0的set_constraint接口。依赖版本以本轮安装清单为准。
