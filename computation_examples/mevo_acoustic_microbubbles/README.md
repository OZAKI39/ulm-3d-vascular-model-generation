# MEVO acoustic–microbubble simulator

## Stage 5A — experimental constrained motion

**FAILED_WITH_EVIDENCE / MOTION_SOLVER_NUMERICAL_GATE_FAILED.** 本轮只开发平移数值算法；原 `lubrication_coupled`、Stage 4L 物理与归档结果保持冻结。新 `constrained_lubrication_v1` 尚未通过连续密集微泡算例，不能用于 production。5 ms 及更长计算未启动。

- [预注册计划](docs/STAGE5A_MOTION_SOLVER_PLAN.md)
- [开发报告与失败原因](docs/STAGE5A_MOTION_SOLVER_DEVELOPMENT.md)
- [图表与机器可读证据](results/stage5/motion_solver/OPEN_RESULTS.html)
- [复核入口和执行范围](scripts/stage5a/README.md)


## Stage 4 — Y-bifurcation branch steering

Stage 4 在 `stage4-bifurcation-control` 分支继续开发，冻结的 Stage 3 标签为 `stage3-minimal-validated`。十个配对种子、40 个纯初级力算例及统计已完成；完整阶段验收尚未生成。

- [2026-10-03 运行记录、已完成结果与次级力适用性问题](docs/STAGE4_RUNTIME_NOTE_20261003.md)
- [进度查询、结果身份和按顺序接续的方法](docs/STAGE4_QUICKSTART.md)
- [预先声明的执行范围与顺序](docs/STAGE4_PLAN.md)

精确的实时阶段以服务器 `results/stage4/CURRENT_EXECUTION_STATE.json` 和 supervisor 状态为准，本地结果为同步快照。条件分支偏向、实际目标到达量和观察截止时未出流比例须同时报告。

## Stage 3 — minimal validated simulator

状态：`MINIMAL_SIMULATOR_VALIDATED_WITH_SCOPE`。四个圆管场景已实跑；73 项生产／回归测试、2 项独立归档测试及 Ruff 检查通过。

本阶段生产链为 **PZT 控制 → 缓存 k-Wave 基底 → APECSS 自适应响应表 → primary / classical secondary Bjerknes → 准稳态 Stokes + FlowField → 几何接触轨迹**。Zhang 已逻辑归档，不再作为下游门禁；生产追踪器不依赖历史力参考模型。

当前状态和定量证据以 [Stage 3 报告](docs/STAGE3_VALIDATION_REPORT.md)及 [PHYSICS_ACCEPTANCE.json](PHYSICS_ACCEPTANCE.json)为准。数值验收与实验物理精度分开记录：历史力敏感性显示明显的输运近似误差，近接触仍只作工程趋势解释。

- [运行方法](docs/STAGE3_QUICKSTART.md)
- [三维运动动画：C3 主动画与 C1–C4 同步对照](docs/STAGE3_MOTION_VISUALIZATION.md)
- [History-force sensitivity：模型、三组 SonoVue 条件、误差与图](docs/HISTORY_FORCE_SENSITIVITY.md)
- [Methods](docs/METHODS_MODEL_DESCRIPTION.md)、[模型映射](docs/MODEL_MAP.md)、[适用范围](docs/KNOWN_LIMITATIONS.md)
- [Zhang 历史证据](validation/deprecated/zhang2011/README.md)、[未来物理](docs/FUTURE_PHYSICS.md)

以下 Stage 1/2/2.5 内容保留为历史基线；其中未运行或受旧门禁阻塞的描述仅适用于当时阶段。

## Stage 2.5 status (2026-10-02)

Branch `stage2.5-pair-physics` adds an independent Zhang forensic reference,
248 matched pair comparisons, measured validity maps, an adaptive APECSS table,
evidence gates and explicit CPU/GPU fallback. The focused numerical response domain
has 790 valid cells and passes 256 independent holdouts at production tolerance.
**Zhang Fig.3 remains unresolved; Lazarus, validated cylinder and sweep were not run.**
See [Stage 2.5 report](docs/STAGE2_5_VALIDATION_REPORT.md) and
[guarded entry points](docs/STAGE2_5_QUICKSTART.md).
Stage1/2 results and their branches/tags remain unchanged.

## Stage2 development status (2026-10-01)

The existing Stage1 demo below is retained. New nonlinear APECSS response tables,
reference validations, strict validity policies and CPU/GPU benchmarks are documented in
[STAGE2_VALIDATION_REPORT](docs/STAGE2_VALIDATION_REPORT.md) and
[STAGE2_QUICKSTART](docs/STAGE2_QUICKSTART.md).
**Stage2 is partially validated:** Zhang2011's distance-induced force sign reversal was not reproduced.
Lazarus clustering and the validated cylinder/frequency sweep remain gated and were not run.
Old linear-surrogate/regularized results remain **demo only**. They are not new SonoVue physics results.


一个可继续扩展的 Python 科研工程：**可控声源 → 3D 线性声场 → 微泡声致平移与空间聚集**。
本目录是现有 BraVa 工作区中的独立声学 demo，不替换或修改 BraVa CFD、正式1565条群体结果或动画。
当前不模拟完整脑血管、颅骨或压电陶瓷内部机电过程。

## 已实现的模型

- 25.6 × 12.8 × 12.8 mm Cartesian 域；默认256 × 128 × 128、0.1 mm间距；外置 PML。
- 轴向 x、长20 mm、半径1.5 mm的圆柱血管。声学介质均匀，血管不引入声阻抗界面。
- 四个独立的3D disc声源，位置、法向、直径、频率、幅值、相位、开关均可配置。
- `kWaveArray` 带限几何源权重，`kspaceFirstOrder` 正式统一 API；GPU实际后端及证据见下面说明。
- CW平滑4周期 ramp，等待两倍域对角线传播时间，末4周期 ROI 录制；每周期取整数时间步，CFL不超过0.25。
- 复声压基底及复梯度的 HDF5 保存、加载、线性组合、有界三线性插值。
- 固定随机种子；单分散、截断 lognormal 或用户给定半径样本。默认直径2.5 µm，**不是完整真实 SonoVue 粒径分布**。
- 线性径向响应替代模型、primary Bjerknes、classical far-field secondary Bjerknes。
- 准稳态 Stokes 阻力、RK2/RK4慢时间积分、位移限制子步、无重叠及圆柱边界投影。
- 每个接受子步的三维位置、速度、力、声压、体积响应、聚类标签、pair-distance histogram；接触事件账本。
- ZeroFlow、Poiseuille，以及实际可用并用 synthetic VTU 验证的静态 SimVascular VTK速度适配器。
- 聚类统计、300 dpi PNG/PDF、ParaView VTR/VTP/PVD。

## 环境与安装

服务器实测 Ubuntu24.04、AMD7800X3D、RTX4090（24564 MiB），driver595.84、CUDA toolkit13.2。
**主机61 GiB RAM不等于容器额度**：本实例内存限额约42.06 GiB、CPU quota7.68核。
项目虚拟环境使用Python3.12.3，不改系统driver或其他项目环境。详细机器信息在 `docs/SERVER_DOCTOR.json`。

```bash
cd /workspace/mevo_acoustic_microbubbles
uv venv --python /usr/bin/python3 .venv
uv pip install --python .venv/bin/python -e '.[dev,gpu]'
# 本次验证版本；已将editable路径规范为当前工程目录
uv pip install --python .venv/bin/python -r requirements.lock.txt
.venv/bin/python -m mevo_mb doctor --output docs/SERVER_DOCTOR.json
```

本实例的 `python` 不在非交互SSH的PATH；请使用 `.venv/bin/python`。
SSH别名 `vast4090` 已配置用户的密钥，指向 `root@50.115.148.16:4159`。本项目不需要新增web服务或开放端口。

## CPU / GPU 验证

```bash
.venv/bin/pytest -q -m 'not gpu'
.venv/bin/python -m mevo_mb acoustic smoke --config configs/smoke.yaml --device cpu
.venv/bin/python -m mevo_mb acoustic smoke --config configs/smoke.yaml --device gpu
.venv/bin/pytest -q -m gpu
.venv/bin/ruff check src tests scripts
```

已实际尝试优先的 `backend="cpp", device="gpu"`：当前binary缺少旧HDF5 ABI，失败原文保留在首次GPU smoke的 `CPP_FAILURE.json`。
安装版0.6.2的C++统一分支还不支持末段起始录制。因此正式demo使用允许的 **`backend="python", device="gpu"`（CuPy）**，不会静默退回CPU。
详见 `docs/BACKEND_DECISION.md`。GPU telemetry记录0.5 s间隔的整设备显存样本，是采样峰值，不能冒称私有分配精确上界。

## 一条命令与模块化CLI

```bash
# 从声场到两种演示；每次创建独立结果目录
.venv/bin/python -m mevo_mb demo
# 已有同身份且叠加验证PASS的声场时复用；不重复求解k-Wave
.venv/bin/python -m mevo_mb demo --reuse-acoustics

.venv/bin/python -m mevo_mb acoustic precompute --config configs/demo_migration.yaml
.venv/bin/python -m mevo_mb acoustic validate-superposition \
  --config configs/demo_migration.yaml --basis results/<basis_run>/acoustic/basis.h5
.venv/bin/python -m mevo_mb bubbles run --config configs/demo_migration.yaml \
  --basis results/<2MHz_basis_run>/acoustic/basis.h5 --pattern baseline
.venv/bin/python -m mevo_mb bubbles run --config configs/demo_migration.yaml \
  --basis results/<2MHz_basis_run>/acoustic/basis.h5 --pattern shifted
.venv/bin/python -m mevo_mb bubbles run --config configs/demo_aggregation.yaml \
  --basis results/<1MHz_basis_run>/acoustic/basis.h5
.venv/bin/python -m mevo_mb visualize --run results/<particle_run>
```

远端长任务使用独立 supervisor program；不修改已有 BraVa program。`scripts/worker.py` 和 `scripts/submit_stage.py`
支持有限声场阶段，失败不自动重试。完整demo同样由独立 `mevo_mb_demo` 管理。本次实际入口和关闭状态写入results内的JSON。
配置支持 `extends`；运行目录保存展开配置和实际控制向量。改变同频率幅值/相位/开关只做基底组合，改变频率或几何需另算基底。
提供2 MHz迁移和1 MHz聚集；没有一次性计算1/1.5/2/2.5/3 MHz全扫频。

## 复数约定、力与单位

全项目采用SI：m、s、Hz、Pa、Pa·s、kg/m³；输出字段后缀明确单位。

`p(t) = Re[P exp(+iωt)]`；控制量 `u_i=A_i exp(iφ_i)`。基底 `H_i` 存储 Pa/source-Pa，故 `P=Σu_i H_i`。
正弦基频用绝对采样时刻作coherent Fourier提取；不会把录制窗口开始当作新的零相位。
梯度对复数实虚部等价地作二阶有限差分，保留空间离散误差；不对物理结果平滑。
`BonA`非空时配置与基底类均拒绝简单线性叠加。

线性替代模型来自径向小扰动振子：

`ξ'' + 2ζω0 ξ' + ω0² ξ = -p/(ρR0)`

`V_hat = -4πR0 P / [ρ(ω0² - ω² + 2iζω0ω)]`，单位m³。

默认 `f0=2.5 MHz, ζ=0.15` **仅为demo假设，不是已验证的SonoVue壳参数**。
默认迁移声源为每元2 kPa，使预测径向幅度保持在小扰动量级；该选择依据线性模型适用性，不是根据轨迹外观调参。
聚集场按初始泡云局部声压中位数75 kPa作线性源标定，源幅值与局部PNP分别报告。
75 kPa下模型不自动具有小振幅或真实SonoVue有效性；每次run报告 `maximum_linear_radial_amplitude_over_R0`。

`F_PB = -0.5 Re[V_hat conj(∇P)]`，m³·Pa/m=N。

`F_SB_ij = ρ/(4πd²) · 0.5 Re[(iωV_i) conj(iωV_j)] · (x_j-x_i)/d`。

同相相同泡吸引，反相排斥；测试直接检查符号、绝对量纲数值、作用反作用及二次幅值缩放。
分母按接触距离/最小距离与正则化长度截断，保留显式cutoff。它是单频、未相互修正径向响应的远场近似，近接触数值正则化不提供近场物理精度。

平移为 `v_free = u_flow + (F_PB+F_SB)/(6πμR)`。边界及接触后的 `total_velocity_m_s` 另有约束修正，
同时保存 `free_total_velocity_m_s`。每个保存状态的速度是**该状态瞬时求速**；RK位移不能直接当作该速度乘前一步长。
固定名义步长25 µs、RK2；为避免跨泡/跨壁，内部步的位移限制为0.2R，全部接受子步都有记录。
几何安全不等于时间步收敛；dt/2对照保留相同物理参数，并报告实际子步数量与最小步长。

## 聚集与边界

聚集定义为空间图连通分量，阈值在YAML；同时报告平均/中位最近邻距离、pair-distance分布、cluster count、最大cluster大小/比例、near-contact比例。
阈值是操作性统计定义，不是物理吸附或融合判据。半径静态保存一次，通过所有时间帧共享，不发生合并、消泡或半径守恒重分配。

壁面处理把球心投影到 `R_vessel-R_bubble` 内并去除外向法向速度；两端也采用显式有限圆柱端盖约束，短时演示不宣称全域出口输运。
接触采用确定性迭代的阻力迁移率加权投影；未收敛时失败，不接受重叠状态。
接触速度通过同一迁移率度量下的整体非负最小二乘投影，同时满足接触对的非穿透约束。
**这些不是润滑、黏附、真实有限面壁面流体动力学或聚并模型**，也没有移植/篡改原BraVa的近壁模型。

迁移演示每种相位150泡，持续10 ms。聚集演示32泡，持续20 ms，初始泡云半径40 µm。
最初100泡开发试跑在密集接触后计算开销显著增大，于1.941884 ms主动中断并保留7935个状态。
32泡用于限定本阶段集成验证规模，没有据此修改力公式、压力、观察时间或接触阈值；不能把它当作100泡结果。
另保留首次接触求解失败和重复调度中断的证据，不将这些部分结果计为成功。

## 数值验证与结果阅读

1. 单源基底与四源直接激励：复压力relative L2、峰值归一化max error、显著幅值区pointwise error、phase RMS/max。
   近压力节点的相位不适定，因此显式报告掩码阈值与排除点数；全ROI复压力L2不删点。
2. 0.1 mm与0.08 mm网格，物理域/源几何/PML厚度相同；在相同血管坐标比较压力、梯度、primary force。
3. 25 µs与12.5 µs粒子步长，比较最终同ID位置、平均位移、最近邻距离和cluster指标。
4. 以上是所测两档差异；即使很小，也不宣称已完成所有尺度的网格无关或真实生理验证。

**实测限制：聚集单泡轨迹尚未证明收敛。** 20 ms的名义dt/2对照末位置RMS差5.72 µm，
平均位移差6.13%；最近邻趋于接触直径，不能仅凭该统计量稳定判断精度。
迁移的平均横向位移在两档声场间相差7.96%，也须继续检查更细网格。
迁移名义dt/2两次均为800内部子步；`scripts/effective_step_check.py`额外检查更小有效步长，
其中聚集补充对照明确只覆盖0–2 ms形成阶段，不替代20 ms密集接触的收敛研究。

`results/<run>/config.yaml`、`metadata.json`、`acoustic/`、`particles/`、`vtk/`、`figures/`、`metrics.csv`保存每次独立结果。
`particles/trajectories.h5`含所有接受子步，`events.jsonl`含接触/壁面事件。
最终粒子运行同时保存 `SOURCE_IDENTITY.json`、`logs/source_snapshot.zip` 和实际使用的
`acoustic/combined_field.h5`。因此细网格对照重绘仍使用原始细网格场。
早期声场开发阶段只记录当时Git提交及运行参数，存在未提交开发修改；不得把该提交号单独当作完整源码快照。
ParaView打开 `vtk/particles.pvd`，按时间播放，使用 `radius_m` glyph球体；叠加 `vtk/acoustic_field.vtr` 查看复声场幅度/相位。
300 dpi PNG和PDF来自原始数值，图文入口为项目 `OPEN_RESULTS.html`，最终定量结果见 `docs/VALIDATION_REPORT.md`。

本次进一步检查的可复现入口（长任务仍应由supervisor运行）：

```bash
.venv/bin/python scripts/convergence_particles.py  # dt/2和已生成细声场，最多2个CPU worker
.venv/bin/python scripts/check_convergence.py      # 比较已关闭运行，输出定量报告
.venv/bin/python scripts/secondary_control.py      # 相同泡云和控制量，关闭secondary force
.venv/bin/python scripts/review_artifacts.py       # 全帧完整性、pair-distance CDF、几何
.venv/bin/python scripts/build_report.py           # HTML/Markdown入口
```

细声场的有限任务由 `scripts/submit_stage.py --help` 所列参数提交；网格物理尺寸不变，步长从0.1 mm到0.08 mm。
含外置PML的实际FFT网格分别为288×160×160和360×200×200，可从每次metadata的物理网格和PML宽度推得。
关闭secondary的控制试验不改变primary force、初始ID或源控制；其结果用于区分泡间作用与primary力本身的汇聚。

## SimVascular与未来扩展

`FlowField.velocity(x_m,t_s)` 隔离上层积分器与流场来源。`SimVascularFlowField`支持VTU/VTK点数据速度，坐标与速度单位缩放显式配置。
域外查询抛错，不用最近邻、零填充或投影掩盖。未来加入时间插值、网格定位缓存及真正患者几何边界。
当前圆柱边界不能直接用于患者血管；更换流场之外还必须更换几何/边界提供器并独立验证。

`MarmottantBubbleModel`与无默认值的参数dataclass已定义，但响应调用明确抛 `NotImplementedError`。
后续校验Keller–Miksis/Marmottant或接入APECSS；预计算 `(R0,f,|P|)→V_hat` 表，而非在每个slow step重解MHz ODE。
`PairInteractionBackend`目前实现NumPy向量化，N<500；预留CuPy、neighbor/cell list扩展，没有不必要的scikit-learn依赖。

未实现：piezoelectric FEM、非线性声传播、声流、颅骨、FSI、气泡散射反馈、两向血流耦合、壳破裂、聚并、变形和高Re拖曳修正。
当前是单向耦合。低Re与小振幅近似的有效性须针对后续真实工况验证。

来源与推导说明：[PHYSICS.md](docs/PHYSICS.md)。所有经验/模型参数均在YAML中；当前不把任何demo参数称为患者或商品材料标定。
