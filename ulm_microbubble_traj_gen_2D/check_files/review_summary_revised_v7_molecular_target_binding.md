# Revised v7 靶区与分子键合实现 Review Summary

日期：2026-07-15

## 1. 本次实现目标

本次实现把 `references/Revised v7.md` 中的真实壁面靶区与确定性分子键合模型接入现有二维/2.5D 微泡连续灌注流程，核心目标如下：

- 靶区只能由物理坐标多边形或带物理坐标轴的 NPZ 掩膜显式给出，不从速度、剪切力或历史键合结果反推。
- 靶区只存在于封闭固体壁面，入口和出口开口不能被误当成可成键壁面。
- 使用连续的期望键数和总切向伸长描述平均场键合，不引入随机离散键、二元“已黏附”状态或失活槽位复用。
- 按局部球冠/反应盘几何计算靶点阳性面积，并在靶区边缘保持连续变化；面积单位显式从 `um^2` 转换为 `m^2` 后再参与表面分子数和成键率计算。
- Bell 滑脱键产生的切向力、法向力和力矩与碰撞载荷一起进入已有迁移率矩阵，使键合通过同一物理通道改变平移和旋转。
- 键状态只在时间步和有限尺寸几何约束均被接受后提交；壁面约束缩短运动时，总伸长使用实际接受的切向位移和旋转推进。
- 在没有实验标定二维成键系数时，提供空血管、无键力接触 pilot 和基于 `Da_on` 的敏感性场景，不伪造“真实”参数。
- 保存靶场、轨迹分子诊断、参数来源和数值有效性信息，同时保持旧配置、旧输出读取和现有可视化路径兼容。

## 2. 修改了哪些文件，每个文件为什么改

### 入口、配置与编排

| 文件 | 状态 | 修改原因 |
| --- | --- | --- |
| `generate_microbubble_trajectories.py` | 修改 | 在生成完成后打印靶场 NPZ 和 contact-pilot 报告路径。 |
| `configs/physics_flow_config.yaml` | 修改 | 新增默认关闭的 `molecular_target`、`molecular_binding` 和 `binding_scenario_sweep` 配置；所有带量纲参数在键名中写明单位。 |
| `utils/config.py` | 修改 | 新增三组配置数据结构、严格解析和交叉验证；拒绝没有物理坐标的 ROI、量纲不匹配的速率、无接触却要求生成二维成键系数等配置。 |
| `utils/runner.py` | 修改 | 在接受 CFD 场后构建靶场、把分子配置传入输运、按需运行 pilot，并保存解析后的配置、靶场和元数据。 |
| `utils/particle_advection.py` | 修改 | 扩展粒子推进调用链，将靶场和分子键合配置显式传递到 mobility 连续灌注实现。 |

### 靶区、键合与无实验敏感性流程

| 文件 | 状态 | 修改原因 |
| --- | --- | --- |
| `utils/molecular_target_field.py` | 新增 | 构建真实封闭壁面靶区；支持物理多边形和坐标化 NPZ 掩膜；保存壁面坐标/法向；以解析区间交集和圆弓积分计算连续靶点阳性反应面积。 |
| `utils/molecular_binding.py` | 新增 | 实现反应盘几何、二维平均场成键率、张力型 Hooke 键、Bell 滑脱率、分子力/力矩，以及指数 Euler/Heun 状态推进与容量/溢出保护。 |
| `utils/molecular_contact_pilot.py` | 新增 | 从无键力轨迹统计逐泡接触时长、接触事件、反应面积和表面滑移，并由观测接触时长生成维数正确的 `Da_on` 敏感性定义。 |
| `utils/molecular_pilot_runner.py` | 新增 | 将配置、轨迹、靶场和 pilot 分析连接起来，为每个捕获距离比例写出独立 YAML 报告；无有效接触时明确报告不可辨识。 |

### 粒子输运与输出结构

| 文件 | 状态 | 修改原因 |
| --- | --- | --- |
| `utils/particle_mobility_transport.py` | 修改 | 在纯 RHS 中计算局部靶区几何和键载荷；将碰撞力与键力/键力矩合并后只通过一次现有迁移率矩阵；输出分子诊断。 |
| `utils/particle_perfusion_transport.py` | 修改 | 为每个永久 ID 增加不复用的 `n/m` 状态；实现 Euler/Heun 预测与接受提交；用几何约束后的实际运动推进伸长；记录容量限幅、Bell 限幅、平均场和局部平面近似诊断。 |
| `utils/types.py` | 修改 | 为轨迹结果增加可选分子状态、力、力矩、速率、反应面积及最终 registry 状态字段；绑定关闭时保持 `None`。 |
| `utils/field_io.py` | 修改 | 保存独立 `molecular_target_field.npz` 和轨迹中的分子数组；分子诊断使用 `float64`，避免有限 double 值在写出时降精度溢出。 |

### 可复现场景、文档与测试

| 文件 | 状态 | 修改原因 |
| --- | --- | --- |
| `configs/molecular_binding_scenarios/literature_psa_pselectin.yaml` | 新增 | 保存 PSA:P-selectin 文献参数的来源型预设；明确文献一阶前向速率不能直接作为二维 `k_on`，因此默认不可直接运行。 |
| `configs/molecular_binding_scenarios/literature_slex_eselectin.yaml` | 新增 | 保存 sLeX:E-selectin 文献参数的来源型预设，并应用相同量纲和捕获距离警告。 |
| `configs/molecular_binding_scenarios/no_experiment_dimensionless_sweep.yaml` | 新增 | 给出不声称实验标定的 `Da_on`、密度和捕获距离比例敏感性轴；只定义场景，不自动批量运行。 |
| `REPRODUCTION_SPEC.md` | 修改 | 更新当前调用链、全部新配置、公式、单位、假设、NPZ schema、pilot 解释、有效性诊断、复现步骤和测试结果。 |
| `test_files/test_molecular_config.py` | 新增 | 测试默认兼容、合法配置和跨 section 的失败条件。 |
| `test_files/test_molecular_scenario_configs.py` | 新增 | 测试三个场景文件的量纲警告、禁用状态和非标定语义。 |
| `test_files/test_molecular_target_field.py` | 新增 | 测试真实壁面交集、开口排除、物理坐标注册、圆盘全覆盖、半覆盖及靶边左右连续性。 |
| `test_files/test_molecular_binding.py` | 新增 | 测试反应几何、单位换算、成键动力学、张力/Bell 公式、符号、容量保护、纯 RHS、Euler/Heun 和实际位移覆盖。 |
| `test_files/test_molecular_contact_pilot.py` | 新增 | 测试接触事件、矩形时间积分、`Da_on` 换算和无接触不可辨识行为。 |
| `test_files/test_molecular_pilot_runner.py` | 新增 | 测试报告生成、场景笛卡尔积、来源上下文和零/无效法向处理。 |
| `test_files/test_molecular_transport_integration.py` | 新增 | 测试分子载荷进入 mobility、键跨越靶区边缘后的持续解离、永久 ID 状态隔离、NumPy/Numba 一致性、几何约束后的伸长和 NPZ schema。 |

## 3. 核心逻辑变化，按调用链说明

1. `load_config` 读取完整物理配置。默认三个分子 section 均关闭；启用键合时必须同时提供真实靶区、正的物理参数和维数正确的二维成键系数。
2. `run_generation` 沿用原有血管加载、栅格化、PhiFlow 求解和接受场检查。只有 CFD 场被接受后，才从粒子壁面距离/法向场建立分子靶场。
3. 靶场构建器将用户 ROI 与封闭固体壁面求交，并排除入口/出口开口及其保护邻域。多边形保留精确边界；掩膜保留原始物理坐标轴和解析后的源路径。
4. 连续灌注仍从空血管开始，仍使用永久唯一 ID、上游等待队列和不复用数组槽位。新 ID 被接纳时，其期望键数与总伸长都严格初始化为零。
5. 每次 mobility RHS 根据粒子中心、半径、壁面间隙、局部内法向和切向构造反应盘。靶点阳性面积由 ROI 区间与有效封闭壁面区间相交后解析积分得到；完整靶区时严格为 `pi*a^2`，跨靶区边缘时连续变化。
6. 平均场模型由当前 `n/m`、局部靶区面积、固定配体/靶点密度和表面滑移计算成键与解离；已有键不会因瞬时反应面积缩小而被直接裁掉，离开靶区后只停止新生并继续 Bell 解离。
7. 分子切向力、向壁法向力和力矩与同步计算的粒子碰撞载荷相加，一起传给同一个自由空间/近壁混合迁移率算子。CFD 背景场不被二次当成力，也不因键合重新求解。
8. Euler 或 Heun 先产生运动候选，再执行有限尺寸壁面、路径和出口约束。Heun 第二阶段使用已约束的预测几何；最终 `m` 的伸长源取实际接受的切向中心位移与旋转，最后才原子式提交位置、旋转和键状态。
9. 每帧保存连续键状态、反应面积、形成/解离率、分子力和力矩；每个永久 ID 另存最终键状态。启用靶区时额外写出独立靶场 NPZ。
10. 若只启用 no-bond pilot，则正式轨迹不施加键力。pilot 从保存帧中测量真实接触，再按 `k_on = Da_on/(rho_T*t_contact)` 生成敏感性报告；没有正接触时不生成任何伪参数。

## 4. 行为变化：改动前 vs 改动后

| 方面 | 改动前 | 改动后 |
| --- | --- | --- |
| 壁面生物靶区 | 不存在 | 可显式配置物理多边形或坐标化 NPZ，并严格限制为封闭固体壁面。 |
| 微泡状态 | 位置、速度、旋转及碰撞等水动力状态 | 在键合启用时增加连续 `n/m` 分子状态；无二元黏附开关。 |
| 粒子附加载荷 | 主要为粒子碰撞 | 碰撞与分子键力/力矩合并后进入同一迁移率矩阵。 |
| 反应面积 | 不存在 | 按有限尺寸球冠接触几何和靶区边界解析积分，靶边连续。 |
| 离开靶区 | 无相关行为 | 停止形成新键，但已有键保留并按受力 Bell 速率解离。 |
| 时间推进 | 原有同步 Euler/Heun 与几何约束 | 原有结构保留；键状态使用预测/接受两阶段，且伸长与最终接受运动一致。 |
| 参数不足 | 无键合功能 | 可先做无键力接触 pilot 和无量纲敏感性扫描；不会把文献 `s^-1` 误作二维 `k_on`。 |
| 输出 | 原有轨迹/场/元数据 | 键合启用时追加 v5 分子轨迹字段和靶场 NPZ；关闭时仍使用原 v4 schema。 |
| 默认运行 | 无分子键合 | 数值行为保持不变，因为新增功能默认全部关闭。 |

## 5. 新增/修改了哪些测试

新增 7 个分子专项测试模块，共 49 项测试：

- 配置与场景：默认禁用兼容、量纲约束、真实 ROI 要求、场景 provenance 和 no-contact guard。
- 靶场几何：多边形/NPZ 注册、开口排除、完整圆盘、精确半圆盘、靶边连续性和源路径输出。
- 分子动力学：形成率、张力、Bell 指数、力/力矩方向、容量限幅、指数 Euler/Heun、纯函数语义和实际位移伸长源。
- pilot：逐泡接触时长、连续事件、面积/滑移加权、单位换算、零法向和无接触不可辨识。
- 完整输运：永久 ID 状态隔离、同步 mobility 耦合、离开靶区后的键衰减、约束运动、NumPy/Numba 一致性以及保存/读取字段。

除新增专项外，完整 `test_files` 回归共发现并通过 148 项测试，覆盖已有心跳脉动、连续灌注、碰撞、粒子补充、壁面几何、PhiFlow 和 PyVista/VTK 路径。

## 6. 已运行的检查命令及结果

使用环境：`D:\anaconda3\envs\pmp\python.exe`

1. 分子专项测试：

   ```powershell
   D:\anaconda3\envs\pmp\python.exe -m unittest ulm_microbubble_traj_gen.test_files.test_molecular_config ulm_microbubble_traj_gen.test_files.test_molecular_target_field ulm_microbubble_traj_gen.test_files.test_molecular_binding ulm_microbubble_traj_gen.test_files.test_molecular_contact_pilot ulm_microbubble_traj_gen.test_files.test_molecular_pilot_runner ulm_microbubble_traj_gen.test_files.test_molecular_scenario_configs ulm_microbubble_traj_gen.test_files.test_molecular_transport_integration -q
   ```

   结果：`Ran 49 tests`，`OK`。

2. 全量回归：

   ```powershell
   D:\anaconda3\envs\pmp\python.exe -m unittest discover -s ulm_microbubble_traj_gen\test_files -q
   ```

   结果：`Ran 148 tests in 9.361s`，`OK`。

3. Python 编译检查：

   ```powershell
   D:\anaconda3\envs\pmp\python.exe -m compileall -q ulm_microbubble_traj_gen
   ```

   结果：退出码 0。

4. 本次所有改动 Python 文件的 Ruff 检查：

   ```powershell
   D:\anaconda3\envs\pmp\python.exe -m ruff check <本次改动的 Python 文件列表>
   ```

   结果：`All checks passed!`。

5. 默认配置与 CLI 冒烟检查：

   - 默认配置成功加载，`molecular_target.enabled=false`、`molecular_binding.enabled=false`、`binding_scenario_sweep.enabled=false`。
   - `python -m ulm_microbubble_traj_gen.generate_microbubble_trajectories --help` 正常退出，退出码 0。

6. 扩大到整个 `utils` 的 Ruff 检查仍报告一个本次范围外的既有问题：`utils/vessel_rasterizer.py:430` 中 `max_radius` 赋值后未使用（F841）。本次没有改动该文件，也没有为了清除无关告警删除原代码。

未执行完整真实血管数据的长时间 PhiFlow 生产仿真；该运行会消耗显著时间，并且真正的分子运行还必须由用户提供真实疾病 ROI 和有明确依据的二维参数。调用链、数学核心、序列化以及合成端到端输运已由上述测试覆盖。

## 7. 可能的回归风险

1. 当前仍是二维/2.5D、局部单平面壁近似；分叉、强曲率或同时靠近两侧壁面时，真实三维接触斑和旋转耦合可能不同。
2. 确定性平均场在反应区内有效配体或靶点数很少时可能不适用；代码只给工程警告，不会自动切换到随机离散键模型。
3. 未经实验标定的有效二维 `k_on`、靶点密度、捕获距离和疾病 ROI 仍是主要生物学不确定性。场景 YAML 是敏感性定义，不是生产真值。
4. pilot 使用保存帧矩形积分，短于输出间隔的接触可能漏计；应同时做输出步长和内部子步收敛研究。
5. 解析面积依赖局部壁面采样和最近壁面归属；网格过粗时，即使积分在靶边连续，几何本身仍可能欠分辨。
6. 高刚度、高成键率或高 Bell 指数会形成刚性问题。容量限幅或指数上限一旦被触发，必须增加内部子步并做收敛比较，而不能把警告当作普通物理结果。
7. 分子诊断改用 `float64` 增加轨迹 NPZ 体积；长时、大浓度运行需评估 `max_particle_frame_records` 和磁盘空间。
8. 键合是单向耦合：它改变粒子运动，但不会反馈到 CFD、壁面剪切、靶点密度或其他粒子表面的化学状态。

## 8. 需要我人工重点看的 5 个位置

1. `configs/physics_flow_config.yaml` 与 `utils/config.py`：确认计划使用的真实 ROI、配体/靶点密度、捕获距离和二维速率单位；不要直接把论文中的一阶 `s^-1` 前向速率填入二维 `k_on`。
2. `utils/molecular_target_field.py` 中 `build_molecular_target_field` 及反应区间构造：把真实疾病区域叠加到 `molecular_target_field.npz` 上人工确认，尤其检查根入口、末端出口和分叉附近。
3. `utils/molecular_binding.py` 中 `evaluate_mean_field_bonds`、`predict_bond_state_exponential_euler` 和 `accept_bond_state_exponential_heun`：核对所选配体体系的单位、力方向、Bell 参数与平均场假设。
4. `utils/particle_mobility_transport.py` 中 `_evaluate_rhs`：确认碰撞载荷和键载荷只合并一次、壁面切向/内法向约定与预期一致，并核对单向耦合是否满足研究目的。
5. `utils/particle_perfusion_transport.py` 中 `_advance_segment` 及分子 metadata/输出：重点检查几何约束后的实际运动如何推进 `m`，并在首个真实算例中核对容量限幅、Bell 限幅、`a/R_curvature`、低分子数及子步收敛诊断。

## 9. 是否有任何超出需求范围的改动

没有引入超出需求的模型改动：

- 未修改 PhiFlow 求解器、现有 Goldman 近壁迁移率系数、碰撞本构、心跳脉动、入口通量算法或微泡直径分布。
- 未修改 `visualize_microbubble_flow.py`、PyVista/VTK 场可视化或 trame viewer；本次只增加数据产物，分子键可视化仍不在范围内。
- 未实现随机离散键、受体耗竭/扩散、靶区自动推断、三维键几何、黏附对 CFD 的双向反馈或自动场景批量运行。
- 未迁移旧 Taichi 路径；现有 NumPy/Numba 后端继续有效，并已做一致性测试。
- 未安装新依赖；当前环境已有 NumPy、SciPy、Numba、PyYAML 等足以完成实现。
- 除按要求更新 `REPRODUCTION_SPEC.md`、新增测试与场景说明外，没有整理或删除工作区中其他既有代码和注释，也没有处理与本次功能无关的 Ruff 告警。
