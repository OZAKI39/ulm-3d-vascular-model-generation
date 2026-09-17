# Revised v16 连续血管壁实施状态（2026-07-20）

本节优先于本文后续仍提到 `mask_npz` 和栅格固壁的旧操作说明。正式粒子输运现已采用 Revised v16 的连续几何—栅格流场混合方案：

- 血管中心线、真实半径及分叉过渡区先生成连续管腔 $\Omega$；规则网格只是该连续管腔在 CFD 上的采样，不再是粒子壁面几何真值。
- 微泡间隙 $g_R$、壁面最近点、法向、有限半径可达性、扫掠路径接触、入口准入和有向出口穿越都直接查询连续边界。
- 单壁面预测接触仍保留 Revised v15 的动力学公式，但其间隙和法向已全部改为连续 v16 几何；真正的双壁面约束仍会以含永久微泡 ID 的诊断异常停止，不会静默投影或穿墙。
- 旧的栅格壁面构造仅保留为旧测试和旧文件读取兼容路径。`runner.py` 与候选靶区正式调用链始终显式构造并传入 `ContinuousVesselGeometry`，不会回退到栅格粒子壁面。

固定靶区已迁移到：

```text
ulm_microbubble_traj_gen/study_inputs/fixed_target_v16_20260720/continuous_wall_target.npz
```

迁移不是把旧靶区像素投影到新壁面，而是保留原来的物理影响圆、随机种子及随机傅里叶实现，在连续壁面弧长上重新评价同一个随机场。旧文件保持不变，其 SHA-256 为：

```text
56CBC83045FA9283BD1D71C1C6426E5AB4CA57DD747B5D39F74652572348DD7E
```

当前连续几何和新靶区共同绑定的几何 SHA-256 为：

```text
5C8A7AEFC83DCEA7E4B629B477457CEA133647FA35769E9B989FF659F9E551BA
```

新靶区在物理影响圆内的连续壁长阳性比例为 `0.8000500789`。该壁长比例不能与旧栅格掩膜的像素面积比例直接比较。

验证状态：247 个自动化测试全部通过。正式 `1.5 um` CFD 网格的 20 秒试验已完整通过，结果保存在 `results/20260720_172800`。该试验创建 671 个永久微泡，最大同时活动 86 个，记录 592 次有向出口事件；成功越过此前约 `1.2715 s`、`2.3637 s` 和 `6.7683 s` 的失败时刻，且穿墙、有限半径可达域违规、接触细化和真实双壁面失败均为 0。1,452,113 条保存状态中有 24 条触发“靠近对侧壁面的单平面水动力近似”警告；这不是同时接触两面壁，也没有触发接触约束。

本次 20 秒试验的三份接触报告均为 0 个靶区接触微泡。正式轨迹的最小连续壁间隙约为 `0.09219 um`，高于当前三组试验的最大捕获距离 `0.058 um`，所以不能再通过调整成键速率来补救：当前首先缺少几何暴露机会。这是新的物理参数筛选结论，不是 Particle transport 运行失败；进入弱、中、强成键参数比较前，需要另行确认是否扩大捕获距离、改变靶区位置/范围或调整输运工况。

---

最适合你的实施顺序是：固定现有靶区不动，先提高接触统计可靠性，再执行弱、中、强三组成键仿真。现阶段不把仿真靶区称为“真实疾病靶区”。

另外，我检查工作区时发现原来的 `20260717_114349` 目录目前已改名为：

[Pipeline validation](</E:/ULM/hatimb-particle_flow_simulator/ulm_microbubble_traj_gen/results/Pipeline%20validation>)

下面均按当前路径说明。

# 第一阶段：永久冻结当前靶区

## 第一步：建立研究输入目录

在 PowerShell 中执行：

```powershell
New-Item `
  -ItemType Directory `
  -Force `
  "ulm_microbubble_traj_gen\study_inputs\fixed_target_20260717"
```

复制已经验证成功的靶区：

```powershell
Copy-Item `
  "ulm_microbubble_traj_gen\results\20260716_153659\selected_molecular_target_mask_fraction_0p8.npz" `
  "ulm_microbubble_traj_gen\study_inputs\fixed_target_20260717\fixed_target_mask.npz"
```

如果原目录中有对应的自动选择报告，也一并复制：

```powershell
Copy-Item `
  "ulm_microbubble_traj_gen\results\20260716_153659\automatic_molecular_target_selection.json" `
  "ulm_microbubble_traj_gen\study_inputs\fixed_target_20260717\automatic_molecular_target_selection.json"
```

记录文件哈希：

```powershell
Get-FileHash `
  "ulm_microbubble_traj_gen\study_inputs\fixed_target_20260717\fixed_target_mask.npz"
```

保存这个哈希值。以后所有固定靶区研究都必须使用同一个 NPZ 和相同哈希。这样可以证明你没有为了得到成键结果而移动或扩大靶区。

需要注意：如果 `fraction_0p8` 是在看到旧的失败结果后，为了获得接触而从较低比例提高到 0.8，那么它应该被称为：

> 合成高表达靶区情景

不能称为真实疾病表达水平。

# 第二阶段：延长无成键接触试验

目前 5 秒试验只有 169 个微泡，其中只有一个微泡接触靶区。因此先增加观测时间，但不修改靶区。

## 第二步：复制一份新的接触试验配置

```powershell
Copy-Item `
  "ulm_microbubble_traj_gen\configs\molecular_contact_pilot.yaml" `
  "ulm_microbubble_traj_gen\configs\molecular_contact_fixed_target_20s.yaml"
```

打开新文件：

[接触试验配置](/E:/ULM/hatimb-particle_flow_simulator/ulm_microbubble_traj_gen/configs/molecular_contact_fixed_target_20s.yaml)

只修改下面几处。

### 模拟时间改为 20 秒

```yaml
particles:
  inlet_number_concentration_mb_per_ml: 20000000.0
  n_steps: 20000
  dt_s: 0.001
  max_particle_frame_records: 5000000
```

因为：

$$
T=n_{\mathrm{steps}}\Delta t
  =20000\times0.001
  =20\ \mathrm{s}
$$

### 靶区改为冻结后的文件

```yaml
molecular_target:
  enabled: true
  region_mode: mask_npz
  mask_npz_path: 'E:/ULM/hatimb-particle_flow_simulator/ulm_microbubble_traj_gen/study_inputs/fixed_target_20260717/fixed_target_mask.npz'
  mask_array_key: target_mask
  x_coordinates_key: x_um
  z_coordinates_key: z_um
  target_density_molecules_per_m2: 0.0
```

### 保持成键关闭、接触分析开启

```yaml
molecular_binding:
  enabled: false
```

```yaml
binding_scenario_sweep:
  enabled: true
```

其余内容保持不变。

## 第三步：运行 20 秒试验

```powershell
D:\anaconda3\envs\pmp\python.exe `
  ulm_microbubble_traj_gen\generate_microbubble_trajectories.py `
  --config "ulm_microbubble_traj_gen\configs\molecular_contact_fixed_target_20s.yaml" `
  --skip-render
```

参数筛选阶段使用 `--skip-render` 只跳过 PyVista/VTK 可视化，数值场、粒子轨迹、
元数据和接触报告仍会保存。需要最终 CFD/壁面剪切场可视化时，去掉该选项重新运行。

当前配置的平均注入率约为：

$$
33.43\ \mathrm{MB/s}
$$

所以20秒大约会创建：

$$
33.43\times20\approx669
$$

个微泡。

按照历史运行中观察到的约 $0.59\%$ 接触率粗略估计，期望值只有约4个接触微泡。因此20秒主要用于检查数值稳定性和接触趋势，不能把“至少10个接触微泡”作为这一步的硬性通过条件。

## 第四步：检查接触报告

在新时间戳文件夹中打开：

```text
molecular_contact_pilot_01_capture_ratio_1p5.yaml
```

重点检查：

```yaml
pilot:
  n_unique_bubbles:
  n_contacting_bubbles:
  contacting_bubble_fraction:
  total_contact_time_s:
  median_positive_bubble_contact_time_s:
```

以及：

```yaml
contact_events:
```

先检查数值正确性；以下项目必须通过：

- 仿真完整运行到20秒，不能因离散网格可达性判断中止；
- `active_outside_lumen_violations: 0`；
- `accepted_negative_gap_count: 0`；
- `contact_nonzero_velocity_zero_progress_count: 0`；
- 出口终止数量和各出口计数一致；
- `discrete_accessibility_disagreement_records` 是离散最近网格点与连续壁面几何之间的审计计数，不是运行失败条件。

数值正确性通过后，再检查统计质量：

- 20秒试验允许只有0到少量接触；若要据此选择成键参数，`n_contacting_bubbles` 应至少达到10，最好达到20；
- 接触不能全部由同一个微泡贡献；
- 不应只有一个长接触事件一直持续到最后一帧；
- 前10秒和后10秒的接触率不应相差数倍。

后四项是统计质量标准，不是20秒数值回归的物理定律。

## 第五步：如果接触样本仍少，增加到50秒

复制20秒配置：

```powershell
Copy-Item `
  "ulm_microbubble_traj_gen\configs\molecular_contact_fixed_target_20s.yaml" `
  "ulm_microbubble_traj_gen\configs\molecular_contact_fixed_target_50s.yaml"
```

修改为：

```yaml
particles:
  n_steps: 50000
  dt_s: 0.001
```

50秒预计注入约：

$$
33.43\times50\approx1672
$$

个微泡。按照 $0.59\%$ 的历史粗略接触率，接触微泡的期望值约为10，但达到至少10个的概率只有约54%。若希望达到至少10个接触微泡的概率约为90%，建议连续运行约72秒；约95%则建议约80秒。

如果20秒运行中接触数为0，不要只靠延长时间来掩盖问题。应先检查目标区域附近的最小间隙、近壁速度、捕获距离以及目标掩膜是否与实际可达壁面重合。

不要通过重复运行完全相同的5秒配置增加样本，因为当前注入序列是确定性的，相同配置通常会重复相同结果。应延长一次连续仿真。

# 第三阶段：从新的接触报告选择成键情景

当前配置在运输开始前固定声明：

$$
t_{\mathrm{ref}}
=
\texttt{da\_on\_reference\_time\_s}
=1.0\ \mathrm{s}
$$

场景中的 $k_{\mathrm{on}}$ 由预先声明的 $Da_{\mathrm{on}}$、靶点密度和 $t_{\mathrm{ref}}$ 推导，不应再由本次运行观察到的接触时间反推。长时间接触报告中的 `median_positive_bubble_contact_time_s` 只用于评估采样质量和接触动力学，不会改变同一配置的 $k_{\mathrm{on}}$。

## 第六步：先只选择三种成键强度

使用新生成的：

```text
molecular_contact_pilot_01_capture_ratio_1p5.yaml
```

选择以下三行：

| 情景     | 报告中的`scenario_index` | $Da_{\mathrm{on}}$ |
| -------- | -------------------------: | -------------------: |
| 弱成键   |                         19 |                  0.1 |
| 中等成键 |                         31 |                    1 |
| 强成键   |                         43 |                   10 |

这三个情景保持以下参数相同：

- 捕获距离与静息长度之比：1.5；
- 靶点密度：$270\ \mathrm{molecule/\mu m^2}$；
- 配体密度：$100\ \mathrm{molecule/\mu m^2}$。

只改变无量纲成键强度。

暂时不要直接运行全部180种组合。

# 第四阶段：建立三个正式成键配置

## 第七步：复制三个配置文件

```powershell
Copy-Item `
  "ulm_microbubble_traj_gen\configs\molecular_contact_fixed_target_20s.yaml" `
  "ulm_microbubble_traj_gen\configs\binding_fixed_target_weak.yaml"

Copy-Item `
  "ulm_microbubble_traj_gen\configs\molecular_contact_fixed_target_20s.yaml" `
  "ulm_microbubble_traj_gen\configs\binding_fixed_target_medium.yaml"

Copy-Item `
  "ulm_microbubble_traj_gen\configs\molecular_contact_fixed_target_20s.yaml" `
  "ulm_microbubble_traj_gen\configs\binding_fixed_target_strong.yaml"
```

如果最终采用50秒接触试验，则从50秒版本复制。

## 第八步：把报告中的情景参数填入配置

以中等情景为例，从新报告的 `scenario_index: 31` 中复制以下数值：

```yaml
molecular_target:
  enabled: true
  region_mode: mask_npz
  mask_npz_path: 'E:/ULM/hatimb-particle_flow_simulator/ulm_microbubble_traj_gen/study_inputs/fixed_target_20260717/fixed_target_mask.npz'
  mask_array_key: target_mask
  x_coordinates_key: x_um
  z_coordinates_key: z_um
  target_density_molecules_per_m2: 270000000000000.0
```

```yaml
molecular_binding:
  enabled: true
  model: deterministic_mean_field

  ligand_density_molecules_per_m2: 100000000000000.0
  capture_distance_um: 0.0435
  rest_length_um: 0.029

  association_rate_m2_per_molecule_s: 从新报告的scenario_index_31复制
  zero_force_dissociation_rate_s: 0.06
  bond_stiffness_pn_per_um: 10000.0
  reactive_compliance_nm: 0.039
  temperature_k: 310.0

  mean_field_warning_count: 10.0
  bell_exponent_limit: 80.0
```

正式成键时关闭接触参数生成器：

```yaml
binding_scenario_sweep:
  enabled: false
  da_on_levels: []
  capture_distance_to_rest_length_ratios: []
  target_density_molecules_per_um2_levels: []
  ligand_density_molecules_per_um2_levels: []
```

弱、中、强三个配置的主要差别是：

```yaml
association_rate_m2_per_molecule_s:
```

分别从新报告的情景19、31、43复制。

## 第九步：分别运行三个正式情景

```powershell
D:\anaconda3\envs\pmp\python.exe `
  ulm_microbubble_traj_gen\generate_microbubble_trajectories.py `
  --config "ulm_microbubble_traj_gen\configs\binding_fixed_target_weak.yaml"
```

```powershell
D:\anaconda3\envs\pmp\python.exe `
  ulm_microbubble_traj_gen\generate_microbubble_trajectories.py `
  --config "ulm_microbubble_traj_gen\configs\binding_fixed_target_medium.yaml"
```

```powershell
D:\anaconda3\envs\pmp\python.exe `
  ulm_microbubble_traj_gen\generate_microbubble_trajectories.py `
  --config "ulm_microbubble_traj_gen\configs\binding_fixed_target_strong.yaml"
```

除了成键参数之外，三个配置的以下内容必须完全相同：

- 靶区文件；
- 微泡浓度；
- 仿真时长；
- 心跳参数；
- CFD参数；
- 粒子积分子步；
- 碰撞参数；
- 微泡直径范围。

# 第五阶段：比较结果

每个正式结果的 `microbubble_field_trajectories.npz` 中都有以下数据：

- `record_bond_count_expected`：平均场意义下的期望键数量；
- `record_bond_force_xz_pn`：分子键总作用力；
- `record_bond_torque_pn_um`：分子键力矩；
- `record_single_bond_tension_pn`：单键张力；
- `record_bond_formation_rate_bonds_s`：成键形成速率；
- `record_bond_dissociation_rate_s_inv`：解离速率；
- `record_target_reaction_area_um2`：有效反应面积；
- `registry_final_bond_count_expected`：每个永久微泡ID结束时的期望键数。

至少比较：

1. 最大期望键数；
2. 期望键数达到或超过1的微泡比例；
3. 第一次达到1个期望键的时间；
4. 累计成键时间；
5. 微泡速度相对无成键试验下降多少；
6. 最终滞留在靶区附近的微泡数量；
7. 正常流出血管的微泡数量；
8. 是否出现容量限幅、Bell指数上限或积分收敛警告。

判断方法：

- 弱、中、强三个情景都不成键：接触机会或参数范围过弱；
- 只有强情景明显成键：系统对动力学参数高度敏感；
- 中等和强情景结果接近：可能进入饱和区；
- 三个情景结论一致：结论对成键强度较稳健。

# 第六阶段：再决定是否扩大敏感性矩阵

完成弱、中、强三组后，再增加以下六组：

| 研究内容        | 报告文件                   | 情景索引 |
| --------------- | -------------------------- | -------: |
| 捕获距离比例1.0 | `capture_ratio_1.yaml`   |       31 |
| 捕获距离比例2.0 | `capture_ratio_2.yaml`   |       31 |
| 低靶点密度      | `capture_ratio_1p5.yaml` |       28 |
| 高靶点密度      | `capture_ratio_1p5.yaml` |       34 |
| 低配体密度      | `capture_ratio_1p5.yaml` |       30 |
| 高配体密度      | `capture_ratio_1p5.yaml` |       32 |

加上基准和弱/强情景，总共约9次正式仿真，已经能完成一套合理的单因素敏感性研究。

# 第七阶段：是否需要多个随机靶区

这一步推荐做，但可以放在固定靶区动力学研究之后。

保持以下参数完全不变：

```text
Influence-region fraction = 0.10
Positive-wall fraction = 0.80
Correlation length = 7.5 µm
Random field modes = 512
```

只改变：

```text
Random seed = 42, 43, 44, 45, 46
```

宏观疾病区域不移动、不扩大，只改变区域内部的微观受体斑块排列。

如果你连微观斑块都不希望改变，那么可以跳过这一步，但最终结论只能写成：

> 结果适用于当前这一张固定的合成靶区掩膜。

不能说明结果对未知靶分布具有普遍性。

# 浓度方面的重要限制

当前验证使用：

$$
C_0=2\times10^7\ \mathrm{MB/mL}
$$

而你此前考虑的物理情景是：

$$
5\times10^5\quad\text{或}\quad10^6\ \mathrm{MB/mL}
$$

当前高浓度适合增加计算样本，但因为代码启用了微泡碰撞，它不等价于单纯“多复制一些互不影响的微泡”。因此建议：

1. 先用 $2\times10^7$ 完成弱、中、强参数筛选；
2. 将结果明确标记为计算富集情景；
3. 对最重要的三个情景，再用 $10^6\ \mathrm{MB/mL}$ 做确认；
4. 高浓度与物理浓度结果分别报告，不能混合。

按照当前约 $0.59\%$ 的接触率，在 $10^6\ \mathrm{MB/mL}$ 下获得10个独立接触微泡可能需要非常长的物理时间。因此，如果后续要大规模研究，最值得新增的是“保存轨迹的离线多靶区重分析/加权试验粒子”功能，以避免每个靶区和参数组合都重新求解整个仿真。

# RBC-induced drift–diffusion reduced-order transport model 正式设置

当前正式运行采用单向耦合的 RBC 诱导漂移—扩散降阶输运闭合：

- 根血管排出血细胞比容固定为 $H_{D,0}=0.35$；
- 红细胞主直径固定为 $d_R=8\,\mu\mathrm m$，且只用于单列流—多列流尺度激活；
- 不进行 $H_{D,0}$ 敏感性对照；
- 局部 $H_D$ 由 Pries–Secomb 2005 二叉分叉相分离关系逐级计算；
- 局部 $H_t$ 由 Fåhræus 关系计算；
- 红细胞项只增加微泡平移速度，不修改 CFD、黏度、迁移率矩阵、力矩或分子参数；
- 已进入平均无细胞层、横向空间无效或尺度开关为零时，RBC 漂移和扩散严格为零；最近壁面不唯一时条件平均漂移与 Fick 修正为零；
- 超出定量血细胞比容或剪切率范围的使用只记录并警告，不截断、不调参；
- 主要靶区暴露量报告为 $\mathcal E_T=\sum_i\int A_{T,i}(t)\,dt$，只有 $\mathcal E_T>0$ 时才解释后续分子成键比较。

正式运行另存 `red_blood_cell_transport.npz`，用于审计每条血管的 $Q$、$D$、$H_D$、$H_t$、$\delta_{CFL}$、尺度激活、分叉相分离参数及红细胞通量守恒误差。粒子轨迹默认保存 RBC 确定性速度分解、$D_\perp$、定量适用性和横向空间有效性；其余兼容诊断仅在 `store_full_diagnostics` 开启时保存。
