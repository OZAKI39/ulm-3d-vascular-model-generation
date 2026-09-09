# 固定参数 SDPD：历史补充分析与恢复准备

状态：RESTART_NOT_VALIDATED。本轮新 GPU 测量 NOT_RUN；CPU 测试 PASS。
长轨迹另有原生能力阻塞：RESTART_NOT_VALIDATED。预算获批也不能把当前粒子 checkpoint 当作完整随机过程恢复。selection=null；人工验收 PENDING。

## 来源与原参数

配置固定引用已交付记录；逐一验证旧包、原始文件、报告和原生库哈希。旧实验 400000 步，t*=0.4，12.662309353678 µs；原结论 STATIONARITY_OR_SAMPLING_INCONCLUSIVE 保留。
旧目录没有 state/XMF/HDF5 checkpoint。诊断快照不是续跑状态；规划独立冷启动，不拼接历史轨迹。
精确冻结参数：{"m_star": 1.0, "n_star": 8.0, "kBT_star": 1.0, "rc_star": 1.0, "id": "sdpd_target_linear", "method": "SDPD", "dt_star": 1e-06, "density_interaction": {"kernel": "WendlandC2", "rc_star": 1.0}, "EOS": "Linear", "sound_speed_star": 120.0, "rho_0_star": 0.0, "viscosity_source": "exact_required_mu_star_from_frozen_mapping", "rationale": "cs=120 predicts max positive pressure demand within 3.5% density change; rho_0=0 keeps pair background pressure positive as in native tests, whose mechanical effect is measured. dt=1e-6 is below acoustic/viscous estimates and is still subject to probe and half-step checks.", "viscosity_mu_star": 3312.4601269222444, "nu_prior_star": 414.0575158652805}
密度 1056 kg/m³；动力黏度 0.00345312 Pa·s；运动黏度 3.27e-6 m²/s；298.15 K；dt*=1e-06，N=4096，盒 [8,8,8]。
单位：{"L0": 5e-07, "M0": 1.6499999999999997e-17, "t0": 3.1655773384195756e-05, "E0": 4.116404993500001e-21, "si_per_star": {"length": 5e-07, "area": 2.5e-13, "volume": 1.2499999999999998e-19, "time": 3.1655773384195756e-05, "velocity": 0.015794907106885803, "acceleration": 498.9581810303033, "particle_mass": 1.6499999999999997e-17, "mass_density": 132.0, "number_density": 8.000000000000001e+18, "force": 8.232809987000002e-15, "pressure": 0.03293123994800001, "dynamic_viscosity": 1.0424638690544628e-06, "kinematic_viscosity": 7.897453553442902e-09, "volume_flow": 3.94872677672145e-15, "mass_flow": 5.212319345272314e-13, "energy": 4.116404993500001e-21}}
观测只复制/归约通道，没有添加力、额外恒温器或速度缩放。

## 旧数据的 CPU 补充证据

正式窗口仍为 (0.30,0.40]，500 个样本。温度描述均值 1.02282627182*；四分窗均值 [1.0307765105794657, 1.0227723754271538, 1.0214265230574044, 1.016329678198956]；两半变化 0.789634%，超过原 0.5% 门槛。这不是可信稳态均值。
All four fixed quarter means decline; negative HAC slope persists across lags 20/40/80. This supports continuing cooling under the diagnostic assumptions, although the old 3-RMS coherent-drift rule did not pass. No steady temperature bias is identified.
温度：四段均值均下降，相关噪声修正下仍支持末段冷却；不能据此确定稳定偏差。HAC 的线性趋势、残差弱相关与渐近区间假设写入 JSON，未替换旧验收规则。
压力：每 200 步一个瞬时值，并非 200 步平均。原生应力包含保守、耗散及随机力。窗口标准差 172.970483 Pa；门槛 0.263449919584 Pa。
压力四分窗均值不单调，多个 HAC 尺度的斜率区间均跨零：尚未证明存在持续漂移。旧同相位应力归约与 CPU 核密度核查已通过；没有删随机应力或用输入 EOS 代替机械压力。
即使独立同分布，500 个样本的平均压力 95% 半宽也约 15.162 Pa；这只是乐观算例，不是正式 CI。条件成本：[{"criterion": "mean_95pct_halfwidth_at_thermal_pressure", "optimistic_iid_samples": 1655999, "duration_star": 331.19980000000004, "physical_duration_s": 0.010484385813690959, "proportional_wall_s_at_old_cadence": 437486.472143911}, {"criterion": "sum_of_two_half_mean_halfwidths_at_thermal_pressure", "optimistic_iid_samples": 13247990, "duration_star": 2649.598, "physical_duration_s": 0.0838750738472183, "proportional_wall_s_at_old_cadence": 3499891.248785664}]
六个保存帧额外重建保守/耗散 virial，随机残差另列；与旧 CPU 重建差异最大 0.00727493*。三个末段帧的保守压力升高约 26.294499 Pa，支持结构仍在松弛、背景被噪声遮盖的解释；三帧不能证明全窗口总压力漂移或给方差归因百分比。
结构：核密度均值/方差持续下降，末段最近邻 q10 依次约 0.34867、0.36336、0.37626；范围 0.02759 rc 超过原 0.02 rc。三帧支持观察到分布变化，但不能给出充分独立的分布稳定置信度；方差缩小不等于结晶。
独立信息：物理最小块长 40 个样本，500 点最多容纳 12 个完整块；这不是 12 个已验证独立块。过渡态 ACF 不被当作稳态 ACF，稳态 CI 留空。更晚 1000 点最多容纳 25 块，不能保证通过。
PROPOSED（未应用）：评估逐步累计完整机械压力、仍每 200 步输出时间平均的开销和相关性；旧瞬时定义、压力门槛、旧结论均保留。

## 各指标分开判断（历史固定窗口）

| 指标 | 状态 | 说明 |
|---|---|---|
| TEMPERATURE_STATIONARITY | INCONCLUSIVE | Frozen point screen exceeded; this alone does not establish persistent drift. |
| TEMPERATURE_TARGET_MATCH | INCONCLUSIVE | Steady temperature and its uncertainty have not been established. |
| PRESSURE_STATIONARITY | INCONCLUSIVE | Frozen point screen exceeded; this alone does not establish persistent drift. |
| STRUCTURE_STATIONARITY | INCONCLUSIVE | Fixed distributions and snapshot quantiles, not frozen particle positions. Three late frames do not provide independent ensemble replication. |
| SAMPLING_SUFFICIENCY | INCONCLUSIVE | Block capacity is an upper bound, not established independence. Transient ACF is not a steady-state correlation time. |
| RESTART_VALIDITY | NOT_TESTED | 未测试；见记录 |
| OVERALL_EQUILIBRIUM | INCONCLUSIVE | Per-observable progress retained; no claim about viscosity, EOS validation, walls or biological components. |

### 保存帧机械压力分解（Pa，仍属于历史补充）

| 保存步数 | 保守 virial | 耗散 virial | 随机项及重建残差 | 原生总 virial |
|---|---|---|---|---|
| 200 | 2758.597577 | -562.352898 | -31.948673 | 2164.296006 |
| 100000 | 3294.692193 | -59.308372 | 48.661137 | 3284.044959 |
| 200000 | 3424.701836 | -14.573968 | 253.577434 | 3663.705303 |
| 300000 | 3475.957678 | -5.857365 | -53.096839 | 3417.003475 |
| 350000 | 3491.031580 | -4.661197 | -45.745782 | 3440.624602 |
| 400000 | 3502.252178 | -3.864358 | -108.430253 | 3389.957567 |

此表为 virial 分量；总机械压力还要加一次 COM 扣除后的动能项。没有删除任何随机贡献。
新实验全部未测量；历史补充趋势不升级为新实验验收。

## 保存/恢复契约与实际限制

原生 SDPD 类及其基类没有 checkpoint/restart override，最终继承 MirObject 空实现；内核虽定义 RNG 序列化，但没有相互作用层调用。启用应力时还有复制内核的 RNG 与 StressManager 时间状态。DPD 的实现不能证明 SDPD 可恢复。当前构建的符号检查也只看到 DPD 与 MirObject 的相关实现。
官方 interactions.py 用 DPD、dt=0；particle_vector.py 是 PV-only、dt=0；均不能代替目标非零 dt 测试。
已实现退出后版本归档、物化软/硬链接、XMF 内部引用与 SHA 校验、完整标记、按 ID 比较、协调器步数与时间验证。文件完整性与随机状态完整性严格分开；缺时间、RNG、文件或哈希错误均拒绝正式恢复。
每次 u.run 重建调度器并将 nExecutions 置零，所以启用 checkpoint 后每个 200 步块起点都会保存，不能把 checkpoint_every=200000 误当成全局周期。原生临时目录用 PingPong 限制体积，退出后物化成不可覆盖版本。每 20 万步段预计 1001 次保存，I/O 开销未实测，18 秒是预算余量而非测量承诺。边界另推进 1 步触发保存；归档以 state.txt 为准，额外步数计费且不混入正式样本。恢复不重抽速度、不扣 COM、不重置时间。
已准备真实对照：A 连续 4000 步；B 保存 2000 步状态后退出（触发保存多推进的 1 步单列），另一个进程恢复并推进 2000 步。比较恢复前、首步及短段的粒子、时间、温度、压力、密度和动量。未运行不称 PASS；源代码 RNG 契约不通过时，即使均值接近也停止。
CPU 测试只验证接口和保护逻辑，不替代真实 GPU 恢复。未修改或重编译原生库。

## 冻结后续设计和完整预算请求

正式窗口 (0.60,0.80]；趋势窗口 (0.40,0.60]；每 0.0002* 采样；最长 0.80*，不自动续到 1.0，也不移动窗口。冷启动 800000 个有效步，4 个 200000 步段，另外 4 步用于边界保存触发，均计费。
旧 400000 步实收 528.365623582999 s；冷启动等比例 1056.731247 s；每段加 15% 和 18 s 未实测开销余量，估计 321.811753 s，每段保留 340 s。
恢复对照 3×30=90 s；条件长实验 4×340=1360 s；一次请求总上限 1450 s。全局原授权 4093 s，已收 4032.996053298019 s；可用于此范围的旧余额 1439.003946701981 s。旧 493 s 扩展已耗尽。
精确缺口 10.996053298019 s；申请整数追加 **11 GPU 秒**。并发 1，单任务不超过 600 s，实收可能更少，失败仍扣实际用量。
该请求是有停止条件的完整上限，不是保证耗尽全部额度：当前原生 RNG 缺口意味着正式 1360 s 不可支出。默认只准备恢复诊断，失败或无法证明 RNG 连续就停止。预算不能单独解除原生能力阻塞；也不自动使用单个超 600 s 的长进程替代。
最大物理范围：25.324618707357 µs；所有新参数与生产 selection 仍为空。
授权通过独立 --register-authorization 命令，绑定真实用户消息、此请求文件哈希与整组 task_ids。预检、查看、执行拒绝分支均不写授权、旧账本或启动 CUDA。

## 可执行命令

```bash
cd /home/lzy/projects/mirheo_starter
.venv/bin/python -B -m py_scripts.run_sdpd_equilibration --config py_scripts/sdpd_equilibration_extended.yaml --preflight-only
.venv/bin/python -B -m py_scripts.run_sdpd_equilibration --config py_scripts/sdpd_equilibration_extended.yaml --analyze-only
.venv/bin/python -B -m test_code.review_sdpd_equilibration --config py_scripts/sdpd_equilibration_extended.yaml --open
# 仅在真实批准并登记本请求后；当前代码会在恢复诊断后按原生缺口停止
.venv/bin/python -B -m py_scripts.run_sdpd_equilibration --config py_scripts/sdpd_equilibration_extended.yaml --execute
```

唯一下一建议：先完成有预算保护的真实非零 dt 恢复对照，确认此原生构建的恢复边界；恢复未验证前不启动新的长观察。
