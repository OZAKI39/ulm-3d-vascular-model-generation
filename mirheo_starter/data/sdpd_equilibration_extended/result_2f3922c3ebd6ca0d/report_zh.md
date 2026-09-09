# 固定参数 SDPD：真实恢复测试失败与 CPU 原因核查

**RESTART_NOT_VALIDATED。三进程恢复诊断已执行，恢复读取失败，长实验按约定停止。**
本轮计费 11.711100826971 GPU 秒；CPU 修复验证 PASS。没有自动重试。selection=null；人工验收 PENDING。

## 实际执行与绝对时间

| 独立任务 | 结果 | 实际步数 | 实际 t* | GPU 秒 |
|---|---|---:|---:|---:|
| restart_A | COMPLETED (exit 0) | 4000 | 0.003999999989900971 | 6.211455132 |
| restart_B_save | COMPLETED (exit 0) | 2001 | 0.0020009999949479607 | 3.780379508 |
| restart_B_restore | STOPPED_OR_FAILED (exit 255) | None | None | 1.719266187 |

A、B 分别从同一初态冷启动，不能拼成连续轨迹；此短对照仅检验恢复，不检验液体平衡。A 达到约 t*=0.004；B 保存的真实边界为第 2000 步、t*=0.002，随后多推进 1 步触发写盘。
B-restore 在 u.restart 内退出，尚未得到恢复后推进前状态；不能把计划再推进的 2000 步记为已完成。三进程均未超时，runner 未发现遗留本任务进程。
本轮正式长轨迹 NOT_RUN，(0.60,0.80] 没有新样本。最近完成的液体长观察仍是旧 400000 步、t*=0.40、12.662309354 µs，旧结论与门槛保留。

## 直接失败原因及独立阻塞

原生日志：`Unrecognised form Other`，定位到 `core/xdmf/xmf_helpers.cpp:103`。实际 `pv.PV.xmf` 中 `saved_forces` 为 `Other / Force`；HDF5 实际形状是 `[4096, 1]`。`type_map.h` 未给 Force 定义 DataForm，回退到 Other；`channel.cpp` 把它当成 1 分量写出，而读取器拒绝 Other。
`saved_stresses` 是受支持的 Tensor6 / Stress，HDF 形状 `[4096, 6]`，不是这次错误源。将 saved_forces 标签改成 Vector 不能补回不存在的三维数据；本轮没有改写原始 checkpoint。
项目 Python 的问题是：启用了无法原生序列化的持久化诊断 Force 通道，却仅检查文件哈希。已先用失败测试复现，再增加 CPU 通道格式检查，在导入 Mirheo/CUDA 和调用 restart 前拒绝该文件。该修复是明确拒绝无效输入，不代表 Force 序列化已修好。
另一项独立限制仍然成立：当前 SDPD 相互作用没有接入内核随机状态的 checkpoint/restart，应力包装器还有独立 RNG 与调度状态。即使解决文件读取，也不能据此建立完整随机过程续跑。未修改、重编译或升级原生库。

## 保存状态与恢复比较

| 核查 | 状态 | 证据 |
|---|---|---|
| 文件完整性 | PASS | manifest、内部引用和每个文件 SHA 验证 |
| 相同冷启动初态 | PASS | 按 int64 粒子 ID 比较，位置/速度最大差均为 0 |
| HDF 与 B 自身第 2000 步 | PASS | 按 ID 比较，位置/速度最大差均为 0；使用原 2e-6 门槛 |
| 独立 A/B 第 2000 步 | FAIL | 最大速度差 5.245208740234375e-06*，原门槛 2e-6*；已在恢复前出现 |
| 恢复后推进前、第一步、短段 | NOT_RUN / NOT_REACHED | 原生读取失败，缺实测状态 |
| 温度/压力/密度/动量恢复跳变 | NOT_MEASURED | 不用均值接近或合成数据替代 |
| SDPD RNG 连续性 | 未验证 | 原生持久化契约缺失 |

A/B 的保存前速度差已超过原先用于 A 对恢复 B 的严格门槛；它不是恢复造成的跳变，也不能忽略。这里只补充描述，未事后放宽门槛。以后验证应单独区分“恢复状态对其自身保存状态”和“独立进程演化差异”，需要在新试验前冻结设计。

## 授权、停止和余额

批准已登记：追加 1379 秒，绑定原计划 a68cad9d265e78ebf5f4606f0d537c7c5f417dba7254bf0144334832e242b3bd 与指定 7 个任务。原 3600 秒基额和旧 493 秒扩展保留，总授权 5472 秒。
全局累计计费 4032.996053298019 秒；可用于本范围的余额 1439.003946701981 秒，其中原基额余额 71.715047528953 秒。
原批准范围上限为 3×30 + 4×340 = 1450 秒；失败后不支出 1360 秒条件长段。未花余额留在原用途账户，不清空、不回填、也不自动授权重试。本次新增申请为 0 秒。
并发 1，单任务上限 600 秒；没有改成单个超上限长进程。正式窗口仍为 (0.60,0.80]，趋势窗口仍为 (0.40,0.60]，最长 0.80*。

## 可复现的只读分析与页面

```bash
cd /home/lzy/projects/mirheo_starter
.venv/bin/python -B -m py_scripts.run_sdpd_equilibration --config py_scripts/sdpd_equilibration_extended.yaml --analyze-only
.venv/bin/python -B -m test_code.review_sdpd_equilibration --config py_scripts/sdpd_equilibration_extended.yaml --open
```
已有 terminal 失败记录时，--execute 也只返回分析包，不重新启动 GPU。执行计划和各进程 frozen_code 保留原哈希；CPU 分析修订另记代码身份，不把修订后的代码冒充此次 GPU 所用版本。
读取失败后自动生成的旧 preparation_779eb32f8aae47f0 包中，NOT_RUN 总状态和重复预算请求是 Python 汇报错误。原包保留，本报告及 restart_execution.json 对其作明确更正。

下一优先项：先解决持久化通道与 SDPD 随机状态的完整恢复契约，再做新的非零 dt 验证。当前约束禁止原生修改，因此本轮在 CPU 证据和拒绝保护处停止，不能宣称可靠恢复或后段平衡已完成。

## 历史物理结果的补充证据（不属于本次短诊断）



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
本轮正式液体长观察未测量；短恢复诊断不作为物性证据，历史补充趋势不升级为新实验验收。

