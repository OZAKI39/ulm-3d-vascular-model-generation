# DPD 与原生 SDPD 纯液体对照

状态：**PARTIAL**；selection=None；用户人工验收 PENDING。

固定目标：ρ=1056.0 kg/m³，ν=3.27e-06 m²/s，μ=0.00345312 Pa·s，T=298.15 K。
原保存单位 L0=4.9999999999999998e-07 m，M0=1.6499999999999997e-17 kg，t0=3.1655773384195756e-05 s；目标 ν*=414.057515865、μ*=3312.46012692。没有重新定义时间单位。

## 材料判定

| 模型/候选 | 实测 ν (m²/s) | 95% CI | 相对目标误差 | 流动态 T/T目标 | 结论 |
|---|---:|---|---:|---:|---|
|DPD/dpd_cost_baseline|未得到有效结果|None|未得到有效结果|未得到有效结果|COST_ONLY|
|DPD/thermal_baseline|8.994092e-09|[8.830366292781555e-09, 9.164004084546772e-09]|0.9972495|1.026059|NOT_QUALIFIED|
|DPD/thermal_high_gamma|2.427875e-06|[2.323885232196549e-06, 2.5416069504261685e-06]|0.2575307|1.05188|NOT_QUALIFIED|
|SDPD/sdpd_target_linear|3.099559e-06|None|0.05212268|1.061956|NOT_QUALIFIED|

### DPD / dpd_cost_baseline
未通过或缺证：density, temperature, stationarity, viscosity, step_and_force_sensitivity, mechanical_pressure_definition, EOS_statistics, EOS_temperature, pressure_range, density_adjustment, measured_EOS_Mach。
平衡温度 T*：[]；步长/驱动：INCONCLUSIVE_OR_FAILED；机械 EOS：INCONCLUSIVE；压力覆盖：UNVERIFIED。

### DPD / thermal_baseline
未通过或缺证：temperature, viscosity, EOS_temperature, pressure_range, density_adjustment。
平衡温度 T*：[1.0260233152160265]；步长/驱动：PASS_PROPOSED；机械 EOS：MEASURED_LOCAL_RESPONSE；压力覆盖：OUTSIDE_MEASURED_EOS_RANGE。

### DPD / thermal_high_gamma
未通过或缺证：temperature, viscosity, step_and_force_sensitivity, EOS_statistics, EOS_temperature, pressure_range, measured_EOS_Mach。
平衡温度 T*：[1.0517656840781202]；步长/驱动：INCONCLUSIVE_OR_FAILED；机械 EOS：INCONCLUSIVE；压力覆盖：UNVERIFIED。

### SDPD / sdpd_target_linear
未通过或缺证：planned_execution_complete, temperature, stationarity, viscosity, step_and_force_sensitivity, EOS_temperature, pressure_range。
平衡温度 T*：[1.1462786934528366]；步长/驱动：INCONCLUSIVE_OR_FAILED；机械 EOS：MEASURED_LOCAL_RESPONSE；压力覆盖：NOMINAL_COVERAGE_UNCERTAIN。

## 实际执行与费用

旧账本已用 1695.861351 s；本轮新增 1599.325158 s；合计 3295.186509 / 3600 s；剩余 304.813491 s。

| 任务 | 来源 | 墙钟 s | 实际物理时长 s | 有效统计 s | 每步总成本 s | 退出状态 |
|---|---|---:|---:|---:|---:|---|
|probe|HISTORICAL_REFERENCE|3.479425|0.0001266231|未得到有效结果|0.001739712|COMPLETED|
|equilibrium|HISTORICAL_REFERENCE|118.0776|0.007913943|0.006331155|0.0009446206|COMPLETED|
|eos_low|HISTORICAL_REFERENCE|118.1342|0.007913943|0.006331155|0.0009450739|COMPLETED|
|eos_high|HISTORICAL_REFERENCE|115.1067|0.007913943|0.006331155|0.0009208537|COMPLETED|
|flow|HISTORICAL_REFERENCE|196.83|0.01266231|0.01012985|0.0009841499|COMPLETED|
|flow_half_force|HISTORICAL_REFERENCE|192.8788|0.01266231|0.01012985|0.0009643941|COMPLETED|
|flow_half_dt|HISTORICAL_REFERENCE|339.8784|0.0109529|0.008762318|0.0009823075|COMPLETED|
|high_gamma_equilibrium|HISTORICAL_REFERENCE|302.4595|9.610693e-05|7.686022e-05|0.0009962435|COMPLETED|
|high_gamma_flow|HISTORICAL_REFERENCE|309.0166|9.610693e-05|7.686022e-05|0.001017841|COMPLETED|
|sdpd_probe|MEASURED|6.132219|1.266231e-07|未得到有效结果|0.001533055|COMPLETED|
|sdpd_equilibrium|MEASURED|210.9335|5.444793e-06|2.722397e-06|0.001226357|COMPLETED|
|sdpd_flow|MEASURED|326.464|8.559721e-06|4.279861e-06|0.001207337|COMPLETED|
|sdpd_half_dt|MEASURED|376.78|4.792684e-06|2.393176e-06|0.00124432|COMPLETED|
|sdpd_half_force|MEASURED|288.4865|6.957939e-06|3.475804e-06|0.001312496|STOPPED_OR_FAILED|
|sdpd_eos_low|MEASURED|184.4083|4.925638e-06|2.462819e-06|0.001185144|COMPLETED|
|sdpd_eos_high|MEASURED|196.8268|4.925638e-06|2.462819e-06|0.001264953|COMPLETED|
|dpd_cost|MEASURED|4.347089|1.266231e-07|未得到有效结果|0.001086772|COMPLETED|
|dpd_cost_density|MEASURED|4.946769|1.266231e-07|未得到有效结果|0.001236692|COMPLETED|

同观测短成本对照：COMPARABLE_EXECUTION_SETTINGS。密度路径增量：{'extra_s_per_step': 0.00017854870579245275, 'relative_to_plain_DPD': 0.25593723116514466, 'scope': 'Matched passive Density path includes intermediate halos, channel preservation and scheduling; not an isolated CUDA-kernel timer. Single short trial, no statistical speed CI.'}。
完整计时、初始化、平衡段、同步计算、采样传输输出、主机可用内存与设备显存位于 comparison_summary.json。设备显存含其他应用；短测试每步成本没有重复实验的置信区间。
目标液体效率排名允许：False。未达到合格材料时，达到合格解的耗时为 null；不能按低黏度历史基线的速度宣布血流加速。

## 测量定义与局限
SDPD 输入的是动力黏度，实测运动黏度来自独立速度剖面拟合。模拟点、目标解析线和拟合线分别保存。
全盒 N/V、N*m/V、核数密度及进入 EOS 的 m*d_i 分开保存；Wendland C2 包含自贡献。少量原始快照由 CPU 独立重算密度，只作观测核查，CPU 耗时单列。
机械压力采用原生配对应力的迹，加一次动能动量通量并除以体积；不再添加 EOS 压力。EOS 局部预测与机械压力的差异保留，不能当作已校准的绝对热力学压力。
旧三个表压差原值保留。每个候选以本身的实测平衡机械压力作为测试零表压参考，不套用旧 LBM 数值偏置或大气压；移动零点不扩大已测压力区间。
压力均值名义覆盖与包含不确定性的覆盖分开。保守核查要求端点压力的向内 CI 边界覆盖“基准压力加表压差”的完整基准 CI；各点是单独 95% CI，不宣称联合 95% 覆盖。
原生 Stats/virial 时间差由该版本两个回调和每次实际 CSV 对应验证。核密度快照取 beforeIntegration，和同相位位置一起保存；其时间与步末速度相差一个实际 dt。
局部材料统计门槛与旧血管短/长窗口分别报告。本轮不是空间收敛测试；相同间距不保证相同空间误差。没有将独立初始化拼成连续轨迹。

## 旧数据重新审查
旧温度偏高同时发生于平衡态和流动态，因此不能只归咎于未扣流速。各流动记录另给基于实测剖面的有限箱宽剪切残余估计；没有通过减去任意速度场强行通过温度检查。
旧高耗散候选的时间步/驱动力与 EOS 证据独立判定，不借用低耗散候选。旧 CSV、配置、报告、UTC 字段和账本全部保留；分析修正只写入本结果。

## 后续能力边界
未运行真实血管、SDF、RBC、微泡、黏附、管壁或开放压力出口。后续需要验证近壁核截断、冻结壁粒子、膜内外液体/膜节点耦合，以及双向压力储液区和有限截面压力测量。普通 DPD/SDPD 的体积黏度未独立匹配。

## 下一步
依据本轮实际失败项先解决初始化/温度和统计充分性，再验证黏度及机械 EOS；保留当前未合格结论，不自动选失败候选为生产默认。
