# 第二阶段真实执行报告

总体：**PARTIAL**。合格选择：`None`。用户人工验收：**PENDING**。

## 已验收来源与固定物理条件

来源：`/home/lzy/projects/ulm_3D_vascular/outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_base_anchor003274_20260901/qc/reference_scaled_base_final.json`，接受检查点 598755。实际生成 Lua、独立 runtime contract、当前物理常数、冻结测量截面哈希及其对应壁面已逐项核对。

ρ=1056.0 kg/m³；ν=3.27e-06 m²/s；μ=ρν=0.00345312 Pa·s。运动体积黏度 2.18e-06 m²/s 保留，未独立匹配。

入口 Q=2.73691323909057e-15 m³/s，质量流量=2.890180380479642e-12 kg/s。出口表压：{'outlet_01': 14.544978101274268, 'outlet_02': 132.20454922317552, 'outlet_03': -13.700626673311461} Pa。边界数值是用户确认的测试值，298.15 K 是用户同意的假设，均非生理实测。

旧 LBM 数值压力偏置 3387510.72 Pa 不作为热力学绝对压力。历史零表压为 `GLOBAL_STRUCTURAL_LEAVES_ZERO_GAUGE`。当前配置指向另一个 Musubi 二进制，实际比较以接受报告的生成快照为准，没有执行旧求解器。

## 单位约束与候选

L0=5e-07 m，M0=1.65e-17 kg，t0=3.16557733842e-05 s，E0=4.1164049935e-21 J，P0=0.032931239948 Pa。密度与 298.15 K 将尺度锁定后，要求 ν*=414.057516。黏度匹配由实测单独判定。

n*=8、m*=1、kBT*=1、rc*=1；a=25。gamma=10 是基线，gamma=3700 是量级探针；两者均非生产默认。驱动力与 dt 比较始终使用同一单位。空间间距 0.25 µm，rc=0.5 µm；真实 SDF 间距未定。

## 实际 GPU 执行

|任务|步数|dt*|墙钟 s|退出码|GPU 显存采样峰值 MiB|状态|
|---|---:|---:|---:|---:|---:|---|
|probe|2000|0.002|3.479|0|3033|COMPLETED|
|equilibrium|125000|0.002|118.078|0|3047|COMPLETED|
|eos_low|125000|0.002|118.134|0|3013|COMPLETED|
|eos_high|125000|0.002|115.107|0|3092|COMPLETED|
|flow|200000|0.002|196.830|0|3063|COMPLETED|
|flow_half_force|200000|0.002|192.879|0|3054|COMPLETED|
|flow_half_dt|346000|0.001|339.878|0|3021|COMPLETED|
|high_gamma_equilibrium|303600|1e-05|302.460|0|3024|COMPLETED|
|high_gamma_flow|303600|1e-05|309.017|0|3031|COMPLETED|

累计 1695.861 / 3600 s；余 1904.139 s，未自动追加批次。包含初始化、失败重试、输出和退出。显存为设备整体采样峰值（包括其他应用），非真实任务峰值。

## 温度、平衡与统计

|任务|热温度 kBT* 均值 ±95%半宽|有效块数|温度状态|排除区间 t*|统计持续 ms|旧短/长窗口|
|---|---|---:|---|---|---:|---|
|probe|无充分统计|0|WINDOW_INSUFFICIENT|—|—|WINDOW_INSUFFICIENT|
|equilibrium|1.026023 ± 0.0020019938695537243|83|INCONCLUSIVE_OR_FAILED|[0, 50.0]|6.33115|AVAILABLE_DURATION_ONLY / AVAILABLE_DURATION_ONLY|
|eos_low|1.025308 ± 0.002070452413530003|83|INCONCLUSIVE_OR_FAILED|[0, 50.0]|6.33115|AVAILABLE_DURATION_ONLY / AVAILABLE_DURATION_ONLY|
|eos_high|1.027795 ± 0.0018189011095997398|83|INCONCLUSIVE_OR_FAILED|[0, 50.0]|6.33115|AVAILABLE_DURATION_ONLY / AVAILABLE_DURATION_ONLY|
|flow|1.026059 ± 0.0015808423555595014|133|INCONCLUSIVE_OR_FAILED|[0, 80.0]|10.1298|AVAILABLE_DURATION_ONLY / AVAILABLE_DURATION_ONLY|
|flow_half_force|1.026276 ± 0.0014812277196188893|133|INCONCLUSIVE_OR_FAILED|[0, 80.0]|10.1298|AVAILABLE_DURATION_ONLY / AVAILABLE_DURATION_ONLY|
|flow_half_dt|1.013693 ± 0.0010894514386045064|125|PASS_PROPOSED|[0, 69.2]|8.76232|AVAILABLE_DURATION_ONLY / AVAILABLE_DURATION_ONLY|
|high_gamma_equilibrium|1.051766 ± 0.0012469080784977703|202|INCONCLUSIVE_OR_FAILED|[0, 0.6080000000000001]|0.0768602|WINDOW_INSUFFICIENT / WINDOW_INSUFFICIENT|
|high_gamma_flow|1.05188 ± 0.001241461934543145|202|INCONCLUSIVE_OR_FAILED|[0, 0.6080000000000001]|0.0768602|WINDOW_INSUFFICIENT / WINDOW_INSUFFICIENT|

平衡温度扣除 COM；流动扣除瞬时分箱均流。分块至少覆盖 5 个积分相关时间、2 个黏性松弛时间，至少 8 块。旧短/长窗口可用仅表示统计跨度足够；本轮未执行真实血管流量/压力验收。

## 黏度与敏感性

|任务|ν*|ν / m² s⁻¹|95% CI / m² s⁻¹|剖面相对 RMS|采样判定|目标匹配|
|---|---:|---:|---|---:|---|---|
|flow|1.1388598|8.9940922e-09|[8.830366292781555e-09, 9.164004084546772e-09]|0.028961|SUFFICIENT|NOT_QUALIFIED|
|flow_half_force|1.1487649|9.0723176e-09|[8.625597038867708e-09, 9.567836817873435e-09]|0.050673|SUFFICIENT|NOT_QUALIFIED|
|flow_half_dt|1.1403049|9.0055053e-09|[8.806733875349614e-09, 9.213456599369185e-09]|0.029018|SUFFICIENT|NOT_QUALIFIED|
|high_gamma_flow|307.42501|2.4278748e-06|[2.323885232196549e-06, 2.5416069504261685e-06]|0.053266|SUFFICIENT|NOT_QUALIFIED|

force / 2：均值相对变化 0.00869742264222051，保守 CI 比值差界 0.08351528131904673，PASS_PROPOSED。

dt / 2：均值相对变化 0.0012689586769398797，保守 CI 比值差界 0.04338328602526831，PASS_PROPOSED。

拟合使用每粒子力 F，而不是加速度；保留 m，体力密度为 nF。解析曲线只用于拟合实测速度点，包含箱宽二次曲线平均修正和自由常数偏置。没有壁面，因此不能证明真实管壁无滑移。

## 局部 EOS 与压力可达范围

状态：MEASURED_LOCAL_RESPONSE。实测范围：[5.070158119893552, 6.197883794960322] Pa。平衡参考：5.6121090780820255 Pa。将旧表压加到该参考后目标为 {'outlet_01': 20.157087179356292, 'outlet_02': 137.81665830125755, 'outlet_03': -8.088517595229435} Pa；OUTSIDE_MEASURED_EOS_RANGE。

dp*/dρ*=43.023990015701486 ± 1.1928793457683569；这是三个实际密度的局部拟合，温度状态 NOT_AT_ACCEPTABLE_TARGET_TEMPERATURE，不能外推成全局 EOS 或已验证的开放出口。

原生输出为 virial 总和，分析补入动能项和体积归一化；配对 0.5 因子已在源码中。Stats 与 virial 同回调采样，但 CSV 标签相差一个 dt，已按源码和实测偏移处理。

## 候选结论

- `thermal_baseline`：NOT_QUALIFIED；EOS_NOT_AT_ACCEPTABLE_TARGET_TEMPERATURE；EQUILIBRIUM_TEMPERATURE_NOT_QUALIFIED；LOCKED_UNIT_VISCOSITY_TARGET_NOT_MATCHED_WITH_NARROW_CI；MEASURED_TEMPERATURE_NOT_WITHIN_PROPOSED_GATE；TARGET_PRESSURES_OUTSIDE_MEASURED_EOS。
- `thermal_high_gamma`：NOT_QUALIFIED；CANDIDATE_SPECIFIC_DT_AND_FORCE_SENSITIVITY_NOT_MEASURED；CANDIDATE_SPECIFIC_EOS_NOT_MEASURED；EQUILIBRIUM_TEMPERATURE_NOT_QUALIFIED；LOCKED_UNIT_VISCOSITY_TARGET_NOT_MATCHED_WITH_NARROW_CI；MEASURED_TEMPERATURE_NOT_WITHIN_PROPOSED_GATE。

没有合格候选时不选择“最不差”的一组，不把仅调 t0 匹配黏度但温度不符的映射当作解决方案。

## 可实施的边界计划

每个端口沿用第一阶段端盖和既有延伸段，保留旧有限测量轮廓及检查截面。控制层拟设在实际端盖内侧 [−rc,0]，横向由有限端口和连通腔体限定；这是设计而非已生成几何。入口采用可控局部供液，以保留测量截面的实测平均 Q 反馈。Q/A 不自动定义速度剖面。

原生 VelocityInlet、DensityControl 和 DensityOutlet 可作为供液、局部反馈力和超额密度移除组件。仍需任意有限截面的有符号通量/牵引力测量、可调供液及双向压力储液插件；现有 DensityOutlet 不补液，PlaneOutlet 不是压力边界且不能限制多分支局部区域。不得删除反向粒子或改压力通过回流门槛。全部端口目前为 REQUIRES_EXTENSION。

固定壁面可复用原生墙面约束与冻结壁粒子的接口设计，但要先在合成通道核验滑移、漏液、壁密度和离散误差。本轮未开发或编译 C++/CUDA 插件。

## 进入③的条件

几何身份、法向、单位和旧测量位置足以讨论 SDF 坐标域与数据组织。正式分辨率与初始化仍受材料黏度、目标温度、压力范围及壁面/局部储液边界限制；本轮空间候选不能直接写成已验证默认值。高 gamma 下一轮优先检查更小 dt、足够统计时长以及同时满足压力刚度的候选，再考虑实际血管。

## 来源与输出保护

Mirheo 提交 `8fa67b9aaa7f04c9de2d74a335c9c8c4665068cf`，精度 `single`，库 SHA256 `1bab2b922b43a33dd0ffc84211234f3cf75f8edc852e30adfcd17f68c8c78c87`。源码、共享库、补丁及第一阶段数据哈希在 provenance 和保护审计中保留；未重编译或替换环境。

本次仅完成旧纯流体工况整理、单位映射、预算内 DPD 标定及入口出口实现方案评估。未运行整段真实血管，未生成真实血管 SDF 或完整粒子模型，未加入 RBC、微泡或黏附。本阶段用户人工验收状态为 PENDING。
