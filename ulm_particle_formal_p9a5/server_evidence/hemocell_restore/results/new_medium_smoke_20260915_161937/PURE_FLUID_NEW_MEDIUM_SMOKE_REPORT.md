# PBS/BSA 开发介质纯流体 smoke 核查

结果：**PASS**。这是新介质的短期数值稳定性与单位转换检查，不是长期收敛或真实实验介质验证。

介质为 1× PBS + 1% BSA，25°C；ρ=1000 kg/m³、μ=0.001 Pa·s、ν=10⁻⁶ m²/s 均为开发阶段假设，未做实验测量。最终实验合同仍为 PENDING。旧整血等效纯流体基线 PASS 保留。

当前 prepare_numerics.py 唯一生成 dt=6.6599470814617197e-09 s，tau=1.0，有效 dx=1.9989918081065344e-07 m。pressure_unit=900908.05926991301 Pa。三出口物理表压不变；生成密度为 1.0000484343922278, 1.0004402376508774, 0.99995437727568648。

压力换算沿用 dt=((tau−0.5)/3) dx²/ν、pressure_unit=ρ(dx/dt)²、rho_LU=1+3 Pg/pressure_unit。求解器只读合同与生成参数文件，并逐值核对；独立 finalizer 用相同物理定义核验单位，核验结果不作为新输入。

生成器完整 JSON 的 SHA256 在 NUMERICS_GENERATION_RECEIPT.json 中；合同内部 output_sha256 哈希指定的 canonical generated_output_digest_payload，避免文件自哈希循环。运行前修正过回执中 input 路径的变量覆盖，旧回执/脚本/合同保存在 provenance/PRE_AUDIT_RECEIPT_FIX_*，数值内容不变。

GPU 核心（batch dispatch、halo、managed memory、Guo、监控归约、kernel layout）与已验证 Stage4 相同。独立 host driver 只适配合同读取、500/5000 步任务边界及本任务指定输出节奏；初始化/Guo/边界数学字节核查见 DRIVER_STATIC_PROOF.json。

## 实际执行与安全检查

500 步与 5000 步各一次，均从 ρ_LU=1、u=0 初始化，ramp=10 步；MPI1，单物理核心绑定。每步保留 nonfinite 与原 GPU 安全归约，50 步记录状态，100 步记录全部 24 个测量面流量，500 步记录控制体质量。完整快照只有各运行 step0 与 5000 步运行的 step5000。没有运行 RBC，也未启动后续任务。

|运行|核查|solver steps/s|端到端 s|初始化 s|VRAM 峰值 MiB|timed GPU中位利用率|
|---|---|---:|---:|---:|---:|---:|
|500|PASS|164.74195462172528|335.8464555340179|330.290893664991|2106.0|52.0%|
|5000|PASS|163.82313256365316|363.9784238199936|330.6619134609937|2106.0|53.0%|

性能只记录，不作为通过门槛。当前更稀疏的流量输出属于用户要求的 smoke 监测节奏，不能据此宣称 GPU 核心优化加速。MPI1 初始化保留原实现，端到端计时包括初始化、必需输出及退出。

全程物理流体 rho 范围 0.9999504812820089–1.000477243174689；流体与 ghost 的最大 Mach=0.000495550987291592；nonfinite=0。既有 density 安全范围 [0.99,1.01]、Mach 硬门 0.05；本 smoke 在运行前设定通过余量 Mach≤0.01。

5000 步末总质量相对初值漂移=0.00014070070071214286，50 步最大质量跳变比例=4.632535680934869e-06，CV 间隔/累计最大闭合误差除以初始 CV 质量分别为 3.8114908974310578e-06 / 3.8114908974310578e-06。

质量门槛、快速单调失控规则、CV 误差与 flux 粗筛均预先冻结于 contracts/SMOKE_SAFETY_CONTRACT.json。CV 闭合采用 storage change + 100 步观测的净向外质量通量梯形积分；这是 transient 简单闭合检查，不是长期 R_flow≤1% 门。

末态主测量面（g2，冻结的4dx面）有符号流量 m³/s：{'Qin': 2.686280389327215e-15, 'Qout01': -1.2964489132866903e-16, 'Qout02': -1.856687265880749e-15, 'Qout03': 5.6906761302776e-16}。入口按向内为正，出口按向外为正；全24面最大 |Q|/Qtarget=4.086220633293507。

短暂回流记录：{"Qout01": {"count": 50, "iterations": [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000]}, "Qout02": {"count": 50, "iterations": [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000]}, "Qout03": {"count": 0, "iterations": []}}。压力驱动 transient 可以出现负出口流量；有限性、符号约定、极端流量门及 snapshot 复算均独立检查，未把负出口流量自动判失败。

原 multiplier=1.1197286861799598 仅沿用数值入口实现，新介质校准 NOT_PERFORMED；不宣称长期 Qin 精度，也未重新校准。未对不同粘度的旧/新介质要求流场相等。

内存 steady-thirds 核查：[{"samples": 17, "vram_median_MiB": 2106.0, "RSS_median_MiB": 7498.11328125}, {"samples": 16, "vram_median_MiB": 2106.0, "RSS_median_MiB": 7498.12890625}, {"samples": 17, "vram_median_MiB": 2106.0, "RSS_median_MiB": 7498.1171875}]；增长判定 {"vram_median_MiB": {"delta": 0.0, "limit": 256.0, "growing": false}, "RSS_median_MiB": {"delta": 0.00390625, "limit": 512.0, "growing": false}}。只说明本短测试未达到持续增长门槛，不证明无限长运行无泄漏。

## 独立终态核查与适用范围

finalizer 独立读取 numerics、每步安全记录、运行诊断和二进制快照，复算最终密度/速度/Mach/总质量/CV质量，以及所有24组多平面插值通量，按预先冻结容差验证。具体检查逐项记录在 verification/SMOKE_*_AUDIT.json。顶层三个历史 CSV 的 run_horizon 区分两个从零运行，避免把两段拼成连续时间。

旧 dt=2.0366810646671923e-9 s；新 dt / 旧 dt=3.2700000000000005。固定 dx 与 tau 时，较小物理运动粘度对应每 LBM 步更长物理时间。5000步物理时间=3.3299735407308598e-05 s。加入 RBC 后有额外计算量与时间尺度约束，不能推导 RBC 必然快3.27倍。

通过仅允许 SUSPENDING_MEDIUM_NUMERICS=PASS、SUSPENDING_MEDIUM_CONTRACT=PROVISIONAL_PASS；EXPERIMENTALLY_MEASURED=NO、FINAL_EXPERIMENTAL_CONTRACT=PENDING、RBC_GPU_CORRECTNESS=PENDING。后续建议 RBC_ONLY_VALIDATION_STAGE_1，但本任务没有启动它。

旧源文件、几何与冻结输入 hash 核查见 verification/SOURCE_AND_INPUT_FINAL_AUDIT.json。远端报告为不可变运行证据；下载后的 WSL SHA256 验证及独立重算记录另见 LOCAL_TRANSFER_VERIFICATION.json 与 LOCAL_NUMERICAL_AUDIT.json。
