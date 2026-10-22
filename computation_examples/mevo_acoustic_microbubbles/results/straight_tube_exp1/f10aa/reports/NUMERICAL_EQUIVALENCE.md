# 完整数值一致性

首轮96个冻结输入后端状态（CPU1/2线程、CuPy、OFF）全部通过原Gate；启动前固定2/4线程又完成64个状态及6条完整轨迹，均通过。32/64/128 固定 ID 的完整窗口 CPU/GPU 最大差 2.4666e-11 nm；不是精度容差被放宽。dt/2 与细缓存的最大逐泡差 1.60274 nm，原门槛 50 nm。

COM、Rg、Q、长短比、XY 伸长、可辨主轴方向、固定 ID 及无ID分布匹配、最近邻、几何边事件序列和时间全部对照；见 full/COMPARISONS 与 analysis/EXTRA_COMPARISONS。主轴近退化时不伪造方向准入。形态数据从全部接受帧计算，不只看末态。

65,259 个接受步独立重算原 resistance 与最终真实法向的力平衡、对偶、真实间隙互补；另按原单位门槛重查全部逐步 residual，all_pass=True。压力逐步范围、OFF=0、源幅相和事件时序沿用原代码，非零压力严格 0.1–40000 Pa。

EXPLORATORY · PHYSICS_UNCALIBRATED · HOMOGENEOUS_WATER_CONTROL_PROXY · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
