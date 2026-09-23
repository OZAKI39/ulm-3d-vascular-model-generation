# Particle-1 小步骤 7 检查说明

三个时间步在首次运行之前由局部最短边长和起点速度计算，固定逐级减半；没有根据结果调 dt。三个重放均从 OUTLET_02 离开，相邻两次的出口时间差为 7.906e-05、7.243e-05 s。这里报告数值趋势，不选择 production dt，也不称作 FEM 时间步研究；轨迹是否有视觉上突然跳跃仍待用户审核。

![检查图](figures/07_validation_timestep_comparison.png)
