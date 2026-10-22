# 原数值要求与独立验证

独立重建阻力矩阵、终点接触法向、力平衡、互补性及非穿透：检查 697500 个接受步，完整轨迹独立审计全部通过=True。原物理及接触门槛未放宽。测试结果：........................                                                 [100%]
24 passed in 1.72s

三个档案候选的 24/24 对照通过；另一个档案外强收缩候选的三个新初态 dt/2 对照出现 104.9、318.9、112.6 nm 差异，超过原 50 nm 门槛，整体跨初态验收 FAIL。完整共享物理时刻的逐泡最大差异、COM、Rg、Q、长宽比、方向及事件时序见 VALIDATION_COMBINED.json 和 ../f10ad_multifrequency/LONG_CONVERGENCE.json。几何连边时序不等于物理融合。数值收敛不等于模型完成实验验证。

`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `training_ready=false`
