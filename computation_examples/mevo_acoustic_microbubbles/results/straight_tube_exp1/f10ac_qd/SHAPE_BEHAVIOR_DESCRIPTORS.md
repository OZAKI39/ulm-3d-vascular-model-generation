# 形态描述量和档案规则

描述量为真正三维的 Δlog(Rg) 与 Δlog(Lmax/Lmin)，范围 [-0.08,0.08) × [-0.16,0.16)，固定 8×8 档案。越界结果保留，不截断、不塞入边界格。质量 = 100/(1+|COM 位移|/1 µm)，仅完整合法数值轨迹可入档。Task A 允许固定粒子集合分裂；Task B 要求终点一个几何连通分量；DT1 的 10 µs（DT05 为 5 µs）采样路径连通另列，不据此声称每个中间步均连通。

保留 Rg、形状张量、主轴长度、方向可辨识性、Kabsch/相似变换残差、非仿射残差、粒子距离、连通分量和接触反力。档案格数不等于独立物理模式数。

`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `training_ready=false`
