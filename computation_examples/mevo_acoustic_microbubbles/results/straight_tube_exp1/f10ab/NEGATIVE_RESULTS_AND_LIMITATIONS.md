# 失败、未运行与不可外推范围

| 轨迹 | 请求/ms | 接受/ms | 原因 |
| --- | --- | --- | --- |
| pilot__S1_r0__A3_X0__500us__DT1 | 0.5 | 0.048 | CONTACT_REACTION_BALANCE_PASS |
| pilot__S2_r0__A3_X0__500us__DT1 | 0.5 | 0.043 | CONTACT_REACTION_BALANCE_PASS |
| pilot__S3_r0__A3_X0__500us__DT1 | 0.5 | 0.033 | CONTACT_REACTION_BALANCE_PASS |
| pilot__S4_r0__A3_X0__500us__DT1 | 0.5 | 0.08 | CONTACT_REACTION_BALANCE_PASS |
| pilot__S5_r0__A3_X0__500us__DT1 | 0.5 | 0.024 | CONTACT_REACTION_BALANCE_PASS |

- 这些5条A3同相高幅轨迹在早期被反力平衡门槛拒绝，保留接受前缀、失败提案和日志；不能补帧称完整。
- A7在不同初态的长宽比变化方向不一致。
- A5方向可重复，但可能裂为多个邻近分量。
- 局部形状响应第二方向低于保守分辨门槛。
- 3%目标涉及连通性失败，不能只展示接近目标的长宽比。
- N128新增形状试验、10ms/1s、可逆控制、Y完整轨迹、PDMS校准、AI/MPPI、硬件均NOT_RUN。N32/64已足以发现当前机制和完整性限制，未为规模展示增加负担。
- 近接触物理、APECSS参数、均匀水空间梯度和电声转换未校准；数值收敛不等于真实装置精度。

适用范围：`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED` · `training_ready=false`。Y 与四源静态结果另加 `HOMOGENEOUS_WATER_PROXY_SPECIFIC`。
