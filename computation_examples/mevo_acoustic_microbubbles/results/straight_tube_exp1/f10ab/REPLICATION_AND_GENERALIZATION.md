# 发现、独立确认与外推边界

| 形状 | 动作 | 三初态长宽比变化 / % | 方向一致 | 全部终点连通 | 全种子dt/cache复核 |
| --- | --- | --- | --- | --- | --- |
| S2 | A5_X180 | 6.275, 5.621, 8.291 | True | False | True |
| S2 | A7_X210 | -3.71, -3.757, 6.572 | False | True | False |
| S4 | A5_X180 | 5.759, 4.961, 4.811 | True | False | False |
| S4 | A7_X210 | 5.403, -1.642, 5.48 | False | True | False |
| S6 | A5_X180 | 0.2097, -0.09372, -0.09834 | False | False | False |
| S6 | A7_X210 | 1.562, 1.106, 0.3313 | True | False | False |

S2/A5三个初态均具有dt/cache完整对照；S4/S6新种子完成原逐步Gate，但未对每条补充dt/cache，故证据层级不同。S2_r0/r1/r2仍只有3个独立配置，不作总体概率或显著性推广。接触网络断开是实际输出，不删除失败种子。目标控制使用同一组三种子，其中r1/r2未参与幅值选择，但它们已用于机制确认，不冒称全新的外部测试集。

适用范围：`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED` · `training_ready=false`。Y 与四源静态结果另加 `HOMOGENEOUS_WATER_PROXY_SPECIFIC`。
