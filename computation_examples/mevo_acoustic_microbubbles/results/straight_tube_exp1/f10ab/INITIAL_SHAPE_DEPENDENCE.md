# 初态依赖

| 形状 | 动作 | 三初态长宽比变化 / % | 方向一致 | 全部终点连通 | 全种子dt/cache复核 |
| --- | --- | --- | --- | --- | --- |
| S2 | A5_X180 | 6.275, 5.621, 8.291 | True | False | True |
| S2 | A7_X210 | -3.71, -3.757, 6.572 | False | True | False |
| S4 | A5_X180 | 5.759, 4.961, 4.811 | True | False | False |
| S4 | A7_X210 | 5.403, -1.642, 5.48 | False | True | False |
| S6 | A5_X180 | 0.2097, -0.09372, -0.09834 | False | False | False |
| S6 | A7_X210 | 1.562, 1.106, 0.3313 | True | False | False |

A5的S2/S4长宽比增大方向重复；A7在S2/S4间甚至同类不同种子中方向变化。S6是分散粒子集，不把其统计形变称为已形成完整泡团变形。

全部发现响应：[矩阵](SHAPE_RESPONSE_MATRIX.md)。

适用范围：`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED` · `training_ready=false`。Y 与四源静态结果另加 `HOMOGENEOUS_WATER_PROXY_SPECIFIC`。
