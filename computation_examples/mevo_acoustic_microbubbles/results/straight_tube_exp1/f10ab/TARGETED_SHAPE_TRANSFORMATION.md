# 预登记目标的实际结果

预登记目标：S2固定64粒子在2ms时长宽比增大3%±0.5个百分点，COM漂移≤20nm，终点仍为1个邻近分量。四个候选600/700/800/900Pa-eq，X相差180°。选择700Pa-eq后独立重算、两个机制确认种子及所有dt/cache对照。目标总体状态：**FAIL**。

| 初态 | 长宽比变化 / % | COM范数 / nm | 终点分量 | 目标量通过 | dt/cache通过 | 完整目标通过 |
| --- | --- | --- | --- | --- | --- | --- |
| S2_r0 | 2.864 | 0.006491 | 3 | False | True | False |
| S2_r1 | 2.481 | 0.003025 | 1 | False | True | False |
| S2_r2 | 3.817 | 0.01153 | 2 | False | True | False |

接近目标长宽比但破坏连通性仍判失败，不修改目标或删除不利粒子。未运行反向/恢复动作：`NOT_RUN`；首先解决团完整性比继续扩展控制更有价值。

3–5ms持续性试验（同一r0、1036.2Pa-eq，不是目标控制）:

| 窗口/ms | COM X/nm | Rg变化/% | 长宽比变化/% | 终点分量 |
| --- | --- | --- | --- | --- |
| 3.0 | -0.02293 | -1.786 | 9.189 | 4 |
| 5.0 | -0.03869 | -2.22 | 13.77 | 5 |

适用范围：`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED` · `training_ready=false`。Y 与四源静态结果另加 `HOMOGENEOUS_WATER_PROXY_SPECIFIC`。
