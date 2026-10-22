# 多频率完整逐泡动力学

| 频率 MHz | 声场/响应 | 完整动力学准入动作 | 所有测试动作通过 |
| --- | --- | --- | --- |
| 2.4 | True | X180, Y90, FOUR | True |
| 2.5 | True | X180, Y90, FOUR | True |
| 2.6 | True | X180, Y90, FOUR | True |


N32/N64 的 0.5 ms OFF、单 X、X180、Y90、四路控制先运行；通过对应 dt/2、细缓存对照后选择性延长到 2 ms。绝不在轨迹中途拼接频率响应。实际轨迹与失败前缀均保留。三次表域拒判发生在新的声力求值阶段，没有积分器拒绝提案；旧 HDF5 中上次成功提案已在派生数据中注明，不能冒称本次失败提案。

| 候选 | 频率 MHz | Δlog Rg | Δlog aspect | 新初态同号复现 | 全部保持连通 | 绝对主轴伸展超过噪声 |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 2.6 | -0.065757 | 0.095171 | True | False | False |
| 1 | 2.4 | -0.001762 | 0.004319 | True | True | False |
| 2 | 2.5 | -0.076256 | 0.079426 | True | False | False |
| OD | 2.5 | -0.286524 | -0.197027 | False | True | False |


`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `training_ready=false`
