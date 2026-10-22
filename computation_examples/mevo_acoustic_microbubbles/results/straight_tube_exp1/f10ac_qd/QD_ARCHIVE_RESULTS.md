# 实际形态档案

| 方法 | 计费评价 | 新增真实调用 | 完整 | 档案外 | 重复建议 | 档案格 | QD score | 计费墙钟 s | 优化器 s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SOBOL | 32 | 24 | 31 | 16 | 0 | 5 | 492.543 | 111.864 | 0.002 |
| BOP | 32 | 24 | 30 | 25 | 1 | 4 | 392.41 | 301.04 | 10.044 |
| CMA_MAE | 32 | 24 | 31 | 9 | 0 | 6 | 593.517 | 88.975 | 1.962 |


| 候选 | 频率 MHz | Δlog Rg | Δlog aspect | 新初态同号复现 | 全部保持连通 | 绝对主轴伸展超过噪声 |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 2.6 | -0.065757 | 0.095171 | True | False | False |
| 1 | 2.4 | -0.001762 | 0.004319 | True | True | False |
| 2 | 2.5 | -0.076256 | 0.079426 | True | False | False |
| OD | 2.5 | -0.286524 | -0.197027 | False | True | False |


候选选择在独立初态复算之前冻结。结果索引见 ELITE_VALIDATION_PREREGISTRATION.json。三条不完整轨迹均为历史低压表域拒判：最后接受态独立 CPU 重建得到某泡约 0.099987 Pa，低于历史下限 0.1 Pa，未裁剪或外推。该重建不是当时失败的 GPU 力记录。

证据：[完整实验索引](ALL_RUNS.json)、[原始逐步独立审计](AUDIT.json)、[数值复算](VALIDATION_COMBINED.json)、[频率准入](INTERPRETATION.json)、[运行预算](PREREGISTRATION.json)。

`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `training_ready=false`
