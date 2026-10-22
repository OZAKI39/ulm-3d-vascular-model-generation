# 三个新初态与完整数值复算

搜索初态为 F10AB 的 `S2_r0`，64 个半径 1.78 µm 的固定 ID 微泡，几何合法、初始邻接连通；它是人工指定的演示初态，不是本轮声学自发聚集得到的初态。三个档案候选的新增验证种子为 301801/301802/301803；档案外候选另用 401801/401802/401803。搜索及各自候选选择均未使用对应验证初态。N32 用于历史回归和新频率短轨迹准入。背景速度严格为零。

| 候选 | 频率 MHz | Δlog Rg | Δlog aspect | 新初态同号复现 | 全部保持连通 | 绝对主轴伸展超过噪声 |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 2.6 | -0.065757 | 0.095171 | True | False | False |
| 1 | 2.4 | -0.001762 | 0.004319 | True | True | False |
| 2 | 2.5 | -0.076256 | 0.079426 | True | False | False |
| OD | 2.5 | -0.286524 | -0.197027 | False | True | False |


三条不完整轨迹均为历史低压表域拒判：最后接受态独立 CPU 重建得到某泡约 0.099987 Pa，低于历史下限 0.1 Pa，未裁剪或外推。该重建不是当时失败的 GPU 力记录。

证据：[完整实验索引](ALL_RUNS.json)、[原始逐步独立审计](AUDIT.json)、[数值复算](VALIDATION_COMBINED.json)、[频率准入](INTERPRETATION.json)、[运行预算](PREREGISTRATION.json)。

各候选在搜索初态和三个新初态均尝试 DT1、DT05、FINE；同号复现需要四个初态均超过 max(1 nm,5×数值残差差异)，描述量亦用 max(1e-4,5×差异) 判断可分辨性。

`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `training_ready=false`
