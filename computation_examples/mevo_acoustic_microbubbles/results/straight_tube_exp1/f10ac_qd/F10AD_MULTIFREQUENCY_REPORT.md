# F10AD 多频率科研报告

采用真实四圆片尺寸/位置/朝向的均匀水 FNM 声场，Catmull–Rom 缓存的梯度取自同一个压力插值多项式。频率只在完整轨迹之间变化。PDMS 反射、弹性波、声流、自由液面和硬件传递函数均未模拟。APECSS 为原 RPAR + 可逆 Marmottant，当前半径相关壳参数、0.036 N/m 初始表面张力假设及 0.002 Pa s 径向有效黏度不作实验标定解释。

| 频率 MHz | 声场/响应 | 完整动力学准入动作 | 所有测试动作通过 |
| --- | --- | --- | --- |
| 2.4 | True | X180, Y90, FOUR | True |
| 2.5 | True | X180, Y90, FOUR | True |
| 2.6 | True | X180, Y90, FOUR | True |


| 候选 | 频率 MHz | Δlog Rg | Δlog aspect | 新初态同号复现 | 全部保持连通 | 绝对主轴伸展超过噪声 |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 2.6 | -0.065757 | 0.095171 | True | False | False |
| 1 | 2.4 | -0.001762 | 0.004319 | True | True | False |
| 2 | 2.5 | -0.076256 | 0.079426 | True | False | False |
| OD | 2.5 | -0.286524 | -0.197027 | False | True | False |


三条不完整轨迹均为历史低压表域拒判：最后接受态独立 CPU 重建得到某泡约 0.099987 Pa，低于历史下限 0.1 Pa，未裁剪或外推。该重建不是当时失败的 GPU 力记录。

证据：[完整实验索引](ALL_RUNS.json)、[原始逐步独立审计](AUDIT.json)、[数值复算](VALIDATION_COMBINED.json)、[频率准入](INTERPRETATION.json)、[运行预算](PREREGISTRATION.json)。

`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `training_ready=false`
