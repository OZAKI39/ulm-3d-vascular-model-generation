# F10AC-QD 科研开发报告

已完成 F10AC-QD + F10AD 的有限真实实验：366 次逐泡尝试、363 条完整轨迹，独立审计 697500 步。最终公平比较选定 CMA_MAE（仅本次预登记范围）。

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


| 频率 MHz | 声场/响应 | 完整动力学准入动作 | 所有测试动作通过 |
| --- | --- | --- | --- |
| 2.4 | True | X180, Y90, FOUR | True |
| 2.5 | True | X180, Y90, FOUR | True |
| 2.6 | True | X180, Y90, FOUR | True |


三条不完整轨迹均为历史低压表域拒判：最后接受态独立 CPU 重建得到某泡约 0.099987 Pa，低于历史下限 0.1 Pa，未裁剪或外推。该重建不是当时失败的 GPU 力记录。

证据：[完整实验索引](ALL_RUNS.json)、[原始逐步独立审计](AUDIT.json)、[数值复算](VALIDATION_COMBINED.json)、[频率准入](INTERPRETATION.json)、[运行预算](PREREGISTRATION.json)。

档案范围始终保持不变；额外对一条完整但越界的负长宽比响应做了单独预登记及三个全新种子复核（OUTSIDE_PREREGISTRATION.json），不计为档案格或优化器收益。初次集成比较保存在 `search/` 和 `COMPARISON.json`：档案外真实结果也被排除出 GP 学习，产生重复探索；它不作为最终算法优劣依据。修复后重新完整运行 `search_v2/`，保持档案范围、动作范围、32 次计费评价和初态不变。前一次计算也计入总成本。

记录阶段总耗时 1831.719 s，包含首次集成审计、物理生成、真实搜索、独立复算和后处理；不把重用共同样本算成新增物理调用。十二张科研图在 figures/，三段真实轨迹视频在 videos/；完整原始逐步数据在服务器，便携真值数据在 dataset/。

完整中文逐项回答见 [FINAL_CHINESE_REPORT.md](FINAL_CHINESE_REPORT.md)。

`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `training_ready=false`
