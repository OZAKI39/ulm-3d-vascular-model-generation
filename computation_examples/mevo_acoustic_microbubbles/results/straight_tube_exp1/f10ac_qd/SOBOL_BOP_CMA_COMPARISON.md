# 统一条件算法比较

| 方法 | 计费评价 | 新增真实调用 | 完整 | 档案外 | 重复建议 | 档案格 | QD score | 计费墙钟 s | 优化器 s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SOBOL | 32 | 24 | 31 | 16 | 0 | 5 | 492.543 | 111.864 | 0.002 |
| BOP | 32 | 24 | 30 | 25 | 1 | 4 | 392.41 | 301.04 | 10.044 |
| CMA_MAE | 32 | 24 | 31 | 9 | 0 | 6 | 593.517 | 88.975 | 1.962 |

各方法 32 次计费评价包含同一组 8 次 Sobol 初始数据；共同样本只物理计算一次，向每种方法计入同样的公共墙钟。每种方法另外实际计算 24 次；失败及越界结果均占预算。

档案范围始终保持不变；额外对一条完整但越界的负长宽比响应做了单独预登记及三个全新种子复核（OUTSIDE_PREREGISTRATION.json），不计为档案格或优化器收益。初次集成比较保存在 `search/` 和 `COMPARISON.json`：档案外真实结果也被排除出 GP 学习，产生重复探索；它不作为最终算法优劣依据。修复后重新完整运行 `search_v2/`，保持档案范围、动作范围、32 次计费评价和初态不变。前一次计算也计入总成本。

直接复用 pyribs 0.12.0 的 `BayesianOptimizationEmitter` / `BayesianOptimizationScheduler`（BOP-Elites）和 `EvolutionStrategyEmitter` / `Scheduler`（CMA-MAE），以及 SciPy 1.15.3 `qmc.Sobol`。未实现自定义 GP、EJIE、CMA 更新或档案阈值更新。唯一档案边界适配禁止 pyribs 默认将越界指标放进边界格，真实完整观测仍反馈给官方优化器。物理失败无完整形态标签：BOP 不学习伪标签，CMA 含失败批次通过官方 emitter 重新初始化；原始轨迹全部保留。

`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `training_ready=false`
