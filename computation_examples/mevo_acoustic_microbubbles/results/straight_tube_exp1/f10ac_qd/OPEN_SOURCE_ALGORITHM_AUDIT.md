# 成熟开源实现核查

直接复用 pyribs 0.12.0 的 `BayesianOptimizationEmitter` / `BayesianOptimizationScheduler`（BOP-Elites）和 `EvolutionStrategyEmitter` / `Scheduler`（CMA-MAE），以及 SciPy 1.15.3 `qmc.Sobol`。未实现自定义 GP、EJIE、CMA 更新或档案阈值更新。唯一档案边界适配禁止 pyribs 默认将越界指标放进边界格，真实完整观测仍反馈给官方优化器。物理失败无完整形态标签：BOP 不学习伪标签，CMA 含失败批次通过官方 emitter 重新初始化；原始轨迹全部保留。

[pyribs BOP 官方示例](https://docs.pyribs.org/en/stable/examples/bop_elites.html)；[CMA-MAE 官方教程](https://docs.pyribs.org/en/stable/tutorials/cma_mae.html)；[SciPy Sobol](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.qmc.Sobol.html)；[APECSS](https://github.com/polycfd/apecss)。版本与许可见 [DEPENDENCIES.json](evidence/DEPENDENCIES.json)，APECSS 固定提交 `a0ee071ca247b5489e997d90e83b5f53c40e4ec7`；构建与依赖冻结日志随交付保存。

档案范围始终保持不变；额外对一条完整但越界的负长宽比响应做了单独预登记及三个全新种子复核（OUTSIDE_PREREGISTRATION.json），不计为档案格或优化器收益。初次集成比较保存在 `search/` 和 `COMPARISON.json`：档案外真实结果也被排除出 GP 学习，产生重复探索；它不作为最终算法优劣依据。修复后重新完整运行 `search_v2/`，保持档案范围、动作范围、32 次计费评价和初态不变。前一次计算也计入总成本。

`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `training_ready=false`
