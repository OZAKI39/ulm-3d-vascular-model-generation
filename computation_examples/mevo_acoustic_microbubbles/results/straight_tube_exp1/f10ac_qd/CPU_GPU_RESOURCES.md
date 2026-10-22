# 资源与真实计算成本

| 阶段 | 墙钟 s |
| --- | --- |
| pilot | 20.724 |
| frequency_build | 41.891 |
| corrected_comparison | 457.39 |
| initial_integration_audit | 487.36 |
| dynamics | 122.531 |
| elite_validation | 98.302 |
| independent_audit | 256.767 |
| extension_winner | 83.066 |
| extension_frequency_2400000 | 48.735 |
| extension_frequency_2600000 | 52.928 |
| initial_interpretation_media_review | 28.862 |
| final_figures_caption_review_tests | 22.17 |
| local_component_postprocessing | 1.106 |
| outside_archive_validation | 109.886 |

已记录阶段合计 1831.719 s；真实逐泡尝试 366，完整 363。并行进程的单轨迹墙钟不能与阶段墙钟相加。

RTX4090 执行 Warp 声学，SciPy CPU 执行原接触求解；预检 1/2/4 worker 墙钟分别 3.969 s, 2.876 s, 2.298 s，选择 4 worker。优化器 GP 在 CPU，未为低维 GP 引入大型 GPU 框架。

| 方法 | 计费评价 | 新增真实调用 | 完整 | 档案外 | 重复建议 | 档案格 | QD score | 计费墙钟 s | 优化器 s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SOBOL | 32 | 24 | 31 | 16 | 0 | 5 | 492.543 | 111.864 | 0.002 |
| BOP | 32 | 24 | 30 | 25 | 1 | 4 | 392.41 | 301.04 | 10.044 |
| CMA_MAE | 32 | 24 | 31 | 9 | 0 | 6 | 593.517 | 88.975 | 1.962 |


`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `training_ready=false`
