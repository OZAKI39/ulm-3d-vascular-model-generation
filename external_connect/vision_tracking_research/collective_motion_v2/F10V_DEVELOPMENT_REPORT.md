# F10V-Vision 科研开发报告

本轮完成了可执行的学习型群体运动观测原型、12项实际训练、验证选模、独立校准、1536样本最终测试和真实录像因果回放。**科研结论为不替换PIV**。代码运行完成与科学准确准入严格分开。

| Method | MAE denominator | Six-parameter MAE (vx, vy, div, vort, s1, s2) | Accepted / 1536 | Failed accepted |
|---|---:|---|---:|---:|
| PIV single pass | 1114 | 1.894, 1.698, 0.020, 0.039, 0.009, 0.010 | 1111 | 85 |
| Direct6 gray | 1320 | 11.420, 7.994, 0.057, 0.327, 0.045, 0.030 | 0 | 0 |
| Direct6 + heat | 1320 | 11.394, 8.011, 0.061, 0.332, 0.040, 0.031 | 0 | 0 |
| Direct6 heat only | 1320 | 13.416, 9.526, 0.059, 0.300, 0.044, 0.031 | 0 | 0 |
| Direct6 + heat robust | 1320 | 11.554, 7.937, 0.057, 0.288, 0.043, 0.026 | 0 | 0 |
| GeoMotion | 1320 | 10.833, 7.597, 0.057, 0.282, 0.040, 0.023 | 0 | 0 |

| Method | Error risk at 25% diagnostic coverage (384/1536) |
|---|---:|
| PIV single pass | 0.52% |
| Direct6 gray | 80.73% |
| Direct6 + heat | 86.98% |
| Direct6 heat only | 83.07% |
| Direct6 + heat robust | 79.43% |
| GeoMotion | 86.46% |

| Method | Pure model p50/p95 ms | Full pipeline p50/p95/max ms | Consumed /120 | Fresh six-parameter Hz |
|---|---:|---:|---:|---:|
| PIV single pass | N/A | 34.35, 44.29, 47.32 | 114 | 25.24 |
| Direct6 gray | 2.45, 3.25 | 19.72, 25.16, 32.23 | 120 | 0.00 |
| Direct6 + heat | 3.24, 4.39 | 16.08, 21.13, 22.98 | 120 | 0.00 |
| GeoMotion | 2.94, 3.42 | 18.92, 23.56, 26.06 | 120 | 0.00 |

训练共258个epoch、11094优化步、11.75分钟；不是仅烟雾测试。最佳权重按冻结验证目标选择，两个GeoMotion主模型最佳点均在epoch1，后续验证恶化的证据没有删除。最终测试没有用于重新训练、调阈值或挑选随机种子。

所有AI因未找到合法验证风险阈值而拒判；零接受不是零错误成功。PIV的1111个接受中85个失败也完整保留。七组图、三段视频展示了实际数据与失败。只有一个录制会话、没有真实六参数真值；EXPLORATORY，PHYSICS_UNCALIBRATED，training_ready=false。相机、PZT、AFG及物理动作世界模型均未操作或训练。

30项逐问回答见[FINAL_RECOMMENDATION.md](FINAL_RECOMMENDATION.md)，分类/条件/种子/素材统计见[results/analysis.json](results/analysis.json)，源码入口与恢复方法见[ACTIVE_CODE_INDEX.md](ACTIVE_CODE_INDEX.md)和[README.md](README.md)。
