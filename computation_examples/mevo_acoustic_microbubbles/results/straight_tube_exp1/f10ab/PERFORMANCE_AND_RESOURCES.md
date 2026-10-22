# 实际计算成本

| 阶段 | 完整/总轨迹 | 物理阶段墙钟/s |
| --- | --- | --- |
| pilot | 30/35 | 33.9 |
| discovery | 25/25 | 105.3 |
| refinement | 14/14 | 72.72 |
| confirmation | 22/22 | 72.21 |
| sensitivity | 16/16 | 102 |
| sensitivity_checks | 16/16 | 156 |
| target_search | 4/4 | 10.25 |
| extension | 6/6 | 48.71 |
| target_confirmation | 9/9 | 30.06 |

静态诊断 38.786s；物理阶段累计 669.987s，低于3600s预算。单条最长32.132s，低于240s；最多两条物理轨迹并行，单任务BLAS线程1。逐条wall之和1199.866s不能当并行墙钟。独立分析阶段之和405.476s；绘图/视频/打包耗时另见evidence/*TIMING.json和videos/MANIFEST.json，传输与环境安装不混入求解速度。峰值worker RSS893.7MiB。

使用RTX4090 Warp声学+CPU接触，未宣称服务器独占；其他已有任务未停止。独立.postvenv只用于绘图/视频，物理.venv未修改。首次后处理安装因Python3.11不满足av19失败，改用现有Python3.12成功，没有更换求解依赖。后处理首次导入另缺PyYAML，已仅在.postvenv补齐6.0.3；失败日志保留。

适用范围：`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED` · `training_ready=false`。Y 与四源静态结果另加 `HOMOGENEOUS_WATER_PROXY_SPECIFIC`。
