# F10AA：CPU/GPU 逐泡求解优化

本轮在新工作树保护 F10S/F10Z 历史代码、输入和结果，新增矩阵自由全粒子对接触、按需工作集线性求解、向量化润滑阻力和可选同流打包读回。没有改变物理公式、原 dt 或数值门槛。

- [中文报告与32项回答](../results/straight_tube_exp1/f10aa/reports/F10AA_DEVELOPMENT_REPORT.md)
- [活动代码与复现](../results/straight_tube_exp1/f10aa/reports/ACTIVE_CODE_INDEX.md)
- [最佳配置与限制](../results/straight_tube_exp1/f10aa/reports/FINAL_RECOMMENDATION.md)
- [机器Gate](../results/straight_tube_exp1/f10aa/FINAL_GATE.json)
- [所有产物SHA](../results/straight_tube_exp1/f10aa/INPUT_OUTPUT_SHA256.json)

64泡2ms：原CPU34.632 s，优化单线程5.760 s，6.012×；128泡1ms：原CPU44.374 s，优化单线程6.561 s，6.763×，启动前固定4线程5.293 s，8.384×。吞吐与单轨迹延迟分开；4个独立进程暖轮吞吐2.706×。

59条轨迹尝试58条完成；128泡2ms原CPU参考在1.659ms因240s预算停止，原始FAIL保留。优化版2ms、dt/2及细缓存完成，但不能宣称完整2ms参考已通过。256泡仅合法短步预检。

主要原始数据保留在本工作树及服务器 `/workspace/mevo_acoustic_microbubbles/f10aa_simulation_20261018/results/straight_tube_exp1/f10aa`；本地完整结果按SHA核对。Git提交保存源代码、冻结小输入、全部小型数值证据、报告、图和视频；环境、大型HDF5与历史外部缓存通过路径及哈希管理。

`F10AA_BLAS_THREADS=4 bash scripts/straight_tube_exp1/f10aa/python.sh ...` 在导入前固定线程。默认1线程；不要在运行中切换到4线程（已复现native故障）。32/64使用1线程，可选WaterPacked；128单任务可用启动前4线程，独立任务并行每进程1线程。

EXPLORATORY / PHYSICS_UNCALIBRATED / HOMOGENEOUS_WATER_CONTROL_PROXY / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED / training_ready=false。
