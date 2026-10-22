# F10Z-Simulation

已完成：102 条实际积分，25 组 dt/cache 连续量对照，最长 5 ms；中等分散态有显著内部收缩与局部新邻近，但完整 XY 位置—形状世界模型尚未准入。完整回答见 [F10Z_DEVELOPMENT_REPORT.md](F10Z_DEVELOPMENT_REPORT.md)。

研究分散32泡在既有四PZT均匀水代理中的平移、相对靠近与形变。频率2.5 MHz，半径1.78 µm，背景流严格为零；幅值是Pa-equivalent，不是电压。

保留 EXPLORATORY / PHYSICS_UNCALIBRATED / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED，training_ready=false。

程序位于 `scripts/straight_tube_exp1/f10z/`，协议在 `configs/straight_tube_exp1/f10z/protocol.json`。`PREREGISTRATION.json`、`INITIAL_STATES.h5` 与 `RUN_SCHEDULES.json` 是运行前冻结输入，不可覆盖。所有物理积分复用 F10S Campaign.run 和原 WaterSequential；新代码只适配输入、调度、导出、分析和图表。

在工程根目录使用独立本地 `.venv`：

```bash
PYTHONPATH=src:. .venv/bin/python -m pytest -q tests/f10z
PYTHONPATH=src:. .venv/bin/python -m scripts.straight_tube_exp1.f10z.analyze
```

远端计算复用现有 mevo 环境，每次使用不可变 Git 快照和全新 Supervisor 标签。请先读协议与当前结果，不能重复运行 prepare 覆盖初态。运行入口：`python -m scripts.straight_tube_exp1.f10z.remote --stage A|B|C --label NEW_UNIQUE_LABEL`；C须通过前级数值和成本门槛。

`dataset/*.h5` 由真实已接受状态导出；失败后不填补时间。`SIMULATION_XY_PROJECTION` 不是已校准的显微图像或 Chronos 检测输出。几何连边不代表融合。


后处理复现顺序（工程根目录）：

```bash
export PYTHONPATH=src:.
.venv/bin/python -m scripts.straight_tube_exp1.f10z.analyze
.venv/bin/python -m scripts.straight_tube_exp1.f10z.enrich
.venv/bin/python -m scripts.straight_tube_exp1.f10z.quality
.venv/bin/python -m scripts.straight_tube_exp1.f10z.figures
.venv/bin/python -m scripts.straight_tube_exp1.f10z.videos
.venv/bin/python -m scripts.straight_tube_exp1.f10z.finalize audit
.venv/bin/python -m scripts.straight_tube_exp1.f10z.deliver
```

图表：figures/；视频：videos/。原始数据：a_formal_v1、b_refinement_v1、c_extension_v1 的 runs/；统一状态：dataset/。数据路径、执行 Git 身份、环境版本和许可证在 evidence/。绘图与仿真运行身份分别保存，不将最终文档提交假称为早先执行源码。

没有操作硬件，没有训练神经世界模型或运行 MPPI，没有 GitHub 推送。原始轨迹、大型派生 HDF5 和虚拟环境保留在本地/原计算服务器；Git 选择纳入复现代码、固定初态/配置、图表、三段小体积科研视频、小日志、报告和哈希。
