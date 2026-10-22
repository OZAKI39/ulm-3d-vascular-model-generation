# F10AA 交付索引

本轮保留原物理和数值门槛，以矩阵自由全粒子对接触＋按需活跃列求解实现完整轨迹加速。32/64推荐单线程CPU力学，可选packed声场；128单任务推荐启动前固定4线程CPU力学，原Warp声学/壁面不变。原引擎及历史结果未覆盖。

| N | 物理时间 | 原CPU完整墙钟/s | 优化CPU/s | Packed/s | CuPy/s | CPU加速 |
|---|---:|---:|---:|---:|---:|---:|
|32|2 ms|7.7872|2.6455|2.4808|4.6503|2.944×|
|64|2 ms|34.6319|5.7603|5.4442|7.7861|6.012×|
|128|1 ms|44.3742|6.5612|6.5523|7.8433|6.763×|


图：`../figures/`，8 组 PNG/PDF；视频：`../videos/`，2 段 MP4，全部帧已解码。数据：`../full/`、`../extension/`、`../parallel/`。审核入口：`F10AA_DEVELOPMENT_REPORT.md`、`FINAL_RECOMMENDATION.md`、`../FINAL_GATE.json`。

新分支 `dev/f10aa-solver-acceleration-20261018`。每次执行的不可变 commit 和 Supervisor 标签见 `../logs/*_EXECUTION.json`，基线 `df76b1e9d311024398bc2bca40b6cf0f6382ba64`。科学计算源码快照与图表源码身份分别保存，不能把后期绘图提交当成全部早期执行身份。

启动器：`F10AA_BLAS_THREADS=4 bash scripts/straight_tube_exp1/f10aa/python.sh ...` 在Python导入前固定线程数，避免运行时resize。32/64或独立任务并行用1线程。

复现：项目根目录先设置 `PYTHONPATH=src:.` 和 `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`。GPU 环境使用 `/workspace/mevo_acoustic_microbubbles/f10aa_simulation_20261018/.venv/bin/python`；执行 `-m scripts.straight_tube_exp1.f10aa.benchmark --output benchmark_replay` 可独立重做固定状态和完整步基准。完整轨迹入口 `full.py`，需在新输出目录/独立 checkout 复现，既有目录存在时主动拒绝覆盖。后续 `parallel.py`、`extension.py`、`analyze.py`、`figures.py`、`videos.py`、`report.py`。所有历史原始缓存按原 SHA 恢复；没有凭空生成替代缓存。

EXPLORATORY · PHYSICS_UNCALIBRATED · HOMOGENEOUS_WATER_CONTROL_PROXY · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
