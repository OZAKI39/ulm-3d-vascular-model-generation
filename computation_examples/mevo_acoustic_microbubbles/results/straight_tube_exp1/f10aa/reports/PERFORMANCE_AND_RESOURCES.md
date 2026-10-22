# 真实资源与时钟

AMD Ryzen 9 7950X3D，16 核/32逻辑线程，约63 GiB RAM；一张 RTX4090 24 GiB，驱动550.120；Python3.12.3，NumPy2.2.6，SciPy1.15.3，CuPy13.6，Warp1.18.0+cu12，实际 BLAS 为两套 OpenBLAS LP64/ILP64。版本/CPU亲和/线程和 GPU 空闲预检保存在 baseline/ENVIRONMENT。

完整正式 campaign 采样 CPU 均值 1.000 核当量，GPU 利用率峰 28.0%，显存峰 550.0 MiB，进程 RSS 峰 1042.2 MiB。1秒采样会漏短峰；CPU 均值包含分析和保存。GPU 开跑无其他计算任务，但服务器始终存在另一 CPU 任务，因此不是整机独占基准。各新后端在相同可用资源下对照，历史 F10Z 探针不作为加速分母。

已完成阶段墙钟（部分固定测试在中断 campaign 中，单列 Supervisor ledger 后才是总数）：

- baseline: 4.483 s
- benchmark_v3: 4.703 s
- full: 216.420 s
- parallel: 14.343 s
- extension: 315.999 s
- profile_extra: 5.219 s
- threaded: 37.543 s

上述完成阶段合计 598.711 s；额外失败安装、native 崩溃及固定状态阶段时间见 EXECUTION_LEDGER.json，不隐匿失败。绘图/视频/下载另列，不混为物理求解时间。

EXPLORATORY · PHYSICS_UNCALIBRATED · HOMOGENEOUS_WATER_CONTROL_PROXY · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
