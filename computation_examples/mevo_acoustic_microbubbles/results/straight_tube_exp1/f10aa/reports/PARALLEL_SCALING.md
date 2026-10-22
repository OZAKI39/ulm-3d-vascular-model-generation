# 独立任务并行吞吐

| 并发 | 批次 | 四条墙钟/s | 物理ms/墙钟s | 平均单轨迹/s |
|---:|---|---:|---:|---:|
|1|cold|3.3510|0.23873|0.5093|
|1|warm|2.6650|0.30019|0.5289|
|2|cold|2.3528|0.34002|0.5366|
|2|warm|1.5328|0.52191|0.5262|
|4|cold|1.8210|0.43932|0.5999|
|4|warm|0.9847|0.81246|0.5981|

相同 4 条 64 泡短窗口（0.2 ms）任务，分散/历史接触初态×OFF/PAIR，spawn 进程，1 BLAS 线程/worker，1/2/4 并发。每条内部原时序严格顺序，输出相同；不合并物理相互作用。

暖轮最高为 4 进程，相对 1 进程吞吐 2.706×。这是当前短任务下的真实结果，含每任务初始化和完整输出；不能宣传为单条轨迹同倍加速。冷轮包括进程/导入/CUDA 初始化；暖轮复用工作进程，但仍真实重建每条任务的引擎。资源监控 parent CPU 不包括 child，子进程 CPU/RSS 从各 RESULT 单列，GPU 为全卡采样。

EXPLORATORY · PHYSICS_UNCALIBRATED · HOMOGENEOUS_WATER_CONTROL_PROXY · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
