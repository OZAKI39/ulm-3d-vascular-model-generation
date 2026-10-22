# 原始瓶颈及计时定义

| 状态 | 稠密接触几何/ms | 阻力/ms | Cholesky/ms | 线性求解/ms | SQP/ms |
|---|---:|---:|---:|---:|---:|
|N32_dispersed|1.5822|0.1896|0.0286|0.0387|1.8018|
|N64_historical_end|11.7532|0.5816|0.1008|12.3022|25.2513|
|N128_dispersed|28.9434|1.7908|0.5356|0.3275|35.3910|

上述为嵌套 inclusive instrumented 时间，不能相加。独立预热单步 p50/95 见原始 PROFILING.json；instrumentation、cold engine、wall、process CPU 分开。排序：反复稠密接触几何；进入接触后的全 B.T 求解；阻力/数组工作；Cholesky。128 泡全部 8128 对、384 列，B 的存储和遍历随 N³ 增长，全右端项求解还会进一步放大开销。

CPU 使用 cProfile、line_profiler、tracemalloc；GPU 为实际 Warp CUDA event kernel 计时，fused 声场/APECSS 无法诚实拆成两个独立内核。F2 pair/reduction 和双墙面 kernel 分别保存。D2H wrapper 计时包含同步等待与复制，不能当成纯 PCIe 带宽。`benchmark_v3/FULL_STEP.json` 补充原始 profiling print=False 未触发回调的缺项；错误未冒充 GPU 测量。Nsight Systems 不可用，未运行内核权限分析器。

Python 分配峰值与进程 RSS、实际线程数、BLAS 线程、数组上传下载、CuPy 额外流量见 `profile_extra`；所有原始 function/line 计时和库配置均保留。

EXPLORATORY · PHYSICS_UNCALIBRATED · HOMOGENEOUS_WATER_CONTROL_PROXY · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
