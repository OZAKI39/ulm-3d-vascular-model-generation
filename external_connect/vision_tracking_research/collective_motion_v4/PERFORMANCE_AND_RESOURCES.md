# 性能和资源

同一服务器AMD Ryzen 9 7950X3D、16物理/32逻辑核、约62GiB RAM；RTX4090 24564MiB，driver550.120。环境独立，PyTorch2.14.0+cu126、CuPy13.6.0、numpy2.2.6、SciPy1.17.1、OpenPIV0.26.1及官方Rust编译扩展。版本声明与真实执行均保存。未修改驱动，未停止其他CPU工程；GPU基准期间没有其他计算GPU任务。

| Run | Frames | p50 / p95 / max ms | Mean valid windows | Fresh group-frame Hz | Drops |
|---|---|---|---|---|---|
| original_offline_regular | 937 | 58.02 / 77.67 / 254.12 | 490.8 | offline | 0 |
| original_offline_sparse1 | 937 | 58.79 / 79.16 / 162.85 | 491.2 | offline | 0 |
| original_offline_sparse3 | 937 | 50.86 / 74.53 / 96.86 | 365.8 | offline | 0 |
| original_offline_sparse6 | 937 | 51.20 / 68.77 / 119.93 | 334.1 | offline | 0 |
| verified_offline_regular | 937 | 51.80 / 70.76 / 209.11 | 490.8 | offline | 0 |
| verified_replay_regular | 778 | 28.98 / 71.43 / 94.30 | 218.8 | 10.42 | 159 |
| verified_replay_sparse3 | 791 | 26.53 / 65.82 / 128.92 | 179.9 | 11.47 | 146 |
| verified_replay_sparse6 | 803 | 30.38 / 63.50 / 114.05 | 174.2 | 12.69 | 134 |


正式完整离线处理含decode、全扫描、PIV、验证、局部群、JSON序列化；gzip科研存档写盘另计wall_s。回放计时还含实际M4本机消费者发送，单帧decode成本加回，不含真实相机链路。完整regular离线p95=70.76ms；paced p95=71.43ms、max=94.30ms。33.3ms和16.7ms目标均未准入。

规则回放CPU核心等价 1.11，进程RSS峰值 867.2MiB；设备显存采样峰值 790MiB，GPU平均利用率 2.2%、采样最大 5%。这是约0.2s采样和设备级显存，不是CUDA allocator瞬态绝对峰值。完整raw资源样本保留。

离线扫描kernel的1000/平均耗时约 111.1Hz，但完整paced实际全通道扫描率仅 24.93Hz。新鲜有效群帧 10.42Hz；空间群状态总计 145.34/s，不是这么多个物理团。主要剩余成本为CPU局部群拟合、候选/序列化；GPUFFT已较快。内部profiler保存FFT/峰/SNR函数耗时，profile样本不混入warm性能统计。

软件环境准备、数据传输、图表编码不能算GPU科学吞吐。各实际job/identity/resources/wall_s和日志在results、evidence/deployment、logs/remote；完整原片SHA固定。
