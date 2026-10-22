# 实际性能与资源

远端科学批次墙钟 1287.478 s（21.458 min），保存的积分时间和 1134.886 s。本地分析 122.252 s；后处理/视频耗时另见 evidence/DELIVERY_TIMING.json。两者不是“科研开发总用时”。
 
RTX4090 24 GB + Ryzen9 7950X3D 32 logical CPU，复用既有环境，未停止其他 CPU 项目。CPU user 1283.683 s，system 5.537 s；有效 B/C 监测 GPU 平均 2.02%、最大 5.0%，显存峰 468.0 MiB，进程 RSS 峰 885.57 MiB。
 
逐步计时总和 1116.685 s，QP 726.390 s、阻力 113.568 s，嵌套时间不能相加声称独立成本。各批 profile.txt/pstats 保留声场查询、F2、矩阵、QP、Python、I/O 累计调用；计时含 profiling 开销，不是生产性能峰值。单条只有 32 泡，CPU 接触和小规模 GPU 查询混合，GPU 利用率低不表示算例错误。
 
A 批监控线程因为旧环境没有 psutil 退出，积分主线程正常完成；A 的资源时间序列缺失，明确保留报错。B/C 改为标准库 /proc 监控，不改物理计算，也未为监控重跑所有物理。环境 lock、版本/许可证清单、Supervisor 标签和不可变执行源码见 evidence/、logs/。
