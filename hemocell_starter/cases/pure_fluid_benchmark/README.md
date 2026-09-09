本例由固定 HemoCell examples/cube 和 examples/pipeflow 的初始化与纯流体路径派生。
保留上游头部及 AGPL-3.0-or-later 许可；Palabos 同为 AGPL-3.0-or-later。
核心求解代码未修改。只调用 HemoCell 管理的 D3Q19 Guo BGK lattice 的 collideAndStream。
没有 cellfield、addCellType、loadParticles、RBC.pos 或 PLT.pos，也没有细胞回调。

坐标为 (i+0.5)dx，周期长度 N dx；两半力严格相反。采样只遍历本 rank unique bulk，
MPI_SUM 汇总计数/矩，MPI_MAX 汇总各 rank 耗时；绝不采样 halo。
computeVelocity 已含外力半步校正。另输出原始动量速度，供独立检查。
开启力时由静止 populations 减去半步动量，保证宏观初速度连续为零。
无力准备阶段和开启体力后的流动阶段在同一个 HemoCell 实例中连续执行。

输出 profiles.csv、timings.csv、completion.json。CSV 是基本流体观测输出，
不是 HemoCell 细胞 HDF5；并行 HDF5 已在 CMake 找到 C 与 HL 并链接。
正式运行由 mirheo_starter 的共同 Python 入口生成独立 config.xml 并调度。
