先读 PASSIVE_MICROBUBBLE_TRANSPORT_V0_REPORT.md 与 FINAL_TERMINAL_SUMMARY.txt。

公开独立复核：`/usr/bin/python3 -B finalize_passive_transport_v0.py --root "$PWD"`。它会重新生成审计与图件；保留封存包时请在副本运行。

CMakeLists.txt 使用只读的上一阶段 LAMMPS build_cpu/build_gpu 静态库，以 ENGINE_ROOT 参数指定原构建目录。当前 CUDA 路径对应本 Vast 的 13.2 环境，迁往其他平台需重新配置。binaries/ 是本次实际构建证据，保留原 build-tree runtime 路径，并非通用可搬移二进制。

正式 frozen contract 位于 contracts/PASSIVE_MICROBUBBLE_TRANSPORT_V0_CONTRACT.json；派生所选步长位于 PASSIVE_TRANSPORT_TIMESTEP_CONTRACT.json。case0/1/2 选中档位的已完成轨迹直接用作 production V0 证据，没有额外重复求解。其它未选中的生产输入保留但未执行。

本任务程序使用 HWLOC_COMPONENTS=-gl、OMP_NUM_THREADS=1。-gl 仅跳过 GL 显示拓扑探测。scripts/run_validation.py 是已完成的一次性执行记录，会拒绝重复运行已有算例。

ParaView 同时打开 visualization/FROZEN_LUMEN.vtp 和 CASE2_TRAJECTORY.vtp / CASE5_TRAJECTORIES.vtp。坐标单位 m。HUMAN_VISUAL_REVIEW=PENDING。
