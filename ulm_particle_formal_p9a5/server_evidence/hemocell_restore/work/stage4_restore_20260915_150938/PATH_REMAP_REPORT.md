# 路径迁移记录

归档保持只读；仅在新工作副本中解析原/workspace命名空间并复制真实目标文件。归档symlink不改写，工作副本中的source、Palabos、输入均逐文件保持SHA256一致。完整OLD_PATH→NEW_PATH见RESTORE_INPUT_MANIFEST.tsv及WORKING_COPY_SHA256.json；前者记录archive源路径，后者记录work路径及相同SHA。

| OLD_PATH | NEW_PATH（相对本work目录） | REASON |
|---|---|---|
| /workspace/hemocell_gpu_poc/20260914_stage4_batch_dispatch/source | source/ | 原冻结C++与CMake原样保存 |
| 原Stage4 upstream/palabos | palabos/ | 固定导出源码，Stage4 patch已应用，不重复打补丁 |
| 原PoC upstream/palabos | official_palabos/ | 官方cavity及200步硬限原样复制 |
| 原PoC frozen_bundle | frozen_bundle/ | 几何/物理/监控原样复制 |
| 原NVHPC toolchains/nvhpc | /workspace/hemocell_restore/toolchains/nvhpc_26_5 | 本实例独立工具链安装 |
| 原run目录 | stage4_200、stage4_1000、stage4_5000 | 每次从零，使用原对应run的合同字节 |
| 旧build_gpu | build/ | 本平台重新编译，原binary只作provenance |

runtime参数文件直接复制原已通过的RTX5090 run，不做文本或数值替换。参数首两行本来就是相对frozen_bundle的路径。历史参数文件max_steps=700000原样保留，实际Stage4 C++外部硬限仅允许200/1000/5000，本任务不会执行长程。

新独立CMake配方只把路径与GPU架构cc120→实测cc89适配；所有C++及原CMake原样封存。OMP1、MPI1、STAGE4_LBM=batch、STAGE2_PROFILE_MODE=off沿用原已验证模式。
