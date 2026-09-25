# 服务器源码与执行证据

来源 `root@50.115.148.16:4159`。本目录保留前一版服务器快照，并纳入当前额外代码、配置、日志和关键科学数据。原始服务器目录与源文件没有搬移。

最新P9-A.4服务器目录为 `/workspace/particle9a4_population_inlet_20260924T233708Z`；其源码、154个已回传输出与本次root `particle_3d/`内容相同时，使用manifest引用而不复制。新增部署manifest位于对应服务器证据目录。原始时间戳和绝对路径属于来源记录。

H0求解、Network配对验证、P8/P9历史运行、shear-lift，以及历史HemoCell RBC、LAMMPS/Palabos资料均已检查。历史模拟的代码/日志保留不等于对其物理模型作新的审核。HemoCell/LAMMPS的原生大状态、完整第三方源码/工具链不全量打包。

本轮逐文件映射、去重引用与排除原因：[server_delta_inventory.json.gz](../sync_metadata/p9a4_flow_particle_rbc_20260925/server_delta_inventory.json.gz)；[汇总](../sync_metadata/p9a4_flow_particle_rbc_20260925/server_delta_summary.json)。`REPRESENTED_BY_EXISTING_ARTIFACT`指相同SHA内容已在仓库其他路径；representation明确普通字节或gzip原始字节。所有此类引用在整理后校验。`EXCLUDED_FROM_DELTA`只表示本次不新增，不删除父提交已经保存的历史证据。

固定svMultiPhysics补丁及实际源文件仍在 `pinned_solver_source/`；完整当前源与许可证见根目录`vendor/`。CUDA/PETSc/MPI安装、二进制和Python环境留在服务器，按构建/版本记录复现。
