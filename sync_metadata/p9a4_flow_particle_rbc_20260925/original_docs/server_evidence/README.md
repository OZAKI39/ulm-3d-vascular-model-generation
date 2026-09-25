# 服务器源码与执行证据

各子目录保留原 `/workspace/<目录名>/` 下当前仍存在的源码、配置和选定日志；Git 会复用相同内容的对象。这些历史证据不作为默认开发代码入口。

- `particle_network_flow_mb_validation_v1_20260924T210047Z/`：当前 30 泡配对验证、源快照和执行记录。
- `flow_mean_2p0_mmps_A_H0_20260924T181140Z/`：当前 H0 求解器入口、输入与实际运行证据。
- `flow_mean_2p0_mmps_A_H0_TaylorHood_20260924T194503Z/`：已停止的 P2/P1 验证代码与资源/失败记录。
- `formal_3D_flow_solver_FEM_SimVascular_lzy/`：保留的 Stage L/N/Q 工具链脚本、配置、日志和版本来源。
- `pinned_solver_source/`：实际运行的 svMultiPhysics 修改文件、上游 commit 和完整 tracked diff。

第三方安装目录、CUDA/PETSc/MPI 二进制、解释器缓存、重复部署包及大中间结果未上传。历史脚本中的服务器路径保持原样；重建和重新运行时应按当前机器配置迁移。
已在代码瘦身中删除的 371 个服务器脚本没有恢复到本分支。
