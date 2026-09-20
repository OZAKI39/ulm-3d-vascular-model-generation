# Frozen Stage Q 源码来源

svMultiPhysics upstream：https://github.com/SimVascular/svMultiPhysics，commit `c3f0bb892b765b718f61069ecd9726dbc6d177fd`。

`vendor/svMultiPhysics_stage_q/` 保存 WSL accepted Stage Q 的可编译源码及许可证。原 source snapshot 包含的 950 个未下载上游示例/安装包指针归档为 JSON provenance，不作为 raw LFS pointer blob 提交。全部 2431 个原始文件按 `reports/sv1_3q/source_patch.json` 的 after hashes 核对；仅发布包装中的两份 `.gitattributes` 用无 LFS 规则替换，原字节归档在 `sync_metadata/upstream_gitattributes.txt`。这不改变 solver source。

继承顺序：

1. Stage N：[PETSc lifecycle adapter](sv1_3n/svmp_petsc325_compat.patch) 与 [来源说明](sv1_3n/PATCH_PROVENANCE.md)，在 MPI_Finalize 前完成 PETSc cleanup。
2. Stage O：[final-output-on-stop](sv1_3o/final_output_on_stop.patch)，native stop 时保留 final VTU/checkpoint。
3. Stage P：[within-timestep ILU reuse](sv1_3p/reuse_within_timestep.patch)。
4. Stage Q：[跨步复用、adaptive recovery 与最终 stop timestamp](sv1_3q/ilu_rebuild_policy.patch)。单独的 `stop_ack_timestamp.patch` 是历史中间补丁，不要重复应用在最终完整补丁上。

从原 upstream 到最终 Q 的累计补丁：[upstream_to_frozen_stage_q.patch](upstream_to_frozen_stage_q.patch)。它包含新增 `sv13q_reuse.h`。恢复源码时选“vendor source + prepare_test_paths 恢复原快照占位指针”或“upstream + 累计补丁”，不要同时应用累计补丁和分阶段补丁。累计补丁已经在临时 pristine upstream 树中应用并逐文件比对，见 `sync_metadata/patch_reconstruction.json`；没有编译或运行 solver。

Stage M 的 PETSc 3.19 CUDA ghost backport 及 repair 历史也保留在 `sv1_3m/`，它属于历史路径，当前 PETSc 3.25.5 不需要把旧 PETSc patch 重新套上去。

实际 build provenance：`reports/sv1_3q/remote/svmp_reuse_build.json`。最终 executable SHA256 `0e509fe21424b5b1f731c4b4f519839b0104bfb2aab0817b01556cc4054f7fc7`，PETSc library SHA256 `b79e98d278e105aa8a1848bd08d66ef36a3e2425458610f97213d3fc8e45c5ef`。未提交 executable、CMake cache 或 toolkit。旧 wrapper/build 脚本保留供审计，移植限制见 `../PORTABILITY_NOTES.md`。
