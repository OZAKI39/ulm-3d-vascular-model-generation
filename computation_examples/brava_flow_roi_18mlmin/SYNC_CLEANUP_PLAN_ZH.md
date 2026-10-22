# BraVa 当前工作树清理与 GitHub 同步计划

本轮仅 filesystem / Git / documentation / 快速 import 与 unit tests，禁止科学重算。

目标分支：`sync/roi-flow-mb1500-wss-stop-analysis-20260928`。远端目标是既有同步 checkout HEAD 的祖先，ahead 11 / behind 0；保留这些已有提交并 fast-forward，不改写历史，不使用 LFS，不更改 .gitattributes。当前算例仓库没有 remote，不能把其独立历史直接推到 GitHub。

删除前全量清单见 SYNC_FILE_INVENTORY.json（两棵工作树，不跟随 symlink，排除 .git 数据库），大文件见 LARGE_FILE_REVIEW.csv。UNKNOWN 一律保留。共享源码仅精确同步 src/tests，禁止 wholesale stage 共享工作树。

## 已核对的候选与保护边界

- 删除已被 V2 取代的 V1 原始 .msh/展示副本、旧 cap/V1 绘图和 CPU monitor；保留 V2 生成/QA 仍引用的原 BREP、INLET、BOUNDARY_TOPOLOGY 和 RAW_MESH。
- 删除 GPU step10–70 的 VTU/checkpoint 及重复 last checkpoint；保留 state71、原 execution.json 中全部区间与318次求解记录。
- 删除旧 CPU 瞬态 VTU、重复 last checkpoint 和原日志，先保存 timeout summary/日志hash。CPU step10 checkpoint 仍由未改动的 GPU prepare 身份检查依赖，因此保留为历史证据，绝不作为 current field。
- 保留 2.82MB GPU solver.log：不属于巨型日志，且原 certifier 直接读取；另提取精简 GPU evidence。保留195KB/5秒 GPU samples，因冻结 identity SHA 与 certifier 直接依赖，删除会使身份链断裂。删除每秒 supervisor samples，保存原hash和汇总。
- 当前图只保留三个最终 PNG，移除重复 PDF/SVG 与退休 mesh/runtime 图。相关展示代码/导航同步改为 PNG，不改数值。
- GitHub 旧完整 svMultiPhysics 第三方 checkout 删除，保留 patch、build recipe、SHA、wrapper、RA header、policy 和 solver_support。
- 删除 pytest 临时目录和 Python cache；不删除兼容链接、正式 flow、accepted states/seeds、共享测试。已在上次清理移除的 MMG/Omega/refine/MeshSim/P2/BC-B 实现不恢复。
- >=100MB 一律不上传；保留本地且以 SHA/size/server/config/recovery 记录。历史原始 population states 和大批逐轨迹派生文件不整体复制进 Git；保留账本、summary、身份、最终图及准确取回路径。
- 服务器仅只读文件/身份与必要 import smoke；不清理服务器、不启动任务。

## 科学状态冻结

REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC + MIXED_FIELD_ERROR + P1_FIELD_TRAJECTORY_IMPACT_LARGE（DIAGNOSTIC_ONLY）。V3、正式 CFD-vs-oracle、dt convergence、新 BraVa validation 均 NOT_RUN。NO PRODUCTION PHYSICS CHANGE。

## 执行前审计

候选删除 1633 项，94,731,213 bytes。所有候选已记录逐文件原因；删除前再次核对 sha256。

最终引用复查修正：历史删除清单中的3个报告文件名已被当前权威报告复用；已恢复本轮开始时的准确SHA，不删除、不采用旧内容。净删除1631项、94,766,770 bytes；逐秒记录另缩减2,414,984 bytes。逐项执行与恢复记录见 production_workflow/validation/sync_cleanup/CLEANUP_ACTIONS.json。
