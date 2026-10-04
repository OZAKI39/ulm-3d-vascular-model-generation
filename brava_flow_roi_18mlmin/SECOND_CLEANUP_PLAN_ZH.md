# SECOND ACTIVE-CONTEXT CLEANUP 删除前计划

本计划和保护清单在任何删除动作前生成。基线提交 `5622c1a9549c19306b28eb5479c87600d76997c0`；仅在独立分支 `cleanup/prune-retired-meshing-experiments` 工作，不改历史。

逐路径审查见 [SECOND_CLEANUP_CANDIDATES.csv](SECOND_CLEANUP_CANDIDATES.csv)，含host/path/type/size/tracked/referenced_by/classification/decision/reason；大于10 MB见 [LARGE_FILE_REVIEW.csv](LARGE_FILE_REVIEW.csv)。保护SHA见 [SECOND_CLEANUP_PROTECTED_PATHS.json](SECOND_CLEANUP_PROTECTED_PATHS.json)。未知项一律KEEP。

依据Python AST import及源文件中的CLI/路径字符串建立依赖闭包；再做git grep和只读smoke。生产scripts/vendor/particle/population、manufactured公式/probes/strict containment、mesh conversion与svMultiPhysics边界adapter保护。源码存在动态引用局限，不确定项保留。依赖图位于 `production_workflow/validation/second_cleanup/DEPENDENCY_GRAPH.json`。

- local: 拟删除2043个路径、302,140,463 bytes；保留所有唯一正式科学结果。
- server: 拟删除4254个路径、1,353,092,229 bytes；保留所有唯一正式科学结果。

删除MMG/Omega_h/refine专用构建、source/install、aborted workspaces、临时输入/metric/log/monitoring；MeshSim临时probe；已退休P1/P2/BC路线脚本、执行副本和中间报告；旧体积metric glue不再被未来Gmsh路径引用。必要日志先提取10–50行（日志不足10行则全长）及SHA到证据JSON，工具版本/关键数字汇入每条路线唯一summary。P0/P1/NW1原始baseline保持。

保留昂贵旧研究CFD/mesh/trajectory及审计数组：复用工具仍引用，且无法据现有证据排除未来科学引用。可能可再缩减的研究大文件标MANUAL_REVIEW_REQUIRED。生产SV/TetGen/PETSc/MPI、环境、系统缓存及兼容链接不删除。

tracked路径逐项git rm；未跟踪仅按冻结删除表unlink，不跟随符号链接；运行前复查size/mtime，变动或意外保护重叠即中止。不存在归档复制。服务器14项已有缺失另列，不计本轮释放字节。

仅执行module import、strict containment和快速核心unit tests、现有正式flow/mesh loader、svMultiPhysics adapter只读smoke。不得调用会运行零源population的旧cleanup_validation.py，不运行任何mesh/CFD/trajectory/population。完成报告、独立commit后停止。

最终删除前复核：去除重复的嵌套历史报告、raw solver logs；正式输入、昂贵研究场/mesh和逐接受状态不动。日志片段及SHA已先提取，未运行任何科学计算。最终逐路径清单数如下。

- local: 2090路径，353,545,022 bytes，590 tracked。
- server: 4303路径，1,404,637,843 bytes，579 tracked。
