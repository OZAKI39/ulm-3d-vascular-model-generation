# Stage P 证据与复查入口

本阶段运行目录已经包含原始证据。运行脚本会拒绝覆盖既有案例，不应为复查而重新启动 CFD。

- `configs/sv1_3p/runplans/`：每次实际血流运行的固定计划。最终生效选项、进程命令、库 SHA256、检查点 SHA256 和单调时钟时间在对应 `reports/sv1_3p/remote/*_execution.json`。
- `logs/sv1_3p/remote/`：完整求解器日志、编译日志和实际 PETSc 参数帮助输出。失败候选的日志与负收敛原因保留。
- `outputs/sv1_3p/`：原生 VTU、检查点、XML 和 PETSc 事件记录。验收脚本核对文件 SHA256，并由单独进程重读最终场。
- `patches/sv1_3p/reuse_within_timestep.patch`：P1 的单文件补丁。`ksp_lifecycle_audit.json` 是修改前审计，`source_patch.json` 保存前后逐文件哈希。
- `external/petsc325_cuda13_hypre/`、`external/petsc325_cuda13_amgx/`：隔离构建记录和原始依赖包归档；远端实际前缀、编译参数、版本和库路径均写入构建报告。
- `analysis_history/` 与 `P1_SMOKE_acceptance_before_parser_fix.json`：保留早期分析结果。更新仅修复证据解析和补充失败运行的参数审计，没有重跑这些 CFD。
- `native_artifact_mirror.json`：远端新增二进制与依赖归档的本机哈希镜像。`delivery_manifest.json` 是最终交付文件清单。

实际血流执行顺序为 P1 smoke、P1 10 步窗口、新 hypre 二进制上的原 ASM2+ILU2 短检查、P2 smoke、P3 smoke、P4 smoke。AmgX 官方构建失败，没有运行 CFD。只有 P1 达到 10% 门槛，故最后只启动一次 `REAL_VASCULAR_GPU_PC_WINNER`。

独立 GPU 小算例只用于证明新增 hypre 库可运行及读取真实帮助参数，不能当作血管性能结果。历史 Stage O 的 10 步基线没有重跑。

只读验收命令：

```bash
.venv/bin/python -B -m pytest -q -p no:cacheprovider tests/test_sv13p_*.py
```

完整测试保留历史失败，不改旧测试，也不因为 pytest 启动任何血流求解。完整输出和逐项统计见 `logs/sv1_3p/pytest_full.log`、`pytest_summary.json`。

比较限制：各候选只运行一次。P 的正式窗口记录 GPU 事件时间，原 A 为历史非诊断运行；其同步开销包含在候选实测时间中。事件有嵌套，不相加；设备级显存采样不是精确进程峰值。CPU/GPU 科学等价验证仍推迟。
