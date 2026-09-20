# 可移植性与历史证据

冻结读取入口、manifest 校验和新测试只使用当前 checkout 相对路径。它们不需要 `/home/lzy`、原 FEM Git 历史、远端 SSH 或 GPU。Python 读取包见 `pyproject.toml` 与 [版本记录](sync_metadata/python_reader_versions.json)。

历史 `scripts/`、`configs/`、`reports/` 和 `tests/` 保留来源字节，不据此承诺原开发机自动化命令可在任意机器直接重跑。尤其 `remote.py`/`invoke.py` 依赖旧 SSH transport 和旧本地 Git blob，build/runner 脚本引用旧 `/workspace/` 安装前缀，部分旧报告依赖 `/home/lzy/` 及未打包历史输出。完整扫描见 [portability_inventory.json](sync_metadata/portability_inventory.json)，每项列路径和行号。不要为运行本阶段校验去执行这些历史命令。

原 Stage Q 20 项测试保留不改，可在当前分支运行：

```bash
python -B scripts/fem_freeze_sync/prepare_test_paths.py
python -B -m pytest -q tests/test_sv13q_*.py
```

准备脚本只在被 Git 忽略的 `external/sv13q/` 生成可删的源码副本，恢复原测试期待的路径及上游 `.gitattributes` 原字节，不执行源码或下载依赖。已提交 vendor 的 `.gitattributes` 仅关闭上游 LFS filter；原文件位置列于 `sync_metadata/upstream_attribute_map.json`（根规则保存在 `sync_metadata/upstream_gitattributes.txt`）。上游 950 个历史示例/安装包本来就是未下载的 LFS pointer；其原始指针文本、OID、逻辑大小和 SHA 全部归档在 `sync_metadata/upstream_lfs_pointer_inventory.json`，不把 raw pointer blob 交给 GitHub 触发未知 LFS 对象检查，也不下载这些示例。准备脚本可在 ignored 源码副本里恢复原指针字节；正式冻结 flow/mesh/boundary/checkpoint 均为实际完整文件，不是 LFS pointer。此分支不启用 LFS。

Stage Q 当时的完整历史 suite：505 passed、8 failed、65 skipped；原 Stage Q 本阶段测试 20 passed、0 failed。下面是仍保留的早期失败证据，不被改写为 PASS，也未转换成 pytest xfail：

- `tests.test_sv11_short_run_mass::test_actual_step_ten_mass_strictly_improves`
- `tests.test_sv13h_petsc_build::test_actual_petsc_cuda12_build_acceptance`
- `tests.test_sv13j_svmp_gpu_smoke::test_actual_official_fluid_gpu_smoke_acceptance`
- `tests.test_sv13l_official_gpu_smoke::test_actual_official_gpu_smoke_acceptance`
- `tests.test_sv13m_official_gpu_smoke::test_official_complete_flow`
- `tests.test_sv1_mass_balance::test_actual_field_mass_conservation`
- `tests.test_sv1_solver_result::test_real_solver_converged`
- `tests.test_sv1_solver_result::test_required_steady_intervals_reached`

本次同步只运行新增完整性测试及可闭合的 Stage Q 测试；不声称精简分支的整个历史 pytest 都能运行。原始 suite 的很多测试需要被排除的依赖树、旧输出或原开发 Git 历史，失败原因不能混同于新的科学失败。[旧完整测试记录](reports/sv1_3q/pytest_summary.json) 原样保留。

`source_manifest.csv` 是本次同步文件清单，不等同于旧阶段 `delivery_manifest.json`：后者完整记录了当时开发机交付范围，包含本次明确不提交的 build/dependencies/intermediates。遗漏清单 [omitted_artifacts.json](sync_metadata/omitted_artifacts.json) 提供原 WSL 路径、size、SHA256；无必要 Particle-0 文件依赖这些遗漏。

Frozen solver XML 保留原相对路径，在 `frozen_reference/run/solver.xml` 下可完整解析到相邻 `SV_MESH/`。`frozen_reference/run/STOP_SIM` 与 `sv13q_stop_ack.json` 是完成时的停止证据；此目录用于只读交接，不是新启动目录，不要拿留下的 STOP_SIM=0 作为初始运行输入。Native checkpoint 保留最终 step71，随仓库提交，不需 pointer。历史运行时 binary/library 的 SHA 保留为 provenance，不上传二进制安装树。
