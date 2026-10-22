# 这些脚本做什么

- `run_small_tests.sh`：运行原有小型测试和本次路径保护测试，输出到全新的目录；不启动 BraVa 正式计算。
- `check_independence.py`：用本工程安装的包导入源码，并在禁止读取旧工作区的条件下读取既有数据。
- `check_real_gpu_field.py`：只在旧 BraVa 场的 2146 个位置比较 CPU/GPU 求值，不积分轨迹。
- `check_mesher_import.py`：由 SimVascular 自带 Python 执行，只加载工具，不生成网格。
- `build_native.sh`：从本工程保留的生产源码编译 BraVa 求解器。依赖目录在 `configs/cfd/native.json` 中明确列出。
- `seal_native.py`：记录本地构建身份，并检查库加载和无输入启动行为。不会进行 CFD。
- `record_migration.py`：一次性迁移审计；按复制前清单只读核对旧文件，更新来源对应表。日常运行不需要它，也不需要原工作区。
- `prefetch_cuda_dependencies.py`：可选安装辅助，按 Torch 的固定版本依赖下载 CUDA wheel，并核对 PyPI SHA256；缓存不属于科研数据。
- `rebuild_python_environment.sh`：在新目录重建 Python 环境，拒绝覆盖已有环境。
- `summarize_validation.py`、`check_current_links.py`：汇总已执行记录、检查当前报告链接；不启动计算。

本次实际执行过的目录变换脚本保存在 `docs/provenance/migration_tools/*.py.txt`，仅是审计文本，不参与运行，也不要作为生产入口重放。普通使用从 `human-flow` 或 `src/human_vessel_flow/main.py` 开始。
