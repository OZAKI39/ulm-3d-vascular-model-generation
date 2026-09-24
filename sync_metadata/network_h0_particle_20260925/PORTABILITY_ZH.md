# 路径迁移与复现边界

本分支将多个工作树合并为一个仓库。科学源文件、冻结数据和原报告保持原字节，历史 JSON 中的绝对路径仍作为原始来源证据保留。逐文件映射见 `local_file_inventory.csv`、`original_path_map.json`；服务器映射见 `server_file_inventory.json`。

## 当前验证

从仓库根目录运行 `python3 sync_metadata/network_h0_particle_20260925/verify_snapshot.py --scientific-inputs`。此检查仅用标准库，核对全部发布文件 SHA，以及 Network 受保护输入、Particle 科学源文件、服务器源码快照和 OLD/NEW 冻结场的原始哈希。完整清单不包含它自身，也不禁止克隆后生成额外文件。

使用 [已验证依赖](requirements-verified.txt) 的环境运行 `run_checks.py --output /tmp/ulm-workflow-checks`。该工具创建临时副本，只迁移以下 5 份 JSON：Network v1 输入哈希与来源、Network v2 算例差异、Particle OLD/NEW 合同与科学保护清单。原清单中的 92 个 `.pyc` 缓存条目仅在临时验证清单中排除；全部 114 个实际科学源文件仍使用原期望哈希检查。源码、数组、测试断言和容差不修改。迁移记录见 `relocation_receipt.json`，最终测试记录见 `summary.json`。

Particle 的默认 `frozen_reference/` 是历史 2 mm/s 场。原 Flow 工作树的早期 Stage Q 输入另存为 `upstream_stage_q_reference/`；工具在临时目录运行流场回归时选用它。不要直接把这两个目录互相覆盖。

## 实际重新计算

当前 Particle 复现命令见 [Particle README](../../particle_3d/README.md)，指定包含固定科学源码和输入的 `server_bundle`，以及新的输出名。该命令会真的积分，不属于本次同步检查。

历史 FEM/Network 提交、构建与收集脚本可能包含原 SSH 别名、`/workspace`、`/home/lzy/projects`、Python 和求解器路径。换机器时先按当前配置和原始路径映射设置工作目录、MPI/PETSc/CUDA 与 svMultiPhysics 安装位置，再执行相应阶段。`export_tools/` 记录本次导出的操作来源，带有原机路径，不是面向新机器的一键运行器。

环境目录和求解器二进制未复制。svMultiPhysics 固定源码在 `formal_3D_flow_solver/FEM_SimVascular/vendor/`；服务器的实际修改文件、补丁和 Git 版本在 `server_evidence/pinned_solver_source/`。Ultraliser 版本见 `ultraliser_dependency.json`，其本地差异只有行尾空白，第三方构建树不入库。完整原始脑数据集未复制；当前 H0 输入哈希闭包已单独保留。

停止的 Taylor-Hood 服务器 `official_smoke` 原有 6 份未下载的上游 LFS 指针。GitHub 因目标对象不存在拒绝这些占位文件（GH008），因此发布时将原始路径、哈希和指针内容完整保存为 `historical_lfs_pointer_evidence.json`，不保留冒充实际网格的占位文件。需要重跑此历史示例时应在对应上游仓库拉取实际 LFS 数据。当前 H0 与 Particle 输入使用完整普通 Git 文件。

原可视化代码和样式保持不变。已有动画、图片及其必要输入按逐文件清单保留；省略渲染逐帧缓存和不参与现有渲染的重复候选池。FFmpeg 来自已记录的 `imageio_ffmpeg` 安装包，环境未必提供 PATH 中的 `ffmpeg`。

115 项通过表示本次选择的当前回归通过，不等同于全部历史测试、重新求解 CFD 或重新积分轨迹；历史环境相关失败见清理记录中的前后比较。
