# BraVa CFD 新服务器迁移与环境验证报告

完成时间：2026-10-09 18:20:02 CEST（Europe/Zurich）。状态：**MIGRATED_AND_ENVIRONMENT_VALIDATED**。

仅迁移 BraVa CFD。数据与源码已落入独立目录，Python 和 native 依赖有各自前缀，文件核验及 CPU/MPI/CUDA 小型运行验证通过。本轮安装/验证作业均已退出，专用 supervisor 已关闭，没有恢复旧生产队列。

## 使用入口

| 项目 | 新服务器位置 |
| --- | --- |
| SSH | `root@95.3.33.46:15283` |
| 独立根目录 | `/workspace/brava_cfd` |
| BraVa 代码与关键数据 | `/workspace/brava_cfd/project` |
| 独立 Python 环境 | `/workspace/brava_cfd/envs/cfd` |
| 激活脚本 | `/workspace/brava_cfd/activate.sh` |
| 求解器入口 | `/workspace/brava_cfd/bin/brava-solver` |
| MPI 入口 | `/workspace/brava_cfd/bin/brava-mpi` |
| 复现配置与脚本 | `/workspace/brava_cfd/runtime` |
| 迁移清单与来源 | `/workspace/brava_cfd/provenance/migration` |

从本机连接可使用本目录的 `connect.sh`；它使用已有 SSH 密钥，没有将私钥复制到服务器。登录后：

```bash
source /workspace/brava_cfd/activate.sh
cd /workspace/brava_cfd/project
python -c "import sys; print(sys.executable)"
brava-mpi -n 2 /bin/hostname
```

后续计算使用新入口和独立的新算例输出目录。`project/gpu_solver_fix/svmultiphysics` 是保留身份的旧 CUDA 13 可执行文件；当前入口指向新构建，未用新文件覆盖旧 binary。历史报告、命令和审计脚本中的旧绝对路径/哈希仍是原执行证据，不应直接冒充新服务器运行身份。

## 迁移范围与完整性

本轮筛选 **20,561 个普通文件、104 个符号链接，普通文件逻辑大小 2,415,693,628 bytes（2.42 GB）**。加上新建环境和构建产物，最终根目录逻辑用量约 4.08 GB；未将本地 45.38 GB 完整备份整包上传。

保留正式网格/端口/边界条件/冻结流场和最终保存状态，当前 CFD 支持代码，固定 PLC 中细档对照场，Gmsh V2 冻结参考，最近八项 CFD/P2/几何/线性代数审查及相关原始证据，以及实际 production svMultiPhysics/PETSc 源码和 MPI/VTK/TetGen 依赖。

未迁移旧微泡群体原始批次、动画、废弃加密批次的大量中间帧、旧环境与重复 CUDA/PETSc 工具链、完整 Git 历史和其它工程。完整历史仍保留在本地备份与本地 Git；本次远端 `project` 是带文件身份清单的 CFD 工作集，并非完整 Git 历史镜像。SimVascular 旧 GUI/嵌入式 Python 环境未迁入，本轮验证范围为已有网格上的 native CFD 与字段处理。

本地当前 Git `a4b238b27437e1760b4f74fc219ede4f0180ee8f` 的相关追踪文件优先，服务器快照补充已执行源码和原始数据。两者有 18 个选定文件差异，逐项保留选择依据，不相互冒称版本。见 [范围](provenance/migration/SELECTION_SUMMARY.json)、[来源身份](provenance/migration/SOURCE_IDENTITY.json)、[逐文件清单](provenance/migration/SELECTED_MANIFEST.jsonl)、[差异记录](provenance/migration/LOCAL_SERVER_DIFFERENCES.json)。

构建前 17,792 个 bootstrap 文件全部逐项核验通过。最终完整检查中 **20,560 个普通文件与选定来源逐字节一致**；另 1 个 `petscvariables` 由 PETSc configure 按新源路径和 `arch-brava-cuda126` 重生成。旧文件、预期新文本及前后 SHA 都保留并验证。16 个项目绝对链接改为目录内相对链接，全部目标可解析。科学网格、场、边界、原记录数组与求解器科学源码没有改变。

[构建前校验](reports/BOOTSTRAP_INTEGRITY.json) · [最终校验](reports/PAYLOAD_INTEGRITY.json) · [链接映射](provenance/migration/RELOCATED_LINKS.json)。

## 独立环境与构建

新机实测 RTX 4090 / SM89，驱动 550.120，CPU 配额 15.36 核、内存上限约 42.07 GiB，Ubuntu 24.04.4。原 CUDA 13 链不能直接用于该驱动；按原源码和数值配置重新构建 PETSc 3.25.5 与 production svMultiPhysics，使用新机已有 CUDA 12.6、GCC 13.3、同一 GPU 架构。没有升级宿主驱动或覆盖系统 CUDA。兼容性依据及 PTX 限制见 [NVIDIA 官方文档](https://docs.nvidia.com/deploy/cuda-compatibility/minor-version-compatibility.html)；最终接受以本机实际运行检查为准。

| 组件 | 配置 |
| --- | --- |
| Python | 3.12.3，`include-system-site-packages=false` |
| NumPy / SciPy | 2.2.6 / 1.15.3 |
| VTK Python / PyVista | 9.7.1 / 0.49.0 |
| CuPy | 14.2.0，CUDA 12 wheel；实测其 runtime API 为 12.9 |
| Native VTK / MPI | 原独立 VTK 9.3 / OpenMPI 4.1.6 |
| Native PETSc | 3.25.5，CUDA 12.6 / SM89，double / 32-bit index |
| Native svMultiPhysics | 原 production host-sync 与 RA 补丁源码重新构建 |

PETSc/native 安装、MPI、VTK、Python 环境均位于 `/workspace/brava_cfd`；缺少的 MPI 系统库用 `.deb` 解压到专用 `deps/sysroot`，没有安装进系统。未向 `/venv/main` 安装包。43 个 Python 包通过依赖检查，精确版本见 [lock 清单](reports/requirements.lock.txt)。构建命令见 [构建记录](reports/NATIVE_BUILD_STEPS.json)。

新 svMultiPhysics SHA-256：`b22f0aa6c757326966b4b41f5c605f0408ee609fe398e8041ebf6947a3b3bc8e`。新 PETSc SHA-256：`66d791cc2f2c074769b7e69e05d3c58cc2e8e741b0ac23fbeecb785c402f40b1`。旧二进制身份保留在原报告中，两套身份明确区分。

## 实际验证结果

- Python 独立性、模块导入、CuPy 实际 GPU 运算、现有 P1 sampler 的单四面体查询均通过。
- 原正式网格 135,852 点 / 692,313 tet、原冻结场与 Gmsh V2 的 36,211 点 / 174,785 tet 及原第 71 步结果读取通过，冻结文件 SHA 未变。
- PETSc 原生 ex19 的 CPU 单进程、CPU 双进程、CUDA 检查均通过。
- 使用已保存的官方小管道夹具，3,107 点 / 11,878 tet，分别执行 CPU1、GPU1、GPU2 各一个时间步。XML、网格和物理参数保持原夹具设置；两后端使用相同 FGMRES/ASM/LU/RA 科学选项，只调整 Mat/Vec backend 和 GPU 计时开关。日志实际显示 `seqaijcusparse/seqcuda` 及 `mpiaijcusparse/mpicuda`，不是仅检查输入参数。

| 比较 | 速度相对 L2 | 压力相对 L2 | 预设 1e-8 门槛 |
| --- | ---: | ---: | --- |
| cpu1 vs gpu1 | 3.3750e-15 | 7.3767e-17 | PASS |
| cpu1 vs gpu2 | 3.4617e-15 | 7.9679e-17 | PASS |

三组输出对原网格的节点坐标匹配误差均为 0。首次收集器错误地要求 native VTK 必须包含 `GlobalNodeID`，CPU 求解本身已成功；修正为对原网格做严格一一坐标对应后复用其已完成结果，没有重跑 CPU 或改数组。初次失败记录及修正后的验证脚本都保留。

[完整验证](reports/ENVIRONMENT_VALIDATION.json) · [初次收集器记录](reports/ENVIRONMENT_VALIDATION_INITIAL_COLLECTOR.json) · [原生运行身份](validation/NATIVE_SMOKE_IDENTITY.json) · [最终状态](reports/FINAL_SERVER_STATE.json)。

这是运行环境和小算例后端一致性验证，不能证明 BraVa 局部流场物理精度。完整 BraVa 新 CFD、P2 CFD、微泡轨迹和群体新批次均为 0。原 P2 资源/输出门槛、Gmsh V2 科学结论、历史正式场和旧 1565 条轨迹身份均保持；后续新 BraVa 边界仍按用户要求使用中心线阻力估算出口压力。

## 收尾与后续

35 个数据传输分段全部成功，临时续传文件和安装下载缓存已清理。本轮 supervisor 的安装、构建、验证作业均退出，专用 daemon 已关闭；配置均为不自动启动，不影响其它项目服务。新机 `/workspace` 未挂载持久卷，本地完整备份和本轮验证交付继续保留。

下一轮开发从本报告的路径、`project/ACTIVE_DEVELOPMENT_SCOPE.md` 与 `MIGRATION_STATUS.json` 进入。没有因换机自动恢复旧全域加密或 P2 大规模任务。旧环境重建脚本是执行证据，不要未经检查重新运行安装脚本覆盖现有环境。
