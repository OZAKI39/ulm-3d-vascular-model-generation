# BraVa CFD 服务器迁移本地保全报告

完成时间：2026-10-09T10:12:37.278991+00:00。状态：**本地备份及完整性核验完成，等待用户给出新服务器指令。**

本轮只迁移 BraVa CFD。没有运行 CFD/粒子任务，没有向新服务器上传或安装环境，没有删除原服务器文件，也没有覆盖本地现有工程。MEVO 声学工程及其数据未纳入备份。

## 本地交付

备份根目录：`/home/lzy/projects/server_migration/brava_vast4090_20261009T094459Z`。

- [备份状态](BACKUP_STATUS.json)：机器可读的完成状态与范围。
- [恢复计划](RESTORE_PLAN_ZH.md)：收到新服务器信息后的路径、环境和验证步骤。
- [范围与排除项](inventory/MIGRATION_PLAN.json)：逐目录的选择依据。
- [源文件清单](inventory/SOURCE_MANIFEST.jsonl)：服务器端 SHA-256、大小、权限、UID/GID 和链接目标。
- [本地校验汇总](inventory/VERIFICATION_SUMMARY.json)：本地独立读取并计算 SHA-256 的结果。
- [逐项校验](inventory/LOCAL_VERIFICATION.jsonl)、[额外文件检查](inventory/PAYLOAD_EXTRAS.json) 和 [原冻结身份核对](inventory/FROZEN_IDENTITY_CHECK.json)。

共保存 **363,344 个普通文件**、**4,554 个符号链接**。普通文件逻辑大小 **45,381,897,009 bytes（45.38 GB / 42.27 GiB）**；这不是去重后的磁盘占用。所有选定普通文件均通过大小及 SHA-256 核验，权限和链接原文也已核对，失败数为 0。备份数据目录不存在清单外文件，只有路径容器目录。源清单内运行时特殊对象计数为 1。

| 内容 | 普通文件逻辑大小 |
| --- | ---: |
| BraVa | 29.918 GB |
| Native CFD stack | 8.765 GB |
| Python and metadata | 1.324 GB |
| CUDA | 5.376 GB |

主体包含正式 BraVa 网格、边界条件、冻结流场、计算输入和 checkpoint、CFD 审计/验证证据、代码与未提交修改、日志和报告；以及 native svMultiPhysics/PETSc/MPI/VTK/TetGen 源码、构建记录、二进制、启动包装器和 CUDA 13.2 工具链。BraVa 内部用于流场审计的历史粒子证据随项目保存，未改写其结果身份。

已排除无关工程、旧 CUDA 12.6 安装包和重复解压工具链、可重建缓存、运行时 socket/device/FIFO 及账户/SSH 凭据。没有复制整个通用 Torch 环境或 MEVO Python 环境；只记录了被 BraVa 使用的共享 Python 依赖版本。

两次早期传输为调整复制策略主动中断，最终四组传输均成功结束。一个 12,189,696-byte 的续传临时文件及其空目录已从本次新建备份目录清理，记录见 [清理记录](inventory/TRANSFER_SCRATCH_CLEANUP.json)。这不涉及源数据删除。已有本地文件仅用于校验后复用复制，备份没有以硬链接绑定现有工程。

初次完整核验发现一个本地 `.git/config` 比源文件多出 `branch.main.vscode-merge-base` 配置项；源文件仍与原清单一致。已重新复制该配置并独立核验大小、SHA-256 和权限；其余文件使用首次完整读取的通过结果。初次失败记录与 [单项修复复核记录](inventory/GIT_CONFIG_REPAIR.json) 均保留，最终核验失败数为 0。

## 服务器状态

最终读取时刻：2026-10-09T10:08:07.094865+00:00，详见 [最终状态](inventory/FINAL_SERVER_STATUS.json)；系统/编译器/软件包与动态库信息见 [初始快照](inventory/SERVER_SNAPSHOT.json) 和 `environment/`。

- SSH 可连接；GPU 为 RTX 4090，驱动 595.84，显存 24,564 MiB，观测使用 1 MiB、GPU 利用率 0%。
- BraVa 的四个系统 supervisor 任务均为 EXITED；未观察到运行中的 BraVa 求解器/MPI 计算。两个历史 supervisor 守护进程仍存在，不代表 CFD 计算正在运行。
- 容器 CPU 配额约 7.68 核，内存上限 42.06 GiB，系统 Ubuntu 24.04.4；原服务器磁盘剩余约 24.70 GB、使用率 93%。
- 备份读取后 cgroup 计费内存接近上限，其中文件缓存约 39.16 GiB、匿名内存约 0.69 GiB；`memory.events` 的 oom/oom_kill 均为 0。不能将这部分文件缓存解读为 CFD 工作内存，详见 [内存记录](inventory/FINAL_MEMORY_DETAILS.json)。
- 其它工程和基础服务保持原状；未停止或迁移它们。

## 代码及数据身份

| 快照 | Git HEAD / 状态 |
| --- | --- |
| 服务器 BraVa | `91ebc02162e7b0c774cf649b5669bcc5c8026e54`；status 输出 43 行，完整工作树与 Git 信息已保存 |
| 当前本地 BraVa | `a4b238b27437e1760b4f74fc219ede4f0180ee8f`；工作树干净，分支 `validation/gmsh-reference-v2` |
| 服务器 native svMultiPhysics source | `c3f0bb892b765b718f61069ecd9726dbc6d177fd`；本地补丁保留 |

备份前后服务器上述 Git HEAD/status 一致。服务器和本地是不同代码状态，不能相互冒称；后续恢复需要对齐。原网格、正式冻结流场、PETSc 库身份与既有记录一致，具体 SHA-256 见 [关键文件](inventory/CRITICAL_IDENTITIES.json)。备份完整性通过不等于局部流场物理精度通过，也不修改已有失败门槛或科学结论。

## 新服务器恢复边界

`payload/` 中保留原服务器绝对路径布局。环境版本、动态库解析、CMake/PETSc 构建信息已保存；这是一份 BraVa CFD 数据与工具链备份，不是整个服务器操作系统镜像。

当前 native 链为 PETSc 3.25.5、OpenMPI 4.1.6、CUDA 13.2，原 GPU 构建架构为 89。新服务器的 GPU、驱动和系统 ABI 尚未知，必要时用已保存源码重建。Python 环境有 3.12/3.13 及共享路径依赖，应按使用阶段分别恢复；不能直接混装所有 requirements。绝对符号链接保留原文，在本地备份目录内不保证可执行。

下一步等用户提供新服务器连接信息及执行指令后，进行硬件/环境盘点、暂存上传、哈希复核、环境恢复和小型运行验证。旧 supervisor 配置仅作参考，不自动启用；新计算另按用户任务执行。

## 复核入口

在备份根目录执行 `sha256sum --check METADATA_SHA256SUMS` 可检查本报告、环境记录、源清单和验证记录等元数据。运行 `python tools/verify_local.py` 可再次逐文件核验 payload；它仅操作本地校验记录，不连接服务器、不执行求解器。重新生成校验记录后，应重新生成元数据校验清单。
