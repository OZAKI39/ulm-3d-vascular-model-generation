# BraVa CFD 迁移恢复计划

本目录只用于 BraVa CFD 的离线保全。新服务器上传、环境安装和任何 CFD 计算都等待用户下一条指令；没有自动恢复队列。

## 数据布局

`payload/` 后面的相对路径就是原服务器绝对路径去掉开头的 `/`。主入口：

- `payload/workspace/brava_flow_roi_18mlmin_20260928/`：正式网格/边界/流场、计算输入、checkpoint、原始日志、生产 GPU 修复、Gmsh V2 reference、P2/P1/common-geometry/边界及线性代数诊断。BraVa 内用于 CFD 审计的历史粒子证据随目录保留，不改变其原始结果身份。
- `payload/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/`：共享 svMultiPhysics、PETSc、MPI、VTK、TetGen 源码、构建文件、已有二进制与启动包装器。旧 CUDA 12.6 安装器和重复解压目录等排除项见 `inventory/MIGRATION_PLAN.json`。
- `payload/root/particle8_2_runs/env/` 和 `payload/.uv/python_install/`：旧 CFD 后处理 Python 环境及其基础解释器副本。
- `payload/usr/local/cuda-13.2/`：当前 native CFD 构建/运行使用的 CUDA 工具链；不是宿主 NVIDIA 驱动。
- `payload/etc/supervisor/conf.d/brava*.conf`：历史任务配置，仅供参考。恢复文件不等于授权启动这些任务。
- `environment/`：Python 精确包版本、解释器信息、系统软件包清单、动态库解析记录。共享 MEVO Python 仅记录被 BraVa 使用的环境信息，没有备份 MEVO 项目代码或声学数据。

`inventory/SOURCE_MANIFEST.jsonl` 保存服务器端逐文件 SHA-256、大小、权限、原 UID/GID 以及符号链接原文。`inventory/LOCAL_VERIFICATION.jsonl` 和 `VERIFICATION_SUMMARY.json` 保存本地独立复核结果。Unix socket 等运行时对象不会迁移；旧 PID 文件即使作为证据保留，也不得用于新服务器的进程控制。

## 收到新服务器信息后执行

1. 读取新机器自己的 AGENTS/实例指南，检测实际 OS、CPU/cgroup、RAM、磁盘、GPU 架构、驱动、CUDA 和现有环境。不要照搬旧机器 worker 数。
2. 先核对本地清单，再上传到新服务器的独立暂存目录并复核 SHA-256。优先保持上述 `/workspace/...` 绝对路径；路径冲突时先制定映射，不批量改写冻结输入和历史 identity。
3. 对齐本地当前 BraVa 代码与服务器快照。两者的 Git HEAD/status/diff 已单独保存；备份过程没有用服务器旧代码覆盖本地工程。不要直接把任一快照标为另一版本。
4. 按新机器条件恢复 native MPI/PETSc/VTK/svMultiPhysics。原生产链使用 PETSc 3.25.5、独立 OpenMPI 4.1.6、CUDA 13.2；原 PETSc 构建指定 GPU 架构 89。新 GPU 或 ABI 不兼容时，使用已保存的源码/configure/CMake 信息重建，并记录新的二进制身份。不要安装或覆盖宿主 NVIDIA 驱动，也不要将 CUDA stubs 当运行库。
5. 建立独立 CFD Python 环境。优先复现相应验证阶段使用的解释器和版本，而不是合并所有 requirements。旧后处理环境为 Python 3.13.15；最新 BraVa P2 诊断也使用过共享路径下 Python 3.12.3 / NumPy 2.2.6 / SciPy 1.15.3。包清单是观测记录，不能当成一个跨环境混装的 lock 文件。
6. BraVa population 环境通过 `parent_torch.pth` 引用 `/venv/main` 的 Torch；本轮没有复制整个通用 Torch 主环境。该依赖已记录。若后续任务需要粒子/GPU 后处理，再按记录建立相容环境或受控的路径映射；仅 native CFD 不要求启动任何粒子批次。
7. 先做环境验证：动态库完整性、MPI 启动、Python 导入、冻结网格/场读取、已有小型数值/CPU-GPU一致性检查。数据完整性通过不代表新主机数值验证通过。不得以生产 CFD 重跑作为默认安装测试。
8. 只有用户明确给出计算任务后，才安装/启用对应 supervisor 配置和新运行目录。所有旧失败门槛、停止状态、原正式场身份保持，不恢复已停止的旧全域加密或旧微泡队列。

## 范围边界

此备份不覆盖 MEVO 声学微泡工程、HemoCell、LAMMPS 或其它旧小鼠算例，也不包含账户凭据、SSH 私钥或宿主驱动。原服务器和本地现有工程均保留，未删除、迁址或覆盖。

恢复时继续区分：正式 BraVa 历史校准场；未来中心线阻力估算出口压力方案；旧 1565 条粒子结果（141 截止、2 数值失败）；尚未产生正式大批量结果的新固定灌注协议。备份不修改任何科学结论。
