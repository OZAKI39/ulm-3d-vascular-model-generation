# 本地软件与原生求解器

本轮使用本地 WSL，没有连接远端服务器。旧服务器迁移记录在 `server_migration_record/`，仅供理解过去的软件与二进制身份，不是当前运行环境。

## Python 环境

- 独立环境：`new_human_vessel_flow/.venv`，Python 3.12.3，由系统 Python 新建。
- 当前自有代码通过 `pip install -e .` 安装。普通 `import` 指向本工程 `src/`；无需旧项目的 `PYTHONPATH`。
- NumPy 2.2.6、SciPy 1.15.3、VTK 9.7.1、PyVista 0.49.0；完整实际版本以根目录 `requirements-lock.txt` 为准。
- GPU 微泡使用 PyTorch 2.6.0+cu126，编译原 CUDA 扩展使用系统 `/usr/local/cuda-12.6` 与 ninja。没有修改 GPU 数值内核。
- 当前设备是 RTX 4060 Laptop，8 GiB 显存、计算能力 8.9。测试编译最多使用 2 个 CPU 作业，普通数值线程限制为 1。
- 历史可选 LAMMPS 接口保留在依赖闭包中，但本轮没有安装或验证 LAMMPS；当前 BraVa 群体的 CPU/GPU 路径不依赖它。它不代表小鼠模型数据已经迁入。

重建环境请使用 `scripts/rebuild_python_environment.sh`；它创建新环境并按锁定版本安装，不复制旧虚拟环境。初次联网下载可能较慢。下载超时、重试和最终安装记录在 `validation/logs/`，未把安装中断当成科学测试失败，也没有隐瞒中断。

## 第三方原生依赖

统一安装目录是 `/home/lzy/projects/new_vessel_flow_solvers`：

| 组件 | 实际位置 | 本轮验证范围 |
|---|---|---|
| SimVascular 2023-05-31 / TetGen | `bin/simvascular` → `simvascular-2023.05.31/` | 启动其 Python 3.5.5，成功导入 sv、VTK 8.1.1 与 TetGen 类；没有生成新网格 |
| BraVa svMultiPhysics | `brava_svmultiphysics/bin/svmultiphysics` | 从新工程的生产源码成功编译、动态库检查、无输入启动检查；没有运行 CFD |
| PETSc | `petsc-3.25.5-cuda12.6/` | 明确记录库位置与 SHA256，由新 BraVa 求解器链接 |
| 原生 VTK | `vtk-9.3/` | 求解器编译/运行依赖，与 Python VTK 分开 |
| MPI | 系统 `/usr/bin/mpiexec`，OpenMPI 4.1.6 | 实际版本检查；本次未发起 MPI CFD |
| CUDA | 系统 `/usr/local/cuda-12.6`；WSL 驱动库 `/usr/lib/wsl/lib` | nvcc 实际版本检查；GPU 微泡验证单独记录 |
| 编译器 | 系统 GCC/G++/GFortran 13 | 生产源码构建日志保留 |

SimVascular 包含自己的 Python 和 VTK，因此通过明确的启动器设置它自身的库路径。这是第三方软件的必要依赖，不是加载旧 BraVa 或小鼠工程 Python。主项目 Python 环境不接收这些设置。

## 为什么新编译了一份 BraVa 求解器

检查发现，统一目录中先前已有的求解器源码缺少 BraVa 正式路径需要的 GPU 主机缓冲同步补丁，预条件器复用头文件名称也不同。直接替代会丢失生产代码身份。

因此保留原安装不动，将 BraVa 自己的生产源码编译到独立的 `brava_svmultiphysics` 子目录。新工程没有引用小鼠 Python 或它的算例。共用的是明确安装的第三方 PETSc、VTK、MPI 和 CUDA 库。

- 上游源码基线：`c3f0bb892b765b718f61069ecd9726dbc6d177fd`。
- 当前源码：`native/svmultiphysics/`。
- 生产补丁与复用支持：`native/production_patch/`。
- 构建脚本：`scripts/build_native.sh`。
- 源码、构建脚本、新二进制和 PETSc 身份：[BUILD_IDENTITY.json](../native/BUILD_IDENTITY.json)。
- 当前配置：[native.json](../configs/cfd/native.json)。
- 实际检查：[NATIVE_ENVIRONMENT.json](../validation/NATIVE_ENVIRONMENT.json) 和 [MESHER_ENVIRONMENT.json](../validation/MESHER_ENVIRONMENT.json)。

P2/P1 实验补丁只保留在诊断目录，没有合并到这份生产源码。新编译二进制具有自己的 SHA，不能冒充过去服务器上的二进制。新编译器、库版本和新环境的加载测试也不能证明完整 CFD 数值等效；本次没有执行新 CFD，这一项为 `NOT_RUN`。

## 软件安装通过不等于科学问题解决

无输入运行检查的正常结果是退出码 1，并提示需要 XML 文件；代码在开始求解前退出。第一次检查器误以为应出现 `Usage`，该记录仍保留为 `NATIVE_ENVIRONMENT_INITIAL_PROBE.json`；修正了检查器对真实启动提示的识别，再次执行通过。没有调整物理验收容差。

同样，CPU/GPU 小测试一致，只证明本次测试范围内实现一致；全局守恒、稳态或动画流畅，均不能替代局部场精度认证。
