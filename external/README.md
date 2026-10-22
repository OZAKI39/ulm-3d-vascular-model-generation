# 当前原生软件安装

安装根目录：`/home/lzy/projects/new_vessel_flow_solvers`。本目录保留源码与原服务器身份，当前二进制、运行环境、编译脚本和小测试证据位于安装根目录。

| 组件 | 当前安装 | 实际小测试 |
|---|---|---|
| svMultiPhysics | 原有效源码 commit c3f0bb8 + 原本地补丁，新编译 | CPU/GPU 各 55 步小管道通过 |
| PETSc | 3.25.5，CUDA 12.6，系统 OpenMPI 4.1.6 | GPU 矩阵和向量运算执行证据已保存 |
| SimVascular / TetGen | 2023.05.31，独立复制安装，嵌入式 Python 3.5 | 454 点、2053 四面体通过 |
| VMTK | 独立环境，1.5.0 / Python 3.11 / VTK 9.2.6 | 小圆管法向与封口通过 |
| Ultraliser | 原 executable 的独立副本，commit 3e4b0ee | 非对称分叉小样例通过；极小特殊输入崩溃记录保留 |
| 粒子 C++ 查询 | 工程自身编译 | 原有字节一致性测试 |

原求解器 SHA `0e509fe21424b5b1f731c4b4f519839b0104bfb2aab0817b01556cc4054f7fc7`、原 PETSc SHA `b79e98d278e105aa8a1848bd08d66ef36a3e2425458610f97213d3fc8e45c5ef` 仅描述历史服务器二进制。新构建不冒充它们；原环境记录已移到 `solver_support/original_runtime/`。

当前入口读取 `configs/environment.json`，校验安装清单自身的 SHA、solver/PETSc/VTK/启动脚本及测试凭据，且要求传入 build JSON 与已核验版本完全一致。物理参数、线性求解选项、稳态停止和 checkpoint 验收未变。安装未修改正式 CFD 或轨迹。

启动入口为安装目录 `bin/svmultiphysics`、`bin/solver-env`、`bin/simvascular`。网格模块 `meshing.launch_embedded()` 使用独立 SimVascular 启动器；不要把 Python 3.5 的 `sv` 扩展导入项目的 Python 3.13 环境。VMTK 同样通过自身环境运行。

运行依赖由安装根目录和系统库提供；系统 CUDA 位于 `/usr/local/cuda-12.6`，MPI、编译器和数学库由 Ubuntu 提供。没有借用旧项目环境。安装说明、已记录限制和重建命令见安装根目录 `README.md`。

LAMMPS 是可选研究桥接工具，不参与当前正式微泡或流场计算，本次未安装。
