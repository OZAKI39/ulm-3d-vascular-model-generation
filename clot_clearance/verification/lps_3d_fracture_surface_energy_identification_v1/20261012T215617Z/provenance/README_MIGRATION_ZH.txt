2026-10-10 更新：本地完整迁移备份及 stage 暂存副本已按用户要求删除；本服务器完整项目保留。当前精选 GitHub 分支为 clot-clearance/y-clot-migration-sync-20261010。删除记录见 migration/LOCAL_BACKUP_DELETION_20261010.json。以下为迁移当日历史报告。

血栓项目迁移与环境验证报告 — 2026-10-09
状态：完成，PASS。

新服务器：root@95.3.33.46，SSH 端口 15283。
独立目录：/workspace/y_clot
代码：/workspace/y_clot/project
独立 Python 环境：/workspace/y_clot/envs/clot
原生依赖：/workspace/y_clot/deps
本次验证输出：/workspace/y_clot/validation

登录（从原本地电脑）：
ssh -p 15283 -i /home/lzy/.ssh/vast_step3b_ed25519 -o IdentitiesOnly=yes root@95.3.33.46
进入开发环境（在服务器执行）：
source /workspace/y_clot/env.sh
cd /workspace/y_clot/project

迁移范围
主载荷 48,307 个文件/链接，5,865,737,992 字节；另补充 Y 网格工作流必需的软件测试状态与证据 39 个文件，3,228,781 字节。合计 5,868,966,773 字节，约 5.87 GB。新编译产物及环境另占空间。
保留当前点云与早期流场耦合实验、前序连续交接与 C1 基线、原始轨迹/断键/逐节点受力审计、原生 Y 几何/网格/稳态流场、最新 4K 视频及绘图代码、参考样式/字体、工程约定、修改后的 LAMMPS 源码、svMultiPhysics/PETSc 源码和已验证的主机同步修补、SimVascular/TetGen 原生命令行运行库、必要 MPI/VTK 依赖。
历史调参运行全集、旧编译树、旧 CUDA13 二进制、重复归档、无关项目和完整 Git 历史未迁移；它们继续保存在原本地完整备份中。原代码工作树和完整备份未改动。
源代码建立了三个新的本地 Git 开发快照；这是迁移快照，不是恢复原仓库完整历史。原 HEAD/状态及源文件 SHA 记录在 migration/original 和迁移清单中。

环境
Ubuntu 24.04，RTX 4090，驱动 550.120。
Python 3.13.11；主要包沿用原环境锁定版本：NumPy 2.5.3、SciPy 1.18.1、VTK 9.6.2、PyVista 0.49.0、mpi4py 4.1.2。
项目私有 OpenMPI 4.1.6。LAMMPS 以 MPI+BPM+MOLECULE+共享库方式重新编译。
PETSc 3.25.5 与 svMultiPhysics 在新机重新编译，CUDA12.6、计算架构89，保留原有 GPU 主机同步修补。
沿用 -use_gpu_aware_mpi 0：GPU 矩阵计算，MPI 经主机内存传输；不是 CPU 后端替代。
未升级系统显卡驱动，未修改其他项目的环境。env.sh 只对当前 shell 生效。
完整版本清单：migration/requirements_new_server.lock.txt、validation/ENVIRONMENT.json。

验证结果
1. 48,307 个迁移项逐项 SHA-256/链接目标核对完成，908 个实验数据/结果文件完全一致；39 个补充依赖文件 SHA-256 核对通过。
   初次核对发现 PETSc configure 自动重写一份架构/路径选择文件，已核实仅为新 PETSC_DIR/PETSC_ARCH；另两份 CFD 工作流文件更新了新路径及新编译库哈希。三项差异逐项登记，原文保留，最终核对无未解释差异。
   记录：validation/TRANSFER_VERIFICATION.json、migration/DOCUMENTED_RUNTIME_ADAPTATIONS.json。
2. 两进程 MPI Allreduce 正确；本地动态库依赖检查无缺失。
3. LAMMPS 两进程 run 0：3,488 节点、157,637 键；初始坐标误差为零；墙接触/早期流力回调、逐节点审计可运行；没有推进断裂仿真。
4. SimVascular 2023.05 原生 TetGen API 读入原 Y 表面并识别 [1,2,3,4] 面；检查范围是命令行原生网格接口，没有重新生成网格或验证 GUI/DICOM 功能。
5. PETSc CPU 单/双进程与 CUDA 自检通过。
6. 原 Y 网格（26,712 节点、146,230 四面体）CPU/GPU 各运行 2 个原生瞬态步，每次使用 2 个 MPI 进程；速度/压力最大相对 L2 差异 7.96446e-15，通过 1e-8 门槛。实际日志确认 GPU mpiaijcusparse/mpicuda 后端。
7. 当前实验全部 501 帧坐标及分量拓扑哈希与原记录一致；15 个原生 Y 资产哈希一致。血栓 1,927 节点、固定锚点 1,561 个，终态脱离 217 节点、208 个脱离分量。
8. 新环境实际生成 4 张 3840×2160 验证图并完成短视频编码/解码；X 轴投影角为 0 度；加载和阻力箭头使用原数据。保留原视频与原科学结果，验证图仅写入 validation/render。
9. CFD 配置入口、依赖导入、网格前置状态、原生网格启动器及新二进制哈希约束均通过。

后续开发入口
当前实验代码：pd_clot/lammps_fragmentation/early_flow_coupling.py
当前实验数据：verification/y_clot_point_cloud_and_early_flow_coupling_v1/
最新完整血管 4K 绘制代码：analysis/y_clot_whole_vessel_4k_v1/render.py
现有结果：results/y_clot_point_cloud_and_early_flow_coupling_v1/whole_vessel_4k_v1/
Y 工作流：/workspace/y_clot/deps/y_workflow/configs/stage4/workflow.yaml
网格启动器：/workspace/y_clot/bin/y-clot-simvascular
环境包装器：/workspace/y_clot/bin/y-clot-run
重编译脚本：migration/build_lammps.sh、migration/build_native.sh。
历史证据 JSON 中旧绝对路径是原始记录；现运行代码通过 pd_clot/portable_paths.py 映射到新目录，不应批量替换科学证据文件。

验证范围
本次仅迁移和验证环境，运行的 CFD 为每后端两步检查，不是重新验证稳态收敛或物理精度。现有流场仍为无血栓背景流场，碎片耦合沿用当前工程模型；没有新增双向流固耦合结论或全量断裂实验。
本地完整备份：/home/lzy/projects/server_migration_backups/y_clot_20261009
本次迁移记录：/home/lzy/projects/server_migrations/y_clot_20261009
