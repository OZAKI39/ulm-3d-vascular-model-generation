# SECOND ACTIVE-CONTEXT CLEANUP 结果

**已完成清理、必要smoke和文档；生产物理零修改，科学结果零重算。完成独立提交后停止，等待下一条指令。**

本地工程和服务器工程同步执行；未启动Gmsh、新mesh、CFD、oracle、trajectory或population。内存中的单元测试几何、静态CUDA代数/状态测试不构成科学仿真。

## Q1：删除了多少

| 指标 | 本地 | 服务器 |
|---|---:|---:|
| 已删除路径（包含链接） | 2,090 | 4,303 |
| 其中tracked删除 | 590 | 579 |
| 已删除源码LOC | 71,918 | 58,839 |
| raw stdout/stderr日志 | 315 | 345 |
| 额外JSONL流 | 18 | 17 |
| 中间mesh/array等数据文件 | 199 | 265 |
| 删除逻辑bytes（毛额） | 353,545,022 | 1,404,637,843 |
| 退休mesher build/install/source bytes | 0 | 1,032,910,993 |

本地毛删除约0.329 GiB，服务器约1.308 GiB；其中服务器退休mesher构建/安装/source约0.962 GiB。两端含重复副本，不能相加称作唯一科学数据量。新增审计清单抵消约10 MB，净变化见下表。原生20M级Omega_h网格只存在于失败进程内存，没有20M最终mesh文件，本轮未虚报删除该文件。

| Checkout规模 | 本地：前 → 后 | 服务器：前 → 后 |
|---|---:|---:|
| 最终Git tree tracked数 | 793 → 232 | 793 → 232 |
| 普通文件数 | 11,856 → 9,791 | 55,569 → 51,387 |
| case checkout一方源码文件数 | 603 → 148 | 1,810 → 1,366 |
| Python LOC | 65,530 → 14,744 | 97,201 → 46,464 |
| C/C++ LOC | 23,416 → 2,336 | 433,439 → 425,389 |
| MD/HTML报告数 | 220 → 66 | 294 → 142 |
| JSON数（含科学身份） | 3,543 → 2,947 | 3,558 → 3,014 |
| 日志/CSV/JSONL数 | 579 → 246 | 651 → 289 |
| working-tree逻辑GiB（排除Git DB） | 22.6949 → 22.3748 | 28.2308 → 26.9318 |
| allocated GiB（排除Git DB） | 22.7212 → 22.3966 | 28.2709 → 26.9626 |

计数口径：不跟随符号链接，排除所有`.git`数据库；源码LOC不包含runtime/env/vendor中的第三方安装。服务器`gpu_solver_fix`中的生产求解器源码构建保留，因此服务器C/C++存量较大。该表的源码数是checkout内一方源码（含保留baseline快照），不能混称为全部生产依赖数。显式活动项目依赖闭包90个文件全部保留，新增1个测试文件；共享粒子330个源码/测试SHA条目保持。日志总数含CSV/JSONL，因此与上表raw stdout/stderr数不同。JSON仍有约3000份，多为保留的正式run/seed/逐粒子身份和科学证据，不以减少计数为由删除。

指标来自清理后快照；最终报告/交付清单的文本封存可能增加数KB。完整口径和时间在`production_workflow/validation/second_cleanup/{BEFORE_METRICS,AFTER_METRICS}.json`。服务器删除前已缺14个tracked文件，其中11个属此次退休路径；交付后其余3个历史缺失继续注明，不冒充本轮删除或新回归。

## Q2：从活动checkout移除的路线

P1 brute-force refinement的旧queue、diagnostic脚本、执行副本、monitoring；local P2/P1与BC-B一次性诊断实现；frozen Jacobian计划此前已停止、无活动实现，本轮索引明确标为NOT_RUN；MMG anisotropic/runtime专用源码、输入、build/install与D0/D1/D2失败workspace；Omega_h适配器、build/source、R1失败workspace；NASA/refine预检与几何转换；MeshSim临时probe；无后续Gmsh引用的TetGen scalar volumetric metric/support-cloud glue。当前生产SimVascular/TetGen/PETSc/MPI未删除。按不删除Git历史的要求，保留原有refine_build/source和omega_h_check/source两份隐藏.git数据库；工作区源码/build已清空，数据库不计active checkout，也未新建archive。

每条路线仅留一个最终中文报告和一个紧凑summary（MeshSim/refine并入Omega_h最终报告，MMG两阶段合并）。历史图保留单一PNG，每报告不超过2张。关键日志保存有行号的短摘录，不复制整份日志。入口：[退休索引](RETIRED_EXPERIMENTS_INDEX.md)、[证据与SHA](RETIRED_EXPERIMENT_EVIDENCE.json)。

## Q3：继续保留哪些历史结果，为什么

- 正式BraVa geometry/mesh/flow、所有原始1565接受状态、source/seed、出口/删失/失败账本、运行身份及固定灌注合同。旧结果仍为458/527/437出流、141观察截止、2数值失败，不能称全成功或新协议结果。
- 高保真R1 PLC、frozen probes、P0/P1原始CFD、合法NW1 mesh、Poiseuille公式、strict containment、u/G/s/curl合同、dt/oracle复用接口、mesh conversion/svMultiPhysics adapter。P0/P1 field均未通过全部门槛；NW1 mesh合法但normal-height未达目标，未计算NW1 CFD。
- 昂贵的旧研究CFD/mesh/singleton状态，以及仍被保留replay helpers读取的FIELD_ARRAYS等审计数据；这些可能用于正式结论复核，不能假定可低成本重建。
- 正式动画、兼容链接、共享粒子源码/全部当前测试、vendor及生产服务器环境。

## Q4：production代码与科学配置有没有变化

**NO PHYSICS CHANGE。** 没有修改任何已有生产/validation Python或C/C++源码；唯一新增源码是5个小型reference interface tests。删除的是退休实现；更新的是范围说明、历史报告归并、索引及当前清理状态。没有改物性、边界条件、mesh/field、source、seed、near-wall/MB–MB/contact/adaptive stepping或默认dt。

删除后哈希核对：本地693项保护项目文件+330项粒子源码/测试全部一致；服务器692项存在的保护项目文件+150项实际部署粒子文件全部一致。服务器1个旧展示HTML和180项旧测试本来就未部署，未声称核验了不存在的文件。所有未获修改授权的保留路径size/mtime无漂移。见`PRESERVATION_LOCAL.json`、`PRESERVATION_SERVER.json`。

必要smoke：本地23 passed、服务器10 passed。particle_3d/population导入、真实GPU模块导入、manufactured分析导入、strict containment（域内/边界/域外及定位器误判）、快速核心阻力/固定合同测试均通过。GPU测试实际使用CUDA，但不推进轨迹。正式135852节点/692313 tetra flow loader及P0 184712 tetra的svMultiPhysics转换输入/GlobalNodeID面片只读smoke通过。未调用旧零源population验证脚本。

AST import/CLI引用闭包和git grep完成；85项宽松字面命中已逐项核对，均是仍保留的同名canonical文件、通用output词或未来运行创建的executed副本，没有活动代码加载已删除实现。`managed_job.py`保留通用监督接口，旧TASK删除；新授权任务必须新建TASK，不能重启旧任务。Markdown/HTML链接核验通过（`:line`本地文件跳转按行号语法解析）。

## Q5：下一轮默认读取

先读[ACTIVE_DEVELOPMENT_SCOPE.md](ACTIVE_DEVELOPMENT_SCOPE.md)。然后只按任务读共享production particle、`production_workflow/scripts/{poiseuille,pipe_probe_geometry,particle_field_inputs,analyze_pipe_calibration,manufactured_boundary,prepare_pipe_cfd,pipe_measurements}.py`，三个命名粒子输入合同/审计/基线文件、manufactured PROTOCOL/PARTICLE_PROBE_CONTRACT、当前tests。未来Gmsh位置预留为`production_workflow/gmsh_reference/`，本轮未开发。不要递归读旧路线或本次大审计清单。

## Q6：MANUAL_REVIEW_REQUIRED

存在，全部KEEP。UNKNOWN共本地229、服务器243个路径；其中>10 MB有本地17、服务器19个，主要是旧mesher probe体网格及wall-query NPZ。另有明确保留的昂贵研究数据可在将来专门审计，但本轮不删除。完整逐文件列表在[LARGE_FILE_REVIEW.csv](LARGE_FILE_REVIEW.csv)。例如18.75µm probe volume约35.7 MB、多个wall-query数组约28.6–31.2 MB；目前未证明它们都是可舍弃的重复证据。

最大存量仍是科学数组/原始结果和生产依赖。这是额外科学数据保护的结果；本轮主要消除了代码、重复报告、构建、monitoring与中间实验上下文。

## 最大剩余20文件（逐端）

### local

| 文件 | MB（10⁶bytes） |
|---|---:|
| `population_v1/cfd_wall_audit/all_cells.csv` | 430.780 |
| `mesh_particle_validation_v2/04_field_audit/M3_fine_h075_frame180/FIELD_ARRAYS.npz` | 365.189 |
| `population_v1/cfd_wall_audit/wall_adjacent_cells.csv` | 357.078 |
| `mesh_particle_validation_v2/04_field_audit/M2_fine_frame180/FIELD_ARRAYS.npz` | 348.422 |
| `mesh_particle_validation_v2/04_field_audit/M3_medium_frame180/FIELD_ARRAYS.npz` | 295.949 |
| `mesh_particle_validation_v2/04_field_audit/M2_medium_frame180/FIELD_ARRAYS.npz` | 289.534 |
| `mesh_particle_validation_v2/04_field_audit/M1_coplanar_fine_h020/FIELD_ARRAYS.npz` | 288.165 |
| `mesh_particle_validation_v2/04_field_audit/M1_coplanar_fine/FIELD_ARRAYS.npz` | 219.235 |
| `population_v1/cfd_wall_audit/recomputed_arrays.npz` | 207.728 |
| `mesh_particle_validation_v2/04_field_audit/M1_fixedplc_fine/FIELD_ARRAYS.npz` | 185.852 |
| `mesh_particle_validation_v2/04_field_audit/M1_coplanar_medium/FIELD_ARRAYS.npz` | 167.411 |
| `mesh_particle_validation_v2/04_field_audit/M1_fine/FIELD_ARRAYS.npz` | 166.545 |
| `mesh_particle_validation_v2/04_field_audit/M1_fixedplc_medium/FIELD_ARRAYS.npz` | 158.868 |
| `mesh_particle_validation_v2/04_field_audit/M1_medium/FIELD_ARRAYS.npz` | 153.376 |
| `mesh_particle_validation_v2/04_field_audit/M0/FIELD_ARRAYS.npz` | 144.855 |
| `mesh_particle_validation_v2/04_field_audit/M3_fine_h075_frame180/AUDIT_VOLUME.vtu` | 123.867 |
| `mesh_particle_validation_v2/04_field_audit/M2_fine_frame180/AUDIT_VOLUME.vtu` | 118.156 |
| `microbubble/data/gpu_mesh_input.npz` | 109.639 |
| `mesh_particle_validation_v2/04_field_audit/M1_coplanar_fine_h020/AUDIT_VOLUME.vtu` | 106.912 |
| `mesh_particle_validation_v2/04_field_audit/M3_medium_frame180/AUDIT_VOLUME.vtu` | 100.496 |

### server

| 文件 | MB（10⁶bytes） |
|---|---:|
| `population_v1/cfd_wall_audit/all_cells.csv` | 430.780 |
| `mesh_particle_validation_v2/04_field_audit/M3_fine_h075_frame180/FIELD_ARRAYS.npz` | 365.189 |
| `population_v1/cfd_wall_audit/wall_adjacent_cells.csv` | 357.078 |
| `mesh_particle_validation_v2/04_field_audit/M2_fine_frame180/FIELD_ARRAYS.npz` | 348.422 |
| `mesh_particle_validation_v2/04_field_audit/M3_medium_frame180/FIELD_ARRAYS.npz` | 295.949 |
| `mesh_particle_validation_v2/04_field_audit/M2_medium_frame180/FIELD_ARRAYS.npz` | 289.534 |
| `mesh_particle_validation_v2/04_field_audit/M1_coplanar_fine_h020/FIELD_ARRAYS.npz` | 288.165 |
| `mesh_particle_validation_v2/runtime/simvascular/svExternals/lib/python3.5/site-packages/tensorflow/python/_pywrap_tensorflow_internal.so` | 250.875 |
| `mesh_particle_validation_v2/04_field_audit/M1_coplanar_fine/FIELD_ARRAYS.npz` | 219.235 |
| `population_v1/cfd_wall_audit/recomputed_arrays.npz` | 207.728 |
| `mesh_particle_validation_v2/04_field_audit/M1_fixedplc_fine/FIELD_ARRAYS.npz` | 185.852 |
| `mesh_particle_validation_v2/04_field_audit/M1_coplanar_medium/FIELD_ARRAYS.npz` | 167.411 |
| `mesh_particle_validation_v2/04_field_audit/M1_fine/FIELD_ARRAYS.npz` | 166.545 |
| `mesh_particle_validation_v2/04_field_audit/M1_fixedplc_medium/FIELD_ARRAYS.npz` | 158.868 |
| `mesh_particle_validation_v2/04_field_audit/M1_medium/FIELD_ARRAYS.npz` | 153.376 |
| `mesh_particle_validation_v2/04_field_audit/M0/FIELD_ARRAYS.npz` | 144.855 |
| `mesh_particle_validation_v2/04_field_audit/M3_fine_h075_frame180/AUDIT_VOLUME.vtu` | 123.867 |
| `mesh_particle_validation_v2/04_field_audit/M2_fine_frame180/AUDIT_VOLUME.vtu` | 118.156 |
| `mesh_particle_validation_v2/runtime/simvascular/svExternals/lib/libQt5WebEngineCore.so.5.11.3` | 113.569 |
| `microbubble/data/gpu_mesh_input.npz` | 109.639 |

## 交付与停止

独立分支与提交标题：`cleanup/prune-retired-meshing-experiments`；清理前恢复点`5622c1a9549c19306b28eb5479c87600d76997c0`。使用git rm删除tracked路径，未创建archive/old/legacy/backup/retired代码目录，未重写历史或force push。提交正文记录removed retired meshing experiments / preserved authoritative reports/results / no production physics changes / no scientific result recomputation。

本地已有未跟踪的P0/P1执行源码副本与vendor继续保持原状，没有混入本次提交；服务器预先缺失项明确记录。外层`BRAVA_CONTEXT_HANDOFF.md`与`brain_report.md`只更新当前入口块，没有改历史正文。服务器`brava_production`为EXITED，旧项目任务STOPPED/EXITED，无自动下一步。

[保留科学图件入口](production_workflow/OPEN_RESULTS.html)。本轮没有生成伪造的新网格或新轨迹图。完成本次提交后停止，等待用户指令。
