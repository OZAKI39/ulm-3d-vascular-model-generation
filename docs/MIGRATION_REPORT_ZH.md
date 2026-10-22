# BraVa 人体脑血管工程整理报告

本次实际建立了 `/home/lzy/projects/new_human_vessel_flow`：可以直接阅读和维护自己的 Python 模块，读取原流场和微泡结果，不再从旧 BraVa、打印或小鼠工程加载运行源码和输入。原文件未删除，服务器和 GitHub 未改动。

需要区分两件事：**工程已整理并进行实际测试；原来尚未通过的流场精度研究仍未通过。** 本轮没有新增正式 CFD、正式微泡批次或轨迹重放，也没有恢复 P2/P1、V3 或历史队列。

## 1. 选了什么版本

调查先于复制。先核对当前 `ACTIVE_DEVELOPMENT_SCOPE.md`、正式身份、Git 工作树、真实 import 依赖，再与 GitHub `sync/brava-flow-project-20261022` 对照，未按修改时间简单挑选。

| 来源 | 核实的身份 | 用途 |
|---|---|---|
| 本地 BraVa 工作区 | `a4b238b27437e1760b4f74fc219ede4f0180ee8f` | 正式输入、结果、原入口和当前诊断 |
| 本地共享粒子工作区 | `37c40d7c08ca6527c8379c0ffa81436ae708f59f`，包含既有未提交变化 | 以真实本地文件和依赖为准，保留当前粒子实现 |
| GitHub 同步快照 | `4e509991fd9c1b094e20318758aeb67a370fcbd2`，远端分支身份已核对 | 对照发布内容、恢复参考、生产原生源码 |
| 历史源码档案 | `3c3fb0498686fe647426bf67d067e90ebf189e52` | 溯源参考，未建立运行时依赖 |

8749 个发布文件中，8747 个与对应本地来源完全一致，2 个差异是发布时调整的 `.gitattributes`，没有被误判为科学实现更新。另审查了 19 个未发布的大文件，共 4.24 GB；没有将“未发布”理解成“可删除”。详细记录在 `docs/provenance/`。

首次选取并复制了 3434 个文件、4,314,176,716 字节，复制前先计算 SHA，再逐个核对。之后加入少量明确诊断内核和协议，并将无用研究启动器从新活动源码中剔除。最终核查的不可变资料为 **3289 个文件、4,313,292,920 字节**，全部逐字节一致；它们包含科研数据、分布契约、原生源码及历史证据。最终共核对 **3439 个原来源文件未变**。这些计数口径不同，不应相加。

## 2. 真实目录结构

完整层级摘要见 [DIRECTORY_TREE.txt](DIRECTORY_TREE.txt)，主要结构如下。未机械建立重复的 `results/frozen_flow`、几何配置或 archive 副本。

```text
new_human_vessel_flow/
├── src/
│   ├── human_vessel_flow/       # 明确的准备、读取、验证和绘图函数
│   │   └── reference_pipe/      # 现有诊断数值内核
│   ├── particle_3d/             # 一套权威粒子实现，含 population_v1/CUDA
│   ├── vascular_validation/     # 必要网格与流场检查
│   ├── flow_solver_support/     # 原求解日志、稳态、checkpoint、WSS 实现
│   └── sonovue_sampler.py       # 原始采样实现，源码 SHA 不变
├── native/                     # 生产原生源码、补丁、网格接口、构建身份
├── configs/                    # CFD 路径、原生依赖、固定灌注和当前契约
├── inputs/                     # BraVa 几何、唯一权威网格、粒径输入
├── cases/                      # 历史正式场、历史校准场，独立保存原身份
├── results/                    # 旧微泡、原图表报告、本次只读总览
├── diagnostics/                # 参考管、事件区域、近壁、P2/P1、局部边界
├── tests/                      # 原小型测试和路径/身份保护测试
├── validation/                 # 本轮实际日志、检查结果和准备示例
├── scripts/                    # 少量明确的检查和构建入口
├── docs/                       # 中文说明、迁移对应表、来源与清理建议
├── .venv/                      # 独立 Python 环境，不是旧环境链接
├── README.md
├── pyproject.toml
└── requirements-lock.txt
```

## 3. 源码与计算流程

从 [main.py](../src/human_vessel_flow/main.py) 开始，按 README 的模块表依次读：几何 → 网格 → 边界 → 原生求解调用 → 冻结场读取 → 微泡 → 检查 → 图表。

1. 几何来自原打印工程的 BraVa 四端口封闭内腔，只迁入必要 STL、SWC、端口和姿态记录。源单位 mm 转为计算单位 m，没有重新塑形。
2. 原 SimVascular/TetGen 调用与网格质量检查整理为明确函数。保留壁面和四端口标记，检查正体积、封闭性与对应关系。本次没有生成新网格。
3. 未来边界仍使用原“中心线阻力估算 → 出口压力”公式。按原算法复算 O1/O2/O3 为 **508.233436510425 / 0 / 24.309360181533066 Pa**。这只是设计值，实际三维分流为 `NOT_RUN`。
4. Python 准备 XML、调用 svMultiPhysics、读取 PETSc/稳态日志，实际 FEM 保留在原生 C++/Fortran。没有重写 CFD 算法。
5. 冻结场读取核对原节点顺序及每个四面体包含的节点。原 VTU 的四面体局部顶点顺序与权威网格不同，但单元对应一致。梯度与 WSS 使用权威连接顺序，不能仅因顺序不同就宣布几何或场不同。
6. 微泡使用原冻结场、SonoVue 分布、有限尺寸几何、阻力、接触、CPU 回退与 FP64 GPU 更新。当前 `particle_3d` 依赖闭包有 77 个模块，未凭旧阶段名称替换算法。
7. 参考管误差分解、事件诊断、P2 几何和资源判断内核保留为普通模块。重型研究启动器不作为当前入口，既有补丁及失败证据单独保存。

模块来源、顶层函数和类分别见 [MODULE_MAPPING.csv](MODULE_MAPPING.csv)、[FUNCTION_CLASS_MAPPING.csv](FUNCTION_CLASS_MAPPING.csv)。当前共 147 行文件/模块对应（含测试与网格脚本）、578 个符号对应；103 个对应文件逐字节不变，32 个有明确差异记录，12 个为新组织入口。完整差异在 [FINAL_SOURCE_ADAPTATIONS.patch](provenance/FINAL_SOURCE_ADAPTATIONS.patch)，未复制多个完整历史代码版本。

新入口还修正了两处结果报告顺序：无原生输出时把已有 FAIL 记录返回给命令行；原分流门槛提前到写入新冻结结果之前检查。门槛值和数值算法未改变，避免新入口在失败时给出成功退出或先写 PASS 的误导记录。

收尾核查曾把已经声明为 `adapted_source` 的 P2 几何 Python 文件误列入不可变资料，产生一次分类 FAIL；已修正检查器的类别处理，并将它纳入源码差异表。原始记录见 `validation/INITIAL_SOURCE_CATEGORY_CORRECTION.json`。没有修改任何科研数据或放宽哈希检查。

## 4. 正式算例、身份与保护

具体位置及关键 SHA 见 [CASE_INDEX.md](CASE_INDEX.md)。历史正式压力为 **556.9706555488448 / 0 / 17.13742715356202 Pa**，保持不变；与未来中心线估算压力分开。等流量校准场也完整保留，没有改称自然分流预测。

原网格为 135852 节点、692313 四面体。正式场出口分流为 33.3334236%、33.3229436%、33.3436328%，全局质量误差约 1.76e-16。原始 WSS 为 101074 个壁面片值，2 个恰为零，面积平均 3.0588756614426456 Pa。零值和近壁风险未被平滑或删除。

旧群体所有 145 个文件完整保留，包含 49,990 次状态快照、27,272,704 条粒子状态记录。源、终态、相互作用账本内容哈希和全部状态块文件哈希通过；ID、出生时间、年龄和状态时间关系通过。仍为 **1565 接受、458/527/437 出流、141 观察截止、2 数值失败**，`scientific_integrity_pass=false` 不变。

新固定协议原文件内容不变，灌注 1.240815660 s、随后观察 10 s、dt=0.5 ms、浓度 4e9 m⁻³。新准备入口默认使用它，旧覆盖停止入口没有成为默认。旧群体没有重新标为新协议成果。

所有新准备操作要求全新的任务目录；当前入口拒绝向冻结输入、旧算例和旧群体目录输出。新算例后处理只产生自己的 `NEW_FLOW_IDENTITY.json`，不覆盖原 `ACTIVE_FLOW.json`。历史 XML、日志和账本中的旧路径保留原文；当前运行通过显式新路径读取。source 源码保持可编辑；冻结资料以 SHA 清单和只读权限防止误写。

## 5. 本轮实际验证及能力边界

最终完整小测试为 **74 项通过、0 失败、0 跳过**，实际用时 81.14 s。独立性检查成功导入 **127 个本工程模块**，禁止旧工作区访问后仍能读取实际输入。真实 BraVa 场 **2146 个位置**的 CPU/GPU 对照通过，速度最大差为 **1.11e-16 m/s**，压力、梯度、涡量及应变率在这些点上的差均为 0；没有积分新轨迹。

最终逐项状态、实际测试数量和日志链接见 [VALIDATION_SUMMARY.md](VALIDATION_SUMMARY.md) 与 [DELIVERY_STATUS.json](../validation/DELIVERY_STATUS.json)，状态使用 `PASS / FAIL / NOT_RUN / BLOCKED / SKIPPED`。

交付文件的完整 SHA 清单见 [DELIVERY_FILES_SHA256.json](../validation/DELIVERY_FILES_SHA256.json)，复查结果见 [DELIVERY_HASH_CHECK.json](../validation/DELIVERY_HASH_CHECK.json)。清单排除虚拟环境、缓存及自身，科研文件和当前源码、配置、报告均在内。

已经执行的只读检查包括：源几何和 mm→m 数组精确对应、端口、网格质量、两套场 SHA、速度压力对应、出口流量、原始 WSS 重新求值与原数据比较、两套最终 MPI checkpoint 结构、全部旧微泡账本和状态身份。WSS 使用原 1e-12 比较容差，没有放宽。保存的 V2 分解数据重新计算后仍得到原分类和误差，没有新增 CFD 或轨迹。

CPU/GPU 测试使用原小型场景与原容差；真实 BraVa GPU 检查只比较少量位置的场求值。与完整正式生产运行、轨迹收敛和高阶精度认证分开记录。未执行内容不能因为环境装好或小测试通过而标成 PASS。

## 6. 独立性、环境与使用方法

日常 import 和读取只需要本工程以及明确安装的第三方工具。检查器用 Python 隔离模式导入每个模块，并在禁止读取旧工作区的条件下加载实际几何、流场、粒径分布和旧群体摘要。完整调用关系在 [SOURCE_IMPORT_GRAPH.json](../validation/SOURCE_IMPORT_GRAPH.json)，实际检查在 [INDEPENDENCE.json](../validation/INDEPENDENCE.json)。

没有添加 `sys.path`、执行旧源码文本或建立跨工程链接。`particle_3d/audit.py` 中两个旧路径只是未使用的历史报告文字；来源表和历史文档也保留旧绝对路径，不能将这些文字误认为运行依赖。一次性迁移工具以 `.py.txt` 保存在 provenance 中，不参加运行。`record_migration.py` 是明确读取旧清单的审计工具，普通使用不需要它。

Python `.venv` 独立新建。第三方 native 统一放在 `new_vessel_flow_solvers`，BraVa 求解器单独编译到其 `brava_svmultiphysics` 子目录，因为已有通用安装缺少 BraVa GPU 同步补丁。没有改动小鼠求解器，也没有混入实验 P2 补丁。详情与构建日志见 [NATIVE_ENVIRONMENT_ZH.md](NATIVE_ENVIRONMENT_ZH.md)。

实际使用命令、输出目录要求、模块阅读顺序和未来配置位置在 [README.md](../README.md)。本次还实际准备了三份不执行的计划：`validation/prepared_case_example`、`prepared_population_example`、`prepared_mesh_example`。可查看其 `PLAN.json`，其中 `executed=false`，不会冒充计算结果。

## 7. 图表说明

打开 [可视化总览](../results/reviews/migration_final/OPEN_RESULTS.html)，或查看 [PNG](../results/reviews/migration_final/frozen_evidence.png) / [PDF](../results/reviews/migration_final/frozen_evidence.pdf)。新图用原数据说明旧微泡分类、三出口流量和原面片 WSS 分布，已检查排版。

四视角原动画保持黑底、蓝泡、英文、无外框、坐标网格和大字号，原视频未重新编码。原显示排除了 2 条数值失败和 58 条风险轨迹，并执行逐帧前后 0.05 s、低于 2 mm/s 的低速筛选；这些都只是展示规则。原科学账本仍包含全部 1565 个接受微泡。动画流畅不能用于证明局部物理精度。

新总览的文件链接单独核查。原样保留的 19 个历史 HTML 中有 54 个旧相对链接目标未迁入，多数是保存网页的资源和已替代展示页面；清单在 `validation/HISTORICAL_HTML_LINK_REVIEW.json`。没有改写这些原始页面来伪造完整性，请从新总览与 [诊断资料索引](../diagnostics/README.md) 进入当前保留的结果。

## 8. 未解决问题与后续事项

- Gmsh V2：`REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC` 保留。
- 场误差分解：`MIXED_FIELD_ERROR` 保留。
- P1 轨迹敏感性：`P1_FIELD_TRAJECTORY_IMPACT_LARGE` 保留；诊断轨迹不是正式收敛结果。
- P2/P1：`STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT` 保留，曲面压力回填仍有失败，没有正式可接受高阶场。
- V3 仅提案；新固定时间窗正式群体没有运行。
- 新编译本地求解器没有做 CFD 数值对照，本轮也不声称新场精度或微泡收敛得到认证。

本轮允许范围内的整理、读取和小测试完成后，后续是否开始新科学计算应由新的任务决定。旧目录处理只输出 [清理审核清单](OLD_WORKSPACE_CLEANUP_REVIEW.md)，**原工程删除数为 0**，没有自动恢复任何旧任务。
