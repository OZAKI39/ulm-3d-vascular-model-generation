# 小鼠脑血管流场与微泡工程整理交付报告

日期：2026-10-22。项目：`/home/lzy/projects/new_mouse_vessel_flow`。

本轮已实际完成有效源码迁移、独立 Python 环境安装、正式数据无损复制、配置归一、运行入口及允许范围内的验证。默认参考为 ROI-only；另外三套正式流场、Stage-Q 与既有 WSS 证据分别保留。没有重新求解正式 CFD、生成新网格、重算正式微泡批次、连接服务器、推送 GitHub 或删除旧工程。

**交付状态：工程整理与本地核验完成；原生 CFD/几何工具环境尚待单独配置，不能据此宣称已完整重跑科研工作流。**

- 31,419 个迁入数据/证据文件，共 4,864,008,448 字节，逐文件 SHA256 核验通过。
- 原位置 34,645 个已盘点文件复查无改动。
- 189 个 Python 模块导入通过；当前源码身份 191 项、原始受保护源码身份 114 项核验通过。191 项包括包标记文件和 C++ 内核，与 189 个导入模块不是同一种计数。
- 测试为 **116 PASS、3 BLOCKED（pytest 显示 skipped）、0 FAIL**。三个阻碍均是迁移前已经缺少的历史 real-scene 测试数据，未补造。
- 五套冻结场/原生 checkpoint、三套共 3500 条正式微泡轨迹、ROI 网络/H0 和 WSS 只读恢复均通过对应检查。

入口：[README](../README.md)；[图表与原媒体入口](../results/refactor_review/index.html)；[本轮数值核验 JSON](../validation/delivery/validation.json)；[测试状态](../validation/TEST_STATUS.json)。

## 1. 本地基线与来源确认

迁移以本地真实调用链为准，没有按修改时间选择模型。先检查发布分支 `sync/mouse-vascular-flow-project-20261022` 的本地已发布快照、README、算例说明及同步清单，再逐文件对照真实来源。快照 commit 为 `987d0638944f09e1acc1f415b26bf5b0479c3f1d`。

对比的 34,645 项中，仅三份 `.gitattributes` 存在发布包装差异；没有发现这批清单内科学源码/数据与本地内容不同。此结论限定于盘点范围，不等于所有 `/home/lzy/projects` 文件均已发布。发布排除清单、用途未确认项和当前选择依据均保留。

| 真实来源 | 采用依据与身份 |
|---|---|
| `_shared/particle_workflow/particle_3d` | 两套 1500 条算例的 `production_entry.py → campaign.py` 及原生几何内核调用、H0 契约、114 文件保护清单；Git 基线 `37c40d7c08ca6527c8379c0ffa81436ae708f59f` |
| `_shared/active_network_development/ulm_3D_vascular` | 当前 H0、ROI-only、best-feasible 实际使用的 21 个网络模块；基线 `69c6af94eabb87ba6db0823134c4b92b30a069db` |
| `_shared/geometry_network` | 用于核对旧几何/网络及路径关系；基线 `fd21a850a16d0ba17ef3d521123badd53864a3a2`；没有把这套较旧网络另复制为运行实现 |
| `_shared/formal_3D_flow_solver/FEM_SimVascular` | 原网格、FEM/WSS 自有支持、配置和 native 来源；外层基线 `a4ff70fa77d49a319cb701617a30c243fc3c37e2` |
| `vascular_printing` 中必要小鼠几何部分 | 实际被调用的 ROI、管腔和表面工具及 fMOST 源样本，记录文件 SHA；不复制打印工程 |

上述旧 Git 工作树原本含大量未提交/历史清理状态，本轮没有提交、重置或整理它们。具体状态见 [git_baseline.json](audit/git_baseline.json)；逐文件本地/发布对比见 [local_vs_published.json](audit/local_vs_published.json)。

## 2. 实际目录结构

```text
new_mouse_vessel_flow/
├── README.md / pyproject.toml / requirements-lock.txt / .gitignore
├── mouse_vessel_flow/                 当前自有代码唯一运行权威
│   ├── main.py / paths.py / provenance.py
│   ├── geometry.py / meshing.py / network.py / boundary.py / fem.py
│   ├── flow_field.py / wss.py / validation.py / visualization.py
│   ├── geometry_tools/               必要几何依赖闭包
│   ├── network_models/               当前 1D/0D 网络实现
│   ├── mesh_tools/                   原生 TetGen 与网格检查支持
│   ├── solver_support/               原生运行监控、checkpoint、WSS
│   ├── particles/                    原源模型、插值、有限尺寸动力学与分析
│   └── presentation/                 原流场/表面/点诊断展示工具
├── configs/
│   ├── geometry/ / flow/ / particles/
│   ├── environment.json
│   └── runtime_source_identity.json
├── inputs/
│   ├── geometry/ / network/
│   └── reference/
│       ├── fem/                      规范网格、原始 2 mm/s 参考和 checkpoint
│       ├── particles/ / sonovue/
│       └── source_identity/          原始源码字节证据 ZIP，只读校验，不执行
├── cases/{roi_only,best_feasible,network_h0,mean_2mmps,stage_q}/
│   └── case.json / README.md          参数、身份和数据入口，不复制求解源码
├── results/
│   ├── roi_only/ / best_feasible/ / network_h0/ / mean_2mmps/ / stage_q/
│   └── refactor_review/              本轮只读统计图、PDF、HTML
├── validation/
│   ├── recorded_cases/ / wss_audit/ / wss_validation_v2/
│   ├── delivery/validation.json / TEST_STATUS.json / pytest_delivery.xml
│   ├── original_small_scene_baseline.json / small_scene_equivalence.json
│   └── prepared_case_smoke/           新输入准备检查，无求解
├── external/
│   ├── README.md
│   ├── svMultiPhysics/{source/,identity.json,local_changes.patch}
│   └── solver_support/
├── tests/                            迁入的原测试与本轮必要回归
├── scripts/                          一次性迁移/审计工具，不参与科学运行
├── docs/{WORKFLOW_ZH.md,SOURCE_SYMBOL_MAP.csv,CLEANUP_REVIEW_ZH.md,audit/}
├── archive/                          精选历史证据，不参与当前计算
└── .venv/ / .cache/ / .pytest_cache/ / mouse_vessel_flow.egg-info/
                                      本轮新建环境与可再生本地文件，Git 忽略
```

没有通过软链接或硬链接维持旧工作区依赖。规范输入集中保存；少数原生算例/成品媒体包内的相同网格副本按原相对路径身份契约保留，详见 [内容重复审查](audit/duplicate_data_review.json)。这不是全历史目录镜像，也不声称数据全局逐字节去重。

## 3. 主要代码与迁移对应

| 旧来源/符号 | 新位置/入口 | 职责与处理 |
|---|---|---|
| 几何 `utils.rodent_vasculature`、sampling、cfd_lumen/preprocess/surface_prepare | `geometry.py`、`geometry_tools/` | 60 个实际依赖模块；读取源样本、ROI、管腔与表面输入，未引入打印模具/支撑 |
| `network_1d0d/roi_only_hydraulics.py: load_roi_cache, regress_h0` | `network_models/roi_only_hydraulics.py`；`network.roi_cache/validate_design` | ROI 图、阻力和 H0/设计回归 |
| `roi_boundary_design.py: design` 与 `parameterized_fem_handoff.py`、`fem_h0_case.py` | 同名 `network_models/` 模块；`boundary.py` 显式导出 | 压力设计、延长段交接、FEM 参数；不混用三种工况 |
| 原 SimVascular 网格生成脚本 | `mesh_tools/tetgen_native.py: generate`；`meshing.launch_embedded` | 保留原 TetGen 生成语句，通过原生嵌入式解释器运行 |
| 原 `solve_a_h0_fem_remote.py` 生产监控入口 | `solver_support/native_run.py: run_native`；`fem.prepare_case/run_case` | 显式传入本地 case/build，保留原停止、checkpoint、验收和求解逻辑 |
| `particle_3d/field.py: FrozenFEMField` | `particles/field.py`；`flow_field.load_case` | 原冻结场插值，增加明确 case 路径及身份入口 |
| 原 SonoVue sampler 与 continuous infusion source | `particles/sonovue_sampler.py`、`sonovue_adapter.py`、`injection.source_for_case` | 单一 sampler；显式 import；原始字节身份和新运行代码身份分别检查 |
| `formal_dynamics_p9a5`、`particle9a_motion`、`particle65_motion` 等 | 同名 `particles/` 模块；`dynamics.py` | 原有限尺寸积分/接触依赖仍然有效，保留阶段名称与算法 |
| 正式 `campaign.py`/`production_entry.py` 中环境和任务绑定 | `particles/dynamics.py`、`trajectories.run_saved_cohort` | 多算例共享一个入口，保存队列和契约保持原样；新任务记录新源码身份 |
| 正式 `native_adapter.py` / `closest_native.cpp` | `particles/` 同名文件 | 保留原 CPU C++ 查询及绑定策略，编译产物写入新工程缓存 |
| 原 WSS facet/traction 实现 | `solver_support/wss.py`、`wss_case.py`；`wss.py` | 原双精度网格、原公式和显示平均定义 |
| 原 flow/surface/streamline 展示脚本 | `presentation/` | 显式新输出目录、当前输入路径；绘图与新点诊断分开 |
| 原散落核验入口 | `validation.py`、`particles/analysis.py`、`visualization.py` | 冻结场、原始事件/轨迹/WSS 只读验证和报告 |

完整 **旧文件/函数/类 → 新文件/函数/类** 对照在 [SOURCE_SYMBOL_MAP.csv](SOURCE_SYMBOL_MAP.csv)，包括 905 个顶层函数/类；文件来源、旧/新 SHA、Git 基线与修改状态在 [source_map.json](audit/source_map.json)。171 个迁入源文件中 141 个字节一致；880 个函数/类 AST 完全相同。剩余差异逐项保存于 `audit/source_diffs_final/`，涉及命名空间、路径、环境/输出接口及身份校验，不能笼统声称全部源码字节不变。

原 114 受保护文件完整字节作为 `inputs/reference/source_identity/original_source_identity.zip` 保存；它们不作为第二套可执行实现。新的运行源码身份单独固定在 `configs/runtime_source_identity.json`。原正式计算契约和日志没有改写成新代码身份。

依赖、I/O、子进程、环境访问见 [dependency_inventory.json](audit/dependency_inventory.json)；导入图为 189 节点、无发现循环，见 [import_graph.json](audit/import_graph.json)。原生适配的有限私有函数绑定仅引用当前包，未动态加载历史文件。完整上游测试资源和旧虚拟环境未迁入。

## 4. 科学流程与保持不变的边界

详细说明见 [WORKFLOW_ZH.md](WORKFLOW_ZH.md)，源码阅读顺序为 `main → geometry → network → boundary → fem → flow_field → particles → validation/visualization`。

1. **几何**：保留 fMOST 样本 `fMOST_0_5_6_0_0_6_0001_02_01`、anchor 003274、真实 ROI/端口与表面元数据。SWC 坐标/半径为 µm，计算转换到 SI；Ultraliser 表面输入沿用 `0.91 × source_radius`，网络仍使用原始半径。本轮未重建几何。
2. **网络与边界**：ROI 109 节点/108 边，原 `R = 8μ/π ∫ds/r⁴`、节点守恒、H0、ROI-only 和 best-feasible 设计保持。真实切口到 CFD cap 的延长段压力转换不变。ROI-only 不用 CFD 反向校准；0D 的目标等分不能替代实际 3D 出口流量。
3. **FEM**：沿用规范四面体和面标签，P1/P1 + VMS、Newtonian `ρ=1056 kg/m³`、`μ=0.00345312 Pa·s`、原广义 α 步进和 PETSc 参数。Python 负责准备和监控，原生 svMultiPhysics 才执行求解。保留原零初值、稳态停止及物理验收逻辑。
4. **冻结场**：按 case 验证精确 SHA、点/单元顺序、GlobalNodeID 与边界映射，再复用原 `FrozenFEMField` 插值。未用平滑/修正场替换原正式结果。
5. **微泡**：保留原 SonoVue 数据、D≤4 µm 条件分布、来源顺序和种子、单次入口准入、有限尺寸阻力、近壁与接触规则、原积分步长、3/6/12 s 观察扩展、出口/静止/失败分类。RBC 和可选 LAMMPS bridge 作为独立研究功能保留依赖边界，未引入 RBC–微泡耦合。
6. **WSS**：原 `τ=(I−nnᵀ)μ(∇u+∇uᵀ)n`，使用原 float64 网格，区分 raw facet 值与节点显示平均。流量守恒与图像平滑不构成局部离散精度证明。
7. **展示**：本轮图表读取原始已保存结果。原动画原样保留；点流线诊断、有限尺寸微泡轨迹和统计图具有不同含义，没有互相替代。

## 5. 正式算例索引与身份

所有算例通过 `cases/<name>/case.json` 显式追踪几何配置、网格/面标签、求解 XML、材料、冻结场、环境、粒子契约与输出根目录。共同规范网格为 `inputs/reference/fem/frozen_reference/SV_MESH/mesh-complete.mesh.vtu`，70,363 节点、371,402 四面体；SHA `1a7a69dcba52f0475b2280775921e4dedcdfefa67b7965096d86cb9dba38bdf9`。

| 配置/说明 | 原工况 | 正式流场与结果入口 | 微泡状态 |
|---|---|---|---|
| [roi_only](../cases/roi_only/README.md) | mean-2p0-mmps-A-ROI-only-balanced-pressure-v1 | `results/roi_only/flow/frozen_flow/`；`results/roi_only/particles/microbubble_roi_only_dt0p5ms_n1500/` | 1500，dt=0.5 ms，默认参考 |
| [best_feasible](../cases/best_feasible/README.md) | mean-2p0-mmps-A-best-feasible-balance-v1 | `results/best_feasible/flow/frozen_flow/`；`results/best_feasible/particles/` | 1500，dt=0.5 ms，独立方案 |
| [network_h0](../cases/network_h0/README.md) | mean-2p0-mmps-A-H0-pressure-v1 | `results/network_h0/flow/frozen_flow/`；`results/network_h0/particles/` | 500，dt=1 ms，历史最终动画交付未完成 |
| [mean_2mmps](../cases/mean_2mmps/README.md) | mean-2p0-mmps | `results/mean_2mmps/flow/frozen_flow/`；`inputs/reference/particles/` | 原始基线与参考验证，未冒充新批次 |
| [stage_q](../cases/stage_q/README.md) | fem_reference_stage_q | `results/stage_q/flow/` | 历史较低入口速度参考，不能称为 2 mm/s |

| 场 | 完整冻结 SHA256 |
|---|---|
| ROI-only | `fb3c6ad0815d156b09b289fc347766486f24ddb521b8e51e8bbc4cf1fd2ace12` |
| best-feasible | `fcc692caa74c70d7ecbf9ae6662d29d45abbae926b04783aa52886d382a4bdb2` |
| Network-H0 | `064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4` |
| 原始 2 mm/s | `129ebb77396550b300b20ce29cde7b03919e6037b7a81de60fe75d6c6efb616d` |
| Stage-Q | `373ff549430cab710d57eb4406f2e7588f7a5bb68ebaa68cf32f8da5662e26f1` |

原生 checkpoint 也逐一读取头信息、节点数、step/time 与状态完整性，均为 step 71、70,363 节点，哈希各不相同。详细 SHA 及原路径→新路径见 [delivery/validation.json](../validation/delivery/validation.json) 和 [data_manifest.json](audit/data_manifest.json)。不能按相似文件名混同不同场。

## 6. 实际验证结果

| 检查 | 状态 | 证据与适用范围 |
|---|---|---|
| 全量迁入文件与原位置保护 | PASS | 31,419 文件 SHA 相同；原 34,645 文件无变化；`data_manifest.json`、`original_files_after.json` |
| 从 `/tmp` 导入/运行 | PASS | 189 模块，旧项目模块 0、项目软链接 0；`validation/delivery/validation.json` |
| 当前/历史源码身份 | PASS | 当前 191 项；原保护清单 114 项；两种身份分开记录 |
| 网络/H0/ROI-only | PASS | 109/108，H0 端口流量最大相对 Qin 差约 `1.75e-14`，压力设计保存量对比最大差 0 |
| 五套冻结场、网格与 checkpoint | PASS | 精确哈希、原 P1 插值、梯度、有限数值与面流量；三套当前场与原验收记录一致 |
| ROI-only WSS | PASS | 45,221 wall facets；WSS、切向牵引、法向、面积、父单元最大绝对差均 0 |
| 三套原微泡数据 | PASS | cohort/事件/初始状态、COMPLETE 清单哈希、轨迹格式、dt、出口/静止/失败分类逐项核查 |
| 原小场景动力学等效 | PASS | wall/pair，dt=0.004 s、时长 0.032 s，状态、接触、接受/拒绝步及诊断规范 JSON 字节完全相同 |
| 正式 C++ 几何内核 | PASS | 本地编译；2048 组最近点查询与原 Python 计算逐位一致；当前算例任务绑定通过 |
| 最小源事件/准入 | PASS | ROI-only/best-feasible 原 ledger 对比；未启动新正式灌注 |
| TetGen 生成语句与 Python 3.5 语法 | PASS | AST 比较与静态兼容检查；不代表 native 网格生成执行通过 |
| 整体 pytest | PASS / BLOCKED | 116 通过，3 因原缺失 fixture 跳过；`pytest_delivery.xml`、`TEST_STATUS.json` |
| CFD 输入准备 | PASS | 只复制输入、校验路径/哈希，`prepared_case_smoke`，solver_executed=false |
| 原生 CFD 环境 | BLOCKED | 未在独立环境安装原 solver/PETSc/CUDA 构建 |
| 原生 TetGen / Ultraliser / VMTK 执行 | BLOCKED | 工具未独立安装；未生成新几何/网格 |
| 新正式 CFD、微泡批次、完整动画重渲染 | NOT_RUN | 本轮不启动；旧冻结结果保持不变 |
| 细档血管/O2 压力敏感性 | NOT_RUN | 原状态 `SKIPPED_BY_USER` 保留，未恢复 |

轻量 wall/pair 基线曾在独立子进程中显式读取旧源码进行一次对照采集，使用的是本轮新 Python 环境。这是有记录的迁移验证，不是隐藏运行依赖。当前测试只读本地保存的 `original_small_scene_baseline.json`，当前运行代码无需原目录。证据见 [small_scene_equivalence.json](../validation/small_scene_equivalence.json)。

三个 BLOCKED 测试为 `test_every_position_is_time_integrated`、`test_no_skipped_or_duplicate_time`、`test_accepted_and_continuous`。所需原 `particle6_5/data/07_real_v1_dt1.json` 在本地与发布快照均已缺失；相关断言未放宽，也未重新生成历史输入。已执行 wall/pair 与真实保存轨迹检查不能替代这三项原测试。

迁移中发现的缺失 import、测试 fixture 路径、H0 步长读取位置已按实际来源修正。插值测试曾误用 float32 质心而产生对照误差，修正的是测试取样精度，没有修改冻结场、插值或放宽容差。保留首次失败和后续日志，不将初次未通过隐去。

## 7. 微泡结果与可视化分析

| 工况 | 接受轨迹 | O1 | O2 | O3 | 总出流 | 接触支持静止 | 执行失败 |
|---|---:|---:|---:|---:|---:|---:|---:|
| ROI-only | 1500 | 787 | 61 | 467 | 1315 | 185 | 0 |
| best-feasible | 1500 | 0 | 133 | 1182 | 1315 | 185 | 0 |
| H0 | 500 | 0 | 139 | 291 | 430 | 70 | 0 |

ROI-only 和 best-feasible 的总出流数相同，但出口分布显著不同：前者以 O1/O3 为主，后者以 O3 为主。H0 使用另一正式队列和 1 ms 步长，不能把三组统计直接归因为仅一个变量的效应。185/70 条是模型判定的接触支持静止，不是全部成功出流，也不是本轮新增数值失败。

已生成 [300 dpi PNG](../results/refactor_review/frozen_campaign_comparison.png)、[PDF](../results/refactor_review/frozen_campaign_comparison.pdf) 和 [HTML](../results/refactor_review/index.html)，内容为出口/静止计数及接受队列粒径经验 CDF。统计来源在 [source_summary.json](../results/refactor_review/source_summary.json)；HTML 连接本地原媒体。CDF 相近或重叠不代表速度场/轨迹相同，图表没有做轨迹平滑或状态修改。

## 8. 已执行启动命令与后续入口

独立 `.venv` 基于 Python 3.13.11 新建，`include-system-site-packages=false`，没有复制旧 venv。主要版本为 NumPy 2.5.3、SciPy 1.18.1、VTK 9.7.0、PyVista 0.49.0、Matplotlib 3.11.2、pytest 9.1.1；其余版本和工具探测见 [environment_verified.json](audit/environment_verified.json)、`requirements-lock.txt`。

以下入口已实际运行，数据检查亦从 `/tmp` 启动通过：

```bash
cd /home/lzy/projects/new_mouse_vessel_flow
.venv/bin/mouse-flow --help
.venv/bin/mouse-flow list
.venv/bin/mouse-flow inspect --case roi_only
.venv/bin/mouse-flow analyze --case roi_only
.venv/bin/python -m pytest tests -q
```

以下是本轮执行使用的输出位置；**目录现已存在，复验时改用不存在的新目录**，程序有覆盖保护：

```bash
.venv/bin/mouse-flow validate --case all --full-integrity --output validation/delivery
.venv/bin/mouse-flow visualize --output results/refactor_review
.venv/bin/mouse-flow prepare-cfd --case roi_only --output validation/prepared_case_smoke
```

CLI 日志保存在 `docs/audit/cli_*.json`、`validate_delivery.log`、`pytest_delivery.log` 和 `prepare_case_smoke.log`。原 `presentation` 三个模块的 `--help` 与导入通过，完整渲染未执行。

后续 `compute-particles --case ... --ids ... --output 新目录` 明确建立新任务，使用原积分实现和新源码身份；`compute-cfd --directory ... --build-json ...` 明确调用原生 runner。它们不是查看/验证命令的副作用。本轮没有执行这些正式计算入口，不能把接口存在当成完整批次验证。

新机器重建 Python 环境：用 Python 3.13.11 创建新的 venv，安装 `requirements-lock.txt`，再 `pip install -e .`。本工程采用带本地数据的 editable 安装，普通 wheel 不含数据；复制工程时必须带上所需输入、结果和身份清单。

## 9. 独立性与外部依赖的明确回答

| 问题 | 结论 |
|---|---|
| 当前运行还需要旧项目 Python 源码吗？ | 不需要。当前运行源码仅在 `mouse_vessel_flow/`；迁移审计脚本保留旧路径用于来源追溯 |
| 已迁入工作流还需要旧目录数据吗？ | 不需要。核验和运行入口使用本工程副本；原生环境未安装与数据缺失是不同问题 |
| 是否存在软/硬链接依赖？ | 没有以链接引用旧工程；正式文件按独立文件复制与哈希验证 |
| 当前包是否含旧工作区硬编码路径、sys.path/runpy/exec 加载？ | 静态检查无此类旧加载；历史日志、原始契约和审计来源中的路径原文保留 |
| archive 是否参与运行？ | 不参与。源码身份原始字节证据在显式 `inputs/reference/source_identity/` 中，只校验、不执行 |
| 哪些已能独立运行？ | 配置与场读取、网络/压力回归、原粒子轻量动力学/源和准入、CPU C++ 内核、轨迹统计、WSS 恢复、验证、报告图生成和 CFD 输入准备 |
| 哪些尚未完整验证？ | 新整批微泡生产/恢复、完整几何→网格→原生 GPU CFD、完整动画重渲染；不得用只读检查代替它们 |
| 哪些工具需另外安装？ | 原 svMultiPhysics/PETSc/MPI/CUDA、SimVascular TetGen、Ultraliser/VMTK；独立研究 bridge 可选 LAMMPS |

粒子动力学实际在 CPU；C++ 内核是 CPU 几何查询加速。历史 GPU 场诊断/EGL 绘图与 PETSc GPU CFD 是不同执行路径。本地 CUDA 或 `nvidia-smi` 存在不表示这些路径已经配置完毕。

原 svMultiPhysics upstream commit 为 `c3f0bb892b765b718f61069ecd9726dbc6d177fd`，有效源、真实 diff 和未跟踪但必需的 `preconditioner_reuse.h` 均保留。原 solver 与 PETSc 库 SHA 门禁未绕过。构建方法、环境记录及其未执行边界见 [external/README.md](../external/README.md)。

原 SimVascular 使用 Python 3.5 嵌入式解释器，不能把其 `sv` ABI 直接安装到当前 Python 3.13。新入口调用原生可执行文件执行本工程内兼容脚本；目前只做原语句 AST/语法检查，未借用旧安装后宣称独立网格生成通过。

## 10. 原始数据保护与旧目录处置

[data_manifest.json](audit/data_manifest.json) 给出每一项原路径、新路径、大小与 SHA；[original_files_after.json](audit/original_files_after.json) 证明盘点原文件无变化。正式流场、轨迹、源队列、COMPLETE 记录、原生 checkpoint、物理验收、既有 WSS 与真实日志均按原内容保留；新的报告单独输出，没有改写原报告中的历史路径。

未迁入旧环境、编译缓存、大量可再生中间帧、未执行细档大网格、无内容的上游测试 LFS 资源和无关工程。原发布排除清单、[本次选择清单](audit/selection_decisions.json) 与 [排除审查](audit/excluded_review.json) 保留依据；这些排除不等于授权删除原件。

具体待清理路径、已迁入范围、共享消费者风险和必须保留项见 [CLEANUP_REVIEW_ZH.md](CLEANUP_REVIEW_ZH.md)。旧 `_shared`、Git worktree、共享 native 安装、所有兼容链接与原科研数据本轮均未删除。后续清理需用户确认该清单后另行执行。
