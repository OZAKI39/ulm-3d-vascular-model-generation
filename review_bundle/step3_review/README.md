# HemoCell Step 2–3 外部审查资料

此 bundle 固定实际执行过的测试源码、冻结合同、日志、关键诊断和来源证据，供外部代码审查。同步状态不代表科学结果升级。

- Step1：整理并冻结真实血管 STL、米制单位、1 inlet/3 outlets 的结构身份、法向/面积、人工延长段及来源。生理身份仍为 ASSUMED。
- Step2：闭合体素化与四端口映射；已人工 ParaView 通过。原格点 493×497×280，有效 dx=1.9989918081065344e-7 m，仅 CANDIDATE_BASELINE。
- Step3：原生 Guo 曲面 no-slip 壁、积分归一化入口体积流量、三个独立表压出口，MPI1 从静止启动至 1000 步，含 0/1/10/100/1000 安全 gate。

Step3 当前原始状态：AUTO_PASS_HUMAN_PENDING。首次运行有 MPI 流量监测器代理缓存错误，原记录单独保存且流量无效；修正监测器后按相同参数重跑。两次各 1000 步，总计 2000；两次物理场哈希一致。当前入口误差约 12.6%，通过事先声明的 25% 短程容差；closure≈1.506，两个正表压出口回流，尚未证明稳态质量闭合。

本次只归档与 Git review 提交，没有修改测试源码、原始报告/CSV/JSON/日志，没有运行任何求解、编译、后处理重生成或旧 CFD。原报告内绝对路径及历史 SHA256SUMS 保留；这些是原工作站证据清单，不是本 bundle 的相对路径校验清单。

## 先读这些文件

1. [Step3 报告](step3/STEP3_REPORT.md)
2. [物理 BC 合同](step3/contracts/physical_bc_contract.json)
3. [BC API 可行性](step3/contracts/bc_api_feasibility.md)
4. [单位换算](step3/contracts/lattice_unit_contract.json) 与 [执行合同](step3/contracts/execution_contract.json)
5. [真实流量历史](step3/diagnostics/flow_history.csv) 与 [自动检查](step3/diagnostics/automatic_acceptance.json)
6. [实际 Step3 源码](../../test_code/step3_vascular_pure_fluid/vascularPureFluid.cpp) 和 [源码说明](../../test_code/step3_vascular_pure_fluid/README.md)
7. [Step2 几何报告](step2/STEP2_REPORT.md) 和 [源码](../../test_code/step2_vascular_voxelization_smoke)
8. [文件清单](review_inventory.md)、[review manifest](manifests/review_manifest.json)、[未上传文件](manifests/omitted_artifacts.json)

## 实际运行代码与日志

两个 test_code 目录的源码保持原样。Step3 当前 C++ 对应 attempt2_pre_execution_freeze 的哈希；首次含错误监测器的源码在 step3/provenance/attempt1_source。Python 文件是实际准备/导出/报告工具，此次只提交其现有内容。

Step3 的 0/1/10/100/1000 阶段来自同一连续运行，因此共同使用 step3/logs/mpi1_run.log，没有虚构独立阶段日志。首次无效监测器日志在 step3/attempts/attempt1_invalid_flux。MPI2 未运行。所有本次选择的运行日志都小于 5 MiB，完整复制。

## 物理与数值来源

物理值来自当前旧工程配置/代码及已接受 downstream runtime/promotion 输出；原始快照见 step3/provenance/source_records。当前旧工程 HEAD、branch、remote、dirty state 和仅限相关文件的差异见 provenance/。19 个物理来源均核对了当前文件与 Step3 快照哈希；未把旧工程改动提交到旧仓库。

dt/nuLU/tau 根据当前 HemoCell/Palabos 公式及 Step2 有效 dx 推导，证据见 lattice_unit_contract；未继承旧 Musubi 数值压力偏置。入口 profile 是明确数值假设，不是实验验证分布。

## 未上传内容与适用范围

冻结 STL 为 3,363,184 bytes，低于 10 MiB，已包含。Step2 小型 VTI/VTP 和 cap ID 数组，以及 Step3 测量面 VTP 对几何/积分审查必要且均小于 5 MiB，作为明确例外包含。完整流场 VTU 系列、68.6 MB 原生网格数组、原始流体二进制快照、编译可执行文件/静态库与大体积审计明细未上传；路径、大小、SHA256、字段/步数与作用见 omitted_artifacts.json。

step3/output/flow.pvd 是未经修改的原始系列引用，所指 flow_step_*.vtu 均留在本地，不是可在本 bundle 中直接打开的完整序列。Step2 小型几何文件可独立审查。

尚未验证：Step3 人工 ParaView、生理角色/实验参数/profile、绝对生理血压、稳态、网格/时间步/积分收敛、cap-wall 精度、MPI2、生产 CFD 与 RBC readiness。不得把本次归档成功解释为这些项目 PASS。

## 工作区与发布状态

HemoCell 既有 cmake/setup_googletest.cmake 为 PREEXISTING_NON_STEP3_CHANGE，没有纳入提交；core、Palabos、官方案例均未加入。原 Step3 记录的旧 Git 索引元数据刷新例外保留，不改写科学证据；此次旧工程 Git 审计使用隔离索引，前后完整性另核对。

用户随后明确目标仓库为 OZAKI39/ulm-3d-vascular-model-generation，以其 codex/cfd-wall-force-numerics-validated-sync-20260830 分支为 review 基线。使用独立 sparse checkout 创建全新 review 分支，只新增两份 HemoCell test_code 与此 bundle；不修改本地旧来源工程、不推送已有开发分支、不改两处现有工作区 remote。

HemoCell 科学源码基线仍为 UvaCsl/HemoCell 的 5a410848bd5c57d5ae1c171112e78eab4a82e650。此目标 GitHub 仓库保存审查用测试代码/证据，不是另一次 HemoCell 安装；若审查者后续要构建，应使用该 HemoCell 基线并将 test_code 放入其源码树，依赖与历史命令见原日志。本任务没有构建或执行这些代码。

发布权限的实际检查结果以 provenance/remote_access_check.json 为准。本次不建 GitHub 仓库、不自动 fork、不建 PR。

本 bundle 的局部 .gitattributes 保留原始证据字节/换行与历史空白，并将已选择的小型 STL/VTP/VTI/NPZ 存为普通 Git blob，避免继承目标仓库全局 STL/VTP 的 LFS 转换。未修改既有 .gitattributes 或启用新的 LFS。两份实现源码仍接受常规 whitespace check。
