# 最新归档：修复版 Mirheo 的本地验证、云端部署与完整 RBC 尝试

本分支 `sync/hemocell-mirheo-cloud-rbc-20260911T073240Z` 从 `sync/hemocell-mirheo-single-rbc-repair-20260910T141046Z`（`26c22ed53df1d867aa39443fd65173b0770409c0`）追加归档。以下是已有实验的记录，不是本次同步新跑的实验。旧血管项目根 README 保持原样。

**最新云端任务没有完成 Γ=4。** 准备阶段发生粗碰撞候选溢出，正式剪切为 0 步、Γ=0。`CLOUD_RBC_RUN_COMPLETE=NOT_COMPLETE`、`CLOUD_RBC_NUMERICAL_SCREEN=FAILED`、`RESULTS_RETURN_VERIFIED=PASS`。这些状态逐字节保留，`qualified_speedup=null`。

- [最终中文报告](cloud_results/cloud-rbc-full-20260911T000225Z-review/report_zh.md) · [离线中文 HTML](cloud_results/cloud-rbc-full-20260911T000225Z-review/rbc_full_review.html) · [最终执行交付记录](cloud_compute/rbc_full/execution_delivery_receipt.json)
- [实际原始数据目录](cloud_results/cloud-rbc-full-20260911T000225Z/simulation) · [实际参数](cloud_results/cloud-rbc-full-20260911T000225Z/actual_parameters.json) · [冻结 spec](cloud_results/cloud-rbc-full-20260911T000225Z/full_spec.json) · [来源链](cloud_results/cloud-rbc-full-20260911T000225Z/provenance.json)
- [当前调度及分析代码](cloud_compute) · [实际执行的冻结代码](cloud_results/cloud-rbc-full-20260911T000225Z/tools) · [云端部署报告](cloud_compute/report_zh.md) · [部署摘要](cloud_compute/deployment_summary.json)
- [本地 A0–A6 修复报告](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/report_zh.md) · [本地完整执行回执](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/delivery_receipt.json) · [本地原始运行数据](mirheo_starter/runs/single_rbc_repair/rbc_repair_20260910T131105Z/gpu) · [本地核查页面](mirheo_starter/test_code/outputs/single_rbc_repair/rbc_repair_20260910T131105Z/comparison_review.html)
- [本次同步说明](sync_reports/20260911T073240Z/README.md) · [文件清单](sync_reports/20260911T073240Z/SYNC_MANIFEST.json) · [路径映射](sync_reports/20260911T073240Z/PATH_MAPPING.json) · [排除说明](sync_reports/20260911T073240Z/EXCLUDED_FILES.json)

本轮新增证据覆盖本地 A0–A6 全部成功/失败输出、隔离构建记录、云端失败构建与成功 sm_120 构建、smoke、唯一正式完整验证及本地 CPU 事后分析。核心 CSV/JSON、全部相关 HDF5/NPZ 与错误快照直接保存，没有删行、降精度或转 LFS。HemoCell 旧结果从基础分支继承，本次云端只验证 Mirheo DPD，没有 HemoCell 新运行或速度排名。

| 结论 / 限制 | 具体原始证据与字段 |
|---|---|
| 实际使用修复库与最新 worker | [compile_identity](cloud_results/cloud-rbc-full-20260911T000225Z/provenance.json)：`identity.runtime_binary_hash`、`identity.original_source_sha256`；[实际 rank 0](cloud_results/cloud-rbc-full-20260911T000225Z/simulation/rank_0.json)、[rank 1](cloud_results/cloud-rbc-full-20260911T000225Z/simulation/rank_1.json)；[spec](cloud_results/cloud-rbc-full-20260911T000225Z/full_spec.json)：顶层 `bouncer_policy=shared`、`continuous_observation=true`、`dt=0.0005` |
| 准备调用失败，未进入剪切 | [phase_relaxation](cloud_results/cloud-rbc-full-20260911T000225Z/simulation/phase_relaxation.json)：`native_returned=false`、`successful_returned_steps=null`；[失败摘要](cloud_results/cloud-rbc-full-20260911T000225Z-review/failure_summary.json)：`last_saved_progress_lower_bound.relaxation=45000`、`returned_shear_steps=0` |
| 粗候选 18074 超过 6400 | [原始碰撞计数](cloud_results/cloud-rbc-full-20260911T000225Z/simulation/bounce_counts_relaxation.csv)：step 45040、outer/local、`coarse_count`、`coarse_capacity`、`fine_count=-1`；[console.log](cloud_results/cloud-rbc-full-20260911T000225Z/simulation/console.log)；细筛未执行不能解释成零碰撞 |
| 46 个定期帧与失败瞬间有不同状态 | [geometry.csv](cloud_results/cloud-rbc-full-20260911T000225Z/geometry.csv)：定期帧 A/V 漂移上限约 1.9626% / 0.7920%；[事后诊断](cloud_results/cloud-rbc-full-20260911T000225Z-review/postmortem.json)：`failure_snapshots[0].membrane.z_wall_violations=7`、`geometry_status=INVALID_PERIODIC_GEOMETRY`；[失败原始膜快照](cloud_results/cloud-rbc-full-20260911T000225Z/simulation/bounce_relaxation_outer_local_45040_membrane.csv)。失败瞬间自交与 A/V 为 null，不拿定期帧零自交替代 |
| 成员只检查初态 192 个探针 | [numerical_checks](cloud_results/cloud-rbc-full-20260911T000225Z-review/numerical_checks.json)：`membership.tested_point_frames=192`、`confirmed_mismatches=0`、`strict_impermeability=NOT_VERIFIED`；故障膜无有效周期展开，未进行强行分类 |
| 速度和流体计数 / 动能统计 | [inner 原始 Stats](cloud_results/cloud-rbc-full-20260911T000225Z/simulation/native_stats_relaxation_inner.csv)、[fluid Stats](cloud_results/cloud-rbc-full-20260911T000225Z/simulation/native_stats_relaxation_fluid.csv)；[定期局部场](cloud_results/cloud-rbc-full-20260911T000225Z/simulation/native_relaxation)。原生 kBT 含整体流动动能，不是扣除流速的温度，也不是 K |
| 实际退出及预算 | [execution](cloud_results/cloud-rbc-full-20260911T000225Z/execution.json)：退出码 255、`solver_process_wall_s=36.9039248029876`；[budget_usage](cloud_results/cloud-rbc-full-20260911T000225Z/budget_usage.json)：一次、1800 秒上限、未复用旧预算、未重试。进程墙钟不是 Vast 账单 |
| MPI rank CPU/RSS 缺测 | [原资源采样](cloud_results/cloud-rbc-full-20260911T000225Z/resources.jsonl)、[事后诊断](cloud_results/cloud-rbc-full-20260911T000225Z-review/postmortem.json)：`resource_sampling_coverage`。只采到了两个管理进程，不能用其 RSS 代表全任务；GPU 是全设备采样 |
| 返回完整、原始文件未被报告修复覆盖 | [云端清单](cloud_results/cloud-rbc-full-20260911T000225Z/results_manifest.json)、[本地回传校验](cloud_results/cloud-rbc-full-20260911T000225Z/LOCAL_ARCHIVE_VERIFIED.json)：970 文件 / 65021041 字节；[本地派生清单](cloud_results/cloud-rbc-full-20260911T000225Z-review/review_manifest.json)。原始封包中的回传 PENDING 是封包前状态，本地最终回执为 PASS |
| 浏览器检查已发生，人工仍待验收 | [实际浏览器记录](cloud_results/cloud-rbc-full-20260911T000225Z-review/browser_final/browser_check.json)：`fixture_only=false`、`status=PASS`、`html_sha256`；[截图](cloud_results/cloud-rbc-full-20260911T000225Z-review/browser_final/review.png)。本次同步仅校对哈希，没有重新打开浏览器 |
| 本地修复仍未解决共同终点和质量 | [本地交付](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/delivery_receipt.json)：`GATE_A_UNRESOLVED_HALF_DT_PREPARATION_COMPLETED_QUALITY_FAILED`、`same_quality_same_endpoint_complete=false`；[A6 准备质量](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/A6_preparation_quality.json)；[各 A0–A6 账本](mirheo_starter/runs/single_rbc_repair/rbc_repair_20260910T131105Z/gpu/budget_ledger.json) |
| 物理验证保持限制 | [最终回执](cloud_compute/rbc_full/execution_delivery_receipt.json)：材料 NOT_MATCHED；材料标定、空间/时间收敛与 HemoCell 比较 NOT_TESTED；`qualified_speedup=null` |

库身份：`mirheo-sm120-20731713865ae510`，SHA-256 `d45b4fd1b4498b6f365a012e758eec35bc58a7844dabc66f22e617a7b8cf84ee`；RTX 5090 / CUDA 12.8 / sm_120 / 原单精度；MPI 两个 rank、一个计算 GPU。原生修复为 `local_before_halo_v2`，当前每阶段一次连续 `u.run`。完整配置及旧 HemoCell 硬件/精度记录保留在对应来源文件，未重新测量。

GitHub 的 HTML 链接用于查看/下载文件，不是已发布网站。下载本分支后，用本地浏览器打开离线页面；它自带 Plotly 与真实膜数据，无外部资源。原始 JSON/HTML 的 WSL 与云端绝对来源路径保留，仓库相对入口及 PATH_MAPPING 提供远程核查位置。运行/部署入口仍含原机器路径和既有单次授权，归档不授予再次运行权限。

完整第三方源码树、`.venv`、编译对象/库、安装缓存、锁文件和 `.env` 未新增上传；[原生源码说明](sync_reports/20260911T073240Z/NATIVE_SOURCE_ARCHIVE.md)保留确切提交、许可证、补丁及未跟踪头文件。全部本次选定的科学观测已归档；无法补齐原实验未测的 rank 资源、连续不可渗透证明和失败瞬间有效 A/V。

归档自身检查（不运行求解器）：

```bash
python3 -B sync_reports/20260911T073240Z/verify_archive.py
```

历史同步清单对应其各自提交中的文件版本，不能拿旧清单校验本分支后来更新的同名报告。下方旧入口作为归档历史原样保留，其阶段性结论不能覆盖上方最新交付。

---

# 最新归档：单红细胞剪切流 CPU 修复与证据

本分支 `sync/hemocell-mirheo-single-rbc-repair-20260910T141046Z` 基于最近的单红细胞归档 `sync/hemocell-mirheo-single-rbc-benchmark-20260910T120318Z`，准确基础提交 `7a0be010fa4cfeee4073bba5bab929216cf4df7b`。本次原样归档 `rbc_repair_20260910T131105Z` 的修复代码、CPU 重分析、失败历史、冻结计划和离线核查页；旧基准的数据与说明继续继承。

比较对象为 HemoCell LBM + 可变形细胞、Mirheo DPD + WLC/Kantor 膜。共同域 24³，剪切率 0.02，目标 Γ=4。**本轮修补版未编译、未运行；材料 NOT_MATCHED，根因 SUPPORTED_NOT_CONFIRMED，qualified_speedup=null，人工 PENDING。** 原有 16 项 CPU 测试与 21 项浏览器检查通过，不代表原生修复验证通过。

- [最终中文报告](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/report_zh.md) · [修复记录](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/fix_log.md) · [结果 JSON](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/comparison_results.json) · [最终交付回执](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/delivery_receipt.json)
- [离线 HTML 原件](mirheo_starter/test_code/outputs/single_rbc_repair/rbc_repair_20260910T131105Z/comparison_review.html) · [浏览器记录与截图](mirheo_starter/test_code/outputs/single_rbc_repair/rbc_repair_20260910T131105Z/browser_cpu_delivery/browser_check.json)
- [新版配置](mirheo_starter/py_scripts/single_rbc_benchmark_repaired.yaml) · [修复入口](mirheo_starter/py_scripts/repair_single_rbc_benchmark.py) · [单位与几何代码](mirheo_starter/py_scripts/single_rbc_benchmark/physics.py) · [冻结短对照计划](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/frozen_benchmark_plan.json)
- [新派生数据、原生补丁与隔离 HemoCell 案例](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z) · [继承的全部原始运行数据](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver) · [本次同步清单和校验](sync_reports/20260910T141046Z/README.md)

GitHub 的 HTML 链接用于读取源码。下载或检出本分支后，保留 `mirheo_starter` 的目录结构，用浏览器打开 `comparison_review.html`。脚本与绘图数据内嵌，9 个新证据链接指向同一归档的数据目录；下方动画与耗时明确属于旧 campaign。旧页面还带有一个浏览器记录相对链接，原修复输出目录缺少该文件；归档在该相对位置补入[旧浏览器记录的原样副本](mirheo_starter/test_code/outputs/single_rbc_repair/rbc_repair_20260910T131105Z/browser_delivery_verified/browser_check.json)，它仍属于旧 campaign，新页面的检查记录在 `browser_cpu_delivery`。HTML 原字节未变。没有新修复轨迹。没有启用 Pages 或 LFS。

本轮新 GPU 求解、CPU 求解与编译用量均为 0。新 GPU 8500 秒、CPU 求解 600 秒、隔离编译 1800 秒仍为未批准申请；旧账本余额保持原用途。旧 HemoCell 两次完成 Γ=4，完整任务为 58.17194、60.78252 秒；旧 Mirheo 四次失败成本与最后存帧见原报告，均不能替代合格解耗时。

硬件与精度依据原环境记录：i7-13700HX、WSL 可见 24 个逻辑 CPU；原 HemoCell 为 CPU 2 MPI rank / double，Mirheo 为 RTX 4060 Laptop 单 GPU、1 compute + 1 postprocess rank / single。新候选尚无原生测量。见[本轮环境快照](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/environment.json)与[原任务启动环境](mirheo_starter/data/single_rbc_benchmark/rbc_shear_20260910/environment_authorized_start_corrected.json)。

| 结论或核查项 | 精确证据 |
|---|---|
| 当前关口与未完成比较 | [结果](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/comparison_results.json)：`status`、`runtime_fix_status`、`material_match`、`qualified_speedup`；[关口](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/candidate_comparability.json)：`gate_A/B/C` |
| 溢出是直接退出原因；准确失败步缺测 | [时间线](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/failure_timeline.json)：`runs[].error.coarse_count/capacity/exact_native_step/failing_step_interval`；[主运行原始日志](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/main_mirheo_1/shear_00000.log) |
| 准备阶段成员不符与长诊断自交 | [主运行逐帧表](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/main_mirheo_1_frame_metrics.csv)、[独立诊断时间线](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/diagnostic_mirheo_full_timeline.json)：`first_confirmed_nonadjacent_intersection`；[原始探针](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/main_mirheo_1/fluid_probes_relaxation_00002500.npz)、[原始膜顶点](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/diagnostic_mirheo_full/vertices.csv) |
| 调度根因仍需运行对照 | [源码证据](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/root_cause_evidence.json)：`hypotheses`、`source_sha256`；[补丁](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/native/local_before_halo_v2.patch)、[完整新头文件](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/native/source/src/mirheo/core/bouncers/repair_trace.h)、[来源与重建说明](sync_reports/20260910T141046Z/NATIVE_SOURCE_ARCHIVE.md) |
| 初始弯曲力审计约 0.14% 相对差 | [CPU 能量核查](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/initial_bending_force_audit.json)：`comparisons[0].relative_force_l2_error`；[逐顶点力 CSV](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/initial_bending_cpu_vs_native.csv) |
| 材料和液体可比性尚未通过 | [材料记录](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/material_matching.json)：`status`、`old_failures`、`new_material_runs`；[隔离材料探针源码](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/hemocell_case/benchmark.cpp) |
| 没有新成本测量或新算力授权 | [结果](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/comparison_results.json)：`new_solver_s/new_compile_s/new_end_to_end_s`；[申请](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/authorization_request.json)：`approved=false`；[旧账本](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/budget_ledger.json) |
| CPU / 浏览器 / 人工状态分开 | [CPU 记录](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/cpu_tests.json)、[原浏览器记录](mirheo_starter/test_code/outputs/single_rbc_repair/rbc_repair_20260910T131105Z/browser_cpu_delivery/browser_check.json)、[交付回执](mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/delivery_receipt.json)：`human_review=PENDING` |
| 原路径、源文件身份和归档完整性 | [路径映射](sync_reports/20260910T141046Z/PATH_MAPPING.json)、[文件清单](sync_reports/20260910T141046Z/SYNC_MANIFEST.json)、[归档校验](sync_reports/20260910T141046Z/SYNC_VALIDATION.json)、[机器可读结论索引](sync_reports/20260910T141046Z/RESULT_INDEX.json) |

完整第三方 checkout（含嵌套 `.git`）、环境、编译库和缓存不新增上传；精确版本、许可证、兼容补丁、修复补丁、所有新头文件与被引用源码均保留。没有省略本轮已存在的数值观测，没有脱敏。旧记录原本缺少准确崩溃步与完整候选状态，归档无法补造；在新机器执行仍需重建原生环境。空的新 runs 目录没有文件，因此 Git 不创建占位结果。

原报告中的“本轮未提交或推送”、结果里的 `git_commit_or_push=false` 是本次同步前的原始交付快照，保持原字节。本次 Git 归档状态见独立同步记录；这不改变实验授权或科学结论。

---

## 历史单红细胞与纯流体归档说明（原文保留）

# 最新归档：单个可变形红细胞简单剪切

本分支 `sync/hemocell-mirheo-single-rbc-benchmark-20260910T120318Z` 从上一轮已发布的纯流体归档分支 `sync/hemocell-mirheo-benchmark-20260909T220959Z`（`bcac1d7cf8b6c596f663d047523febee82082b0a`）派生。本次增加单细胞基准的实际代码和全部成功/失败证据，保留下方原纯流体归档说明。没有重跑求解、编译或科学测试。

本轮是 **HemoCell 原生 LBM + High Order RBC + IBM** 与 **Mirheo DPD + WLC/Kantor 膜 + 双向作用**。共同有效域 24³、X 流向/Z 梯度、同一参考细胞和 Γ=4 目标。HemoCell 为 CPU 2 MPI rank / double；Mirheo 为 RTX 4060 Laptop 单 GPU、1 compute + 1 postprocess rank / single。

原报告结论：**工作流 PARTIAL、模型可比性 PARTIAL、BENCHMARK_SCREEN FAILED、qualified_speedup=null**；研究适用性 NOT_VALIDATED，人工 PENDING。HemoCell 两次完整成本为 58.172、60.783 秒（含必要分析）；Mirheo 两次主运行只存到 Γ=0.6、2.9 后失败，半步到 Γ=0.7，受限队列长诊断到 Γ=3.4。失败进程秒不当作完成 Γ=4 的端到端秒。

- [最终中文说明及实际命令](mirheo_starter/test_code/README_single_rbc_benchmark.md)
- [原始离线 HTML](mirheo_starter/test_code/outputs/single_rbc_benchmark/rbc_shear_20260910/single_rbc_review.html) · [原始结果 JSON](mirheo_starter/test_code/outputs/single_rbc_benchmark/rbc_shear_20260910/results.json) · [便于网页阅读的结果字段索引](sync_reports/20260910T120318Z/RESULT_INDEX.json)
- [共同配置](mirheo_starter/py_scripts/single_rbc_benchmark.yaml) · [单位、几何与指标代码](mirheo_starter/py_scripts/single_rbc_benchmark/physics.py) · [新增 HemoCell 案例](hemocell_starter/cases/single_rbc_shear_benchmark/benchmark.cpp)
- [全部原始运行目录](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910) · [计划、授权及数值证据](mirheo_starter/data/single_rbc_benchmark/rbc_shear_20260910) · [本次归档校验与排除说明](sync_reports/20260910T120318Z/README.md)

GitHub 的 HTML 链接是源文件入口。请下载本分支，保留对应输出目录后在本地浏览器打开 HTML；Plotly 和数据已内嵌，浏览器证据链接使用同目录文件。原结果 JSON 为 57,568,143 字节，已单独审查并原样保留；若网页无法预览，请下载或使用 Git 读取。轻量字段索引不替代完整数据。

| 原报告结论或核查项 | 精确证据入口与字段 |
|---|---|
| 双方没有完成共同终点比较 | [结果字段索引](sync_reports/20260910T120318Z/RESULT_INDEX.json)：`formal_cold_completed=2`、`workflow=PARTIAL`、`qualified_speedup=null` |
| HemoCell 两次完整成本及散布 | [完整结果](mirheo_starter/test_code/outputs/single_rbc_benchmark/rbc_shear_20260910/results.json)：`repeat_statistics.hemocell`；[第 1 次原生计时](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/main_hemocell_1/timings.csv)、[第 2 次原生计时](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/main_hemocell_2/timings.csv) |
| D、倾角、面积、体积来自实际膜 | [HemoCell 原始顶点](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/main_hemocell_1/vertices.csv)、[Mirheo 原始顶点](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/main_mirheo_2/vertices.csv)；`runs[].frames[].metrics` |
| DPD 黏度半步未通过原门槛 | [材料测量](mirheo_starter/data/single_rbc_benchmark/rbc_shear_20260910/fluid_measured.json)：`primary.nu`、`strict.nu`、`strict_difference`；[原始剖面](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/material_dpd/profiles.csv)、[半步剖面](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/material_dpd_strict/profiles.csv) |
| 原生膜响应及准备后形状未充分匹配 | [膜测量](mirheo_starter/data/single_rbc_benchmark/rbc_shear_20260910/membrane_measured.json)：`comparison`、`max_response_relative_difference`；完整结果 `prepared_shapes` |
| 空通道与细胞对流体的实际影响 | [HemoCell 空流剖面](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/empty_hemocell/profiles.csv)、[DPD 空流剖面](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/empty_dpd/profiles.csv)、[有细胞局部流场](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/main_mirheo_2/local_flow.csv)；完整结果 `feedback` 保留不同可用窗口与热噪声限制 |
| Mirheo 主运行和半步碰撞溢出失败 | [停止说明](mirheo_starter/data/single_rbc_benchmark/rbc_shear_20260910/native_failure_final.json)、[第 1 次日志](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/main_mirheo_1/console.log)、[第 2 次日志](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/main_mirheo_2/console.log)、[半步日志](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/strict_mirheo/console.log) |
| 受限队列诊断未修复长任务 | [冻结诊断计划](mirheo_starter/data/single_rbc_benchmark/rbc_shear_20260910/runtime_candidate_frozen.json)、[长诊断日志](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/diagnostic_mirheo_full/console.log)、[原始探针所在目录](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/diagnostic_mirheo_full)；完整结果 `runtime_diagnostics` |
| 26 次求解尝试、34 条执行记录与全部费用 | [求解账本](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/budget_ledger.json)、[构建账本](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/preparation/budget_ledger.json)、[执行顺序核查](mirheo_starter/data/single_rbc_benchmark/rbc_shear_20260910/execution_integrity.json) |
| 原有 CPU / 浏览器检查的身份 | [26 项 CPU 测试](mirheo_starter/data/single_rbc_benchmark/rbc_shear_20260910/cpu_tests_final.json)、[19 项最终浏览器检查](mirheo_starter/test_code/outputs/single_rbc_benchmark/rbc_shear_20260910/browser_delivery_verified/browser_check.json)、[原最终交付记录](mirheo_starter/data/single_rbc_benchmark/rbc_shear_20260910/final_delivery.json)；此前浏览器 FAIL 记录也保留 |
| 精度、库身份、资源及 Palabos 来源更正 | [更正后的启动环境](mirheo_starter/data/single_rbc_benchmark/rbc_shear_20260910/environment_authorized_start_corrected.json)、[现有构建文本](hemocell_starter/metadata/single_rbc_benchmark/rbc_shear_20260910/build_evidence/CMakeCache.txt)、[具体任务资源](mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/solver/diagnostic_mirheo_full/host_resources.json) |
| 原始 WSL 路径如何对应本仓库 | [PATH_MAPPING.json](sync_reports/20260910T120318Z/PATH_MAPPING.json)；原配置、provenance 和 HTML 均未改写 |

本轮数值原始数据、NPZ 探针、分析缓存、冻结 worker/案例源码以及失败历史均已保留。排除临时锁、编译程序和对象、完整运行环境与第三方源码树；准确版本、已有补丁及许可证由上轮归档继承，新增构建文本另存。已有数据的主要结论可以离线核查；在另一台机器重新执行原生求解仍需要相应环境。没有对原文件脱敏或改变科学状态。

---

## 上轮纯流体归档说明（原文保留）

# HemoCell / Mirheo 基准远程核查入口

本分支仅归档已经完成的无壁周期双向 Poiseuille 纯流体基准，供远程读取代码和实测数据。
本次同步没有重新运行仿真、编译、标定、物理测试或浏览器测试，没有重新判断哪套方法更好。
基础分支为 `sync/mirheo-starter-20260908`，base commit 为 `48d53e9c86963a0d3520b8d262b79ec21fe76980`；本分支为 `sync/hemocell-mirheo-benchmark-20260909T220959Z`。

## 直接查看

- [最终中文报告（保持原始字节）](mirheo_starter/test_code/outputs/solver_benchmark/review_7b2382c0e17324e6/结果说明.md)
- [最终离线 HTML（保持原始字节）](mirheo_starter/test_code/outputs/solver_benchmark/review_7b2382c0e17324e6/benchmark_review.html)
- [完整 results.json](mirheo_starter/test_code/outputs/solver_benchmark/review_7b2382c0e17324e6/results.json) 与 [原交付记录](mirheo_starter/test_code/outputs/solver_benchmark/review_7b2382c0e17324e6/delivery_record.json)
- [共同配置](mirheo_starter/py_scripts/solver_benchmark.yaml) 与 [冻结物理定义、任务和精确单位映射](mirheo_starter/data/solver_benchmark/solver_benchmark_20260909/frozen_plan.json)
- [原始 CPU 运行目录](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu) 与 [原始 GPU 运行目录](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/gpu)
- [HemoCell C++ 案例](hemocell_starter/cases/pure_fluid_benchmark/benchmark.cpp)、[CMake](hemocell_starter/cases/pure_fluid_benchmark/CMakeLists.txt)、[安装脚本](hemocell_starter/scripts/setup_hemocell.sh)
- [共同 Python 调度入口](mirheo_starter/py_scripts/benchmark_mirheo_hemocell.py)、[分析实现](mirheo_starter/py_scripts/solver_benchmark/analysis.py)、[原操作说明](mirheo_starter/test_code/README_solver_benchmark.md)
- [本次同步清单与校验说明](sync_reports/20260909T220959Z/README.md)；[原路径到仓库路径映射](sync_reports/20260909T220959Z/PATH_MAPPING.json)

GitHub 中的 HTML 链接通常显示源码。下载该 HTML 后，用 Edge/Chrome 直接打开即可；页面内嵌数据与 Plotly JavaScript，无外部资源依赖，无需 Web 服务。
读取原始 CSV/JSON 可通过克隆本分支或下载所需目录完成；本次未启用 GitHub Pages、Git LFS 或部署服务。
原始配置保留 WSL 绝对路径，查看入口的链接使用仓库相对路径。远程读取不要求重跑 GPU；原调度器不是经过路径迁移的新执行环境。

## 原报告所记录的范围和状态

比较双方为 **HemoCell/Palabos 原生 D3Q19 Guo BGK LBM 纯流体路径**与 **Mirheo SDPD（Wendland C2、Linear EOS）**。
共同盒为 4 µm 周期立方盒，物性和体力见冻结配置；本轮未运行 DPD 新基准，也未运行圆管、RBC、真实血管或压力出口实验。
HemoCell 使用 i7-13700HX / WSL2 CPU、双精度、1/2/4 MPI rank，OMP_NUM_THREADS=1；WSL 可见 24 个逻辑 CPU。
Mirheo 使用 RTX 4060 Laptop GPU、现有单精度构建，1 个计算 rank 加 1 个后处理 rank。

原结果总体状态为 **`PARTIAL`**，人工验收为 **`PENDING`**。
HemoCell 的 8 次主基准通过本轮 `PROPOSED` 局部平均流场筛选（含细网格核查）。
SDPD 完成一次成功重试，但流场误差、温度、漂移及独立块数未满足原门槛；它的合格解耗时和全部合格速度比均为 `null`。
本分支保留这一结论，不用执行成本重新计算“赢家”。圆管保持 `PIPE_BENCHMARK_DEFERRED`。

共有 19 次 CPU 数值运行，计费 119.365107859951 秒；GPU 首次启动失败（零流体步）加一次成功重试，合计 320.162673518993 秒。
下表直接列出原 results.json 的记录；资源审计、低 I/O 和冒烟中的 `qualified=false` 不代表重新判定科学失败，而是原记录不用于合格解比较。

| 任务与冻结参数 | 原 role | 原 status | 实际步数 | 原 qualified |
|---|---|---|---:|---|
| [deployment_smoke](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/deployment_smoke/task.json) | smoke | COMPLETED | 64 | False |
| [lbm_N16_r1_main_rep1](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N16_r1_main_rep1/task.json) | main | COMPLETED | 2511 | True |
| [lbm_N16_r1_main_rep2](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N16_r1_main_rep2/task.json) | main | COMPLETED | 2511 | True |
| [lbm_N32_r1_main_rep1](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N32_r1_main_rep1/task.json) | main | COMPLETED | 10045 | True |
| [lbm_N32_r1_main_rep2](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N32_r1_main_rep2/task.json) | main | COMPLETED | 10045 | True |
| [lbm_N16_r1_lowio_rep1](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N16_r1_lowio_rep1/task.json) | low_io_cost | COMPLETED | 12000 | False |
| [lbm_N16_r1_lowio_rep2](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N16_r1_lowio_rep2/task.json) | low_io_cost | COMPLETED | 12000 | False |
| [lbm_N16_r2_main_rep1](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N16_r2_main_rep1/task.json) | main | COMPLETED | 2511 | True |
| [lbm_N16_r2_main_rep2](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N16_r2_main_rep2/task.json) | main | COMPLETED | 2511 | True |
| [lbm_N16_r2_lowio_rep1](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N16_r2_lowio_rep1/task.json) | low_io_cost | COMPLETED | 12000 | False |
| [lbm_N16_r2_lowio_rep2](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N16_r2_lowio_rep2/task.json) | low_io_cost | COMPLETED | 12000 | False |
| [lbm_N16_r4_main_rep1](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N16_r4_main_rep1/task.json) | main | COMPLETED | 2511 | True |
| [lbm_N16_r4_main_rep2](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N16_r4_main_rep2/task.json) | main | COMPLETED | 2511 | True |
| [lbm_N16_r4_lowio_rep1](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N16_r4_lowio_rep1/task.json) | low_io_cost | COMPLETED | 12000 | False |
| [lbm_N16_r4_lowio_rep2](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N16_r4_lowio_rep2/task.json) | low_io_cost | COMPLETED | 12000 | False |
| [lbm_N16_r1_memory_audit](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N16_r1_memory_audit/task.json) | memory_audit | COMPLETED | 2511 | False |
| [lbm_N32_r1_memory_audit](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N32_r1_memory_audit/task.json) | memory_audit | COMPLETED | 10045 | False |
| [lbm_N16_r2_memory_audit](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N16_r2_memory_audit/task.json) | memory_audit | COMPLETED | 2511 | False |
| [lbm_N16_r4_memory_audit](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N16_r4_memory_audit/task.json) | memory_audit | COMPLETED | 2511 | False |
| [sdpd_main](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/gpu/sdpd_main/task.json) | main | FAILED_OR_STOPPED | 0 | False |
| [sdpd_main_attempt_02](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/gpu/sdpd_main_attempt_02/task.json) | main | COMPLETED | 252719 | False |

## 结论—证据索引

下列 results 字段按 `comparison.results` 数组中的 `task_id` 定位；不依赖数组顺序。

| 原记录中的结论或限制 | 具体证据与字段 |
|---|---|
| HemoCell 已安装并推进零细胞原生流体 | [installed.json](hemocell_starter/metadata/installed.json)：`hemocell_commit`、`palabos_commit`、`binary_sha256`；[冒烟 completion](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/deployment_smoke/completion.json)：`cell_count=0`、`actual_steps=64`、`completed`；[execution](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/deployment_smoke/execution.json)：`exit_code` |
| 物性、体力及各后端 dt 来自共同定义 | [冻结计划](mirheo_starter/data/solver_benchmark/solver_benchmark_20260909/frozen_plan.json)：`physics`、`gpu_plan.dt_si`、`cpu_tasks[].dt_si`；[原单位映射](mirheo_starter/data/fluid_physics/result_1b008103675eb6ff/unit_mapping.json)：`L0`、`M0`、`t0`、`required_mu_star`、`required_nu_star` |
| LBM 16³/32³ 测量与细网格核查 | [16³ 原始 profiles.csv](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N16_r1_main_rep1/profiles.csv)、[32³ profiles.csv](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N32_r1_main_rep1/profiles.csv)：`step,time_si,bin,count,ux_si,raw_ux_si,half_force_si`；[结果](mirheo_starter/test_code/outputs/solver_benchmark/review_7b2382c0e17324e6/results.json)：`accuracy.profile_relative_l2`、`accuracy.half_flow_relative_error`、`accuracy.apparent_nu_relative_error`、`comparison.refinement_relative_l2` |
| SDPD 实际完成 252719 步，非计划值替代 | [成功 completion](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/gpu/sdpd_main_attempt_02/completion.json)：`actual_steps`、`actual_time_si`；[实际执行](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/gpu/sdpd_main_attempt_02/execution.json)：`exit_code,elapsed_monotonic_s,timeout`；[当时 worker](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/gpu/sdpd_main_attempt_02/mirheo_worker.py) 与 [native_mpi](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/gpu/sdpd_main_attempt_02/native_mpi.py) |
| SDPD 剖面、通量与拟合依赖完整原始序列 | [SDPD profiles.csv](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/gpu/sdpd_main_attempt_02/profiles.csv)：`time_si,bin,ux_si,count`；[结果](mirheo_starter/test_code/outputs/solver_benchmark/review_7b2382c0e17324e6/results.json)：`accuracy.measured_profile_si`、`accuracy.full_window_descriptive_metrics`、`accuracy.half_flow_si`、`accuracy.apparent_nu_si`；[分析实现](mirheo_starter/py_scripts/solver_benchmark/analysis.py) 与 [拟合和分箱函数](mirheo_starter/py_scripts/solver_benchmark/physics.py) |
| SDPD 温度与独立采样不足，无有效 CI | [moments.csv](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/gpu/sdpd_main_attempt_02/moments.csv)：`time_si,temperature_K`；[结果](mirheo_starter/test_code/outputs/solver_benchmark/review_7b2382c0e17324e6/results.json)：`temperature_statistics`、`accuracy.block_statistics`、`accuracy.profile_ci95_si`、`accuracy.stationarity_relative_drift` |
| 执行成本与物理时间成本，保留全部 rank/重复 | [SDPD 计时块](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/gpu/sdpd_main_attempt_02/timings.csv) 与 [LBM 示例计时块](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N16_r2_main_rep1/timings.csv)：`chunk_steps,compute_max_rank_s,sampling_max_rank_s,time_si`；各任务 `execution.json.elapsed_monotonic_s`；[结果](mirheo_starter/test_code/outputs/solver_benchmark/review_7b2382c0e17324e6/results.json)：`timing.compute_ms_per_step`、`timing.workflow_wall_s_per_us` |
| 首次启动失败计费；只获准重试一次 | [失败 console.log](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/gpu/sdpd_main/console.log)、[失败 execution](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/gpu/sdpd_main/execution.json)；[修正与重试请求](mirheo_starter/data/solver_benchmark/solver_benchmark_20260909/explicit_retry_request.json)、[结构化重试确认](mirheo_starter/data/solver_benchmark/solver_benchmark_20260909/explicit_retry_confirmation.json)；[GPU ledger](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/gpu/budget_ledger.json)：`attempts` |
| CPU/GPU 累计成本与缓存复用 | [CPU ledger](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/budget_ledger.json)、[GPU ledger](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/gpu/budget_ledger.json)：`attempts[].charged_s`；[缓存重执行核查](mirheo_starter/data/solver_benchmark/solver_benchmark_20260909/cache_reexecute_check.json)：`before,after,solver_budget_files_unchanged`；[当前共享池单独快照](mirheo_starter/data/solver_benchmark/solver_benchmark_20260909/archive_evidence/shared_budget_pool_at_sync.json) |
| 原主机采样漏掉 ranks，CPU 后补审计，Mirheo 主机 RSS 仍缺测 | [CPU 4 rank 审计](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/cpu/lbm_N16_r4_memory_audit/host_memory.json)：`ranks,sampled_peak_tree_rss_bytes,minimum_WSL_available_bytes`；[结果](mirheo_starter/test_code/outputs/solver_benchmark/review_7b2382c0e17324e6/results.json)：`memory.host.coverage`、`memory.host.sampled_peak_tree_rss_bytes`；[GPU execution](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/gpu/sdpd_main_attempt_02/execution.json)：`resource_samples,device_sampled_peak_used_MiB` |
| 核对冻结脚本与原始输出身份 | 每次运行的 `output_sha256.json`；[SDPD 示例清单](mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909/gpu/sdpd_main_attempt_02/output_sha256.json)；[原 worker 修正记录](mirheo_starter/data/solver_benchmark/solver_benchmark_20260909/implementation_revision_integrator.json) |
| 原 129 项 CPU 回归与浏览器检查，不是本次重跑 | [原 CPU 回归日志](mirheo_starter/data/solver_benchmark/solver_benchmark_20260909/CPU_delivery_regression.log)；[原浏览器检查](mirheo_starter/test_code/outputs/solver_benchmark/review_7b2382c0e17324e6/browser_local_profile/browser_check.json)：`status,html_sha256,checks`；[截图](mirheo_starter/test_code/outputs/solver_benchmark/review_7b2382c0e17324e6/browser_local_profile/overview.png) |
| 历史数据仅作背景参考 | [结果](mirheo_starter/test_code/outputs/solver_benchmark/review_7b2382c0e17324e6/results.json)：`historical_reference.source_category=HISTORICAL_REFERENCE`、`source_sha256`、`limitations`；[对应历史源文件](mirheo_starter/data/sdpd_diagnostics/result_4ad869c67164f241/corrected_tasks.json)，由基础分支继承 |
| 编译与环境身份，没有为归档重新构建 | [HemoCell 环境](hemocell_starter/metadata/environment_before.json)、[源和补丁](hemocell_starter/metadata/sources.json)、[构建证据](hemocell_starter/metadata/build_evidence)、[实际构建日志](hemocell_starter/logs/build_library.log)；[Mirheo 编译器](mirheo_starter/metadata/compiler.txt)、[CUDA](mirheo_starter/metadata/cuda_toolkit.txt)、[原 Python 清单](mirheo_starter/metadata/python_packages.txt)、[同步时只读依赖补充](mirheo_starter/data/solver_benchmark/solver_benchmark_20260909/archive_evidence/python_packages_at_sync.json) |
| 合格解速度比保持 null，圆管及人工状态未提升 | [结果](mirheo_starter/test_code/outputs/solver_benchmark/review_7b2382c0e17324e6/results.json)：`comparison.qualified_speedups`、`comparison.pipe_benchmark`、`human_review`；[最终中文原报告](mirheo_starter/test_code/outputs/solver_benchmark/review_7b2382c0e17324e6/结果说明.md) |

## 保留的阶段差异与归档限制

原 YAML 的 `policy.GPU_authorization=PENDING_SCOPE_APPROVAL` 是运行前冻结状态；实际授权和用量见结构化授权及账本，最终 results 的预算为已授权状态。
原始主机采样文件中的早期“tree RSS”未覆盖 MPI ranks，最终分析已明确降级为未测；四次 CPU 补充审计单列。此次同步未改写任何上述原件。
SDPD 初始 `moments.csv` 的 `kernel_rho_mean_si` 和 `kernel_rho_cv` 为未测的空字段，按原 CSV 保留；JSON 中未知值保持 `null`，没有补值或删行。

完整 vendor 源码、`.venv`、编译对象、可执行文件、共享库/静态库、下载缓存、临时浏览器 profile、重复 ZIP 和两份独立聊天答复文本不新增上传。
仅保留与预算相关的结构化授权；源码准确提交、官方补丁、许可和实际库哈希均已保留。
本轮用于平均流场、温度、拟合、计时与资源结论的原始序列已逐字节归档，排除项不影响这些主要结果复核。
Mirheo 主机 rank RSS 是原实验缺测；该缺口无法通过归档补齐。被排除的授权原文本和本地 ZIP 无法远程逐字节重算其源哈希，相关记录明确保留。
完整执行环境未打包；这不会阻止无 GPU 的 CSV/JSON/源码审阅。各项路径、大小、哈希及影响见 [EXCLUDED_FILES.json](sync_reports/20260909T220959Z/EXCLUDED_FILES.json)。

HemoCell、Palabos 及派生 C++ 案例为 AGPL-3.0-or-later，Mirheo 为 MIT；参见 [案例 COPYING](hemocell_starter/cases/pure_fluid_benchmark/COPYING) 和 [官方脚本/补丁/许可](hemocell_starter/metadata/upstream_evidence/HemoCell)。
本分支继承的历史文件不清理、不改写；原仓库根 README 和血管源码保持不变。
