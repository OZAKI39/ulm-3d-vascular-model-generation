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
