# F10K：科学适用范围、持续形态变化与有限规模测试

日期：2026-10-12。基线为本地完整 F10J `91b09993167625fc103ebd6a026773753a12d3f4`，工作分支 `dev/fast-pzt-f10k-science-scope-scaling-20261012`，独立 worktree `mevo_fast_pzt_f10k_20261012`。科研结果仍为 **EXPLORATORY / PHYSICS_UNCALIBRATED / training_ready=false**。

本轮得到一个 B2 持续图聚合的正例和有限规模性能证据；没有得到完整拓扑收敛、指定目标导航成功或物理校准。A/B/C 是主要任务，D 在独立预算内完成，没有扩展为生产求解器。

## 结果与证据

| 层次 | 实际结论 | 机器可读证据（相对 results/fast_pzt/f10k） |
|---|---|---|
| 历史保护 | F10J 的 71 项结果和 13 项原始输入 SHA-256 通过 | BASELINE_AUDIT.json |
| 模型范围 | 保留 F5–F10J 基线证据，另加 F10K 有限新证据 | analysis/SCIENTIFIC_SCOPE_FINAL.json |
| F5 候选 | 417.5 ms 状态几何无可分辨穿壁，但近壁 B2 不准入 | analysis/F5_ADMISSION.json |
| 替代状态 | F10D prepare_1101 第 42 帧、绝对时间 2.1 ms，远壁真实状态 | analysis/ALTERNATIVE_SCREEN.json |
| 状态恢复 | 新缓存 M0 重算原始下一个 50 µs 保存间隔，坐标逐位一致 | topology_dt10_v1/RESTART_CHECK.json |
| 新 B2 对照 | 10/5 µs，各推进 3 ms；同成员持续聚合通过 | analysis/B2_TOPOLOGY_COMPARISON.json |
| 瞬时连边 | 新算例 46/47 条，最大配对时间差 60 µs，FAIL | 同上；历史 12/29 FAIL 未改变 |
| 指定目标 | 实际下游 20 µm 截面在约 2.024 ms 通过，额外 Z 位移约 20.3 nm，导航未建立 | analysis/NAVIGATION.json |
| 固定状态性能 | N=32/128/512 两端收敛；1024 因预设 CPU 数组内存预算 NOT_RUN | analysis/SCALING.json |
| 数据准入 | 不具备 AI 科学数据生产准入 | analysis/FINAL_GATE.json |

## 新 B2 科学实验

F5 保存帧 SHA-256 为 `20d3f3f56a4771a3a3fa8078f83556f2b4674f9e8cef5cc8bfac1e0c76e003c7`。独立原三角面查询最小泡面壁距为 −4.02e−21 m，仅为零点舍入，不能称为实际穿壁。其壁面接触不满足本轮内部管腔 B2 的适用范围，因此明确记录 `B2_TOPOLOGY_CASE_NOT_ADMITTED`；没有平移泡团，也没有拿历史 M0 的 418.54 ms 分裂当作拟合目标。

有界筛选检查 4 个现存替代状态，复用 F10J 的 98 条历史账本清单及 6 个不完整账本排除记录，没有随机搜索。选中 `results/fast_pzt/f10d/campaign_v2_diversity/prepare_1101/trajectory.h5`，SHA-256 `372fcfed6a93719693c69fd891521f7c389b1fd9809db5538af784f91b2de7dd`。第 42 帧为 31+1 连通分量，最小泡面壁距 1288.687 µm，原始声场和原四面体 P1 CFD 查询均有效。原控制为 O、CFD 开启、F2 二次力开启，全部 ID 活跃，没有出口事件。

原坐标、半径、控制和时钟可恢复。M0 新缓存重放 50 µs 后坐标逐位一致，但这只支持可重建初始化，不构成 B2 精确 checkpoint。新 B2 使用局部时间 0 和明确绝对时间偏移 2.1 ms；B2 连通历史重新开始，不继承 M0 已发生的持续性。这是**跨模型初始化的探索性反事实实验**。

采用原 F10H 几何单边弹簧 B2、局部润滑阻力及刚度 0.04641319518378658 N/m；两步长共 900 个接受步。没有放宽原 17.8 nm 重叠、178 nm 单次修正或 1e−13 m 轨迹接触求解容差。先完成 10 µs 窗口并检查质量和持续事件，才执行 5 µs 窗口。

| 指标 | 10 µs | 5 µs |
|---|---:|---:|
| 物理窗口 | 3 ms | 3 ms |
| 接受步数 | 300 | 600 |
| 持续图聚合时刻（局部） | 0.290 ms | 0.290 ms |
| 观察持续时间 | 2.710 ms | 2.710 ms |
| 最大重叠 | 11.419 nm | 11.533 nm |
| 最大单次修正 | 1.938 nm | 0.964 nm |
| 独立实际方程残差 | 9.971e−14 m | 6.571e−14 m |
| 原始连边事件数 | 46 | 47 |
| 积分及质量计算墙钟 | 3.201 s | 4.336 s |
| 含诊断输出的进程墙钟 | 12.219 s | 17.926 s |

ID 5 与其余 31 泡连通后保持至终点；无右删失的持续聚合、无持续分裂。301 个共同真实时刻上成员集合一致；最大质心差 0.039854 nm，固定 ID 位移差 27.569 nm，Kabsch RMSD 12.153 nm，Rg 相对差 0.01555%，长短比 0.11935%，形状张量 0.19329%。这些通过原门限。原始连边仍 FAIL：配对时间差最大 60 µs，5 µs 多出 (12,17) 在 2.945 ms 的断边。总体仍是 `MORPHOLOGY_NOT_CONVERGED`。持续图聚合不等于物理泡融合；缺少持续分裂正例和实验校准。

## 工程修改、实际失败与恢复

新增 F10K 脚本、协议和 opt-in `ScalingB2` 子类。唯一修改历史求解器的位置是 `f10i_gpu.FixedStateB2` 的可覆盖规模上限，默认仍 32；原 `OptimizedB2(128)` 继续拒绝，N=32 三次调用与 F10J 位置、力逐位一致。未修改原默认 M0/B2 方程、接触刚度、PCG/非线性迭代上限或历史结果。

10 µs 科学轨迹完成并落盘后，公共运行器写配置时缺少 `protocol_sha256` 引发 KeyError，原进程 exit 1。修复身份字段后单独执行元数据恢复，逐项验证完成终点、步数、质量和原文件哈希，写 `topology_recovery_v1/RECOVERED_CONFIG.json`；没有重跑轨迹或覆盖失败文件。原 10 µs HDF5 没有新增绝对时刻数组，偏移在恢复元数据中；5 µs 文件直接保存绝对时刻。

规模测试 v1/v2 在导入阶段缺少 threadpoolctl，未执行数值求解。现有 venv 无 pip，保留失败安装记录，通过 uv 仅增加 `threadpoolctl==3.6.0`（官方 wheel 来源记录在安装日志和版本配置），不升级其他包或驱动。v3 完成正式测试。早期导航离线分析因浮点相等比较选行失败，改用整数步长构造，原轨迹未修改。所有失败 stderr、计时和不可变源码均归档，详见 `analysis/COST_AND_FAILURES.json`。

两条轨迹的诊断输出分别约 7.57/12.46 s，明显高于单纯积分，应将全步审计视为本轮研究成本，不能用简化计时声称生产吞吐。六个服务器进程（含失败启动和恢复）累计约 41.247 s；不含本地分析、绘图、传输、安装及研发时间。

## 执行与归档身份

服务器独立目录 `/workspace/mevo_acoustic_microbubbles/f10k`，复用父目录 venv；RTX 4090 / Ryzen 9 7950X3D。实际环境版本见 `ENVIRONMENT.json`。执行前核查 Supervisor 与 GPU，未打断其他科学项目；归档时 F10K 进程均已退出，GPU compute 列表为空。

| 执行 | 不可变 Git 提交 |
|---|---|
| 10 µs 科学窗口（元数据尾部失败） | 49f8109dd5b1ab607133e73ea4c6124f667f4277 |
| 仅元数据恢复 | b3a9d7e271c13cc05c78077906ba0049430af72d |
| 5 µs 窗口 | 1b3514acc3a6aa6aac25d3dcc83451725bdb6ca0 |
| 两个缺包启动失败 | 9ba0b3e97abecfec0e2bc19a055979b91bbb9210 |
| 正式规模测试 v3 | 8e1e9aaec50fd11f5716051e7c301c19da31fdca |

77 项服务器结果逐文件回传验 SHA-256；五个源码快照分别 481/482/483/485/486 项文件与对应 `git archive` 核对，包括启动失败快照。13 项输入重新核验。输入和旧轨迹沿用已验证归档，不复制无关大文件。最终本地 `ARTIFACT_MANIFEST.json` 包含科学日志、图件、分析、文档变更证据及源文件索引。

四组英文 PNG/PDF 位于 `results/fast_pzt/f10k/figures/`。绘图使用原壁面、实际逐泡中心和真实球径；局部坐标仅平移原点，不平移科学初态，不聚类简化或平滑。PDF 保留文字和曲线，三维表面栅格化。实际打开 PNG 及 PDF 渲染页检查，记录 `VISUAL_REVIEW.json` 与 `FIGURE_PROVENANCE.json`。

本地科研进度文档采用新增短段落，备份、before/after、diff、SHA 和标题结构检查见 `progress_document/`。本轮没有 GitHub 推送任务，交付保存在独立本地分支。

## 复核与复现

在本 worktree 项目目录，使用已有共享 Python 环境，`PYTHONPATH=src:. OPENBLAS_NUM_THREADS=1`。离线执行 `python -m scripts.fast_pzt_f10k.verify` 核查，不启动服务器计算。需要重建分析时，依次运行 `analyze_topology`、`navigation`、`analyze_scaling`、`delivery`、`figures`；重绘后必须重新人工打开图件并更新视觉审查哈希。原始执行身份与当前分析身份分开，不把重分析提交当成原科学运行提交。

实际数值运行由 `remote.py` 从干净 Git 提交创建不可变快照并使用独立 Supervisor 标签；原有标签和结果不得覆盖。D 是 B=1 单个接触问题随 N 扩展，不是 B 条轨迹并行。详细范围及下一步协议见同目录 `FAST_PZT_F10K_SCIENTIFIC_SCOPE.md`，性能口径见 `FAST_PZT_F10K_SCALING.md`。
