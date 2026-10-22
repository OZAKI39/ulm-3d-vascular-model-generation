# F10D：完整逐泡后处理、多控制家族与独立初态检查

本轮工程后处理与性能验收通过，真实不同初态的三控制覆盖未完成：四个种子中两个准备失败；另两个种子的 O/S 分支被原接触预算停止，只有 OFF 分支完整通过。不能将 `campaign_v2/COMPLETE.json` 的流程 PASS 解读为所有科学轨迹成功。保持 **EXPLORATORY / PHYSICS_UNCALIBRATED，training_ready=false**。

## 基线、环境与实验身份

- 从完整 F10C `2fc90fd4ad57d16addfe5bed2878c7c3085efd61` 建立独立 worktree 和分支 `dev/fast-pzt-f10d-batch-data-pipeline-20261011`，未回退到 F5，未覆盖旧工作区。
- WSL：`/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10d_20261011/computation_examples/mevo_acoustic_microbubbles`。
- 服务器：`/workspace/mevo_acoustic_microbubbles/f10d`；使用原独立 `.venv`，Warp 1.18.0、NumPy 2.2.6、SciPy 1.15.3、h5py 3.16.0、Python 3.12.3；RTX 4090，Ryzen 9 7950X3D。未升级驱动或安装新计算依赖。
- preregistration/profile：`9905de6ee1308bf5ca4c15ab6cb6ca54e1f3c560`；formal regression：`040f8e915ddde0cfb5f903005f09a8aae7435e2b`；正式 GPU campaign_v2：`deb2821a976c5d0205270851d1bf692061bcc60d`。各结果保存 IDENTITY、输入与执行源码 SHA-256。
- F10A/F10B/F10C 分别 228/90/915 项历史文件及 13 个 F5 原输入核验；最终再次核验见 `HISTORY_PRESERVATION.json`。
- GPU 科学任务由独立 Supervisor 管理，采样全部 compute PID 和其他科学 Supervisor 任务。17:46:21–17:47:23 UTC 共105次采样未观察到冲突；不能保证采样间隙没有极短外部活动。没有停止其它项目。

## 真实数据路径审计

1. 原 `scripts/fast_pzt_f10c/runner.py::snapshot` 每个保存点读取14类数组：位置、声力、自由速度、复压力、接触诊断、active、保守壁距、壁面统计、CFD速度与梯度、CFD有效性与单元ID、近场诊断、自由位移。初始力尚未求值的语义原样保留，不能把它当作 t=0 的已求解力。
2. GPU 执行原逐泡声力、F2反馈、局部润滑、12轮PBD、CFD、壁面及每10 µs接受步连边。形态统计、Kabsch和持续拓扑事件都在CPU；GPU没有新的形态替代模型。
3. 201保存帧产生多次同步下载，最后质量摘要按成员读取诊断数组，确有重复读回；本轮保留以保证原质量语义，单独计时。没有把全部原组合成本归咎于下载。
4. 保存帧数组集中转换；旧路径逐帧逐成员调用 metrics、shape、SVD和SciPy稀疏连通图，重复算部分形状量。本轮将固定形状操作批量化，保持float64。
5. 原HDF5已按完整数组写gzip数据集，并非逐帧追加；source/medium/flow/bubbles/control/quality存JSON，部分身份是硬链接。输出末尾重开文件更新时间属性，并未重新积分或重新跑全套物理检查。
6. 形态组合计时还包括数组整理、事件下载/排序、质量摘要和metadata准备。本轮区分纯离线后处理与完整运行，不把组合计时全部称为Kabsch。

## 热点证据与实际优化

先读取原F10C B64的201×64×32轨迹，仅做CPU分析。未经profiling插桩的原后处理首次5.657 s、warm5.282 s；另一次正式新旧配对的warm原路径5.654 s、新路径1.299 s。不同运行不混成严格对照。

独立cProfile运行总9.433 s：12,864次frame_metrics累计4.500 s，22,272次SciPy图验证累计2.936 s，66,816次稀疏构造累计1.479 s；Kabsch函数累计0.662 s，SVD累计0.172 s。**主要成本是逐帧Python调用、重复稀疏图验证/转换和重复形状计算；Kabsch不是唯一或首要瓶颈。累计时间相互包含，不能相加作wall-clock分解。**原质心/Rg在frame_metrics/shape内部，未独立隔离的旧子项明确记为未测。

| 文件 | 本轮职责 |
|---|---|
| `src/mevo_mb/fast_sim/f10d_metrics.py` | 独立可选的NumPy批量形态、连通标签和接受步事件处理 |
| `scripts/fast_pzt_f10c/runner.py` | 显式postprocess/身份/存储参数钩子；默认保留原路径，原GPU循环未变 |
| `scripts/fast_pzt_f10d/postprocess_reference.py` | 原metrics/morphodynamics可重复回归入口 |
| `scripts/fast_pzt_f10d/pipeline.py` | 完整同步控制家族、完整输出、身份、质量与时间口径 |
| `scripts/fast_pzt_f10d/profile.py`, `validate.py`, `storage.py` | 实测热点、严格回归、无损输出候选测量 |
| `scripts/fast_pzt_f10d/benchmark.py`, `families.py`, `diversity.py`, `campaign.py` | 有限GPU对照、20条真实串行对照、四种子及资源监控 |
| `scripts/fast_pzt_f10d/analyze.py`, `figures.py`, `verify.py` | 数据清单、拒绝提案审计、四组图与最终准入 |

优化连通算法按最小粒子ID生成与SciPy完全相同的component标签，最多512帧一块；用batched float64归约/eigh/SVD和相同determinant修正计算形态，退化Kabsch使用原实现。沿用原特征值间隙判据屏蔽不可辨识主轴。当前正式路径为32个固定成员，未声称任意多分散/任意成员子集都已验证。

每个接受步的原始边事件仍保存并按真实时刻处理。只有成员标签完全没变时省去重复的集合比较，不用保存帧事件代替接受步事件。2 ms持续性及右删失定义未改。使用CPU向量化已足够，本轮没有增加Numba、GPU后处理、16核进程池或新物理近似；BLAS/OMP/MKL按原启动脚本固定单线程，避免小矩阵并行开销。

HDF5 byte shuffle候选60个数据集和属性/硬链接回归通过，但写入0.548→0.608 s，文件26.13→27.57 MB，未获收益，因此正式管线保留原gzip1。`storage_selection.json`记录拒绝。campaign_v1在候选选择尚未落实时启动，读取结果后主动停止了本轮自己的任务，全部残留文件和日志保留为`INTERRUPTED_IMPLEMENTATION_SELECTION`，不混入v2速度。它不是物理失败或一次完整性能测量。

## 定义、误差与F5真实释放

在正式服务器同环境对照中，B64和F5各39类检查有38类完全一致；只有主轴角末位舍入不同。B64最大2.01549e-7度；F5最大7.96196e-13度；新增不同初态的接受前缀检查最大2.71361e-7度。预登记角度绝对门槛1e-5度，其余长度1e-18 m、张量1e-28 m²、无量纲1e-10（相对门槛见protocol）。全部离散标签、成员掩码、原始边、时刻、持续事件、控制、原始逐泡坐标/力及质量保持原值。

独立软件检查包括40种连通图、空/全图、重复特征值、共线/零点云、单点、proper rotation、持续性边界；空点云显式拒绝。合成拓扑只作软件边界测试，不伪装成科学分裂算例。

F5输入来自`results/fast_pzt/f5/science_v1/main/trajectory.h5`，SHA-256 `20d3f3f56a4771a3a3fa8078f83556f2b4674f9e8cef5cc8bfac1e0c76e003c7`。没有重算600 ms物理轨迹：1201帧、保存间隔0.5 ms，674条原接受步边事件；33次persistent merge、11次persistent split、1次right-censored split全部一致。300 ms关声时最大连通团32，600 ms为7。事件表示几何连通变化，不是物理破裂/融合。

图4在WSL用同一保存轨迹重新计算两个实现以绘制叠线。WSL与服务器的NumPy/LAPACK环境不同，长短比有约2.13e-14的跨环境舍入差；按原登记门槛核验，保留`WSL_PLOT_RECOMPUTATION.json`。同服务器零差异结论不等于跨平台逐位承诺。

## 控制来源、独立初态与错误审计

OFF/S/O均取F10A `extended_10_v1/{OFF,S_1000,O_1000}/CONFIG.json`，完整幅相/开关/时刻经控制家族哈希分组。1 kPa表示原初态局部声压，使用原`EXPLORATORY_HIGH_DRIVE`源配置：S各源约1489.499 Pa；O约[1219.572,1221.522,1042.077,2192.951] Pa。不能改称291.2377 Pa历史源上限内配置。没有新增3 kPa扫描；一个batch仍仅共享一个同步日程。

原历史成熟轨迹的共同祖先、平移与控制分支不能算四个独立家族。因此按预登记种子1101–1104，在原F10A内部中心、半径12 µm球内随机顺序放置32个半径1.78 µm球，中心距至少2.05R。所有种子、生成尝试、初态坐标和哈希保存；六组配对Kabsch RMSD为10.31–11.94 µm，不是相同团的刚体平移/旋转。四组均由原O控制完整逐泡准备10 ms，不按成功挑种子或重试。

| 任务 | 结果 | 保存的首次错误时间区间（局部ms） | 拒绝提案中超限泡对（固定ID） |
|---|---|---|---|
| seed1101准备 | 完成，最大组32 | — | — |
| seed1102准备 | 完成，最大组32 | — | — |
| seed1103准备 | 失败2048 | (9.45,9.50] | (1,6)，重叠17.8056 nm |
| seed1104准备 | 失败2048 | (2.20,2.25] | (24,28)，重叠17.8313 nm |
| O分支，两个准备后初态 | 整批中止 | (4.85,4.90] | seed1102 (28,31)，重叠17.8313 nm |
| S分支，两个准备后初态 | 整批中止 | (5.00,5.05] | seed1101 (12,31)，重叠18.2601 nm |
| OFF分支 | 两条都完成10 ms | — | — |

2048来自原`vessel_kernels.py`最终提案检查：中心距<1.99R，即重叠超过17.8 nm。没有放宽。原`failure_members`诊断摘要未捕获该最终提案的局部原因；新离线审计读取保存的`rejected_interval/proposal_m`，补充具体sample_id和泡对证据，不改原quality或状态。GPU未保存首先置位的线程或精确子步，**只能给50 µs检查区间，不能伪造精确首错时刻**。

O、S各有一个另外成员没有在保存提案中观察到超限，但整批停止，确有潜在吞吐损失；它们的未来是否能独立成功未知，仍标ABORTED_BATCH，不算成功。6条失败/中止记录、接受前缀、原始错误和拒绝提案都在清单；未进行隔离重试。两个准备后状态只是10 ms时连通32泡，不据此宣称稳态成熟或实验真实密团。

## 数据生产语义与交付

`analysis/SAMPLE_CATALOG.json`逐成员保存路径/哈希、sample_id、batch_id、control_family_id、原初态家族、split_group、初态来源、原物理模式、数值状态、全部计时字段。370条目录记录包含配对/串行/重复基准，不是370个独立初态；同一seed所有准备和控制分支同组，无虚构train/validation/test拆分。

297个远程文件合计466,315,268 bytes已回传WSL并逐项SHA-256核验；包括失败及中止候选。`local_validate_v1`仅是源码未提交时的本地开发证据，明确不作为正式commit回归；正式服务器执行快照独立验证。最终目录另含图表、文档快照、数据清单与总manifest。

四组英文PNG/PDF位于`results/fast_pzt/f10d/figures/`，分别为后处理瓶颈、完整吞吐、多控制家族、形态拓扑验证；坐标和泡半径真实，图4只投影、不放大群内间距。绘图代码、输入SHA、输出SHA、PDF检查和人工视觉复核分别保存。总门禁`analysis/FINAL_GATE.json`区分工程通过、完整三控制多样化数据未通过及物理未校准。

进度文档`/home/lzy/projects/microbubble_aggregation.md`只在F10C后插入F10D叙述，原F1–F10C和标题结构不变；备份`.bak_F10D_20261011`，修改前后快照、diff及SHA见`progress_document/`。

## 复现入口与后续重点

在项目根目录用`PYTHONPATH=src:.`，远程固定commit由`python -m scripts.fast_pzt_f10d.remote <module> --label <new_unique_label>`部署并运行；历史label不可覆盖。纯CPU入口为profile、validate、storage；正式GPU的campaign要求离线回归已通过、设备空闲，包含B16/B64、三家族与四种子；archive回传核验，analyze生成清单，figures出图，progress更新文档，verify最终核验。不要为了复现图表重跑整个GPU campaign。

F10E优先在这些真实失败窗口补齐每成员错误原因与首错子步记录，然后验证**不改变原数值预算**的失败隔离/可追溯重试，并检查原接触迭代在独立密团上的稳健性。必须记录额外时间、失败尝试和对照，不能只报重试成功者。多控制分组本轮已有收益，无证据支持先上完全异步kernel。接受步拓扑仍是优化后最大CPU项，但不应先微优化它而忽略多样化数据失败率。

本轮没有新CFD/k-Wave、G1冻结团、多速率Q1/Q2/Q3、全局RPY或AI训练。快、数值一致、初态不同与真实物理校准是四件不同的事。

实现参考：NumPy官方[stacked SVD](https://numpy.org/doc/stable/reference/generated/numpy.linalg.svd.html)、[eigh](https://numpy.org/doc/stable/reference/generated/numpy.linalg.eigh.html)，h5py官方[dataset/filter pipeline](https://docs.h5py.org/en/stable/high/dataset.html)。本轮性能结论来自上述本项目实测，不来自文档推测。
