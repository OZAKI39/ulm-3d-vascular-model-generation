# GPU_POC_STAGE_3 — 优化后正确性与持续性能审计

本阶段状态为 **PASS**；持续正确性 **PASS**，数值差异趋势 **STABLE**，持续性能 **STABLE**。是否具备后续长程纯流体验证条件：**YES**；建议平台：**GPU**。本报告没有宣布 Step3C 正式长验证完成，也没有运行 RBC。

10000 步 CPU/GPU 持续速度分别为 **17.643724 / 52.336876 steps/s**，solver 倍率 **2.966317×**；包含启动、初始化、输出和正常退出的实测端到端倍率为 **1.534809×**。判断依据为本轮新测 CPU MPI12 与 GPU MPI1，没有沿用旧 17.618650 作为分母。

**监控中断与样本完整性**

原gpu_5000_r2的solver输出实际到达5000，安全检查通过；外部监控只记录到约3014步，正常退出码、精确E2E与后段资源证据缺失。原调度器因此退出FAIL，未启动10000。完整原始文件散列和真实限制见 verification/SESSION_INTERRUPTION_AUDIT.json；不以文件mtime推造精确退出时间。

该中断样本的内部solver速度仍有直接文件证据：52.348399 steps/s；其完整E2E和晚段资源为UNVERIFIED，因此不作为具备完整审计的正式重复样本。

原冻结比较器对中断样本的数值量比较均通过，但依赖正常退出/保存绑定证据的两项审计检查不能通过；详见 verification/INTERRUPTED_OUTPUT_ASSESSMENT.json 和 provenance/interruption_output_comparison/。这些证据没有被改写为完整PASS样本，也没有混入正式E2E统计。

用户随后明确授权一次额外GPU5000执行，原失败调度回执与整组输出保存在 provenance/ORIGINAL_SCHEDULE_TERMINAL.json 和 interrupted_runs/。替代样本从零初始化，沿用原冻结二进制、比较器和止损门；没有自动重试。GPU5000实际执行数为3，其中2个具备完整审计证据的样本用于正式统计。新调度器通过scratch私有supervisor写持久日志，关闭自动重启。

**科学与实现冻结**

仅使用 Stage2 已通过的 CPU-preferred AL metadata + GPU accessed-by 优化。GPU adapter 和 MPI/monitor 数学实现保持原 SHA256；Stage3 仅增加受限 horizon、相同的额外快照、初始化/完整迭代/持续迭代计时及 100 步时间窗口。没有新性能优化，没有修改正式 HemoCell、Palabos baseline、Step3B 或 Step3C。

| 冻结记录 | 内容 |
| --- | --- |
| GPU Palabos commit | 4127697e90169bbef982295f1d1c933cf6e90caa |
| CPU Palabos commit | 05712164d940a42e06afdd705249912fa0c49f14 |
| 编译器 / GPU flags | nvc++ 26.5-0 64-bit target on x86-64 Linux -tp znver5 ; NVIDIA Compilers and Tools; Copyright (c) 2026, NVIDIA CORPORATION & AFFILIATES.  All rights reserved. / -stdpar -Kieee -gpu=cc120,nofma -fopenmp -O3 |
| GPU adapter SHA256 | 4b9dcfa512acdb3e8853d0086585903750883848a29ed96c99768515f0fecf4c |
| 优化 patch SHA256 | 9f46077323886448c947cf1eea279c7c9a0ad458459cbe604fe4a5ae4a72df92 |
| 原数值合同 SHA256 | 5e78c61974be7ca2a52b3add288858330804c8d74d47ae909fced717d1f02456 |
| 冻结输入映射 SHA256 | da8b9a5f75acee8323078050dd0b633f7b0da6456e6c4c8c3c43eab87ead12b6 |
| Stage2 归档 SHA256 | 032c94b84412403e52cc67faf1300c91a9fe6c250babb55a4e2f9c64667a485b |

Palabos 为不带 .git 的源码归档；commit 来自已冻结的 Stage1/2 来源记录，源码完整性由既有源树哈希核查。逐个 source hash、编译记录及本轮 harness diff 见 STAGE3_FROZEN_GPU_IMPLEMENTATION.json 与 provenance/。

geometry/STL、rho=1056 kg/m³、nu=3.27e-6 m²/s、dx=1.9989918081065344e-7 m、dt=2.0366810646671923e-9 s、tau=1、Qtarget=2.7369132390905703e-15 m³/s 均保持。Step3C multiplier=1.1197286861799598；三个出口压力为 14.544978101274268、132.20454922317552、-13.700626673311461 Pa。Guo、入口/出口算法、double 精度、原初始条件、10 步 ramp、control volume、24 组测量平面及 monitor 公式保持。

原比较 JSON 全文保留。它所带的旧 benchmark/run-limit 元数据由用户本轮明确要求的独立 STAGE3_EXECUTION_AND_ANALYSIS_CONTRACT.json 覆盖；没有修改比较阈值、科学参数或正式文件。执行/分析规则与源码在首个 run 前封存，散列见 provenance/execution_hashes.json。

**从零开始、公平性与计时口径**

每次都启动新的 MPI/solver 进程，建立几何、lattice、边界与 GPU 状态，从 iteration 0 开始；没有 restart 或跨 run 复用初始化状态。CPU 使用 12 个独立物理核心，GPU 使用 MPI1。相同 launcher 选项和线程限制保留，内存绑定警告不被误写为 NUMA 绑定成功；实际核心亲和性单独核查。文件系统/CUDA 缓存未强制冷却，“从零”指模拟状态重新建立。

所有必要 safety 和全部 24 flux 分组逐步执行，flow/safety CSV 逐步记录，控制台日志频率一致。两侧均保存覆盖到的 0、100、200、1000、5000、10000 步 rho/u 与 sampled-field 快照：保留 Stage2 的 200 步快照，再加入本轮要求的长 horizon checkpoints。没有通过减少 GPU 输出取得优势。不同 horizon 的稀疏快照成本会被不同步数摊薄，所以跨 horizon 的速度差不能全部解释为 kernel 变快。

| 指标 | 定义 |
| --- | --- |
| INITIALIZATION_SECONDS | MPI 启动前的 parent monotonic 时间至 step0 必需检查/输出完成、准备开始 step1；包含 launcher 等待、MPI_Init、几何/lattice、GPU setup、边界与首个必需同步 |
| 初始化明细 | 另记 launcher_to_main、MPI_initialization、main_entry_to_ready；不同层级存在包含关系，不重复相加 |
| 持续 solver | 1000 case 为完成 step100 后至 step1000；5000 为500后至5000；10000为1000后至10000。分别900/4500/9000步 |
| full-iteration average | step1至终点的全部迭代时间，含 warmup、全部检查、相同快照与窗口记录 |
| END_TO_END_SECONDS | Popen 前至 MPI/solver 正常退出；包含初始化、所有迭代、必需诊断/输出、正常 shutdown |
| 未计入 solver E2E | 启动前冻结校验/输入路径物化、运行结束后的 CPU/GPU 离线比较及无损存储去重；两侧同口径 |

普通 benchmark 关闭阶段 profiler；没有额外每阶段 device synchronize。持续计时边界保留原必要 GPU sync/MPI barrier，输出策略两侧一致。低开销 100 步 wall-time 窗口用于检查运行内是否退化。所有精确命令、环境、二进制 SHA256、资源轨迹与正常/异常退出记录均保存在 runs/<name>/。

**各 horizon 的已收齐样本中位数（未满规定重复数时为临时值）**

| 步数 | CPU/GPU 重复数 | CPU steps/s | GPU steps/s | solver 倍率 | CPU E2E秒 | GPU E2E秒 | E2E倍率 | 状态 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1000 | 3/3 | 17.721833 | 51.688852 | 2.916676 | 241.933992 | 319.362893 | 0.757552 | COMPLETE |
| 5000 | 2/2 | 17.692465 | 52.198655 | 2.950333 | 469.639619 | 395.423571 | 1.187687 | COMPLETE |
| 10000 | 2/2 | 17.643724 | 52.336876 | 2.966317 | 753.288572 | 490.802705 | 1.534809 | COMPLETE |

速度中位数与端到端时间中位数分别计算；solver speedup=GPU速度中位数/CPU速度中位数，E2E speedup=CPU总时间中位数/GPU总时间中位数。样本数按用户要求：1000各3次，5000和10000各2次；10000两次GPU差异超过均值的5%才追加一次。原始逐次值、min/max、完整迭代平均和初始化时间见 CPU_RUNS.csv、GPU_RUNS.csv、PERFORMANCE_BY_HORIZON.csv。

10000 步首次两次 GPU 的相对差异为 0.3451%；是否触发第三次：False。不存在自动 solver 重试。

5000 与10000 的 GPU 持续速度相对变化为 0.2648%，按事前 ±5% 规则分类为 STABLE。原短区间 46.589479765553 steps/s 仅作参考，持续与端到端结论都来自上述新测数据。

**初始化固定成本、盈亏平衡与投影**

| 后端 | launcher至main中位秒 | MPI_Init中位秒 | main至可迭代中位秒 |
| --- | --- | --- | --- |
| CPU | 120.445605 | 60.437071 | 65.535099 |
| GPU | 120.641112 | 60.415315 | 176.791791 |

上表 MPI_Init 已包含在 main至可迭代 中，不能再加一次。两侧都保留本实例实测的 launcher/MPI 等待；本轮没有为某一侧移除这类固定成本。

直接测得初始化中位数：CPU **185.971242 秒**，GPU **297.432879 秒**。该值包含本实例的进程启动/MPI固定等待，不能将所有时间归给 GPU device allocation。拟合中的 A 还包含正常退出和其他固定/平均输出成本，故不必与直接初始化中位数完全相等。

| 后端 | 拟合 A 秒 | 拟合 sustained S | 拟合残差 RMS 秒 | 最大绝对残差秒 |
| --- | --- | --- | --- | --- |
| CPU | 185.502179 | 17.609503 | 0.357288 | 0.522212 |
| GPU | 299.822551 | 52.350403 | 0.645034 | 0.994062 |

模型 T_CPU(N)=A_CPU+N/S_CPU、T_GPU(N)=A_GPU+N/S_GPU 用本轮各次有效、无 profiler 的 E2E 数据最小二乘拟合。拟合交点为 **3033.54 步**；这是模型值，不是恰好在该步数实际运行测得。

实测约束：Observed median GPU advantage first at 5000; prior losing tested horizon 1000。只报告测试点所支持的区间，不把模型交点当作直接测量。

没有增加额外 targeted horizon。只有收齐规定样本且测试点能夹住交点时才报告正式实测区间；未完成时不推测缺失值。拟合残差、全部输入点与整数模型阈值保存在 BREAK_EVEN_ANALYSIS.json。重复次数有限，未宣称精确到单步的硬性性能保证。

长时间投影状态为 **PROJECTED**。只有5000/10000速度基本稳定、正确性通过且无实质差异增长时才投影；使用当前10000步持续速度与本轮拟合固定开销，未使用旧CPU/GPU速度。

| 步数 | CPU小时 | GPU小时 | 预计E2E倍率 | 性质 |
| --- | --- | --- | --- | --- |
| 100000 | 1.625900 | 0.614034 | 2.647900 | PROJECTED |
| 240000 | 3.830020 | 1.357083 | 2.822244 | PROJECTED |
| 360000 | 5.719266 | 1.993983 | 2.868263 | PROJECTED |
| 700000 | 11.072129 | 3.798531 | 2.914845 | PROJECTED |

所有100000/240000/360000/700000结果都只是 PROJECTED，并未执行。这些投影假设当前数值负载、monitor与输出成本保持；正式长验证可能有不同诊断/输出要求，不能把投影视为已完成验证或确定收费时长。

**数值一致性、误差增长与止损**

沿用原阈值：rho normalized L2与rho min/max/mean绝对差≤1e-8；velocity normalized L2及最大速度归一化差≤1e-5；全部24组flux绝对差/Qtarget≤1e-4；总质量和CV质量相对差≤1e-8；NaN/Inf为0；rho∈[0.99,1.01]、物理/ghost Mach≤0.05，ghost rho为正。每个全场文件必须具有完全一致的182694个整数节点ID、排序、有限rho/u。

| GPU run | 步数 | 检查数 | 失败数 | 数值结果 | 趋势 |
| --- | --- | --- | --- | --- | --- |
| gpu_10000_r1 | 10000 | 370068 | 0 | PASS | STABLE |
| gpu_10000_r2 | 10000 | 370068 | 0 | PASS | STABLE |
| gpu_1000_r1 | 1000 | 37062 | 0 | PASS | STABLE |
| gpu_1000_r2 | 1000 | 37062 | 0 | PASS | STABLE |
| gpu_1000_r3 | 1000 | 37062 | 0 | PASS | STABLE |
| gpu_5000_r1 | 5000 | 185065 | 0 | PASS | STABLE |
| gpu_5000_r2 | 5000 | 185065 | 0 | PASS | STABLE |
| gpu_profile | 200 | 7459 | 0 | PASS | STABLE |

每步标量与 safety 记录均逐项比较，完整检查明细保存在 verification/<run>/ALL_NUMERICAL_CHECKS.csv.gz；关键 checkpoints 在 SUSTAINED_CORRECTNESS.csv，所有量的全程最大误差在各 MAXIMUM_ERRORS.csv。比较的是相同 CPU transient，不要求10000步达到稳态，也没有重新校准Qtarget或出口压力。

| step | rho L2 | velocity L2 | 最大速度差 | mass relative | CV mass relative | max flux/Qtarget |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 |
| 100 | 0.00000000e+00 | 6.28628962e-16 | 0.00000000e+00 | 7.86266869e-16 | 1.06084613e-15 | 6.34106323e-15 |
| 200 | 0.00000000e+00 | 6.89571867e-16 | 0.00000000e+00 | 2.88297642e-15 | 2.91732477e-15 | 3.17053162e-15 |
| 1000 | 8.99779162e-19 | 8.42368372e-16 | 1.63942550e-16 | 6.55218842e-16 | 7.95630313e-16 | 3.60287684e-15 |
| 5000 | 1.03895750e-18 | 1.67653375e-15 | 9.58831112e-16 | 1.17937350e-15 | 1.06082209e-15 | 4.03522206e-15 |
| 10000 | 2.20393450e-18 | 1.88216695e-15 | 0.00000000e+00 | 2.88287746e-15 | 3.18242706e-15 | 5.04402757e-15 |

上表取同一 checkpoint 所有正常 benchmark 对照中的最大误差，不混入 profile run。事前趋势规则：任何原 tolerance 失败即 FAIL；≥5000 checkpoint误差同时超过其step1000值4倍、且超过原tolerance的1%时标记 GROWING；否则为 STABLE（在本轮合同与测试范围内无实质增长）。原始误差大小仍完整显示，STABLE 不等于数学上保证无限时间无误差。

5000步止损门状态：**PASS**。在进入10000前核查数值趋势、相对1000的速度、5000运行内前后1000步速度、steady thirds的VRAM/RSS变化及新短trace的迁移量。异常理由为 []。门的阈值在结果前冻结，没有看到结果后调宽。

晚期5000步的managed迁移字节未单独trace，为 UNVERIFIED；其数值、速度与内存趋势单独实测。短trace不能证明任意长时间都没有迁移。凡观察到止损条件即不再启动10000 case；失败不会通过补参数或自动重跑抹去。

**资源与优化后新 profile**

| GPU run | 首个预热后1000步速度 | 最后1000步速度 | 后/前速度比 | VRAM后/前中位差MiB | RSS后/前中位差GiB |
| --- | --- | --- | --- | --- | --- |
| gpu_10000_r1 | 52.943931 | 51.991087 | 0.982003 | 0.000000 | 0.000000 |
| gpu_10000_r2 | 52.835111 | 51.717200 | 0.978842 | 0.000000 | 0.000000 |
| gpu_5000_r1 | 52.115005 | 51.736867 | 0.992744 | 0.000000 | 0.000000 |
| gpu_5000_r2 | 52.161782 | 51.776887 | 0.992621 | 0.000000 | 0.000000 |

以上为运行内实测描述，不新增或调整判定阈值。最后1000步包括规定的终点快照。CPU对应值、steady thirds资源中位数与样本数见 SUSTAINED_WINDOW_DIAGNOSTICS.csv。

| 采样范围 | 样本数 | GPU util median | mean | peak | VRAM peak MiB |
| --- | --- | --- | --- | --- | --- |
| warmup后正常GPU迭代 | 557 | 82.00 | 80.78 | 82.00 | 2129.00 |
| 正常GPU全端到端（含启动） | 2672 | 0.00 | 18.70 | 86.00 | 2129.00 |

整个运行每秒记录GPU utilization、VRAM、求解进程CPU使用率、进程RSS及可用RAM，见 GPU_RESOURCE_TRACE.csv。最终利用率键采用warmup后迭代样本；显存峰值采用正常GPU全程峰值。完整起停区间统计另列，避免用短时高利用率掩盖初始化期间的等待。进程RSS求和可能重复计算共享页，不等同于唯一物理内存占用。

只另运行一次200步短profile，采集101–200共100个完整timestep，没有trace10000步。新GPU kernel活动占比 **71.0977%**，无kernel活动占比 **28.9023%**；每步launch **360.020**，中位kernel持续 **34.368 µs**。

该trace区间elapsed / 三次普通1000运行的匹配101–200区间elapsed中位数 = **1.030607**。trace会影响性能；它用于归因，未混入正式速度中位数。GPU active包含kernel内部等待/stall，不是有效SM compute utilization。

| 新测应用阶段 | 占timestep |
| --- | --- |
| LBM | 62.0934% |
| halo_sync | 10.0625% |
| Guo_fused_wall_inlet_outlets | 1.7321% |
| MPI | 0.1163% |
| monitor | 16.7871% |
| output_logging | 8.9308% |
| other | 0.2778% |

上表为新profile的互斥应用范围，各行合计100%；memory migration、host gap、kernel/API等待都是这些范围内部的另一个视角，不重复相加。Guo wall与入口/三个出口在同一原kernel内融合，生产独占份额为UNVERIFIED，表中只给新测parent总时间，不复制Stage2分类计时冒充新测值。

新测 LBM 内非kernel主机间隙占整个trace **12.9791%**；launch API自身仅占 **3.4653%**。主机间隙不能全部称为launch API耗时。完整kernel名称、grid/block形状和持续时间见 POST_OPT_KERNEL_SUMMARY.csv。原Nsight Systems trace与SQLite均保留。

短profile回执的timed_seconds包含cudaProfilerStart/Stop控制成本，不能用其派生steps/s代表生产速度；阶段占比与扰动比使用100个NVTX timestep的实际capture wall区间。profile run在GPU_RUNS.csv明确标记normal_benchmark=False，未参加任何性能中位数或盈亏平衡拟合。

主机在CUDA阻塞同步API中的等待占比为 73.6661%，它大量覆盖kernel运行时间，不能解释为同等比例可消除的同步损耗。

Nsight Compute宿主硬件计数权限在上一阶段已明确不可用，本轮没有重试或修改宿主权限；occupancy、硬件带宽/compute利用率仍为UNVERIFIED。

新trace可见的managed DMA活动区间并集为 1.870489 ms，占capture wall的 0.0874%；其中 1.206284 ms 与kernel重叠。这只是可见复制活动，不代表全部page-fault服务或远程metadata访问等待，不能直接与应用阶段占比相加。

| 迁移口径 | Stage2 before bytes/step | Stage3 after bytes/step | reduction factor |
| --- | --- | --- | --- |
| 含各自规定快照的采集区间 | 8074280.96 | 331939.84 | 24.324531 |
| 去掉末端快照的普通步 | 7385004.41 | 8853.98 | 834.088690 |

旧trace为50步含1次快照，新trace为100步含1次快照，故同时提供去掉末端快照的比较；新step101紧接step100快照，还额外给出102–199区间。数值均来自实际CUDA memory事件，不通过速度猜迁移量。CPU常驻metadata的GPU直接读取不必表现为managed migration；未测PCIe读取总流量，不能宣称所有host/device总线传输消失。

显式H2D=0.00、显式D2H=512.00 bytes/step；FULL_FIELD_TRANSFER_EVERY_STEP=NO，完整lattice writeBack调用为0。离散checkpoint输出保留。

| 新瓶颈排名 | 组件 | 新测占比 |
| --- | --- | --- |
| 1 | LBM | 62.0934% |
| 2 | monitor | 16.7871% |
| 3 | halo_sync | 10.0625% |

**第二次优化投资判断**

SECOND_OPTIMIZATION_WORTHWHILE = **YES**；候选：Batch LBM block dispatch while reusing original collision mathematics。本轮没有实现第二次优化。新首要组件被完全消除时的Amdahl上界为 **2.638063×**；该上界不意味着可以低成本完全消除。

候选成本为 MEDIUM，科学风险 MEDIUM，现有Palabos能力复用 PARTIAL。提案为 Expose a scratch-only immutable block job view and end-of-batch buffer-swap hook; reuse original collision functor and exact pull/collide/push ordering. Preserve original block partition, halo sequencing, framework processors/time counters, Guo and all monitors.。

| 条件情景，非实测预测 | 整体倍率 |
| --- | --- |
| 消除LBM非kernel间隙的50% | 1.069399 |
| 消除LBM非kernel间隙的80% | 1.115863 |

设备只读属性为 170 个SM，单个LBM launch只有17–19个CTA，同时可覆盖SM的上界约 11.1765%。这是grid覆盖上界，不是实测occupancy。LBM kernel自身占比 49.1141%，因此机会还包括缩短kernel时间。现有population/tmp buffers为private，需scratch batch view与swap hook，并保持框架processor/time-counter语义；不存在现成drop-in开关。

| 跨block batch的条件情景，尚未实施 | 整体倍率 |
| --- | --- |
| LBM kernel加速1.5倍，同时减少25%的LBM host gap | 1.244031 |
| LBM kernel加速2.0倍，同时减少50%的LBM host gap | 1.450254 |

这些约24–45%的情景支持一项有界原型的投资判断，不构成实际收益预测。真实收益仍为UNVERIFIED；没有把170/19的SM数量比直接当作速度倍率。

小grid已由trace直接证实，但它不能单独确定kernel内部stall来源。CPU常驻metadata的GPU远程读取延迟也可能贡献时间，当前没有硬件计数器证据将二者分离。因此kernel加速1.5/2倍是待验证假设，不能作为承诺或成本预算依据。

情景只量化值得进一步验证的时间余量，不保证实现能够获得相应收益。若trace相对匹配普通窗口的扰动超过20%，不据此支持生产收益预测。建议同时受用户20%新瓶颈门槛与LOW/MEDIUM成本门槛约束；如果收益上界很小，应优先停止抠纯流体性能并进行正确性/正式验证。具体理由与边界见 SECOND_OPTIMIZATION_ASSESSMENT.md。

**最终边界与证据**

GPU_READY_FOR_LONG_PURE_FLUID_VALIDATION=YES；RECOMMEND_STEP3C_LONG_VALIDATION_ON=GPU；READY_FOR_RBC_GPU_ARCHITECTURE_AUDIT=YES。RBC_GPU_FEASIBILITY=UNVERIFIED，RBC_INCLUDED=NO。架构审计准备度不代表允许直接加入RBC。

所有结束运行的原始输出先正常写出并计入E2E，随后才对SHA256完全相同且属性一致的Stage3快照做硬链接去重；所有路径与字节保留。对应证明在 provenance/*_snapshot_dedup.json。没有修改旧Stage2归档或通过丢弃输出节省计时成本。

最终 source/binary、所有run清单、1541个冻结输入、原Stage2及其受保护旧源树一并校验；结果记录在 verification/。SHA256SUMS封存全部归档文件。目录使用既有实例内依赖链接，不宣称是独立portable bundle；运行结果与正式源项目的远程发布/Git同步不在本次任务范围。

| 文件 | 用途 |
| --- | --- |
| STAGE3_FROZEN_GPU_IMPLEMENTATION.json | 固定GPU实现、commit、compiler、patch与合同散列 |
| STAGE3_EXECUTION_AND_ANALYSIS_CONTRACT.json | 运行前步数/重复/输出/比较趋势与止损规则 |
| SUSTAINED_CORRECTNESS.csv | 关键checkpoint的rho/u/flux/mass与安全比较 |
| CPU_RUNS.csv / GPU_RUNS.csv / PERFORMANCE_BY_HORIZON.csv | 逐次与中位性能、初始化、完整迭代和E2E |
| BREAK_EVEN_ANALYSIS.json / LONG_RUN_PROJECTIONS.json | 拟合、残差、实测夹逼与明确标记的投影 |
| POST_OPT_PROFILE_SUMMARY.json / POST_OPT_KERNEL_SUMMARY.csv | 新profile的阶段、活动区间和kernel证据 |
| POST_OPT_MEMORY_MIGRATION.json / POST_OPT_BOTTLENECK_RANKING.json | 新迁移量与新瓶颈排名 |
| SECOND_OPTIMIZATION_ASSESSMENT.md | 成本、风险、Amdahl与投资判断；无新实现 |
| GPU_RESOURCE_TRACE.csv / GPU_RESOURCE_SUMMARY.json | 全程资源原始样本与清晰分域的统计 |
| runs/ / verification/ / profiling/ | 原始终态、场文件、全检查压缩CSV与原始trace |
| FINAL_TERMINAL_SUMMARY.txt / SHA256SUMS | 终端摘要与完整封存校验 |
