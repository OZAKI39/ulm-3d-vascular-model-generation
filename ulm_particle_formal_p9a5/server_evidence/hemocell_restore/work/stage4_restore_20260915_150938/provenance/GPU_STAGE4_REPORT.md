# GPU Stage4：有界 LBM 批量派发验证

最终生产候选：**STAGE4**；Stage4 接受：**YES**；收益等级：**ACCEPT_MAJOR**。
决定依据：ALL_FROZEN_GATES_PASS。本轮结束后停止纯流体优化；下一步为 Step3C GPU 正式长验证。本轮没有启动该正式长验证。

## 实施与冻结

已先保存 Stage3 完整算例源码、GPU Palabos src/externalLibraries 及 SHA256，再在独立 scratch 副本中变更。已有静态库、冻结输入和工具链只读复用；独立补丁为 STAGE4_PATCH.diff。新增非虚方法及嵌套描述类型不改变原生对象布局，编译时确认 AtomicAcceleratedLattice3D 仍为504字节。原逐单元 pull/collide/push 主体按原文提取；CollisionKernel<BGK,NoDynamics>、readBit 表达式、精度、streaming 索引及边界分支不变。框架 internal processors、atomic/lattice 时间计数及调用次序不变。Guo、halo、monitor、MPI reductions 均未优化。

LBM_DISPATCH_INVENTORY.csv 在优化前生成。305次 LBM + 2次halo + 1次Guo + 4次fluid/ghost reduction + 48次flux reduction = 360次普通step launch。每100步trace中的终点快照另有2次packing launch，故基线360.02。每个block的运行时兼容性和独立路径见 runs/*/diagnostics/lbm_runtime_compatibility.csv。persistent descriptor和cell routing只创建一次；两套原生population指针通过parity选择，完成批量kernel后按原逻辑交换原生buffers。

## 预先冻结的门槛

先200步正确性，在0/1/10/100/200比较全流体rho/velocity及每步24组flux、total/CV mass、rho范围、Mach和nonfinite。数值阈值逐字节复用Stage3，不因结果修改。随后各3次200步；提升<5%停止。达到5%才各2次1000步；提升<10%停止。达到10%才各2次5000步，重复相对差≤5%、相对1000吞吐下降≤5%，且5000提升≥10%、无主要迁移退化方可保留。没有10000步，也没有自动重试。

## 实测性能

| 区间 | baseline steps/s | batch steps/s | 提升 | 结论 |
|---|---:|---:|---:|---|
| 200步，预热100 | 48.009064411 | 103.218381441 | 114.9977% | PASS |
| 1000步，预热100 | 51.615828782 | 121.606383214 | 135.5990% | PASS |
| 5000步，预热500 | 52.368973448 | 125.880863075 | 140.3730% | PASS |

同口径 Stage3 CPU MPI12 5000步中位数为 17.692464626 steps/s。最终候选为 125.880863075 steps/s，相对CPU为 7.114942193×。CPU短sanity只作同环境核查，不把200步CPU结果冒充5000步长期成绩。

5000步两次batch吞吐的相对差为 **0.467244%**（门限5%）；5000步中位吞吐较1000步中位数提升 **3.515013%**，没有持续降速。稳定性门PASS。

## 正确性证据

汇总状态：PASS；共 **496,473** 项检查通过。全部逐检查证据为verification/*/ALL_NUMERICAL_CHECKS.csv.gz；里程碑汇总为STAGE4_CORRECTNESS.csv。明细沿用旧列名cpu_value/gpu_value；各comparison路径和汇总指定其实际reference/variant，多数是Stage3-equivalent GPU baseline对batch GPU，并非重新跑CPU。

| 比较 | 检查数 | 状态 | 全场逐字节一致 |
|---|---:|---|---|
| cpu_sanity_comparison | 7459 | PASS | False |
| initial_strong_correctness | 7465 | PASS | True |
| pair_1000_r1 | 37062 | PASS | True |
| pair_1000_r2 | 37062 | PASS | True |
| pair_200_r1 | 7459 | PASS | True |
| pair_200_r2 | 7459 | PASS | True |
| pair_200_r3 | 7459 | PASS | True |
| pair_5000_r1 | 185065 | PASS | True |
| pair_5000_r2 | 185065 | PASS | True |
| profile_correctness | 7459 | PASS | True |
| stage3_baseline_equivalence | 7459 | PASS | True |

## GPU机制证据

| 指标 | Stage3 pre | Stage4 post |
|---|---:|---:|
| launch/step | 360.020000000 | 56.020000000 |
| median kernel µs | 34.368000000 | 103.423000000 |
| active fraction | 0.710976558 | 0.560084445 |
| idle fraction | 0.289023442 | 0.439915555 |

迁移：pre 331939.84 bytes/step；post 337182.72 bytes/step。两者trace均为101–200且包含相同step200快照。去掉快照的101–199、以及102–199细分见PRE_POST_MEMORY_MIGRATION.json。GPU active是kernel占用的时间并集，包含kernel内部stall，不能解释为SM计算利用率。GPU active/idle、managed DMA、CPU同步等待与应用阶段范围重叠，不能相加。

预热后batch GPU资源采样：median=71.0%、mean=71.1304347826087%、peak=72.0%，样本69。这里实际使用两次5000步运行的预热后区间，完整scope字段保留在GPU_RESOURCE_TRACE.csv。

| 应用阶段（可加的父范围） | pre | post |
|---|---:|---:|
| LBM | 0.620934001 | 0.180228700 |
| halo_sync | 0.100625131 | 0.216524149 |
| Guo_fused_wall_inlet_outlets | 0.017321101 | 0.037929902 |
| MPI | 0.001162666 | 0.002328556 |
| monitor | 0.167871149 | 0.362376531 |
| output_logging | 0.089308158 | 0.195403395 |
| other | 0.002777794 | 0.005208766 |

Trace相对同一101–200无profiler窗口中位数的耗时比为 1.015214906385856。trace不作为生产steps/s。Guo墙/入口/出口属于未修改的融合kernel，不能从本trace拆出独占壁面占比。Occupancy/带宽硬件计数仍UNVERIFIED，未绕过宿主机计数器限制。

## 保留、拒绝与外推

依据真实5000步中位数作**PROJECTED**外推，仅含solver迭代，不含初始化/I/O变化/未来物理工况风险：

- 240000步：0.529601小时（PROJECTED，未实际执行）。
- 360000步：0.794402小时（PROJECTED，未实际执行）。
- 700000步：1.544670小时（PROJECTED，未实际执行）。

下一步：STEP3C_LONG_VALIDATION_ON_GPU；之后才是RBC_GPU_ARCHITECTURE_AUDIT。本轮RBC/IBM均未包含，自动Stage5纯流体优化禁用。

## 来源、完整性与限制

执行模式、重复次数和gate原文见STAGE4_EXECUTION_CONTRACT.json。durable supervisor无自动重启，终态见provenance/SCHEDULE_TERMINAL.json。每次均从零初始化，MPI1 GPU/CPU MPI12，正常计时无逐阶段新增同步；descriptor prefetch同步仅发生于初始化。

文件空间使用仅涉及本轮Stage4：运行成功并记录end-to-end后，快照gzip无损压缩与按raw SHA256去重。每个原文件raw SHA、压缩SHA及还原检查见LOSSLESS_STORAGE.json；不会将压缩耗时算入某个variant的solver成绩，也未修改旧档案。原始nsys-rep及sqlite按实际生成保留。

## 批量派发的直接机制证据

LBM每步kernel数从305降至1，所有2621304个原分配单元仍执行相同原生计算。原kernel每次17–19个CTA，新batch为5120个CTA，均为256 threads/block；这证明dispatch规模改变，不等同于测得occupancy。单次LBM kernel中位时间从34.368µs增至814.405µs，但每步LBM kernel总时间从10.509428ms降至0.815304ms，约12.890倍。

整体launch减少6.426633倍。idle占比从28.90%升至43.99%，但绝对idle从6.184522ms/step降至4.326786ms/step，下降30.038%；kernel active绝对时间下降更多，因此不能用idle比例上升否定吞吐改进。

| 组件 | pre ms/step | post ms/step |
|---|---:|---:|
| LBM | 13.286469 | 1.772582 |
| halo | 2.153132 | 2.129554 |
| Guo | 0.370629 | 0.373047 |
| MPI | 0.024878 | 0.022902 |
| monitor | 3.592032 | 3.564039 |
| output | 1.910976 | 1.921828 |
| other | 0.059438 | 0.051229 |

monitor、halo、Guo绝对时间接近原值，而占比上升来自LBM明显缩短。本轮按要求保留它们，停止进一步纯流体优化。输出占比来自每100步1次终点快照的短trace；5000步计时区间快照密度不同，不能把此占比当作5000步常态。

最终选择的是经过验证的GPU实现；本轮benchmark可执行仍保留5000步硬上限。正式Step3C长验证是下一阶段，本轮没有移除这个上限。

最终Stage3保护审核已PASS，完整记录见verification/STAGE3_FINAL_PROTECTED_AUDIT.log。本轮全量源码、运行文件及快照无损还原检查见verification/FINAL_INTEGRITY_AUDIT.json；归档校验见SHA256SUMS。
