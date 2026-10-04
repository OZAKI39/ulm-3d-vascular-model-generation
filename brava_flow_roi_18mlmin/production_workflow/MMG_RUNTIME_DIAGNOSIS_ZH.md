# MMG R1 运行成本诊断

**当前状态：MMG_R1_NOT_PRACTICAL_WITHIN_BOUNDED_RUNTIME。** 本轮已停止，无自动下一任务。尚无通过全部门槛的 R1 网格，不具备进入 manufactured CFD field gate 的条件。

本轮只处理 mesh。P0/P1 CFD、TetGen NW1、legacy particle baseline、旧1565 population 与旧 wall audit 均未重跑。冻结 A1、已安装 binary 及 surface/refs 保持原身份。

## 资源和构建

CPU：AMD Ryzen 7 7800X3D 8-Core Processor；8物理核/16逻辑核，NUMA=1。cgroup CPU `cpu.max=768000 100000`，即7.68核额度；allowed cpuset 0–15。MMG固定在允许核6，防止迁移；没有宣称获得宿主机排他预留。内存上限42.0557 GiB；审计时宿主可用RAM=56.85 GiB，不等于容器可分配上限。active/inactive文件缓存与不可回收估计分开记录，没有清缓存。

RTX4090，24564MiB，driver595.84，已装PyTorch CUDA13.2可用；详情见能力审计（历史引用 `MMG_RUNTIME_CAPABILITY.json`；结论见本报告和最终摘要，旧路径已清理）。MMG只链接libm/libc，没有OpenMP/MPI/CUDA运行时或对应符号；官方CMake配置与实际运行1线程一致。Release实际编译flags为`-O3 -DNDEBUG`，没有fast-math，**无需重建**。`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`。没有伪造GPU或多核MMG加速。[固定版本官方构建](https://github.com/MmgTools/mmg/blob/v5.8.0/CMakeLists.txt)

MMG binary SHA256：`8b0df6c26c9b0ab4d5a72ed550db7a9c4dad0e3bbe7c7169f52f5dbc72f9b59e`。已有GPUWall是粒子BVH/有限三角形查询接口，不是tetra topology/metric mesher。已实现的新网格QA路径使用NumPy批量、稀疏图和两个并行NumPy任务；本轮R1无最终输出，这条完整R1 QA路径未执行。box fixture和合成拓扑拒绝测试已执行。没有为GPU利用率开发新mesher。GPU计算时间/峰值新增内存均0；nvidia-smi既有占用不是本轮计算。

## D0：原 A1 的 10 分钟证据

input.mesh SHA=`4bc843481d9de15ac661f3f8c4f97d4fdbfcb73cfe78f478054345bdabdd6957`；input.sol SHA=`fdcba868eef55ed303cd2123eabceca2142d184ba97cfb5c33a9f709f1532bbc`。唯一MMG选项变化是verbose4→6；`-nosurf -nosizreq -hgrad -1 -hgradreq -1 -m 16000`全部保持。tensor、全部required实体、输入坐标/连接均保持。

native wall=600.340s；user=599.761300s，system=0.572001s，平均CPU=99.999%，线程数[1]。native child maxRSS=1.6366GiB（wait4峰值含exec启动继承页，另有运行中RSS采样）；voluntary/involuntary ctx switches=6/839；末次proc采样iowait=0.0s、read_bytes=102400、write_bytes=0。无OOM增量，无最终mesh。SIGTERM后MMG自身返回0且打印Program killed，所以exit0不算完成。

| major batch结束时间s | split | collapse | swap | move |
|---:|---:|---:|---:|---:|
| 14.084 | 57188 | 1965 | 247062 | 26024 |
| 65.913 | 286761 | 5089 | 1377575 | 132976 |
| 301.289 | 993531 | 31151 | 5635575 | 478160 |

上表是完整批次计数，不把同一批次的子日志和汇总重复相加。D0最后一个完整主要批次在301.289s，后约299s未有新批次完成日志；CPU持续运行不等于这段时间已经有可量化网格收益。总体明确做了适配操作，允许D1；没有最终质量轨迹或逐实体历史，不能据此确认重复refinement/coarsening cycle或质量停滞。因此gradation不改。

## D1：仅取消人工 required

保留`-nosurf`，只从ASCII输入删除RequiredVertices/RequiredTriangles两个块，前面vertices/tetra/triangles/refs字节完全一致；A1 metric逐字节相同。官方实现将nosurf边标为MG_REQ+MG_NOSURF，并区分真实required标记；选项含义不能代替输出表面检查。[官方实现](https://github.com/MmgTools/mmg/blob/v5.8.0/src/mmg3d/analys_3d.c)

已有1mm box新配置测试PASS：uniform、anisotropic输出points/tetra/triangles/refs与旧显式required输出数组完全相同。细方向比0.164202，切向比0.850165/0.894587，bulk比1.067103；表面坐标差0。仅执行新锁定配置，没有重跑旧fixture或覆盖旧输出。fixture证据（历史引用 `MMG_D1_FIXTURE.json`；结论见本报告和最终摘要，旧路径已清理） · D1计划（历史引用 `MMG_D1_PLAN.json`；结论见本报告和最终摘要，旧路径已清理）

首段1800s。到达门槛时，必须已有至少3个非零已完成操作日志、最近600s至少2条且最后一条距门槛≤300s，才允许同一进程延长一次至总3600s；否则停止。此规则在执行前记录，CPU/RSS增长不作为延长依据。外层监督额外时间仅供QA/收尾，不能增加native预算。

## D2：最后一个受限配置

D1因预算停止后才启动D2。唯一变化：volume h_theta在0–10µm由7.323538平滑增至25µm，10–25µm保持25µm，25–75µm恢复183.512362µm，之后bulk。hn=25µm及其原transition、hz=238.711713µm及原transition不变；hgrad/hgradreq=-1保持。冻结时间2026-10-03T13:04:09.134919+00:00，metric SHA `e6880dc8499b5023e5608afb01693bb1f45a09b850da1720145b9ce24fcad566`。15872个真实WALL节点保留A1 tensor逐字节相同；继承Float32壁面坐标的最大解析径向残差2.76855e-11m按拓扑WALL的d=0处理，没有移动坐标。轴向metric逐字节一致，法向相对差仅1.19e-15。全部37625个原P0节点确定性采样，SPD通过；MMG仍插值所提供tensor，并不连续调用解析尺寸律。分段采样数量{'d_le_0': 7316, 'd_0_10': 9178, 'd_10_25': 754, 'd_25_75': 4448, 'bulk': 15929}。因此只能验收实际输出高度，不能将连续曲线当成已实现的网格高度。冻结D2计划（历史引用 `MMG_D2_METRIC_PLAN.json`；结论见本报告和最终摘要，旧路径已清理）

## 结果及质量验收

| Case | native wall s | CPU % | child peak RSS GiB | 最终mesh | wall median µm | bulk/P0 |
|---|---:|---:|---:|---|---:|---:|
| D0 | 600.340 | 100.00 | 1.6366 | 无 | 未评估 | 未评估 |
| D1 | 1800.356 | 100.00 | 2.1594 | 无 | 未评估 | 未评估 |
| D2 | 1800.365 | 100.00 | 2.1604 | 无 | 未评估 | 未评估 |

### D1

状态 `BOUNDED_RUNTIME_NO_FINAL_MESH`；停止原因 `FIRST_GATE_NO_RECENT_OPERATION_EVIDENCE`。没有最终输出，surface/高度/bulk/质量全部未评估，不写成质量失败。延长决定见原始 EXTENSION_DECISION.json；这里只报告真实运行到的预算。
### D2

状态 `BOUNDED_RUNTIME_NO_FINAL_MESH`；停止原因 `FIRST_GATE_NO_RECENT_OPERATION_EVIDENCE`。没有最终输出，surface/高度/bulk/质量全部未评估，不写成质量失败。延长决定见原始 EXTENSION_DECISION.json；这里只报告真实运行到的预算。


已实现、待最终R1输出触发的完整QA包括全表面坐标、三角形/ref、真实tetra exterior、free/non-manifold edges、单连通体、zero/negative/数值near-zero、重复cell和all-WALL；本轮没有执行这套R1 QA。near-zero以`abs(det)/max_edge^3≤64*eps64`定义；大Euclidean aspect单独不判失败。图和报告不改变科学数据。MMG初始quality不作为最终metric-space quality；若有最终输出，原生质量原文与上下文保留。[官方输出质量计算](https://github.com/MmgTools/mmg/blob/v5.8.0/src/mmg3d/quality_3d.c)

本轮运行结论按预先记录的保守近期计数规则得出。D1/D2均在30min首段停止，没有实跑60min；因此不能声称已证明它们必然无法在60min内完成。两次在22.8min左右仍有交换计数，停止前仍占满CPU；这也不是进程挂死的证明。没有证据将失败归咎于OOM、I/O或MMG数学能力。

D1/D2前四个主要批次的完成时刻和计数在[比较JSON](MMG_RUNTIME_COMPARISON.json)中逐项列出；第四批约902.65/900.57s，差约0.2%，不是显著加速证据。只有一次每配置的有界观察，没有统计重复试验，不能推断所有网格/尺寸律的性能。

保留检查：本地（历史引用 `MMG_RUNTIME_PRESERVATION_LOCAL.json`；结论见本报告和最终摘要，旧路径已清理） · 保留检查：服务器（历史引用 `MMG_RUNTIME_PRESERVATION_SERVER.json`；结论见本报告和最终摘要，旧路径已清理） · 下一步工具比较（仅建议、未安装运行）（历史引用 `MESH_TOOL_FALLBACK_DECISION_ZH.md`；结论见本报告和最终摘要，旧路径已清理）

## 六个问题

1. D0 显示真实计算负载：CPU≈100%、单线程、I/O wait=0、连续完成多个操作批次；不能把计算时间归因于 I/O。 D1 与 D0 前3个可对照主要批次的操作数全部相同。未观察到去人工 required 后适配工作量减少。
2. 1 mm box 已实测足够：坐标差0、connectivity/refs通过，uniform及anisotropic输出数组与旧required版本完全一致。 R1 D1 表面检查 NOT_EVALUATED_NO_FINAL_MESH。
3. D1 未在本次按证据允许的预算内生成最终网格；实际耗时 1800.356 s，不把30min停止写成实际运行了60min。
4. D2 实际耗时 1800.365 s，最终网格 未生成；无法从截尾时间声称获得完成加速比。
5. 见上表；没有输出的候选全部未评估，不能以目标metric代替真实altitude。
6. 本轮已停止，无自动下一任务。尚无通过全部门槛的 R1 网格，不具备进入 manufactured CFD field gate 的条件。

## 实现、可复核性和运行问题

监控使用userspace `/proc`/psutil、wait4/getrusage和带UTC/elapsed的stdout/stderr。实例禁止kernel profiler，本轮未用perf/eBPF，也未申请root或修改cgroup。第一次启动缺少/usr/bin/time，在MMG启动前失败；改用上述等价记账。D0 wrapper最初只查子进程导致采样空，已附加只读sampler到同一MMG进程；native详细采样从约79s起，完整CPU/context-switch/RSS总账覆盖全程。原始错误和原始PROFILE保留，根目录D0报告明确给出修正来源，没有重启D0或补造早期样本。

D1首次box启动中，监控字段status重名触发Python异常，没有完整输出；原始目录保留为D1_fixture_instrumentation_failure。修复字段名及异常时子进程组清理后，新box配置通过。预算/近期进展判据、D2只改周向metric、实际拓扑损伤拒绝、输入块删除的4项测试通过。修复上述运行问题未改科学配置或MMG算法；D2的周向metric改变按独立冻结计划执行。

监督程序的exit0也只代表收尾脚本完成，不代表网格成功；以native stop_reason、最终输出存在与完整QA共同判断。

保留核对：本地8项正式输入、205项冻结证据、330项粒子源码/测试全部相符，3896项历史文件SHA与6936项size/mtime也相符。服务器正式8项/冻结205项/MMG与SimVascular二进制全部相符；实际部署的150项粒子源码/测试相符，另180项旧测试不在服务器部署中，未声称远端核验了330项。首次检查将未部署测试计为失败，原始结果保留；后续明确报告覆盖缺口，没有补拷测试或掩盖已有文件的差异。

## 下一阶段GPU执行计划（本轮未执行）

只有后续mesh winner通过manufactured CFD的u/G/s/curl及strict containment field gate后，才运行已有RTX4090 production GPUPopulationEngine。按冻结2个半径×5个gap组成10 canonical singleton cases；每档dt批量处理10例，dt=0.5/0.25/0.125ms，合计30个参数组合。网格/field/BVH一次上卡、状态数组常驻、同档统一批量推进，减少kernel launch和host-device搬运；保留既有float64、精确壁面和contact/event规则。若合批多个独立算例，必须明确隔离scenario，避免引入原本不存在的MB–MB作用，不能直接把10个oracle泡混作一个物理群体。少量最小gap/中间gap/远壁CPU reference负责parity；完整性能需CUDA events+同步边界、wall/CPU时间、峰值VRAM和传输统计。具体半径/gap按冻结合同读取，不重新抽样。本轮0 particle/oracle/trajectory/population，正式默认dt未改。

[300 dpi 图文入口](OPEN_RESULTS.html) · [机器可读比较](MMG_RUNTIME_COMPARISON.json) · D0完整profile（历史引用 `MMG_D0_PROFILE.json`；结论见本报告和最终摘要，旧路径已清理）

## 清理后的历史定位

本路线CLOSED_FOR_NOW。原A1单次1800.628秒无最终网格、tree RSS峰值2.4355 GiB；D0/D1/D2结果及MMG 5.8.0源码commit、binary SHA在比较JSON。这里的“本轮”指历史运行，所有专用构建、适配器、原生日志及逐秒监控已退休。保留2张PNG，不恢复计算。
