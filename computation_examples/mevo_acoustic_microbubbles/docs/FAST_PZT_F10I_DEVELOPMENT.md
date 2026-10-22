# F10I：B2连边敏感性、误差分离与GPU可行性

日期：2026-10-12。独立分支：`dev/fast-pzt-f10i-b2-edge-events-gpu-20261012`。

**结论：当前O/seed1102、32泡、5 ms内，微小间隙差被缓慢阈值跨越放大为130 µs事件差；收紧非线性接触求解容差不足以消除它。团始终连通，位置与整体形态有独立数值稳定证据，但原严格事件验收继续失败。GPU固定状态原型实际运行，32泡比较误差很小，内部收敛仍失败，按预登记停止扩张。**

这不是新的接触物理校准，也不是全部形态收敛或生产通过。`training_ready=false`。

## 基线、范围与执行身份

- 从F10H最终提交 `9419837e0ae1ac0a2c6629c63f9a69027b2fef8f` 建独立worktree，未回退F5。
- WSL：`/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10i_20261012/computation_examples/mevo_acoustic_microbubbles`。
- 服务器：`/workspace/mevo_acoustic_microbubbles/f10i`，既有独立 `.venv`；RTX 4090、Ryzen 9 7950X3D。
- 13份原始输入SHA与F10H的213个清单文件逐一核验。`BASELINE_AUDIT.json`记录HEAD、worktree、旧严格FAIL、真实接受轨迹与ID、适用AGENTS；服务器指南SHA为 `1d264a8e24863e481852969f4bedd2834abfed79c68dab66d13e9d61a2a7360b`。
- 原物理刚度0.04641319518378658 N/m、R=1.78 µm、原声力/CFD/控制/近场模型不变。原重叠17.8 nm、单次修正178 nm、事件10 µs门槛不变。默认生产求解器未改。
- 仅新增一条从原始初态与完整控制开始的2 ms、10 µs、容差1e−14 m轨迹（200个接受步）；其余为旧数据分析、固定状态重解。2.5 µs、新10 ms、新seed、B64、CFD/k-Wave/APECSS、AI训练均未运行。

|执行|不可变提交|Supervisor标签|结果|
|---|---|---|---|
|唯一短容差控制|`fc53eadf11f4f95957c1d8750a09e6edbfcc3c95`|`mevo_fast_pzt_f10i_tolerance_v1`|200/200接受，2 ms完成|
|GPU原型v1|`f663a2dc454f1e6760dde85c561e52775dbb2ee6`|`mevo_fast_pzt_f10i_gpu_v1`|缓冲容量属性错误，0个科学算例完成|
|GPU原型v2|`458dabada891b218d8dc91b256a6de6c3d7d496d`|`mevo_fast_pzt_f10i_gpu_v2`|18小算例匹配，32泡内部收敛失败|

`REMOTE_VERIFICATION.json`：29份远程文件回传SHA通过；上述三提交分别446、446、442份执行快照文件对Git archive验证（按提交排序见JSON）；13份输入一致。运行前整卡空闲检查保存在启动日志，归档时GPU无计算进程，未干预其他项目。

## 科学结果分层

|层次|实际证据|状态|
|---|---|---|
|逐泡事件身份|29/29事件按pair、开关与同pair次序一一匹配；旧GPU事件逐条重建一致|EDGE_ID_SEQUENCE_PASS|
|严格事件时间|12/29超过10 µs，最大130 µs|RAW_EDGE_TIMING_FAIL|
|跨泡对全局顺序|一对事件相对顺序反转，时间分离超过原容差|FAIL|
|团级位置|质心0.048867 nm；固定ID10.411515 nm；Kabsch 4.637528 nm|GROUP_POSITION_PASS（限定窗口）|
|团级形状|Rg 0.006359%、长短比0.037649%、张量0.059827%|GROUP_SHAPE_PASS（限定窗口）|
|成员与连通分量|501/1001个含初态接受帧均为同一32-ID分量|GROUP_LEVEL_NUMERICAL_EVIDENCE|
|持续分裂/合并|没有多分量、持续事件或右删失事件|NO_PERSISTENT_SPLIT_OR_MERGE_OBSERVED|
|持续事件预测精度|窗口内没有可检验的持续分裂/合并事件|PERSISTENT_TOPOLOGY_NOT_EVALUATED|
|GPU候选|32泡输出符合CPU/独立参考比较限，但内部停止判据失败|GPU_B2_NOT_ADMITTED|
|物理/生产/AI|未实验校准；没有批量或长窗证据|EXPLORATORY / PHYSICS_UNCALIBRATED；training_ready=false|

历史`MORPHOLOGY_NOT_CONVERGED`保留，团级新增标签不能覆盖旧严格门槛。

## 16个问题的实际回答

1. **130 µs来自哪些因素？**接受坐标独立重建排除了时钟错位。10/5 µs离散路径产生nm级间隙差，ID7–24在阈值附近有浅谷和慢速离开；微小空间差被放大为事件时间差。线性跨越估计仍相差135.447 µs，说明不只是5/10 µs采样量化。固定状态和唯一短容差控制支持“接触停止容差不是当前主要来源”。尚未分离所有声力更新、阻力冻结和时间积分子误差，更不能从两条轨迹断言绝对时间阶数。
2. **ID7–24如何过178 nm？**最后断边时，10 µs轨迹在1.650→1.660 ms由177.999758→178.047093 nm，5 µs在1.785→1.790 ms由177.992241→178.070134 nm。斜率分别4.7335/15.5787 µm/s。两条轨迹各先断、再连、再断，共3次变化；这是有限的近阈值回连，不是大量无规律抖动。这些间隙为正，泡对几何弹簧力为零。
3. **29条能否追溯？**能。`analysis/EVENT_MATCHES.json`包含pair、开关、各步真实时间和跨越括区、前后间隙、速度、阈值差、重复开关及反转索引；`ACCEPTED_GEOMETRY_DT*.npz`保存全496泡对全时序。
4. **阈值变化多大？**原2.100R为29/29、最大130 µs；仅离线改2.095R为28/28、最大30 µs；2.105R为27/27、最大30 µs。三个阈值均未通过原10 µs全部事件门槛；没有选取新“正确阈值”。全程最大团比例均为1。
5. **影响持续拓扑吗？**当前窗口没有影响：32-ID成员集合和唯一分量一直相同。2 ms持续性规则下没有事件可比较，右删失0；这不是对未来分裂/合并预测能力的通过。
6. **求解容差影响多大？**14状态的1e−14 m全部收敛；1e−15 m为13/14，DT10_0141停在1.002579e−15 m并保留失败。全部紧容差末次结果最大位置变化0.000097524 nm、接触力2.4198e−16 N，活动集合/边相同。失败末次迭代不是合格参考。独立方程根残差最大1.039e−22 m。
7. **是否分离了误差？**有界地区分了当前证据：同dt、收紧10倍容差重放2 ms，固定ID最大变化0.005588 nm，23条前缀事件时间全不变；与5 µs轨迹相比仍130 µs。没有将局部残差直接当成全局误差界。
8. **2.5 µs是否运行？**没有，`NOT_RUN`。当前诊断已能解释任务相关敏感性；不继续为消除130 µs无界加精度。本轮唯一短补充已用于同dt容差传播对照。
9. **旧严格验收？**仍FAIL，12/29事件时间超限且存在跨pair顺序失败。
10. **独立团级证据？**有，位置/形状均过原对应工程限、每帧成员相同，单列`GROUP_LEVEL_NUMERICAL_EVIDENCE`；不自动设训练可用。
11. **GPU是否真正开展？**是。Warp双精度、矩阵自由阻力、真实球面接触、块Jacobi-PCG、复用缓冲和两份CG Graph，B=1固定状态；条件启动前A–C分析和`GPU_START_GATE.json`已保存。
12. **误差和耗时？**32泡相对CPU位置2.328e−19 m、力1.072e−20 N；相对独立根位置2.721e−15 m、力1.376e−18 N；实际方程残差2.607e−15 m，集合完全一致。CPU接触调用12.799 ms（不含阻力装配），GPU单次完整尝试35.457 ms（不含首次图构建40.560 ms），其中35次标量同步22.242 ms。没有完成公平的5次重复性能基准，不能报加速倍数。
13. **值得生产集成吗？**当前不值得直接集成。GPU输出接近CPU并不等于求解器通过：保守位移界1.850e−14 m未达1e−14 m，线搜索20次耗尽，基础PCG也没达自身递推残差目标。成本中CPU同步占比高，按协议停止并保留`GPU_B2_NOT_ADMITTED`。
14. **图在哪里？**见下方四组PNG/PDF，全部英文、白底、300 dpi，实际打开PNG与PDF检查。只重画图，不重算轨迹。
15. **进度文档？**对`/home/lzy/projects/microbubble_aggregation.md`最小追加F10I，备份`.bak_F10I_20261012`；before/after、diff、SHA和标题检查在`progress_document/`。
16. **F10J先做什么？**优先用同一固定状态检查GPU基础/外层真实残差和停止语义，避免每次Armijo全主机同步，再做有限公平计时。只有这些通过且端到端成本合理才讨论集成。另行预登记任务级持续拓扑验证；不要把当前“未出现事件”泛化成分裂/合并正确，也不要仅为消除130 µs扩大轨迹扫描。

## 图与证据索引

所有图目录：`results/fast_pzt/f10i/figures/`。

|图|PNG / 同名PDF|内容|
|---|---|---|
|1|`Fig_F10I_1_edge_crossing.png`|间隙、真实接受步、速度、29事件、阈值敏感性|
|2|`Fig_F10I_2_topology_significance.png`|27变化pair边状态、成员、持续性、Rg和长短比|
|3|`Fig_F10I_3_solver_vs_timestep.png`|旧残差、14状态重解、1项非收敛、唯一短控制|
|4|`Fig_F10I_4_gpu_feasibility_admission.png`|实际GPU误差、残差、单次成本、失败准入|

`analysis/FINAL_GATE.json`是分层最终结果；`GPU_FAILURE_ANALYSIS.json`是只读失败复核。`protocol.json`先于离线分析登记，`gpu_protocol.json`先于GPU科学算例登记；四泡负例识别在v1零科学完成后、v2之前记录。`DEVELOPMENT_DECISIONS.json`给出唯一短控制必要性，`GPU_START_GATE.json`给出GPU启动依据。

## 真实资源与复现

离线全事件/形态分析主体5.640 s；14固定状态、3容差级及独立根主体1.448 s（WSL CPU）。唯一短控制远程进程8.626 s，积分/质量2.278 s，详细observer输出4.942 s。GPU v1/v2进程3.313/1.563 s，合计三个进程13.502 s；这不是全研发时间、纯GPU时间或生产吞吐。

原通用计时器在短控制记录内输出`valid_sample_count=1`及诊断前缀每小时率；它只指该短窗记录完成，**不是本轮新增独立科学有效AI样本或完整10 ms样本**。本轮这些正式产出数为0，不据此前缀率推算生产吞吐。

GPU v2进程500 ms抽样观察到498 MiB，短控制482 MiB，v1为418 MiB；不是精确峰值。32泡显式缓冲78,376字节还不含CG私有工作区、CUDA上下文和图。环境记录见`ENVIRONMENT.json`及CPU结果自身环境：服务器Warp 1.18.0、NumPy2.2.6、SciPy1.15.3；WSL参考NumPy2.5.3、SciPy1.18.1。没有升级驱动。仅单独临时工具环境固定安装ruff0.14.0用于语法/风格检查。

本地离线入口（需`PYTHONPATH=src:.`）：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src:. /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.fast_pzt_f10i.verify
```

`events`、`solver`、`propagation_analysis`、`final_analysis`和`figures`分别负责离线阶段。原运行输出已带当时提交及源码SHA，不为刷新身份而重算。`remote`部署Git archive不可变快照、检查GPU并使用独立Supervisor，`timed`保存冷进程计时、退出状态和资源抽样，`archive`逐项验证回传及Git执行快照。最终本地`ARTIFACT_MANIFEST.json`覆盖所有F10I证据（清单自身除外），`DELIVERY_VERIFICATION.json`记录历史不变、结果一致、文档与图完整性。Git同步未作为本轮任务，交付保存在本地独立分支。
