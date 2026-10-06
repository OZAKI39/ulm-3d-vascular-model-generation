# P2/P1 线性代数与资源审查：重新执行

日期：2026-10-06。本轮授权：“重新执行任务”。独立目录 `p2_linear_algebra_rerun_20261006`，上一轮结果和禁止重试标记原样保留。

**P2_PIPE_LINEAR_ALGEBRA_FEASIBLE = NO**

**NEXT_RECOMMENDATION = STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT**

这里的NO限定为：本轮预声明候选在当前资源保护范围内未能完成setup、证明安全可行；不是数学上证明所有PETSc方法都不可行。本轮新候选通过P1回归，完整P2真实系统在CPU PC setup期间被监督器按资源阈值停止。未发生OOM；没有放宽上限或重试P2。

| 问题 | 结果 |
|---|---|
| 223,119个额外storage slots是什么？ | P2边节点上的非数学压力占位，逐层mapping继承验证。 |
| phantom pressure占多少rows/NNZ/内存？ | 21.509%行；compact理论去掉36,344,511项、即32.51045%结构NNZ；单份CSR节省0.4070 GiB，非整机峰值。 |
| 111.79M结构NNZ的主要来源？ | 6,987,082个TET10节点邻接块统一预分配4×4槽，含压力边占位及数值零项。 |
| ASM overlap=2实际重复多少？ | 冻结测量：2,024,020子域行（1.95120倍）、207,653,552子矩阵项（1.85748倍）。 |
| CPU LU fill预计/实测多少？ | P2 LU未运行；继承P1 fill21.2544估算overlap0/1/2的value/index为25.50/31.33/49.33 GiB。不是P2实测。 |
| 本轮最多2个候选是什么？ | 仅1个：GPU主矩阵/向量 + 同矩阵CPU副本上的原生full Schur selfp/GAMG，避免原CUDA PC稀疏乘积。 |
| 候选在P1是否保持同一科学解？ | YES，9次非线性迭代，最大Δu=1.858e−16 m/s、Δp=7.816e−14 Pa。 |
| 完整P2单系统setup实际峰值？ | setup未完成。被停时实测过程树RSS峰值25.217 GiB，是完成setup所需峰值的下界。 |
| 是否保留≥20%安全余量？ | **未通过完成setup后的门槛**。停前瞬时保守余量9.903 GiB大于所需8.411 GiB，但内存仍增长、先触发RSS上限。 |
| 是否实际运行P2 pipe CFD？ | NO：完整CFD=0，线性求解完成数=0。仅真实装配和一次PC setup尝试。 |
| actual P2场是否改善？ | NOT_EVALUABLE，没有实际P2场。 |
| 是否值得下一轮BraVa P2 pilot？ | NOT_EVALUABLE，本轮不支持启动。 |

[三张300dpi图](OPEN_RESULTS.html) · [决策JSON](P2_LINEAR_ALGEBRA_DECISION_GATE.json) · [执行收据](EXECUTION_RECEIPT.json) · [停止原因](P2_STOP_DIAGNOSIS.json)

## 实际改动与回归边界

保留8 MPI、CUDA Amat/outer Vec、FGMRES、right/unpreconditioned norm、rtol=1e−10、atol=1e−24、max_it=2000、restart=100。PC仍为PETSc full Schur，selfp使用已有的A11−A10 diag(A00)⁻¹ A01，两场均用已安装GAMG；无HYPRE或新压力算子。

validation-only适配用 `MatConvert(..., MATMPIAIJ)` 建立CPU Pmat，再用 `KSPSetOperators(A,P)`。每次复制后 `MatEqual` 检查通过。RA仅增加可选Pmat参数，重建判断、最大age5、1.5倍迭代增长逻辑及原RHS恢复逻辑不变。[适配diff](CPU_PC_AND_OBSERVATION.diff)、[RA两处改动](RA_BACKEND_ONLY.diff)和[构建身份](BUILD_IDENTITY.json)可审阅。

该路线是GPU主系统配合CPU预条件器。没有CPU CFD benchmark，没有改变P2/TET10速度、P1/TET4压力、common-curved补丁、rho/mu/Umean、BC、网格、时间步、generalized-alpha、stabilization、非线性容差或探针；没有更换production binary/PETSc/CUDA。P1受控段只从原step10继续到step11。

成熟API依据：[PETSc KSPSetOperators](https://petsc.gitlab.io/petsc/main/manualpages/KSP/KSPSetOperators/)、[PCFieldSplitSetOffDiagUseAmat](https://petsc.org/main/manualpages/PC/PCFieldSplitSetOffDiagUseAmat/)。实际执行语义以本机3.25.5源码为准，[源码摘录](CPU_PC_NATIVE_SOURCE_EVIDENCE.txt)与SHA记录对应本机库。

## P1同检查点验证

原baseline复用冻结证据，先重新验证输入/源码/库身份。本轮新候选声明在运行前写入 [candidate plan](LINEAR_SOLVER_CANDIDATE_PLAN.json)。CSR结构与RHS逐位相同；矩阵最大绝对差1.05879e−22，逐行尺度最大差3.34971e−16（1.509个double epsilon），低于预声明128 epsilon门槛。字节SHA不同如实保留；没有放宽科学容差。

| 指标 | 原baseline | 本轮CPU Schur |
|---|---:|---:|
| 同一step10→11耗时/s |249.056|124.416|
| 自身任务RSS峰值/GiB |9.251|5.562|
| 总Krylov迭代 |242|2264|
| 非线性迭代 |9|9|
| Q/m³·s⁻¹ |1.1175570752082609e−7|1.1175570752082611e−7|
| 中央拟合dp/dz /Pa·m⁻¹ |−9050.309514439|−9050.309514439004|

相对状态差3.420e−16，Q和dp/dz相对差约1.18e−16和6.03e−16。9次CPU Pmat检查通过，各级实际预条件矩阵均CPU AIJ，GPU主矩阵和向量身份保持。速度与压力逐项比较见 [parity JSON](P1_LINEAR_SOLVER_PARITY.json)。这证明受控段解等价，不是稳态精度认证或一般性能保证；Krylov次数明显增加，仍在预声明10倍baseline门内。

首次P1尝试在iteration-0监测中停止：检查器误把pressure finest-level隐式Schur operator要求为AIJ；已改为检查其真实Pmat，算法和容差未变。首次日志、[修正diff](OBSERVER_CORRECTION.diff)及所有尝试收据保留。准备阶段也修复过JSON列表继承类型错误，无科学计算影响。

## P2真实系统的资源结果

本轮实装配保持1,037,320行、111,793,312结构项；精确非零61,418,083，结构存在但数值零50,375,229。223,119个dummy pressure rows仍只有正对角，RHS/非对角为0。与上轮精确非零项相差3，不据此声称两个P2矩阵逐位相同；本轮A→P复制的MatEqual一致性已直接检查。

完整节点布尔图G²实际计数97,504,968，为selfp pressure product的结构上界；原生日志也实际出现相同97,504,968个结构项的压力矩阵（相对原压力对角块6,987,082项达到13.96倍），对应CSR payload约1.091 GiB。[原生日志计数](P2_PRESSURE_MATRIX_OBSERVATIONS.json)中的重复打印不能当成多个同时存活副本累加。[图预算](CPU_SCHUR_GRAPH_BUDGET.json)明确它不约束所有MPI通信副本、symbolic临时数组、AMG或全进程峰值。原17.932 GiB任务预估被实测突破，不能用于下一次放行。

| 实测项 | 本轮值 |
|---|---:|
| cgroup总上限 |42.055664 GiB|
| P2启动时保守非回收占用（含其他项目） |10.887859 GiB|
| 自身任务RSS保护上限 |24 GiB|
| setup中断时采样自身任务峰值 |25.217113 GiB|
| cgroup raw峰值（含文件缓存） |42.055664 GiB|
| 保守非回收内存峰值 |32.152225 GiB|
| 对应保守余量 |9.903439 GiB|
| 20%所需余量 |8.411133 GiB|
| 末20样本非回收内存增长拟合 |+0.292864 GiB/s|
| 最小采样GPU free |12,293.8125 MiB|
| 主机OOM / OOM kill增量 |0 / 0|
| CUDA allocation failure / admission denial |0 / 0|
| setup监督耗时 |51.253 s|

监督器在50.245秒采样触发 `PROCESS_TREE_RSS_LIMIT`，约51.253秒完成退出。**不是“余量已经小于20%”，也不是“GPU和CPU同时OOM”。** 余量门须在setup/solve完成且不持续增长时评价，本轮这些前提不成立。

日志出现velocity hierarchy至6×6矩阵，随后是259,330×259,330压力乘积/子矩阵摘要；没有进入after_PC_setup回调或P2线性迭代。可定位到首个CPU PC setup阶段，但没有捕获停止瞬间的确切分配函数，不能断言某个具体AMG层或LU是唯一原因。各rank的getrusage只记录到setup前，因此图中的rank柱为此前高水位下界，不能当成完整rank峰值。

CPU副本本身的CSR约1.253 GiB，并不足以解释整个增长；MPI分布稀疏乘积的中间结构、重复数据和多重网格对象仍需内存。P2实测过程树峰值已比规划量高7.285 GiB且持续增长。没有为获取更大峰值而重新运行。

本轮CUDA分配保护拦截cudaMalloc/cudaMallocAsync/cudaMallocManaged，在共享锁内先查free并留2 GiB，4 KiB正常分配及4 TiB请求“分配前拒绝”smoke均通过。它不是全部driver内部路径的数学保证。本次P2捕获92次允许的调用，最大约1.05 MiB；还未进入完整Krylov工作区，不据此推断后续GPU峰值。GPU每秒NVML v1采样的used含device reservation，与上一轮nvidia-smi used口径不同，跨轮比较用free或明确口径。

## 冻结证据、compact与科学停止门

A/B/C存储、native preallocation和ASM overlap审计复用带SHA的冻结证据，不重复这些科学前置任务。本轮P2再次实际构造同一原生图，旧图/本轮读数来源分别标记；`MatGetInfo.memory=0`仍表示API不可得，而非矩阵不占内存。结构计数推导的字节数不冒充RSS或PC内部独占分配。

compact mixed graph保留为reference；单份CSR可省0.4070 GiB，但没有可靠归因证明完整PC setup峰值改善≥15%。因此不重构编号、不删除数学DOF。[compact gate](COMPACT_PROTOTYPE_GATE.json)。P2 local LU仍未执行，symbolic例程会分配数值数组，不能用无边界symbolic试跑冒险。

上一轮独立复现的 `Integrator::corrector_taylor_hood` 曲面边压力回填异常仍为FAIL：用四角点仿射TET4几何反解曲面边点。本轮仅允许线性代数改动，未修复或重跑此科学/输出路径。[原生输出gate](HIGH_ORDER_OUTPUT_GATE.json)。即使内存gate通过，完整CFD也不能跳过此门。

完整P2 CFD、BraVa CFD、trajectory、population、dt sweep均为0，三个P2科学结果JSON全部NOT_RUN；不生成04/05科学精度图、不推断actual P2误差改善。旧1565条结果仍为458/527/437出流、141截止、2失败；不改成新灌注协议结果。

资源限制下仅保留一个下一方向：**更高内存服务器上的受控验证**。当前证据不支持指定一个“保证够用”的内存数值，也不授权自动换机或重跑；曲面输出gate仍需独立处理后才能进入完整CFD。

## 交付与检查

图01区分全cgroup、过程树RSS、计数推导payload及rank下界；图02保留实测overlap与估计LU fill的区别；图03比较同段P1回归和中断P2资源，不能用其推断物理场精度。三图均300dpi，[图数据](FIGURE_DATA.json)可复查。

[TEST_RESULTS.json](TEST_RESULTS.json)记录本轮工程检查，[FINAL_INTEGRITY.json](FINAL_INTEGRITY.json)在收尾再次核对保护输入、2431个生产源码文件、binary与PETSc；[执行收据](EXECUTION_RECEIPT.json)含11个已退出监督作业（累计约200.9秒作业时间，不是全部开发耗时）。所有自身TASK移除；未触碰其他项目GPU进程、兼容链接或系统设置。仅本轮目录和三个入口文件独立commit，不push。
