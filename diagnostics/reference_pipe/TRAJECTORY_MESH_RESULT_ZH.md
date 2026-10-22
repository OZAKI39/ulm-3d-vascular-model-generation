# 唯一 R1 trajectory mesh 结果

**STOPPED_TRAJECTORY_MESH_RUNTIME_LIMIT。没有最终网格，不能进入 manufactured CFD 或 particle validation。** 本轮一个 R1 尝试、零 correction；未启动其他 mesher candidate。停止原因：`TIME_LIMIT_1800S`。

## 实际运行

| 项目 | 实测 |
|---|---:|
| 原生运行 wall time | 1800.715 s |
| 最后采样累计原生 CPU time | 8371.81 CPU·s |
| 最后采样 CPU time / elapsed | 4.65 cores |
| 最大原生线程数 | 7 |
| 原生 RSS 峰值 | 18.071 GiB |
| 受管进程树 RSS 峰值 | 17.860 GiB |
| 不可回收内存保守估计峰值 | 21.021 GiB |
| 新增 OOM / OOM kill | 0 / 0 |
| 输入 tetra | 184,712 |
| 日志中最大 tetra | 26,001,463（输入的 140.77 倍） |
| 日志中最后 tetra | 19,974,733 |
| 完成的目标 metric 步数 | 0 |
| output.msh 文件状态 | 不存在 |

原生初始记录：`NATIVE_INITIAL vertices=37625 tetra=184712 quality=0.00602143 max_length=14.3465 omp_max_threads=7`。默认准入条件为 native quality≥0.2 且 metric edge length≤2.82843；当前 P0 在这一原生 implied metric 下不能直接进入 `approach_metric`，薄适配器调用官方 `fix`。它先从输入单元得到 implied metric，进行原生质量修复。日志中该阶段开始时 quality 最小值约 0.00140746、最大 metric edge length 74.0853，随后大量 refine/coarsen。

本次膨胀发生在目标25 µm metric正式推进之前，不能归因为该目标尺寸本身导致全域加密；也不能据此判断一个已成功生成的25 µm网格能否通过u/G/s/curl。现有证据仅把本次执行瓶颈定位到默认初始质量修复及其后续几何位移处理，没有继续修改工具内部算法追究更多机制。

以上是当前默认集成的执行路径和观测结果，不是 P0 物理体积非法的证明，也不是 Omega_h 普遍能力失败的证明。日志中的 `[0.00,1.00]` 经四舍五入，不等于检出了零/负体积。目标尺寸约25 µm也不是实际壁邻单元高度。没有为获得外观良好的网格而修改输出、放宽 gate 或重试参数。

CPU 路径实际使用 OpenMP；本作业没有 CUDA 分配，原生 GPU VRAM 为0。资源记录保留 device-wide `nvidia-smi` 观察值，不能把其他进程或显示服务占用记作本作业 GPU 利用。CPU累计值止于最后一次活进程采样；受管时间包括停止和收尾。没有因达到20分钟目标预算而自动增加参数搜索，硬上限为1800s。

## Mesh gates

| 验收项 | 结果 |
|---|---|
| 最终壁面/端面几何、中央 finite-facet deviation≤0.05 µm | NOT_EVALUABLE_NO_FINAL_MESH |
| WALL/INLET/OUTLET、闭合/流体连通/非流形、体积合法 | NOT_EVALUABLE_NO_FINAL_MESH |
| all-WALL tetra=0、明显病态单元 | NOT_EVALUABLE_NO_FINAL_MESH |
| 实际 wall altitude median20–30 µm及P5/P95 | NOT_EVALUABLE_NO_FINAL_MESH |
| bulk/P0=0.75–1.25 | NOT_EVALUABLE_NO_FINAL_MESH |

能力检查中的坐标/标签往返和解析点投影 PASS 不能填入此表冒充最终 mesh PASS。没有合法最终网格，也没有明确的 metric 参数→实际 wall altitude 单调关系，因此不满足一次 deterministic correction 的准入条件。本轮停止，不自动启动 BraVa、其他工具或旧任务。

## 图与证据

[HTML入口](OPEN_RESULTS.html) · （历史图件已退休；结论保留于本报告/JSON，tracked文件可从清理前算例Git提交91ebc021恢复。） · PDF（历史引用 `figures/trajectory_mesh_execution.pdf`；结论见本报告和最终摘要，旧路径已清理） · [结构化结果](TRAJECTORY_MESH_RESULT.json) · 冻结计划（历史引用 `TRAJECTORY_MESH_PLAN.json`；结论见本报告和最终摘要，旧路径已清理） · 原生日志（历史引用 `manufactured_pipe/R1_TRAJECTORY_MESH/mesh/native.log`；结论见本报告和最终摘要，旧路径已清理） · 逐秒资源（历史引用 `manufactured_pipe/R1_TRAJECTORY_MESH/mesh/resources.jsonl`；结论见本报告和最终摘要，旧路径已清理）。

图中 CPU/RSS 横轴是真实秒数；单元数横轴是原生日志统计序号，没有伪造每条日志的时间戳。native metric quality 仅用于解释原生算法行为，不能代替体积/几何/流场 QA。没有生成伪造的新流场或粒子轨迹图。

## 后端预检的最终结论

## MeshSim

现有 `sv.meshing.MeshSim()` 报告 interface is not defined；`create_mesher(Kernel.MESHSIM)` 也不能创建内核。安装中找到 Python 配置定义不等于找到可工作的 MeshSim kernel。`license_status=NOT_ESTABLISHED`，不能把接口缺失解释成已证实 license 无效。volume、boundary layer 和 boundary IDs 均没有可用实例供实际验证。

检查在 10 分钟预算内结束，没有申请 license、安装新的 SimVascular 或修改 production runtime。原始检查（历史引用 `TRAJECTORY_MESHSIM_CAPABILITY.json`；结论见本报告和最终摘要，旧路径已清理）及 SimVascular 启动日志（历史引用 `TRAJECTORY_MESHSIM_PROBE.log`；结论见本报告和最终摘要，旧路径已清理）保留。JSON 中极小的 `elapsed_s` 仅为 Python API 调用部分，不包括应用启动和环境检查。

## NASA/refine

固定源码提交 `1563f71f224d59c8712e2181369ff93440957b0f`，Release、7 个并行编译 workers。原生 tensor `.solb` 与 MPI/EGADS 接口来自 [NASA/refine 源码](https://github.com/nasa/refine/tree/1563f71f224d59c8712e2181369ff93440957b0f)。几何依赖使用 [MIT 官方 ESP 预编译发行包](https://acdl.mit.edu/ESP/PreBuilts/linREADME.txt)，隔离提取到项目目录；没有系统升级，也没有改 production PETSc/SimVascular。

| 能力 | 证据与边界 |
|---|---|
| 读取当前 tetra | 37,625 vertices、184,712 tetra；二进制转换坐标、连接、标签往返精确 |
| anisotropic metric | 原生 6 分量对称 tensor API 已核对；没有实际 R1 adaptation |
| 并行 | 使用既有 production MPI 4.1.6，8 rank help 成功；系统 MPI 启动挂起后仅清理本轮相关进程 |
| 明确 geometry | 自动关联后原生 volume validation 失败，未导出有效关联网格 |
| 边界角色 | 输入 CAD face→WALL/INLET/OUTLET 映射已建立；没有最终网格的标签保证 |

这里的失败属于**当前集成路径的关联预检**，不能推广为 refine 普遍不支持圆柱或 anisotropic adaptation。SI 与毫米两个单位的关联尝试都出现相邻节点重合；没有调整 R1 metric，也没有恢复 MMG。直接投影探针中的 7.32 µm 是坐标位移，不能当成径向几何误差。Omega_h 后续探针也证实了“点位移”和“投影后径向误差”不同。

检查从 17:13:43 UTC 开始，最后一个探针约 17:28:56.8 UTC 结束，约 913.8 s，**超过用户 15 分钟要求约 13.8 s**；决策记录写于 17:30:19 UTC。不宣称严格满足这一检查时间限制。构建记录还保留一次系统 MPI help 挂起，成功构建与 help 启动耗时分开解读。完整能力记录（历史引用 `TRAJECTORY_REFINE_CAPABILITY.json`；结论见本报告和最终摘要，旧路径已清理）。

`CYLINDER.json` 是早期错误抽取了日志末行的原始失败记录，不作为 JSON 几何证据；`CYLINDER_MM.json` 的旧字段名带 `_m`，数值实际为 mm，不能直接当 SI 使用。这两份原始记录未重写。最终 Omega_h 几何记录明确写 `model_units=mm`。

## Omega_h

固定 [Omega_h 9.34.13 源码](https://github.com/sandialabs/omega_h/tree/dc235450cf15b2b9549fabe8472269e7f8618cae)，提交 `dc235450cf15b2b9549fabe8472269e7f8618cae`。独立 Release/OpenMP/EGADS build；没有安装、调用 Gmsh mesher，`.msh` 只是由 Omega_h 自带解析器读写的交换格式。

构建进行了两项小修复：EGADS 调用点的旧 `filter_graph` 名称改接现有原生 `filter_graph_edges`；链接时补齐 `-fopenmp`。前者保留逐字 diff 和源文件 SHA，后者保留 CMake 命令。没有改 split/collapse/swap/warp 等算法。安装没有导出 `Omega_h_egads.hpp`，薄适配器编译时显式引用同一固定版本的源码头文件。

圆柱采用 EGADS 原生 `-CYLINDER`，不人为拆成两个半圆柱面；物理 R、L 不变。原生投影/边界关联使用 [Omega_h EGADS 实现](https://github.com/sandialabs/omega_h/blob/dc235450cf15b2b9549fabe8472269e7f8618cae/src/Omega_h_egads.cpp)，适配调用沿用官方 cylinder 示例和 native quality repair。该几何路径在 host 上执行，本轮选 7 线程 OpenMP，未尝试为 GPU 重写 mesher，也未把 CUDA 的可配置性当成已验证的端到端 GPU 几何支持。

原生预检实际结果：

- P0 坐标逐值保持，tetra incidence 和带标签 boundary incidence 往返一致。
- CAD face 1=WALL、3=INLET、2=OUTLET；对应 30,720／8,192／8,192 个边界三角形。
- OpenMP 实际团队大小 7。读取、分类、投影、导出的进程约 0.432 s；adaptation 调用数为 0。
- 原生投影后的 wall radius 最大误差 `1.11e-13 µm`，cap plane 最大误差 0；最大点位移 `14.6468 µm`，两者不可混用。
- 诊断性地把完整 warp 一次加到所有节点，会得到 174 个非正体积 tetra。这个数组不是 native adaptation 的输出；Omega_h 的原生 `warp_to_limit` 会限制位移并检查质量。没有把这个未经限制的数组送给 CFD，也没有据此宣称最终网格合法。

证据：能力记录（历史引用 `TRAJECTORY_OMEGA_CAPABILITY.json`；结论见本报告和最终摘要，旧路径已清理）、原生预检 QA（历史引用 `manufactured_pipe/trajectory_backend_capability/OMEGA_NATIVE_QA.json`；结论见本报告和最终摘要，旧路径已清理）、唯一冻结计划（历史引用 `TRAJECTORY_MESH_PLAN.json`；结论见本报告和最终摘要，旧路径已清理）。预检没有评估最终 finite-facet deviation、真实 wall altitude 或 bulk scale，必须由真正完成的候选独立验收。
