# Stage SV0 — SimVascular 网格与有限元血流求解可行性验证

**STAGE SV0 STATUS: BLOCKED — ENVIRONMENT**

**可行性分类：SOLVER_ENVIRONMENT_FAIL；具体原因：CONTAINER_RUNTIME_UNAVAILABLE。**

本次执行于 2026-09-18，完成 Step 1–4，在 Step 5 的正式 solver container 环境门槛停止。
官方 SimVascular 无头 Python API 已运行成功；实际血管导入、重新网格与 svMultiPhysics 求解均未执行。
这里的 `SOLVER_ENVIRONMENT_FAIL` 表示容器环境不可用，不能解释为已经运行的 official smoke 或血流方程求解失败。
没有达到 CONDITIONAL PASS，也没有获得 production CFD 推荐依据。

状态证据：[final_status.json](final_status.json)、[execution_status.json](execution_status.json)。
WSL 新工程是 source of truth；远端探测与只读环境审计已 fetch 回本工程。

## 为什么尝试 SimVascular？

旧 Stage 3 已定位到具体的 P2/P1 网格与函数空间相容性问题：两个压力自由度在壁面速度约束后失去有效耦合。
本阶段建立独立工具链，保留这一既有诊断，评估成熟工具重新三角化表面并重建体网格后的可行性。
旧 DOLFINx solver 没有被修改、复制或求解。

旧 FEM 本机审计覆盖 **7440 项**，包括 source、config、tests、reports、outputs、logs、环境和 `.git`；
文件 SHA256、大小、模式、mtime 与符号链接均核对，读取产生的 atime 不作为变更。
最终新增、删除、修改均为 0，HEAD 保持 `e41cc9c56852fcf88fe4640b2b598a08a082370d`。
既有用户改动 `reports/stage01/REPORT.md` 的字节与 Git diff 指纹保持不变。
见 [old_fem_readonly_audit.json](old_fem_readonly_audit.json) 和保留的 baseline/final 全量文件指纹。
远端另核对旧 Python 可执行文件及 Conda package metadata 共 198 项，指纹一致；该检查的范围不等同于远端全目录审计。
见 [old_remote_environment_audit.json](old_remote_environment_audit.json)。

## 使用了哪个血管模型？

正式 reference 是 Stage 1.7 selected exterior surface。

| 输入 | 原始路径 | SHA256 |
| --- | --- | --- |
| Exterior surface | `/home/lzy/projects/formal_3D_flow_solver/FEM/outputs/stage01_7/selected/surface/tagged_surface_si.npz` | `31746b2c2ef95a0b9044d96ef4b470a83592a91ae9132c5d6810736880734ac6` |
| Selected volume | `/home/lzy/projects/formal_3D_flow_solver/FEM/outputs/stage01_7/selected/mesh/volume_mesh.npz` | `fc036eb0eaa0b365b06a5033f5f9eb31addfe464af4beb53c4cd96a25ed594a3` |
| Frozen physical config | `/home/lzy/projects/formal_3D_flow_solver/FEM/configs/stage03_reference_vascular.yaml` | `ca7b9037295f05f50a36e6b46bee0a608ab01b6a5025b348ffd4e69e19db9539` |

总共复制 9 个最小正式输入，合计 5,077,324 bytes；每个 source SHA 均等于 copied SHA。
源路径、逐文件 SHA、大小、复制时间与旧 Git 状态见
[MANIFEST.json](../../inputs/fem_reference/MANIFEST.json)。没有复制整个 FEM 项目。
源表面当前为保留标签的 NPZ；尚未进入 Step 8 的正式 PolyData/VTP 适配。

| Original semantic role | Original ID | SimVascular face ID |
| --- | ---: | --- |
| WALL | 1 | 未导入，未建立 |
| OUTLET_03 | 2 | 未导入，未建立 |
| OUTLET_01 | 3 | 未导入，未建立 |
| INLET | 4 | 未导入，未建立 |
| OUTLET_02 | 5 | 未导入，未建立 |

没有使用 Stage 1 fan-cap 或 Stage 1.8 MMG/TetGen 候选。
`source_geometry_and_faces.png` 位于 Step 7，因 Step 5 已阻塞而未生成。

## SimVascular 环境是否真的可用？

**已证明的范围是官方 headless Python API；实际 meshing 能力仍未测试。**
从[官方 GitHub release](https://github.com/SimVascular/SimVascular/releases/tag/2023-05)
下载 Ubuntu 20 Linux distribution，包版本 **2023.05.31**。
探测时官方 GitHub 的更新 release 未提供 Linux binary asset；官方网页指向的 SimTK 下载页返回 HTTP 403。
因此选择这一可公开获取的官方 Linux binary，而非宣称它是最新 SimVascular 软件版本。
仅在新项目 `outputs/sv0/simvascular_distribution/` 解包，没有修改系统安装或旧 FEM 环境。

发行包 SHA256：`38ad92410629042246394b5eb4ee703a3062a15e130b54884ef0a20081a3cff1`，大小 840,795,966 bytes。
精确 URL 与选择依据见 [simvascular_distribution.json](simvascular_distribution.json)。
官方 launcher 实际执行 `simvascular -python -- scripts/sv_api_probe.py`，退出码 0；
嵌入 Python **3.5.5**、VTK **8.1.1**，成功创建 `sv.meshing.TetGen` 和 `sv.modeling.PolyData`。
`surface_mesh_flag=true`、`volume_mesh_flag=true`、`use_mmg=false`；构造探测没有加载血管，没有执行 `generate_mesh`。
部分未使用的 MITK DICOM 模块提示缺少 `libicuuc.so.66`，但上述 API 探测确已成功。
官方包自带内部 meshing libraries；没有另外安装或执行独立 MMG/TetGen。
见 [simvascular_api_probe.json](simvascular_api_probe.json) 与
[完整 stdout/stderr](evidence/simvascular_headless_attempt_01.json)。

| Host | OS / architecture | 初始 Python | 容器实测 |
| --- | --- | --- | --- |
| WSL | Ubuntu 24.04.2 / x86_64 | 3.13.11 | Docker Desktop bridge 存在，但 version/info/pull 均退出 1，要求启用 WSL integration；无 Docker socket |
| Remote | Ubuntu 24.04.4 / x86_64 | 3.12.3 | Docker/Podman/Apptainer/Singularity 均未安装；`unshare -Ur` 返回 EPERM，无 Docker socket 或 `/dev/fuse` |

远端实例说明明确其运行于非特权 Docker 容器，不能再运行嵌套容器引擎。
该结论仅针对当前实例，不外推为所有远端机器均不可用。见
[environment_probe.json](environment_probe.json)、[实例说明](evidence/instance_guide.txt)。

按[官方 svMultiPhysics README](https://github.com/SimVascular/svMultiPhysics)选择 `simvascular/solver`。
已实际尝试两端 `docker pull simvascular/solver:latest`：WSL 返回 1，远端 Docker executable 不存在。
另外读取并 SHA 校验了官方公共 registry 的 manifest/config，解析到 linux/amd64 不可变引用：

```text
simvascular/solver@sha256:87928224ff90506d8b84553cf201a8f5a163baa1c6dc5b4c5cb5ccad405c12eb
```

Registry 创建时间：`2026-09-17T21:26:12.535606993Z`。
**这只是 registry metadata，未拉取 image layers，没有 local image ID、docker inspect 或已运行的 container。**
本次没有任何容器挂载。以后实际 pull/run 应使用此不可变引用，不能把解析 digest 当成环境成功。
原始 metadata 和两端失败输出见 [solver_container_manifest.json](solver_container_manifest.json)。

官方 `fluid/pipe_RCR_3d` smoke：**NOT_EXECUTED_ENVIRONMENT_BLOCKED**。
尚未选择兼容 revision、获取 LFS case 或执行 solver；exit code、MPI ranks、运行时间、VTU SHA、finite u/p、官方 expectation 全部没有实际值。
见 [official_smoke.json](official_smoke.json)。因此没有进入真实模型调试。

## SimVascular 如何重新生成网格？

尚未重新生成。用户指定流程允许重新排列表面三角形，同时用几何距离、体积和端口检查约束真实血管形状。
预定通过正式 `sv.meshing.TetGen` 启用 surface 与 volume meshing。
按顺序，必须先通过官方 smoke，再完成几何导入，自动计算 Stage 1.7 唯一体网格边长中位数 `h_ref`、wall 尺度及端口半径，冻结 policy 后生成 M0。
只有 M0 生成失败或 hard geometry/validity 失败，才允许 `0.80*h_ref` 的 M1。
本次 `meshing_policy.json`、M0、M1、face map 均未产生；未做缩放、边界层、各向异性或参数 sweep。

## 血管几何改变了吗？

没有新候选可比较，结论是 **NOT_EVALUATED**，不能把“未 remesh”报告为 geometry PASS。
未来检查须包含双向 wall 距离 max/P95/RMS、体积相对差，以及四端口的面积、中心、法向和尺度相关平面误差。
`surface_before_after.png` 与 `surface_geometry_error.png` 均未生成。

## 新网格质量怎样？

以下仅引用复制的 [Stage 1.7 baseline quality](../../inputs/fem_reference/baseline_quality.json)，没有在 SV0 重算旧质量，也没有新 SimVascular 结果：

| 指标 | Stage 1.7 baseline | SimVascular |
| --- | ---: | --- |
| vertices | 43178 | 未生成 |
| tetra | 147948 | 未生成 |
| q_min / Gmsh minSICN | 0.06238017030953989 | 未评估 |
| P1 | 0.3513432334394028 | 未评估 |
| P5 | 0.4856705783671268 | 未评估 |
| median | 0.7218556061894749 | 未评估 |
| P95 | 0.9293127193055747 | 未评估 |
| N(minSICN < 0.1) | 3 | 未评估 |
| low-quality fraction | 2.02773947603e-05 | 未评估 |

未来的 invalid tetra hard gate 与 baseline-relative production quality gate 必须分开。
本次无 `simvascular_mesh_cutaway.png` 或 `mesh_quality_comparison.png`，不能推荐新网格质量。

## 原来两个 pressure island 还在吗？

旧 Stage 3 frozen diagnosis：unsupported pressure DOF = **2**。
新 SimVascular mesh：**未生成，实际数量未知，不是 0**。
旧 audit 的只读调用接口已保留，但实际新网格 compatibility audit 尚未执行；没有组装或 PDE solve。
见 [baseline_pressure_support.json](baseline_pressure_support.json)。
这是 P2 velocity / P1 pressure / wall no-slip 下的 DOLFINx compatibility diagnostic，
不是 svMultiPhysics 自身的必要条件，也不能替代真实 solver PASS。
`pressure_support_comparison.png` 未生成。

## SimVascular 求解的边界条件是什么？

**以下是用户指定的后续目标，不是已经运行的 case。**
Wall no-slip；三个出口均为 zero-traction Neumann；入口为明确指定的 plug velocity profile，按实际离散入口面积归一化到 `Q_reference`。
写入的速度场还必须独立积分，满足 `|Q_written-Q_reference|/Q_reference <= 1e-12`。
真实模型不允许 Resistance 或 RCR outlet；官方 smoke 自身的 RCR case 与真实模型目标不同。

指定入口速度剖面不等价于旧 FEM 的 integral-Q-only condition。
没有创建 real solver XML 或 `sv_real_reference.yaml`；当前 image 内 binary 的 equation mode、BC syntax 和 pressure reference 尚未核验。
不能声称已实现 steady Stokes。若后续只能用恒定边界的 transient incompressible Navier–Stokes，必须如实标记 transient-to-steady。
Zero traction 也不能被解释为所有出口节点逐点 `p=0`。

## 用的 Q 是实验流量吗？

**NO. REFERENCE NUMERICAL CONDITION / NOT EXPERIMENTAL PUMP FLOW。**
物理值来自 Stage 3 frozen YAML，并与其记录的 JSON 执行快照核对；没有读取当前 LBM 配置替代 frozen reference。
新 [reference_condition.yaml](../../configs/reference_condition.yaml) 明确 `experimental: false`。

| Quantity | Frozen reference | Unit |
| --- | ---: | --- |
| rho | 1056 | kg/m³ |
| nu | 3.27e-06 | m²/s |
| mu = rho × nu | 0.0034531200000000001 | Pa·s |
| Q_reference | 2.7369132390905703e-15 | m³/s |

整个参考条件为 SI：m、kg、s、m/s、Pa、kg/m³、Pa·s、m³/s。
没有 SI→mm adapter 或实际新 mesh 的 A、Dh、Re；新网格 Reynolds number 不能报告为已验证。

## solver 是否稳定运行？

**未运行。** Official smoke 和 SV_REAL_PLUG 均无 solver 进程、MPI ranks、timesteps、convergence history 或 final state。
真实模型默认 4 MPI ranks 只是用户规定的未来设置，实际 MPI ranks 为 N/A。
Mesh generation time/RSS、solver elapsed/peak memory、result VTU size 均 N/A；真实 result count = 0；GPU used = false。
无 `solver_resource_usage.png`。

## 是否达到稳态？

**NOT_EVALUATED**。没有保存状态或收敛图。
未来标准是最后连续 5 个保存区间 `E_u <= 1e-5` 且 `E_Q <= 1e-6`，因此需要至少 6 个状态。
本次只对验收函数做合成数据测试，没有真实场的 L2 时序积分。
未计算或冻结 dt，也没有第一 block 或 restart extension。
`steady_convergence.png` 未生成。

## 流量守恒吗？

**NOT_EVALUATED**。没有 B_INT，也没有真实结果 VTU/VTP 的独立积分。
Q_target 是上述参考值；Q_in、Q_out_total、epsilon_Q、epsilon_mass 都是未知，不写成 0。
未来四端口的两套积分相对差与实际 inlet/mass 误差均须满足各自 `1e-6` 门槛。
`flux_balance.png` 未生成。

## 三个出口分别流多少？

OUTLET_01、OUTLET_02、OUTLET_03 的 signed outward flow、flow fraction 和面积平均压力全部 **N/A**。
没有预先分配流量或从参考输入填充假结果。
`outlet_flow_split.png` 未生成。

## 速度和压力场是什么样？

没有真实解，因此 velocity/pressure finite、pressure range、velocity norm 都未评估。
`velocity_global.png`、`velocity_slices.png`、`pressure_global.png`、`pressure_sections.png` 均未生成。
Synthetic VTU 新进程读写测试通过只能说明初步 parser/I/O 工作，不能当成实际 solver solution reload PASS。
实际结果 reload：**NOT_EXECUTED**。

## inlet profile 会影响结果吗？

本次 **NOT_AVAILABLE**：plug 尚未运行，也未进入官方 parabolic preprocessing capability 检查。
同网格 plug/parabolic 的速度/压力 L2 差异、出口比例、面积平均压力及距离入口大于 `2*Dh` 的区域敏感性全部未知。
没有 `profile_sensitivity.png`；其条件性生成要求没有满足。

## 和旧 DOLFINx 能直接比较吗？

当前不能宣称严格同 BC 比较：两者入口数学条件不同，并且没有新流场。
已核对的共同参考是来源几何、frozen 物性和总 Q；本次不能比较实际 gross flow behavior。
正式 FEM-vs-SimVascular same-BC comparison 需要另立阶段，本次没有自动开始。

## 是否值得继续使用 SimVascular？

**分类：SOLVER_ENVIRONMENT_FAIL。阶段状态：BLOCKED / ENVIRONMENT。**
Headless API 可用是有限的环境正面证据，尚不足以支持 SV_FEASIBILITY_ONLY 或 SV_PRODUCTION_CANDIDATE。
继续的具体前提是可运行官方 OCI/Singularity solver image 的远端环境；然后从已解析 digest 实际 pull、inspect，
使用兼容官方 revision 和所需 LFS 数据完成 fluid smoke，再按原顺序恢复 Step 7。
本次没有绕过容器使用 standalone solver，也没有更换或购买远端实例。

软件与证据测试：**48 passed、0 failed、12 skipped**，保留全部 23 个要求的测试模块。
10 类指定恶意 regression 均覆盖，包括角色交换/缺失、几何漂移、混合单位、入口未归一化、错误 Resistance、NaN、退出 0 但缺少 VTU、质量守恒失败、旧 source 改动检测。
对旧 source 的恶意修改仅在临时 fixture 上模拟，没有修改旧 FEM。
另验证 NaN outlet 不能被聚合隐藏，且修改状态 JSON 不能让尚未实现的真实验收假通过。
真实科学验收的 12 项 skip 明确由环境阻塞触发；后续这些集成断言仍需实现，不能将初步验收库当成完整 solver pipeline。
见 [test_summary.json](test_summary.json)、[pytest_results.xml](pytest_results.xml)。

## 尚未证明什么？

本次未证明真实血管 import、健康 tetra remesh、pressure-support 消除、边界功能、真实 solver、稳态、守恒、场有限性或真实结果 reload。
没有 experimental pump Q、formal same-BC comparison、mesh convergence、WSS validation、RBC、microbubble 或 particle interpolation。
16 张要求的图片均未生成：流程止于 Step 5，连 Step 7 的源几何图也未执行。
[visualization_manifest.json](visualization_manifest.json) 逐项说明缺失；没有用状态卡片、合成图或旧图冒充本次科学可视化。
没有待用户审核的真实 SV0 场图，因此不申请 CONDITIONAL PASS 的图像审核。

本次阻塞交付已保留 reference inputs、全量本机只读指纹、真实环境 stdout/stderr、registry metadata、远端回传、依赖锁与测试。
关键原始记录有逐文件 SHA：[evidence_manifest.json](evidence_manifest.json)；
终端汇总见 [terminal_summary.txt](terminal_summary.txt)。此处停止，不启动下一阶段。
