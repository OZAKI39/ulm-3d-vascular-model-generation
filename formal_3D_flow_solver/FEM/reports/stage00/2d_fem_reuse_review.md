# 可以基本复用的设计

本次只阅读，没有导入、运行或复制旧求解器。以下路径均相对于只读工程
`/home/lzy/projects/ulm_microbubble_traj_gen_2D`；逐文件 SHA256 见同目录
`reuse_provenance.json`。新代码为独立实现。

| 检查内容 | 来源 | 可以保留的原则 |
|---|---|---|
| Taylor–Hood | `utils/flow/dolfinx_gmsh_solver.py:897`、`:1087` | 使用稳定的 P2 速度/P1 压力配对；黏性项使用对称速度梯度，保留不可压缩耦合。3D 要使用三分量速度和动力黏度 Pa·s。 |
| PETSc 求解和诊断 | 同文件 `:1247` | 将组装、求解、导出分开；小算例用直接解作对照；保留收敛原因、迭代次数、实际求解器配置。 |
| 按边界积分流量 | 同文件 `:1178` | 在真实 tagged facets 上积分 `u·n`，随后做 MPI SUM；入口外向流量为负。不能用采样格点平均速度代替边界积分。 |
| 物理验收与输出诊断分离 | 同文件 `:430`–`:577`、`utils/io/field_io.py:401` | 线性系统成功与守恒验收分别记录；笛卡尔网格采样散度只作导出诊断。保留失败即停止、几何哈希和有效性标志。 |
| 有针对性的回归测试 | `test_files/test_dolfinx_gmsh_backend.py` | 保留标签分离、流量逐口检查、无效配置拒绝、有限值检查、剪切流/刚体旋转应力检查的测试思想。 |

这些是可复用的设计，不代表把旧源文件原样接入即可。旧边界条件与本项目不同。
`DOLFINX_BACKEND.md` 与 `configs/physics_flow_config.yaml` 已交叉核对实现：
旧主模型是 X–Z 平面、所有 opening 给速度、压力只设数值基准。

# 需要重写但思想可复用的设计

| 检查内容 | 来源 | 3D 改写要求 |
|---|---|---|
| 变分式和代数布局 | `dolfinx_gmsh_solver.py:925`–`:1094` | 旧代码把两种普通单元放进 mixed element；新模型有全局 Real 约束，必须用 block/MixedFunctionSpace 或显式分块。普通 P2/P1 可以混合；不能把 Real 一并塞入普通 `basix.ufl.mixed_element`。全局 multiplier 只有一个自由度，不能每个 rank 一个互不相连的约束。 |
| PETSc / SciPy 路由 | 同文件 `:1247`、`:1305`；`environment-dolfinx.yml` | 旧 PETSc 为 `preonly+LU+MUMPS`，SciPy 为单 rank CSR `spsolve`，后者显式拒绝 MPI。新 3D 优先 CPU PETSc/MPI；SciPy最多保留小型代数对照用途，不能作为大型并行替代。旧依赖版本 `<0.11` 不沿用。 |
| 残差与物理阈值 | 同文件 `:1297`、`:1392`、`:1224` | 旧 PETSc 分支保存的相对残差是 NaN，只看 KSP 原因；新实现应独立计算 `||Ax-b||/||b||` 并检查有限。旧 `max(Q, machine_epsilon)` 用在 μm 单位下；若直接搬到 SI，约 1e-15 m³/s 的流量接近绝对机器 epsilon，阈值会失真。用有量纲参考尺度和相对/绝对容差，不能对 SI 流量随意加无量纲 epsilon。 |
| Velocity-gradient export | 同文件 `:1145`–`:1162` | 旧做法实际为 `Expression(grad(u))` 到 DG(Pk−1) 的插值，虽然日志称 project。仿射 P2 单元中梯度为 DG1；3D 要输出 3×3 梯度（s⁻¹），处理并行所有权。曲边映射下不能默认梯度仍是精确一次多项式。 |
| WSS 与采样 | 同文件 `:1647`–`:1703` | 保留真实壁法向和 FEM 梯度的思想。改写为三维切向牵引和 3D cell 定位，记录失配/未定位点；不能静默填零。 |
| FEM-to-particle field export | 同文件 `:1588`；`utils/flow/hybrid_velocity.py` | 旧导出通过显式单元内求值重建多项式，避免依赖 DoF 顺序，此思想可用。旧数组仅支持仿射三角形、两分量 X–Z 坐标；以后须重写为四面体、三分量、SI 单位与分布式所有权。Stage 0/第一版不实现粒子。 |

新出口流量是求解结果，不再要求每个出口等于旧的固定份额；应检验总流量守恒、
入口积分约束、壁面无滑移、有限值及压力/牵引约定。出口设置零外部表压的自然牵引，
入口总 Q 用全局乘子约束，入口法向牵引乘子作为结果输出。

版本依据：[DOLFINx 0.11 官方 Real element 说明](https://docs.fenicsproject.org/dolfinx/v0.11.0.post0/python/release_notes.html#the-real-element)
推荐独立 Real 空间配合 `ufl.MixedFunctionSpace` 与 `ufl.extract_blocks`。
这只是下一阶段的设计约束，本阶段没有组装 Stokes/Real 系统。

# 3D 中不应该继续使用的设计

| 旧设计 | 具体证据 | 为什么不能成为新实验边界条件 |
|---|---|---|
| inlet velocity Dirichlet | `dolfinx_gmsh_solver.py:1019`–`:1056` | 旧代码构造抛物线并插值固定边界速度；新入口只允许给总 Q，速度剖面必须由方程和全局约束共同求出。 |
| outlet velocity Dirichlet | 同一循环同时处理所有 terminal outlets；metadata `outlet_boundary_condition` | 旧方法预先固定出口流量/速度分布。新出口是大气表压的自然牵引 `σn=-p_ext n`，其中 `p_ext=0 Pa`；不能预设分流。 |
| pressure pin | 同文件 `:1058`–`:1073` | 旧全部速度边界留下压力常数自由度，钉一个压力点只为消除该自由度。新压力基准由出口牵引给出，不应再照搬 pressure pin，更不能把 pin 当成出口压力边界。 |
| effective 2D thickness | 同文件 `:982`；`DOLFINX_BACKEND.md`；`configs/physics_flow_config.yaml:28` | 旧用 `Q3D/h_eq` 得到平面通量，`h_eq=ΣπR²/Σ2R`。真实 3D 应直接在面积上积分得到 m³/s，根本没有等效厚度参数。 |
| μm 内部计算和运动学压力 | 同文件 `:1084`、`:530`；配置 `kinematic_viscosity_um2_s` | 新内部统一 m、s、m/s、m³/s、Pa、Pa·s、kg/m³；输入适配层才允许 μm→m。不得先求运动学压力再把 mmHg 导出惯例当内部物理单位。 |
| 旧额外物理和固定分流模型 | `configs/physics_flow_config.yaml`；求解 metadata `fixed_total_inlet_equal_terminal_shares` | 不引入 RBC、微泡输运、黏附、脉动、非牛顿黏度；不照搬旧等份 terminal flow。 |

参考测试多数是依赖较轻的单元测试，并非完整 3D 求解验证；旧测试还依赖另一个
vascular generator 包，因此本阶段没有在只读参考目录执行它们，也不宣称它们通过。
