# STEP 3 — 真实血管纯流体物理边界与短程验证

## 1. Executive summary

**STEP3_STATUS=AUTO_PASS_HUMAN_PENDING; STEP3_AUTO_CHECK=PASS.**

当前结果为第 2 次 MPI1 验证，完成 1000 步（2.03668106e-06 s）。第一次也完成 1000 步，因流量监测器缓存 MPI 代理的缺陷而归档，不能使用其流量列。两次共 2000 个流体时间步，没有物理参数调整。

最终入口 Qin=2.39300257e-15 m³/s，相对目标误差 12.566%，满足运行前设定的短程 25% 容差。该容差不是最终准确性标准。两个正表压出口仍回流；closure=1.50637，尚不能证明稳态质量闭合。失败项：无。Step3 人工检查仍待完成。

## 2. Frozen Step 1 / Step 2 inputs

旧工程 HEAD：fd21a850a16d0ba17ef3d521123badd53864a3a2；HemoCell HEAD：5a410848bd5c57d5ae1c171112e78eab4a82e650。冻结 STL SHA256：840da5e1c43ec31ac70ec781cbf75b32940bc73538b5dba83c58af1ddf5e36fb。

Step2 人工检查已通过。物理坐标、dx、493×497×280 格点、四端口中心/法向、191 cap 三角形与全部 497 原标签复用。原生重建只为提供 Guo 所需 TriangleHash，逐项检查原生 bulk flags 与 closed lumen，没有重新识别端口。物理 closed lumen=182694 格点；opened 标记包络=183191。

结束工作文件完整性检查：PASS。Step1/Step2 全文件 SHA256、HemoCell 已有 tracked files SHA256 与 Step2 代码均核验；现有 GoogleTest compatibility patch 保留。只有 Step3 新目录允许新增源码。

审计例外必须保留：首次结束审计的 Git status/diff 调用刷新了旧工程 .git/index 及 .git 目录的元数据，尽管设置 GIT_OPTIONAL_LOCKS=0。旧工程工作文件的前后属性与完整 git status 列表没有变化；此前没有索引二进制 SHA256，不能宣称索引逐字节未变。这属于只读流程的 Git 元数据例外，不隐藏时间戳或回滚索引。后续审计使用 RUN 内复制的 GIT_INDEX_FILE，并禁用 diff 自动刷新。详见 provenance/frozen_integrity_after.json。

## 3. Physical BC provenance

contracts/physical_bc_contract.json 为各物理量保存值、单位、来源、上下文、状态与 SHA256。链路为 s3/1D 配置、s4 延长段修正 → 当前 cfd_flow.yaml 与 validated_contract.py → 已接受 runtime contract → 下游 promotion 输出。

早期 s3/s4 的 Q=7.693508475538942e-16 m³/s 与当前目标不同，采用当前配置、常量、实际输出三者一致的目标。出口则保留当前合同仍采用的延长段压力修正结果。旧工程只读，未执行 s3/s4、Seeder、Musubi 或旧 CFD。旧 adaptive_flux_pressure、pressure_eq、wall_libb、数值压力偏置、tau/dt 分列为 OLD_SOLVER_NUMERICAL_METHOD，不自动迁移。

## 4. Physical BC contract

密度 1056.0 kg/m³，运动黏度 3.27e-06 m²/s，动力黏度 0.00345312 Pa·s。

目标 Q=2.7369132390905703e-15 m³/s；目标质量流量=2.8901803804796421e-12 kg/s。表压参考 GLOBAL_STRUCTURAL_LEAVES_ZERO_GAUGE；绝对生理血压 UNVERIFIED。所有生理角色为 ASSUMED，结构端口身份已确认。壁静止、刚性、no-slip。

| Port | Role | Gauge pressure (Pa) | rhoLU target | Final outward Q (m³/s) |
|---|---|---:|---:|---:|
| outlet_01 | ASSUMED | 14.5449781012743 | 1.0000042893791 | -9.59800583e-17 |
| outlet_02 | ASSUMED | 132.204549223176 | 1.00003898771286 | -1.29538256e-15 |
| outlet_03 | ASSUMED | -13.7006266733115 | 0.999995959623917 | 1.79622994e-16 |


## 5. HemoCell / Palabos BC API survey

实际读取 pipeflow、pipeflow_with_preinlet、curvedflow_with_preinlet、stl_preinlet、AR2、preinlet_shear 的相关实现及 on/off-lattice public API。周期体力管流不匹配当前物理合同。0P/0N/1P/1N/2P/2N 依赖 Cartesian population 方位，不能用大 Box3D 冒充斜 cap。

当前 bundled Palabos aneurysm 示例已使用 Guo 曲面速度入口与恒定密度压力出口。逐端口结论、源码上下文和 SHA256 在 contracts/bc_api_feasibility.md/json。

## 6. Arbitrary-port BC feasibility

四个 cap 均由 frozen triangle IDs 标记。public tagDomain 只接受其三个顶点都属于 cap 的三角形，随后验证整张网格的标签恰好对应 191 个冻结 cap 三角形，其他面均为 wall tag 0。无 bbox 端口、新圆形 cap、重定位或新 population reconstruction。

原生 Guo 在 outerBorder 的 dry/ghost 节点完成边界；这些是数值支持，不是新增物理管腔。Step2 的 opened 标签完整保存，物理体积与 native 支持分开统计。width-2 envelope 只用于第二邻点存储。真实 completion 中四个端口与 wall profile 均被调用。

DensityNeumann 固定密度并外推法向动量，代表恒定压力及相应速度外推，并非把指定压力替换成任意零梯度压力。曲面 cap/wall 接缝有有限网格误差，仍需后续独立精度研究。

## 7. Lattice-unit conversion

当前 mechanics/constantConversion.cpp:36–54 给出 νLU=(τ−0.5)/3、dt=νLU·dx²/νphys、dm=ρphys·dx³、df=dm·dx/dt²。D3Q19 cs²=1/3 与 BGK 二阶 equilibrium 来自当前 Palabos；其二阶矩为 ρcs²I+ρuu。

固定参考表压 p=(ρLU−1)cs²ρphys(dx/dt)²；参考物理表压 0 Pa 对应 rhoLU=1。不得使用会减去瞬时平均密度的便捷 OffLatticeBoundaryCondition.computePressure 来输出本合同的固定参考表压。u物理=uLU·dx/dt；Q物理=QLU·dx³/dt；质量单位=ρphys·dx³。

全部数值、单位与源码证据在 contracts/lattice_unit_contract.json。时间步独立推导，未复制旧 Musubi 的 dt 或数值压力偏置。

## 8. dx/dt/tau preflight

requested dx=1.9999999999999999e-07 m；effective dx=1.9989918081065344e-07 m，只称 CANDIDATE_BASELINE。STEP3_CANDIDATE_TAU=1.0，νLU=0.16666666666666666，dt=2.0366810646671923e-09 s。

入口预期均速 0.00035 m/s，Re=0.000337732，入口 Mach=6.17647e-06。保守启动筛查假设完整压力范围跨一个 dx 作用 1000dt 且无黏性衰减，上界 Mach=0.0248485<0.05；这是筛查估计，不是数学稳定性证明。

压力目标密度区间 [0.999995959624,1.00003898771]。运行前冻结的流体密度限值 [0.99,1.01] 和 Ma≤0.05 为本候选设定，没有机械复制旧 CFD 阈值。

## 9. Wall BC implementation

HEMOCELL_WALL_METHOD=原生 GuoOffLatticeModel3D + NoSlipProfile3D；regularized、second-order、all directions。BGK 流体，outside 设 NoDynamics，outerBorder 保留原生边界支持，与官方 aneurysm 设置一致。

test-local 包装仅返回原生 profile 数据。没有自行编写 population reconstruction，没有复制 wall_libb，没有改核心或官方案例。近壁格点中心非零速度本身不能当作壁穿透证据。

## 10. Inlet flow implementation

物理合同是 TARGET VOLUMETRIC FLOW。数值实现为 u=−n·Q/Aproj，其中 Aproj 是原生已 inflate 的真实 tagged cap 有向投影面积。所有 cap 三角面上 ∑u·n_triangle A_triangle=−Q，相对积分误差=0。

原生有效面积 7.8213689069616978e-12 m²；冻结 cap 物理面积 7.8197530070341531e-12 m² 保持不变。二者差异来自 Step2 已使用的官方 0.001 LU inflate，归一化的是数值积分，不是修改物理 cap 面积。

uniform plug 明确标为 STEP3_NUMERICAL_INLET_PROFILE_ASSUMPTION。当前物理合同没有实验 profile；这是最小、可积分核查的原生速度实现，人工 extension 内测量其实际响应。旧历史 PARABOLIC 不当作生理事实。

实际 lattice 流量另行积分，最终误差 12.566%。预先冻结的 25% 容差只用于短程一致性，不是最终 CFD 准确性结论。无 controller、无拟合。

## 11. Three outlet pressure implementations

| Port | Role | Gauge pressure (Pa) | rhoLU target | Final outward Q (m³/s) |
|---|---|---:|---:|---:|
| outlet_01 | ASSUMED | 14.5449781012743 | 1.0000042893791 | -9.59800583e-17 |
| outlet_02 | ASSUMED | 132.204549223176 | 1.00003898771286 | -1.29538256e-15 |
| outlet_03 | ASSUMED | -13.7006266733115 | 0.999995959623917 | 1.79622994e-16 |

每个 cap 的 rhoLU=1+p_gauge/[cs²ρphys(dx/dt)²]，分别传给原生 DensityNeumannBoundaryProfile3D。三者不平均、不全部归零，负表压保留。NUMERICAL_STARTUP_RAMP：Q 与 gauge pressure 从 0 至第 10 步线性达到完整目标，此后固定。ramp 是数值启动安排，不是物理参数。

真实 reconstruction 由原生 Guo 根据边界固定密度与流体邻居动量完成，test-local 层只组合 native profile 数据。

## 12. Flux measurement method

详见 contracts/flux_measurement_contract.md/json 和 flux_quadrature.tsv。各面位于 cap 内侧 4dx，处于已知人工 extension 首四分之一内。与全冻结 STL 相交后选取包含投影 cap 中心的唯一闭合局部轮廓，并检查没有到达邻支。没有新识别端口或修改几何。

真实轮廓三角剖分，积分三角形最大边≤dx，每三角形三个求积点。每点用 8 个实时格点三线性速度插值，必要时包含 native ghost 支持。保留法向投影与实际面积权重，非 mean velocity × contract area。Qoutward=∑(u·n)dA；Qin=−Qoutward(inlet)，出口保留符号。

首次缺陷：MPI MultiBlockLattice::get 返回复用 distributedCell，缓存其地址导致所有采样点别名。修正为 public sparse locate + getComponent 的 atomic 格点；1767 个坐标对应 1767 个独立地址，第 10 步逐点与即时 MPI get 读取的速度差为 0。原始错误源码、程序、日志与数据保留。

内移截面截去端部短段，瞬时 closure 不等价于全体积精确离散质量残差。物理质量另按 closed lumen 统计。积分加密与网格精度尚未验证。

## 13. Runtime staging

| Step | Safety | Time (s) | Qin (m³/s) | Target error | Closure |
|---:|---|---:|---:|---:|---:|
| 0 | PASS | 0 | -0 | 1 | 0 |
| 1 | PASS | 2.036681e-09 | -0 | 1 | 0 |
| 10 | PASS | 2.036681e-08 | 5.837275e-16 | 0.78672 | 2.87527 |
| 100 | PASS | 2.036681e-07 | 2.330666e-15 | 0.148433 | 2.35879 |
| 1000 | PASS | 2.036681e-06 | 2.393003e-15 | 0.125656 | 1.50637 |

当前第 2 次执行的连续链为 0→1→10→100→1000。每步安全检查正常才推进，重要阶段输出。两个执行各 1000 步，总计 2000 步；第一次 flow 列标记 INVALID，不并入当前 CSV。所有物理参数、BC、ramp、单位与 quadrature 相同，修正只针对监测器。

两次 0/10/100/1000 的物理原始场 SHA256 比对在 diagnostics/physical_repeat_check.json，可核查监测器修正未改变求解轨迹。MPI1、OMP 1，只用系统编译器/MPI；MPI2=UNVERIFIED。没有长稳态计算。

## 14. Stability diagnostics

NaN/Inf、MPI_ABORT、segfault 均无。最终 rhoLU=[0.99999589665461,1.000039547467]，u_max=0.000602064014 m/s，Ma_max=1.06246665e-05。

每个时间步检查物理流体与 ghost 的有限值、正密度与 Mach；物理流体另检查 [0.99,1.01] 密度范围。CSV 每 10 步及关键阶段采样；逐步硬 gate 全程完成与退出码为 0 是执行证据。

## 15. Flow diagnostics

最终 Qin=2.393002566e-15 m³/s，ΣQout=-1.211739624e-15 m³/s。closure=|Qin−ΣQout|/max(|Qin|,Qtarget·1e−12)=1.506368，入口目标误差=12.56564%。

物理流体总质量 1.541075268e-12 kg；相对初始质量漂移 5.997932e-06。outlet_01 与 outlet_02 有向内回流，outlet_03 向外，符号没有用绝对值掩盖。完整回流时间点见 diagnostics/automatic_acceptance.json。

均匀 0 表压静止初始化后，正表压出口可先向内推动流体，负表压出口可先向外抽流。当前约 2 μs，压力传播和储质量过程尚在继续。这个解释是结合初始化、压力差和实际流向的数值推断，并非实验生理验证。入口目标接近独立验收；短程不设稳态 closure 门槛。小相对总质量漂移不能单独证明出入口稳态守恒。

## 16. ParaView outputs

output/flow.pvd 包含 0/10/100/1000 四帧，VTU 是原 Step2 体素的稀疏子集，不重采样。字段包括 FluidMask、PortLabel、速度矢量/大小、rhoLU、物理密度、固定参考表压、GlobalLatticeId。

物理场先 Threshold FluidMask=1。原 497 个端口标签在边界支持格点上原值保存，PhysicalFieldValid=0；其速度、密度等零值仅占位，不是实际 ghost 测量。四个 measurement_*.vtp 标明积分面。VTK reader 对四帧逐字段与原标签、格点中心位置回读校验 PASS。单位在 output/field_units.json。

## 17. Limitations

本轮只验证短程数值安全、目标入口的定量接近、明确的三出口压力施加、瞬态可解释性与输出一致性。25% inlet tolerance 是运行前声明的 smoke 容差；12.6% 误差不能称高精度。

出口回流和 closure≈1.51 明确显示尚无稳态质量闭合证据。profile 假设、cap-wall 接缝、native ghost 和截面积分都存在待验证的离散误差。没有实验、稳态、网格、时间步或 MPI2 验证。候选 dx/tau 不能升级为最终值。

## 18. Unverified items

- Step3 人工 ParaView
- 生理入口/出口身份
- 实验流量/黏度/profile
- 绝对生理血压
- 稳态收敛
- 网格/时间步收敛
- 壁与 cap 接缝离散误差
- 流量积分收敛/精确离散质量残差
- MPI2
- 生产 CFD / RBC readiness

## 19. Step 3 status

**AUTO_PASS_HUMAN_PENDING**；自动检查=PASS；Step3 人工检查=PENDING。Step2 人工通过不传递到新的流场。

NEXT_RECOMMENDED_ACTION=manual ParaView review of Step 3 pure-fluid short run。

本轮停止在人工审查前，不添加 RBC/particles，不自动做长稳态、网格收敛，不 git add/commit/push。
