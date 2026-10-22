# 实际 production stabilization 审计

状态：**SOURCE_CONFIRMED**。审计对象为服务器生产二进制，不以 upstream 自动代替。完整链见 [源码身份](PRODUCTION_SOURCE_IDENTITY.json) 与 [2431 文件清单](EXACT_DEPLOYED_SOURCE_MANIFEST.json)。二进制 SHA256 `43cd7d7b9ae9bde91275982225c6bf084859141cdaf2ec1030e1b83c3b94691d`，基础源码 commit `c3f0bb892b765b718f61069ecd9726dbc6d177fd`。基础树仅 PETSc/preconditioner reuse、stop acknowledgment 等执行侧改动；正式 host-sync 补丁只改 PETSc device→host 同步。已核对构建 receipt、目标文件依赖、source manifest 和当前文件/二进制 SHA；未声称做过可重复二进制重建。

## 实际调用链和参数

`Integrator::assemble_equations → eq_assem::global_eq_assem (Equation_fluid) → fluid::construct_fluid → fs::get_thood_fs → nn::gnn → fluid_3d_m / fluid_3d_c → PetscLinearAlgebra::assemble`。当前 TET4、`nFs=1`，故 `vmsStab=true`；momentum 与 continuity 分别调用 iOpt=1/2，但此配置两个空间均为 P1。PETSc 包装器复用 FSILS 的 CPU 元素装配，再用既有 PETSc CUDA 线性求解；“GPU solver”不等于全部装配在 GPU 上。

实际 XML：CFD dt=0.01 s（不是 particle dt=0.0005 s），rho=1056 kg/m³，mu=0.00345312 Pa·s，nu=3.27e−6 m²/s；constant viscosity；spectral radius=0.5；无 ALE、Brinkman、URIS、体力。Taylor–Hood 参数省略，默认 false。入口 Flat/Impose_flux/Zero_out_perimeter，外向 Q=−3e−7 m³/s；历史出口 traction 556.97065554884477/0/17.137427153562019 Pa，WALL=0，backflow coefficient=0。没有改用新边界设计，也没有重标历史校准场。

## 精确公式与单位

输入四面体若 `((x1−x0)×(x2−x1))·(x3−x2)>0`，reader 交换节点0/1；节点3不变。`nn` 的 TET4 solver→basis 排列为 `{1,2,3,0}`，所以参考形函数导数为 `[I | −1]`，**参考原点是节点3**：

```text
J = [x0−x3, x1−x3, x2−x3]
K = J^(-T) J^(-1)                              [m^-2]
kT = 4(ctM/dt)^2 + (mu*Brinkman/rho)^2 + URIS^2 [s^-2]
kU = u^T K u                                   [s^-2]
kS = ctC (mu/rho)^2 (K:K)                       [s^-2]
ctM=1, ctC=36; 当前 kT=40000
τM = 1 / (rho sqrt(kT+kU+kS))                   [m^3 s/kg]
τC = 1 / (τM tr K)                             [Pa s]
```

未替换成 h_min，也没有添加未在实际源码中的 tetra metric 对称化。原0/1交换只交换参考轴，K不变；任意换原点通常会改变K。known cells 的四积分点使用原 TET4 对称规则 s=(5+3√5)/20、t=(5−√5)/20，各权1/24，另保留单元中心。

`u' = −τM [rho(a − f + (u·∇)u) + ∇p − ∇·(mu(G+Gᵀ)) + Brinkman/URIS项]`。当前后两类阻力项关闭。P1、常黏度的单元内部强形式黏性二阶导数为零；这**不表示全局弱形式黏性力为零**，更不能忽略跨面应力。

VMS 开启时 `ua=u+u'`，`pa=p−τC div u`；另有 `τB=rho/√(u'ᵀKu')`，源码 is_zero 时对分母量置 machine epsilon。`pa` 在 momentum 弱式中产生 `+τC div(u) div(test)`，是有限的 grad-div/LSIC-like 项。continuity 残差为 `Σg w[Nq div u − u'·∇Nq]`，pressure-pressure tangent 含 `w αf γ dt τM ∇Na·∇Nb`。它要求连续压力测试空间下的加权残差为零，**不等价于每个四面体 div u=0**。因此 nodal residual 收敛和非零 cell Qnet 可以同时出现。

Generalized-alpha：αm=5/6、αf=2/3、γ=2/3；`uα=(1−αf)u_old+αf u_new`，`aα=(1−αm)a_old+αm a_new`；tangent 使用 `T1=αf γ dt`、`amd=αm/T1`。稳态加速度趋零，τ 中 `4/dt²` **仍保留**。本轮离线量级使用已有 frozen steady 速度，不冒称读到了原 Newton 中间态积分点载荷。

## 与论文的对应及边界

[Jia–Esmaily 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10625993/) Eq5/6 对应 `τSUPG=rho τM`、`νC=τC/rho`；其 CI=3，当前代码 ctC=36，完整 VMS 弱式也有额外项。故 `JIA_ESMAILY_TRADITIONAL_FORM_MATCH=true` 仅指时间尺度公式家族。

论文 Eq8/9 用 `ω²` 替代 `(2/dt)²`，`ω=||∂t u||L2(Ω)/||u||L2(Ω)`：全域欧拉时间导数，不能换成局部 `Gu` 或材料加速度。首步 ω=2/dt；随后用上一接受步状态，每时间步更新一次，在非线性迭代中固定。ALE 形式按 Eq10 使用网格坐标时间导数和相对速度。本文只求 ω→0 的稳态极限，保留当前 ctC=36、metric 和速度，标为 **OFFLINE_HYPOTHETICAL_TIME_CONSISTENT_EVALUATION**；未测全域 ω，永久删 kT 也不是论文完整算法。

论文 pipe 与复杂血流 benchmark 说明该机制在部分条件下重要，不能据此推定 BraVa 局部 u/G 改善。论文给出的 [MUPFES 作者页](https://sites.google.com/site/memt63/tools/MUPFES) 本次直接读取返回 Google 登录页，未获取并核对其公开实现；不声称作者未公开代码。本项目 Stop C 已独立成立，无需猜实现继续试算。

## 可定位源码证据

所有行号对应实际基础源码树 `Code/Source/solver/`，文件 SHA 逐项核对生产源清单。

| 项目 | 文件/行号 | SHA256 |
|---|---|---|
| assembly_dispatch | `eq_assem.cpp:403` | `b07ac95bede7dd7ca10e35778c3f72a3bfdd691166e1a72c000d26fc2270b30d` |
| fluid_entry | `fluid.cpp:480` | `df04669e8393de6f1eef8fcc18ce883dd39369bad64c3888f665b5525ef2ee64` |
| vms_switch | `fluid.cpp:496` | `df04669e8393de6f1eef8fcc18ce883dd39369bad64c3888f665b5525ef2ee64` |
| momentum | `fluid.cpp:1768` | `df04669e8393de6f1eef8fcc18ce883dd39369bad64c3888f665b5525ef2ee64` |
| continuity | `fluid.cpp:1443` | `df04669e8393de6f1eef8fcc18ce883dd39369bad64c3888f665b5525ef2ee64` |
| tau_momentum | `fluid.cpp:991,1279,1655,2001` | `df04669e8393de6f1eef8fcc18ce883dd39369bad64c3888f665b5525ef2ee64` |
| tau_continuity | `fluid.cpp:1000,1290,1673,2018` | `df04669e8393de6f1eef8fcc18ce883dd39369bad64c3888f665b5525ef2ee64` |
| time_parameters | `initialize.cpp:464,468` | `3b8084bb7ccbe4707127878dbef81d838cf28deb80de09b422da70f25c90f55f` |
| metric | `nn.cpp:814` | `ba138fc9aaba3424ffbcbb1f2ef72cfd2fff11298466100958d1bbb96dd543b8` |
| metric_order | `nn.cpp:174` | `ba138fc9aaba3424ffbcbb1f2ef72cfd2fff11298466100958d1bbb96dd543b8` |
| mesh_orientation | `read_msh.cpp:657` | `2e9c4d3e92c673d1ca5cb82482dbce255b9b6dd6df8777f1e675a2301aa8be8a` |
| TH_reader | `read_files.cpp:1741` | `57ab378dc3f7e8442b28b99a908b91a546db2f63a769244a1a9e2510d7e624b4` |
| TH_pressure | `fs.cpp:365` | `470029137da2ff1c245f0a12135752303594b7a0c37bd3a947fcb675dfc65e45` |
| TH_edge_pressure | `fs.cpp:394` | `470029137da2ff1c245f0a12135752303594b7a0c37bd3a947fcb675dfc65e45` |
| TH_integrator | `Integrator.cpp:140` | `45e31f4ee1706de5eacf4f5b4bcff8e4395158ca34cb3f21c4f38c0941f15a7f` |
| TH_pressure_corrector | `Integrator.cpp:1087` | `45e31f4ee1706de5eacf4f5b4bcff8e4395158ca34cb3f21c4f38c0941f15a7f` |
| PETSc_assembly | `PetscLinearAlgebra.cpp:58` | `fef89d109493d420ee731fa40ecd74845fcde026debd7c572ef3368d3d83e8bc` |
| VTK_TET10_write | `VtkData.cpp:737` | `888b1826e32d1a63318e0ed6f3f829e6c7586ab2b457766e79d5c9808eab711a` |
| VTK_TET10_read | `vtk_xml_parser.cpp:51` | `e72a7608ac340feca357519d543838226bee532fa113acbb5742e2577a1f966c` |
| restart | `output.cpp:298` | `da2546ac340479b5e1fdb1dbd51dab497baadffc367650272fb93ea8d639f1e3` |
| arbitrary_basis | `FE/Basis/BasisFactory.cpp:51` | `8cdaef5238bb22d770358bdcccc5bf24f6e235d57238f3b661d5fa1f2d3972bf` |
| TET4_quadrature | `nn_elem_gip.h:214,474` | `870f77c32ddb2f4a71eacd59e853916edf607f763d7cc3a27dba34c14b3e0082` |
| TET4_quadrature_default | `Parameters.cpp:3126` | `bc37aad1b96c009ed452955c725c30bb68d6b37e3a0db3fc3760a2d7e4c51f62` |
