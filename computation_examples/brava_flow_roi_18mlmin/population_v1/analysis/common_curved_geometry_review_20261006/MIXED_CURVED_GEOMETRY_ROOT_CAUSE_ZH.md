# Mixed P2/P1 共同曲面几何：根因与补丁范围

当前结论：**体装配修复通过，边界弱 Dirichlet 分支仍阻断**。补丁只存在于本轮独立 validation 构建，未修改生产源码、binary、P1 结果或 stabilization 公式。

## 正确的数学语义

唯一物理映射为 `x(ξ)=Σ_A N_A^P2(ξ) X_A`，`J=∂x/∂ξ`。速度和压力使用不同参考基函数，但都位于这个映射定义的物理单元中：

- `u=Σ_A N_A^P2 U_A`，`∇x u=(∇ξ u) J⁻¹`（此处速度分量按行排列）。
- `p_h=Σ_i λ_i P_i`，`∇x p_h=J⁻ᵀ Σ_i P_i∇ξ λ_i`。
- `Kxi=J⁻ᵀJ⁻¹`，`w=w_ref detJ`。压力基函数不得另行产生几何坐标或积分权重。
- 二阶导数满足 `H_x N=J⁻ᵀ[H_ξ N−Σ_i (∂N/∂x_i)H_ξ x_i]J⁻¹`。漏掉曲率项或沿用其它积分点的几何都不是这个映射的 Hessian。

**参考 P1 压力不是物理坐标中的完整线性多项式空间。** 任意物理线性压力 `p=c0+c·x` 拉回曲面 P2 几何通常是二次函数，四角点 P1 无法精确表示。对任意物理二次速度也有类似限制：其拉回通常达四次，不能把 affine P2 的精确表示结论泛化到任意 curved P2。

本轮在查看修复结果前冻结测试：两原单元的投影只改变 x/y，z 映射保持仿射，因此选择非恒定压力 `p=0.3+(z−z0)/H` 和速度 `u_z=(z−z0)²` 作可精确表示的物理场测试；同时保留一般三分量物理线性压力的角点插值作为表示误差对照。后一组只要求正确计算 `p_h/∇p_h`，不冒称精确恢复任意物理线性压力。详细尺度见 [冻结协议](COMMON_GEOMETRY_HARNESS_PROTOCOL.json)。

## 六个根因问题

1. **为什么覆盖 Jacobian？** `nn.cpp:814–899` 的 `gnn` 同时接收某个场的参考导数和该场节点坐标，构造几何、变换场导数，并通过引用覆盖 `Jac/ks`。`construct_fluid` 分别传入10节点 `xwl` 和4角点 `xql`，却共用一组 `Jac/ksix`。
2. **API 耦合还是顺序错误？** 两者都有。仅调换调用顺序只能改变最后留下的积分 Jacobian，不能让压力梯度使用高阶几何。`fs::get_thood_fs` 已在相同参考积分坐标上生成两组基函数；问题发生在物理变换阶段，不是 quadrature 坐标生成阶段。
3. **为什么连续性首点与后续点不同？** `fluid.cpp:705–723` 先调用速度 `gnn`，再调用压力 `gnn`；TET4 压力的 `lShpF=true` 使后者只在 g=0 执行。首点 `Jac/ksix` 留下角点仿射值，后续点留下速度曲面值。这个缓存 flag 只描述参考基函数导数恒定，不能证明曲面映射后的物理导数恒定。
4. **动量有同类问题吗？** `634–651` 中压力先、速度后，故积分权重和 metric 最后是速度几何；但 `Nqx` 一直来自压力角点映射。没有首点权重异常不代表动量 mixed 语义正确。
5. **pressure-gradient / PSPG / VMS / tau 使用什么几何？** `fluid_3d_m/c` 从 `Nqx` 得到压力梯度；momentum 的 tauM 来自传入 Kxi，原动量循环的 Kxi 是速度几何，原连续性首点则被压力几何覆盖。原生 `nFs=2` 会关闭连续性 PSPG、动量 tauC/tauB，但动量仍有 tauM/up 项；本轮保持这组既有开关及全部公式/系数。新补丁统一传入高阶 metric、压力变换和权重，没有以“修几何”为名更换 stabilization。
6. **boundary traction/Cauchy 有同类问题吗？** 分路径回答。`set_bc_trac_l:2008–2027` 通过 TRI6 几何求法向模和权重，插值外部给定的向量牵引，不调用 TET4 pressure geometry；真实 native 夹具通过。强 Dirichlet 遍历 face.gN，包括三个边中点。**`set_bc_dir_wl:1445–1457` 则确有同类缺陷**：高阶速度反解后，重新以 `xql` 四角点坐标独立反解压力。曲面 wall 的点可落在旧平面单元外，触发原生异常；cap 上即使未越界，压力参考坐标仍不一致。它与 prescribed full Cauchy traction 是两个不同路径，不能相互代证。

## 最小验证补丁

[COMMON_GEOMETRY_PATCH.diff](COMMON_GEOMETRY_PATCH.diff) 只修改独立副本的 `fluid.cpp`，新增局部 helper header。`SVMP_VALIDATION_COMMON_GEOMETRY` 为显式编译门，仅在 3D、nFs=2、TET10 分支生效；未定义宏或进入 P1/P1 时仍执行旧装配块。

每个积分点用 TET10 参考基函数计算一次共同几何，中心化节点坐标以减少相消；用原库 `mat_det/mat_inv` 得到 detJ/invJ。P2 和 P1 参考导数分别乘同一 inverse；Jac、metric 和积分权重均来自该几何。基函数/积分规则、材料、generalized-alpha 和 fluid_3d_m/c 的物理装配代码均不变。

初版复用旧 `gn_nxx`，新增严格 Hessian 检查出现20项失败（同一几何检查在5个场上重复）：H²归一化最大误差3.20e−11至4.91e−11，高于固定1e−11。修正版本采用同一逆矩阵、中心化几何曲率及完整 Hessian chain rule；未放宽容差。最终该误差最大3.97e−12。初版 patch、身份和失败记录均保留在 `ATTEMPT_01_*`，不混入最终通过计数。

实际生产对象链接的原始装配、仅观测重编译、补丁无观测/有观测四个独立 native binary 相互对照，避免以复制装配算法冒充 native。common reference 的几何、导数和 Hessian 由独立 NumPy 计算，再传给**未修改**的 native physics kernels；连续性和零速度压力载荷另有直接弱形式积分对照。reference 不声称独立验证整个 Navier–Stokes kernel。

压力空间物理表示误差、原生边界反解的停止容差和未执行的 CFD 都没有被这个体装配补丁消除。按本轮 Phase F 停止规则，未继续 patch 弱边界，也未自动运行 pipe CFD。
