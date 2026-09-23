# Stage 2 数学与符号约定

本阶段只开发并验证通用 Stokes 核心。几何为独立合成算例；Stage 1 的真实血管网格不参与求解。

## 物理量与边界

所有坐标和网格尺寸为 m。速度 u 为 m/s，压力 p 和全局乘子 λ 为 Pa，动力黏度 μ 为 Pa·s。密度 ρ 只用于 Reynolds 数，不进入 Stokes 矩阵或右端。

n 总是流体域的外法向。正式 API 的 Q_target > 0 表示流入的体积流量大小（m³/s），所以入口有符号积分为负：

\[
\int_{\Gamma_{in}}u\cdot n\,dS=-Q_{target}.
\]

壁面 u=0；入口不固定任何速度自由度，仅施加上述全局约束；入口自然牵引为 σn=−λn。出口为 σn=0，采用大气表压对应的零牵引。λ 是未知的入口法向推动强度，正式语义为 **inlet normal traction multiplier**，一般不等于任一局部静压。

\[
\sigma=-pI+2\mu\epsilon(u),\quad
\epsilon(u)=\tfrac12(\nabla u+\nabla u^T),\quad
-\operatorname{div}\sigma=0,\quad\operatorname{div}u=0.
\]

没有压力 Dirichlet 边界、压力 DOF 固定或常数压力 nullspace。出口零牵引已经提供压力基准：把 p 增加常数会改变出口牵引，不能仍满足原边界条件。

## 弱形式和代数符号

对壁面为零的速度测试 v、压力测试 q 和全局标量测试 η：

\[
2\mu(\epsilon(u),\epsilon(v))_\Omega
-(p,\operatorname{div}v)_\Omega
+(q,\operatorname{div}u)_\Omega
+\lambda\int_{\Gamma_{in}}v\cdot n\,dS
+\eta\int_{\Gamma_{in}}u\cdot n\,dS
=-\eta Q_{target}.
\]

若 D 是离散散度，则这条弱形式的原始矩阵为
`[A, -D^T, C; D, 0, 0; C^T, 0, 0]`。
如为对称矩阵装配把整条连续性方程乘 −1，并定义 B=−D，得到
`[A, B^T, C; B, 0, 0; C^T, 0, 0]`；物理压力的符号不变。
不能同时声称原始的 +q div(u) 行与 −p div(v) 列本身互为无符号转置。

速度采用 [P2]^3，压力 P1，约束采用 Basix 原生 Real 元素的独立 FunctionSpace。
使用三个独立空间组成 block/MixedFunctionSpace；不将 Real 塞入普通 mixed_element。
先在当前安装的 DOLFINx 0.11.0 / Basix 0.11.0 上独立验证 1、2 rank 的全局 DOF 数和块装配，再实现完整求解。

## 圆管的解析公式及边界兼容性检查

合成基准 R=5e−6 m、L=50e−6 m、μ=3e−3 Pa·s、ρ=1050 kg/m³、Q₀=1e−14 m³/s。
这些数字是 SYNTHETIC NUMERICAL BENCHMARK，NOT EXPERIMENTAL CONDITION。

无限长/充分发展圆管的 Poiseuille 场为
\[
u=w(r)e_s,\quad w=2\frac{Q}{\pi R^2}(1-r^2/R^2),\quad
p(s)=\Delta p(1-s/L),\quad\Delta p=\frac{8\mu LQ}{\pi R^4}.
\]

**先做解析代入，不预设有限管的端面边界与该解相容。**
在法向 n=±e_s 的端面，
\[
\sigma n=-pn\;\pm\;\mu\,w'(r)e_r.
\]
除中心线外 w′(r)≠0。Poiseuille 的端面需要切向黏性牵引；入口纯法向牵引、出口完整零牵引都没有该切向项。
因此本阶段指定的有限长圆管边值问题会有端部效应，λ 不必等于充分发展公式 Δp，压力也不必是从出口零值开始的整管直线。
不能通过网格加密消除模型边界条件的差别，不能把动力黏度应力偷偷改为向量 Laplacian 的另一种自然边界来追求该解析答案。

用户已明确选择：**保留正式边界条件；另增牵引一致的 Poiseuille 解析验证案例，并单独报告原圆管的端部效应。**

据此分为两组，使用相同几何、离散空间和全局流量约束：

1. `natural_traction`：正式边界不变，入口 σn=−λn、出口 σn=0。报告端部效应及对本边值问题的网格收敛，不把它误标为整管精确 Poiseuille。
2. `poiseuille_traction_reference`：仅在解析回归接口补充已知的端面切向牵引 t_tan=(n·e_s) μ∇⊥w。入口 σn=−λn+t_tan，出口 σn=t_tan。它不规定任何入口速度值，λ 仍由流量约束求得；此附加案例的精确解才是上述 Poiseuille 场和 λ=Δp。

补充牵引对速度测试的作用出现在右端 `∫ t_tan·v ds`。生产接口不开放此参考模式。
解析速度/压力/λ 及其网格收敛的 hard gates 应用到第二组；正式组仍必须严格满足 Q、质量闭合、唯一压力基准、符号、线性、黏度/密度和 MPI 回归。任何数学 gate 失败均如实报告 FAIL。

## 输出与限制

λ 原值保存，不取绝对值。入口与出口流量均由 FEM 解重新积分并做 MPI reduction。
入口相对约束误差和质量闭合均要求 ≤1e−10。数学反向 RHS 只在内部回归接口使用，正式 API 拒绝非正 Q。
使用 CPU PETSc direct LU/MUMPS，无 pressure pin、无 nullspace 注入，不使用 GPU，不计算 WSS。
Stage 2 完成后停止，不进入 Stage 3。

API 依据：[Basix 0.11 real_element](https://docs.fenicsproject.org/basix/v0.11.0/python/_autosummary/basix.ufl.html#basix.ufl.real_element)、[DOLFINx 0.11 官方源码](https://github.com/FEniCS/dolfinx/tree/v0.11.0)。边界兼容性结论由上面的应力直接代入获得。
