# 开源实现能力核查

结论：Pecnut 上游能执行本轮的自由空间、多固定硬球、一个活动球的完整力/力矩混合问题。它是固定硬球墙，不是 regularized blob 墙。工具能力已实测通过；墙体能否逼近连续 no-slip wall 仍需另外判断。

|候选|固定/规定运动|平动、转动、力矩|应力偶与润滑|任意坐标|本轮使用|
|---|---|---|---|---|---|
|Pecnut/stokesian-dynamics|YES，原生 UFTE 模式和完整 resistance 矩阵|YES|YES，11 DOF/球；原生两球润滑 excess|YES|主实现，实际测试|
|RigidMultiblobsWall|YES，prescribed kinematics 与刚体约束|YES，刚体 K 矩阵|原生 RPY 刚体；不能将 RPY 等同硬球精确润滑|YES|只读能力审查；未作第二个运行后端|
|RigidBodyIB 该独立仓库|不提供完整规定刚体运动求解流程|marker mobility 为基础，非现成完整 6×6 球响应|经验 marker mobility / RPY；非本轮完整混合求解器|marker pair 可输入坐标|只读能力审查|

- [Pecnut 原始仓库](https://github.com/Pecnut/stokesian-dynamics)，提交 `6b9117d1601c5e50dfea0475fd4f10d53447174a`，MIT。README 的固定球例子与 `functions/generate_grand_resistance_matrix.py`、`generate_Minfinity.py`、`generate_R2Bexact.py` 已读。原生方法为 `R = inv(M_infinity) + R_2B_excess`，含 force/torque/stresslet；可对全部墙 U/Ω/E 置零，取目标球六个响应列。没有使用平墙 Green 函数。
- [RigidMultiblobsWall](https://github.com/stochasticHydroTools/RigidMultiblobsWall)，提交 `8d41e464d7de9b6514a85a741dd227f1219e3a49`，GPL-3.0。`body/body.py` 的 K 矩阵、`multi_bodies/multi_bodies.py` 的 prescribed_kinematics 分支，以及 no_wall mobility 路径确认可构造固定刚体 blob 体系。但目标球的离散自由空间转动精度需要单独验证，本轮未据此声称已有等价后端。
- [RigidBodyIB](https://github.com/stochasticHydroTools/RigidBodyIB)，提交 `8f1abe43d0a1ac55ffca10a6f20b0293fcc1a747`，GPL-3.0。该小仓库提供 `getEmpiricalMobilityMatrix`、`getRPYMobilityMatrix`，不是完整 IBAMR 应用。独立矩阵代码本身不需实际 Eulerian 网格；经验模型有 DX、IB kernel 参数。不能从仓库名称推定已具备所有 mixed rigid-body 功能。

## 主实现与薄适配层

源码副本保留原始水动力函数。Numba 启用为运行配置；不同 scalar 数据 SHA 使用不同 JIT cache，避免编译常量污染。分块 Schur 消元和 SciPy Cholesky 只是标准矩阵运算：固定墙 11N×11N M∞ 因子分解可复用；活动球保留完整 11 自由度再取六个力/力矩响应。没有重写水动力核，没有使用伪逆、对角校准或新经验润滑公式。

墙—墙润滑项处于完整 resistance 的墙—墙列，乘以严格为零的墙运动。因此只计算目标响应时可跳过这部分；这不等于删去活动球—墙润滑。本轮以 7、13、7 个墙球的原生完整矩阵（含接触墙球案例）检查，缓存与原生目标矩阵最大相对差约 `3.98e-16`。接触墙球案例仅诊断；正式间距比为 1.01，无物理重叠。

上游自带 scalar 表只有 λ=.01,.1,1 及倒数；本轮 β=.5、.25 使用上游原始 Fortran/midfield 和 nearfield 生成器生成，未在尺寸比间插值、未用最近值替代。保持原始 38 个间距节点、Fortran `globalerror=1e-8`。上游合并脚本的人类可读输出有 list.shape 错误，薄 IO 层只做数组拼接和列合并，随后调用原始 subtract-R2Binfinity 脚本。

科学局限：上游两球 scalar 间距插值和渐近/中场分界仍是来源方法误差；没有把离散硬球方法当成 Stokes 方程精确解。原生小间距插值在最低 s'=2.00001 以下夹取；仅接触诊断涉及该正则化。目标—墙正式测试不靠移动墙面或逐间隙拟合来消除误差。RMBW 表的有限间隙 RR/TR 参考资格限制保留。

## 实测证据

`validation/TOOL_FREE_SPACE_GATE.json`、`raw/FREE_SPACE_FIRST.npz`、三个 `raw/NATIVE_CACHE_*_FIRST.npz`、`logs/PECNUT_UPSTREAM_TESTS_CLEAN_ARGV.log`。上游测试以清空应用参数的 Python 包装运行，单个 pytest 测试内部的 88 个两球 fixture 比较通过；首次 pytest 命令行被上游 settings 当算例编号解析的失败日志也保留。

CPU 后端为 NumPy/SciPy/Numba，远端 AMD Ryzen 7 7800X3D、BLAS 4 线程；RTX 4090 存在但未使用。没有为 GPU 改写算法。
