# Particle-4：有限尺寸粒子接触

本目录新增 P4 模块；P0–P3 的源码、测试、科学参数和输入均不改变。P3 当前人工验收见
[PARTICLE3_MANUAL_ACCEPTANCE.json](reports/particle3/PARTICLE3_MANUAL_ACCEPTANCE.json)。
P3 原验证 JSON 保留原测试时刻的快照，因此其人工审核字段仍记录当时的 PENDING；新的独立验收记录为 PASS_WITH_CARRY_FORWARD_LIMITATION。

## 模型范围

支持 sphere / FREE_OBLATE ellipsoid / CAPILLARY_DEFORMED capsule 的全部六类配对。
粒子间硬接触无摩擦、无反弹，不使用质量、惯性矩、弹簧、阻尼、润滑或碰撞诱导变形。
输出的是速度修正和 KINEMATIC_CONTACT_MULTIPLIER，不是 force 或 impulse。
P3 的真实 RBC 通行仍为 NOT_ESTABLISHED；代理不可行不能解释为生理堵塞。

## 几何

`pair_geometry.py` 直接使用 P3 不可变 shape support。计算 Minkowski 差，GJK 得到分离方向，
再对真实光滑 support 求平稳方向，避免用近似多面体的面距离充当光滑粒子的间隙。
胶囊的直段 support 含绝对值：枚举所有端点符号区域、单个 kink 大圆和两个 kink 的交线。
椭球采用解析 support 梯度、Hessian 和切平面 Newton 修正；内部穿透采用固定 26 个球面种子及 GJK 方向搜索平稳极值。
这是经过独立数值参考验证的 support-stationary penetration 方法，不声称为任意凸体的通用 EPA 或形式化全局优化证明。
独立参考包括有限线段解析式、凸优化原始/对偶界，以及独立全局角度优化。未闭合的 witness 恒等式明确抛出 PairGeometryError。

令 n 从 j 指向 i，签名间隙为最大 support 分离平面间隙。外部分离时，point_i - point_j = gap*n；
内部相交时返回负的最小平移分离深度，测试涵盖负值，绝不伪造负数。
稳定 ID 决定算术顺序；交换调用输入时仅交换两个 witness 并反转法向。ID 不按空间位置重编号。
geometry roundoff 沿用 P3：512×float64 eps×实际坐标/形状尺度。没有物理 padding、固定比例 overlap 或扩大壁面 tolerance。

## 同时约束

`kinematic_contact.py` 将所有 WALL 与 PAIR 行一起按 canonical ID 排序。
对 FREE_OBLATE 使用变量 z=[V, ell*Omega]，ell=max(a,b,c)；球和胶囊只提供平移校正自由度。
Jacobian 的旋转列为 (r×n)/ell。目标为最小 1/2 ||delta z||²，满足全部接触法向速度非负。
确定性 dual active-set 解非负乘子，并检查 primal、dual、complementarity 和实际接触点速度。
数值界只来自 eps、矩阵尺寸、保留奇异值的条件数和实际速度尺度。
每个 PAIR 的平移贡献为 +lambda*n 与 -lambda*n；WALL 不要求粒子的修正总和平衡。
约束退化不相容、循环或 KKT 失败均报告 MULTI_CONTACT_INFEASIBLE。
CAPILLARY_DEFORMED 不接受角速度校正；球也不因法向接触改变保存的 spin。

## 真实时间积分

`particle4_motion.py` 复用 P3 的真实二分时间递归，左右区间都处理，失败尝试不推进时间。
每个接受状态检查全部配对和实际 WALL。初始配对重叠返回 INITIAL_PARTICLE_OVERLAP 并保留 IDs、接触点、法向和 gap。
一般固定模式区间采用双粒子的 surface-motion bound 与连续 separating-plane 证明。
转动椭球支撑函数 h=||diag(axes) R(t)^T n|| 满足
|h''| ≤ (amax²/amin+amax)|Omega|²，故端点弦减 B*dt²/8 给出全区间保守下界。
不能认证时才真实细分。改变形状模式或胶囊轴/尺寸时明确要求 P3 提供重新计算的路径；旧固定形状证书不能跨改变使用。
本轮真实 FEM 只运行双球，P3 的 RBC 变形更新和通行限制不被绕过。

孤立双球在持续擦边接触时，直线弦推进会造成不断离开/返回接触面的极短自由段。
为避免这种数值抖动，同一模型有解析积分分支：w=V_i_free-V_j_free，
n'=(w-(w·n)n)/R，c=n·w_hat 满足 c'=|w|(1-c²)/R。
解析积分保留实际初始中心间距，积分质心平移；c 达到 0 后恢复自由运动。它不将中心投影回接触面。
这一分支仅用于无 WALL、无其他粒子的孤立双球，并通过独立 solve_ivp 参考验证。
其状态速度是区间末的瞬时速度，区间平均速度另列 integration.interval_average_velocities_m_s；
其他分支记录该区间应用的常速度。projection 记录均针对区间起点。

## 候选筛选和验证场景

`pair_broadphase.py` 使用实际 support AABB 的稳定 sweep；swept box 只加运动上界和数值舍入。
保留 all-pairs 精确参考，所有最终接受状态仍检查全配对。耗时只作 INFORMATIONAL SMOKE。
五个 RBC 分位样本来自原 100000 个 P2 population，seed=2026092002；RBC 尺寸与面积预算没有更改。
人工四粒子场景使用两球及两个原 P2 椭球，中心线接触、零自由角速度、dt=0.2/0.1/0.05 秒，均为 VALIDATION_ONLY。
该场景姿态差为零不能用于宣称真实 FEM 姿态收敛。
真实双球使用原 SonoVue sampler 的固定种子 20260920 和 20260921。
在原 P3 接受路径上按索引确定性寻找正 WALL gap、正 pair gap 的初值；要求后续 256 个参考路径点也能容纳原较大 MB，第二个粒子不缩小。
最初只检查瞬时可行的候选随后发生大球壁面接触，P3 的时间细分产生大量子步，开发尝试被中止并保留日志。
正式 smoke 使用上述足够宽的参考窗口，不把它解释为大球持续沿曲壁接触的通用验证。
烟雾轨迹覆盖 256 个 P3 最细验证步，遇首个出口可提前结束。未到出口和无自然配对接触都是明确记录的有效结果。

## 复现

在仓库根目录运行，解释器为当前 WSL 已有验证环境：

```bash
P4_PY=/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv/bin/python
$P4_PY -B particle_3d/scripts/run_particle4_validation.py --stage geometry
$P4_PY -B particle_3d/scripts/run_particle4_validation.py --stage contact
$P4_PY -B particle_3d/scripts/run_particle4_validation.py --stage mixed
$P4_PY -B particle_3d/scripts/run_particle4_validation.py --stage real
$P4_PY -B particle_3d/scripts/generate_particle4_report.py
$P4_PY -B particle_3d/scripts/finalize_particle4.py --check
```

最终报告保持 PENDING_USER_REVIEW，不启动 Particle-5，不执行 CFD，不 push 或 merge。
