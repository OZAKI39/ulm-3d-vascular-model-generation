# Particle-3：有限尺寸 WALL、硬接触与 RBC 胶囊代理

所有代码和证据限于 Particle-3；没有 CFD 求解、Particle-4、push 或 merge。
本阶段沿用 P0 场插值、P1 SonoVue 尺寸与入口起点、P2 的 100000 个 RBC 分布和 quaternion/Jeffery 实现。
P2 用户聊天验收单独记录在 `705dd0c`，不修改 P2 科学数据和实现。

## 数据合同

正式输入从 `formal_3D_flow_solver/FEM_SimVascular/frozen_reference/mesh_manifest.json`
及其 boundary manifest 解析。只有 `WALL.vtp`（ID 1，45221 个原三角形）是实体壁；入口与三个出口是开口。
WALL 刚性、静止、物理涂层厚度为 0。原三角形绕序给出外法向，取负得到内法向；18 个不同位置、方向、面积样本均对照唯一所属腔内四面体检查。
几何核心使用 m、s、m/s、1/s、m³、m²；图中才转换为 µm 等展示单位。

## 有限三角形间隙算法

`convex_triangle.py` 实现 `SUPPORT_FEATURE_STATIONARY_SIGNED_DISTANCE_V0`，是针对三种凸形状的有限特征距离算法，不依赖 GJK/EPA 库，也不以表面点云代替几何。
对单个三角形，gap > 0 为凸粒子表面与有限三角形的欧氏距离，gap < 0 为与该三角形相交时的最小平移脱离深度取负。
WALL 查询取候选三角形 gap 的最小值。它不声称给出同时退出多个重叠壁面的整体最短平移。
P0 的 center-inside-lumen 独立保留；三角形两侧的正间隙不能独自证明中心在腔内。

- 球：精确有限面投影和三条有限线段最近点，再减真实球半径。
- 椭球：令 `Q = rotation diag(a²,b²,c²) rotationᵀ`。支撑点为 `c + Q d / sqrt(dᵀ Q d)`。
  面的支撑方向解析求解；边投影成二维椭圆，顶点转成三维椭球最近点问题。
  枚举世俗方程 `sum(alpha_i*y_i²/(alpha_i+lambda)²)=1` 的全部实平稳根，包括不活跃极点的奇异轴解。
  在相邻活跃极点之间，该方程左侧严格凸，用导数零点隔离两个可能的根，随后括区间求根。
  所有候选都以实际有限三角形上的点核验，取最大的分离支撑值，得到外部距离或内部最小脱离深度。
- 胶囊：原三角形与负轴段的 Minkowski 和是最多六顶点的棱柱，精确查询其有限面，再减胶囊半径。
  棱柱退化为平面、胶囊轴段为零及面/边/顶点都覆盖。批量内核与逐三角形内核互相回归。
- VTK 仅做 BVH 候选选择。中心最近点的最终距离也由自己的有限三角形算法复算。
  BVH 距离搜索半径来自粒子包围半径和当前已认证距离上界，不是物理涂层。

固定舍入规则为 `512 * float64_eps * 实际坐标或几何尺度`，不设置经验接触距离。
细长三角形的有限特征判定在这个长度预算内验证最近点，避免把固定的无量纲重心坐标阈值误当空间精度。
半径搜索接受阈值使用更严格的形状尺度预算，为独立连续路径认证保留坐标尺度舍入余量。
任何阈值都不是 glycocalyx、润滑厚度或生物物理参数。

## 运动与接触

自由运动保持 `V_free = u_inf`。MB 自由角速度为半涡量；FREE_OBLATE 使用原 P2 Jeffery 角速度。
当接触点 `v_c = V_free + Omega_free × r_c` 向墙运动，只补偿中心法向速度；角速度不变。
球和 CAPILLARY_DEFORMED 使用平移约束。多个接触点用三维最小速度改变量的单边法向约束求解，最多三个独立活跃约束。
没有接触力、反弹系数、切向摩擦、弹簧或速度惩罚。

轨迹记录的 `free_velocity_m_s` 是当前接受位置的自由速度；`corrected_velocity_m_s` 是到达该行所用的区间速度。
因此验证位移用 `x[i] = x[i-1] + (t[i]-t[i-1])*corrected_velocity[i]`，约束的自由速度来自前一行。
初始行的 corrected 值是尚未施加时间步的自由速度。`normal_constraint_error_m_s` 对应该行刚完成区间起点的接触约束。
它是求解器诊断，已满足数值判据而未作修正的分支返回 0。另用 `audit_particle3_contact_velocity.py`
对每个真实接触区间的实际速度和所有接触法向独立复算，完整原始残差保存在 `04_real_contact_residuals.csv`。
最终机器记录采用完整复算的最大值。预算包含原有 QP 界 `1024*eps*velocity_scale`，以及原有
`256*eps` 近重合法向合并规则所产生的显式传播界 `||n-n_kept||*||V_applied||+|b-b_kept|`；两项逐区间单列。
这不是增加阈值：QP、法向合并和 WALL 的原阈值均不变，轨迹与速度也不变。首次审计只用 QP 界而漏计法向表示误差，失败证据保留在日志中。

## 真实物理时间与连续路径

失败区间 `[t0,t1]` 递归为 `[t0,tmid]` 和 `[tmid,t1]`；右半段从已完成的左半段状态继续。
只有通过连续认证的区间才能进入 ledger。不能移动中心来修复穿透，也不能原地消耗失败步时间。
最大深度 48 是实现保护，时间 ULP 下限同样显式抛出 `PhysicalTimeRefinementError`，保存区间、gap、三角形、模式和状态。

刚体用支撑分离平面或表面运动保守界认证整个区间。胶囊变化重新计算形状路径，不沿用旧刚体界。
在端点自由流方向之间，轴取归一化线性插值。定体积胶囊的连续支撑极值由不超过 12 次多项式的所有实平稳点、端点及绝对值折点共同检查；另有密集独立路径与中途最大支撑反例测试。

准静态模型在每个时刻选择最大的可行半径，不意味着线性插值两个端点最大半径。
若端点半径线性路径无法认证，算法寻找一个**仅用于证明存在性**的连续辅助胶囊路径：保持同一 RBC 体积、满足同一面积预算、使用同一轴路径，并在整个时间区间与 WALL 分离。
该路径证明每个时刻的可行半径集合非空；紧闭集合上的最大可行半径因此存在，并按定义满足墙约束。
辅助半径绝不写入接受状态；接受状态仍由完整有界搜索选择最大可行半径。
若找不到充分证明，则真正细分物理时间，而不把“证明失败”猜成“不可行”。
瞬时椭球/胶囊模式替换仅是准静态几何选择；两端替换形状和完整胶囊运动区间重新认证，不模拟膜展开的过程。

## 体积、模型面积预算与全区间搜索

原始 RBC 的 `a=b>c`、D、V、r 均继承 P2。`A_budget` 是原 oblate 的解析光滑面积，不是直接实测膜面积。
胶囊轴段长度 `L=V0/(pi*R²)-4R/3`，这里 L 是两个半球球心之间的距离。
面积 `A=2V0/R+(4/3)pi R²` 不超过预算；不强迫相等，未使用面积只代表模型未展开部分，不解皱褶。
在 `(0,R_equal_volume_sphere]` 中确定性求出面积允许的半径区间。

先检查原椭球的实际 WALL 冲突；只有冲突才启用胶囊。冲突只需一个认证的相交三角形，因此记录的
`original_oblate_conflict_gap_m` 是一个实际冲突证据，不保证是所有壁面中的最负值。
将有限三角形仿射映射到单位球坐标，可快速寻找相交候选；最终冲突符号仍由欧氏间隙核认证。
胶囊轴是归一化局部自由流速，数值不可定义时显式 `DEFORMATION_AXIS_UNRESOLVED`。

最大可行 R 使用最大半径优先的有界搜索；**不假设 R 越小越容易通过**。
对 `R∈[lo,hi]`，轴段长度 `L(R)>=L(hi)`，可用包含关系排除上端必穿透半径带；另一排除界来自形状 Hausdorff 变化的 Lipschitz 界。
剩余区间继续分解，直到全局上界与已有可行半径差在固定半径舍入预算内。
不能因为某个小区间自己的宽度小就丢弃它：若它仍高于当前最优解的允许误差范围，继续分解到排除成功或真正的可表示半径下限。
20000 次查询是实现保护；不能认证时返回 `DEFORMATION_SEARCH_UNRESOLVED`。
全区间被排除才返回 `DEFORMATION_SURROGATE_INFEASIBLE`，停止该验证案例。

搜索加速复用每个实际三角形的已计算 gap：固定中心/轴时，`|ΔR|+|Δ(L/2)|` 是形状 Hausdorff 变化上界，
旧 gap 减去这个界（再扣原坐标舍入预算）可认证该三角形在新半径下仍然分离。
未被认证的候选全部实际重算；一个真正相交的有限三角形足以排除当前半径。
这只加速可行性搜索，最终接受前仍重新查询整幅 WALL 的精确最小 gap 和最近点。
永久测试对照关闭缓存的完整搜索，以及直接遍历全部 45221 个三角形的结果。

CAPILLARY_DEFORMED 保留并冻结 P2 quaternion，胶囊几何只用局部流向。
`JEFFERY_ORIENTATION_NOT_INTERPRETED_WHILE_DEFORMED`；不能用冻结 q 的零差宣称真实姿态收敛。

## 验证与复现

Python 使用既有环境 `/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv/bin/python`，不修改那个历史工作树。
从项目根目录执行，建议设置 `PYTHONDONTWRITEBYTECODE=1`，使用 `-B` 和 pytest `-p no:cacheprovider`。

```bash
python -B particle_3d/scripts/run_particle3_validation.py --stage geometry
python -B particle_3d/scripts/run_particle3_validation.py --stage contact
python -B particle_3d/scripts/run_particle3_validation.py --stage deformation
python -B particle_3d/scripts/run_particle3_real_validation.py --stage controlled
python -B -m pytest -q -p no:cacheprovider particle_3d/tests/particle3 --ignore=particle_3d/tests/particle3/test_real_fem_evidence.py --ignore=particle_3d/tests/particle3/test_particle3_artifacts.py
python -B particle_3d/scripts/run_particle3_real_validation.py --stage mb
python -B particle_3d/scripts/run_particle3_real_validation.py --stage rbc
python -B particle_3d/scripts/audit_particle3_contact_velocity.py
python -B particle_3d/scripts/generate_particle3_report.py
python -B particle_3d/scripts/finalize_particle3.py --check
```

最终检查器运行 P0/P1/P2/P3 全部永久测试、原冻结工作树的 18 项 handoff 测试、冻结校验与 manifest 校验，并保存日志、JUnit、命令和 SHA。
各图的 JSON source manifest 保存数据路径与 SHA，图可用 `--only 2 --output /tmp/p3-plots` 单独重作。
真实 WALL 的局部控制试验保留全局间隙，明确 PATCH_ONLY；真实 FEM 轨迹始终用整幅 WALL。
128 个固定分布样本、9 个 48 边形管径、控制试验间隙比例 0.25、指定速度及 1 秒时长都是 VALIDATION_ONLY。
真实 dt 由 P2 已有最细验证步长及两次减半确定；真实终止预算沿用 P1。

## 证据解释

完整结果见 `reports/particle3/PARTICLE3_REVIEW.md` 和 `PARTICLE3_VALIDATION.json`。
开发失败日志也保留；只有最终完整回归结果可作为自动验收。
人工图像审核始终 `PENDING_USER_REVIEW`。
wall lubrication、near-field resistance、transit-time penalty 留给 Particle-5；RBC→MB displacement 留给 Particle-4/5。
不输出真实膜应力、接触力或压降，不重新求因堵塞改变的流场。生产步长未冻结。
