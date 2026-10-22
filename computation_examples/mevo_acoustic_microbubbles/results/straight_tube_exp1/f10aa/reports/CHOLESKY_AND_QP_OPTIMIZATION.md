# Cholesky、QP 和 SQP

原求解器一进入有约束路径就求解 A⁻¹Bᵀ 的全部 M 列；新实现仅对当前工作集新增列解线性系统，缓存限制在同一 QP 的同一 A/B 中，SQP 重线性化后重新建缓存。原目标、工作集排序、96 次上限、循环检测、gelsy 最小二乘/cond=1e-12、负对偶移除、违例添加和失败返回均不变。

Cholesky 没有跨物理步偷用旧因子，也没有显式矩阵逆。最大收益来自不再求解不需要的全约束右端项，以及不再每次构建稠密 B，不是把 Cholesky 本身宣称优化数倍。

FCC 532 活跃约束存在冗余/反力非唯一可能，物理验收采用净反力、KKT、真法向力平衡及位移，不要求任意数学对偶向量唯一。矛盾等式测试返回 active_set_inconsistent_equalities，非法几何、非 SPD、非法 SQP 上限均拒绝。

EXPLORATORY · PHYSICS_UNCALIBRATED · HOMOGENEOUS_WATER_CONTROL_PROXY · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
