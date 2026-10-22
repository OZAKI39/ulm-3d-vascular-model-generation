# 矩阵自由接触表示

每个约束行保留 i,j,n 三项；Bv=n·(vj−vi)，Bᵀλ 在 i/j 上累加 −λn/+λn。所有 M=N(N−1)/2 对都保存和检查，未裁剪远距对、未删除 F2、没有空间筛选假设。Bᵀ 是相同算子的真实转置；独立稠密乘积/伴随测试覆盖全部 N。

128 泡稠密 B=8128×384，24,969,216 bytes；PairMatrix 索引和法向为 325,120 bytes，另有 gap、pair mapping 和活动小矩阵。这是 **B 表示存储** 对比，不是整个仿真器 RSS 的 76.8× 缩减。阻力 A 仍为稠密 SPD；全系统力学并未变成稀疏近似。

工作集选取、全部非活跃不等式重查、原 primal/dual、真实最终法向、互补与力平衡均保留。

EXPLORATORY · PHYSICS_UNCALIBRATED · HOMOGENEOUS_WATER_CONTROL_PROXY · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
