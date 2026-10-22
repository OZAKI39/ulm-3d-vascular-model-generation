# 三维与XY形状指标

固定ID、不重新匹配、不删除粒子。Q=YᵀY/N，Rg=√tr(Q)，主轴RMS长度=√λ，长宽比=√(λmax/λmin)。XY仅为投影，不是Chronos图像。邻近阈值d≤2.1R；硬接触几何统计d≤2R+1nm，二者均不代表融合。

主轴相对特征值间距≤0.03时标记 `ORIENTATION_UNIDENTIFIABLE`；全局F要求满秩、条件数<100、det(F)>0、非仿射相对残差≤10%。无效时保存NaN，仍保留协方差和成对距离。距离分位数为0/10/50/90/100%。

[分析定义源码](../../../scripts/straight_tube_exp1/f10ab/metrics.py)；[测试证据](evidence/TESTS.json)。

适用范围：`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED` · `training_ready=false`。Y 与四源静态结果另加 `HOMOGENEOUS_WATER_PROXY_SPECIFIC`。
