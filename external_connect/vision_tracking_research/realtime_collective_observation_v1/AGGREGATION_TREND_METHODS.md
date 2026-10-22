# 聚集趋势方法与适用边界

这是二维可见分布观测，不是接触/融合识别，也不是 PZT 声场因果推断。

固定全通道及 300×300 px 图像支持区域，使用约 0.8 s 时间窗，至少 8 次、0.5 s 支持。每帧通过 cKDTree 计算最近邻和 kNN 密度，协方差给出 Rg、长短比、方向与形状张量。DBSCAN eps=55 / min_samples=4 标记候选区域；邻近中心的短时持续性只表示空间持续，不是物质身份。

收缩候选：NN 与 Rg 对数斜率均 < −0.025 s⁻¹。运动证据：去掉平均平移后的径向速度中位数 < −1 px/s，二维速度梯度迹 < −0.04 s⁻¹。两类证据同时成立才显示 VISIBLE_CONVERGENCE。纯平移被消去；旋转的径向与散度应接近零。

拒判：候选数时间窗极差/均值 >15%，焦点指标变动 >25%，亮度极差 >12，边界候选比例 >15%，候选或有效运动支持 <8。阈值是本轮冻结工程规则，不是已标定生物物理门槛。短时窗口跨越 onset 会有延迟，不能把延迟当作物理聚集时间。

分别保存 visible_distribution_change、local_convergence_evidence、aggregation_trend_status、observation_quality。低质量仍保留原连续指标，但不发布物理聚集概率。数量变化可以被观测，但不能因此自动认定聚集。

残余风险：恒定数量的候选替换、密集泡合并后面积仍落在筛选范围、光学叠加、出入视野、固定反射、非仿射 PIV 假向量，仍可能导致错误。当前真实录像没有独立二维聚集真值，更没有声学归因真值。

实现依据：[SciPy cKDTree](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html)、[scikit-learn DBSCAN](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html)。
