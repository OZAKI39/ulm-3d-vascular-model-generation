# 后续模型扩展（本轮未实现）

先解决源点和终端身份，再讨论 Pries–Secomb apparent viscosity、Fåhræus–Lindqvist 效应，以及 hematocrit/phase separation。它们分别改变有效黏度、管径相关阻力和分叉血细胞分配，会使半径与流量反馈耦合，需要新的模型参数、适用范围和独立测试。

不能把这些模型和 outlet boundary 环境同时改变后，将分流差异全部归因于 ROI 裁剪。本轮一直固定 μ=0.00345312 Pa·s、ρ=1056 kg/m³，没有实现任何血细胞比容输运、非牛顿黏度或微泡模型。

±5%/±10% 是声明的测量敏感性幅度，不是从本 sample 实测推得的误差分布。没有找到可以用于逐节点/逐分支误差模型的重复半径测量或校准不确定度，因此标记 NO_EMPIRICAL_RADIUS_UNCERTAINTY_AVAILABLE。
