# Particle-2 文献与分布合同

已知的直接统计：Moss 等的 WT C57BL/6Case 成熟 RBC 共有 1,156,720 个，
major diameter 均值 6.79 µm、SD 0.93 µm、众数及中位数 6.67 µm、IQR 1.33 µm。
同文 Table 4 的 47.9±2.6 fL 来自 4 只鼠的 MCV；2.6 不是单细胞体积 SD。
该来源是 **PREPRINT，非同行评审定论**。[Moss 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12132290/)

MCV 量级交叉核查：C57BL/6J 雌鼠 47.8、雄鼠 48.4 fL；另一个 C57BL6 control 为
49.3±0.9 fL，RDW 13.2±0.9%。这些研究不是同一个亚品系、样本或仪器。
[Rivera 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3656420/)，
[De Franceschi 2005](https://pmc.ncbi.nlm.nih.gov/articles/PMC1895196/)

V0 假设：直径用上述均值和 SD 的正态近似，体积用均值 47.9 fL、CV=14.8% 的正态模型。
14.8% 是用户明确选定、参考 C57BL/6 RDW 文献量级的模型选择，
不是直径数据集直接测得的单细胞体积 CV。原文核查发现 Moss Table 4 的 WT RDW 为
17.3±0.5%；本次如实记录这个差异，仍按已授权的 V0 选择 14.8%，不擅自更改模型。
我们没有可靠的 matched single-cell D/V 数据，因此先独立抽样，不能说真实 D/V 已证明独立。

三个标准差范围 D=[4.00,9.58] µm、V=[26.6324,69.1676] fL 是模型保护范围，
不是生物学 min/max。半轴 a=b=D/2、c=3V/(πD²)，厚度 2c、r=c/a 和 Jeffery λ
都由模型推导，不是文献直接测量。c≥a 的组合整对拒绝，不裁切、不修改 D 或 V。
拒绝筛选会改变最终分布，甚至引入 D/V 相关性；最终统计单独报告，不冒充原始文献值。

机器合同见 [C57BL6_RBC_GEOMETRY_DISTRIBUTION_V0.json](../../contracts/C57BL6_RBC_GEOMETRY_DISTRIBUTION_V0.json)。
文献元数据见 [sources.json](literature/sources.json)。没有在仓库复制整篇论文。
