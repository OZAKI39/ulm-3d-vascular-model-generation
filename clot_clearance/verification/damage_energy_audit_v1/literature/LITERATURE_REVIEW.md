# 原始文献与固定版本源码对照

检索日期：2026-10-01 UTC。本地 35 次 replay 和 106 个独立审计任务完成后才开始文献调查，顺序依据见 `../comparisons/LOCAL_AUDIT_COMPLETE.json`。这里只比较方法，不实现新模型。

## Peridigm：必须区分具体 material 路径

仓库：<https://github.com/peridigm/peridigm>。固定 commit：`f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903`。许可证：BSD-3-Clause，原文保存在 [LICENSE.md](LICENSE.md)。仅保存六个相关源码文件，不引入整个仓库；原始文件头保持完整，下载 URL 和 SHA-256 见 `FETCH_BATCH1.json`、`FETCH_BATCH2.json`。

| 检查项 | 文件与行范围（固定 commit） | 实際结论 |
|---|---|---|
| 经典 K 与 F | [correspondence.cxx:810–920](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/materials/correspondence.cxx#L810-L920) | 有 K 与 numerator 的积分及右乘逆矩阵，但此函数不接收 bond-damage 数组；第 843–844 行是零值占位。不能把文本中的 `(1-bondDamage)` 误称为活动损伤权重。 |
| 经典 force、hourglass | [CorrespondenceMaterial.cpp:288–396](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/materials/Peridigm_CorrespondenceMaterial.cpp#L288-L396)，[correspondence.cxx:1252–1349](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/materials/correspondence.cxx#L1252-L1349) | Piola stress 经 K 逆矩阵转为力，另加 hourglass force；后者的 bondDamage 也为零占位。其公式不等于本项目 r/trace(K) 的能量项。 |
| bond-associated 重构更新 | [BondAssociatedCorrespondenceMaterial.cpp:1066–1110](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/materials/Peridigm_BondAssociatedCorrespondenceMaterial.cpp#L1066-L1110)，[correspondence.cxx:3040–3250](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/materials/correspondence.cxx#L3040-L3250) | damage 增长触发 gradient-weight 重算；有效 influence 含 1−bondDamage，参与 moment matrix 和 gradient weight，随后重算 F。与本项目共同点是损伤会影响重构；算法与本构能量并不相同。 |
| weighted volume | [correspondence.cxx:1378–1475](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/materials/correspondence.cxx#L1378-L1475) | 同时存在 damaged / undamaged 两种函数；不能因为存在某函数就断定所有 material 都调用它，也不能把它等同本项目不含 influence 的 remaining。 |
| 不可逆断键示例 | [CriticalStretchDamageModel.cpp:118–207](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/damage/Peridigm_CriticalStretchDamageModel.cpp#L118-L207) | 保留旧 damage，超临界伸长时设为 1，只有更大的 trial 值才更新。这是独立 damage model；不说明它与任意 correspondence material 的搭配成立。 |

所查经典和 bond-associated 路径未提供本项目形式 `sum(V*r*W(F))+Ustab` 的 fixed-x damage 单调性证明或检查。没有在这些路径找到保证，不等于证明整个 Peridigm/所有公开方法都没有。代码差异本身也不足以认定本项目有实现 bug。

[Littlewood et al., The Peridigm Meshfree Peridynamics Code](https://doi.org/10.1007/s42102-023-00100-0) 于 2023 年在线发表，卷期为 2024。已阅读出版社全文 §2.3–2.4：correspondence 使用传统 stress law 并需处理零能模态；state-based bond energy 与整个 family 有关，因此单键 stretch 不是通用的能量释放度量。全文通过浏览工具读取；直接下载的 `peridigm_2023.html` 是访问检查页，不能当论文全文。此论文用于方法背景，以上源码判断均来自固定 commit。

## Silling 原始 correspondence 与稳定化

[Silling et al. 2007, Peridynamic States and Constitutive Modeling](https://doi.org/10.1007/s10659-007-9125-1)：已读作者 M. Epton 上传的[全文展示](https://www.researchgate.net/publication/225366931_Peridynamic_States_and_Constitutive_Modeling)，§3 的式 (16)–(19)，§18 的式 (134)–(142)，印刷页 180–182。K 是影响函数加权参考二阶矩；F 为加权当前/参考混合矩乘 K 逆；再由本构应力生成力状态。这解释了本项目重构的理论来源，但这些式子没有证明任意改变影响权重会降低 W。UNL 官方库确认身份，[下载路径](https://digitalcommons.unl.edu/cgi/viewcontent.cgi?article=1094&context=usdoepub) 返回 403，OSTI PDF 返回 404；这里依据已实际可读的作者全文，不把下载失败标成成功，也未重构缺失公式。

[Silling 2017, Stability of peridynamic correspondence material models and their particle discretizations](https://doi.org/10.1016/j.cma.2017.03.043)：出版社文章信息已核对。全文公式对照使用同作者、同标题的完整 [SAND2016-6710 官方报告](https://www.osti.gov/servlets/purl/1457611)，报告 DOI [10.2172/1457611](https://doi.org/10.2172/1457611)。下载 32 页原文，已读 §4–6，尤其报告页 17–20；此处的公式编号属于报告，不冒充期刊最终版编号。

报告在对应能量之外加入非仿射残差 `z=η−Fξ` 的正二次能量，并由该能量求力；其特例 β=γω 同时使断键影响稳定化。报告页 20 还讨论了该选择在大损伤/离散化中的限制。证明的是针对形变扰动的材料稳定性；不能把它解释成对损伤变量的任意增量保证非负耗散。本项目 residual energy 与其思想相关，但 r、trace(K) 的损伤依赖和系数形式不同，不能直接借用其定理。

## Rausch 2021：直接核对全文，而非移植参数

[Rausch et al., Hyper-viscoelastic damage modeling of whole blood clot under large deformation](https://doi.org/10.1007/s10237-021-01467-z)，[NSF 官方全文](https://par.nsf.gov/servlets/purl/10317555)。本地 `rausch_2021.pdf` 13 页；已核对印刷页 1646–1648，§2.1.1–2.1.3，式 (1)、(11)–(14)，并目视核对第 3 页公式。详细演化依据由论文指向 Waffenschmidt/Ostwald 等工作；没有凭记忆补出他们的演化方程。

论文将运动学 F 与非递减局部内部变量 κ、独立非局部场 φ 区分。式 (11) 用退化函数乘平衡及非平衡能量；式 (13) 另外包含梯度项及局部/非局部场差异的惩罚项。不可逆性是内部变量演化的要求，不等于只将某个 r 降低。

**以下仅是本审计对所见式 (11)、(13) 的偏导，不是论文另一个被直接引用的演化定律**：固定 C、Fv、φ、∇φ，令未退化能量和为 Ψ0，则局部部分的共轭力为 `−f'_d(κ) Ψ0`；把惩罚能计入后为 `−f'_d(κ) Ψ0 + βd(φ−κ)`。所以任意手动增大 κ，也不能仅凭局部退化因子就宣称全耦合自由能必然下降。需要完整自由能、驱动力、不可逆约束和耦合演化的一致性。本项目 D 尚未具有这些依据，且 D 会改变重构 F。

这些文献足以完成本轮机制对照；未实现任何新方程、未引入 Ogden、未移植材料参数。后续至多两类候选见报告第 17 节，其可行性尚未验证。
