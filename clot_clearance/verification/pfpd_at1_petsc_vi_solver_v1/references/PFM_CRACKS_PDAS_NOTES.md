# pfm-cracks / 2015 PDAS 源码对照

审计对象：[tjhei/cracks，paper-2015](https://github.com/tjhei/cracks/tree/paper-2015)，固定 commit `342678db32e1736b3f1eae601f55322907f584c9`。仅保存约 2.5 MB 的浅层源码快照供复核；未引入、编译或调用 deal.II。`cracks.cc:1–6` 明确 GPL version 2 or later；LICENSE 为 GPL v2 文本，故源码标识为 GPL-2.0-or-later。本项目没有复制其求解实现。

对应论文：[Heister–Wheeler–Wick (2015)](https://doi.org/10.1016/j.cma.2015.03.009)。以下实现结论直接来自固定版本源码，而非仅根据论文标题推断。

| 检查点 | 固定源码位置 | 观察 |
|---|---|---|
| 迭代初始化 | `cracks.cc:2285–2317` | 设置边界，装配全残差，保存旧解与上一次集合 |
| 分类与不可逆性 | `cracks.cc:2320–2376` | 相场节点由 residual/mass + c·(new−old) 的符号分类；受约束节点回到 old 值，Newton 增量置零 |
| 集合监控 | `cracks.cc:2377–2395` | MPI 汇总受约束数和集合是否变化 |
| Newton 和线搜索 | `cracks.cc:2397–2445` | 当前约束下装配并解线性系统；残差下降检查、回退和阻尼 |
| 收敛与循环限制 | `cracks.cc:2448–2463` | 必须残差足够小且集合未变化；最大迭代数超限抛出不收敛。未发现专门的周期检测或反循环记忆 |
| 输出 | `cracks.cc:2676–2678` | 输出 active_set 场 |

相同的结构：未知集合由迭代确定；不可逆性由历史解约束；约束自由度从 Newton 更新中固定；最终必须同时满足方程与约束。其 intactness 相场 φ 随破坏减小，本项目 damage s 随破坏增加，转换 s=1−φ 后历史上界对应历史下界。源码中的 active 指“约束生效”，本项目旧诊断的 active 曾指“正在演化的内部相场键”；报告与图分别使用 bound-fixed 和 free/interior，不能仅对照名称。

不能照搬的内容：FE 节点/单元遍历、质量对角归一化、hanging-node 约束、网格细化、FE 块预条件以及采用的本构近似；本项目是键相场和 BAQP 力。也不能把其 c 值直接用作本项目阈值。

[PETSc RSLS 官方说明](https://petsc.org/main/manualpages/SNES/SNESVINEWTONRSLS/) 与匹配本地版本的 [vi.c](https://github.com/petsc/petsc/blob/v3.25.3/src/snes/impls/vi/vi.c) 已提供边界符号分类、自由子系统 Newton 和线搜索，不必复制 PDAS。该版本 `vi.c:120–122,293–300` 确认下界 F>0、上界 F<0 属于固定集合；`vi.c:395` 及 RSLS 使用自由分量；`vi.c:600–604` 定义三个监控选项。`SNESGetVIInactiveSet` 返回约束不生效的自由集合。

# 2023 年效率参考的有限用途

[Kolditz–Mang–Wick, arXiv:2306.03501v2，§3.4](https://arxiv.org/html/2306.03501v2)：文中针对 Newton 残差已下降而少数约束自由度仍反复切换的成本，讨论调整集合判定及组合迭代，并另行迭代退化项线性化以接近单体极限。本轮没有实现这些改动。只有未来 RSLS 出现大量集合切换、重复集合循环或过高 Newton 成本，才值得做单独的效率研究；不能由本轮少量 C5 算例推断所有网格均无循环。
