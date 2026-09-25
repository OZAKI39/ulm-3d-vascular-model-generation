# 冻结 stable lubrication 源码审计

目标版本固定为 22Jul2025 Update 6，commit `9c5ab448c78a14fd534619622162ba418d6a1fb1`。本审计没有修改上游源码，也没有使用 2Sep2026 candidate。原始源码、片段提取行号与 SHA256 均位于 `provenance/`。

| 模式 | lubricate/poly | lubricateU/poly | 独立数学核验及生产处理 |
|---|---|---|---|
| normal squeeze | IMPLEMENTED | IMPLEMENTED | XA 的奇异标量通过；保留相同标量数学 |
| tangential shear | IMPLEMENTED | IMPLEMENTED | YA 的奇异标量通过；保留相同标量数学 |
| transverse pump | IMPLEMENTED | IMPLEMENTED | 原生项不能直接作为耦合分解后的 pump；独立修复 |
| translation–rotation | PARTIAL | PARTIAL | 半径乘 YA 不能给出异径球的正确 YB；独立修复 |
| axial twist | NOT_IMPLEMENTED | NOT_IMPLEMENTED | PENDING_BLOCKER，未声称通过 |

`lubricateU/poly` 内存在多段阻力计算。首次提取器抓到背景应变路径而非完整 R_FU 路径；后续已改为定位声明 `a_pu` 后的完整系数段，并保留提取器编译失败日志。不能把应变路径没有 pump 的事实当成整个 pair style 没有 pump。

修正后的源片段 probe 使用 C++ 原文，不将重写后的 Python 公式冒充上游实现。100,000 组状态覆盖等径、异径、最小–最大径、中位–90% 径和多个间隙。结果见 `validation/COEFFICIENT_AUDIT.json`：

- stable squeeze/shear 对独立 XA/YA 标量的误差为浮点舍入量级。
- 原生 `a_pu` 对 i/j 交换不具备本任务要求的球对耗散对称性。
- 原生 `a_pu` 是单球 C11 类系数；shear 已经产生转动力矩，直接相加不能等同于独立耦合矩阵的 Schur 补 pump。
- `lubricateU/poly` 完整旋转段与 `lubricate/poly` 还存在系数差异：前者的 beta 因子实现必须单独核查，不能用另一文件的结果代替。
- `a*YA` 与独立 `|YB11|` 的异径球差异不是改变随机种子或调整容差可以消除的。

生产模块采用 project-owned `rigid_math.cpp`、`rigid_lammps_main.cpp`。normal/shear 标量保留；shear 的转动臂和 transverse pump 改为独立矩阵平方配方得到的系数，并使用明确记录的有限间隙客观性补全。完整公式、来源和近似范围见 `THEORY_DERIVATION.md`。

生产采用原冻结 CPU/GPU 静态库，`pair_style zero ... full` 仅请求 LAMMPS 邻居表。新增 embedding 与 custom storage fix 处理六自由度过阻尼解、RK2 更新和硬约束。原生 `lubricate/poly` / `lubricateU/poly` 不直接驱动正式近场算例。

## Twist 审计

冻结文档明确指出未包含 twist。Ball–Melrose 的作者仓储无全文、出版社返回 403，本次没有声称读到完整原文；其公开摘要、冻结文档中的模式表达式，以及 Radhakrishnan 独立推导已读取。Townsend 式 (26–27) 的轴向 XC 项包含非奇异集体阻力和远场贡献。与当前单球 Stokes 转动阻力如何匹配尚未通过独立验证，因此不能直接拼接，也不能编造一个 `A_tw`。

`TWIST_ROTATION=PENDING_BLOCKER`；`CASE_C3=BLOCKED_TWIST_NOT_IMPLEMENTED`；`FULL_ROTATIONAL_LUBRICATION` 不得标 PASS。其余必需项只有在系数、交换、耗散、矩阵、MPI、轨迹和真实旧 8 球重放全部通过之后，才允许声明 `PASS_WITH_TWIST_PENDING`。
