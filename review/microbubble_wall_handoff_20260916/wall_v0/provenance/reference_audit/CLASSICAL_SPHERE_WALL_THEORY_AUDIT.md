# 经典刚性球—平面壁面理论核查

本文件区分可读取的理论、软件实现和数值比较。壁面为 z=0，流体 z>0，球心高度 z_c=a+h，epsilon=h/a。所有转矩均绕球心，角速度采用右手系。只研究刚性 no-slip 球与 no-slip 平面；没有真实 STL、可变形气泡或壁面生产模块。

## 法向平移：Brenner 级数

Brenner 1961 原始 DOI 的全文未取得，记为 **UNVERIFIED_SOURCE_ACCESS（原刊全文）**：[原刊](https://doi.org/10.1016/0009-2509(61)80035-3)。可取得并逐页核对的独立来源是 [Ascoli 1988 Caltech 学位论文](https://thesis.caltech.edu/4428/3/Ascoli_EP_1988.pdf) 的 Eq.(29)，PDF 第 23 页；该研究直接复现并使用 Brenner 级数。归档保留 PDF、提取文本和公式页图像。

令 alpha=acosh(1+epsilon)，法向阻力/自由 Stokes 阻力为：

```math
lambda = (4/3) sinh(alpha) sum_{n=1}^infinity [n(n+1)/((2n-1)(2n+3))]
 * [(2 sinh((2n+1)alpha)+(2n+1)sinh(2alpha))
 / (4 sinh((n+1/2)alpha)^2-(2n+1)^2 sinh(alpha)^2)-1].
```

独立 `classical_wall_reference.py` 使用 mpmath 的 60 位和 80 位精度分别求和。每次连续 8 项低于当前和的 1e-32 才停止；14 个 gap 的两次结果转换 float64 后一致，epsilon=0.001 需要 802 项。该计算完全不调用两个软件包。近壁主导行为为 lambda~1/epsilon，即 M_perp/MT0~epsilon；不能以“单调下降但趋于非零常数”代替此渐近要求。

## 平移、旋转及耦合的可追溯来源

[Goldman–Cox–Brenner 1967 I 原刊](https://doi.org/10.1016/0009-2509(67)80047-2) 能核对书目信息和摘要，原刊完整公式 **UNVERIFIED_SOURCE_ACCESS**，没有凭记忆抄写完整双极坐标解。

可实际读取的 [O’Neill 1964 原论文](https://doi.org/10.1112/S0025579300003508) 讨论不旋转球的平行平移，PDF 第 7 页给出 F*、G* 及另一旋转问题的 F、G 表格。第 5 页定义 G*=G_y/(8*pi*mu*U*a²)，这里 G_y 是流体作用转矩。与外加负载阻力矩阵比较时，应使用 G*=-R[Omega_y,Vx]/(8*pi*mu*a²)，不能丢失这个负号。

另取表中 8 个 alpha 作为明确标注的次级诊断，z_c/a=cosh(alpha)，不把印刷表格中四位小数的高度当作精确输入。F*、G* 数值经公式页和表格图像核对；所有原始数值和输出存于 `raw/*_ONEILL_TABLE.json`。只比较平移列 F*、G*。Dean–O’Neill 原旋转表曾被 GCB 修正，因此没有把另外两列 F、G 用作准确性裁判。印刷表中远场 G*=0.00001 有明显有效位数限制，不能当作高精度耦合基准。

[Sprinkle 等 2020，Appendix A2 与 Table I](https://arxiv.org/pdf/2005.06002) 可完整读取，PDF 第 35–37 页明确列出各壁面阻力系数的渐近来源及切换范围。令 q'=(V,a*Omega)，阻力以 6*pi*mu*a 归一化，则表中给出：

| 模式 | 已核对的表达式 | 表中渐近区间 |
|---|---|---|
| Xtt（法向） | 1/epsilon - (1/5)ln(epsilon) + 0.9713 | epsilon<0.1 |
| Ytt（平行） | -(8/15)ln(epsilon) + 0.9588 | epsilon<0.01 |
| Ytr（交叉） | (4/3)[(1/10)ln(epsilon)+0.1895-0.4576*epsilon] | epsilon<0.1 |
| Xrr（绕壁法向旋转） | (4/3)[1.2021-3(pi²/6-1)epsilon] | epsilon<0.01 |
| Yrr（绕壁平行轴旋转） | (4/3)[-(2/5)ln(epsilon)+0.3817+1.4578*epsilon] | epsilon<0.1 |

论文将 Xtt 指向 Cooley–O’Neill Eq.(5.13)，将 Ytt、Ytr、Yrr 的近壁项指向 GCB Eqs.(2.65a,b)、(3.13b)。**Ytr 和 Yrr 的线性项是该论文作者对 multiblob 数据增加的拟合项，并非原始 GCB 的完整精确解。** 本次仅运行其冻结实现并核查出处，没有新增拟合。表中常数和源码常数保留不同印刷精度，例如 0.9713/0.971280、1.2021/1.2020569；不能拿印刷四位数设置机器精度 identity 门。

[Delong 等 2015 Appendix D](https://arxiv.org/pdf/1506.08868) 解释旧 `sphere/` 混合方法的组成。其 RR/TR 使用 162-blob 数据插值；这种描述不能证明接触极限精度。PyStokes 的 TT 自项与该附录 D1 的低阶壁面镜像近似一致：当 epsilon->0，法向归一化迁移率趋于 1/4，平行为 1/2，因此这类有限阶公式不提供润滑极限。

## 绕壁法向自转与绕壁平行轴旋转

[Liu–Prosperetti 2010](https://pages.jh.edu/aprospe1/publications/PapersPublished/Rotating/LiuJFM_2010.pdf) 的完整 PDF 可通过网页阅读器访问，已核查第 10 页 Eq.(4.1) 和第 13 页 Eqs.(5.1–5.4)。本地直接下载返回 HTTP403，原失败记录保留；**网页全文已读取**与**本地 PDF 未下载**分别记录。

Eq.(4.1) 给出绕壁法向旋转的阻力/8*pi*mu*a³ 为 zeta(3)-3(pi²/6-1)*epsilon+…，接触极限有限，zeta(3)≈1.2020569。绕壁平行轴旋转则具有对数发散。论文中轴命名与本任务不同，映射依据“轴平行/垂直于平面”，不是机械照搬其 x/z 标签。本审计只使用其 Stokes 极限内容，不使用有限 Re 数据。

另有一处必须保留的文献差异：Liu–Prosperetti Eq.(5.1) 的平行轴旋转常数为 **0.3709**，而 Sprinkle Table I / RMBW 使用 **0.3817**。两者共同支持 -(2/5)ln(epsilon) 主导项，但常数不同；本轮未取得 GCB 原刊完整推导来裁定差异来源，记为 **UNRESOLVED_CONSTANT_DISCREPANCY**。没有改写 reference 常数，也不把对 Sprinkle 表达式的一致性称为独立验证。有限 gap 的 RR_parallel 精度仍有此限制。

绕壁法向自转是**单球—平面壁面问题**；上一阶段未完成的 bubble–bubble center-line twist 是另一个问题，本阶段不改变其 pending 状态。

## 能支持与不能支持的结论

- 法向：独立 Brenner 级数可在全部 14 个 gap 作为数值参考；RMBW lubrication 最大相对阻力差约 0.8295%，epsilon=0.001 为 4.95e-7。
- 平行平移与交叉转矩：有 GCB 近壁公式出处、O’Neill 原始表格和本轮独立对照；不把软件之间的一致性当作唯一验证。
- RR_parallel：近壁主导对数项与作者公式可以追溯，但整个有限 gap 区间的高精度 RR 误差尚无独立数值解界定。
- RR_normal：有限接触极限有独立理论支持；中间 gap 仍依赖原有 2562-blob 表格。
- 远场：RMBW lubrication 在 epsilon>9.018296 使用低阶 bulk fallback，交叉项为零、RR 为 bulk。不能把这一区间的输出当作高阶 TR/RR reference。
- 正确的互易性、正性和量纲缩放是必要条件，不是任意 gap 精度认证。后续生产模块应从公开理论独立实现，并分别对照经典理论和合法使用的 reference 数值。
