# PFPD_CANDIDATE_V1 人工审核报告

本轮结论：**DEVELOPMENT VERIFICATION PASSED**，限定为 BAQP 方法复现、规定损伤空间图案下的运动学/弹性能分离检验。完整 PFPD 损伤演化未实现；生产求解器与历史结果均未改变。

## 1. 为什么做这一轮

上一轮发现：即使把所有粒子的位置固定住，轻微增加当前侵蚀变量 D，保存的弹性能也可能增加。本轮检查一种具体解释及其修正原则：软化材料时，能否先保持“从邻居位置估算变形”的方法不变。工作是在已有工程中新增独立候选和诊断，并不产生新的侵蚀时间轨迹。

实际仓库根目录为 `/home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930`；项目在其中的 `clot_clearance/`。分支 `clot-clearance/energy-regularization-20260930`，已有 HEAD 为 `f0e891686b9b107510ab3292aa27ce83638b3339`。适用上级目录和项目内未发现 AGENTS.md。开始时建立了包含未跟踪成果的 5307 文件保护清单，而非只保护 Git HEAD。

## 2. 当前旧模型哪里有问题

旧 NOSB 模型的 `1−D` 同时影响材料软化、邻域形状矩阵 K 和重构的变形梯度 F。于是同一组实际坐标也能在 D 变化后被“读成”另一种变形；这部分储能增加会超过软化减少的储能。

本轮再现了全部 35 个旧跳步的能量、K/F 变化及上一轮的重构分解。旧模型每一跳的弹性能都增加，最大增加 `4.21051296176331e-18 J`；最大相对 F 变化为 `8.811990076664277e-9`。这不是仅凭动画外观或整体守恒作出的解释。上一轮已确认的稳定化项在这些跳步中是下降项，不能改说成它造成了增加。

## 3. PFPD 用什么思路解决

`g(s)=(1−s)²` 控制能储存多少弹性能；`h(s)` 控制邻居数据如何进入变形重构。`sc` 是两者开始分开的界限：当 `s≤sc` 时 h 恰好等于 1，只有超过 sc，运动学权重才下降。论文示例的 sc 为 0.95；本轮直接使用该值，并保留 sc=1（禁用运动学退化）和 sc=0（h=g）的预定诊断，不按结果挑选参数。

BAQP 把能量和应力放到每一条键的中点上计算。先平均两个端点估计的变形，再纠正该键方向，使估计与这条键的实际位移完全一致。一个物理键仅有一份状态；两个有向索引只用于端点邻域贡献和体积积分。

## 4. 论文有没有公式歧义

**有，而且损伤演化歧义没有解决。** 已用 pdftotext 提取本地原文，并检查 PDF 页 4、5、6、7、8、9、13、14、19 的渲染原图。完整逐式记录见 [公式审计](resources/PFPD_EQUATION_AUDIT.md)，身份见 [PAPER_IDENTITY.json](resources/PAPER_IDENTITY.json)。论文 SHA256：

`a77d29c88bd45bbd1cba78f1e823c09a9157b3be436b5eb6b3855c97a260a157`

γ=s、g=(1−s)² 和 Eq3.17 给出的固定变形约束极小值是 `max(0,1−Yc/Y)`；Eq3.25 却写成 `Y/(Y+Yc)`，低于阈值时也非零。两者不是文字解析差异，Algorithm 2 第16行也采用后一形式。仅画两式的比较图，不选择其中任何一个作为规范演化。[arXiv 当前页面](https://arxiv.org/abs/2604.09215)仅有 v1；原始网页快照已保存。

另外，Eq2.19 和 Algorithm 2 第20行把直接力项分母印成长度 r，作者公开 BAQP 源码则使用 r²。这里采用作者代码中的 r²，并通过作者可执行实现和独立能量有限差分核对。它解决的是力学实现依据，不代表损伤演化门已解决。h 改变 K/F 时，局部损伤导数还涉及邻域重构依赖，不能将固定 F 的驻值论证直接扩展到这一阶段。

## 5. 有没有找到作者公开 PFPD implementation

**PUBLIC AUTHOR PFPD CODE: NOT_FOUND。** 原文第19–20页称 PFPD 位于私有开发分支；本轮检查公开 80 个分支、16 个标签、16 个 release 和全部 156 个 PR，并检索 97 个本地远端/标签 refs 的源码、测试、文档和示例。所要求关键词中，只有 sc 命中接触刚度、表面修正等无关符号，没有找到 PFPD 演化代码。这是本次公开检索结果，不是对未公开代码的判断。

实际运行的外部库是作者的 [Peridynamics.jl](https://github.com/kaipartmann/Peridynamics.jl)，固定 commit [`b6e86563845a83b5cb06c805af19a817bde6b106`](https://github.com/kaipartmann/Peridynamics.jl/tree/b6e86563845a83b5cb06c805af19a817bde6b106)，版本 `0.6.0-DEV`，该提交无 release tag，MIT 许可证 SHA256 为 `0dc42bc9160a936d5fa5ec6f622e94d2baa1e02d13a675e9080795c156da483e`。可用的 BAQP 是 `RKCMaterial(monomial=:C1)`，不能把另一种 Chen–Spencer `BACMaterial` 混作同一算法。

Julia 1.10.12 仅装在本轮 external scratch；项目和 depot 独立，未修改系统。锁定 Manifest SHA256 为 `bb60d4987b8f59b7cc45eec6261f2471f92274f678871a0ae6763335bc06c808`。版本、源码定位、LICENSE、API 原始快照和每个 ref 检索证据在 resources/；命令见 [REPRODUCE.sh](provenance/REPRODUCE.sh)，Project/Manifest 小文件另存在 comparisons/，不需要保留庞大的安装缓存作为交付数据。

## 6. BAQP intact mechanics verification

完整无损 BAQP 先通过力学门，之后才做规定损伤诊断。Saint-Venant–Kirchhoff 用于论文方法及 Julia 比较；当前项目的 Neo-Hookean 另作兼容性检验。同一个比较使用相同材料参数和本构，不引入 Ogden。

|检查|结果|
|---|---|
|拉伸、压缩、剪切、等向伸长，内部与截断边界|Fbar/Ftilde 最大绝对误差 `1.5543122344752192e-15`|
|刚体平移、旋转|能量在约 `1e-36–1e-37 J`，最大节点力约 `1.4e-18 N`|
|非仿射变形旋转客观性|最大相对能量差 `8.83e-15`|
|全局节点力平衡|最大相对残差 `3.91e-17`|
|全局力矩平衡|最大相对残差 `2.54e-17`|
|完整坐标能量有限差分，五个步长、两种本构|全部步长中的最大相对 L2 误差 `2.363143867921908e-7`，最佳约 `1.40e-10`|
|与作者 Julia 的仿射和非仿射参考|最大相对能量差 `1.01e-15`，最大相对力差 `4.08e-15`，Fbar/Ftilde 绝对差 ≤`4.44e-16`|

有限差分步长为 `1.25e-8、3.75e-9、1.25e-9、3.75e-10、1.25e-10 m`，阈值和连续步长验收规则在运行前声明。另有不等体积、体积缩放、键方向精确性测试。完整记录见 [MECHANICS_GATE.json](baqp_verification/MECHANICS_GATE.json) 与 [OFFICIAL_COMPARISON.json](comparisons/OFFICIAL_COMPARISON.json)。

作者点能量采用单端归一化，论文使用端点对称核；二者全局能量一致，逐点能量的定义不同，因此比较全局能量。作者当前 weighted_volume 在运动学权重改变后也变化，本轮候选按论文 Eq2.14 保留完好 ω0；Julia 的等价性检验严格限定为无损状态。秩不足的重构明确报错，不替它补造碎片运动学。

## 7. c0 verification

当前二次核 `(1−ρ)²` 的分子是 `1/60`，分母是 `1/15`，因此 **c0=0.25**。手算和 scipy.integrate.quad 一致；常数核为 0.375，线性核为 0.300。见 [NORMALIZATION_CONSTANTS.json](normalization/NORMALIZATION_CONSTANTS.json) 和可选符号交叉检查记录。

这一步确定“核函数怎样把裂纹体积积分换成面积能量”的比例。它没有测出血栓的 Gc；所有裂纹能数值测试仅使用 `Gc_demo=1 J/m²`。

## 8. Griffith discrete normalization

在边长 3 mm 的规则立方 coupon 中放置一张中心平面，跨面键设 s=1，其余设 s=0。计算严格采用有向双积分：无向键求和时乘 2，即 `Σ(i<j) 2 wbar Vi Vj Gc γ/(c0δ)`，不额外除以 2。显式物理键求和与按晶格偏移量加速的计数实现经过独立小试样核对。

|每方向点数|spacing (μm)|δ/spacing=3.01|5.01|8.01|
|---:|---:|---:|---:|---:|
|20|150.00|0.984489|0.990889|1.000659|
|32|93.75|0.980393|0.987940|0.996534|
|48|62.50|0.978143|0.986316|0.994760|

表中为 `EΓ/(Gc A)`，目标是 1。九个有限试样最大误差 **2.185684%**，没有表面修正。固定 δ/spacing 加密时，自由边界影响减小，但晶格求积误差并不消失；某些有限试样更接近 1 是两种误差部分抵消，不能宣称它们最好。

额外的无限内部晶格计数中，δ/spacing 从 3.01 增至 5.01、8.01、16.01、32.01，比例依次为 `0.973699、0.983101、0.991227、0.997304、0.999258`；误差从 **2.630119%** 降至 **0.074222%**。因此本轮 **GRIFFITH NORMALIZATION: VERIFIED** 指解析常数、双计数约定、离散实现和已报告的收敛趋势通过验证，绝不表示所有粗网格都等于目标，或血栓网格客观性已解决。完整数据见 [GRIFFITH_STUDY.json](normalization/GRIFFITH_STUDY.json)。

## 9. historical 35-state compatibility

只读使用 `results/streaming_erosion_equilibrium_v1/` 的 20+10+5 个 fixed-x 跳步。坐标固定，不再求平衡、不推进时间。历史 D 仅作为 **prescribed spatial pattern proxy**，最大值 `0.0032794785449227204`，不称为 PFPD 历史 s。

所有能量变化在本轮采用 **new−old**；上一轮审计的 old−new 已显式换号。F/K 指标是在每个粒子上计算矩阵 Frobenius 范数变化，再相对旧值取最大值。

|模型|最大相对 ΔK|最大相对 ΔF|35 跳中弹性能增加次数|
|---|---:|---:|---:|
|A：原 NOSB，1−D 同时软化及重构|1.2014373e-4|8.8119901e-9|35|
|B：BAQP，g=(1−s_proxy)²，h=1|0|0|0|
|C：BAQP，同一 g，sc=0.95|0|0|0|
|D：BAQP，sc=0，h=g 诊断|2.3998920e-4|1.7623156e-8|35|

B/C 的 K、Fbar、Ftilde 不仅在误差范围内相等，而是逐位相同；二者的所有弹性能也一致。每跳弹性能变化范围为 `−1.11965286e-17` 至 `−2.79778410e-18 J`。D 的范围为 `+6.33323127e-19` 至 `+2.54352000e-18 J`，再现重构敏感性。A 与保存审计的能量差为 **0 J**，重构分项和 K/F 也完整再现。

A 保留原 NOSB 稳定化项，B/C/D 用 BAQP 键能量，所以 A 与 B 的绝对能量差不能单独归因于 h。更严格的受控证据是 B/C/D 使用同一 BAQP、本构、坐标和 g，只改变 h。每跳数组、绝对/相对 K/F 变化、软化和重构贡献均保存在 historical_compatibility/；[35 跳 CSV](historical_compatibility/ALL_35_JUMPS.csv) 可直接查看。

## 10. 是否消除了 reconstruction anomaly

**在本轮审计的固定坐标、规定损伤状态中，h=1 下已消除。** 当 h 不变，K/F 不再随微小损伤变化，未退化的键能量也不变；g 随损伤增加而下降，因此弹性能下降。sc=0 的重构敏感性及能量增加再次出现，是一个有控制变量的机制检验。

状态为 **REMOVED_FOR_AUDITED_STATES**。不能扩展成任意大损伤、自由碎片、不同邻域或动态裂纹中的普遍保证，也不能仅凭这组结果把旧侵蚀率升级为热力学损伤律。

## 11. 是否已经证明 thermodynamic PFPD

**没有。** 完整演化门为 **UNRESOLVED_IN_V1**；未实现规范驱动、历史推进和相场演化。两条相场关系仅保存在公式比较图/CSV，不作为候选 damage API。

历史 s_proxy 是外部规定图案，不能要求“退化弹性能+裂纹表面能”自动下降，更不能拿这种任意更新代替变分演化论证。c0 与平面裂纹积分验证也不等于整个 PFPD 的热力学证明。

## 12. 与 current production solver 的区别

新增 6 个 `pd_clot/pfpd_candidate/` 文件，仅提供 family、BAQP 固定状态力学、两种明确本构、g/h 与归一化。生产 `mechanics.py`、`fragment_mechanics.py`、`regularization/erosion.py`、`erosion_equilibrium_runner.py` 均未更改，也没有删原稳定化或把 D 换成 s。

所有新诊断数据和报告集中在本目录；代码、测试和执行脚本仅为新增文件。保护清单逐一核对 5307 个既有文件，生产修改 **0**、历史修改 **0**；证据见 [PRESERVATION_CHECK.json](provenance/PRESERVATION_CHECK.json)。未触碰 BraVa CFD、旧轨迹、兼容路径链接或旧停止任务。

算例使用已有直管项目的规则点云和保存状态；没有加入真实微泡 CFD、材料校准、墙面内聚律、疲劳或侵蚀/PFPD 耦合。约 8.18 GB RAM 下按至少 20% 余量限制为 3 个 Python worker；BLAS/OpenMP 各 1 线程，Julia 1 线程完成后才启动历史检验进程池。硬件与运行日志在 provenance/。

## 13. 测试结果

|阶段|passed|failed|skipped|日志|
|---|---:|---:|---:|---|
|本轮编码前原有基线|130|0|0|[BASELINE.log](tests/BASELINE.log)|
|新增独立候选测试|39|0|0|[NEW_TESTS.log](tests/NEW_TESTS.log)|
|最终合并测试集|169|0|0|[ALL_TESTS.log](tests/ALL_TESTS.log)|

原 60 条 warnings 仍存在，未通过隐藏 warnings 或删测试改变结论。A–J 验收证据由 [FINAL_VERIFICATION.json](provenance/FINAL_VERIFICATION.json) 汇总；本轮 Python/Julia 比较、历史重放、归一化与输出检查均已完成。未对图形效果调整任何历史科学数据。

## 14. 图怎么看

全部图同时提供 PNG 与 PDF；保持黑底、Arial 常规字重、英文、无箭头。横轴中的 35 个跳步是三组独立历史试验的排列，不是拼接的一条新轨迹。

1. [Energetic vs kinematic degradation](figures/energetic_vs_kinematic_degradation.png)：左图显示全部 s 范围，右图放大当前历史的小损伤区域。绿色阴影范围中 g 已下降，而 B/C 的 h 始终等于 1，解释了为什么早期损伤不改变重构。
2. [Normalization constant](figures/normalization_constant.png)：每种核的解析值和数值积分成对显示。二次核的两根柱都落在 0.25，说明常数核和项目核不能共用一个未经换算的 c0。
3. [BAQP affine reproduction](figures/baqp_affine_reproduction.png)：比较规定仿射变形与点/键重构的误差，分别列出内部和边界。所有结果远低于预先声明的容差，说明无损方法先通过基本力学检验。
4. [Historical K/F change](figures/historical_KF_change.png)：A 和 h=g 诊断有非零重构变化，B/C 都在零线上。纵轴用含零的对称对数刻度，零值没有替换成假小数。
5. [Historical elastic energy change](figures/historical_elastic_energy_change.png)：上方代表弹性能增加，下方代表减少，B/C 完全重合。图中明确标注 KINEMATIC / ELASTIC COMPATIBILITY TEST ONLY，不能读成总断裂能或完整演化验证。
6. [Griffith normalization convergence](figures/griffith_normalization_convergence.png)：左图是有限试样及内部晶格极限虚线，右图是增加邻域点数后的内部求积收敛。固定邻域比例的加密与增加邻域比例是不同操作，不能只挑最接近 1 的点下结论。
7. [Formula consistency](figures/formula_consistency.png)：AT1 驻值关系在阈值之前保持零，Eq3.25 曲线则立即上升。图用于呈现原文矛盾，明确标注 NOT A MODEL SELECTION，不用于选出“最好看”的损伤律。

## 15. 当前状态标签

|标签|状态|
|---|---|
|PFPD_CANDIDATE_V1|DEVELOPMENT VERIFICATION PASSED|
|PFPD PAPER AUDIT|COMPLETE|
|PFPD FORMULA CONSISTENCY|UNRESOLVED_IN_V1|
|PUBLIC AUTHOR PFPD CODE|NOT_FOUND|
|BAQP INTACT MECHANICS|VERIFIED|
|CURRENT KERNEL c0|VERIFIED_0.25|
|GRIFFITH NORMALIZATION|VERIFIED（限本节说明的数值归一化范围）|
|ENERGETIC-KINEMATIC SEPARATION|TESTED|
|CURRENT RECONSTRUCTION ANOMALY UNDER h=1|REMOVED_FOR_AUDITED_STATES|
|FULL PFPD DAMAGE EVOLUTION|NOT_IMPLEMENTED_DUE_TO_FORMULA_AMBIGUITY|
|PRODUCTION MODEL|UNCHANGED|
|PHYSICAL CLOT CALIBRATION|NOT PERFORMED|
|REAL MICROBUBBLE CFD|NOT USED|

既有状态继续继承：CORE SOFTWARE VERIFIED；MECHANICAL EQUILIBRIUM VERIFIED_FOR_CURRENT_DEVELOPMENT_CASE；ADR VERIFIED AGAINST SCIPY SMALL-CASE REFERENCE；STREAMING_EROSION RATE UPDATE SOFTWARE VERIFIED；STREAMING_EROSION SLOW-TIME STEP PASSED DEVELOPMENT CONSISTENCY STUDY；DAMAGE ENERGY AUDIT COMPLETED；ENERGY INCREASE SOURCE IDENTIFIED FOR AUDITED STATES。当前 D 的 THERMODYNAMIC DISSIPATION NOT ESTABLISHED；FAILURE/DETACHMENT NOT OBSERVED；MESH OBJECTIVITY 与 WALL ATTACHMENT STILL OPEN。

## 16. 下一步建议

唯一最小下一步是补齐完整演化的权威依据。若下一轮需要推进完整 PFPD 演化，能否提供作者后续源码、更新版论文，或明确允许按哪一个变分模型继续？在该选择明确前，本轮候选保持固定状态方法与兼容性诊断的定位，等待人工审核。
