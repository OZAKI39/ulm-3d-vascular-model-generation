# F10Q 开发与科研交付

已在独立分支 `dev/fast-pzt-f10q-sequential-qp-xy-navigation-20261014` 完成两条独立路线。基线 `1c785f4803ca7dc0ea2fe91756b02487c7c75a59`；原F10P工作区未覆盖。

**主要结论：X原短目标达到；Y质心数值通过但目标未达到、形态与ID未收敛；旧强控制仍在0.884 ms拒绝。**

1. 开发前有界引用审计未发现新增可证明安全删除的对象，因此没有删除历史文件；清理审计与活跃索引单独提交。后续仅清除新绘图/分析代码中未使用的导入和变量。

2. F1–F10P报告、M0/PBD、B2、原硬球、R3壁面、F2/APECSS、四PZT缓存、失败状态、历史入口和不可变执行快照全部保留；关键源码SHA复核一致。

3. 实际0.464 ms拒绝提案支持有限切向球面曲率和冻结法向机制；并非仅凭公式猜测。后续0.884 ms的新失败属于活动集循环，不能混为同一原因。

4. ID 19/21的切向相对滑移9.261203 nm，法向位移−0.011488 nm，法向转动0.149052°；真实曲率间隙0.012046312 nm，二阶近似0.012046333 nm；原反力归一化后复现0.010056962 nm。

5. 每个物理步固定原x_n、F和阻力A，在候选位置重建真实距离与法向、重新线性化同一位移QP并更新反力；每次都从同一x_n形成x_new。独立检查局部KKT、真实法向受力平衡、互补、几何、修正及壁面，无重复推进时间或事后PBD。

6. 先保留三次QP仍不合格的冻结复核；依据残差连续下降，封顶增至五次。原冻结失败在第四次通过；新强轨迹每步1–5次，X完整轨迹1–2次，Y为1–4次。

7. 保持原F10G成对润滑阻力与过阻尼平衡，原声学/F2力仍在步初显式评估，步内只更新接触几何和反力。没有改水动力方程或引入弹簧。

8. 原17.8 nm重叠、178 nm修正、0.01 nm真实互补、局部数学残差、1 nm硬球诊断和原压力/壁面门槛全部不变。

9. 已越过0.464 ms，原初态重新积分的0.5 ms强控制通过；进一步的合法前缀为0.883 ms。

10. 没有完成强控制2 ms。0.884 ms提案出现active_set_cycle并拒绝，按登记顺序不再做强完整半步/细缓存；路线B照常运行。

11. 新强轨迹已接受前缀：最大重叠1.6547212e-08 nm，真实互补8.2430715e-09 nm，真实法向迁移率残差3.1093491e-15 m，局部QP残差8.2472954e-22 m，最大力残差1.1415031e-16 N。拒绝提案单独保留，不能计入合法轨迹。

12. 强控制未做完整2 ms步长/缓存对照，逐泡、质心及形态完整一致性均NOT_RUN/未准入；不能用局部一步质量或旧前缀比较替代。

13. 精确共同保存时刻为0与0.1 ms，无插值。M0/B2、M0/硬球、B2/硬球最大COM差分别约0.000012427、0.000005803、0.000018229 nm。

14. 形态不一致：对应Rg差约13.255825、2.509145、15.764970 nm，最大相对Q差约0.5453%；同ID跨模型差可达18.445099 nm。跨模型差异与各自dt/cache差分开记录。

15. 只能支持该初态、该强控制和两个共同保存时刻下团级平移很接近；不能证明一般模型不敏感，更不能外推2或30 ms。

16. 实际独立运行了8个新的中等X/Y候选，均完成0.5 ms；选择X_P1200和Y_P1800各做完整2 ms、半步与细缓存三条独立轨迹。

17. OFF严格零位移；X质心位移(0.200896324,-1.09787383e-06,≈0) µm，目标误差0.000896325125 µm；Y质心位移(≈0,0.051776285,≈0) µm，目标误差0.148223715 µm。目标均未改变。

18. X真正达到原0.2 µm目标：完整窗口、基础质量和COM准入通过，误差≤0.075 µm。Y未达到。成功范围仅为当前未校准模型中的短距XY数值原型，不是实验导航校准。

19. 有：Y是真实的团级位置通过而逐泡失败案例；dt/cache最大同ID差395.608440/764.640172 nm，超过原50 nm。

20. 有：同一Y案例位置通过而形状失败；Rg差96.417259/306.729687 nm，相对Q差4.0428%/12.1646%，匹配分布也超限。X各层数值检查均通过，X/Y各自初始团始终连通。

21. 未尝试5 ms或两段控制。X短距目标已达到，Y的形态/ID误差尚未解决；先解决这些问题比延长轨迹更有价值。5/10/30 ms均NOT_RUN。

22. 核心科学数值执行共139.215017 s：路线A 13.742350 s，路线B 125.472667 s；保留额度未耗用。远端四进程138.613354 s，峰值RSS 660.375 MiB。额外只读最终核验时间独立见VERIFICATION.json；不把研发/绘图耗时或CPU传输余量称为纯GPU时间。

23. 四组PNG/PDF在本独立工作区results/straight_tube_exp1/f10q/figures/，附绘图源码、FIGURE_PROVENANCE.json、逐图SHA、PDF渲染及VISUAL_REVIEW.json。

24. 进度文档备份为/home/lzy/projects/microbubble_aggregation.md.bak_F10Q_20261014；仅在原F10P段后追加中文F10Q段，保留全部历史、标题、公式和引用，before/after/diff/SHA及结构审计在results/straight_tube_exp1/f10q/progress/。

25. 下一轮先用可重放冻结状态解决活动集循环，并针对Y快速内部重排检查步长与场梯度敏感性；导航优化应加入实际接触/形态风险，再有限扩大X距离或探索较低压力Y控制。不先扩大长轨迹或声称任意形态控制。

## 分层状态

| 项目 | 强控制2 ms | X控制2 ms | Y控制2 ms |
|---|---|---|---|
| BASIC_NUMERICS_PASS | FAIL | PASS | PASS |
| GROUP_POSITION_NUMERICALLY_ADMITTED | NOT_RUN | PASS | PASS |
| GROUP_SHAPE_NUMERICALLY_ADMITTED | NOT_RUN | PASS | FAIL |
| GROUP_TOPOLOGY_NUMERICALLY_ADMITTED | NOT_RUN | PASS | PASS |
| PARTICLE_ID_NUMERICALLY_ADMITTED | NOT_RUN | PASS | FAIL |
| TARGET_REACHED | false | true (X) | false |

ENGINEERING_PASS由最终独立VERIFICATION.json确认。CONTACT_MODEL_DEPENDENT：历史形态对照成立，新控制跨模型依赖尚NOT_VERIFIED。EXPERIMENTALLY_CALIBRATED=NOT_VERIFIED；TRAINING_READY=false。持续保留EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED。

具体方法与证据见同目录TRUE_GAP_AUDIT、SEQUENTIAL_QP、GROUP_COMPARISON和XY_NAVIGATION文档；原始轨迹、每步QP反力/位置、质量表、源快照、运行日志与哈希均位于results/straight_tube_exp1/f10q/。所有图实际打开检查后才完成最终验收。
