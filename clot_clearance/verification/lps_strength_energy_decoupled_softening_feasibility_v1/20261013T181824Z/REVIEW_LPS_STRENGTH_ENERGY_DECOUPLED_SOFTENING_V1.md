# LPS_STRENGTH_ENERGY_DECOUPLED_SOFTENING_FEASIBILITY_V1 科学审查

run_id：`20261013T181824Z`。服务器 `/workspace/y_clot/project`；分支 `migration-20261009`，实际HEAD `cc5ee77149f55e58d715bb5a4660c6467231a826`。本机选择性GitHub快照仅作为交付工作树，科研源以服务器为准。

**结论：已建立一个相对于原瞬断LPS，在保留未损伤弹性与固定软化起点时，可独立调节后续不可逆耗能的离线数学与数值候选。宏观强度、材料Gc、真实血栓实验均未认证，未移植生产C++。**

候选选择自由能路线A。完整数学定义见 [SOFTENING_FREE_ENERGY_DEFINITION.md](SOFTENING_FREE_ENERGY_DEFINITION.md)，历史演化见 [DAMAGE_EVOLUTION_DEFINITION.md](DAMAGE_EVOLUTION_DEFINITION.md)，论文逐式审查见 [YANG_2019_EQUATION_AND_MODEL_AUDIT.md](YANG_2019_EQUATION_AND_MODEL_AUDIT.md)。

| 341节点共同运动学协议 | 初始广义刚度 N/m | 首软化阈值 | 完全损伤参数 | 累积耗散 nJ | 连续解析 nJ |
|---|---:|---:|---:|---:|---:|
| instant_reference | 5115.000000 | 0.15 | 0.15 | 143.859377 | 143.859375 |
| soft_sc018 | 5115.000000 | 0.15 | 0.18 | 173.156623 | 172.631250 |
| soft_sc024 | 5115.000000 | 0.15 | 0.24 | 230.776817 | 230.175000 |

上表完整经历加载—卸载—压缩—再加载；所有节点位移可实现且真实保存，无整体平衡求解。广义反力并非单轴材料应力。两软化轨迹的分裂积分误差约0.3043%与0.2615%；全部路径的功账本相对闭合残差最大2.35e−12，两种误差不能混用。理论sc处的浮点剩余反力分别约4.43e−17、5.60e−18 N，首个逐位d=1帧晚0.0005；该微小差异明确记录，不篡改门槛或实际帧。

初稿失败与修正：共用两端内部伸长产生3.49525646e−7 N的隔断后远侧响应。正式路径前在同一变分柔度机制内改为独立端点邻域状态，远侧响应变为0，原失败及代码保留。另有一次数值结果全部写出后的控制台numpy.bool JSON打印失败；只修正输出转换。交付阶段一次SSH文件传输中断导致一次报告封存缺少PAPER_SOURCE.json，补传后仅恢复元数据封存；无正式数值路径重跑。上述失败均保留记录。

**1. 实际读取的Yang原文**

用户指定路径的正式出版14页PDF，路径 `/home/lzy/projects/An improved ordinary state-based peridynamic model for cohesive crack growth in quasi-brittle materials..pdf`；与此前已读桌面副本SHA256一致：`e4a99edd7f2b02cf41ce3c36074b0b5667f567097e194e749a0c3bf284bc75d7`。第3–6、11–13页逐页检查原始公式排版，全14页文字已阅读。

**2. Yang与现有LPS的关键差异**

Yang选用Madenci状态本构、不同径向归一化；软化起始后移除该作用对其他邻居状态的影响并转入cohesive力。当前采用固定离散m、1/r影响函数、部分体积taper与完整势能耦合。

**3. 为何不能直接乘损伤系数**

一般非均匀损伤下，旧扩展空间刚度H与对角衰减矩阵Q不对易，QH不再对称，不是标量势能的Hessian。即便另写加权势能，也必须对全部θ耦合项求导，不能遗漏交叉项。

**4. 候选机制**

一个按端点邻域定义的变分串联柔度机制，内部弹性状态局部极小化，无向键共享不可逆损伤。初稿内部状态分配错误已保存并修正，未开展多个本构竞赛。

**5. 完整损伤自由能**

U=Σ_i min_v_i[ v_iᵀH_i v_i/2 +Σ_b k_ib(1−d_b)(e_b−v_ib)²/(2d_b)]，d=0用v=e约束，d=1去掉惩罚。当前K−5G/3>0保证完整二次项非负。

**6. 对应内力**

由包络定理求各端t_ib，合成无向键t_b，再成对沿当前键方向组装F=−∂U/∂u；代码与固定损伤差分核对。

**7. 未损伤响应保持**

是。当前最终两个状态对实际原生库最大相对L2力差6.196e-16，与既有独立能量相容。保存历史5键断裂后的固定拓扑参考也恢复。

**8. 软化开始条件**

实际几何键拉伸历史κ=max(0,历史最大e/r)超过s_onset后d连续增加。不是先验裂尖、宏观峰值或自然裂纹起裂认证。

**9. 是否保留0.15**

是；K、G、rho、dx、horizon、alpha及0.15的历史工程身份保留。没有以6.78kPa换算Yang的ft或s0。

**10. 新增耗能参数**

sc，两个预注册值0.18、0.24；k_ib由冻结完整H_i和参考r派生，不是额外材料调参。

**11. 对软化区间的影响**

κ从0.15到sc为渐进损伤区间。在本轮均匀等比例点阵路径上可解析得到线性下降，连续总耗散U0(s0)sc/s0；一般非仿射状态仍耦合，不保证每条键独立双线性。

**12. 不可逆历史变量**

明确保存每条原键κ及d，原键身份和原始参考m不重建。所有实际坐标、内力、θ、d、κ及功账本均FP64保存。

**13. 卸载愈合**

没有。四条路径最大d回退为0；零位移可回收能为0，但已累积损伤耗散不回退。

**14. 再加载历史**

遵循原κ和d，未越过此前历史峰值时恢复降刚度路径，不恢复完整刚度。压缩不新增拉伸损伤，原损伤保留。

**15. 固定损伤力与梯度**

候选独立梯度通过，最大绝对误差8.415397e-16 N、相对全局力尺度1.703861e-09，最大位置atom 269、axis 0，差分步长5e−10及1e−10 m；该最大误差来自5.00e-10 m。原生只作为未损伤和旧二值拓扑参考，原生库尚无部分损伤，不能声称原生软化也已通过。

**16. 损伤增加是否升能**

没有。固定坐标的泛函随d单调不增加，极小值保序；混合变形、受压、受拉有限增量和损伤导数均检查。

**17. 负耗散**

四条路径无负损伤增量，原始带符号最小值0；没有绝对值、裁零或丢弃事件。

**18. 完全失效残余拉力**

完全失效键实际力严格0；全体损伤后力与能均0。2470条跨平面键人为完全损伤诊断中跨力0，单侧变形对另一断块的内力影响0。实际部分平面路径仅断49条，不能称平面已完全分离。

**19. θ与固定m处理**

m、原拓扑及taper永久冻结；θ_el,i=L_i v_i，从各端弹性状态求得，随损伤连续松弛。完全断键可能保留内部松弛变量，但不参与实际跨面传力且断块能量可分离。

**20. 固定起点能否改变耗能**

是；瞬断/0.18/0.24的实测总耗散分别143.859377、173.156623、230.776817 nJ；记录ENERGY_CONTROL_DEGREE_OF_FREEDOM_DEMONSTRATED。

**21. 局部峰值如何变化**

本轮均匀扩张协议的峰值状态最大真实键力均3.893799733 µN，总广义反力均38.3625 mN，初始广义刚度5115 N/m。Q=−ΣF·(X/dx)，不是单轴边界拉力或材料应力。一般非均匀变形不承诺峰值不变。

**22. 宏观强度是否保持**

未认证。无缺口或缺口边值试样的整体峰值还需独立低成本验证；本轮只冻结局部起点并证明特殊受控路径可辨识。

**23. 能否达到6–8J/m²**

本轮有限平面开口范围不能；对应单调外功上界仅0.02129468 J/m²。没有把均匀全体损伤总能量任意除以裂面面积，也未计算其他路径的6–8J/m²参数。

**24. 具体限制**

当前线性LPS来自小应变映射，0.15起点本身未认证；最大0.245只是宽松数值诊断上限。有限horizon、离散几何和可控开口功预算均限制本协议耗能；缺少可信大变形本构与自然裂纹面积。不得外推极大sc。

**25. 新数字是否材料Gc**

否。均匀扩张仅报告nJ；附加部分平面开口0.00275878556 J/m²仅名为PROTOCOL_DISSIPATION_PER_AREA。材料Gc=null。

**26. Yang三维GF能否直接用**

不能。Eq.11–13依赖论文力密度、软化转移、下降段能量与完整邻域几何，不是当前状态能量和实际离散协议的现成映射。

**27. 仍缺哪些条件**

宏观强度、自然裂纹能量/面积、加载步长与网格收敛、非仿射实际裂纹验证、接触及原生MPI/排序/restart/损伤提交都尚缺。Eq.15命名约定疑点已明确；论文关键公式无缺失。

**28. 实际新增代码**

analysis/softening_feasibility_v1/：common、prepare、model、checks、freeze、run、postprocess、figures、finalize及config；完整源码快照、CANDIDATE_SOFTENING_IMPLEMENTATION.py和source_additions.diff交付。

**29. 运行数量**

4条正式离线路径（瞬断+两个sc的同一路径，另同几何部分平面循环）；初稿与最终各1轮门槛，合计2个小型原生静态实例、4个小变形原生状态；1个初稿断块独立性诊断、1项离线能量预算。无正式路径重跑。

**30. 实际资源**

单线程；正式路径墙钟7.261410 s、CPU计时7.587506 s、峰值221.191 MiB；含初稿与原生检查的已记录峰值666.027 MiB。已测阶段墙钟总和23.412 s，不包含开发阅读、SSH传输及未计时短探针。资源细项与文件体积见RUN_RESOURCES和manifest。

**31. 必要可视化**

四组PNG/SVG：实际反力与耗散、真实点阵和损伤、完整功账本及历史、耗能控制与限制；已实际打开检查并修正一处注释和失败键对比度。静态OPEN_RESULTS.html可直接浏览。

**32. 生产库与历史保护**

现场再次核查：最终诊断库SHA 00b246ebe4224e648c27b74afb64fc22b9cea4ba6ffb23f7c40d96c37ca22ae9；生产库SHA d1b2f53d6b942d1ca9277a8b8cd9c195ce75ea6f616681ec21c311e4faf3b41a。全部21891个既有verification/results文件及30个源码/库/环境文件SHA未变。无commit、push、reset、stash、清理或环境永久修改。

**33. WSL文档是否真更新**

是，/home/lzy/projects/clot_rupture.md已真正写入并逐字读回；保留全部旧正文和章节标题，新增本轮结果和旧大试样路线暂停说明。

**34. 备份diff与SHA**

本目录wsl_document/含before、after、逐行diff、DOCUMENT_UPDATE_AUDIT.json。修改前SHA c8653d48a343bb95c75319369fe9348530485415600978cfc89846ef1755ddc8；修改后SHA a2cdba711f63a3d775d1d22e836e3033b88e2a161d42f38377e8d4f42c9ebbe8；写入前再次核对并发修改。

**35. 是否值得进入隔离C++研究**

值得进入下一轮隔离原生开发的科研候选；先做同态力与能量、按原键ID损伤、排序/MPI迁移/ghost/restart、提交时序和接触验证。当前不能替换生产库，也没有自动启动下一轮。


| 独立验收项 | 状态 |
|---|---|
| PAPER_MODEL_INTERPRETATION | REVIEWED_USER_FULL_14_PAGE_PDF_WITH_EQ15_CONVENTION_CAVEAT |
| CURRENT_LPS_CONSTITUTIVE_COMPATIBILITY | RECOVERS_BASELINE_NEW_VARIATIONAL_DAMAGE_NOT_YANG_DIRECT_COPY |
| BASELINE_ELASTIC_RECOVERY | PASS_NATIVE_AND_INDEPENDENT |
| SOFTENING_FREE_ENERGY_DEFINED | PASS_ENDPOINT_STAR_VARIATIONAL_MINIMUM |
| FORCE_FROM_ENERGY_GRADIENT | PASS_FIXED_DAMAGE_FD |
| SOFTENING_ONSET_DEFINED | S_ONSET_0P15_ENGINEERING_THRESHOLD |
| DAMAGE_HISTORY_IRREVERSIBLE | PASS |
| NO_UNPHYSICAL_HEALING | PASS |
| NONNEGATIVE_DAMAGE_DISSIPATION | PROVEN_AND_NUMERICALLY_CHECKED |
| LOADING_UNLOADING_RELOADING | PASS_CONTROLLED_KINEMATICS_ONLY |
| FULL_SEPARATION_FORCE_BEHAVIOR | PASS_FULL_DAMAGE_AND_DIAGNOSTIC_CUT;PARTIAL_PLANE_PATH_NOT_FULLY_SEPARATED |
| DISCONNECTED_COMPONENT_LOCALITY | PASS_FINAL_ENDPOINT_STAR;FAIL_INITIAL_SHARED_STATE_DRAFT |
| GLOBAL_FORCE_MOMENT_CONSISTENCY | PASS |
| STRENGTH_ENERGY_DECOUPLING_FEASIBILITY | ENERGY_CONTROL_DEGREE_OF_FREEDOM_DEMONSTRATED |
| ENERGY_TARGET_ADMISSIBILITY | HIGH_ENERGY_TARGET_NOT_YET_ADMISSIBLE_FOR_TESTED_SCOPE |
| MACROSCOPIC_STRENGTH_NOT_YET_CERTIFIED | True |
| MATERIAL_GC_NOT_YET_CERTIFIED | True |
| EXPERIMENTAL_CLOT_FRACTURE_VALIDATED | False |
| NATIVE_CPP_UNCHANGED | PASS_SHA256 |
| PRODUCTION_LIBRARY_UNCHANGED | PASS_SHA256 |
| HISTORICAL_RESULTS_PROTECTED | PASS_21891_FILES_SHA256 |
| WSL_DOCUMENT_UPDATE | PASS_READBACK_AND_CONCURRENT_HASH |
| eligible_for_next_isolated_native_research | True |
| production_deployable | False |
| Gc_J_m2 | None |
| training_ready | False |
| native_fracture_runs | 0 |
| V3_runs | 0 |
| large_mode_I_runs | 0 |
| candidate_count | 1 |
| formal_path_count | 4 |
| scientific_design_correction_count | 1 |
| logging_failure_count | 1 |
| conclusion | YES: one offline, energy-consistent, irreversible endpoint-star compliance model demonstrates energy control at fixed elastic response and onset. No macro strength, material Gc or experimental certification. |
| VISUAL_REVIEW | PASS_FOUR_PNG_AND_SVG |

本轮完成一个候选与必要验证后停止，等待后续科研审查。
