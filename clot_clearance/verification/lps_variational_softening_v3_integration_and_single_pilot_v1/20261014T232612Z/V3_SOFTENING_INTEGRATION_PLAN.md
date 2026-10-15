# V3 渐进软化最小集成计划

定位：EXPLORATORY_APPLICATION；单张工程材料卡、一次从真实初态开始的正式 C1_GAUSSIAN_A200 物理工况。

实际3488节点=1927材料点+1561固定锚点，157637原始键=82270本体(type1)+22040界面(type2)+53327锚点间(type3)。全部原始图、ID、m、r0、taper保留；锚点不积分、不计脱离质量、不受背景阻力。

隔离修改：
- bond_bpm_peri_variational.cpp/.h：原代码仅允许一个type。增添按不可变原始键类型索引的s0/sc、严格共用K/G/horizon检查及typed restart；endpoint-star弹性能和局部消元保持原式。type1=0.15/0.18；type2=0.140625/0.16875（EXPLORATORY_UNCALIBRATED）；固定锚点间type3保留起点1.0，尾部1.2仅占位，所有锚点固定使其κ应保持零。
- fix_lps_variational_dynamic.cpp：维持POST_INTEGRATE提交损伤、拓扑删除、重新邻接及POST_FORCE统一SI转换。仅允许在转换之后的单个命名V3外力callback；禁止其他未审核修改力的fix。
- native.py / external.py：仅clot组原生nve；锚点不清零原始反力但不积分。复用真实Y壁查询/penalty、graph_components及drag_law，外力按当前物理time与当前x计算，采用已有COM阻力并按材料点分配。试步callback不提交不可逆脱离标签；标签仅在接受后登记，回滚重建callback状态。
- integrator.py / mechanics.py / journal.py：沿用Trial/Reject/Subdivide，按实际混合阈值预测Δd；保留真实原生restart与累计能量状态。账本分开人工载荷功、背景阻力功、附着数值阻尼、壁面弹性/阻尼和材料D。每步标量/稀疏损伤事件，关键事件完整FP64。
- preflight.py：相同V3初态以及一个小幅非均匀扰动，逐节点对比旧双修复LPS、独立完整参考LPS、新d=0力学；一个少量真实键的混合损伤、锚点、接触/拓扑/restart核查。必要短步只用于评估成本及接口。
- analyze.py / figures.py / report.py：真实附着/游离拓扑与离散体积直径，6组PNG/SVG、28项中文回答、单案例AI原型以及实际WSL更新。

复用：已通过的C++低秩核、1.0dx粒子接触、事件分辨积分、历史Y几何/无栓背景CFD、固定参考面面积/切向和C1载荷。历史输入精确幅值780.2322247780305Pa、预算200Pa、σ=0.28mm、T=1.1605ms、观察到2.33ms均不改；在变dt下必须使用真实累计time，不能用step×旧dt。

历史界面起点精确值为0.140624859375，与本轮明确指定0.140625有1ppm差异，必须声明，不能声称材料参数逐位相同。旧217/208结果含历史原生缺陷；另列双修复199/190结果作为更接近的参照，仅只读使用。

成本：现场约50GiB可用RAM、CPU负载约2；单rank/单线程共享原始图较合理，构建可用多核。延续用户取消墙钟/RSS固定停止上限的授权，不重设1800s/4GiB；正式前按短步成本估计。保留输出及Δd≤0.10、归一化能量≤2%、身份/运动/接触安全门槛。允许事件步细分，不扫描载荷或材料。

Gc、宏观强度及真实微泡streaming未认证；CRITICAL_LOAD=UNAVAILABLE。完成一次正式工况和交付后停止。
