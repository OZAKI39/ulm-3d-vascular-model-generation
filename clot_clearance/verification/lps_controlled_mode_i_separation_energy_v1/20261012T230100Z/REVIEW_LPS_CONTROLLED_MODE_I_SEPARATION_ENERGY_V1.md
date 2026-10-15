# LPS_CONTROLLED_MODE_I_SEPARATION_ENERGY_V1

Run：`20261012T230100Z`。

**本轮完成一次真实独立弹性预检与最终原生库静态核对；按冻结门槛停止，没有启动正式Mode-I断裂。实际等效分离阻力仍为UNAVAILABLE。**

上下均布20 Pa候选使保存端点的裂尖伸长领先于加载端，但最大伸长0.138705超过预检查0.03上限，最大残余力1.954×10⁻¹² N超过1×10⁻¹³ N。不得把优化器的success标记当成物理验收；也未放宽门槛或追加试样。零载初态通过，失败发生在已加载预检查。

| 指标 | 实际值 | 解释 |
|---|---:|---|
| 原始节点 / 唯一参考键 | 16000 / 909392 | 历史身份一致 |
| 初始预切 / 存活跨面键 | 16416 / 38992 | 均为参考拓扑，不是新断裂 |
| 预检外牵引 / 单面合力 | 20 Pa / 40 µN | 全节点自由，外边界加载 |
| 裂尖 / 加载端 s/s00 | 0.924699 / 0.824168 | 未接受弹性端点趋势 |
| 预检平均开度 | 140.609525 µm | 非规定开度、非物理时程 |
| 独立端点U | 2.508694009 nJ | 不是路径外功或裂面耗能 |
| 正式原生断裂工况 | 0 | 设计准入停止 |
| G_SEPARATION_PROTOCOL / GC_EFFECTIVE | UNAVAILABLE | 不能填0 |

本轮回答是：**尚不能测得与标称裂面能量相容的实际等效分离阻力。直接阻碍为选定牵引候选的预检幅值与静态平衡准入；没有证据将原因归结为材料太脆、实际非目标断裂或能量失守。**

## 34项科研审查

### 1. 实际使用哪个库？

实际创建1个原生静态核对实例，进程映射为 `/workspace/y_clot/project/verification/lps_native_state_identity_and_fracture_bugfix_v1/20261011T203015Z/build_owner_fix/liblammps.so.0`；SHA `00b246ebe4224e648c27b74afb64fc22b9cea4ba6ffb23f7c40d96c37ca22ae9`。最终双修复包含逐ID保存完整参考m及跨rank owned/ghost保护，源码未修改。本轮1 MPI rank静态检查，不冒充新的多rank断裂认证。

### 2. 原始点、键、m是否保持？

是。16000个实际节点、909392条完整唯一参考键、131196条跨Y周期键。原始NPZ SHA `ef145b6d6db6c17584dce6cbea1fa1b9d95847d71a69265a73df6334d1695c03` 已复核，逐ID m 相对差2.1921×10⁻¹⁵。先在完整拓扑初始化m，再明确预切16416键，剩892976键；没有从删键后的邻域重算m。

### 3. Mode-I几何？

X×Y×Z=5×0.4×1 mm；dx=50 µm，δ=150.5 µm；Y周期厚度0.4 mm，X/Z有限边界。z=0预裂纹从X负面进入1.50 mm，初始前缘X=−1 mm，剩余韧带3.50 mm。K=100000/3 Pa、G=940 Pa、ρ=1000 kg/m³、s00=0.15、α=0保持不变。这是理想化周期层，不是完整血栓实验。

### 4. 采用什么外部边界？

仅选一次对称均匀dead traction候选：上下各800个最外层节点，z方向各±20 Pa×dx²=±50 nN，单面合力±40 µN。所有材料位移自由度自由，Y盒固定周期；没有内部裂尖强制位移、额外夹具实体或局部材料强度保护。对称零位移初猜和自平衡载荷保留刚体平移自由度，没有新增物理点钉扎。

### 5. 如何处理夹具角点风险？

取消原外端面统一轴向位移，换成均布外力以允许轴向局部响应。历史未接受高载状态裂尖/夹具最大利用率约0.6364/0.9168；本次未接受端点变为0.9247/0.8242。只能说弹性伸长排序改善，不能声称已避免首断竞争：两组载荷类型和端点幅度不同，也都不是认证临界态。本轮不追加第二夹具或正式试验。

### 6. 初始态可信平衡吗？

零载A0可信：原生力、theta、U、K、速度和接触储能均为0，预切、存活拓扑和m一致。独立接触3.865×10⁻³⁷ J及4.32×10⁻²⁰ N是浮点距离舍入。20 Pa预检端点则未接受：最大残余力1.954014×10⁻¹² N，超过1×10⁻¹³ N；L2相对残差5.234×10⁻⁶虽通过另一项门槛，不能覆盖绝对门槛失败。最大伸长0.138705又超过预注册0.03的小变形预检上限。

### 7. 首次真实断键时间、位置、ID？

NOT_MEASURED。没有正式断裂时程，也没有真实新断键；静态原生核对使用break no且无积分器。参考键269977（原子4690–4851）只是保存端点最大伸长键，参考中点约(−1,−0.075,0) mm，不能写成首断键。

### 8. 是否从预裂纹前缘扩展？

没有这项观测。裂尖最大弹性伸长只给位置趋势，未进行自然断键及连续推进。

### 9. 非目标损伤是否发生？

静态固定拓扑未产生任何新损伤；动态竞争未测。加载端最大伸长0.123625，利用率0.824168，仍不可忽略。不能用它证明正式动态过程中会发生或不会发生非目标损伤。

### 10. 是否形成新生裂面？

没有新生裂面的静态证据。当前全部材料仍为一个16000点连通体；本轮没有动态新生裂面测量。有限裂纹可以绕前缘保留全局连接，因此单一全局连通分量不是唯一局部判据。

### 11. 多少跨面参考键断裂？

55408条原始跨z=0键中16416条是明确初始预切。静态检查新增拓扑变化0；正式原生断键计数为NOT_MEASURED，不能将预切键算成本轮自然断裂。

### 12. 多少跨面键仍活着？

38992条。原参考ID、跨Y最小影像、全部长键与存活掩码保存在完整精度NPZ；按100个X列及8个Y周期条带计数。

### 13. 新增面积如何定义？

候选X列宽dx=50 µm、完整Y周期厚度B=0.4 mm，每列Bdx=2×10⁻⁸ m²。必须从初始前缘连续推进，且该列各周期Y条带内全部原始跨面相互作用失去；另审查跨面长键的完整X投影覆盖。若交点判据与长键覆盖判据不同，保留前缘不确定带，不认证唯一面积。初始1.50 mm预切面积6×10⁻⁷ m²排除。

### 14. 是否只有部分损伤？

本轮连部分新损伤也没有计算；静态拓扑保持初始状态。交点分箱有30个初始空列；完整X投影筛查有28个空列，差异只是预切前缘的非局部影响范围，不能反向改写预切长度或制造新增面积。当前静态新增完全分离面积确为0；正式动态面积仍为null/NOT_MEASURED。

### 15. 真实载荷—开度曲线？

没有物理加载曲线。只有零载初态与20 Pa未接受静态端点，平均开度140.609525082 µm、单面外力40 µN。Fig4仅画分离的端点符号，不连接成曲线；1132次优化函数评估是求解迭代，不能用作时间。

### 16. 边界广义外功？

NOT_MEASURED。单时刻收缩ΣFext·u=5.6243810033e-09 J不是从零加载的路径外功。由于所有加载点法向力相同、方向相反，数学上全DOF功率等于Ftop×平均开度速度；此恒等已核对。仍须真实时程积分，不能拿终点40 µN乘开度当作实际输入功。

### 17. 独立LPS储能？

零载U=0，保存未接受端点U=2.50869400880696e-09 J。用已有力相容展开式和完整固定m；原生力与独立梯度相对差1.7302×10⁻¹³、最大绝对差1.484×10⁻¹⁹ N。原生bond PE×V=1.39692678308621e-09 J未当成完整U。没有物理U(t)时程。

### 18. 带符号拓扑能量？

NOT_MEASURED。保留定义ΔU_topology=U(x*,B_after)−U(x*,B_before)，但本轮没有自然事件，不伪造累积数值或取绝对值。初始预切不属于新裂纹能。

### 19. 动能、阻尼、接触？

两次原生静态坐标核对中速度/K/接触U均为0；没有积分器或物理时程。动态K、阻尼耗散、接触做功均为NOT_MEASURED。拟议动态候选无阻尼，并未执行；L-BFGS迭代不能称为物理阻尼过程。

### 20. 能量是否闭合？

未测，不能给PASS。已明确计划余额R=Δ(U_LPS+U_contact+K)+D_num−Wext−ΣΔU_topology，以及lagged-theta和时间积分误差；因为没有物理路径与Wext，任何断裂能余额都未被测定。状态ENERGY_COMPONENT_UNRESOLVED只表示缺少路径分项，不表示已经观察到能量不守恒。

### 21. 是否准静态或低惯性？

未认证。候选dt=0.5 µs、24 ms仅保留为未执行设计草案。剪切波速0.969536 m/s，H剪切穿越1.031421 ms，W穿越5.157106 ms；24 ms约为后者4.654倍。裂纹臂的柔顺弯曲时间尺度尚未核查，不能仅凭穿越倍数宣称准静态。

### 22. 可以报告分离协议能量吗？

不能。G_SEPARATION_PROTOCOL=UNAVAILABLE（JSON null），不是0。没有真实分离加载、面积及动力能量账本。

### 23. 可以报告实际扩展阻力吗？

不能。GC_EFFECTIVE_LOCAL和G_PROPAGATION_EFFECTIVE均UNAVAILABLE；没有把理想标称G0或静态储能直接改名。

### 24. 与0.06969623734234有什么关系？

上一轮G0_DISCRETE=0.06969623734234 J/m²及G0_THEORY=0.07170070833333 J/m²原样保留。这一轮没有实际单位裂面能量，不能计算比值、相容性或“偏高/偏低”。

### 25. 已证实的阻碍与推测？

已证实：20 Pa候选并非冻结的小变形预检状态；优化器因目标下降达浮点条件退出，但绝对残余力未达要求。源码/参考m/原生力与独立梯度身份核对通过。边界加载让裂尖伸长排序领先是端点证据；是否能自然优先扩展、动态惯性、面积唯一性、能量闭合及材料机制不足仍未测试。不能把全部失败归为材料太脆。

### 26. 支持单临界伸长耗能不足吗？

本轮不能支持新的实际耗能不足结论。理想名义能量低于Liu参考的数量级限制仍存在，但真实分离阻力未测。Liu人体总能约5.90、牛血总撕裂能约18.27、特定模型内禀反演约6–8 J/m²口径不同，不构成当前实验误差。

### 27. 足以要求软化升级吗？

不足。本轮未改变s00、G/K、δ、dx、接触或theta时层，也未实施Yeoh/Mullins/CZM/软化。先解决选定试样的静态准入，不能用更改本构绕过本轮失败门槛。

### 28. 实际新增多少算例？

1个独立弹性预检（1132函数评估）；1个最终库静态身份实例，含完整参考run0和2个无积分器坐标评估；正式原生断裂工况0、动力接受步0。没有第二夹具、第二低载、残差修正重跑或V3重算。

### 29. 资源和输出规模？

预检墙钟24.896568 s、峰值RSS 233.238 MiB；静态原生核对6.043792 s、峰值667.688 MiB。阶段合计与CPU/绘图/输出快照见RUN_RESOURCES.json；审查时98文件、52963758字节，不包括历史保护资产。完整最后规模由ARTIFACT_MANIFEST.json给出。服务器15.36核配额、42.07 GiB cgroup RAM；本轮单进程及单线程，未用GPU。没有把预检秒数当全任务耗时。

### 30. 库与历史是否保护？

最终双修复库、生产库、env.sh及历史资产均按HISTORICAL_PROTECTION_BEFORE/AFTER逐文件SHA复核。没有reset/clean/stash/push/rsync删除，未终止别的任务。监督器仅新增本任务预检服务，退出后不自动重启。

### 31. 生成哪些图？

5组必要PNG/SVG加1组明确未接受预检的补充图，共12图片文件：点阵边界、唯一初始接受态与缺失事件说明、跨面存活拓扑、静态载荷/储能与缺失路径账本、标称值对不可用实际值；补充图示真实未放大预检变形及区域伸长。没有虚构断裂帧、连续外壳、载荷曲线或MP4。每个显示子集的真实ID及数据SHA可追溯。

### 32. WSL文档已更新吗？

是。实际 /home/lzy/projects/clot_rupture.md 在⑤节追加本轮负面科研结论，并追加交付链接。①–⑨及全部历史文本保持；并发SHA检查、原子写入、实际读回通过。

### 33. 备份、diff和SHA？

本run的wsl_document/保存before/after/逐行diff及DOCUMENT_UPDATE_AUDIT.json，入口WSL_DOCUMENT_UPDATE.json。before `ba8b5875e71be49b06ddfc2b1783d4733a262d08d4244ce241c9fae156785780`；after/readback `151dd328a4a9ec629f7667c2d73ccf164a65ab96c6730b8b2f4c19929437df52`。文件清单与本地回收身份见ARTIFACT_MANIFEST.json、WSL_DELIVERY_VERIFICATION.json。

### 34. 下一轮唯一优先任务？

对已固定的同一分布牵引试样开展独立的无损伤静态准入任务：预先选择真正低载的小变形检查及可靠残差恢复方法，验证裂尖/加载端响应和最慢变形时间尺度；不迭代夹具、不变材料，通过后再单独申请自然断裂试验。本轮到此停止，不自动实施下一轮。

## 分项状态

| 项目 | 状态 | 范围 |
|---|---|---|
| SPECIMEN_IDENTITY | PASS | 原16000点、909392参考键及完整m，预切16416键逐ID核对 |
| PRISTINE_REFERENCE_M | PASS | 原16000点、909392参考键及完整m，预切16416键逐ID核对 |
| PRECRACK_IDENTITY | PASS | 原16000点、909392参考键及完整m，预切16416键逐ID核对 |
| MODE_I_LOADING_VALIDITY | CANDIDATE_NOT_ADMITTED | 对称外表面牵引候选；预检超出小变形范围且绝对平衡不达标 |
| INITIAL_EQUILIBRIUM | PASS_ZERO_LOAD_ONLY | 原生A0零载平衡可信；不等于20 Pa端点被接受 |
| NATIVE_LPS_LIBRARY_IDENTITY | PASS_ACTUAL_STATIC_INSTANCE | 最终双修复库实际进程映射及SHA |
| NATIVE_BOND_BREAK_VALIDITY | NOT_EXERCISED | 本轮仅静态检查，未开展真实断键 |
| CRACK_FRONT_LOCALIZATION | UNACCEPTED_ELASTIC_TREND_ONLY | 端点裂尖伸长高于加载端，不是实际起裂位置 |
| GRIP_FAILURE_COMPETITION | NOT_CERTIFIED | 伸长排序改善；未认证临界邻域或排除动态竞争 |
| CROSS_PLANE_SURVIVING_BONDS | PASS_STATIC_TOPOLOGY | 38992条跨面键存活，初始预切与新断裂分开 |
| NEW_CRACK_AREA_IDENTIFICATION | NOT_MEASURED_DYNAMIC | 静态拓扑未变，动态新增面积未测 |
| EXTERNAL_WORK_VALIDITY | NOT_MEASURED | 只有静态力-位移收缩，不冒充路径外功 |
| INDEPENDENT_LPS_ENERGY_VALIDITY | PASS_STATIC_ENERGY_FORCE_IDENTITY | 保留完整固定m的力相容展开式 |
| SIGNED_TOPOLOGY_RELEASE | NOT_MEASURED | 没有自然断键事件及其同坐标能量差 |
| KINETIC_ENERGY_AUDIT | STATIC_ONLY_DYNAMIC_NOT_MEASURED | 静态速度及原生接触为零；无动力轨迹耗散测量 |
| NUMERICAL_DAMPING_AUDIT | STATIC_ONLY_DYNAMIC_NOT_MEASURED | 静态速度及原生接触为零；无动力轨迹耗散测量 |
| CONTACT_ENERGY_AUDIT | STATIC_ONLY_DYNAMIC_NOT_MEASURED | 静态速度及原生接触为零；无动力轨迹耗散测量 |
| ENERGY_BALANCE | NOT_MEASURED_NO_PATH | 不能用端点势能推造外功与断裂能余额 |
| QUASISTATIC_SCOPE | NOT_CERTIFIED | 无动力K/U和时程认证 |
| G0_NOMINAL_REFERENCE | ACCEPTED_HISTORICAL_NOMINAL | 0.06969623734234 J/m²保持，非实际Gc |
| G_SEPARATION_PROTOCOL | UNAVAILABLE | 没有通过准入的正式分离及新增物理面积 |
| GC_EFFECTIVE_IDENTIFICATION | UNAVAILABLE | 没有通过准入的正式分离及新增物理面积 |
| MATERIAL_GC_EXPERIMENTAL_VALIDATION | NOT_VALIDATED | 没有材料固有Gc或实验认证 |
| PRODUCTION_LIBRARY_UNCHANGED | PASS | 10526个受保护文件SHA一致，缺失0、变化0 |
| HISTORICAL_RESULTS_PROTECTED | PASS | 10526个受保护文件SHA一致，缺失0、变化0 |
| WSL_DOCUMENT_UPDATE | PASS_ACTUALLY_UPDATED_AND_READ_BACK | 真实WSL文件并发检查、原子写入、完整旧文保留及读回通过 |

## 实施边界与复核说明

`PRECHECK_FREEZE.json`先于预检写入，包含唯一边界、幅值、求解方法、预算和门槛。阶段一失败后才写出的`MODE_I_DESIGN_FREEZE.json`明确标记“正式设计未完成/未启动”，不伪造先前不存在的动力时程或幅值。动态字段为null；候选0.5 µs / 24 ms仅作未执行草案。

只读后处理发现预检报告的有限X表面分区使用了节点中心±2.475 mm，而继承定义使用材料外界面±2.5 mm。已保存原标签并单独按原意重新分区，见SPATIAL_REPORTING_AUDIT.json；裂尖/加载端的领先键、全部力学状态与失败裁定不变，没有重算或后改验收门槛。初态原生/独立相对力差为1是分母仅浮点接触舍入的结果，绝对4.32×10⁻²⁰ N才是适当比较量。

当前库实际源码是语义依据。公开参考：[LAMMPS BPM/PERI](https://docs.lammps.org/bond_bpm_peri.html)说明state-based力和伸长断键；[fix addforce](https://docs.lammps.org/fix_addforce.html)给出外加节点力语义。此轮未在动力积分中使用addforce；外力只进入独立静态总势。Zhang & Qiao 2018与Liu原文身份、公式及口径继承上一轮经过SHA认证的资料，不重复声称本轮全文重读。Taylor DOI页面未能取得正文，没有据其摘要或未知内容增加结论。

结果入口：[OPEN_RESULTS.html](OPEN_RESULTS.html)。代码位于`analysis/mode_i_v1`；执行方法见[REPRODUCE.txt](REPRODUCE.txt)，完整浮点状态和原输入在本run目录。交付阶段的清单与哈希核验不启动任何新的求解。
