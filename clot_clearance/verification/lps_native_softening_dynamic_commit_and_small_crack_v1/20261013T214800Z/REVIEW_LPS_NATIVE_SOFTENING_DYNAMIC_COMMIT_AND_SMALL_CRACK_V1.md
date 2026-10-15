# LPS_NATIVE_SOFTENING_DYNAMIC_COMMIT_AND_SMALL_CRACK_V1

运行目录：`/workspace/y_clot/project/verification/lps_native_softening_dynamic_commit_and_small_crack_v1/20261013T214800Z`。本轮严格先P0、再P1/P2、一个341点真实动力学工况，最后一次由实际风险触发的短dt/2窗口；P4未准入。

**核心结论：在0–122µs实际窗口内，原生Velocity-Verlet已产生不可逆自然软化和完全失效，自由节点真实运动，力、历史及动量一致，能量余额被独立重建。但主时程第244步的最大Δd=0.102554超过冻结0.10门槛，已停止；不能判为完整动力学验收通过。没有足够证据进入更大裂纹、宏观强度或材料Gc标定。**

## 运行前冻结与实际边界

G=940Pa，K=100000/3Pa，ρ=1000kg/m³，dx=50µm，horizon=150.5µm，s0=.15、主sc=.18，α=0。真实物理质量ρV=1.25e-10kg；原生BPM力密度由新fix在第二半步前转换为N。245自由点、左右49/47夹持点。完整加载预注册700µs，实际122µs后严格停止，没有临时增幅。图和表只用真实已保存数据。

短窗口105–122µs从同一接受restart出发，0.5/0.25µs对应34/68步。最大Δd由0.1025542483降为0.0516505045；能量余额−6.051331549pJ降为−3.002602740pJ，比值0.496188767；同一终点失效ID集合完全相同。最大坐标分量差1.226177nm、速度分量差0.000123467m/s。仅认证局部比较，不外推全局时间收敛。

## 完整能量账本

$W=\Delta U+\Delta K+D+\Delta U_c-J_c+R$。初态U=K=Uc=D=0，数值阻尼不存在。

| 项目 | 停止值 / J | 定义与计量 |
|---|---:|---|
| 实际边界功W | 5.0132696030981718e-09 | 边界实际位移×梯形材料反力，加精确边界ΔK；数值求积估计 |
| 完整变分U | 2.0399269792234265e-09 | 原生完整端点邻域自由能，独立重建 |
| 实际K | 2.5189077249078215e-09 | 全节点物理质量、接受整步速度 |
| 不可逆D | 4.6169616489359898e-10 | 同新坐标的U_old−U_new累计 |
| 数值阻尼 | 0 | 未设置阻尼 |
| 接触Uc / 接触拓扑注能Jc | 0 / 0 | 本窗口没有活动接触 |
| 带符号材料拓扑变化 | -4.6169616489359898e-10 | −D，仅展示，不重复累计 |
| 残差R | -7.2612659266755342e-12 | 真实未闭合部分，不改写为耗散 |

相对冻结累计能量尺度最大残差0.144631437%。为独立识别时间层误差，令d_n下的末态力为F_old(x_(n+1))：

- `Verlet defect = ΔK_free − Σ_free ((F_n+F_(n+1))/2)·Δx`。
- `Spatial defect = U(x_(n+1),d_n)−U(x_n,d_n) + Σ_all ((F_n+F_old(x_(n+1)))/2)·Δx`。
- `Switch work = 0.5 Σ_all (F_(n+1)−F_old(x_(n+1)))·Δx`。
- 接触缺陷按相同求积定义，显式扣Jc；其本轮只剩舍入量。

三个非接触项累计分别为+0.046242415809pJ、−0.011885524449pJ、+7.226909035287pJ；负和重建R，差2.8192e-23J。主要来源为损伤后的力切换，dt/2窗口约减半；这支持有限分裂误差的解释，不表示持续时间/速率/尺度独立的物理认证。

## 分项状态

| 项目 | 本轮状态 |
|---|---|
| DEPRECATED_CODE_AUDITED | PASS_28_CANDIDATES |
| PROVEN_DEAD_CODE_REMOVED | NO_PROVEN_SAFE_DELETION_0_FILES |
| CLEANUP_REVERSIBLE | PASS_NO_DELETION_RESTORE_UNNEEDED |
| ACTIVE_CODE_ENTRYPOINT_CREATED | PASS |
| HISTORICAL_EVIDENCE_PRESERVED | PASS_31001_FILES |
| PRODUCTION_LIBRARY_UNCHANGED | PASS |
| FROZEN_SOFTENING_MODEL_PRESERVED | PASS |
| DYNAMIC_INTEGRATOR_CONNECTED | PASS_ISOLATED_NATIVE_VERLET |
| ACTUAL_FREE_NODE_MOTION | PASS_245_FREE_NODES |
| PHYSICAL_TIMESTEP_IDENTITY | PASS_OBSERVED_244_MAIN_STEPS |
| ONE_COMMIT_PER_ACCEPTED_STEP | PASS_OBSERVED_WINDOW |
| RUN0_AND_TRIAL_READ_ONLY | PASS_AFTER_PYTHON_ASSERTION_REPAIR |
| DAMAGE_HISTORY_IRREVERSIBLE | PASS_OBSERVED_WINDOW_UNLOAD_NOT_REACHED |
| FIRST_SOFTENING_EVENT | OBSERVED_STEP170_85US_BOND3393 |
| FIRST_COMPLETE_FAILURE_EVENT | OBSERVED_STEP185_92_5US_7_BONDS |
| ORIGINAL_BOND_IDENTITY | PASS_SINGLE_RANK_DYNAMIC_AND_RESTART |
| DYNAMIC_FORCE_CONSISTENCY | PASS_ALL_245_SAVED_STATES |
| DYNAMIC_ENERGY_LEDGER | PASS_WITHIN_STOPPED_WINDOW_RESIDUAL_0_144631_PERCENT |
| MOMENTUM_CONSISTENCY | PASS_OBSERVED_WINDOW |
| TIMESTEP_SAFETY | STOPPED_MAX_DD_0_102554_GREATER_THAN_0_10 |
| POSTBREAK_CONTACT_CONSISTENCY | ZERO_FAILED_BOND_TENSION_VERIFIED_ACTIVE_CLOSURE_UNOBSERVED |
| SMALL_NONAFFINE_CRACK_VALIDITY | NOT_RUN_P3_ADMISSION_NOT_PASSED |
| CRACK_AREA_IDENTIFIABILITY | UNAVAILABLE |
| MATERIAL_GC_NOT_CERTIFIED | NOT_CERTIFIED_GC_UNAVAILABLE |
| MACROSCOPIC_STRENGTH_NOT_CERTIFIED | NOT_CERTIFIED |
| FULL_V3_NOT_RUN | NOT_RUN |
| PRODUCTION_SCALABILITY_NOT_CERTIFIED | NOT_CERTIFIED |
| HISTORICAL_RESULTS_PROTECTED | PASS |
| WSL_DOCUMENT_UPDATE | PASS_ACTUAL_FILE_UPDATED |

## 51项逐项回答

### 1. P0之前实际存在多少候选废弃代码？

受限入口审计共28个候选源码文件；“候选”不是已经证实废弃。包括当前两套analysis模块、两份冻结原生源码及CMake编译探针。详见CODE_CLEANUP_INVENTORY.json。

### 2. 哪些代码被确认已经无效？

没有文件同时满足“确定不被调用且无复现/历史价值”的安全删除条件。正式源和复现代码继续保留，编译探针用途不确定。

### 3. 实际删除了哪些文件？

删除0个；DELETED_CODE_MANIFEST.json路径列表为空。未通过改名、移动或删除历史数据制造清理成果。

### 4. 删除依据是什么？

白名单为空；实际根据限定Python imports、C++ includes、冻结manifest、构建依赖、Git状态及科研来源判断。无充分依据的文件不删。

### 5. 替代代码在哪里？

没有需要替代的已删文件。本轮新活跃入口analysis/dynamic_softening_v1/，独立原生源码在本轮source/src/BPM；冻结离线参考仍为analysis/softening_feasibility_v1/model.py。

### 6. 是否仍有不确定的历史代码？

有，清单中的UNCERTAIN记录保留，尤其CMake编译探针。未遍历历史实验全文作过度清理。

### 7. 删除是否可恢复？

没有删除，无需恢复。原始候选SHA、Git before/after及所有受保护历史SHA均记录；旧科学源快照完整保留。

### 8. 是否建立了精简的活跃入口？

已建立根目录ACTIVE_DEVELOPMENT_ENTRYPOINTS.md。P0版27行单独保留；最终入口补充停止状态与唯一后续优先项，总行数远低于150。

### 9. 是否减少了下一轮需要读取的代码范围？

入口明确当前四个C++文件、冻结model.py、动力学驱动和关键审计。后续可按需读取，不必展开旧PROJECT_ACTIVE_CONTEXT或全部verification。未虚构量化的上下文节省比例。

### 10. 是否保护了历史科研证据？

实际逐文件复核31001个历史文件和34项受保护资产，冲突0；已有tracked Git diff逐字相同。P0先于任何新C++实现。

### 11. 实际采用什么Velocity-Verlet调用路径？

锁定LAMMPS Verlet::run → modify.initial_integrate（free组fix nve两次半步中的第一半步与漂移、grip组fix move）→新fix POST_INTEGRATE→邻居重建→pair/bond→reverse_comm→POST_FORCE力密度乘V→final_integrate第二半步。245个自由节点三分量真实积分。

### 12. 损伤在哪个物理时间层提交？

在x_(n+1)已漂移、第二次速度半步尚未执行时。先用d_n计算同坐标U_old，再提交κ/d，计算U_new与带符号D增量；新力使用d_(n+1)。这是明确的损伤算子分裂。

### 13. 每个接受时间步是否只提交一次？

主时程真实step1–244与commit计数、token逐步相同；从170步restart续算计数不重置。短窗口另有68个真实步，原生整数step211–278，对应物理105–122µs。无force调用即接受的混淆。

### 14. run0和重复force evaluation是否保持只读？

是，实际重复run0、改变试探坐标求力和恢复均不改变κ/d/D/事件/提交token；有动态fix时显式外部commit被拒绝。第170步的Python断言曾错误比较不同坐标的U，失败记录保留；修正后同接受态续算，未更改C++或物理协议。

### 15. 是否发生真实自由节点运动？

发生。245个free点位置和速度由原生力更新，停止态相对参考最大自由位移23.355614862µm。独立核对每一步两次半更新方程，位置最大误差1.35637e-20m，速度误差5.46438e-17m/s。

### 16. 实际dt是多少？

主dt=0.5µs；初始完整切线给出dt*omega_max=0.07190393、线性临界dt=13.90744591µs。实际损伤门槛触发后，仅105–122µs一个窗口使用dt=0.25µs。没有dt/4或全程dt/2重跑。

### 17. 首次软化时间与原始bond ID是什么？

第170步、85µs，原始键3393。全部历史重建与原生κ/d逐值相同。

### 18. 首次完全失效时间与bond ID是什么？

第185步、92.5µs，同步首失效7条：2180、2194、3368、3393、3904、9197、10264。所有新增失效时的实际拉伸均达到sc=0.18。停止时520个唯一原始键各失效一次。

### 19. κ/d是否有任何回退？

已保存主时程及短窗口无回退，固定坐标D增量非负。卸载段没有达到，不能把此范围描述成动态全循环不可逆性验证。

### 20. 完全失效后是否仍有不合理非接触传力？

已失效边的两个端点标量t均为0，原始边非接触拉力为0；固定停止态梯度最大相对误差1.3468e-9。整个试样仍连通，不宣称本轮形成了完全分离的大块。

### 21. 接触何时启用？

d首次到1的POST_INTEGRATE内删除相应owned键和1–2屏蔽，随后同物理步邻居重建、pair计算可启用1.0dx接触。当前保存的自然张开时程没有实际闭合接触，活动接触步数0，因此动态闭合覆盖状态为UNOBSERVED。

### 22. 动力学速度与步间位移是否合理？

最大速度0.4443383054m/s，步位移/dx=0.0044333927，相对端点步位移/dx=0.0078199333，均低于冻结数值门槛。但最大Δd=0.1025542483超过0.10，故停止；这不是已认证生物加载速率或全程步长。

### 23. 是否存在数值异常冲量？

已观察步的实际总动量与含夹持惯性反力的外部冲量余额最大1.088211e-25N·s，未发现额外异常冲量。损伤后的正常力重分配和有限时间分裂误差另计。

### 24. 是否执行P4？

NOT_RUN。主P3损伤增量门槛未过，严格停止，没有创建或运行第二种试样。

### 25. 为什么选择该试样？

本轮唯一主试样复用历史341点/12047键，能直接对比冻结材料和原始身份。由于P3停止，没有对P4几何作虚构选择。

### 26. 实际预裂纹与加载边界是什么？

P4缺口与加载均未创建。P3没有新增预切键；原12047条键初始完整，左右x端49/47点规定对称位移，其余245点自由，夹持点y/z固定。幅值每侧60µm，升载300µs、保持50µs、卸载300µs、末保持50µs，冻结700µs；实际只到122µs。

### 27. 损伤是否从预期区域发生？

P4 NOT_RUN，没有预期裂尖可验。P3实际损伤分布见Fig3和FP64数据；该341点阶段只做动力学机制诊断，不追认它为缺口定位实验。

### 28. 是否出现真实裂纹扩展？

观察到真实非均匀损伤和原始键自然失效；没有证据将其认证为从预定缺口前缘推进的裂纹。P4没有运行。

### 29. 是否形成可辨识完整裂面？

没有建立可辨识新裂面，停止态为一个341点连通体。部分d、520条失效键或节点颜色不能代替裂面识别。

### 30. 是否发生非目标损伤？

P4未执行，目标/非目标带指标UNAVAILABLE。P3只保存实际分散损伤，不能直接用P4的目标裂尖语义分类。

### 31. 是否仍有跨面存活键？

有。作为明确标注的只读拓扑诊断，x=20µm平面仍有2041条存活跨面键；这不是预注册裂纹面积测量面。完整存活键11527条。

### 32. 是否能够报告面积？若不能，原因是什么？

新增裂面面积UNAVAILABLE，不能以0替代。P4没有准入、无已识别连续完全分离面；不用断键数乘dx²。

### 33. 实际边界外功是多少？

到122µs为5.0132696030981718e-09J，即5.013269603nJ。使用真实边界位移与反力的梯形功，加规定边界动能增量；是数值路径功估计，不是终点F·u。

### 34. 完整变分储能如何变化？

从0到2.0399269792234265e-09J。来自完整端点邻域自由能；逐态独立重建最大误差1.0960e-23J，未使用旧瞬断不完整PE替代。

### 35. 实际动能如何变化？

从0到2.5189077249078215e-09J，全部材料点含规定边界；独立sum(0.5*rho*V*v²)核对native KE，最大差6.2039e-24J。free点动能也独立存表。

### 36. 不可逆损伤耗散是多少？

累计4.6169616489359898e-10J，即0.461696165nJ；每步由同一新坐标的U_old−U_new实际计算，独立增量最大差6.2039e-24J。未按残差反推D。

### 37. 数值阻尼和接触贡献是多少？

没有阻尼，阻尼功恰为0；实际活动接触为0，Uc和接触拓扑注能为0。提取Ftotal−Fbond的舍入残量所积接触功约8.30e-28J，不把它称为实际闭合接触。带符号材料拓扑变化−0.461696165nJ与D只计一次。

### 38. 动态能量账本是否闭合？

在停止窗口按W=ΔU+ΔK+D+ΔUc−Jc+R核对，R=-7.2612659266755342e-12J；最大相对余额0.144631437%，低于冻结2%门槛。它是有有限误差的窗口通过，不能称整个700µs协议已通过。

### 39. 新模型是否出现未经解释的能量注入？

未发现未解释来源。余额主要来自损伤后力切换功+7.226909pJ，另有Verlet求积+0.046242pJ和固定d空间求积−0.011886pJ；其负和重建R，差2.82e-23J。这仍是有限时间误差，不宣称精确守恒。

### 40. 当前是否可以报告材料Gc？

不可以，UNAVAILABLE。旧G0_THEORY=0.07170070833333J/m²和G0_DISCRETE=0.06969623734234J/m²保留为旧瞬断模型的标称量；不等于新软化材料Gc。历史173.156623/230.776817nJ也不直接换算成本轮Gc。

### 41. 新模型是否已证明宏观强度不变？

没有。冻结未损伤刚度与s0不等于宏观峰值强度已认证，341点短动力学窗口不能替代材料实验和尺度验证。

### 42. 本轮是否有必要修改sc或进一步升级本构？

没有支持这种修改的证据。main sc=0.18保持；sc=0.24仅一个固定状态回归，未另跑动态。历史新旧二值力L2差8.814180703270e-6N原样保留；本次停止坐标同二值掩码差1.239512740840495e-5N属于两种自由能的本构差异。当前优先问题是时间分辨率。

### 43. 实际修改了哪些C++源码？

本轮独立副本中的bond_bpm_peri_variational.h/.cpp，加动态提交入口、同坐标缓存及完整全失效快速返回；新增fix_lps_variational_dynamic.h/.cpp负责hook、物理力换算和运行约束。pair、Verlet、NVE及生产源均不改。src/.gitignore仅补包文件登记。

### 44. 新隔离共享库SHA是多少？

866065394c48ebc8a6d32f147bde992a325c2a503a49d50392243eee10fea9bd。本轮全部正式步实际加载此SHA；FINAL_LIBRARY_IDENTITY.json还记录四个源码、最终补丁及旧生产库SHA。

### 45. 编译和运行总成本是多少？

P0 5.365s；首次编译与一次生命周期析构防护修正的增量编译合计40.101s，12-job上限。主244步3.138s，短窗口68步0.724s，正式动力合计3.863s、CPU 3.861s，峰值657.605MiB。主平均12.862ms/步包含观测与保存，不是裸内核性能。零步表达式失败、静态回归和图形费用在资源表分开；没有虚报未计时编排耗时。

### 46. 是否发生任何历史文件意外变化？

31001个历史文件及34项关键资产SHA无冲突，已有tracked diff未变。只有本轮新路径、活跃入口和明确授权的真实WSL文档发生预期修改。

### 47. 是否保留生产库？

保留。生产库SHA d1b2f53d6b942d1ca9277a8b8cd9c195ce75ea6f616681ec21c311e4faf3b41a；旧双修复库00b246...和上轮b80bfc...亦完全保持，env.sh未换库。

### 48. 是否生成必要科研图？

实际生成5组PNG+SVG及OPEN_RESULTS.html：代码结构、物理时间/提交/损伤、P3实际停止点云、完整能量、时间步与分项准入。P4未运行，因此第三组明确为P3停止证据，未伪造缺口推进或医学连续外壳。没有动画。

### 49. WSL文档是否真正更新？

实际修改/home/lzy/projects/clot_rupture.md，读回一致；SHA由61e0fe9947692ffea38272645758f38a994d900a56af00a177759e419bed476c变为b13d3f1a78874fb4ddf18633f674bcf5a97ad3941756c3729d8633e68a1666b8。保留原①–⑨结构、中文文风、材料参数和全部旧结果。

### 50. 备份、diff、删除清单和SHA在哪里？

本轮verification目录内：clot_rupture.before.md / after.md / diff、WSL_DOCUMENT_UPDATE.json、DELETED_CODE_MANIFEST.json、NATIVE_SOURCE_PATCH.diff、FINAL_LIBRARY_IDENTITY.json、HISTORICAL_PROTECTION_AFTER.json、PAYLOAD_MANIFEST.json。source_snapshot保留驱动/四个C++与读取的核心依赖，binary保存实际库。

### 51. 下一轮唯一最合理的开发任务是什么？

同一341点冻结工况的损伤时间步准入与全时程能量验证。先预注册完整较小dt协议，保持材料和夹具；不能把一个17µs窗口当作全程通过。通过后再审查P4，本轮不自动续算或开展更大裂纹/宏观标定。

## 复现、中断与停止证据

固定状态回归只使用intact、sc=.18部分、sc=.24部分、全失效四个状态，未形成第二动态族。C++最初编译后检查到LAMMPS析构时先释放force、后释放fix，增量修正析构防护，正式科学步均使用最终同一SHA。初次变量min/max语法失败有0接受步；改为数学等价ternary。第170步Python只读断言错误并非材料历史错误，INSTRUMENTATION_RECOVERY逐值证明同一接受态恢复；没有重放前170步。第244步Δd越界为真实停止条件。短窗口没有继续到更晚物理时间，没有第二窗口。

科学局限：本轮1MPI，无新增跨rank动力学迁移认证；旧2MPI静态/运动学证据保留。当前小图复制实现不证明大规模伸缩性。邻居列表本工况每步重建，不能把它写成任意邻居刷新策略都已覆盖。活动闭合接触、卸载段及P4均未覆盖。新旧二值本构差异保留，不把旧0.069696J/m²标称量用于新材料认证。

关键证据：DYNAMIC_DAMAGE_COMMIT_DESIGN_FREEZE.json、DYNAMIC_INTEGRATOR_IMPLEMENTATION_AUDIT.md、DYNAMIC_341_NODE_AUDIT.json、DYNAMIC_FORCE_ENERGY_AUDIT.json、STOP_STATE_TOPOLOGY_GRADIENT_AUDIT.json、DYNAMIC_ENERGY_LEDGER.json、DT_HALF_WINDOW_COMPARISON.json、RUN_RESOURCES.json、REPRODUCE.txt及PAYLOAD_MANIFEST.json。五组图见OPEN_RESULTS.html。

本轮开发及科学交付结束，等待下一轮科研审查；不部署生产、不运行完整V3、CFD、微泡或AI。
