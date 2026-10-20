# 当前可达裂面与多涡旋剪切入口（20261017T175456Z）

当前任务 clot_shear_multivortex_fragmentation_v2；analysis/clot_shear_multivortex_fragmentation_v2/：shear_field_library → current_surface → native → run_cases → analyze。

科学标签 SHEAR_MECHANICAL_REQUIREMENT_ONLY。支撑单元边界、实际开口、外部通道连通和当前遮挡的几何代理已接入；仍有未表示表面，不能称完整真实流体润湿。材料/native内核不改。

原生状态=COMPLETED；等功={'counter_rotating': 'WORK_MATCH_NOT_ACHIEVED', 'checkerboard': 'WORK_MATCH_NOT_ACHIEVED'}；完整碎块目标={'counter_rotating': False, 'checkerboard': False}；优化准入=NO_AUTOMATIC_OPTIMIZATION_ADMISSION。

先读 verification/clot_shear_multivortex_fragmentation_v2/20261017T175456Z/FINAL_STATUS.json、REVIEW_CLOT_SHEAR_MULTIVORTEX_FRAGMENTATION_V2.md、DELIVERY_VALIDATION.json；图表 results/clot_shear_multivortex_fragmentation_v2/20261017T175456Z/OPEN_RESULTS.html。

历史A仅为预算参考；不追加正式工况或自动优化。所有接受证据与失败初始化记录均保留，新实验使用新RUN。

## 历史入口原文（保留）

# 当前随动剪切开发入口（20261017T165129Z）

当前任务 `clot_shear_load_inverse_design_v1`，入口 `analysis/clot_shear_load_inverse_design_v1/`：shear_field_library.py → current_surface.py → native.py → run_cases.py → analyze.py。

两组3.6ms小试样原生完成。SHEAR_MECHANICAL_REQUIREMENT_ONLY；182原始顶面随动、正1/3载荷映射和相对速度拖曳。新裂面润湿、未表示碎片表面流体力及旧V3局部映射仍未解决。材料/native内核不改，无压头/弱格子/人为删键。

两个冻结工况均未形成20–200节点凝聚性游离块。单向剪切历史最大块988节点，双涡旋2节点；本次双涡旋没有展示多目标碎块优势，尚未识别成功剪切需求。表面支撑、空间覆盖、相对速度和存活三维桥接限制推断；不支持盲目扩大BoTorch搜索。

准入：NO_USEFUL_TARGET_OR_NEAR_SEPARATION_SIGNAL；没有优化trial或完整V3。

先读 `verification/clot_shear_load_inverse_design_v1/20261017T165129Z/FINAL_STATUS.json`、`REVIEW_CLOT_SHEAR_LOAD_INVERSE_DESIGN_V1.md`、`DELIVERY_VALIDATION.json`；图表 `results/clot_shear_load_inverse_design_v1/20261017T165129Z/OPEN_RESULTS.html`。

所有新计算使用新RUN，已封存结果只读。

## 历史入口原文（保留）

# 当前声流／剪切开发入口（20261017T132706Z）

当前任务：`clot_acoustic_streaming_shear_multifragment_v1`。两组3.6ms原生小试样已完成，材料与native库保持不变。

**当前主线为流体剪切，不继续压头迭代。** 真实微泡Streaming尚未确定；可计算路径为明确标注的连续反向Stokes剪切模板及独立机械需求牵引。

- 入口：`analysis/clot_acoustic_streaming_shear_multifragment_v1/acoustic_streaming_shear.py`、`streaming_interface.py`、`shear_loading.py`、`native.py`、`run_shear_cases.py`。表面牵引、N单位原始ID传递、不可逆损伤与真实存活连接均实际执行。

- 试样1792材料点+336锚点；dx50µm，原95666键/2736预切保持。实际顶面非负三角权重；不新增预切、不改G/K/s0/sc。

- 模式：sin(2πy/700µm)沿x反向剪切，25µm间隙、μ1mPa·s、4mm/s规定上边界；2ms升、2.4ms降、2.8ms停驱、3.6ms结束。B只放大牵引至1000Pa，不代表微泡流体能力。

- stokes：峰值0.156890216Pa；末态游离0点，历史最大块0点，历史最多20–200节点块0，完整目标=False。

- inverse：峰值1000Pa；末态游离111点，历史最大块15点，历史最多20–200节点块0，完整目标=False。

- 声学/气泡输入已核对，但缺少当前血栓配准、间隙及近壁一阶/二阶流场。APECSS缓存黏度2mPa·s，与上述独立模板不同；声压、R(t)、质心速度均不能直接当作Streaming剪切。

- V3 3611实际CFD三角面仅作映射审计。旧带符号映射存在明显局部差异；非负核牺牲力矩。不能仅凭全局守恒认证局部传递。未运行新完整V3、FSI或优化trial。

- 参考面冻结，移动/脱离后仍沿原始ID加载；大变形阶段只支持死牵引机械探索。新裂面应力和当前表面重算尚未完成。

- 先读 `verification/clot_acoustic_streaming_shear_multifragment_v1/20261017T132706Z/FINAL_STATUS.json`、`ACTIVE_STREAMING_SHEAR_PIPELINE.json`、20问 `REVIEW_CLOT_ACOUSTIC_STREAMING_SHEAR_MULTIFRAGMENT_V1.md`、`DELIVERY_VALIDATION.json`；图表入口 `results/clot_acoustic_streaming_shear_multifragment_v1/20261017T132706Z/OPEN_RESULTS.html`。

- 复用库 `verification/clot_multifragment_breakup_v1/20261017T110849Z/build/liblammps.so`，SHA `465f8144c6da5c9de3aec4137b9c325f53c26737b199c536fb4a20845511ef77`。本轮没有更改生产本构、MPI或restart；新实验必须新RUN，保留本轮封存输入与结果。

- 后续先补齐微泡近壁流体来源与局部表面耦合条件，再讨论微泡参数反向优化和V3。不得把损伤前兆、单大块或细粉当成多个凝聚性碎块。

## 历史入口原文（保留，不代表当前加载主线）

# 当前clot原生软化开发入口

当前：lps_variational_softening_v3_integration_and_single_pilot_v1 / 20261014T232612Z。

- 正式性质EXPLORATORY_APPLICATION；状态**FULL_PROTOCOL_NUMERICAL_PASS**，实际到2.33ms，计划2.33ms。一次正式case，从真实V3初态开始；没有第二材料卡/载荷，没有重跑旧模型。
- 1927材料点+1561固定锚点、157637参考键，本体82270/界面22040/固定支架53327。保留完整m/r0/V/taper及原始atom/bond/type身份；材料质量守恒，锚点不计血栓质量。
- 本体G940Pa/K100000/3Pa/s0.15/sc.18；界面s0.140625/sc.16875，EXPLORATORY_UNCALIBRATED。完整endpoint-star变分力，局部精确消元，typed restart；固定支架无损伤。
- 新隔离库SHA da705cdd2358eab259d9f99ab402e743f6a64c829959eb094c148f8492552ae7；冻结SHA 49a534fc47ebe6827cc647f21855f0a55bd523549ef98ad5db52022a1108b65c。生产库/旧库不改。正式前生成工厂sizeof失配已补编译并复验，失败证据保留。
- C1_GAUSSIAN_A200：453历史表面点，200Pa面积预算、系数780.232224778Pa、sigma0.28mm、1.1605ms脉冲，当前位置/原面积切向。未叠加580点含栓CFD，人工Pa不是超声功率。
- 事件分辨VV：10502接受步、1182拒绝、2307接受细分步；最大Δd0.0999920876，最大归一化R0.0814340371%。没有RSS或墙钟上限；原数值门槛/输出2GiB保留。
- 末态本体/界面失效9567/420；真实脱离180点/174分量，质量2.25e-08kg；最大3点，≥200µm分量0。1–2点碎屑不能称临床微栓子。
- 同步真实存活图与锚连接；游离后即时复用无栓Y背景COM拖曳，附着仅原数值阻尼，非双向FSI。1.0dx粒子断后接触与实际Y壁罚接触分开；交接只查真实FP64状态，不追加长时输运。
- 单rank×单线程、无GPU；主段墙钟3130.263s，峰值RSS0.995697GiB。special_bonds初始化顺序优化避免无用三跳表，图/力学不改。
- 原历史217/208与修正199/190只读探索对照；新旧损伤律、界面取值1ppm、dt及MPI不同，不作单一因果改善归因。
- 单案例AI schema已经输出；EXPLORATORY_UNCALIBRATED，CRITICAL_LOAD/GC_EFFECTIVE=UNAVAILABLE，宏观强度未认证，formal_training_data_eligible=false；真实微泡streaming未接入。
- 先读verification/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/FINAL_STATUS.json、28项REVIEW、V3_ENERGY_LEDGER.json、V3_FRAGMENT_TOPOLOGY.json与V3_PILOT_FREEZE.json。代码analysis/v3_softening_pilot_v1/；仅在新授权目录重现，禁止重跑/覆盖本归档。
- 六组PNG/SVG：results/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/OPEN_RESULTS.html，真实点云/事件/曲线，没有连续医学外壳。WSL实际clot_rupture.md已更新、备份、diff、SHA并读回。
- 上轮3080点2.5ms已接受基线20261014T223001Z及全部历史保留。旧G0标称值不是新V3实际Gc。
- 本轮计算和交付完成后停止；下一项调查见28项报告第25项，不自动加大载荷、改材料、运行CFD/微泡/正式AI或追加正式case。
