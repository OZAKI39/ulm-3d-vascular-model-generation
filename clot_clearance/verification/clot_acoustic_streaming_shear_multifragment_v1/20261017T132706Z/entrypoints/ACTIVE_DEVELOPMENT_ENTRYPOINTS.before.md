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
