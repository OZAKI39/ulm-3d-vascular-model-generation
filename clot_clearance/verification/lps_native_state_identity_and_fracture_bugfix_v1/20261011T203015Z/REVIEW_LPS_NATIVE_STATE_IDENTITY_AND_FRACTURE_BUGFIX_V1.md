# LPS_NATIVE_STATE_IDENTITY_AND_FRACTURE_BUGFIX_V1 — 中文科研审查

本轮确认并在**独立诊断构建**中修复两项实际原生缺陷：参考m身份/持久化，以及跨MPI断键访问ghost的本地键数组。最终小试样通过。**唯一V3对照以中间m修复版本运行，发生原生段错误；最终双修复版未追加V3。不能宣称完整生产验证通过，也没有大碎块改善证据。**

运行标识：20261011T203015Z。WSL工作树：`/home/lzy/projects/temp_storage/clot_github_sync_20261011/clot_clearance`；服务器项目平铺目录：`/workspace/y_clot/project`。基于GitHub32e95bf快照；原生真实迁移提交519a5a，bond源码逐字节吻合锁定上游8de817dd，pair既有补丁一致。原生产库SHA仍为d1b2f53d6b942d1ca9277a8b8cd9c195ce75ea6f616681ec21c311e4faf3b41a。

## 25项明确回答

1. **wvolume是否错配？** 是。341点、12047初始键的同一三维非均匀邻域试样，参考m独立按初始ID/r0/taper/V求和。原生T0误差约1e−15；真实排序后153节点m错误，两rank迁移后85节点m错误，远超预注册2e−12门槛。
2. **历史V3默认排序是否发生？** 输入只设map array，没有关闭默认1000步排序。Verlet setup会调用Atom::sort；历史第1045步以后仍有明确ownership变化，结合Verlet源码支持后续排序调用。历史未保存本地顺序，**无法直接证明某一次具体非恒等排列或精确首次排序步**；本轮小试样保存了实际非恒等排列。
3. **历史是否有MPI迁移？** 是。壁回调按真实owned ID计数，rank0/rank1从1120/807变为1105/822，首次计数变化在第455步。小试样明确记录49个原子换rank及ID清单。
4. **m排序/迁移/restart是否正确？** 原版失败，断后扩容/重新run/restart亦有6节点m改变。最终版本真实排序、49点迁移、ghost、扩容、重新run及1/2rank新格式restart均保持初始m，最大相对误差约1.10e−15。旧损伤restart缺m信息不能自动补全。
5. **修了哪些C++？** bond_bpm_peri.cpp/.h：专用property/atom参考m和ready，跟随copy/grow/exchange/border/restart，刷新别名，不再销毁初始m。bond_bpm.cpp：process_broken两处端点访问限制为0≤index<nlocal。只在新source_owner_fix/build_owner_fix中实施。补丁没有材料、阈值、contact或theta公式修改。
6. **m/theta/F是否符合独立参考？** 是，最终无当步断键的静态、排序、迁移与断后固定拓扑状态均约1e−15。复用既有energy.cpp独立评估，所有方向差分最大相对误差5.69e-09，低于2e−6。
7. **断前弹性是否被改变？** 缺陷尚未触发时，三个配对T0的坐标、m、theta、力、存活键全部逐位一致。修复应改变身份已错乱后的结果，这不属于材料参数变更。
8. **断键同一步theta不一致？** 有。5条实际断键，最大Δtheta=0.06208、ΔF=4.41074 µN，约31.8526%独立节点力尺度。同坐标下一次力计算消除该差异。
9. **是明确bug还是时间层限制？** 测得力与源码明确的“断前theta＋断后存活键”预测一致至2.21e−16，判为EXPLICIT_SCHEME_LIMITATION。差异不小，不能宣称时间步收敛或忽略冲量；也不能把明确显式顺序直接改成另一算法。
10. **是否修theta顺序？** 没有。保留原失效判断和force pass。拓扑势能下降约5.9133e−11 J单独记录；ΔFΔt最大2.2054e−12 N·s，仅为估算，不是无积分测试中的实测速度跳变。
11. **参数和判据是否不变？** 是。G=940 Pa、K=33333.333333 Pa、bulk s00=.15、alpha0、dx50 µm、delta150.5 µm不变；V3原interface阈值、200 Pa Gaussian预算、0.5 µs、1.0dx接触、壁参数和拖曳完全保留。
12. **哪些字段语义有误？** largest_free_fraction/FI_free为全部分量最大值/其补量，名称不能解释成游离。debris_fraction及resolved_fragment_count是全部≤2/≥3节点分量，可能包含附着点；保留其原定义。新增8个明确字段并标记ATTACHMENT_AWARE_V2。旧classifier函数AST不变，新增failure_mode_detached_v2独立保存。
13. **真实游离质量和数量变化？** 没变。末态217点、2.7125e−8 kg、208个游离分量；201单点、6双点、1个4点。最大真实游离质量占比0.00207577；旧0.887390指1710点附着主体。旧≥3计数2→真正游离≥3计数1。213游离点为未解析碎屑。≥200/300 µm始终为零。
14. **统计修复是否改拓扑？** 没有。234保存帧重建与原component_id/attached逐项吻合，旧4项统计逐帧数值一致。5个人工图测试＋200/300µm边界与非等质量公式检查通过。独立新分类与旧泛化分类在221帧标签不同，显式标记POSTPROCESS_CLASSIFICATION_SEMANTICS_CHANGED；不是重写历史AI标签。
15. **跑了V3吗？** 跑了唯一一次正式尝试，中间m修复版，原C1完整输入字节一致、2MPI；发生SIGSEGV退出139，未完成4660步。最终owned guard加入后不再追加V3。
16. **V3比较结果？** 17完整帧至第320步（0.160 ms）坐标/速度/拓扑完全相同；第320步均30断键、vmax0.665175 m/s、KE约1.21625e−9 J、脱离0点。43完整力审计时刻至327步也完全相同。首次bulk均241；新轨迹未到旧首次脱离593或interface621。崩溃确切步未落盘，不能伪造。没有末态尺寸/锚反力改善结论。
17. **能说旧碎屑主要由bug造成吗？** 不能。P0确有风险，但完整修复版的V3因本轮预算未运行；部分对照发生在历史首次可证实迁移455步前且相同。旧217点/208分量负面结果保留。
18. **仍有哪些物理/数值限制？** 材料Gc未识别，瞬断不是耗散标定；显式theta lag未作时间步收敛；原生PE×V偏能计量不符合独立储能，已有文档限制保留；人工Gaussian非泡声场反演；背景无栓CFD、单向COM阻力、Gaussian与拖曳可能重叠为独立载荷问题。没有修改新表面、黏弹/疲劳/孔弹或微泡机制。
19. **必要原生与独立验证通过吗？** 最终小试样通过，包含额外证据驱动的原版两MPI断键崩溃复现和最终4键跨rank断裂/restart。完整V3最终版本未验证，不能总PASS。粒子pair源码字节未改；壁罚力离线公式检查通过。小试样pair zero不直接验证新断键接触，部分V3未有接触或拖曳事件，应标NOT_EXERCISED。
20. **计算数量/资源？** 12个小诊断场景执行（10成功，首次排序harness配置错误和原版跨MPI断裂2次失败），包含4次必要restart，共16次小型LAMMPS实例初始化；另1次V3，0次旧基线重跑。10成功小场景墙钟合计8.818s；两次独立编译合计78.880s，各8核，最大单子进程RSS约541144 KiB。V3 Supervisor记录墙钟7s，崩溃导致CPU/峰值RSS未留存，标NOT_RECORDED，非0。CPU cgroup15.36核，内存cgroup45176848384字节；无需GPU。输出分项规模见IMPLEMENTATION_CHECKS，原生两套src/build留服务器，约260 MB，失败V3约40 MB。
21. **原始数据保护？** 2283个既有结果、生产源码与库/可执行文件SHA复查0变化。保留失败日志、原版诊断及中间构建。只按授权改后处理代码，旧分类函数数值保持。未push、未删历史数据、未停其他科研任务。
22. **科研图？** 四组PNG/SVG：实际ID m审查、断键时间层DIAGNOSTIC_ONLY、1927点历史真实附着/游离图、唯一V3部分失败对照。第四组明确失败边界，不画假定末态或改善。
23. **WSL文档是否真实更新？** 是，实际`/home/lzy/projects/clot_rupture.md`。完整读取、before备份、写前SHA并发校验、after备份、diff和读回均完成；保留①–⑨、旧表格、参数和Gc未识别结论。审计见DOCUMENT_UPDATE_AUDIT.json。
24. **证据在哪里？** 本目录是完整审查入口；BUGFIX_SOURCE_DIFF.patch为最终3文件补丁，BUGFIX_WVOLUME_ONLY.patch为中间修复；provenance/native_source、python_before及clot_rupture.before.md为备份；所有小状态reference/T0/事件/断后/增长/restart均全精度NPZ。SOURCE_AND_BINARY_IDENTITY、FINAL_PATCH_IDENTITY、复现说明和SHA清单随附。最终服务器构建在本run的build_owner_fix，生产入口env.sh未切换。
25. **下一轮唯一优先项？** 下一次授权后，用最终m＋owned guard诊断库，对原C1做**一次**完整V3回归，逐步落盘ID/排序/迁移/断键/力证据，首先确认跨rank稳定与首次有意义差异；不调材料追求大块。本轮停止。

## 分项状态

| 项目 | 状态 |
|---|---|
| `WVOLUME_IDENTITY` | `CONFIRMED_BUG_FIXED_MINIMAL_TEST_PASS` |
| `ATOM_SORT_SAFETY` | `PASS_FINAL_REAL_NONIDENTITY_SORT` |
| `MPI_MIGRATION_SAFETY` | `PASS_FINAL_49_ACTUAL_MIGRATIONS` |
| `RESTART_SAFETY` | `PASS_FINAL_NEW_FORMAT_1_AND_2_RANK` |
| `LPS_FORCE_ENERGY_CONSISTENCY` | `PASS_INDEPENDENT_GRADIENT_NATIVE_PE_TALLY_LIMITATION` |
| `FRACTURE_TIME_LAYER` | `EXPLICIT_SCHEME_LIMITATION` |
| `POSTPROCESS_FRAGMENT_SEMANTICS` | `PASS_ATTACHMENT_AWARE_V2_LEGACY_RETAINED` |
| `NATIVE_BUGFIX_VALIDATION` | `PASS_SMALL_TESTS_FULL_V3_NOT_VALIDATED` |
| `V3_REPLAY_STATUS` | `FAILED_PARTIAL_SINGLE_ATTEMPT_LAST_SAVED_STEP_320` |
| `MATERIAL_GC` | `NOT_IDENTIFIED` |
| `PRODUCTION_LIBRARY_UNCHANGED` | `PASS` |
| `HISTORICAL_DATA_PROTECTED` | `PASS_2283_FILES_SHA256` |
| `MPI_FRACTURE_OWNERSHIP` | `CONFIRMED_ORIGINAL_CRASH_FIXED_MINIMAL_TEST_PASS` |
| `DOCUMENT_UPDATE` | `PASS_ACTUAL_WSL_UPDATED` |

## 误差、使用和可移植性边界

m单位m⁴（本核omega=1/r），原生力乘V成为N。最大绝对误差、显式参考尺度、出错ID及邻居数、所有rank顺序/归属、ghost值在WVOLUME_IDENTITY_TEST/MPI_SORT_MIGRATION_AUDIT及NPZ。初始零力不使用虚构分母。参考几何/X/r0从独立初始模型计算，绝不以错误原生数组自身作真值。

最终修复只在当前等体积LPS、固定原子集合、锁定MPI库和源版本验证；未声称覆盖所有BPM派生模型、动态增删原子和跨MPI数restart。旧损伤restart丢失参考m的迁移须另有初始键资料。原生PE的单独差异早于本轮，已有clot_rupture说明；本轮force-energy PASS指独立变分势能，不把原生PE表述为PASS。

本轮两次额外定位是由真实错误触发的有界回归，未扫描材料或dt。所有停止/失败保留。原始生产代码和历史轨迹仍是原版本；最终诊断库尚未部署。

[实际源码路径审查](WVOLUME_SOURCE_AUDIT.md) · [各项状态](FINAL_STATUS.json) · [原生力与独立势能](NATIVE_FORCE_ENERGY_AUDIT.json) · [单次V3失败对照](V3_REPLAY_COMPARISON.json) · [统计语义](FRAGMENT_STATISTICS_SEMANTIC_AUDIT.json) · [复现](REPRODUCE.txt)
