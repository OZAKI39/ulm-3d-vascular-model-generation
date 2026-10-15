# 原生参考 m、断键与通信源码审查

本表行号对应 provenance/native_source 中保存的实际生产源码，不是在线最新版。服务器原生仓库的迁移快照提交为 519a5a049285af2bb6b99e95506e3df6a6ebbb3c；BPM bond 三个文件与记录的上游 8de817dd79bfe4525d5d39246a212d833e6dee07 字节一致。pair 文件保留既有 1.0dx 接触补丁，SHA 与 SOURCE.json 一致；见 SOURCE_DIFFERENCE_REPORT.json。没有用上游覆盖服务器。

| 路径 | 原始实现与具体位置 | 实测结论 / 最终处理 |
|---|---|---|
| 分配 | bond_bpm_peri.cpp:206 compute，215–225 destroy/create wvolume；头文件内部裸数组 | 不是 Atom 管理属性。移入独立 property/atom 持久状态，theta/winv 保持每步 scratch |
| 初始化 | bond_bpm_peri.cpp:505 compute_wvolume，524–528 对存活键 r0 V h 求和，531 标志置1 | 初次求和正确，T0 m/theta/F 约1e−15 |
| SORT_IDENTITY | atom.cpp:2366 Atom::sort，2447–2453 AtomVec::copy；atom_vec.cpp:348–350 仅调用已注册 fix copy_arrays | 裸 m 未注册，真实非恒等排序后153节点错配；最终 formal fix copy 随ID |
| MPI_MIGRATION_IDENTITY | comm_brick.cpp:599 exchange，657 pack、702 unpack；atom_vec.cpp:1159/1250 和1240/1332 fix回调 | 49节点真实换rank，85节点m错误；最终 property/atom exchange 保留m |
| GHOST_COMMUNICATION | bond_bpm_peri.cpp:263–270，1065/1080 forward pack/unpack | 只把当前owner数组值发给ghost，不能纠正owner身份错配；修复后真实ghost逐ID核对 |
| ARRAY_GROWTH | bond_bpm_peri.cpp:214–225 扩容销毁m并重置标志 | 断后5键试验扩容把6节点m改为存活邻域；最终 property/atom grow保持旧元素 |
| 重新run | bond_bpm_peri.cpp:825 init_style，862 wvolume_setup=0 | run pre yes重建存活邻域，违背固定初始m；最终初始化标记随原子保存，仅未初始化节点求和 |
| RESTART_IDENTITY | bond_bpm_peri.cpp:1009/1029 只保存类型材料系数；裸m没有逐原子restart | 原版断后重启丢失原始m；最终 fix pack/unpack保存m与ready，单rank和2rank均验证 |
| REFERENCE_M_AFTER_BREAK | compute先theta后break，allow_breaks在setup关闭；普通连续计算不重建m | 连续、未重排时保留m本来正确，不是缺陷；修复保持此定义，绝不按存活键覆盖初始m |

正式存储机制：fix_property_atom.cpp:217–219 注册 GROW/RESTART/BORDER；587 grow_arrays；636 copy_arrays；797/833 exchange；869/909 restart。独立保留 ID `BPM_PERI_REFERENCE_M`，字段 `d_bpm_peri_reference_m` 和 `d_bpm_peri_reference_ready`。额外ready区分“已初始化但m=0”与新存储。名称冲突明确报错，不复用或覆盖 d_vfrac/d_s0/d_smin。每次compute刷新m指针别名。原生生产lib与源码未改。

旧版受损restart没有保存原始m，不能仅凭当前剩余键恢复丢失的参考邻域；本轮通过的restart是最终修复版本生成/读取的检查点。未验证动态添加/删除原子、跨MPI数目restart、混合BPM/PMB/VES/EPS、多种horizon及非等体积模型，不能外推认证。

## 断键显式时间层

bond_bpm_peri.cpp:270 compute_dilatation → 316起force loop → 345–351 阈值判断与process_broken → 363起存活键力。首个受控事件5条键，原生力与“断前theta、断后存活键”的预测误差2.21e−16；同坐标下下一次计算恢复变分力。这个明确的显式次序解释了31.85%的局部力尺度差，绝不称其微小或物理已收敛。本轮没有改theta算法、失效阈值或积分器。独立势能拓扑下降与 ΔF Δt 冲量估计分开报告。

## 跨MPI断键内存缺陷（由唯一V3失败触发的限定补查）

BondBPM::process_broken 位于 bond_bpm.cpp:422。代码把ghost经atom->map映射后，仅用 `id != -1`（460）和 `jd != -1`（477）检查，即进入只对owned atoms有效的 num_bond/bond_atom/bond_type/history。映射可以仍是ghost（>=nlocal）。V3与原生产库同一341节点小试样的两MPI断键均在该函数段错误；反汇编+0x267是 setflag[bond_type[jd][m]] 访问。参见两个launcher原始日志和 process_broken.disassembly.txt。

最终仅将两处条件改为 `0 <= index && index < nlocal`。跨rank每端由其owner处理；低ID owner的 nbroken/event规则、UPDATE_SPECIAL_BONDS调度和owned bond history换尾删除保持原样。最终同一小试样4条实际断键通过，同坐标下一次力计算和两rankrestart均符合独立参考。没有把此内存修复称为theta时间层修复。

最终补丁含3个文件：bond_bpm_peri.cpp、bond_bpm_peri.h、bond_bpm.cpp。中间m修复库的V3在第320步仍与历史逐坐标/速度/拓扑一致，最后完整力审计到327步亦相同；随后崩溃。最终两处owned guard只经过小试样复验，未再次V3验证。

## 数值门槛与试验实现边界

DESIGN_FREEZE在改原生源码前保存。m误差 / max独立m ≤2e−12；theta ≤2e−10；节点力最大范数误差 / 最大独立节点力 ≤2e−9；独立储能梯度方向差分 ≤2e−6。零尺度只报绝对误差和null相对值，无隐藏分母下限。原始m的完整ID、原始键、r0、taper、V、degree在各reference.npz。力观察为原生力密度乘V；所有比较坐标相同。

原生ABI观察器不改生产C++，只读取m/theta/winv；扩容是明确调用AtomVec::grow的测试刺激。第一次排序harness运行错误地在pre no路径才打开排序，未初始化排序bins而崩溃；保留 failed_harness_sort_initialization，未用它证明材料缺陷。随后在初始setup预置sort2，用刚体90°旋转真正触发非恒等排序；该更正不放宽门槛。V3导致增加一次同几何跨rank断键定位及最终回归，未进行材料、网格、dt、载荷扫描。

参考机制：[官方BPM/PERI](https://docs.lammps.org/bond_bpm_peri.html)、[atom_modify](https://docs.lammps.org/atom_modify.html)、[property/atom](https://docs.lammps.org/fix_property_atom.html)、[restart](https://docs.lammps.org/read_restart.html)。最终数值结论以保存的服务器锁定源码和实际运行记录为准。
