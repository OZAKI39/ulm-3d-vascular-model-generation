# LPS_VARIATIONAL_SOFTENING_NATIVE_ISOLATED_V1 科研验收报告

Run ID：`20261013T201800Z`。服务器项目：`/workspace/y_clot/project`。现场分支 `migration-20261009`，HEAD `cc5ee77149f55e58d715bb5a4660c6467231a826`。

**核心结论：已经在真正的 LAMMPS 原生 C++ 中复现上一轮冻结模型的节点力、完整储能、损伤历史及不可逆耗能。该结论限于 341 节点 / 12,047 原始键的规定运动学研究与本轮有限身份/接触诊断，不能推广为 V3、宏观强度或实验材料 Gc 认证。**

## 主要数值证据

| 同一规定加载协议 | sc=0.18 | sc=0.24 |
|---|---:|---:|
| 实际 native 接受状态数 | 354 | 354 |
| 最终不可逆 D / nJ | 173.156623298536 | 230.776817064078 |
| 与冻结离线 D 的绝对差 / J | 6.511565782e-21 | 8.099752558e-21 |
| 全路径节点力最大分量差 / N | 5.421010862e-20 | 5.421010862e-20 |
| 全路径完整 U 最大差 / J | 4.116050728e-21 | 7.358598729e-21 |
| κ/d 回退 | 无 | 无 |
| 最小有符号耗散增量 / J | 0 | 0 |
| 功账本最大相对残差 | 3.974518887e-14 | 1.021964594e-13 |
| 最終有效原生键 / 首次失效事件 | 0 / 12047 | 0 / 12047 |

关键独立证据：固定损伤梯度相对差 1.756198696e-09；完全失效跨块响应严格 0；真实排序 341 槽；真实 2MPI 迁移 49 atom；restart 后 262 帧继续一致；闭合产生 49 对原生接触。

## 冻结模型与实施边界

完整推导及变量对应见 [C++ 审查](VARIATIONAL_CPP_IMPLEMENTATION_AUDIT.md)、[模型冻结](NATIVE_SOFTENING_DESIGN_FREEZE.json)、[逐文件身份](SOURCE_AND_LIBRARY_IDENTITY.json)。新 style 是对最终双修复源码的两文件独立扩展。使用低秩精确局部消元、独立端点内部状态、共享不可逆 κ/d、固定参考 m 和积分权重。损伤在显式接受 token、固定真实坐标处一次性提交；trial/run0/差分不改变历史。完整 U 与力出自相同变分能量。

所有新力由真实 LAMMPS atom force 按已固定体积换算 N，完整能量导出 J；旧 BPM force-density 和 bond tally 约定明确区分。全节点规定坐标协议不是自由动力学或准静态平衡，无物理 dt、动能和阻尼项。没有使用优化器求值次数定义物理时间。

历史五键断裂：新模型对自身冻结参考正确，但相对旧删键 LPS 的节点力向量 L2 差为 **8.814180703270e-06 N**。旧保存原生力与独立旧二值参考差 1.312812572e-20 N。这是原始冻结新自由能对失效端点内部约束松弛的本构差异，不是 sc 调整或移植误差；未强行宣称二者相等。

状态按稳定原始键身份存储；每 rank 复制整个小图、汇集 owned 坐标并仅向 owned atom 施力。因此通过 2MPI 的实际迁移，不代表大规模 MPI 伸缩性。未实际触发原子 nmax 扩容，该项保留 NOT_EXERCISED；多体积、多类型、周期/重复端点、任意复杂 topology 及动态积分未认证。restart 使用独立 LPSVAR_RESTART_V1，错误版本实测拒绝。

## 40 项明确回答

### 1. 是否完全复用了上一轮正式离线软化模型？

是。先核验上一轮 127 项封存资产，manifest SHA `3265476d2e6b464f0e2b771fc52731310f97bdc746e54287654725a2abff6985`；正式 model.py SHA `3eb1396b376cf30e40586c26ba649909ab9382d9eac9607666359258273acf74`。完整自由能、独立端点状态、双线性 κ→d、固定 m/V/r/taper 逐项继承。公式来自完整冻结源码，不来自摘要。

### 2. 是否避免恢复旧版两端共享内部状态的错误？

是。每一无向原始键共享 κ/d，但两个端点有各自的邻域内部变量 v。完全断开 A 块移动后，B 块非接触力变化严格为 0；未恢复上一轮初稿的共享 v。

### 3. 新原生 C++ 修改了哪些文件？

仅新增 `src/BPM/bond_bpm_peri_variational.h` 和 `.cpp`，新样式 `bpm/peri/variational`；原 4,729 个 src 文件 SHA 相同。真实新增源码在 `source_snapshot/`，完整新增补丁 `NATIVE_SOURCE_PATCH.diff`。Python 驱动在 `analysis/native_softening_v1/`，不属于生产路径。

### 4. 实际编译出的诊断库 SHA 是多少？

最终隔离库 SHA256：`b80bfc0a99f4160d62567c0e2d9ddc611004ad40f79c8e49370293c681560773`，11,745,088 字节。服务器 `build/liblammps.so`，交付副本 `isolated_library/liblammps.so`。修正前库另存 implementation_history，不能混作最终版本。

### 5. 旧生产库是否保持不变？

是。结束重新计算 SHA：`d1b2f53d6b942d1ca9277a8b8cd9c195ce75ea6f616681ec21c311e4faf3b41a`。旧双修复诊断库仍为 `00b246ebe4224e648c27b74afb64fc22b9cea4ba6ffb23f7c40d96c37ca22ae9`；env.sh 与 Python 绑定未改。

### 6. 未损伤原生力与旧 LPS 是否一致？

是。实际运行新库中原样保留的最终双修复 bpm/peri 与新样式，对照逐 atom ID 三维力；旧样式源文件逐字节与基线一致。参考态力 0；小仿射相对 L2 差 6.375870502e-16，最大分量差 5.293955920e-22 N；非仿射相对差 6.256769681e-16。

### 7. 未损伤完整储能是否一致？

是，与独立冻结完整能量比较。小仿射绝对差 1.615587134e-27 J，非仿射 4.796274304e-28 J；参考态 U=0。本轮新 style 的完整 native tally 乘 V 后一致。旧 bpm/peri 不完整 PE 没有被用作完整储能参考。

### 8. 当前局部内部状态采用何种求解方法？

精确对角加秩一消元，等价于每个端点正定二次能量极小化；不是迭代优化、伪逆或经验缩放。6 个端点另以稠密 FP64 系统独立核对，最大 v 差 1.0587911841e−21 m。

### 9. 是否利用了低秩结构提高计算效率？

是。每个端点以标量 p 消去秩一耦合，每次求力仅对有向边做线性遍历。341 个 endpoint solves / 次 / rank；每条正式路径 2,131 次能量力求值、726,671 个端点求解。每次内部求解约 0.43 ms；这不包括全部 Python、通信和输出开销。

### 10. 部分损伤的原生力与离线参考最大差是多少？

同一实际坐标和历史：sc=.18 固定状态最大分量差 2.371692252e-20 N、完整 U 差 9.264422861e-22 J；sc=.24 分别 5.421010862e-20 N、1.323488980e-21 J。两个状态 κ/d 逐元素一致；另核对非均匀损伤。相对 L2 差分别 3.519702540e-16、5.073106526e-16。

### 11. 固定损伤下力是否等于自由能负梯度？

在已测试自由度内通过。8 个自由度、2 个步长、正负中心扰动共 32 次 native 能量计算；最大绝对差 8.673891077e-16 N，相对实际全局分量尺度差 1.756198696e-09，真实分母 4.939014645e-07 N。未设置人为分母下限，扰动期间提交次数 0。

### 12. 新模型是否保留固定参考 m 和原始积分权重？

是。341 个参考 m、V=dx³、原始 r 与 taper 均来自冻结点阵，损伤不会重新归一化 m。排序、MPI、restart 均核查固定 m；全程未换点阵或 horizon。

### 13. 是否正确处理 d=0、0<d<1 和 d=1？

是。d=0 恢复完整弹性，部分损伤使用冻结串联柔度，d=1 释放端点约束且该键非接触力严格为零；不加 epsilon 残余刚度。注意新模型二值极限不等于旧直接删键法，另见下方单列差异。

### 14. κ 与 d 是否具有稳定原始键身份？

是。从 body.data 恢复真实 12,047 个原始 bond ID，以 bond ID + 两端 immutable atom tags + type + 零镜像标识；不使用当前 bond loop 下标。即使物理 bond slot 被删除，完整历史仍在版本化原始图中。重复端点对、周期/三斜几何明确不支持，重复身份已实测拒绝。

### 15. 是否有任何卸载愈合？

没有。已测试加载、卸载、压缩、再加载全路径 κ/d 不下降；完全失效 restart 与后续闭合也没有愈合。只对本轮已执行状态作此结论。

### 16. 是否出现负损伤耗散？

没有。记录原始有符号 U_old−U_new，不取绝对值、不截负值。两条正式路径最小损伤增量为 0 J，所有已测增量非负。

### 17. 原生 sc=0.18 和 0.24 的耗能分别是多少？

sc=.18：1.7315662329853586e-07 J = 173.156623298536 nJ；sc=.24：2.3077681706407765e-07 J = 230.776817064078 nJ。为同一运动学协议的离散累计耗能，不是单位面积材料 Gc。

### 18. 是否复现已接受离线耗能和广义反力？

是。累计 D 绝对差分别 6.511565782e-21、8.099752558e-21 J。两个 sc 的 Q 峰值均 38.3625 mN，初始广义刚度 5115 N/m。旧瞬断 143.859377 nJ 仅复用已接受离线对照，未虚构原生部分损伤旧曲线。

### 19. 实际力学与功账本是否闭合？

是。沿冻结的分段仿射路径，以旧 d 两点 Gauss 外功积分，再在接受坐标提交损伤；W−ΔU−D 最大相对残差分别 3.974518887e-14、1.021964594e-13。无动能/阻尼项，没有把 run0 次数称为物理时间。连续解析积分与离散更新的历史差异没有被当作 C++ 错误。

### 20. 完全失效后断块之间是否仍存在虚假弹性耦合？

本轮未发现：2,470 个跨面键完全失效，A 块独立位移时对应非接触跨键力=0，B 块力变化=0。整体路径最终 12,047 个键全失效，原生 bonds=0、最终节点力和 U 均严格为 0。

### 21. 部分损伤键是否与接触重复施力？

未出现。原生 1–2 special exclusion 在 0<d<1 时保留。部分损伤闭合的附加 pair 力最大 2.117582368e-22 N，为 FP64 相减舍入，不能写成数学上的严格零。

### 22. 完全失效后的原生接触是否实际经过验证？

是。完全失效提交移除原生键与对应 1–2 屏蔽，正常邻居重建后闭合出现 49 对真实接触；保留原 1.0dx 罚力。逐物理对仅计一次的参考与 native 最大分量差 4.235164736e-22 N，相对差 1.826333786e-16；闭合 κ/d 不变。大开口只为触发诊断断面，不代表有效材料应变区间，未接入壁面/V3。

### 23. 是否真实触发原子排序？

是。标准 atom_modify sort + run0 触发 341 个 owned 数组槽非恒等改变，保存排序前后实际 ID 序列。未使用手动直调 Atom::sort。

### 24. 排序前后损伤历史是否一致？

是。原始 ID、κ、d、m、节点力、完整 U 的排序差均为 0。排序后新接受状态产生预期 217 个失效事件且只提交一次。

### 25. 是否真实触发 MPI ownership 迁移？

是。最终 2MPI 分区实际迁移 49 个原子，ghost 数从 [145,196] 变为 [194,147]。第一次分区位于同一空隙，没有迁移，明确保留 MPI_MIGRATION_NOT_EXERCISED，不计通过。最终分区 −40 µm 按原坐标层位置确定。

### 26. 跨 rank 损伤更新是否正确？

本轮 2MPI 用例通过。2,470 条跨 rank 原始键，更新前 49 条部分损伤，更新后 217 条跨 rank 失效；所有 rank 历史/账本一致，一个 accepted commit 对应 217 个唯一事件。实际单 rank 与两 rank 的力及 U 差为 0，迁移前后 κ/d/m/F/U 也均相同。每 rank 复制小图并只对 owned atom 施力，不能外推生产伸缩性。

### 27. restart 后是否保持完整 κ、d 及参考 m？

是。frame91 的 12,047 个部分损伤键写实际 LAMMPS restart，重新构造实例读回后 κ/d/m/IDs/F/U/D 均一致。继续后续 262 个接受状态，与不中断原生路径的力、U、d、D 差均为 0。完全失效 restart 也实测无愈合。

### 28. 是否重复累计历史耗能或失效事件？

没有。两个主路径最终各 354 次 accepted commit、12,047 次首次失效；每条路径 4 个重复 token 检查无新增 D/事件。restart 重放已保存 token 无新增，排序/MPI/接触重复提交也不增加。相同 token 配不同坐标返回 −5 并保持历史；旧 style 接口返回 −99。

### 29. 是否发现新的原生 Bug？

代码审查发现第一版一个未被科学路径调用的辅助 lpsvar_sort 入口不安全：直接 Atom::sort 绕过 exchange/borders 生命周期且缺少 sortfreq 条件。已用唯一一次最小 C++ 修正删除，数学核未改，最终库重做必要门槛与两路径。未发现力学核差错。首次 MPI 未覆盖迁移是诊断分区问题，修正并保留记录。拷贝符号链接和缺少 /usr/bin/time 的准备错误另存 SETUP_RECOVERY，不伪装成力学失败。

### 30. 本轮是否修改材料 G/K/s0/sc？

没有。G=940 Pa、K=100000/3 Pa、ρ=1000 kg/m³、dx=50 µm、horizon=150.5 µm、alpha=0、s0=.15、sc=.18/.24 固定。未扫描、拟合或扩大 sc。

### 31. 当前是否能保证宏观峰值强度不变？

不能。只证明同一仿射受控协议下广义 Q 峰值相同；38.3625 mN 不是宏观抗拉强度，5115 N/m 不是材料模量。

### 32. 是否能够报告实际材料 Gc？

不能。本轮 D 单位是 J，尚无自然裂纹传播、有效裂面及几何/尺度收敛认证。material_Gc_J_m2=null。

### 33. 能否说明 6–8 J/m² 目标已经实现？

不能。没有实现或认证实验断裂阻力，也未通过调参追求该值。历史 G0_THEORY=0.07170070833333 J/m²、G0_DISCRETE=0.06969623734234 J/m² 保留其历史语境，不能替代本轮材料 Gc。

### 34. 实际新增多少正式小型测试？

只有一个真实 341 节点主图、两个预注册 sc 的正式 354 状态路径。机械门槛组包括 3 未损伤状态、历史 5 键冻结状态、2 部分损伤状态、1 非均匀固定损伤状态、6 个端点稠密核对和 32 次差分能量求值；随后各一个断块、排序、2MPI、restart、接触场景，另 4 项接口拒绝检查。因唯一安全修正，机械门槛/两条路径/断块/排序各执行两批，完整主路径实际共 4 条；MPI 执行 2 次，其中仅 1 次真正迁移通过。restart 1 场景含 3 个原生实例。未把这些求值计成不同材料试验。

### 35. CPU、内存、墙钟和文件成本是多少？

首次配置+编译 36.969028 s（编译段 34.097954 s），最终增量 4.338934 s。编译限 12 jobs，采样峰值 9.388 核、进程树 RSS 1928.69 MiB。最终单 rank 门槛 4.029501 s、两路径 42.679517 s（CPU 42.734104 s、RSS 747.30 MiB）；2MPI 1.380293 s，每 rank 峰值 RSS [1027984, 1031372] KiB；restart 14.654854 s，接触 0.739370 s。GPU 未使用。逐段完整成本、首次复验、服务器源码/build 磁盘占用在 RUN_RESOURCES.json；精确交付字节数在 PAYLOAD_MANIFEST.json。编译 CPU/RSS 为每 0.5 s 采样，不声称捕获短瞬态；MPI 各 rank 峰值不能相加称同时峰值；开发、传输耗时未计入上述求解墙钟。

### 36. 生成哪些实际科研可视化？

5 组各 PNG/SVG：原生/离线逐 ID 力与 U；6 状态原生点阵损伤；Q/U/W/D/残差；原始跨面键、部分损伤、真实分离、实际闭合接触；排序/MPI/restart 历史、力/能量误差及成本。OPEN_RESULTS.html 静态浏览。所有新曲线源自保存的实际 C++ 状态，旧瞬断离线对照显式标注。5 张 PNG 已实际打开审图，调整标签重叠后复核。

### 37. 生产原生库及历史数据是否受到保护？

是。结束核验 22,013 个历史文件、5,399,515,016 字节，SHA/字节数零变化；另 30 个源/库/环境资产与 4,729 个原源码文件一致。无 reset/clean/stash/commit/push、无 rsync --delete、无生产部署，最终 Git 分支/HEAD 仍单列记录。

### 38. WSL clot_rupture.md 是否实际更新？

是，实际 `/home/lzy/projects/clot_rupture.md` 已读取、备份、并发 SHA 核验、插入科研结果并读回一致。保留所有旧文字、章节、参数及离线历史结论。

### 39. 备份、diff 和 SHA 在哪里？

本次结果目录 `wsl_document/`：clot_rupture.before.md、clot_rupture.after.md、clot_rupture.diff、DOCUMENT_UPDATE_AUDIT.json。before SHA `a2cdba711f63a3d775d1d22e836e3033b88e2a161d42f38377e8d4f42c9ebbe8`；after SHA `61e0fe9947692ffea38272645758f38a994d900a56af00a177759e419bed476c`。包含实际 WSL 路径和读回检查，不只是待应用补丁。

### 40. 是否具备进入下一轮非均匀裂纹验证的条件？

具备进入隔离的小型非均匀裂纹研究的条件，非生产部署准入。下一轮需另行研究非仿射裂纹演化、宏观强度、真实分离能与尺度/时间层收敛。本轮到此停止，未自动启动下一轮、V3、CFD、微泡或 AI。

## 分项验收状态

不以单一 PASS 替代边界。下表与 FINAL_STATUS.json 一致：

| 验收项 | 结论 |
|---|---|
| `OFFLINE_MODEL_IDENTITY` | `PASS_FROZEN_MANIFEST_AND_SOURCE_SHA` |
| `BASELINE_LPS_FORCE_PARITY` | `PASS_INTACT_PER_ATOM_ID` |
| `BASELINE_LPS_ENERGY_PARITY` | `PASS_COMPLETE_VARIATIONAL_U;LEGACY_PE_NOT_REFERENCE` |
| `ENDPOINT_STAR_STATE_SEPARATION` | `PASS_INDEPENDENT_DIRECTED_STARS` |
| `NATIVE_SOFTENING_CPP_COMPILED` | `PASS_ISOLATED_ONLY` |
| `FIXED_DAMAGE_FORCE_GRADIENT` | `PASS_NATIVE_FIXED_HISTORY` |
| `PARTIAL_DAMAGE_NATIVE_PARITY` | `PASS_BOTH_SC_AND_NONUNIFORM_DAMAGE` |
| `DAMAGE_HISTORY_IRREVERSIBILITY` | `PASS` |
| `NONNEGATIVE_DISSIPATION` | `PASS_RAW_SIGNED_LEDGER` |
| `LOADING_UNLOADING_RELOADING` | `PASS_KINEMATIC_PROTOCOL` |
| `FULL_DAMAGE_ZERO_TENSILE_TRANSFER` | `PASS_EXACT_ZERO` |
| `DISCONNECTED_COMPONENT_INDEPENDENCE` | `PASS_EXACT_ZERO_REMOTE_RESPONSE` |
| `STABLE_BOND_IDENTITY` | `PASS_ORIGINAL_IDS_AND_ENDPOINT_TAGS` |
| `ATOM_SORT_DAMAGE_SAFETY` | `PASS_341_CHANGED_SLOTS` |
| `MPI_DAMAGE_HISTORY_SAFETY` | `PASS_2_RANKS_49_MIGRATED_ATOMS` |
| `RESTART_DAMAGE_HISTORY_SAFETY` | `PASS_PARTIAL_CONTINUATION_AND_FULL_FAILURE` |
| `DAMAGE_COMMIT_TIME_LAYER` | `PASS_EXPLICIT_ACCEPTED_TOKEN;TRIAL_READ_ONLY` |
| `POSTBREAK_CONTACT_CONSISTENCY` | `PASS_1DX_TWO_COMPONENT_CLOSURE` |
| `SOURCE_AND_BINARY_REPRODUCIBILITY` | `PASS_ISOLATED_BUILD_RECIPE_AND_SHA` |
| `MACROSCOPIC_STRENGTH_NOT_CERTIFIED` | `True` |
| `MATERIAL_GC_NOT_CERTIFIED` | `True` |
| `EXPERIMENTAL_FRACTURE_NOT_VALIDATED` | `True` |
| `FULL_V3_NOT_RUN` | `True` |
| `PRODUCTION_LIBRARY_UNCHANGED` | `PASS_SHA256_UNCHANGED` |
| `HISTORICAL_RESULTS_PROTECTED` | `PASS_22013_FILE_SHA256` |
| `WSL_DOCUMENT_UPDATE` | `PASS_ACTUAL_WSL_READBACK_WITH_BACKUP_DIFF_SHA` |
| `OLD_BINARY_DELETION_EQUIVALENCE` | `NOT_EQUIVALENT_BY_FROZEN_MODEL_DESIGN` |
| `NATIVE_ARRAY_CAPACITY_GROWTH` | `NOT_EXERCISED;KEYED_GRAPH_INDEPENDENT_OF_ATOM_NMAX_BY_SOURCE_AUDIT` |
| `PERIODIC_OR_DUPLICATE_ENDPOINT_GEOMETRY` | `REJECTED_UNSUPPORTED` |
| `PRODUCTION_SCALABILITY` | `NOT_CERTIFIED_REPLICATED_SMALL_GRAPH` |
| `eligible_for_nonaffine_small_crack_research` | `True` |
| `production_deployable` | `False` |
| `material_Gc_J_m2` | `None` |
| `material_parameters_changed` | `False` |
| `full_V3_runs` | `0` |
| `native_main_nodes` | `341` |
| `FIGURE_VISUAL_REVIEW` | `PASS_5_PNG_REVIEWED_SVG_VALIDATED` |

## 复现与审计入口

- [五组图](OPEN_RESULTS.html)、[FP64 数值比较](NATIVE_SOFTENING_COMPARISON.csv)、[资源记录](RUN_RESOURCES.json)。
- [复现步骤](REPRODUCE.txt)、[最终库身份](ISOLATED_BINARY_IDENTITY.json)、[新增源码补丁](NATIVE_SOURCE_PATCH.diff)。
- [唯一 C++ 安全修正](MINIMAL_REPAIR_RECORD.json)、[补丁](MINIMAL_REPAIR.diff)、`implementation_history/pre_safety_cleanup/` 保留初版库、源码、真实计算与首次未触发迁移的结果。
- [历史保护结束检查](HISTORICAL_PROTECTION_AFTER.json)、[文档更新](wsl_document/DOCUMENT_UPDATE_AUDIT.json)、[图像审阅](FIGURE_VISUAL_REVIEW.json)。
- `PAYLOAD_MANIFEST.json` 封存可交付文件；`DELIVERY_VERIFICATION.json` 在 WSL 校验实际下载结果后生成。服务器保留完整隔离 source/build，本地交付实际新增源码、最终库副本和必要证据，不同步整个冗余构建树。

本轮不自动进入下一轮研究。宏观峰值强度、材料 Gc、6–8 J/m² 目标、实验断裂验证、完整 V3 和 AI 材料卡均未认证。
