# LPS_MODE_I_LOW_LOAD_ADMISSION_AND_SOLVER_RECOVERY_V1

Run `20261013T122000Z`。

**尚未建立可信的小变形静态弹性基线。** 2 Pa的实际载荷、原始试样、原生力/独立梯度及小伸长范围通过；主路径240秒和一次600步CG恢复后，绝对/相对残余力及冻结力矩门槛仍未通过。没有放宽门槛、重试另一优化器或启动断裂。

| 核心量 | 实际值 | 裁定 |
|---|---:|---|
| 最大伸长 | 0.01537935146 | ≤0.03，通过范围检查 |
| 原生最大残余力 | 1.187074951e-12 N | >1e-13，失败 |
| 原生相对残余 | 2.390067922e-5 | >1e-5，失败 |
| 原生/独立梯度差 | 7.239375529e-12 | ≤1e-9，通过算子核对 |
| 正式动力断裂 | 0 | 未运行 |
| 可信柔度 / 实际Gc | UNAVAILABLE | 不填0 |

裂尖与加载端的保存候选伸长比约1.03421，但该状态未验收，**不能据此确认可信低载态中裂尖仍领先，也不能预测真实首断位置**。

## 1. 是否使用原始16000点？

是。逐ID坐标及参考NPZ SHA ef145b6d6db6c17584dce6cbea1fa1b9d95847d71a69265a73df6334d1695c03一致；没有重建近似试样。

## 2. 原参考键、预切和m是否一致？

909392完整参考键，初始预切16416、存活892976；131196条参考跨Y缝键、38992存活跨z=0键。先完整初始化再预切，原生m逐位保持，相对独立值最大差2.1921e-15。

## 3. 2 Pa节点力与面积正确吗？

上下各800真实原ID，每点面积dx²=2.5e-9m²，力±5e-9N；每面2e-6m²、合力±4e-6N。原生已施加，不只是写在独立势能里。

## 4. 合力、力矩自平衡吗？

参考外载合力/力矩为零。最终总残余力归一化1.80478e-13通过；当前残余力矩=[1.9177535194461277e-16, -4.712251804915775e-18, 9.835794135828044e-18] N·m，归一化9.60426e-09高于1e-10。固定周期胞X/Z取矩有边界力矩限制；相容Y分量约−4.712e-18N·m也超过冻结2e-18N·m尺度，不能忽略。

## 5. 20 Pa提前停止直接原因？

旧L-BFGS-B满足ftol目标下降条件，而非绝对力收敛。旧1132函数评估日志不是接受步序列，线搜索及条件数不能补造。见专门只读诊断。

## 6. 平移零模或其他病态方向？

3个平移HVP精确为零。固定Y胞仅Y转动保持内部长度，但dead-load势不保持，未投影。保存2Pa未接受态的随机方向与Y转动切向Rayleigh商约0.1361855与1.11644e-7N/m，差约1.22e6倍；支持弱方向尺度差异，不是精确谱条件数或已定位全部病态原因。

## 7. 优化数值尺度？

q=u/(1e-6m)，Pi/(1e-11J)，梯度=-P(F)×1e5；ftol=0、gtol=1e-11，对应单分量1e-16N。缩放不自动改变Hessian条件数。

## 8. 是否规约，有无改变约束？

仅均匀质心平移规约，无固定点、支撑弹簧或转动投影。平移探针保留伸长、储能与开度；最终使用未投影原生节点力，规约反力为零。

## 9. 求解方法？

一次L-BFGS-B，从真实零载A0开始；Pi=完整固定m独立LPS U+现有保守接触U−ΣF·u。接触每次评估，未用不完整原生bond PE。

## 10. 需要恢复吗？

是。预留的一次解析HVP Newton-CG同状态恢复，共601 HVP调用；600次CG未返回收敛，因此没有线搜索接受修正，立即停止。

## 11. 新增几次静态求解？

1个物理2Pa工况，1次主路径、1次同状态恢复调用；1个原生静态实例含完整参考run0和3次无积分器坐标/力评价。正式动力断裂0。

## 12. 优化器为何停止？

本轮主路径在240s预算检查停止，没有SciPy成功返回；恢复因CG_STOP_600停止。原驱动错误地排除所有异常后的恢复，已只修正调度启用预先保留的120s额度；没有增加总预算或第二恢复。旧停止记录与所有坐标保留。

## 13. 最大自由力？

原生1.18707495053e-12N=1.187075pN，是1e-13N门槛的约11.87倍，FAIL。

## 14. 相对残余通过吗？

2.39006792224e-05，高于1e-5，FAIL。

## 15. 能量和位移平台？

末3个不同接受态：相对能量跨度9.91963698425e-09≤1e-8，位置跨度1.60125913168e-12m≤1e-10m。窗口3在运行前已冻结，上轮20Pa无显式平台门槛；旧wide用10。额外报告末10能量跨度5.52117946271e-08未通过同一阈值；末10位置4.35009039195e-12m通过。不能用短窗口平台掩盖力学失败。

## 16. 最大伸长低于0.03吗？

是，0.0153793514577。本轮失败来自数值收敛，非已证实载荷超小应变范围。

## 17. 原生力与梯度一致？

相对差7.2393755293e-12≤1e-9；最大分量差6.79871611865e-19N。HVP与1nm方向差分相对差1.25453519189e-08≤1e-6。通过表示算子一致，不表示残差收敛。

## 18. 是否保留原始键与裂纹？

是，参考全键保留，固定拓扑诊断只屏蔽原16416预切键；活键892976准确，未从删键邻域重算m。

## 19. 新增断键？

0。实际所有活键类型1，s00=0.15不变；break no用于静态诊断，最大s同时低于0.03。无自然损伤事件。

## 20. 裂尖最大伸长？

仅未接受候选：0.0153793514577，s/s00=0.10252901；原键258238，atom [4490, 4650]，参考中点[-0.0010500000000000002, -0.000175, -2.500000000000001e-05]m。

## 21. 加载端最大伸长？

仅未接受候选：0.0148705690865，s/s00=0.0991371272；原键274983，atom [4779, 4780]，参考中点[-0.001025, 0.00012500000000000003, 0.00045000000000000004]m。

## 22. 可信低载态裂尖仍领先吗？

未建立可信低载平衡，不能回答为已认证领先。保存候选比值1.03421405、约高3.42%，仅空间诊断；不预测真实动态首断。

## 23. 最大位移与平均开度？

仅未接受候选：最大63.681003434µm，atom 100；平均端面开度18.546818396µm；平均位移模长9.759454234µm。全部16000点u已保存，不认证为弹性基线。

## 24. 单点割线柔度？

UNAVAILABLE/null。由于受载态未验收，本轮不报告可信割线柔度，也不画已测连续载荷—开度曲线。

## 25. 独立储能增量？

零载U0=0；候选U及增量3.66934372689e-11J=36.693437269pJ，仅未接受静态坐标的可回收能。

## 26. 终点外力位移收缩？

候选ΣF·u=7.41872735821e-11J=74.187273582pJ；Pi=-3.74938363132e-11J。不是路径外功。

## 27. 二者合理对应吗？

不作已验收线性弹性一致性结论。状态未平衡，不能把近似半F·u关系当作能量/力验收或材料实验验证。

## 28. 剪切传播尺度？

cs=sqrt940/1000=0.969535971m/s；H穿越1.031421ms，W穿越5.157106ms。

## 29. 24ms已准静态吗？

没有证据。24ms约4.654个W穿越，裂臂慢弯曲及裂纹传播尺度未认证；没有对未接受位移做已认证模态或静态柔度推断，也没有运行24ms动力学。

## 30. 为何不能报告Gc？

无自然扩展、新生面积和物理加载路径外功，且静态准入失败。G0_THEORY=0.07170070833333、G0_DISCRETE=0.06969623734234J/m²只保留为历史标称参考。实际Gc、分离协议能和带符号拓扑释放均null。

## 31. 需要改原生LPS或材料吗？

本轮没有证据要求修改。生产库、G/K/s00/alpha、horizon、点距、接触及显式theta语义均不变。失败不能直接归咎于材料太脆。

## 32. 资源和墙钟？

现场CPU cgroup={'/sys/fs/cgroup/cpu/cpu.cfs_quota_us': '1536000\n', '/sys/fs/cgroup/cpu/cpu.cfs_period_us': '100000\n', '/sys/fs/cgroup/memory/memory.limit_in_bytes': '45176848384\n'}；单进程/单线程，GPU未用。主路径240.161748s、峰值232.316MiB；恢复10.629846s；原生5.399946s、峰值654.434MiB。主路径记录1036次储能/力评估、988接受迭代；恢复1次额外状态评估和601HVP。末次显示评估计数从接受点续接少计2个旧试探，原日志保留，按真实1037行修正资源统计。全部989接受状态坐标/残余保留，最终规模以manifest为准。

## 33. 哪些科研图？

5组PNG/SVG：实际点阵/存活键/2Pa载荷，主求解与一次失败恢复，候选×1变形，分区伸长及未认证柔度，分项科学准入。所有数据来自实际状态，标STATIC_STATE_NOT_ACCEPTED；无医学外壳、虚构断裂或动画。

## 34. 保护历史及生产吗？

PASS，逐文件复核10632个保护资产；变化0、缺失0。未reset/clean/stash/commit/push/替换生产库或env.sh。

## 35. WSL文档实际更新？

PASS_ACTUALLY_UPDATED_AND_READ_BACK。真实/home/lzy/projects/clot_rupture.md仅追加⑤负面证据和结果链接，①–⑨及原有历史文本全部保留。

## 36. 备份diff SHA在哪？

本run的wsl_document/保存before、after、逐行diff和审计。before 151dd328a4a9ec629f7667c2d73ccf164a65ab96c6730b8b2f4c19929437df52；after/readback c8653d48a343bb95c75319369fe9348530485415600978cfc89846ef1755ddc8。WSL_DOCUMENT_UPDATE.json为入口。

## 37. 下一轮唯一优先项？

对同一保存2 Pa状态预注册一次弱方向/预条件静态收敛专项，保持物理载荷和真实自由度，先通过原生残差验收；不再追加更大载荷或夹具试错。 本轮已停止。

## 逐项科学状态

| 项目 | 状态 | 范围 |
|---|---|---|
| SPECIMEN_IDENTITY | PASS | Exact original 16000 nodes, 909392 bonds, 16416 cuts, fixed full m and periodic images |
| PRISTINE_REFERENCE_M | PASS | Exact original 16000 nodes, 909392 bonds, 16416 cuts, fixed full m and periodic images |
| PRECRACK_IDENTITY | PASS | Exact original 16000 nodes, 909392 bonds, 16416 cuts, fixed full m and periodic images |
| PERIODIC_REFERENCE_TOPOLOGY | PASS | Exact original 16000 nodes, 909392 bonds, 16416 cuts, fixed full m and periodic images |
| LOW_LOAD_PHYSICAL_INPUT | PASS | Actual 2 Pa native force added once after one V conversion |
| LOAD_FORCE_MOMENT_CONSISTENCY | REFERENCE_LOAD_PASS_CURRENT_EQUILIBRIUM_FAIL | Reference total load/moment zero. Current residual torque fails; full-vector periodic scope and Y component reported. |
| SOLVER_STOP_CAUSE | DIAGNOSED | Old20Pa ftol; new2Pa main wall budget, reserved recovery CG600 cap |
| RIGID_MODE_HANDLING | PASS_SCOPE_JUSTIFIED | Only 3 translations removed. Weak Y rotation is not removed. Final unprojected forces used. |
| SCALING_AND_GRADIENT_VALIDITY | PASS_NUMERICAL_IDENTITY_ONLY | Consistent q/Pi scaling; native gradient and read-only HVP FD pass, not convergence |
| ABSOLUTE_FORCE_CONVERGENCE | FAIL | Pre-frozen gates; plateau refers to last3 accepted states. Last10 retained separately. |
| RELATIVE_FORCE_CONVERGENCE | FAIL | Pre-frozen gates; plateau refers to last3 accepted states. Last10 retained separately. |
| ENERGY_PLATEAU | PASS | Pre-frozen gates; plateau refers to last3 accepted states. Last10 retained separately. |
| DISPLACEMENT_PLATEAU | PASS | Pre-frozen gates; plateau refers to last3 accepted states. Last10 retained separately. |
| SMALL_STRAIN_ADMISSION | PASS | Pre-frozen gates; plateau refers to last3 accepted states. Last10 retained separately. |
| NATIVE_FORCE_ENERGY_CONSISTENCY | PASS_FORCE_GRADIENT | Native force matches full independent gradient; native bond PE is incomplete and not used as U |
| NOTCH_LOADING_END_ELASTIC_RESPONSE | NOT_CERTIFIED_UNACCEPTED_STATE | Candidate ratio only; does not establish credible elastic ranking or first failure |
| LOW_LOAD_SECANT_COMPLIANCE | UNAVAILABLE | No accepted loaded endpoint |
| LOADING_TIMESCALE_SCOPE | WAVE_SCALES_ONLY_QUASISTATIC_NOT_CERTIFIED | H/W crossing 1.031/5.157ms; no accepted static-shape slow mode, no24ms dynamics |
| NEW_BOND_BREAKS_ZERO | PASS_STATIC_FIXED_TOPOLOGY | No new cuts; fixed diagnostic break no, actual s00 remains0.15 |
| STATIC_EQUILIBRIUM_ACCEPTED | FAIL | Absolute/relative force and torque not admitted; optimizer success not substituted |
| GC_EFFECTIVE_NOT_MEASURED | NOT_MEASURED | No natural fracture, new area or loading-work path; null rather than0 |
| PRODUCTION_LIBRARY_UNCHANGED | PASS | 10632 protected files checked, no historical replacement |
| HISTORICAL_RESULTS_PROTECTED | PASS | 10632 protected files checked, no historical replacement |
| WSL_DOCUMENT_UPDATE | PASS_ACTUALLY_UPDATED_AND_READ_BACK | Actual WSL section⑤ negative evidence update, all old text retained, SHA readback |

## 复现与限制

结果入口[OPEN_RESULTS.html](OPEN_RESULTS.html)，完整轨迹和残差见STATIC_SOLVER_TRACE.csv、accepted_optimizer_states/；没有物理时间轴。REPRODUCE.txt区分只读复核与需要下一轮授权的新求解。

[SciPy官方L-BFGS-B停止条件](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html)支持旧ftol/gtol解释，实际算法和环境身份以本run源码与PREFLIGHT为准。未重新进行大规模论文检索，也未把公开文档当作本地原生实现的替代依据。
