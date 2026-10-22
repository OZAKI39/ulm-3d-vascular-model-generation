"""Build the review package from immutable measured artifacts, without simulation."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,html
ROOT=Path(__file__).resolve().parents[4];R=ROOT/'particle_3d/reports/particle9a4_population_inlet';D=R/'data'
def read(p):return json.loads(p.read_text())
def write(name,s):
 with (R/name).open('x') as f:f.write(s.strip()+'\n')
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def table(headers,rows):
 return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,row))+' |' for row in rows])
b=read(D/'baseline_git.json');s=read(D/'inlet100k/summary.json');q=read(D/'qacc_diagnostic.json');m=read(D/'smoke30_summary.json');c=read(ROOT/'particle_3d/contracts/P9A4_CONTINUOUS_INFUSION_V1.json');t=read(D/'final_regression.json');a=read(D/'code_logic_audit.json');g=read(D/'inlet_audit_gate_evidence.json');syn=read(D/'synthetic_poiseuille_audit.json')
assert m['passed'] and all(g['checks'].values()) and t['tests']==115 and t['failures']==0 and t['additional_tests']['tests']==50 and t['additional_tests']['errors']==0
assert 'PASS: 46670 published file hashes' in (R/'logs/final_protection_manifest.txt').read_text()
stage='P9A4_POPULATION_INLET_READY_FOR_LARGE_SAMPLE_REVIEW'
qtable=table(['D (µm)','Q_acc/Q_in','点态 95% Wilson CI'],[[f"{x['diameter_um']:g}",f"{x['fraction']:.5f}",f"[{x['ci95'][0]:.5f}, {x['ci95'][1]:.5f}]"] for x in q['rows']])
dtable=table(['分布；单位 µm','mean','p10','p50','p90'],[[name]+[f'{v[k]:.6f}' for k in ['mean','p10','p50','p90']] for name,v in [('source',s['source_diameter_um']),('entering',s['enter_diameter_um'])]])
summary=dict(baseline_git_commit=b['ACTUAL_HEAD'],expected_git_commit=b['EXPECTED_HEAD'],new_branch=b['branch'],new_worktree=str(ROOT),
 flow_sha256=c['flow_SHA'],Q_in_m3_s=c['Q_in_m3_s'],source_contract_sha=digest(ROOT/'particle_3d/contracts/P9A4_CONTINUOUS_INFUSION_V1.json'),
 source_distribution=c['source_distribution'],source_D_min_um=c['D_min_m']*1e6,source_D_max_um=4.,concentration_semantics=c['concentration_semantics'],
 C_total_m3=c['C_total_m3'],F_original_4um=c['F_original_4um'],C_modeled_m3=c['C_modeled_D_le_4um_m3'],
 concentration_is_ground_truth=False,concentration_semantics_resolved=True,absolute_clock_is_model_assumption=True,
 lambda_source_s_inv=s['lambda_source_s_inv'],N_source_proposals=s['N_source_proposals'],N_accepted=s['N_accepted'],N_rejected=s['N_rejected'],
 acceptance_fraction=s['acceptance_fraction'],rejection_fraction=s['rejection_fraction'],lambda_enter_measured=s['lambda_enter_measured'],lambda_enter_expected=s['lambda_enter_expected'],
 lambda_enter_expected_ci95=q['rate_estimate']['ci95_s_inv'],physical_source_time_span_s=s['physical_source_time_span_s'],
 source_mean_interarrival_s=s['source_mean_interarrival_s'],accepted_mean_interarrival_s=s['accepted_mean_interarrival_s'],
 diameter_unit='um',**{'source_D_'+k:v for k,v in s['source_diameter_um'].items()},**{'enter_D_'+k:v for k,v in s['enter_diameter_um'].items()},
 accepted_clearance_m=s['accepted_clearance_m'],Qacc_by_D=q['rows'],rejection_reasons=s['rejection_reasons'],
 point_basin_counts=m['point_basin_counts'],point_basin_scope='FIRST_30_ACCEPTED_EVENTS_DESCRIPTIVE_ONLY',smoke30_status=m,
 confirmed_bugs_found=0,confirmed_bugs_fixed=0,existing_scientific_files_modified=[],code_change_ledger=[],
 legacy_tests_passed=130,new_tests_passed=35,total_final_tests_passed=165,final_tests_failed=0,final_tests_skipped=0,
 baseline_regression=read(D/'baseline_regression.json'),final_regression=t,old_results_changed=0,protected_published_files=46670,
 worker_determinism=dict(workers=[1,3,6],proposals=100000,accepted_births=19221,bitwise_identical=True,ledger_sha256=read(D/'worker_replay_1.json')['proposal_ledger_sha256']),
 open_inlet_confirmed=True,artificial_inlet_cap_size_exclusion=False,
 runtime_s=dict(baseline=read(D/'baseline_regression.json')['runtime_s'],inlet_100k_6workers=s['runtime_s'],inlet_replay_1worker=read(D/'worker_replay_1.json')['runtime_s'],
  inlet_replay_3workers=read(D/'worker_replay_3.json')['runtime_s'],smoke30_6workers=m['smoke_runtime_s'],synthetic_40k=syn['runtime_s'],
  final_portable_test_time=t['runtime_test_s'],final_50_test_time=t['additional_tests']['runtime_test_s']),
 large_sample_status='READY_FOR_P9A4_LARGE_SAMPLE_REVIEW',final_status=stage,formal_500_executed=False,pushed=False,merged=False,
 server_root=read(D/'server_deployment.json')['server_root'],report_built_utc=datetime.now(timezone.utc).isoformat())
write('data/final_summary.json',json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False))

write('P9A4_POPULATION_INLET_REVIEW_ZH.md',f'''
# 一句话结论

新 P9-A.4 direct Poisson thinning 已通过入口联合分布、全量 worker 重放和 30 泡 smoke 审核，达到 **{stage}**。绝对时间使用明确的名义稳态浓度假设；本轮没有启动正式 500，也没有 push 或 merge。

实际基线 `{b['ACTUAL_HEAD']}` 与预期相同。新 worktree 为 `{ROOT}`，分支 `{b['branch']}`。旧科学源码、流场、合同、报告和输出均未修改。

# 一秒钟里微泡是怎么来的

`lambda = C Q`：入口每秒通过的体积 Q 乘以每立方米中的泡数 C，给出平均 source proposal 数。本例 C 为 D≤4 µm 群体的 `{c['C_modeled_D_le_4um_m3']:.12g} m⁻³`，Q 为 `{c['Q_in_m3_s']:.16g} m³/s`，所以平均每秒产生 `{s['lambda_source_s_inv']:.9f}` 个候选泡。它是期望值，不是每秒强行放入固定数量。

# 为什么时间不是等间隔

每次间隔独立取 `-log1p(-U)/lambda_source`，然后按 source ID 顺序累加。100,000 次实测平均间隔 {s['source_mean_interarrival_s']:.6f} s，方差 {s['source_variance_interarrival_s2']:.6f} s²。固定窗口 source Fano={s['fano_source_audit']:.6f}，entering Fano={s['fano_accepted_audit']:.6f}，仅作描述性检查。旧 deterministic scheduler 完全保留。

# 一颗 source bubble 怎么决定大小

唯一正式 source 为 `SONOVUE_D_LE_4UM_CONDITIONAL`；使用冻结 SonoVue 直方图和原分箱内均匀分布，条件化而非把超过 4 µm 的泡裁到 4 µm。每个 source ID 只抽一次 D，范围 {c['D_min_m']*1e6:g}–4 µm。没有 4 µm 堆积。新合同 SHA：`{summary['source_contract_sha']}`。

# 为什么入口位置按血流量抽

直接复用原 `InletFluxSampler`，在整个真实 FEM INLET cap 上按 `(u·n)+/Q_in` 抽一次。原 P1 正向通量裁剪与三角采样不改，不预先限制到某一尺寸的可行区域，不用未来出口标签。

# 为什么只检查一次

现实中的泡不会因为这个位置放不下就瞬间传送到另一个位置。每个候选泡的时间、D、位置生成后，调用原 `FiniteSizeAdmission.check(active={{}})` 一次。失败就保留一次拒绝；不重抽位置或尺寸，也不延后重新插入。这里只代表独立单泡轨迹群体，不包含 simultaneous many-body pair exclusion。

# rejected source proposal 是什么

它是到达入口、但不满足真实 WALL 或原 lower-handoff 条件的物理候选事件，**不是 solver failure**。全部 100,000 行记录保存在 `data/inlet100k/proposal_ledger.jsonl`：19,221 accepted，80,779 rejected；拒绝率 80.779%。其中 WALL_REJECTED=80,672，WALL_NEARFIELD_REJECTED=107。每行都有 source ID、时间、D、位置和原因；只有 accepted 才有连续 particle ID。

# source distribution 和 entering distribution 为什么不同

大泡需要更大的空间，所以实际进入群体向较小直径偏移。必须分别解释两种 PDF。

{dtable}

独立控制样本预测的 entering CDF 与实测最大差为 {g['entering_size_cdf_error']:.6f}；16 个尺寸×位置联合单元的最大差为 {max(x['absolute_error'] for x in g['joint_cells']):.6f}。这些统计验证支持目标联合分布，不等于对任意分辨率分布的一致性证明。

# finite-size accessible flux 是什么

Q_acc(D) 是尺寸 D 的球在真实壁距与 handoff 规则下能进入的那部分流量。它影响进入概率，不参与预先改变 source 位置分布。独立 20,000 个全通量位置估计如下，800 个位置/尺寸组合与原 checker 的判断完全一致。

{qtable}

零次接受的 MC 点不单独证明数学上严格为零；保留其上置信限。表中 CI 为逐尺寸点态区间，使用共同位置，因此不同尺寸估计相关。

# 和 Method B 有什么区别

B 固定 anchor 后反复抽原完整 SonoVue 尺寸，accepted PDF 是位置条件化分布。旧 replay 保持。本次 500 个入口诊断请求得到 306 个 accepted，共 108,361 次尺寸抽样；这是 legacy 尝试预算下的统计，不是 physical Poisson acceptance。

# 和 Method C 有什么区别

C 对全入口不可能的尺寸重抽，随后为可行尺寸重抽位置；其时钟定义为 deterministic entering `C_MB Q`。本次 500 个入口诊断请求得到 500 个 accepted，542 次 source 尺寸抽样和 859 次位置抽样。P9-A.4 保留 source rejection，进入率自然降低。B/C 的 500 只是入口统计，未运行 500 条轨迹。

# open inlet 为什么不是墙

真实 INLET 的 179 个三角不在 WALL BVH 中，cap/WALL 面重叠数为 0。球心位于 cap，球可以向上游跨出；不要求整个球进入体积。开放/封闭 cap 对照、rim overlap、真实 cap 场查询及边/顶点 owner tests 均通过。没有 inlet-plane 距离≥radius 的人工限制。[几何证据](OPEN_INLET_AUDIT_ZH.md)。

# absolute concentration 是不是 ground truth

不是。历史 8.5e12 m⁻³ 来自全泡 bolus dose/blood-volume 名义锚点，不是 D≤4 µm 浓度，也不是连续灌注实测值。新 C_modeled = 8.5e12 × F_original(4 µm)，F={c['F_original_4um']:.12f}。已澄清 total 与 conditional 语义；把这一锚点维持为恒定稳态属于 `MODEL_ASSUMPTION`。[追溯与文献](CONCENTRATION_PROVENANCE_ZH.md)。

# 当前发现了哪些真实代码 bug

在本轮指定科学路径和永久测试覆盖范围内，CONFIRMED_BUG=0，修复=0，existing scientific files 修改=0。deterministic timing、retry、旧 flow contract 都按 legacy/model scope 处理，没有把新模型偏好包装成旧 bug。新开发测试/分析脚本的执行问题及一个旧回归工具的相对路径限制另列于 [ledger](BUG_FIX_LEDGER_ZH.md)，没有隐去失败日志。未对所有 legacy 任意输入作无 bug 保证。

# NEW Network-H0 inlet audit

实际 NEW SHA `{c['flow_SHA']}`；OLD SHA 被新入口拒绝。100k source 的物理时钟跨度 {s['physical_source_time_span_s']:.6f} s（约 {s['physical_source_time_span_s']/86400:.4f} 天），这是恒定模型的大样本统计窗口，不是建议一次真实 bolus 实验运行这么久。

实测进入率 {s['lambda_enter_measured']:.9f} s⁻¹，独立预测 {s['lambda_enter_expected']:.9f} s⁻¹，预测 MC 95% CI [{q['rate_estimate']['ci95_s_inv'][0]:.9f}, {q['rate_estimate']['ci95_s_inv'][1]:.9f}]。两者差为合并标准误 {g['rate_discrepancy_in_combined_standard_errors']:.3f} 倍。accepted 平均间隔 {s['accepted_mean_interarrival_s']:.6f} s。acceptance 95% CI [{s['acceptance_ci95'][0]:.6f}, {s['acceptance_ci95'][1]:.6f}]。

accepted clearance mean/p10/p50/p90 为 {', '.join(f'{s["accepted_clearance_m"][k]*1e9:.3f}' for k in ['mean','p10','p50','p90'])} nm。按直径和入口三角的拒绝统计见 `data/rejections_by_diameter.csv`、`data/rejections_by_inlet_triangle.csv`。1/3/6 worker 的完整 100k ledger 和 19,221 个 births 逐字节一致；拒绝事件 checkpoint/restart、重排合并和重复 ID 拒绝测试通过。

# 30 smoke

首 30 个 accepted births 在独立 Vast.ai 目录、6 CPU worker、每库单线程下运行，P9-A.1 动力学源码冻结。结果 **28 completed / 2 stationary / 0 solver fail**。stationary ID 4、8 有连续 32 步 rank-3 接触与 KKT 支持，不解释成生理滞留。穿透、handoff 违规、入口逃逸、NaN/Inf 均为 0（按原 roundoff 容差）。

MB 出口 O1/O2/O3=0/6/22；同一已接受样本的 point basin 为 1/12/14，noexit=3（30 s point horizon）。只用于描述，既不证明 population split，也不用于 quota 或接受门控。[完整 smoke 审核](SMOKE30_REVIEW_ZH.md)。

# P1 local conservation limitation

NEW P1/P1+VMS 合法 root section max≈3.454212%、RMS≈1.789924% 的局部守恒偏差保留。本轮只使用真实 FEM INLET，不修改 CFD、Taylor-Hood、grad-div 或 P9-A.3 internal section。

# 是否可以进入正式 500

**建议人工审核 source/entering PDF、80.779% rejection 的物理含义、名义浓度与进入率、两条 supported stationary 及 bug audit 后，再决定启动正式 500。** 本阶段仅 `READY_FOR_P9A4_LARGE_SAMPLE_REVIEW`，没有自动授权或执行 500。最终 165 tests pass、0 fail、0 skip；46,670 个旧发布文件 SHA 不变。

[图件总览](OPEN_RESULTS.html) · [数学合同](P9A4_MATHEMATICAL_CONTRACT_ZH.md) · [代码审核](CURRENT_CODE_LOGIC_AUDIT_ZH.md) · [复现](REPRODUCE.md) · [机器摘要](data/final_summary.json)
''')

write('P9A4_MATHEMATICAL_CONTRACT_ZH.md',r'''
# P9-A.4 数学与计算合同

设真实入口 S，n 指向血管内，q(x)=(u·n)₊，Q=∫S q(x)dA。尺寸 PDF f 是冻结 SonoVue 在 D≤4 µm 的条件 PDF，C 是该条件群体的 number/m³。稳态源是带独立标记的齐次 Poisson 点过程，强度测度为 `C f(D) q(x) dt dD dA`。

每个 source k：Δtₖ=−log1p(−Uₖ)/(CQ)，tₖ=ΣΔtᵢ；独立取 Dₖ~f、xₖ~q/Q。使用原真实 WALL 和 lower handoff，定义固定可行指示 A(D,x)。active={}，因此 A 不依赖其他粒子或过去事件。在这个假设下，确定性标记筛选仍为 Poisson thinning：

`Q_acc(D)=∫S q(x) A(D,x)dA`

`p_acc=(1/Q)∫f(D) Q_acc(D)dD`

`lambda_enter=C∫f(D)Q_acc(D)dD=lambda_source p_acc`

`p_enter(D,x)=f(D)q(x)A(D,x)/∫f(d)Q_acc(d)dd`

其尺寸边缘正比于 f(D)Q_acc(D)，不是原 f(D)。在失败后重抽位置或尺寸会重新条件化分布；若仍沿用 source CQ 时钟，则改变进入过程。因此一次 source 只调用一次 check，不调用 attempt 或 MethodCSource.event，不删除拒绝行。

正式生成不需要先知道 Q_acc。独立 MC 诊断只解释结果，不改变 source RNG、位置分布或 acceptance。无未来出口信息；无 outlet quota；不把有限尺寸 population split 强行等同 fluid split。

## 圆管解析 benchmark

Poiseuille u(r)=u_max(1−r²/R²)，Q_total=πu_maxR²/2。设 a=D/2、h=2 nm、R_c=R−a−h，x=max(0,min(1,R_c/R))。对半径 0…R_c 积分 2πr u(r)dr，得到 `Q_acc/Q=2x²−x⁴`。

全源径向 CDF `G(r)=2(r/R)²−(r/R)⁴`。固定 D 的 accepted 径向 CDF 在 r≤R_c 时为 G(r)/G(R_c)，随后为 1。离散尺寸概率 pⱼ 的进入边缘为 pⱼG(R_c,j)/ΣpᵢG(R_c,i)。tests 使用 R=2 µm、尺寸 0.6/2.4 µm、概率 0.4/0.6、lambda_source=5 s⁻¹、固定 seed 2594 和 40,000 个 proposals。

数值结果见 [synthetic receipt](data/synthetic_poiseuille_audit.json)：source 径向 CDF 最大误差 0.00495；两尺寸 accepted 比例 0.920753、0.284066，解析值 0.922048、0.293057；条件径向最大误差 0.001420、0.007332；小尺寸 entering 概率实测 0.683363、解析 0.677163；进入率实测 2.682858、解析 2.723267 s⁻¹。固定宽容差 tests 全通过。Fano 只作描述，不作为严格随机 gate。

## 数值实现与可恢复性

RNG 使用 PCG64/SeedSequence(master_seed,source_event_id,role)，独立角色 940/941/942/943 对应 arrival/D/x/orientation；原 STREAMS 顺序不变。worker 只生成独立 marks，主进程按 source ID 顺序 float64 累加时间并分配 accepted particle ID。身份绑定 NEW flow、浓度、源合同、原 114 源码 manifest 与两份新源码 SHA。

完整100k在1/3/6 worker下 ledger SHA 均为 `b89b019874b8d142f31f227ce45ed20f2d0b5301595b1b8dc3dd92cdad3c8cd7`；accepted births 也逐字节一致。这是相同软件/平台数值合同下的保证，不宣称跨任意 NumPy/VTK/CPU 版本必然 bitwise 相同。

checkpoint 保存全部 proposals、下一个 source ID、下一个 accepted ID、累计时间和身份哈希；完成 manifest 最后写入。恢复检查文件哈希并重放有序 reducer，不抽新 RNG。partial、篡改、丢失拒绝、身份不符、重复/缺失 ID、已有输出目录均被拒绝。这里只提供 population checkpoint；不迁移旧轨迹 checkpoint 的源码身份。

负/非有限浓度、非正/非有限正向 Q、rate under/overflow 显式错误。C=0 表示没有下一个 source event；极小正 rate 若超出有限 float64 时间范围则停止报错，不写 Inf。CDF 使用原冻结条件化实现，没有 4 µm 端点裁切堆积。时间精度耗尽会报错而不静默重复正间隔。

## 适用边界

恒定流、恒定条件源浓度、球形、独立单泡、既有 WALL+handoff 模型。P9-A.1 动力学及场插值不变。真实 continuous infusion 浓度、PK/clearance、多泡排斥或随时间变化流场不在本合同内；这些因素可能破坏当前齐次独立 thinning 前提。
''')

write('OPEN_INLET_AUDIT_ZH.md',f'''
# Open inlet 审核：PASS

WALL 只包含真实 WALL。INLET cap 179 个面与 WALL 的 GlobalNodeID 面集合交集为 0；perimeter 顶点属于真实 WALL，保留 rim 碰撞。面积 7.756795804427715e-12 m²，等效直径 3.1426516127 µm。cap 略非平面，最大平面残差 1.48136e-11 m；查询使用实际三角和 owner tetra，不把拟合平面当碰撞面。

永久 tests 已运行并通过：

| 验证 | 结果 |
| --- | --- |
| sphere_straddles_open_cap_accepts | 开放圆管入口允许同一球向上游跨出 |
| same_sphere_with_solid_cap_rejects | 加入实体 cap 的对照拒绝 |
| rim_wall_overlap_rejects | 真实 rim/WALL overlap 拒绝 |
| inlet_cap_not_in_wall_bvh | 实际 cap/WALL 面交集 0 |
| center_on_inlet_field_sample_is_valid | 179 cap 面心均有效 |
| edge_and_vertex_boundary_ownership_stable | 重复及反序查询 owner 稳定 |

上述测试分布于 `test_open_inlet_straddling_sphere_accepts.py`、`test_real_wall_overlap_rejects.py`、`test_inlet_cap_not_solid.py` 和旧 open-cap 测试。真实100k proposals 没有 CENTER_OUTSIDE 拒绝；WALL 80672、lower handoff 107。

birth center 位于真实 cap，只要求 center 有效及球不碰真实 WALL/不越过 lower handoff；不要求整球 inside，也没有 distance_to_inlet_plane≥radius。较大尺寸接近零 acceptance 来源于狭窄真实壁面，不是把 cap 封死。

30 smoke 全部保存线段按向外方向与真实 INLET 相交事后复判：逃逸 0。初始 cap 上向内运动不误记为逃逸。此观察只证明本次保存轨迹与原连续 WALL 证书的结果；不构成任意上游球扫掠几何的证明，也没有新增 upstream-to-cap 动力学。

[几何机器证据](data/open_inlet_geometry.json) · [最终测试日志](logs/final_new_and_legacy_pass.txt) · [smoke 逐轨迹指标](data/smoke30_metrics.json)
''')

write('CONCENTRATION_PROVENANCE_ZH.md',f'''
# 浓度语义追溯：total / conditional 已解决，稳态绝对值为模型假设

旧 `PARTICLE7_INLET_POPULATION_V0.json` 将 8.5e12 m⁻³ 明确命名为 `NOMINAL_WELL_MIXED_POST_BOLUS_ANCHOR`。原论文的动物实验部分给出 32.6 g 平均体重、2×10⁷ 个 microbubbles 的 bolus，未将该剂量限定为 D≤4 µm。该文没有测量本血管的连续灌注稳态浓度。[原论文 2.6.1](https://pmc.ncbi.nlm.nih.gov/articles/PMC12445601/)

血容量取历史合同的 72 mL/kg；独立机构指南也给出小鼠这一估计值。它是估算而非本实验逐鼠测量。[University of Michigan guideline](https://az.research.umich.edu/animalcare/guidelines/guidelines-blood-collection)

计算：V_blood=0.0326×72=2.3472 mL，C_total_nominal=2e7/2.3472×1e6={c['concentration_provenance']['unrounded_nominal_total_m3']:.12g} m⁻³。继承原舍入锚点 C_total=8.5e12 m⁻³。不能把这个全部 SonoVue 浓度直接写成 conditional modeled 浓度。

冻结分布的 F_original(4 µm)={c['F_original_4um']:.15f}。所以 **C_modeled_D≤4={c['C_modeled_D_le_4um_m3']:.15g} m⁻³**。源分布条件化与源强度扣除被排除尺寸质量必须同时发生。数学等价于从全体源删除 D>4 的 marks，再对真实入口做有限尺寸 thinning。

该文 bolus dose 支持“全部泡的名义剂量锚点”这一语义；本文进一步假设连续灌注把此名义浓度维持为常量。这一步是显式 `MODEL_ASSUMPTION`，不是从论文推得的生理 ground truth。未新增 pharmacokinetics、clearance、泡寿命或循环再入模型。若实际浓度变化，当前绝对时间需重新建模；固定稳态 C 的缩放不改变 normalized source/entering joint law，但按比例改变 lambda 和时间尺度。

100k 模型时钟约 8.85 天只用于统计精度，不表示 bolus 在活体持续这一时间。100k 入口生成无轨迹推进；30 smoke 为独立事件的局部 residence integration，非连续多日 many-body 模拟。

原始分布 SHA `{c['source_original_distribution_sha256']}`；sampler SHA `{c['source_original_sampler_sha256']}`；源合同 SHA `{c['source_original_contract_sha256']}`。原文全文通过 Europe PMC fullTextXML 核验，保存核验字段和 XML SHA，未将论文全文拷入仓库：[核验记录](data/concentration_literature_verification.json)。
''')

write('SMOKE30_REVIEW_ZH.md',f'''
# 首 30 个 accepted births smoke：PASS WITH SUPPORTED STATIONARY

在 inlet gate、115 portable、35 new 和15 legacy tests 全部通过后，取100k群体的前30个 accepted births；不按出口筛选。cohort SHA `aebee7de17212e461b50b8a47f20ee79452161d4d3253268bfc3328316af06d2`。服务器独立目录 `{summary['server_root']}`，Python `/root/particle8_2_runs/env/bin/python`，6 CPU workers、OMP/BLAS/MKL/NUMEXPR=1。没有 GPU 轨迹内核。

结果 completed=28、stationary=2（ID 4、8）、solver fail=0、time-limit=0。MB 出口 O1/O2/O3=0/6/22。两条 stationary 的最后32个已接受求解均有 contact rank=3、平移速度在原 KKT budget 内、非负 multipliers，因此按原审核规则是受约束 stationary；不把其解释成生理 trapping，也没有更改接触模型以消除 stationary。

保存样本 {m['saved_sample_count']}，连续 handoff 证书叶节点 {m['continuous_certificate_count']}。原诊断及独立真实 WALL 最近距离复算均给 penetration=0、handoff violation=0。逐线段入口向外相交复查给 escape=0。NaN/Inf=0，出口分类复判全部一致。

最小真实 wall gap={m['minimum_gap_m']:.17g} m；最小 g_nf={m['minimum_g_nf_m']:.17g} m。后者微小负值约 −1.88e−20 m，小于冻结 roundoff 2.0094114631e−17 m，报告原始值而非裁为0；因此无超容差违规。连续证书验证贯穿原 accepted path；额外独立最近距离检查覆盖保存状态。

point tracer 在同一已接受的30个初始中心事后运行：O1=1、O2=12、O3=14、noexit=3（ID18、25、28达到30 s point horizon）。这不是100k accepted群体的精确 basin split，更不能以 fluid O1/O2/O3 比例作 quota。finite-size MB 和 point 的路径模型不同。

MB运行 {m['smoke_runtime_s']:.3f} s；保存数据审核 {m['analysis_runtime_s']:.3f} s。部署204个输入/源码的 SHA，回传154个输出逐一 SHA 验证。P9-A.1、积分入口、初始速度规则、dt=0.00025 s、MB horizon=1.5 s 均复用原 Network runner；新输入只有入口群体/身份。原 integrate_admitted 对既定 birth 还有启动一致性验证，不重抽、不重记 source acceptance、不改变 ledger。

[逐轨迹指标](data/smoke30_metrics.json) · [接触与连续证书证据](data/smoke30_contact_audit.json) · [汇总](data/smoke30_summary.json) · [图08](figures/Figure_08_smoke30_trajectories.png) · [服务器日志](logs/server_smoke30.txt)

结论为可供人工 large-sample review；本轮没有运行正式500。
''')

write('BUG_FIX_LEDGER_ZH.md','''
# Bug / change ledger

本轮指定科学逻辑审核中 CONFIRMED_BUG=0、fixed=0；existing scientific files 修改0。[code_change_ledger.json](data/code_change_ledger.json) 是空数组。115 baseline 全部通过才开始新实现。旧 deterministic clock、attempt retry、Method B/C 条件化属于 legacy 模型行为；P9-A.4 通过新增两个模块建立独立合同，没有重写旧模块。

新增 `continuous_infusion.py` 提供counter RNG、Poisson间隔、单次proposal、ordered ledger/checkpoint；新增 `population_inlet_p9a4.py` 显式加载NEW输入、浓度合同与CPU worker。仅复用 TruncatedSonoVue 类，没有迁移或改写其文件，保护旧 imports/SHA。新 tests 22个指定文件共35项，另有conftest及包命名文件。

开发阶段失败和处理（全部保留原日志，不计为已有科学代码 bug）：

| 事件 | 原因 | 最小处理与验证 |
| --- | --- | --- |
| new_initial | 新测试写错 BridgeParticle.radius_m 属性 | 改为原 radius；新35通过 |
| inlet_gate | 新gate将numpy.bool_交给标准JSON编码 | 新辅助脚本显式bool；成功gate与失败partial均保留 |
| smoke_analysis | 新分析脚本误把带identity的cohort字典当事件列表 | 读取events字段；同一已保存结果审核通过，未重跑动力学 |
| final_new_and_legacy | 新tests从全局conftest导入，与其他测试目录重名 | 新目录加包名并用相对导入；组合50项通过 |
| final_full_regression | 原portable工具在子进程临时cwd中解释相对--output，找不到JUnit | 使用绝对输出路径；完整115通过；未改原工具 |
| build_reports | 新报告汇总脚本误读deployment的remote_root键 | 改为实际server_root键；未影响任何计算或数据 |

最后一项是复现工具调用/相对路径可移植性限制，非Particle科学算法bug。本轮将该限制写入复现命令，不放宽科学身份，也不修改受保护旧工具。未执行所有历史模块的任意无效输入、多写者破坏性竞争或跨版本checkpoint迁移；这些边界不宣称已证明安全。

最终165项通过、0失败/0跳过。46,670个发布旧文件、40个原network输入、114个原科学源文件、114个server-bundle源文件及OLD/NEW场哈希全部匹配。92个生成pyc只在portable临时副本的旧源码合同中排除，原文件/合同没有编辑。保护日志：[final_protection_manifest.txt](logs/final_protection_manifest.txt)。
''')

# Update only this stage's preliminary audit with exact execution evidence.
tests_by_index={
 1:'test_poisson_inverse_transform.py; test_poisson_fixed_seed_replay.py; source-level legacy clock audit',
 2:'test_rejected_event_never_retries_position.py; test_rejected_event_never_retries_size.py; rejection call-count tests',
 3:'test_old_method_c_replay_unchanged.py; test_truncated_sonovue_no_4um_spike.py; legacy particle9a2_inlet',
 4:'test_old_method_b_replay_unchanged.py; 500-request inlet-only diagnostic',
 5:'test_flux_position_sampler_unchanged.py; complete 100k worker replay; baseline Network flux checks',
 6:'legacy particle9a2_inlet capacity tests; source inspection confirms absent from new production path',
 7:'conditional CDF parity / source SHA / exact D-radius tests',
 8:'test_inlet_cap_not_solid.py: actual cap centroids, shared vertices and edges',
 9:'test_inlet_cap_not_solid.py; test_open_inlet_straddling_sphere_accepts.py; test_real_wall_overlap_rejects.py',
 10:'unchanged original analyze.one + analyze_smoke.py all 30 saved paths; outlet matches and 0 outward INLET crossings',
 11:'source-level original checkpoint audit; new test_restart_preserves_next_source_event.py; old checkpoint migration intentionally not performed',
 12:'source-level original saved-interval restore audit; no claim of new cross-version trajectory checkpoint replay',
 13:'30 isolated outputs complete +154 returned file hashes; source-level writer review; new population partial/overwrite tests',
 14:'portable P9-A.1 tests + 30 smoke / 14592 continuous certificates; no equations changed',
 15:'test_old_flow_rejected_for_p9a4.py; test_new_contract_uses_network_h0_flow.py',
 16:'test_worker_independent_proposal_ledger.py; full 100k 1/3/6 byte equality',
 17:'test_old_method_c_replay_unchanged.py; ordered merge and duplicate-ID tests',
 18:'portable Network paired regression + unchanged runner used for smoke; original source SHA match',
 19:'test_new_contract_uses_network_h0_flow.py; concentration_literature_verification.json; original contracts SHA unchanged',
 20:'test_restart_preserves_next_source_event.py / duplicate-ID tests; simultaneous external writers NOT TESTED, risk retained'}
audit=['# 当前科学代码逻辑审核：最终证据','',
 '源级审查在新实现前于 '+a['audit_before_implementation_utc']+' 完成；原始20条发现和源码SHA保存在 data/code_logic_audit.json。指定科学路径未确认真实bug，existing科学源码修改0。以下区分源码审查、执行过的测试及尚未覆盖的风险，不把legacy行为归为bug。','',
 'baseline115通过；最终portable115 + 新35 + legacy/open-cap15 =165通过，0失败/0跳过。测试命令与开发执行问题见 REPRODUCE.md、BUG_FIX_LEDGER_ZH.md。']
for i,item in enumerate(a['findings'],1):
 audit+=['',f'## {i}. {item["file"]}','']
 for key,label in [('function_or_class','function/class'),('classification','分类'),('current_behavior','当前行为'),('expected_behavior','预期行为'),('minimal_reproducer','最小复现/检查方法'),('severity','严重性'),('modified','是否修改')]:audit.append(f'- **{label}**：{item[key]}')
 audit.append('- **永久测试/执行证据**：'+tests_by_index[i])
audit+=['','## 共同单位、边界与限制','',
 'SI: m、m/s、Pa、m³/s、number/m³；直径与半径转换在source标记处明确完成。INLET inward、OUTLET outward，GlobalNodeID→index减1，float64稳定owner；q=0/backflow按原positive clipping。新流场metadata及实文件双重SHA，禁止OLD替代与默认旧路径fallback。',
 '永久测试覆盖inverse-transform的NaN/Inf、负/零/极小rate、CDF端点4 µm、正负通量裁剪、ID排序和拒绝事件恢复；负/非有限浓度及非正/非有限Q的构造器拒绝路径做了源级审查，未另计为已运行参数化测试。真实rim/WALL与open cap对照测试通过。active={}刻意排除many-body pair conflict；旧pair算法只做源级审查，未宣称本轮验证simultaneous人口模型。',
 '输出单ID单worker及独占新目录是当前合同；外部同时写同目录/同ID没有锁的风险保留。旧轨迹checkpoint应在旧版本身份下重放；本轮不放宽旧身份来迁移。数学采样测试覆盖固定seed和宽松预设容差，不以出口配额判断correctness。',
 '浓度语义已解析为full anchor×F4，但恒定稳态绝对浓度仍是假设。P1局部守恒、局部平壁动力学近似、supported stationary的生理解释属于模型限制。既有科学代码0个CONFIRMED_BUG只描述本轮证据范围，不能推导全历史代码无缺陷。']
(R/'CURRENT_CODE_LOGIC_AUDIT_ZH.md').write_text('\n'.join(audit)+'\n')
write('data/code_logic_audit_execution.json',json.dumps(dict(original_audit_sha256=digest(D/'code_logic_audit.json'),evidence=tests_by_index,confirmed_scientific_bugs=0,scope='Specified active scientific paths and selected permanent regressions; untested risks retained'),indent=2,ensure_ascii=False))

write('REPRODUCE.md',r'''
# P9-A.4 reproduce / review

本轮只在新worktree开发，没有push/merge。原worktrees只读。所有下列命令从 `/home/lzy/projects/ulm_particle_population_inlet_p9a4` 执行。Python 3.13.11、NumPy 2.5.3、SciPy 1.18.1、PyVista 0.49.0；服务器Python3.13.15、相同NumPy/SciPy/PyVista。记录源/输入SHA而不假设跨任意软件版本bitwise一致。

```bash
cd /home/lzy/projects/ulm_particle_population_inlet_p9a4
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONPATH="$PWD/particle_3d/src"
export SONOVUE_ROOT="$PWD/sonovue_size_distribution_v0"
P9A4_PY=/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python
P9A4_CHECKDIR=$(mktemp -d /tmp/p9a4-review-XXXXXX)
"$P9A4_PY" sync_metadata/network_h0_particle_20260925/run_checks.py --output "$P9A4_CHECKDIR/portable"
"$P9A4_PY" -m pytest -q -p no:cacheprovider --junitxml="$P9A4_CHECKDIR/new_legacy.xml" particle_3d/tests/particle9a4_population_inlet particle_3d/tests/particle9a2_inlet particle_3d/tests/particle8_2a/test_open_inlet_not_solid_wall.py
"$P9A4_PY" sync_metadata/network_h0_particle_20260925/verify_snapshot.py --scientific-inputs
```

`run_checks.py --output` 必须绝对路径；其子进程工作目录为临时副本。该旧工具将5份路径manifest只在临时副本重定位，并只在临时源码保护合同中排除92条生成pyc，所有原科学源SHA仍强制匹配。新tests在独立包命名空间，避免多conftest重名。

## 入口生成与checkpoint

本次已经完成：`scripts/run_inlet_audit.py prepare`（独立20k Qacc + immutable contract）、`generate --workers 6 --count 100000 --label inlet100k`、`replay --workers 1`、`replay --workers 3`。不要在已有合同上重跑prepare，不要覆盖现有receipt/ledger。重新生成入口审计使用新label，例如：

```bash
"$P9A4_PY" particle_3d/reports/particle9a4_population_inlet/scripts/run_inlet_audit.py generate --workers 6 --count 100000 --label inlet100k_review_repeat
```

只生成proposals与accepted births，不跑轨迹。新label如果存在会显式拒绝。需要检验1/3/6而不写既有receipt，可在Python中加载同一冻结contract并调用 `make_source(load_new_environment(root), contract)`、`generate(source,100000,workers=n)`，对每行 `canonical_bytes` 的串接计算SHA，与 `data/worker_replay_1.json` 比较。source ledger SHA预期 `b89b019874b8d142f31f227ce45ed20f2d0b5301595b1b8dc3dd92cdad3c8cd7`。生产源身份包含两份新源码SHA；源码若改变不能假装原ledger仍是同一科学版本。

`PopulationLedger.restore(folder, source.identity)` 读取完整manifest、ledger和state，验证后得到同样next IDs/time。`generate(source, extra_count, ledger=restored, workers=n)` 从下一个source ID延续，保存到全新checkpoint目录。永久测试覆盖53个事件（包括拒绝）checkpoint后75个事件续跑，对照128个不中断事件。缺失/篡改/identity不同的checkpoint拒绝，不能删除拒绝行来减小文件。

## 实际输入与保护

NEW输入：`particle_3d/reports/network_derived_flow_mb_validation_v1/server_bundle/inputs/NEW.vtu`。
SHA：`064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4`。
source合同：`particle_3d/contracts/P9A4_CONTINUOUS_INFUSION_V1.json`。
原科学源码保护：`particle_3d/reports/network_derived_flow_mb_validation_v1/data/code_snapshot_manifest.json`。

两份新增科学模块直接使用原 full-positive-flux sampler、原 conditional SonoVue sample、原 admission.check。原 legacy B/C、STREAMS、FrozenFEMField、WALL builder、P9-A.1 equations、OLD contract均未编辑。旧轨迹checkpoint由原版本重放，不为新增module放宽其科学身份。

## 审核脚本与输出

`gate_inlet.py` 读取data中的既定宽容差、100k及独立控制样本。真实审计容差在gate执行前写入，但不是在所有初步aggregate计算前盲定；synthetic永久测试容差预先固定。`compare_legacy.py` 的500表示入口请求，不是trajectory；比较accepted分布时应保留full vs conditional source与分母差异。

`run_smoke30.py` 只允许 `/workspace/particle9a4_population_inlet_<UTC>`，核对deployment manifest、gate、固定30cohort与NEW身份后创建独占输出目录。它复用原Network mb_job；原积分器启动时再次验证既定birth一致性，不重抽source。`deploy_smoke.py` 本轮创建的remote_root见data/server_deployment.json。已完成30 smoke，不要为审核图表重复运行或把它扩大为正式500。服务器操作遵循已保存的logs/server_agent_guide.txt。

`analyze_smoke.py` 只读保存结果，复用原analyze.one并额外独立壁距/入口相交审核，不推进轨迹。`render_figures.py inlet/smoke` 只读数据，英文白底300dpi PNG+PDF；会拒绝覆盖已有图件。`build_reports.py` 从审计产物构建报告，也拒绝覆盖最终交付；重建草稿需另存目录/版本。原阶段输出永不覆盖。

完整100k ledger约88 MiB，是重要科学证据；全部reject行保留。服务器原始30轨迹、metadata、压缩audit、point结果均已回传并保存SHA。报告支持的结论来自这些完整数据。图件目录figures/，主报告和机器摘要见OPEN_RESULTS.html；不需要旧可视化目录参与本轮。

## 本轮执行结果

baseline115 passed（87.810 s）；smoke前portable115、新35、legacy15通过；最终portable115 + new35 + legacy15=165通过、0fail/skip。完整100k worker1/3/6逐字节一致；30 smoke28completed/2supported stationary/0solverfail。46,670旧文件SHA不变。失败的开发测试/辅助脚本日志保留于logs，原因见BUG_FIX_LEDGER_ZH.md。没有正式500、没有push或merge。
''')

reports=['P9A4_POPULATION_INLET_REVIEW_ZH.md','CURRENT_CODE_LOGIC_AUDIT_ZH.md','P9A4_MATHEMATICAL_CONTRACT_ZH.md','OPEN_INLET_AUDIT_ZH.md','CONCENTRATION_PROVENANCE_ZH.md','BUG_FIX_LEDGER_ZH.md','SMOKE30_REVIEW_ZH.md','REPRODUCE.md','data/final_summary.json']
cards=[]
for p in sorted((R/'figures').glob('Figure_*.png')):
 name=p.stem;cards.append(f'<section><h2>{html.escape(name)}</h2><a href="figures/{name}.png"><img src="figures/{name}.png" alt="{html.escape(name)}"></a><p><a href="figures/{name}.pdf">PDF</a> · <a href="figures/{name}.png">300 dpi PNG</a></p></section>')
write('OPEN_RESULTS.html','<!doctype html><html lang="zh"><meta charset="utf-8"><title>P9-A.4 review</title><style>body{background:white;color:#17212b;font-family:system-ui;max-width:1200px;margin:40px auto;padding:0 20px}img{width:100%;height:auto}section{margin:45px 0}a{color:#1767a6}</style><h1>P9-A.4 population inlet review</h1><p>100,000 source → 19,221 accepted；30 smoke：28 completed / 2 supported stationary / 0 solver fail。</p><p>READY_FOR_P9A4_LARGE_SAMPLE_REVIEW；未运行正式500。名义连续灌注浓度为模型假设。</p><ul>'+''.join(f'<li><a href="{x}">{x}</a></li>' for x in reports)+'</ul>'+''.join(cards)+'</html>')
print('REPORTS_COMPLETE',stage)
