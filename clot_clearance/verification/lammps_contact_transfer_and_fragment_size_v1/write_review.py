"""Write review from completed actual cases; no solver launch."""
from pathlib import Path
import sys,json,csv
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from pd_clot.lammps_fragmentation.targeted_contact_transfer import O,R
from pd_clot.lammps_fragmentation.config import dump

s=json.loads((O/'ALL_CASE_SUMMARIES.json').read_text());cross=json.loads((O/'CROSS_PATTERN_CONTACT_DECISION.json').read_text());conditional=json.loads((O/'CONDITIONAL_CASE_DECISION.json').read_text());perf=json.loads((O/'PERFORMANCE_REAL_CASES.json').read_text());viz=json.loads((O/'VISUALIZATION_QA.json').read_text())
main=['GAUSSIAN_100_MODIFIED','OPPOSING_100_D2','GAUSSIAN_200_D2'];modified=[v for v in s.values() if v['contact_model_id']=='MODIFIED_1P0DX'];formal=[k for k in s if (O/'CASES'/k).exists()]
status=dict(CONTACT_1P0_TRANSFER=cross['status'],CONTACT_1P0_PRODUCTION_DEFAULT='CONDITIONAL_EXPLORATORY_CANDIDATE' if cross['status']=='SUPPORTED_FOR_TESTED_PATTERNS' else 'NOT_PROMOTED',SECOND_PATTERN_TEST='PASS',STRONGER_LOAD_TEST='PASS',OLD_MAP_REFRESH='TARGETED_POINTS_ONLY' if conditional['selected'] else 'DEFERRED',DETACHED_SIZE_200UM='OBSERVED' if any(q['thresholds']['200']['EVER_OBSERVED'] for q in modified) else 'NOT_OBSERVED',DETACHED_SIZE_300UM='OBSERVED' if any(q['thresholds']['300']['EVER_OBSERVED'] for q in modified) else 'NOT_OBSERVED',PERSISTENT_LARGE_DETACHED_FRAGMENT='OBSERVED' if any(q['thresholds']['300']['persistent_membership_episodes'] for q in modified) else 'NOT_OBSERVED',WHOLE_CLOT_DETACHMENT='OBSERVED' if any('WHOLE_CLOT_DETACHMENT' in q['dynamic_failure_classification'] for q in modified) else 'NOT_OBSERVED',PULVERIZATION_AT_200PA='OBSERVED' if 'PULVERIZATION' in s['GAUSSIAN_200_D2']['dynamic_failure_classification'] else 'NOT_OBSERVED',LPS_PROGRESSIVE_SOFTENING='NOT_IMPLEMENTED',VES='DEFERRED',B1_B2='TERMINATED_BY_PROJECT_POLICY_CHANGE',PD2_FRACTURE_DEVELOPMENT='STOPPED',REAL_BUBBLE_STREAMING='NOT_CONNECTED',BIOLOGICAL_QUANTITATIVE_VALIDATION='NOT_ESTABLISHED',PATIENT_SPECIFIC='NOT_ESTABLISHED',PRIMARY_USE='QUALITATIVE_FRAGMENTATION_SCREENING',EXPLORATORY_ONLY=True,PERSISTENT_200UM_REFERENCE_EPISODES='OBSERVED',PERSISTENT_LARGE_DEFINITION='DETACHED d_eq >=300 um with unchanged membership for >=3 consecutive saved frames; exact durations separately reported',NEW_FORMAL_SOLVES=len(formal),HISTORICAL_SOLVER_RERUNS=0,MATERIAL_RECALIBRATIONS=0,OFFICIAL_CONTACT_DEFAULT_UNCHANGED=True,OLD_MAP_UNCHANGED=True,STOP='COMPLETE_STOPPED_FOR_HUMAN_REVIEW')
dump(O/'FINAL_STATUS.json',status)
def pct(x):return f'{100*x:.4f}%'
def ms(x):return '未出现' if x is None else f'{x*1000:.4f}'
def title(k):return {'GAUSSIAN_100_MODIFIED':'Gaussian 100 / 1.0（复用）','GAUSSIAN_100_NATIVE':'Gaussian 100 / 原生1.35','OPPOSING_100_NATIVE':'Opposing 100 / 原生1.35','OPPOSING_100_D2':'Opposing 100 / 1.0','GAUSSIAN_200_D2':'Gaussian 200 / 1.0','GAUSSIAN_31P622_D2':'Gaussian 31.6228 / 1.0','OPPOSING_31P622_D2':'Opposing 31.6228 / 1.0'}.get(k,k)
metrics='\n'.join(f"| {title(k)} | {pct(s[k]['broken_bulk_fraction'])} | {pct(s[k]['broken_interface_fraction'])} | {pct(s[k]['debris_fraction'])} | {pct(s[k]['largest_attached_fraction'])} | {s[k]['final']['failure_mode']} |" for k in ['OPPOSING_100_NATIVE',*main,*conditional['selected']])
sizes='\n'.join(f"| {title(k)} | {d} | {z['MAX_SIMULTANEOUS_COUNT']} | {100*z['PEAK_DETACHED_MASS_FRACTION']:.4f}% | {z['AT_FINAL_FRAME']} | {ms(z['FIRST_APPEARANCE_TIME'])}–{ms(z['LAST_SEEN_TIME'])} | {z['FRAME_PERSISTENCE_COUNT']} | {1000*z['SAVED_GRID_PRESENCE_DURATION_S']:.3f} | {z['persistent_membership_episodes']} |" for k in main for d,z in s[k]['thresholds'].items())
contacts='\n'.join(f"| {title(k)} | {ms(s[k]['first_bulk_break_time_s'])} | {ms(s[k]['first_interface_break_time_s'])} | {s[k]['max_active_contact_pairs']} | {s[k]['min_nonbonded_distance_over_dx']:.5f} | {s[k]['peak_clot_KE_J']*1e9:.3f} |" for k in ['OPPOSING_100_NATIVE',*main])
tractions='\n'.join(f"| {title(k)} | {s[k]['tau_budget_Pa']:.3f} | {s[k]['reference_peak_local_traction_Pa']:.3f} | {s[k]['actual_peak_local_traction_Pa']:.3f} | {s[k]['peak_actual_mean_traction_budget_Pa']:.3f} |" for k in main)
runtimes='\n'.join(f"| {title(k)} | {s[k]['wall_seconds']:.3f} | {s[k]['original_postprocess_wall_s']:.3f} | {s[k]['size_postprocess_wall_s']:.3f} | {perf['cases'][k]['solver_performance']['LAMMPS_sections']['Output']['avg_s']:.4f} |" for k in formal)
chronology='\n'.join(f"| {title(k)} | {mode} | {ms(z['first_time_s'])} | {1000*z['saved_grid_duration_s']:.3f} | {'是' if z['at_final'] else '否'} |" for k in main for mode,z in s[k]['dynamic_failure_classification'].items())
review=f'''# LAMMPS_CONTACT_TRANSFER_AND_FRAGMENT_SIZE_V1 中文审查

本轮已完成，停止等待人工审查。**1.0dx 接触修正在已测 Gaussian / Opposing 100 Pa 模式上获得迁移支持；Gaussian 200 Pa 未再次粉碎。出现过短暂的 ≥200 µm 真正游离块，但没有 ≥300 µm 游离块，也没有任何 ≥200 µm 游离块持续到终态。** 所有尺度来自参考节点体积，结果保持 EXPLORATORY_ONLY。

严格执行 2 个必做 + 2 个条件补点，共 4 次新增 LAMMPS 半球求解。历史原生 1.35dx、Gaussian 100 / 1.0dx 均未重算；没有强度校准、材料调参、接触因子扫描或新 C++ 开发。使用上轮独立变体，其可执行文件、patch 和源码身份完全匹配。

固定配置为 2971 clot + 1203 三层 anchors，140555 bulk / 16983 interface / 40395 wall bonds；R=0.5625 mm、dx=50 µm、horizon=3.01dx。K=33333.333333333336 Pa、G=940 Pa、rho=1000 kg/m³、数值阻尼 500/s、s00_bulk=0.15、s00_interface=0.140624859375、alpha=0 均不变。实际脉宽为 1.1605 ms（2.000260435389031 T_cross），dt=0.5 µs、4660 步、终时 2.33 ms、每 20 步保存、MPI=2，全部按已有取整输入保留。近边界预算精确复用 31.622776601683793 Pa。

## 1. Opposing 100 Pa 是否仍抑制异常粉碎？

是，CONTACT_1P0_TRANSFER=SUPPORTED_FOR_TESTED_PATTERNS。Opposing 首次 bulk 断键在原生和修改版均为 **0.1660 ms**；此前 17 个保存帧的位置、速度逐位一致。界面首断由 0.3935 ms 变为 0.4860 ms。末态 debris 从原生约 99.5624% 降至 **4.9815%**，减少 **94.5809 个百分点**；修改版仍有真实活跃接触，保留 **94.2444% 附着主体**。这不是只让碎屑飞得更近，而是 surviving-bond graph 中 bulk 断键大幅减少。

| 工况（预算 Pa） | bulk 断键 | interface 断键 | 1–2 节点 debris | 最大附着主体 | 终态旧分类 |
|---|---:|---:|---:|---:|---|
{metrics}

原 ≥3 节点 RESOLVED_FRAGMENTATION 分类保持不变，只表示离散图上的可辨认分量，不能据它称出现了大游离碎块。

| 工况 | bulk 首断 ms | interface 首断 ms | 峰值活跃接触对 | 全程最小非键合间距/dx | 峰值物理节点动能 nJ |
|---|---:|---:|---:|---:|---:|
{contacts}

三例主比较的修改版均无保存帧距离 <0.5dx 的非键合深度重叠，因此影响比例和持续时间均为零；没有 lost atoms、missing bonds、NaN 或提前终止。活跃对按实际保存位置、原生排斥距离和 surviving-bond 排除规则重建，包含真实输入中仍启用的 pair 力；这是保存帧诊断，不是逐次 pair force 日志，不能据此推断闭合接触功。有限刚度接触的 0.72–0.84dx 局部压近也不能称绝对不可穿透。所有采样门槛限于工程筛查，非生物精度标准。

## 2. Gaussian 200 Pa 是更多大块、整块脱落，还是再次粉碎？

是更早起裂、较多小碎屑并保留主体，**没有全局粉碎或整块脱粘**。首次 bulk 断键 0.1325 ms，interface 首断 0.4690 ms；末态 bulk 断键 8.1605%、debris 6.3952%、最大附着主体 93.2683%。最终最大游离分量只有 4 节点，d_eq≈98.47 µm（离散体积统计，非可信真实微粒尺寸）。

该强载荷下最大曾出现的游离块为 **48 节点、225.45 µm**；≥200 µm 出现在 0.430–0.480 ms，随后继续碎裂。Gaussian 100 的最大值为 49 节点、227.01 µm，≥200 窗口反而更长。因而这一个加倍载荷探针没有增加工程大块输出，不能据少量点声称载荷越大就产生越大的游离块。200 Pa 峰值动能 82.671 nJ，出现在约 0.580 ms；没有粉碎时间点，故没有“接触修正失败导致再次粉碎”的证据。

| 工况 | 动态类别 | 首次保存时刻 ms | 保存网格累计持续 ms | 末帧仍存在 |
|---|---|---:|---:|---|
{chronology}

## 3. 49-node 例子的尺寸和生命周期是多少？

独立读取原生 broken-bond events，逐帧重建 bulk graph，再检查 surviving interface 到 wall 的连接，确实重现 **0.500 ms 的 49 节点 DETACHED 分量**。d_eq=(6×49×dx³/π)^(1/3)=**227.005 µm**。相同 49 个原节点的完整集合在 0.500–0.540 ms 共 5 帧存在，观测跨度 **40 µs**；采用保存网格左保持计时至下一帧为 **50 µs**。0.550 ms 的下一帧分成 47+1+1 节点；47 节点子块不再冒充原来完整的 49 节点块。

整个 ≥200 µm 游离尺寸区间延续到 0.680 ms，但其中包含后续不同成员集合；“49 节点块持续 0.19 ms”会混淆原块和子块。身份采用相邻帧粒子集合交叠并严格以成员集合不变为一个 episode；每次分裂所有子块都另开 episode，父子关系、节点 ID 和交叠数完整保存，未发现多父合并歧义。

1/2/3/5/34/49/113/114 节点及 attachment/wall 排除检查均通过。3 节点的等效球直径约 89.47 µm，不达到两档尺度，也缺乏充分几何分辨率；≥3 节点只能保持旧图分类含义。34 才达到 200 µm；**114 才达到 300 µm，113 不够**。单节点体积为 125000 µm³，单节点的球等效直径也不等于网格间距 50 µm，不能称“真实 50 µm 血栓微粒”。

## 4. 两个尺度的峰值、终值和持续情况？

下面仅统计**游离 clot-only 分量**，不把大的附着主体计入游离质量；质量比例分母为初始 clot 质量 3.71375e-7 kg。200 µm 是接近 [Chueh 等取栓实验 >200 µm 分组](https://pubmed.ncbi.nlm.nih.gov/26919972/) 的参考量，本轮用 ≥200；300 µm 是本项目 dx 相关工程定义，均非普适临床危险阈值。

| 工况 | 门槛 µm | 最大同时块数 | 峰值游离质量比例 | 末态块数 | 首次–末次出现 ms | 出现帧数 | 左保持累计存在 ms | ≥3 帧相同成员 episodes |
|---|---:|---:|---:|---:|---|---:|---:|---:|
{sizes}

所有上述终态对应尺度质量比例均为零；两项 31.6228 Pa 补点从未出现这两档游离块。≥200 的相同成员 episodes 总数分别为 Gaussian100 的 10、Opposing100 的 20、Gaussian200 的 4；其中连续 ≥3 帧者为 3、7、1。最长相同成员游离观测跨度分别为 40、60、20 µs。这些 episode 数包含分裂后子块，**不是独立栓子的总数**；真实同时块数最多为 1、2、1。

PERSISTENT_LARGE_DETACHED_FRAGMENT 按预声明的 ≥300 µm + 相同成员连续 ≥3 保存帧定义，结果 NOT_OBSERVED；≥200 参考尺度下确实有短暂持续 episodes，另行保存，不能把它们隐藏成“从未出现大于 200 的块”。每帧 COM、速度、相对参考与相对首次拓扑脱离的 COM 位移均保存；拓扑 DETACHED 不代表已远离壁面。

## 5. 有没有细长结构被等效球直径误认为立体大块？

本轮三例主比较及补点没有 ≥114 节点的游离候选，因此没有该尺度下需要排除的细长游离假阳性，也不能反过来证明所有较小碎块都是立体的。代码对所有 ≥114 节点候选保留参考位置 PCA 三轴、节点中心跨度、包含原 dx 立方节点体积后的包围盒跨度，并对最短轴 <2dx 标 GEOMETRICALLY_THIN_CANDIDATE，未过滤任何分量。微型线状节点检查确认标记能生效，未运行任何力学 synthetic benchmark。

等效体积直径与最长粒径不同。1–2 节点碎屑明确为 NUMERICALLY_UNRESOLVED_DEBRIS；[microtripsy 原始实验](https://pmc.ncbi.nlm.nih.gov/articles/PMC4528908/) 的主要细碎屑尺度远低于本模型网格，不能由这里的离散颗粒反推真实粒径分布或临床安全性。

## 6. 哪些旧 map 边界仍可解释，哪些不兼容？

37 个旧结果、类别和图均保留，仍可作为 **NATIVE_1P35DX 历史模型自身**的响应记录。不能把这些标签直接移给 MODIFIED_1P0DX。本轮门槛在两例必做完成后通过，先写 CONDITIONAL_CASE_DECISION.json 冻结输入，才各运行一次两个 31.6228 Pa 补点。

- Gaussian / 31.6228 Pa / 2T：旧 PULVERIZATION → 新 NO_FAILURE。新模型仍有约 27 条 bulk、22 条 interface 键断开，属于旧类别容许的低比例断键，**不是零损伤**；无游离块、无 debris。
- Opposing / 31.6228 Pa / 2T：旧、新均 NO_FAILURE，本轮新结果零断键。
- 两模式 / 100 Pa / 2T：旧 PULVERIZATION → 新 RESOLVED_FRAGMENTATION；新类别仍是旧 ≥3 节点定义。
- Gaussian / 200 Pa / 2T：新模型 RESOLVED_FRAGMENTATION；不存在同预算旧 native 数据，不推断或补跑。

原 Gaussian 2T 下 [10,31.6228] 的 NO_FAILURE→PULVERIZATION 区间不适用于新模型；原 Opposing [31.6228,100] 的上端粉碎标签也不适用。新模型两模式只有 31.6228 的 NO_FAILURE 与 100 的 RESOLVED 两个相邻已测点，提示分类变化在其间，但**不拟合阈值、不宣称全区间单调，也不把它称为粉碎边界**。Gaussian 200 仍 resolved。新模型其它预算/脉宽/涡旋均 UNKNOWN。第三张图按接触模型分栏，只画真实点。

## 7. 1.0dx 能否升级为有限范围的探索筛查默认候选？

可以标记 **CONTACT_1P0_PRODUCTION_DEFAULT=CONDITIONAL_EXPLORATORY_CANDIDATE**：两种已测模式在 100 Pa 下均保留真实 pair 排斥、自然断键/运动，断键前可比，严重重叠未出现且异常级联显著减轻；参数全部冻结。200 Pa 的正常结果提供一个额外强载荷点，仍不是对全部载荷的保证。

这是候选资格，不是已测真实 clot 接触半径，也没有覆盖官方默认 1.35 或把旧 map 重标为新模型。前一轮 bulk/interface 强度兼容性证据在相同 executable 下沿用，本轮没有重新校准。Liu 整血 bulk / 猪皮剥离界面来源不变，界面仍是 PORCINE_SKIN_PEELING_MODEL_PROXY，非真实人体血管壁；小应变 LPS 映射不等于完整材料复现。

## 8. 真实 case 和后处理、渲染成本？

| 本轮新 case | solver wall s | 原图/分类后处理 s | 独立尺寸/生命周期后处理 s | LAMMPS Output 平均 s |
|---|---:|---:|---:|---:|
{runtimes}

新增 solver wall 合计 **{perf['new_solver_wall_sum_s']:.3f} s**，全部 2 MPI ranks、4660 步。已有 Gaussian100 修改版 solver 的历史 25.769 s 不计入本轮新增成本。完整原生事件/轨迹独立重分析（含只读 native 对照）合计 **{perf['size_analysis_wall_sum_s']:.3f} s**；部分与新求解重叠，不能与各项墙钟简单相加称总历时。三张图 {viz['figure_wall_s']:.3f} s，短片渲染 {viz['movie_render_wall_s']:.3f} s；本轮没有 build、安装新包或额外性能实验。每例真实 CPU、RSS、I/O 输出段和数据体积见 PERFORMANCE_REAL_CASES.json。

同预算不等于相同局部峰值；预算仅为初始参考构型总绝对切向力除以原 exposed-node 权重面积。当前位置 provider 随离场衰减，运动中不重新归一化。

| 工况 | 参考平均预算 Pa | 参考局部峰值 Pa | 保存帧实际局部峰值 Pa | 保存帧实际平均预算峰值 Pa |
|---|---:|---:|---:|---:|
{tractions}

200 Pa 完全属于工程探索预算，非实测微泡载荷或临床破裂阈值。rho 和数值阻尼也仍属工程假设。

## 9. 是否已有证据要求 softening、VES 或 PeriLab？

没有。接触修正已让两个模式从全局粉碎变为有限损伤，强载荷点也正常。本轮未形成 ≥300 µm 游离块是诚实的尺寸结果，不能据此就加入软化、加大载荷、重调 s00、借用 VES 参数或换求解器。LPS_PROGRESSIVE_SOFTENING=NOT_IMPLEMENTED，VES=DEFERRED；B1/B2 路线终止、PD2 fracture 停止。停止等待人工审查，不自动展开下一轮。

## 10. 将来怎样分别优化碎裂程度与大游离块的风险相关指标？

未来连接真实 streaming 场后，应分别记录目标效果与风险代理：

- 碎裂/移除效果：最大附着主体的减少、总附着质量下降、脱离质量和其离场运动；断键比例作为机制量，不单独当疗效。需要真实流场/输运输入后才谈恢复流通，不把当前节点运动冒充再通。
- 大块相关惩罚：≥200/≥300 游离块的峰值与末态质量比例、最大等效尺寸、持续时间及离开局部区域的质量/通量；后者需要流向和下游截面信息。两档尺度是报告/工程代理，不能直接转换为患者栓塞概率。
- 分辨率约束：单独输出 NUMERICALLY_UNRESOLVED_DEBRIS 与接触数值异常标志。不要通过奖励 1–2 节点颗粒化来假装获得安全微碎屑，也不要简单最大化最大碎块或全部游离体积。

因而适合保留分开的多目标输出或约束关系，先明确希望移除多少附着主体，同时限制大游离块质量与持续时间。当前真实 microbubble/CFD 还未连接，PZT/surrogate 本轮不实施；输入场单位、坐标、时间、有效域和输运证据缺口不能靠这些代理统计填补。

## 交付和停止

三张英文科学图的 PNG/PDF 与一个同步 MP4/GIF 位于 results/lammps_contact_transfer_and_fragment_size_v1。短片包含 Gaussian100 历史修改版与本轮两例主 case，同一物理时刻、同一相机，全部节点和 surviving bonds 保留，无插值、无位移或碎屑放大；灰色 anchors、蓝色附着主体、按真实游离尺寸着色的分量分开。上排是固定原邻域，下排是固定完整粒子云范围。播放 fps 不代表实验/超声时间。

PER_FRAME_DETACHED_COMPONENTS.csv、FRAGMENT_LIFETIMES.csv、FRAGMENT_SPLIT_RELATIONS.csv、FRAGMENT_SIZE_SUMMARY.csv 和每例 membership_and_shapes.json 保存完整统计；图中生命周期只展示最大八个 episodes，完整数据不删。HASH_MANIFEST.json 与轻量 evidence.zip 区分入包证据与原位保存的大文件。旧两轮 verification/results 全目录与原脚本、可执行文件均按启动时清单校验；不 commit/push/PR，不改写历史 map。

**最终状态：COMPLETE_STOPPED_FOR_HUMAN_REVIEW / EXPLORATORY_ONLY。**
'''
(O/'REVIEW_LAMMPS_CONTACT_TRANSFER_AND_FRAGMENT_SIZE_V1.md').write_text(review)
perf['figure_wall_s']=viz['figure_wall_s'];perf['movie_render_wall_s']=viz['movie_render_wall_s'];dump(O/'PERFORMANCE_REAL_CASES.json',perf)
print('REVIEW_AND_FINAL_STATUS_WRITTEN',status)
