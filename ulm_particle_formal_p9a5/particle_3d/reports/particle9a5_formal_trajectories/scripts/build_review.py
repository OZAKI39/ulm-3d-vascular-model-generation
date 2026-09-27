"""Plain-language review and an HTML entry that starts with scientific outcomes."""
from pathlib import Path
import json,html,sys
ROOT=Path(__file__).resolve().parents[4];R=ROOT/'particle_3d/reports/particle9a5_formal_trajectories';D=R/'data'


def number(x,digits=4):return '未观察到 / 不适用' if x is None else f'{x:.{digits}g}'
def yes(x):return '是' if x else '否'
def span(d,scale=1,unit=''):
 if not d['n']:return '没有可统计的样本'
 return f"均值 {d['mean']*scale:.4g}{unit}，中位数 {d['median']*scale:.4g}{unit}，P10–P90 {d['p10']*scale:.4g}–{d['p90']*scale:.4g}{unit}，范围 {d['minimum']*scale:.4g}–{d['maximum']*scale:.4g}{unit}"


def main():
 s=json.loads((D/'final_summary.json').read_text());n=s['final_formal_N'];core=s['core500'];full=s['full_formal'];bench=s['worker_scaling_results'];points=json.loads((D/'point_complete.json').read_text())
 dt_comparison=json.loads((D/'dt_comparison_24.json').read_text())
 missing=[o for o in ['O1','O2','O3'] if core[o]==0];final_missing=[o for o in ['O1','O2','O3'] if full[o]==0]
 if n==500:
  why='前500颗已经自然覆盖三个出口，因此在满足最低500颗要求后停止扩样。' if s['all_outlets_naturally_observed'] else '已完成固定CORE500，但仍有出口未观察到；按生产前声明的资源上限停止，未更改入口或挑选轨迹。'
 else:
  why=f"CORE500 中 {'、'.join(missing)} 尚未出现，因此继续使用后续连续 accepted births，最终累计到 {n} 颗。"
  why+= '此时三个出口已全部自然出现，按预定规则停止。' if s['all_outlets_naturally_observed'] else f"此时达到生产前声明的资源上限，{'、'.join(final_missing)} 仍未观察到。"
 first='；'.join(f"{o}：{number(s['first_'+o+'_at_N'],8)}" for o in ['O1','O2','O3'])
 safety='、'.join(f'{label}={s[k]}' for k,label in [('penetration_count','穿透'),('handoff_violation_count','handoff违规'),('inlet_escape_count','入口逃逸'),('nan_inf_count','NaN/Inf'),('unclassified_corruption_count','未分类损坏')])
 outcome_table='| 类别 | CORE500 数量 | CORE500 比例 | FULL 数量 | FULL 比例 | 比例差（百分点） |\n| --- | ---: | ---: | ---: | ---: | ---: |\n'
 for key in ['O1','O2','O3','stationary','censored','failed']:
  outcome_table+=f"| {key} | {core[key]} | {100*core['fractions'][key]['fraction']:.3f}% | {full[key]} | {100*full['fractions'][key]['fraction']:.3f}% | {100*s['core_full_fraction_difference'][key]:+.3f} |\n"
 uncertainty='| 出口 | CORE500 95% 区间 | FULL 同时区间 |\n| --- | --- | --- |\n'
 for o in ['O1','O2','O3']:
  a=core['fractions'][o];b=full['fractions'][o];uncertainty+=f"| {o} | {100*a['lower']:.3f}%–{100*a['upper']:.3f}% | {100*b['lower']:.3f}%–{100*b['upper']:.3f}% |\n"
 expansion='| 检查点 N | 本轮新增 | O1 | O2 | O3 | 三出口都有 |\n| ---: | ---: | ---: | ---: | ---: | --- |\n'
 for c in s['population_checkpoints']:expansion+=f"| {c['N']} | {c['added']} | {c['O1']} | {c['O2']} | {c['O3']} | {yes(c['all_outlets_naturally_observed'])} |\n"
 scaling='| workers | 秒 / 固定24条 | tracks/s | CPU核等效利用 | 峰值PSS GiB | 失败 |\n| ---: | ---: | ---: | ---: | ---: | ---: |\n'
 for c in bench['configurations']:scaling+=f"| {c['workers']} | {c['wall_seconds']:.2f} | {c['tracks_per_second']:.5f} | {c['cpu_utilization_core_equivalents']:.2f} | {c['peak_tree_pss_bytes']/2**30:.3f} | {c['solver_failures']} |\n"
 diameter='\n'.join(f"- {k}（n={v['n']}）：{span(v,1,' µm')}。" for k,v in s['diameter_by_outcome'].items())
 final_tests=s['final_tests'] or {'tests':'待最终回归','failures':'待核验','errors':'待核验','skipped':'待核验'}
 report=f'''# 一句话结果

完成 {n:,} 条按 P9-A.4 accepted 顺序确定的正式微泡轨迹；O1/O2/O3 分别为 {full['O1']}/{full['O2']}/{full['O3']}，有接触支持的静止 {full['stationary']} 条，长驻留删失 {full['censored']} 条，求解失败 {full['failed']} 条。

**最终状态：`{s['final_status']}`。** 三出口是否全部自然观察到：**{yes(s['all_outlets_naturally_observed'])}**。

# 这次为什么至少要跑500颗

30颗smoke用于确认流程能够安全运行，不能替代正式样本库。CORE500是在任何新轨迹计算前冻结的前500个accepted births，成员、直径、位置、方向、出生时刻与source_event_id都保留；以后版本之间可以直接使用同一个队列比较。

# 最后实际跑了多少颗

正式 FULL_FORMAL_N = **{n}**，最低500颗要求：**已满足**。CORE500始终单独保留。{why}

# 为什么500之后还继续增加

{why} 扩样只追加原P9-A.4已审核列表中紧接着的accepted IDs，不重新开始source population。

# 500颗时与最终N时三个出口分别有多少

{outcome_table}

比例的分母是该队列的全部正式成员，包括静止、删失和失败，不只计算成功到出口的微泡。CORE和FULL是嵌套样本，表中的差值是描述性比较，不是两个独立随机组的显著性检验。

{uncertainty}

CORE500使用固定样本、逐类别Wilson 95%区间。FULL的停止数量依赖出口覆盖，所以采用Clopper–Pearson区间并对6个类别×最多13个预定检查点作Bonferroni调整，整体覆盖率至少95%，允许按既定覆盖规则停止；它会更保守。没有观察到某出口并不等于其概率为零，至少出现1颗也不代表概率已估计准确。

# 有没有为了让某个出口出现而挑微泡

**没有。** 没有按future outlet、point basin或O1友好区域选择，没有出口配额、重抽位置、重抽直径，也没有删除不想要的结果。Point tracer是在最终MB cohort冻结后运行的事后解释工具。

# 微泡数量是怎样顺序增加的

{expansion}

首次出现的accepted前缀位置：{first}。首次三个出口同时存在的最小前缀为 {number(s['all_outlets_observed_at_N'],8)}；按预定批次正式作出停止决定的检查点为 {number(s['all_outlets_observed_at_checkpoint_N'],8)}。这两个量不同于并行任务的先后完成顺序。

# 模拟时间为什么从旧值延长

从历史P8.2A/PPT500的500行trajectory_catalog实际核实：21条达到上限的记录均为**1.5 s**，并与PPT_REVIEW/PPT_VALIDATION中的队列数量和状态核对。本轮初始maximum age为**3 s**，名义dt按用户后续要求调整为**{s['dt_s']:.3f} s（1.0 ms）**；接触处理仍可使用更小的内部子步。

本轮实际最大使用horizon为**{s['maximum_horizon_used_s']:g} s**；3→6延长 {s['horizon_extension_3_to_6_count']} 条，6→12延长 {s['horizon_extension_6_to_12_count']} 条。已经到出口或已获接触支持的静止轨迹不会为了填满3秒而继续运动。若正常运动的轨迹达到边界，使用同一个存活stepper继续；延长前的数组前缀另存并逐点核对。总provider计算预算随授权horizon按比例放宽；dt、力学、contact、handoff及拒绝/进展保护逻辑均保持原样。

# 哪些微泡顺利走出去了

共 {full['completed']} 条实际穿过正式出口。保存的最后一段位移再次与原始出口三角形独立核验。通行时间仅为真正完成者定义：{span(s['transit_stats'],1,' s')}。

# 哪些微泡没有出去

有接触支持的静止 {full['stationary']} 条；达到12秒仍正常运动而被删失 {full['censored']} 条；求解失败 {full['failed']} 条。所有成员的记录都保留在trajectory_catalog及原始轨迹目录中；失败也保留最后有效状态、debug证据和错误原因，不从分母剔除。

# stationary和时间不够有什么区别

“静止”要求原P9-A.1停滞信号及最后32个接受的求解步满足三方向接触支持、非负乘子和KKT速度预算。单纯时间到3/6/12秒不会被改叫静止。12秒仍未完成且仍正常运动的记录归为LONG_RESIDENCE_CENSORED，不能据此声称发生生理性堵塞。

# 不同大小的微泡有什么差别

{diameter}

这些是同一预先确定队列中的结果分组，不能把分组差异直接解释成独立的粒径因果效应；出生位置也会影响路径。

# 有限尺寸微泡与point tracer有什么不同

Point出口/未完成数量为 `{s['point_basin_counts']}`。完整转移表见Figure10和final_summary。Point使用原有CPU原生P1诊断实现，保留30秒/2毫米诊断预算；它与MB的3/6/12秒分类不同，不能用point是否走出代替MB是否完成。它只用于事后说明相同出生中心对应的流体路径。

# 近壁过程是否安全

{safety}。最小原始壁间隙统计：{span(s['wall_gap_stats'],1e9,' nm')}。近壁暴露时间统计：{span(s['nearwall_stats'],1,' s')}，近壁定义为原始壁间隙/半径≤0.1。所有保存位置还独立检查了原始壁面和入口外逃；连续handoff证书与accepted诊断保留。

# 服务器实际用了多少计算资源

服务器为原Vast实例，独立目录 `{s['server_root']}`。CPU可见16个逻辑核，容器配额**{bench['CPU_quota_cores']:.2f}核**。固定同一24条轨迹比较1/2/4/6/8/12/16 workers，科学结果SHA一致，并与同一服务器上原P9-A.1积分器以相同1.0 ms步长重新计算的24条参考轨迹逐点一致。旧0.25 ms结果单独归档，仅用于步长变化比较。

{scaling}

正式使用**{s['production_workers']} workers**，各数学库线程数为1。12和16 workers作为超配对照；正式选择配额范围内、没有内存压力且吞吐量最高的实测配置。当前轨迹kernel没有经过验证的GPU实现，因此**微泡轨迹和point轨迹都在CPU计算**，未为使用GPU改写动力学。

生产批次墙钟合计 {s['runtime_seconds']:.2f}秒（{s['runtime_seconds']/60:.2f}分钟）；worker benchmark合计 {s['benchmark_runtime_seconds']:.2f}秒；point诊断 {s['point_runtime_seconds']:.2f}秒。24条已完成、源身份和SHA匹配的benchmark轨迹在正式CORE500中直接复用。

生产前预估5000条全部走到12秒的保守存储成本约 {bench['worst_case_5000_storage_bytes']/2**30:.2f} GiB、计算成本约 {bench['worst_case_5000_runtime_seconds']/3600:.2f}小时；存储按实测每单位物理时间的最大数据量外推，包含自适应子步和被拒绝尝试的诊断；计算按实测每次provider调用成本外推，并非严格解析上界。生产前声明的RESOURCE_LIMITED_N_MAX为**{bench['resource_limited_N_max']}**，原因：`{bench['hard_limit_reason']}`。实际短轨迹成本通常更低，但上限不会在看到想要的出口后临时改写。

# 1.0 ms 步长修订与差异

名义步长由用户明确指定为**1.0 ms**；旧0.25 ms批次在指令到达后停止并独立归档，其结果不计入本轮正式N。CORE500 SHA和流场SHA保持不变，内部安全细分仍可使用更小子步。固定前24个相同出生样本的新旧步长比较中，出口标签变化为**{dt_comparison['outlet_changes']} / 24**；逐条时间与空间差异见 `data/dt_comparison_24.csv`。

同一1.0 ms步长下，不同worker和原积分器参考的逐点一致性用于验证实现与并行调度；这**不等于时间离散误差已经收敛**。本轮没有把1.0 ms结果宣称为0.25 ms结果的等价替代。步长变更、早期验证PNG/PDF和原始归档位置均可追溯。

# 当前结果可以说明什么

它给出固定NEW Network-H0流场和当前P9-A.1有限尺寸模型下，可复现、可审核的独立微泡路径库，以及固定CORE500和顺序扩样后的结果分布。三出口是否出现是自然观测的事实，不是入口目标配额。

# 当前结果仍然不能说明什么

不能将稀有出口仅出现1颗当作准确概率估计；不能将未观察到当作绝对不可能；不能把接触支持静止直接称作真实血管堵塞。流场仍是已接受的P1/P1 VMS Network-H0参考场，已有局部守恒限制没有在本轮重算。模型没有MB–MB或RBC耦合，恒定绝对浓度平台仍是P9-A.4的模型假设。轨迹库中的出生时刻并不意味着动画中高亮的微泡同时注入。

# 验证、SHA与复现附录

- baseline commit：`{s['baseline_commit']}`。
- worktree：`{s['new_worktree']}`；branch：`{s['new_branch']}`。
- NEW flow SHA：`{s['flow_sha256']}`。
- P9-A.4 source contract SHA：`{s['source_contract_sha']}`。
- CORE500 cohort SHA：`{s['core_cohort_sha']}`。
- FULL cohort SHA：`{s['final_formal_cohort_sha']}`。
- baseline tests：{s['baseline_tests']}。
- 新永久测试：{s['new_tests']}，保留于`particle_3d/tests/particle9a5_formal_trajectories/`。
- final tests：{final_tests}。
- 当前保留科学输入/结果的变化：`{s['old_files_changed']}`；详见 retained_science.sha256 与最终保护核验。历史快照已按用户要求从 WSL 清理，不再要求其完整保留。
- CORE500与FULL原始cohort、batch manifest、完成标志、逐文件SHA、原始npz和压缩audit均保留。

[HTML总览](OPEN_RESULTS.html) · [完整摘要](data/final_summary.json) · [逐轨迹CSV](data/trajectory_catalog.csv) · [原始轨迹](outputs/formal/) · [图件](figures/) · [动画](animations/) · [测试日志](logs/)

# 为什么最后跑的是{n}颗，而不是500颗？

{why} 500是最低正式样本量，后续只能按accepted顺序追加。全过程没有为了某个出口改变入口population或挑选微泡。
'''
 (R/'P9A5_FORMAL_TRAJECTORY_REVIEW_ZH.md').write_text(report)
 cards=[('最终正式数量',str(n)),('500最低要求','已满足'),('三出口自然出现',yes(s['all_outlets_naturally_observed'])),('O1',str(full['O1'])),('O2',str(full['O2'])),('O3',str(full['O3'])),('Completed',str(full['completed'])),('有支持的 stationary',str(full['stationary'])),('Long-residence',str(full['censored'])),('Solver failure',str(full['failed'])),('名义 dt','1.0 ms'),('实际最大 horizon',f"{s['maximum_horizon_used_s']:g} s"),('穿透',str(s['penetration_count'])),('Handoff violation',str(s['handoff_violation_count'])),('入口逃逸',str(s['inlet_escape_count'])),('NaN / Inf',str(s['nan_inf_count'])),('Production workers',str(s['production_workers'])),('正式计算墙钟',f"{s['runtime_seconds']/60:.1f} min"),('计算设备','CPU')]
 cardhtml=''.join(f'<div class="card"><span>{html.escape(k)}</span><strong>{html.escape(v)}</strong></div>' for k,v in cards)
 dt_html='<section><h2>步长调整验证</h2><p>名义 dt = 1.0 ms，按用户后续指令调整。原0.25 ms结果独立归档，不计入本轮正式N。相同步长的实现一致性不等于时间收敛证明。</p><img loading="lazy" src="figures/Stage_B_dt1ms_validation.png"><p><a href="figures/Stage_B_dt1ms_validation.pdf">早期验证 PDF</a> · <a href="data/dt_comparison_24.csv">24条配对步长差异</a> · <a href="data/dt_revision.json">修订与归档记录</a></p></section>'
 images=''.join(f'<section><h2>{p.stem.replace("_"," ")}</h2><a href="figures/{p.name}"><img loading="lazy" src="figures/{p.name}"></a><p><a href="figures/{p.stem}.pdf">下载 PDF</a> · <a href="figures/{p.name}">300 dpi PNG</a></p></section>' for p in sorted((R/'figures').glob('Figure_*.png')))
 videos=''.join(f'<section><h2>{p.stem.replace("_"," ")}</h2><video controls preload="metadata" src="animations/{p.name}"></video><p>使用保存轨迹；出口代表为该出口 completed 中 particle_id 最小者。总览按归一化轨迹进度展示，不代表同时注入。</p></section>' for p in sorted((R/'animations').glob('*.mp4')))
 links='<a href="P9A5_FORMAL_TRAJECTORY_REVIEW_ZH.md">中文完整审核报告</a> · <a href="data/trajectory_catalog.csv">逐轨迹 CSV</a> · <a href="data/final_summary.json">完整 JSON 摘要</a> · <a href="data/CORE500_COHORT.json">CORE500</a> · <a href="data/FINAL_FORMAL_COHORT.json">FULL cohort</a> · <a href="outputs/formal/">原始轨迹</a> · <a href="logs/final_tests.txt">最终测试日志</a>'
 page=f'''<!doctype html><html lang="zh"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>P9-A.5 · {n}条正式微泡轨迹审核</title><style>body{{font-family:system-ui,sans-serif;max-width:1380px;margin:36px auto;padding:0 18px;background:#f2f5f8;color:#203044;line-height:1.7}}h1{{font-size:30px}}h2{{font-size:20px}}.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px}}.card,section{{background:white;padding:20px;border:1px solid #dde5ec;border-radius:10px;margin-bottom:18px}}.card{{margin:0}}.card span{{display:block;color:#617387;font-size:14px}}.card strong{{font-size:26px}}img,video{{width:100%;height:auto;border-radius:5px}}a{{color:#1767a6}}.status{{font-family:monospace;font-size:13px;word-break:break-all}}</style></head><body><h1>P9-A.5 · 正式 {n:,} 条微泡轨迹</h1><p>{html.escape(why)}</p><div class="cards">{cardhtml}</div><section><h2>先看结论与比较</h2><p>CORE500：O1 / O2 / O3 = {core['O1']} / {core['O2']} / {core['O3']}；FULL {n}：{full['O1']} / {full['O2']} / {full['O3']}。</p><p>没有按出口或point basin选样。静止必须有接触/KKT支持，时间不够会单独标为删失。至少1颗只表示观察到了该出口；概率估计仍需看不确定区间。</p><p class="status">{html.escape(s['final_status'])}</p><p>{links}</p></section>{images}{videos}{dt_html}<section><h2>来源与复现</h2><p>{links}</p><p>旧回归 {s['baseline_tests']['tests']} 项；新永久测试 {s['new_tests']['tests']} 项；最终回归 {final_tests['tests']} 项。旧文件变化：{html.escape(str(s['old_files_changed']))}。</p><p>本页的结论来自实际正式manifest及完整轨迹库。轨迹为固定流场中的独立单泡计算，不包含MB–MB或RBC耦合。</p></section></body></html>'''
 (R/'OPEN_RESULTS.html').write_text(page)
 print('REVIEW_AND_HTML_READY',n)
if __name__=='__main__':main()
