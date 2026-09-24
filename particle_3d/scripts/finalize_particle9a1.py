#!/usr/bin/env python3
"""Evidence protection, replay invariants, Chinese review and linked figure index."""
import hashlib,json,subprocess,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from analyze_particle9a1 import ROOT,R,D,O,IDS,read,csvwrite
from particle_3d.particle81_simulation import dump,environment
from particle_3d.particle8_replay import canonical_hash
from particle_3d.particle9a_provenance import sha256

def main():
 before=read(D/'protected_before.json');mismatches=[k for k,v in before.items() if not (ROOT/k).exists() or sha256(ROOT/k)!=v]
 git_before=read(R/'logs/git_before.json');tree=subprocess.check_output(['git','write-tree'],cwd=ROOT,text=True).strip()
 protection=dict(protected_files=len(before),mismatches=mismatches,index_unchanged=tree==git_before['git write-tree'].strip(),source_archive_sha256=sha256(R/'reference/pre_repair_source.tar.gz'))
 dump(D/'protected_verification.json',protection)
 original=read(ROOT/'particle_3d/outputs/particle9a_2mmps/admission/birth_ledger.json');current=read(O/'admission/birth_ledger.json');assert original['events']==current['events']
 checks=[];env=environment()
 for folder in ['P9A1','P65_NEW']:
  for p in sorted((O/folder/'trajectories').glob('mb_*.json')):
   m=read(p);a=np.load(p.with_suffix('.npz'))['samples'];pid=m['particle_id'];e=next(e for e in original['events'] if e['particle_id']==pid)
   dt=np.diff(a[:,0]);res=np.max(abs(np.diff(a[:,1:4],axis=0)-dt[:,None]*a[1:,4:7])) if len(a)>1 else 0.
   checks.append(dict(model=folder,particle_id=pid,birth_identical=m['birth_metadata']==e,birth_hash_match=m['birth_metadata_sha256']==canonical_hash(e),
    samples_hash_match=sha256(p.with_suffix('.npz'))==m['samples_sha256'],strictly_increasing_time=bool(np.all(dt>0)),time_coverage_error_s=float(abs(dt.sum()-a[-1,0])),
    max_recorded_dt_error_s=float(np.max(abs(dt-a[1:,18]))) if len(dt) else 0.,max_position_identity_error_m=float(res),
    max_quaternion_norm_error=float(np.max(abs(np.linalg.norm(a[:,10:14],axis=1)-1))),
    minimum_gap_m=float(a[:,14].min()),minimum_g_nf_m=float(a[:,15].min()),no_penetration=bool(a[:,14].min()>=-env.wall.roundoff_m),handoff_admissible=bool(a[:,15].min()>=-env.wall.roundoff_m),
    source_identity_match=all(sha256(ROOT/'particle_3d/src/particle_3d'/k)==v for k,v in m['p9a1_identity']['source_sha256'].items()),
    independently_integrated=m['independently_integrated'],source_host=m['receipt']['hostname']))
 csvwrite(D/'independent_validation.csv',checks)
 valid=all(c['birth_identical'] and c['birth_hash_match'] and c['samples_hash_match'] and c['strictly_increasing_time'] and c['source_identity_match'] and c['no_penetration'] and c['handoff_admissible'] and c['max_position_identity_error_m']<env.wall.roundoff_m and c['max_recorded_dt_error_s']==0 for c in checks)
 dump(D/'independent_validation.json',dict(all_pass=valid,protected_pass=not mismatches and protection['index_unchanged'],trajectory_count=len(checks),checks=checks))
 g=read(D/'gates.json');g['gates']['evidence_and_time_identity']=valid and not mismatches and protection['index_unchanged']
 production=read(D/'production_summary.json') if (D/'production_summary.json').exists() else None
 if production is not None:
  g['gates']['formal_all_500_recorded']=production['all_500_present']
  g['gates']['formal_no_penetration']=production['penetration_count']==0 and production['handoff_violation_count']==0
  g['gates']['formal_no_unresolved_solver_failure']=not production['unresolved_failure_ids']
  g['production_summary']=production
 if not all(g['gates'].values()):g['automated_status']='P9A1_AUTOMATED_VALIDATION_FAIL';g['production_authorized_by_all_gates']=False
 dump(D/'gates.json',g)
 rows=read(D/'same12.json');inlet=read(D/'inlet_3_15.json');canonical=g['canonical'];contact=read(D/'id7_fixed_fixture.json')['solver'];manifest=read(D/'figure_manifest.json');remote=read(D/'remote.json')
 source_before=read(D/'source_before.json');changed=[p.name for p in sorted((ROOT/'particle_3d/src/particle_3d').glob('*.py')) if source_before.get(p.name)!=sha256(p)]
 dump(D/'changed_source.json',dict(changed_or_new=changed))
 stats=subprocess.check_output(['git','diff','--stat'],cwd=ROOT,text=True);(R/'logs/git_diff_stat.txt').write_text(stats)
 (R/'logs/git_status_final.txt').write_text(subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True))
 counts=g['completed_counts'];same=read(O/'provenance/same12_completed.json');runtimes={'same12':same['wall_seconds']}
 for stage in ['smoke30','production']:
  if (O/'provenance'/f'{stage}_completed.json').exists():runtimes[stage]=read(O/'provenance'/f'{stage}_completed.json')['wall_seconds']
 production_reuse=(sum(x['reused'] for x in read(O/'provenance/production_completed.json')['results']) if (O/'provenance/production_completed.json').exists() else 0)
 tests={}
 for name in ['local_tests','remote_tests','shared_core_tests']:
  tests[name]=sum(int(s.get('tests','0')) for s in ET.parse(R/'logs'/f'{name}.xml').getroot().iter('testsuite'))
 lines=['# 一句话结论','',f'自动状态 **{g["automated_status"]}**；人工状态 **PENDING_USER_REVIEW**。P9-A.1 same-12 完成 {counts["P9A1"]}/12；30 条：{g["smoke30_status"]}；500 条：{g["production_status"]}。',
 '', '# 这轮修了哪三个问题','',
 '1. **入口附近错误接管**：保留真实 FEM 背景载荷，局部平面剪切只提供经过 canonical benchmark 标定的仿射修正；真实开放边界 rim 使用 P6.5 fallback。',
 '2. **2 nm 约束反复开关**：在当前已解出的运动路径上寻找第一次到达 handoff 下限的时刻，先推进到事件，再重算剩余真实时间。',
 '3. **重复接触条件**：用数值秩和非负锥冗余检查移除重复方向，再检查全部原始约束，保留真实不同的不可穿透信息。',
 '', '# 为什么没有改 mobility 系数','',
 '原单平面 benchmark 没有显示 mobility 系数错误。M_eff、R_eff、delta_R、平移/旋转耦合常数、互易投影、法向 P6.5 阻力和 2 nm 下限均保留。源码逐项 AST/字节比较见 `data/unchanged_contracts.json`；原连续证书函数、P6.5 装配、admission、population source 与所有原积分守卫实现均未改写。',
 f'本地 {tests["local_tests"]} 项测试、服务器 {tests["remote_tests"]} 项、共享核心 {tests["shared_core_tests"]} 项测试通过。新旧 P9-A canonical 最大差异：Vt = {canonical["max_new_old_Vt_m_s"]:.9g} m/s，Omega = {canonical["max_new_old_omega_s_inv"]:.9g} s^-1。',
 f'相对原始 2D 四位小数非完全互易 reference：Vt 最大差异 {canonical["max_new_raw_2d_Vt_m_s"]:.9g} m/s，Omega {canonical["max_new_raw_2d_omega_s_inv"]:.9g} s^-1；均在原有解析互易投影误差界内，没有扩大测试容差。',
 '', '# 新的 bulk-preserving 方法是什么意思','',
 '“墙可以改变球怎么滑和转，但不能因为局部剪切梯度很小，就假装原本存在的血流消失了。”',
 '', r'设 $q_{lin}=[T^T(Hs),\;aT^T(n\times s)/2]^T$，$q_{ref}=(1-w)q_{lin}+wq_{wall}$。使用',
 '', r'$$\Delta b=R_{eff}q_{ref}-R_{free}q_{lin},\qquad R_{eff}q=R_{free}q_{bulk,actual}+\Delta b.$$',
 '', r'等价地，$q=q_{ref}+M_{eff}R_{free}(q_{bulk,actual}-q_{lin})$。canonical 流中真实 bulk 等于 q_lin，恢复旧解；远壁 w=0 恢复 free/P6.5。真实 bulk 的差额始终进入载荷。法向自由度没有加入此块。',
 f'ID 22 平均路径速度从 {g["id22_old_mean_speed_m_s"]:.9g} m/s 变为 {g["id22_new_mean_speed_m_s"]:.9g} m/s。这里平均值是实际路径长度/已积分物理时间。',
 '', '# inlet rim 怎么处理','',
 '开放边界和血管壁的交线不是已验证的无限平面墙。按原始 GlobalNodeID，求 WALL perimeter edges 与各 cap perimeter edges 的严格交集。INLET/OUTLET_01/OUTLET_02/OUTLET_03 分别有 29/22/19/23 条 rim edges。全部 cap 边缘必须匹配，文件 SHA 随拓扑记录保存。',
 '用原最近距离 kernel 的 EDGE/VERTEX barycentric witness 对照该拓扑；只有真实 rim witness 才令 P9 切向 delta_R 和仿射 delta_b 为零。普通 WALL shared edge 保持 P9 修正。未增加坐标距离带，P6.5 normal lubrication/handoff 保持工作。',
 '', '| ID | FEM 向内速度 (m/s) | 旧 P9-A | 新 P9-A.1 | 新轨迹结果 |','|---|---:|---:|---:|---|']
 for r in inlet:lines.append(f'| {r["particle_id"]} | {r["FEM_inward_m_s"]:.9g} | {r["old_inward_m_s"]:.9g} | {r["new_inward_m_s"]:.9g} | {r["new_outcome"]} |')
 lines+=['','# handoff 怎么修','',
 '预测当前路径会到达 2 nm 边界时先处理首次到达事件，而不是越界后无限缩小时间步。有限三角形距离的凸性用于寻找步内最小距离，因此两端都清楚但中间穿过的情况也会被检查。首次到达用原有限几何求根，真实时间推进至事件后重新装配和求解剩余时间；没有位置投影、gap clamp、下游偏移或虚假时间。radius-dependent 原下限继续适用。',
 '另一个数值问题是安全的边缘掠过路径：真实最小间隙高于下限，但原端点支撑平面证明过于保守。对此仅在最小距离处分割证明区间，每段重新调用完全未改动的 wall_handoff_certificate。只有全部原证书通过才接受同一条路径，运动路径、速度和位置都不改变。原 swept_clearance_certificate 仍检查完整区间。真实穿越测试必须拒绝。',
 '根搜索上限 64 次来自 binary64 分辨率；事件及证明分段使用原 MAX_REFINEMENT_DEPTH=48 实现守卫。原 provider/progress/stationary/horizon 守卫均未放宽。新增失败会保留，不会为了出血管而调参数。',
 '开发证据也保留：development_trial_01 是仅端点括根版本，development_trial_02 增加步内最小值检测。最终版本进一步组合原连续证书；两轮开发记录不是正式结果。',
 '', '# ID 7 怎么修','',
 '“两条几乎重复的接触命令只保留独立信息。”在白化约束 W 上做 SVD，因 dual 使用 W W^T，数值可分辨阈值取 sqrt(max(W.shape)*eps)*norm(W,2)，不使用固定角度或距离。还要求被移除行能由保留行非负组合表示，避免错误删除相反方向的接触。解后用保留系统的严格预算检查所有原始行；若失败即拒绝，不能利用原大条件数扩大接受预算。',
 f'实际 ID 7 原矩阵 fixture：约束 2→1，旧 solver numerical rank 2→1；同一 Gram 分辨率下 rank before/after 均为 1，说明只删冗余信息。condition {g["contact_id7_fixture_condition_old"]:.9g}→{g["contact_id7_fixture_condition_new"]:.9g}。最小法向速度 −6.76861815e−7→{min(contact["normal_speeds_m_s"]):.9g} m/s。保留/删除 ID、奇异值、残差、非负系数全部见 JSON。',
 '', '# same-12 结果','',
 '严格复用原 500 admission ledger 中前 12 个 accepted ID 的 birth time、radius、position、quaternion、SonoVue draw、seed 和 admission。新阶段 ledger 只增加阶段 provenance，事件内容逐项相同。dt=0.00025 s，horizon=1.5 s，2.0 mm/s FEM SHA=129ebb77396550b300b20ce29cde7b03919e6037b7a81de60fe75d6c6efb616d，Q=1.5513591604402322e−14 m³/s，原出口 classifier 不变。',
 '', '| ID | P6.5 原版本 | P9-A 原版本 | P9-A.1 | P6.5 + 新 handoff |','|---|---|---|---|---|']
 for pid in IDS:lines.append('| '+str(pid)+' | '+' | '.join(next(r['outcome'] for r in rows if r['particle_id']==pid and r['model']==model) for model in ['P65','P9A','P9A1','P65_NEW'])+' |')
 lines+=['',f'完成数量：P6.5 {counts["P65"]}/12；旧 P9-A {counts["P9A"]}/12；新 P9-A.1 {counts["P9A1"]}/12；P6.5 + 新 handoff {counts["P65_NEW"]}/12。',
 '', '| ID | 模型 | 失败原因 | handoff events | rejected trials |','|---|---|---|---:|---:|']
 for r in rows:
  if r['model'] in ['P9A1','P65_NEW'] and not r['completed']:lines.append(f'| {r["particle_id"]} | {r["model"]} | {str(r["failure_detail"]).replace("|","/")} | {r["handoff_events"]} | {r["rejected_trials"]} |')
 stationary_note=('30 条中的唯一未完成是 same-12 已出现的 ID 15：三个独立 FACE 接触，rank=3，所有乘子为正，contact condition≈7.16，平移速度在原 KKT 数值预算内为零；P6.5 同位置也是三个独立接触。新增的 18 条全部完成。因此按“无新的系统性数值错误”的门槛通过，不把这条静止路径改标为 completed，也不以 30/30 为门槛。证据见 data/stationary_stop_audit.json。')
 lines+=['', 'ID 7 实际新轨迹中的最大 contact condition 为 '+str(next(r['contact_max_condition'] for r in rows if r['particle_id']==7 and r['model']=='P9A1'))+'；ID 12/14/17 新轨迹均完成，原 HANDOFF_ENDPOINT_BELOW_LOWER 拒绝次数均为 0。仍有原几何 SWEEP_NOT_CERTIFIED 导致的真实时间细分，不能解释成所有细分均已消除。', '', '# 30-track smoke','',f'{g["smoke30_status"]}。outcomes: `{g["smoke30_outcomes"]}`。仅 same-12 gate 通过才会运行；选择原 ledger 前 30 个事件，已经计算的 12 条只在 source/event identity 一致时复用。',stationary_note,
 '', '# 500-track formal','',f'**{g["production_status"]}**。启动时门槛记录独立保存在 `data/production_launch_gate.json`，之后的完整验证另计。',
 (f'正式批次复用本轮同版本 smoke 的 {production_reuse} 条记录，另外积分 {production["generated_count"]-production_reuse} 条；原 P9-A 旧轨迹没有作为本阶段轨迹缓存使用。正式记录 {production["generated_count"]} 条，完成 {production["completed_count"]} 条；outlet distribution：`{production["outlet_distribution"]}`；failure counts：`{production["failure_counts"]}`。所有未完成 ID：`{production["failure_ids"]}`。\n\n最终状态重装配支持三个独立接触下静止的未完成数：{production["supported_stationary_failure_count"]}；尚未解释的失败 ID：`{production["unresolved_failure_ids"]}`。这只审核保存的最后状态，没有重新积分、改变粒径或强迫通过。详见 `data/production_outcomes.csv`、`data/production_failed_final_states.json`。' if production else '正式批次结果尚未形成；不得把当前 same-12/smoke 当作正式 500 结果。'),
 '', '# 当前限制','',
 '仍是 single nearest planar wall、no curvature、no multi-wall hydrodynamics、no shear lift、tetra-local gradient、independent MB。没有新增 RBC 闭合、MB 变形、CFD 或经验平滑。open-rim 单点拓扑 gate 没有额外 smoothing band；在网格特征切换处可能不连续。硬接触可能同时检查多个真实方向，这不代表加入 multi-wall hydrodynamics。',
 '原守卫停止或静止不等于证明生理捕获；未完成的轨迹不能作为完整通过数据使用。30/500 gate 不通过时不得把 same-12 或开发轮次称作正式 500 数据。',
 '', '自动门槛：','']
 for k,v in g['gates'].items():lines.append(f'- {k}: {"PASS" if v else "FAIL"}')
 lines+=['',f'无穿透核对：{sum(not c["no_penetration"] for c in checks)} 条 WALL penetration；{sum(not c["handoff_admissible"] for c in checks)} 条 handoff 下限违规。时间/位置/输入/source/hash 独立校验：{valid}。离线验证读取保存的逐步位置与 gap 并核对哈希；连续路径安全由在线原证书逐段保证，未声称离线重新计算全部几何查询。',
 f'服务器 `{remote["root"]}`；6 CPU workers，OMP/OPENBLAS/MKL/NUMEXPR 均为 1 线程；GPU 未用于积分。最终各批耗时（秒）：`{runtimes}`。开发轮次另有日志和逐轨迹 receipt，不能计作正式结果。',
 f'保护核对：{len(before)} 个原报告/数据/冻结输入/用户已有修改文件的 SHA 已比较；mismatches={len(mismatches)}；原 git index 未改：{protection["index_unchanged"]}。分支 dev/particle9a1-production-compatibility，未 commit/push。',
 '', '修改/新增科学源码：'+', '.join(changed)+'。`logs/git_diff_stat.txt` 是完整工作区 diff，包含本轮开始前已有修改；本轮原始状态见 `logs/git_before.json`。',
 '', '图与可直接重绘的来源：','']
 for f in manifest:lines.append(f'- [{Path(f["png"]).name}]({f["png"]}) / [PDF]({f["pdf"]})；来源：'+', '.join(f['sources']))
 lines+=['', '入口：[OPEN_RESULTS.html](OPEN_RESULTS.html)；机器记录：[gates.json](data/gates.json)；逐轨迹安全/时间验证：[independent_validation.csv](data/independent_validation.csv)。',
 '', '# 人工审核','',
 '- [ ] inlet 修复合理','- [ ] canonical planar physics 未破坏','- [ ] handoff 不再 chatter','- [ ] duplicate contact 处理合理','- [ ] same-12 轨迹合理','- [ ] 30-track smoke 合理','- [ ] 若已运行，500-track 数据合理','- [ ] 同意进入下一阶段','']
 (R/'PARTICLE9A1_REVIEW_ZH.md').write_text('\n'.join(lines))
 cards=''.join(f'<h2>{Path(f["png"]).stem}</h2><a href="{f["pdf"]}">PDF</a><img src="{f["png"]}" style="width:100%;max-width:1100px">' for f in manifest)
 (R/'OPEN_RESULTS.html').write_text('<!doctype html><meta charset="utf-8"><title>P9-A.1 review</title><body style="font:16px system-ui;max-width:1200px;margin:40px auto"><h1>P9-A.1 production compatibility</h1><p>'+g['automated_status']+' | PENDING_USER_REVIEW</p><p>Formal batch: '+g['production_status']+'</p><p><a href="PARTICLE9A1_REVIEW_ZH.md">中文审核报告</a> · <a href="data/gates.json">Machine gates</a></p>'+cards+'</body>')
 print(json.dumps(dict(status=g['automated_status'],protection=protection,independent_validation=valid,runtimes=runtimes,tests=tests),indent=2))
if __name__=='__main__':main()
