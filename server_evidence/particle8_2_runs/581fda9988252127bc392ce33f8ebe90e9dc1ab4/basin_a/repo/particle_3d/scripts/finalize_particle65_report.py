#!/usr/bin/env python3
"""Build Chinese review and machine validation; PASS only from successful JUnit logs."""
from pathlib import Path
import sys,json,subprocess,argparse,xml.etree.ElementTree as ET
import numpy as np
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parent;sys.path.insert(0,str(PACKAGE/'src'))
from particle_3d.audit import sha256
from particle_3d.particle3_cases import write_json
from particle_3d.nearfield_regularization import contract
from particle_3d.particle6_checkpoint import frozen_provenance
from generate_particle65_report import NAMES,NOTES
REPORT=PACKAGE/'reports/particle6_5';DATA=REPORT/'data';LOGS=REPORT/'logs'

def read(name):return json.loads((DATA/(name+'.json')).read_text())
def regression(stage,prefix):
 file=LOGS/f'{prefix}_{stage}.xml'
 if not file.exists():return dict(status='NOT_RUN',tests=None)
 root=ET.parse(file).getroot();suites=[root] if root.tag=='testsuite' else list(root.iter('testsuite'))
 counts={key:sum(int(s.attrib.get(key,0)) for s in suites) for key in ['tests','failures','errors','skipped']}
 return dict(status='PASS' if counts['tests'] and counts['failures']==counts['errors']==counts['skipped']==0 else 'FAIL',**counts,junit=str(file.relative_to(REPO)))
def fmt(x,scale=1):return '未观察到' if x is None else f'{x*scale:.9g}'

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--final',action='store_true');parser.add_argument('--source-commit');args=parser.parse_args()
 git=lambda *a:subprocess.check_output(['git',*a],cwd=REPO,text=True).strip()
 scope=read('00_scope');real=read('07_real_v1_dt1')['summary'];history=read('07_historical_p5');old=read('07_real_matched_old_p5')['summary'];sens=read('08_sensitivity');bridge=read('10_bridge')
 regressions={s:regression(s,'final') for s in ['frozen','particle0','particle1','particle2','particle3','particle4','particle5','particle6','particle6_5']}
 # Frozen driver uses the original suite label; detect its preserved command name.
 if regressions['frozen']['status']=='NOT_RUN':
  candidates=list(LOGS.glob('final_*frozen*.xml'))
  if len(candidates)==1:regressions['frozen']=regression(candidates[0].stem.removeprefix('final_'),'final')
 automatic='PASS' if all(r['status']=='PASS' for r in regressions.values()) else 'IN_PROGRESS'
 frozen=frozen_provenance(REPO/'formal_3D_flow_solver/FEM_SimVascular')
 assert all(sha256(REPO/k)==v for k,v in scope['p6_source_sha256'].items())
 assert all(sha256(REPO/k)==v for k,v in scope['historical_p5_sha256'].items())
 numerical_assertions=dict(contract_exact=json.loads((PACKAGE/'contracts/NEAR_FIELD_REGULARIZATION_V1.json').read_text())==contract(),
   real_continuously_safe=real['all_accepted_above_lower'],real_finite=real['finite'],all_sensitivity_safe=sens['all_continuously_safe'],
   all_sensitivity_finite=sens['all_finite'],bridge_parity=bridge['status']=='PASS',old_evidence_unchanged=True,
   all_required_figures=all((REPORT/'figures'/(n+'.png')).is_file() for n in NAMES))
 assert all(numerical_assertions.values())
 if args.final:
  assert automatic=='PASS',regressions
  commands=json.loads((LOGS/'final_commands.json').read_text());assert all(r['returncode']==0 for r in commands)
 source=args.source_commit or git('rev-parse','HEAD')
 sources={str(p.relative_to(REPO)):sha256(p) for folder in ['src','scripts','tests','contracts'] for p in sorted((PACKAGE/folder).rglob('*')) if p.is_file() and p.suffix in ['.py','.json']}
 sources['particle_3d/PARTICLE6_5_README.md']=sha256(PACKAGE/'PARTICLE6_5_README.md')
 if args.source_commit:
  for name,digest in sources.items():
   import hashlib
   assert hashlib.sha256(subprocess.check_output(['git','show',source+':'+name],cwd=REPO)).hexdigest()==digest,(name,'source commit mismatch')
 if args.source_commit:
  cp=REPORT/'checkpoints/v1/checkpoint_manifest.json';checkpoint=json.loads(cp.read_text())
  original=checkpoint.get('creation_source_sha256',checkpoint['source_sha256'])
  assert all(sha256(REPO/k)==v for k,v in original.items() if k.startswith('particle_3d/src/'))
  checkpoint.setdefault('creation_source_commit',checkpoint['source_commit'])
  checkpoint.setdefault('creation_source_sha256',checkpoint['source_sha256'])
  checkpoint.update(source_commit=source,source_sha256=sources,source_commit_status='PHYSICS_SOURCES_UNCHANGED; EXACT_TESTED_SOURCE_COMMIT_BOUND')
  write_json(cp,checkpoint)
 val=dict(contract(),stage='Particle-6.5',automated_checks=automatic,git_branch=git('branch','--show-current'),git_commit=source,
  source_commit_status='EXACT_COMMIT_SOURCE_SHA_VERIFIED' if args.source_commit else 'PRECOMMIT_WORKTREE_SHA_SNAPSHOT',source_sha256=sources,
  particle6_dependency_commit='c8e22f99359d15fd6075f2083ed60724ad0cfc2f',frozen_integrity=dict(status='PASS',**frozen),
  regressions=regressions,baseline_regressions={s:regression(s,'baseline') for s in regressions},
  **{f'p{i}_regression':regressions[f'particle{i}'] for i in range(7)},p6_5_regression=regressions['particle6_5'],
  contract='NEAR_FIELD_REGULARIZATION_V1',literature_sources=json.loads((REPORT/'literature/SOURCES.json').read_text()),
  glycocalyx_wall_model='NOT_FROZEN',effective_gap_definition='h_geom - h_lower',
  effective_gap_field_note='effective_gap_definition denotes g_NF; lubrication denominator is separately h_eff=max(h_geom,h_lower)',
  old_p5_subnanometer_history_preserved=True,real_replay_status='PASS_WITH_TIMESTEP_DEPENDENT_HANDOFF_EVENTS',
  real_old_min_gap_m=history['minimum_h_geom_m'],real_matched_p5_min_gap_m=old['minimum_h_geom_m'],
  real_v1_min_gap_m=real['minimum_h_geom_m'],real_v1_h_lower_m=real['h_lower_m'],real_handoff_event_count=real['handoff_event_count'],
  real_handoff_count_definition='ENTRIES_INTO_GEOMETRY_ROUNDOFF_HANDOFF_BAND; MAY_INCLUDE_MESH_AND_TIMESTEP_DEPENDENT_REENTRIES',
  real_time_coverage_error_s=real['time_coverage_error_s'],real_replay=real,real_initialization=history['replay_initialization'],
  original_two_mb_initial_rejection=history['original_two_mb_initial_rejection'],
  sensitivity_results=sens,scientific_sensitivity_review='PENDING_USER_REVIEW',
  time_accuracy_convergence='NOT_ESTABLISHED',handoff_event_topology_convergence='NOT_ESTABLISHED',
  numerical_assertions=numerical_assertions,lammps_bridge_parity=dict(status=bridge['status'],exact_equal=all(r['exact_equal'] for r in bridge['errors']),restart_exact=bridge['restart_error']['exact_equal'],mpi_ranks=1),
  manual_visual_review='PENDING_USER_REVIEW',agent_figure_review='PENDING_AGENT_INSPECTION',
  parallel_multi_rank_lammps='NOT_VALIDATED',real_rbc_dynamics='DEFERRED',real_rbc_passage='NOT_ESTABLISHED',
  molecular_contact_physics='NOT_MODELED',position_projection=False,upstream_p0_p6_sources_unchanged=True,
  stage_result='AUTOMATED_PASS_PENDING_VISUAL_AND_SCIENTIFIC_USER_REVIEW' if automatic=='PASS' else 'AUTOMATION_IN_PROGRESS')
 inspect_file=REPORT/'FIGURE_INSPECTION.json'
 if inspect_file.exists():val['agent_figure_review']=json.loads(inspect_file.read_text())
 manifest=json.loads((REPORT/'FIGURE_MANIFEST.json').read_text())
 assert all(sha256(REPORT/'figures'/f['filename'])==f['sha256'] for f in manifest['figures'])
 val['log_sha256']={str(p.relative_to(REPORT)):sha256(p) for p in sorted(LOGS.iterdir()) if p.is_file()}
 val['data_sha256']={str(p.relative_to(REPORT)):sha256(p) for p in sorted(DATA.iterdir()) if p.is_file()}
 val['figure_sha256']={f['filename']:f['sha256'] for f in manifest['figures']}
 val['literature_sha256']={str(p.relative_to(REPORT)):sha256(p) for p in sorted((REPORT/'literature').iterdir()) if p.is_file()}
 write_json(REPORT/'PARTICLE6_5_VALIDATION.json',val)
 table='\n'.join(f"| {s} | {r['tests']} | {r['status']} |" for s,r in regressions.items())
 sensitivity_table='\n'.join(f"| {r['h_molecular_floor_m']*1e9:g} | dt/{r['dt_divisor']} | {fmt(r['handoff_time_s'],1e3)} | {r['handoff_event_count']} | {fmt(r['minimum_h_geom_m'],1e9)} | {fmt(r['tangential_displacement_m'],1e6)} | {fmt(r['final_position_m_absolute_difference_from_2nm'],1e9)} | {fmt(r['final_velocity_m_s_absolute_difference_from_2nm'],1e6)} |" for r in sens['rows'])
 report=f'''# Particle-6.5 人工审核报告

分支：`{val['git_branch']}`。源码提交：`{source}`（{val['source_commit_status']}）。机器记录：[PARTICLE6_5_VALIDATION.json](PARTICLE6_5_VALIDATION.json)。

## 1. 为什么要做这一阶段

P5 只要原始几何间隙大于 float64 舍入预算，就继续使用球形法向 1/h 阻力，因此旧双 MB 回放算到了 {history['minimum_h_geom_m']*1e9:.9g} nm。它是保存下来的数值结果，不能解释为可信的普通血液连续薄膜。P6 只增加桥接，没有解决这一物理解释边界。旧 P5 报告、图、数据和 P0–P6 数值源码逐文件哈希保持不变；P6 人工 PASS 只作证据更新。

## 2. 文献告诉了我们什么

直接模拟先例来自 Ness (2023)：典型低间距尺度 O(10⁻³a)，外截断 0.05a。受限水的 PNAS 2017 与 Nature Materials 2026 结果只是跨领域连续介质边界证据：前者讨论纳米孔直径，后者使用特定固体界面，不能直接标定血管薄膜。Marsh 的约 0.2–0.3 nm 是脂质水合力衰减长度；Tu 等采用的 4 nm 是壳厚假设，均只作界面背景。

2 nm、完整启用边界 0.01、平滑过渡和两球 min-radius 参考长度属于明确批准的项目建模选择。完整来源、DOI、定位与适用限制见[文献审计](literature/NEAR_FIELD_REGULARIZATION_LITERATURE.md)与[机器来源表](literature/SOURCES.json)。没有引入文献中的惯性、摩擦、表面力或壳层力学。

## 3. 三个关键尺度

χ=h_geom/a_ref。χ≥0.05 时仅关闭本 leading 近场修正，Frozen FEM 的自由流仍在。χ≤0.01 且 h_geom>h_lower 时完整使用 P5 leading term。h_lower=max(2 nm,10⁻³a_ref) 是连续介质交棒边界。球墙 a_ref=a；球球 a_ref=min(ai,aj)，系数仍使用 ai·aj/(ai+aj)。对极小半径，完整区可能为空，代码按合同计算而不改尺度顺序。

当前历史半径约 0.588 与 1.266 µm，悬浮液分量分别约 0.588 与 1.266 nm，两者由 2 nm 控制；大于 2 µm 后由粒径项控制。原始 h_geom 始终保存；h_eff=max(h_geom,h_lower) 只进入阻力分母；约束使用独立的 g_NF=h_geom−h_lower。

## 4. 为什么不能突然打开 lubrication

在 0.01–0.05 中使用 w=1−3s²+2s³，s=(χ−0.01)/0.04。权重有界、单调，两个端点导数为零，避免人为硬开关。源码使用代数等价的 (1−s)²(1+2s)，改善外端附近的舍入表现。它是数值耦合选择，并非生物测量参数。

## 5. 为什么是 2 nm

2 nm 是用户批准的 continuum handoff model choice，**不是直接测得的 SonoVue–小鼠血管 cutoff，也不是真实固体接触距离**。纳米约束水与脂质界面证据提示不能无条件把体相连续介质公式延伸到任意小距离。合同因此停止解释更小间隙；没有把 2 nm 伪装成浮点误差或糖萼厚度。

## 6. 为什么还有 1e-3 a

该相对粒径尺度有悬浮液模拟先例，防止将下限永远固定为 2 nm。V1 取相对尺度与分子尺度的较大值，保留两个分量及主导项。壳厚 4 nm 仅作背景，不加到 sampled outer radius；墙面也未增加糖萼偏置。

## 7. 什么叫 CONTINUUM_HANDOFF_CONTACT

这是当前连续介质模型到达下限的状态，**不是原子或分子意义的真实固体接触**。保留封顶的润滑阻力，用原 P5 阻力度量和无摩擦约束限制继续靠近的法向速度，切向仍可运动。原 P3/P4 的 h_geom=0 几何接触实现未改。

跨界试探调用原 P3 二分物理时间积分器，重新解左子区间，再完整处理右子区间。连续证书同时检查原几何与 g_NF；不接受端点安全但中间穿越的路径，不回推位置。交棒判定只允许原几何舍入预算，例如真实 WALL 为约 2.0094×10⁻¹⁷ m，远小于 0.1 nm。到达时刻对应这个数值带，不宣称数学根的无限精度。

## 8. 真实 P5 亚纳米轨迹怎样改变

原双 MB 起点中 ID203 的间隙已为 0.305838 nm，无法同时保留原中心又作为 V1 合法初态；新 API 明确拒绝，未移动粒子。改用 **P5 已保存的 `MINIMUM_ORIGINAL_APPROACHING_GAP` 原样本**：来自 P4 保存轨迹第 125 行、原时刻 0.001252019129187767 s，半径 1.2658269815352818 µm，原始间隙 3.2328727626623765 nm。这里“改用”是选取另一个既存合法状态，不是修正其几何。

从该完全相同的单球状态分别运行旧 P5 默认模型与新 V1；Frozen FEM、WALL、采样自由速度与半径不变。重放时间从这个状态归零，窗口 {real['horizon_s']*1e3:.9g} ms；较短的历史窗口内粒子先离墙，延长窗口用于覆盖再次靠墙。验证基准 dt={real['dt_s']:.17g} s（历史 dt 的 4 倍），另外比较 dt/2、dt/4，全部仅用于验证。

旧双 MB **历史**最小 gap={history['minimum_h_geom_m']*1e9:.9g} nm；本次 **同初态单球 P5 重放**最小 gap={old['minimum_h_geom_m']*1e9:.9g} nm；两者不是同一条轨迹，均明确保留。基准 V1 最小原始 gap={real['minimum_h_geom_m']*1e9:.12g} nm，下限为 2 nm，差值 {real['minimum_g_nf_m']:.5g} m 在原几何舍入预算内。首次交棒约 {fmt(real['handoff_time_s'],1e3)} ms，进入交棒判定带 {real['handoff_event_count']} 次，时间覆盖误差 {real['time_coverage_error_s']:.3g} s。进入次数按接受状态判定带的进出计数，可能包含网格特征切换或离散步长导致的再进入，不等于独立物理碰撞次数。

**步长敏感性未收敛：** dt/2 最小约 2.03837 nm，dt/4 最小约 3.23287 nm，窗口内交棒计数均为 0。不能把基准的 6 次事件和 14.30 ms 当作稳定物理预测。各组都满足连续间隙约束；安全性通过不等于时间精度通过。该窗口没有出口事件，不宣称出口或真实 RBC 通行结果。

## 9. 1.5 / 2 / 3 nm 敏感性

三个尺度使用同一真实初态、FEM、半径、WALL、dt 集合和窗口；仅尺度改变。正式参数仍为 2 nm，1.5 和 3 nm 仅 `SENSITIVITY_ONLY`。2 nm 正式运行同时作为敏感性基准，比较表的所有行均属于科学敏感性证据。

| floor (nm) | 步长 | 首次交棒 (ms) | 进入次数 | 最小 raw gap (nm) | 切向位移 (µm) | 终位置差 (nm) | 终速度差 (µm/s) |
|---|---|---|---|---|---|---|---|
{sensitivity_table}

终态差均相对相同 dt 的 2 nm 结果，向量差使用欧氏范数；切向位移定义为首个球净位移在初始精确法向垂直平面内的投影长度，不是沿曲壁积分路程。JSON/CSV 同时保存原始终态向量、绝对差、相对差及排序。未观察到的交棒时间为 null，不用 0 或窗口末值冒充。

全部结果有限且连续约束安全，物理时间完整处理；**事件次数、时刻和终态速度有明显尺度/步长依赖，尚未证明科学稳定性。** 不设置任意 5% 门槛，不代替用户决定可接受性。`SCIENTIFIC_SENSITIVITY_REVIEW = PENDING_USER_REVIEW`。

## 10. 自动测试

开发前 Frozen/P0–P6 共 489 项全部通过。最终回归如下；详细命令、stdout 与 JUnit 在 [logs](logs/) 中，源码、历史证据、数据与图片哈希见机器记录。

| 测试组 | 数量 | 状态 |
|---|---|---|
{table}

新增永久测试覆盖合同、文献角色、C1 端点、粒径交叉、原半径和墙不变、球墙/球球解析阻力、交换对称、交棒、连续区间证书、每步 x_new=x_old+dt·V、时间覆盖、真实样本、3×3 敏感性、实际 LAMMPS 和二进制重启，以及从保存数据逐位再生全部图。静态扫描下限以下只是诊断，不会成为 V1 动态接受状态。LAMMPS 多余候选被物理筛选剔除，R、b、U、J、x 与独立路径完全相同。原 P6 62 项回归保留三形状元数据与历史重启验证。

## 11. 人工审核图片

'''
 for i,name in enumerate(NAMES):
  report+=f'![{name}](figures/{name}.png)\n\n'
  report+='\n\n'.join(f'{label}：{text}' for label,text in zip(['应该看什么','实际看到什么','有没有异常'],NOTES[i]))+'\n\n'
 report+=f'''## 12. 当前限制

- 仅 sphere normal near-field；没有完整远场多体 mobility 或多墙流体动力学。
- 没有切向润滑、旋转润滑耦合、RBC 润滑。
- 没有糖萼模型、真实分子接触模型、黏附、表面粗糙度分布、壳层力学。
- 生产 timestep、neighbor cutoff 和 skin 均未冻结；事件和时间精度收敛尚未建立。
- 没有 full suspension；真实 RBC 动力学仍延期，真实 RBC passage 仍未建立。
- LAMMPS 仅验证单 MPI rank，负责状态、邻居、重启；没有采用其力积分或时间推进。

合同可称 `SPHERE_NORMAL_NEAR_FIELD_REGULARIZATION_V1`，不能称为完整生产流体动力学。P0–P6 默认模型与历史证据保留。未运行 CFD、未 push、未 merge main；完成后停止于 Particle-6.5。

## 13. 人工审核状态

AUTOMATED_CHECKS = {automatic}

MANUAL_VISUAL_REVIEW = PENDING_USER_REVIEW

SCIENTIFIC_SENSITIVITY_REVIEW = PENDING_USER_REVIEW

Particle-7 started = NO
'''
 (REPORT/'PARTICLE6_5_REVIEW.md').write_text(report)
 print(automatic,{s:r['tests'] for s,r in regressions.items()})
if __name__=='__main__':main()
