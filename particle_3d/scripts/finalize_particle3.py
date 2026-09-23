#!/usr/bin/env python3
"""Persist full regression evidence and a Chinese human-review report; no next stage."""
from pathlib import Path
import argparse,json,subprocess,sys,xml.etree.ElementTree as ET
from collections import Counter
import numpy as np
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parent
sys.path.insert(0,str(PACKAGE/'src'))
from particle_3d.audit import sha256,read_frozen
from particle_3d.particle3_cases import write_json,P3_BRANCH,P2_COMMIT
from particle_3d.particle2_audit import check_dependencies
REPORT=PACKAGE/'reports/particle3';DATA=REPORT/'data';LOGS=REPORT/'logs'
FEM=REPO/'formal_3D_flow_solver/FEM_SimVascular'
HANDOFF=Path('/home/lzy/projects/ulm-3d-vascular-model-generation/formal_3D_flow_solver/FEM_SimVascular')
FROZEN='c83dccea0f1fe5cd9fd3f5939e2eba883aeba9c2'
P0='7cfe5141382600e28582f05bff712a6f09c38a39';P1='6e59c605cb8251b9fa2e80d2dbed0cfa6daaf03c'


def read(name):return json.loads((DATA/name).read_text())
def git(*args,cwd=REPO):return subprocess.check_output(['git',*args],cwd=cwd,text=True).strip()


def source_hashes():
    paths=[]
    for directory in ['src','tests','scripts','contracts']:
        paths.extend(p for p in (PACKAGE/directory).rglob('*') if p.is_file() and p.suffix in ['.py','.json','.md','.toml'])
    paths.extend([PACKAGE/'PARTICLE3_README.md',PACKAGE/'pyproject.toml'])
    return {str(p.relative_to(REPO)):sha256(p) for p in sorted(paths)}


def run(name,command,cwd):
    path=LOGS/(name+'.log')
    with path.open('w') as stream:
        stream.write('$ '+' '.join(map(str,command))+'\n');stream.flush()
        result=subprocess.run(command,cwd=cwd,stdout=stream,stderr=subprocess.STDOUT,text=True)
    return dict(name=name,command=list(map(str,command)),cwd=str(cwd),returncode=result.returncode,log=str(path.relative_to(REPO)))


def junit(path):
    root=ET.parse(path).getroot();suites=[root] if root.tag=='testsuite' else list(root.findall('testsuite'))
    result={k:sum(int(s.attrib.get(k,0)) for s in suites) for k in ['tests','errors','failures','skipped']}
    result['passed']=result['tests']-result['errors']-result['failures']-result['skipped']
    result['status']='PASS' if result['errors']==result['failures']==0 else 'FAIL'
    return result


def collect():
    p0,p1=check_dependencies(REPO)
    p2=json.loads((PACKAGE/'reports/particle2/PARTICLE2_VALIDATION.json').read_text())
    assert p2['git_commit']==P2_COMMIT and p2['manual_visual_review']=='PASS'
    for group in ['source_sha256','data_sha256','figure_sha256']:
        for path,digest in p2[group].items():assert sha256(REPO/path)==digest,path
    scope=read('00_particle3_scope_and_wall_contract.json');sphere=read('02_sphere_gap.json');ellipsoid=read('03_ellipsoid_gap.json');caps=read('03_capsule_gap.json')
    contact=read('04_contact_velocity.json');physical=read('05_time_summary.json');population=read('07_feasibility_map.json');examples=read('06_deformation_examples.json');controlled=read('09_controlled_summary.json')
    comparison=read('12_timestep_comparison.json');cases=comparison['cases'];assert len(cases)==18
    states={p.name:json.loads(p.read_text()) for p in sorted(DATA.glob('*_real_*_states.json'))};assert len(states)==18
    cap_rows=[r for rows in states.values() for r in rows if r['shape_mode']=='CAPILLARY_DEFORMED']
    v_errors=[r['volume_relative_error'] for r in population+examples if 'volume_relative_error' in r]+[abs(r['capsule_volume_m3']/r['original_volume_m3']-1) for r in cap_rows]
    a_errors=[r['area_violation_m2'] for r in population+examples if 'area_violation_m2' in r]+[max(0.,r['capsule_area_m2']-r['area_budget_m2']) for r in cap_rows]
    normals=[r['normal_constraint_error_m_s'] for r in contact]+[r['max_normal_error_m_s'] for r in controlled]+[r['normal_constraint_error_m_s'] for rows in states.values() for r in rows]
    normal_audit=read('04_real_contact_velocity_audit.json')
    normals.append(normal_audit['max_actual_normal_residual_m_s'])
    tangent=max([r['tangential_velocity_error_m_s'] for r in contact]+[r['max_tangential_error_m_s'] for r in controlled])
    counts=Counter(r['status'] for r in population);wall=scope['wall']
    result=dict(stage='Particle-3',automated_checks='PENDING_FINAL_CHECKS',git_branch=git('branch','--show-current'),git_commit=git('rev-parse','HEAD'),
        source_commit_status='PENDING_SOURCE_COMMIT',particle0_dependency_commit=P0,particle1_dependency_commit=P1,particle2_dependency_commit=P2_COMMIT,
        particle2_acceptance_evidence_commit='705dd0c',frozen_base_commit=FROZEN,
        wall_path=wall['path'],wall_sha256=wall['sha256'],wall_triangle_count=wall['triangle_count'],wall_boundary_ids=wall['boundary_ids'],
        only_wall_is_solid=True,wall_rigid=True,wall_stationary=True,physical_coating_thickness_m=0.,
        wall_gap_engine=scope['wall_gap_engine'],sphere_gap_tests='PENDING_FINAL_CHECKS',ellipsoid_gap_tests='PENDING_FINAL_CHECKS',capsule_gap_tests='PENDING_FINAL_CHECKS',
        max_synthetic_sphere_gap_error_m=max(r['error_m'] for r in sphere),max_synthetic_ellipsoid_gap_error_m=max(r['error_m'] for r in ellipsoid),max_synthetic_capsule_gap_error_m=max(r['error_m'] for r in caps),
        contact_model='HARD_FRICTIONLESS_KINEMATIC_V0',normal_velocity_constraint_error=max(normals),tangential_velocity_error=tangent,
        real_contact_velocity_audit=normal_audit,
        physical_time_refinement=dict(maximum_depth=48,guard_role='IMPLEMENTATION_SAFETY_ONLY',failed_step_consumes_time=False,center_projection=False,
            capsule_continuity='QUASISTATIC_MAX_FEASIBLE_RADIUS; continuous volume/area-valid auxiliary witness proves nonempty feasible set; accepted radius remains maximal'),
        no_tunneling_test='PENDING_FINAL_CHECKS',time_coverage_error=max([r['time_coverage_error'] for r in physical]+[r['time_coverage_error_s'] for r in cases+controlled]),
        rbc_deformation_model='REDUCED_ORDER_CAPSULE_V0',rbc_volume_conservation_error=max(v_errors),rbc_volume_conservation_error_unit='RELATIVE',
        surface_area_budget_role='MODEL_DERIVED_FROM_PARTICLE2_OBLATE_GEOMETRY',max_area_budget_violation=max(a_errors),area_violation_unit='m^2',
        deformation_validation_count=len(population),free_oblate_count=counts['FREE_OBLATE'],capillary_deformed_count=counts['CAPILLARY_DEFORMED'],deformation_infeasible_count=counts['DEFORMATION_SURROGATE_INFEASIBLE'],
        additional_deformation_examples_count=len(examples),additional_deformation_example_counts=dict(Counter(r['status'] for r in examples)),
        distribution_modified=False,distribution_array_sha256=read('07_selection.json')['population_array_sha256'],
        numerical_roundoff_rule=scope['roundoff_rule'],numerical_wall_tolerance_expanded=False,
        deformation_radius_roundoff_rule='512*float64_eps*equal_volume_sphere_radius; unchanged after search-floor fix',
        real_wall_tests=dict(status='PENDING_FINAL_CHECKS',normal_samples=18,static_samples=15,controlled_cases=len(controlled),regions=5,
            controlled_scope='ACTUAL_WALL_FINITE_TRIANGLE_PATCHES; full-WALL static gaps retained separately; real FEM uses entire WALL'),
        real_fem_mb_min_gap_m=min(s['minimum_gap_m'] for s in cases if s['particle_type']=='MB'),
        real_fem_mb_results=[s for s in cases if s['particle_type']=='MB'],real_fem_rbc_results=[s for s in cases if s['particle_type']=='RBC'],
        real_fem_deformation_events=[dict(case=s['case'],events=s['shape_mode_events'],termination=s['status']) for s in cases if s['particle_type']=='RBC'],
        real_fem_upstream_orientation_comparison=comparison['upstream_orientation'],real_orientation_timestep_convergence='NOT_ESTABLISHED',
        deformed_orientation_status='JEFFERY_ORIENTATION_NOT_INTERPRETED_WHILE_DEFORMED',
        wall_lubrication='DEFERRED_PARTICLE5',near_field_resistance='DEFERRED_PARTICLE5',transit_time_penalty='DEFERRED_PARTICLE5',rbc_mb_lateral_displacement='DEFERRED_PARTICLE4_5',
        adhesion='OFF',particle_particle_contact=False,production_particle_timestep_frozen=False,production_wall_model_status='V0_VALIDATION_ONLY_PENDING_USER_REVIEW',
        manual_visual_review='PENDING_USER_REVIEW',particle4_started=False,no_cfd_executed=True,
        trajectory_column_roles=dict(free_velocity_m_s='current accepted position free velocity',corrected_velocity_m_s='velocity applied over the interval ending at this row; initial row stores free velocity',normal_constraint_error_m_s='constraint at preceding interval start'),
        source_sha256=source_hashes(),data_sha256={str(p.relative_to(REPO)):sha256(p) for p in sorted(DATA.iterdir()) if p.is_file()},
        figure_sha256={str(p.relative_to(REPO)):sha256(p) for p in sorted((REPORT/'figures').glob('*.png'))})
    return result


def write_review(v):
    status=v['automated_checks'];tests=v.get('test_results',{})
    table='\n'.join(f"| {name} | {r['status']}：{r['passed']} passed | 0 failures / errors | logs/{name}_pytest.log |" for name,r in tests.items())
    rbc='\n'.join(f"| {s['geometry_index']} / {s['r']:.3f} | {s['validation_dt_s']*1e6:.6f} | {s['status']} | {s['rows']} | {s['last_accepted_time_s']:.9g} | {'无接受状态' if s['minimum_gap_m'] is None else format(s['minimum_gap_m'],'.5e')} |" for s in v['real_fem_rbc_results'])
    figures='\n\n'.join((REPORT/f'{i:02d}_step_notes.md').read_text() for i in range(13))
    text=f'''# Particle-3 人工审核报告

分支：`{v['git_branch']}`。被测源码提交：`{v['git_commit']}`（{v['source_commit_status']}）。
机器记录：[PARTICLE3_VALIDATION.json](PARTICLE3_VALIDATION.json)。复现与算法说明：[PARTICLE3_README.md](../../PARTICLE3_README.md)。

## 1. 这一阶段做了什么

现在第一次真正检查完整粒子有没有碰墙。球形 MB、按当前姿态放置的 RBC 椭球，以及窄管中的定体积胶囊代理，都与原始有限 WALL 三角形比较。
只读正式 WALL，共 {v['wall_triangle_count']} 个三角形；WALL 刚性、静止，物理涂层厚度为零。入口和三个出口不是实体壁。
WALL SHA256：`{v['wall_sha256']}`。P0/P1/P2 实现、科学参数、分布和 Frozen FEM 均保持原样。
用户对 P2 的聊天验收已先单独提交 `705dd0c`，只改 P2 两份审核证据；P2 分布仍为 PROVISIONAL_PASS，真实姿态步长收敛仍未建立。

## 2. wall gap 是什么

中心没有出血管，不等于整个粒子没碰墙。gap 是完整粒子表面到有限三角形的欧氏间隙：正值分离，零为接触，负值表示相交深度。
对每个三角形计算面、边、顶点的有限几何，再对 WALL 取最小值；中心在不在腔内仍由 P0 独立判断。
负值是该单个三角形的最小平移脱离深度取负，不声称等于同时脱离所有重叠壁面的整体最短位移。
数值舍入界固定为 512×float64 eps×实际坐标/几何尺度；真实 WALL 的量级约 2.0094e-17 m。
这是浮点误差预算，不是额外的物理间隙；开发中没有扩大它。

## 3. MB 怎么算

使用冻结 SonoVue sampler 给出的真实球半径；正式重放半径为 5.88015722765549e-7 m，起点完全沿用 P1。
球心到有限三角形最近点的距离减去半径，得到球面间隙，面、边和顶点均单独验证。
三个真实步长均记录每个接受状态。MB 最小间隙为 {v['real_fem_mb_min_gap_m']:.9e} m。
三个 MB 重放终止事件：{', '.join(s['status'] for s in v['real_fem_mb_results'])}；接触次数：{[s['contact_count'] for s in v['real_fem_mb_results']]}。
没有为展示碰墙而挪动入口起点。

## 4. RBC 怎么算

FREE_OBLATE 使用 P2 原始 a、b、c 和当前 quaternion。旋转会改变支撑高度，因此不能用中心距离减去一个固定 RBC 半径。
算法按椭球有限面/边/顶点的全部平稳特征求距离；独立凸优化的原始/对偶夹界验证外部距离，独立椭圆边界积分/优化及平面解析式验证其他分支。
真实重放使用原 5%、25%、50%、75%、95% r 的五个样本。所有轨迹始终查询整幅 WALL。

| 几何序号 / r | dt (µs) | 终止事件 | 接受状态数 | 最后接受时间 (s) | 最小 gap (m) |
|---|---:|---|---:|---:|---:|
{rbc}

无接受状态的初始不可行案例用 null gap 表示，不虚构零间隙或轨迹。它们的初始尺寸、面积、体积和冲突原因仍完整保存。
对于有变形的案例，只能解释为：当前 rigid surrogate 与 WALL 几何冲突，V0 reduced-order model 找到了满足 volume / area-budget / wall-clearance 的 capsule surrogate。
这不表示真实 RBC 在那里一定如此变形。

## 5. 硬接触怎么处理

不反弹、不摩擦、不穿墙，切向照常走。接触点向墙运动时，只调整中心法向速度来抵消自由平移与转动造成的入墙分量，Omega 不变。
球和胶囊采用平移约束；没有接触或正在离墙时不干预。多个同时接触点通过最小必要法向速度修正共同处理。
最大法向残差 {v['normal_velocity_constraint_error']:.6e} m/s；单法向受控试验最大切向变化 {v['tangential_velocity_error']:.6e} m/s。
真实残差逐区间用实际速度和所有实际接触法向独立复算。原有速度求解器舍入界为 1024×float64 eps×速度尺度；原有近重合法向合并阈值为 256×float64 eps。
合并法向与原始法向的微小差异按“法向差的长度×实际速度大小”单独传播到残差预算中，另加转动项差异。最大实际残差与完整预算之比为 {v['real_contact_velocity_audit']['maximum_residual_to_bound_ratio']:.9g}。
共有 {v['real_contact_velocity_audit']['intervals_above_qp_only_bound']} 个区间只看速度求解器预算时会超界；它们均受原有法向合并误差的传播界覆盖。首次漏计该项的审计失败记录已保存；任何 WALL、速度和法向合并阈值均未改变。
这是运动学约束，不计算真实接触力。
CSV 的自由速度是当前行位置的自由速度；corrected_velocity 是到达当前行所用的区间速度，须与前一行自由速度配对解释。
状态文件中的 normal_constraint_error 是求解器诊断；独立复算的完整残差另存 04_real_contact_residuals.csv，最终最大误差采用完整复算结果。

## 6. 为什么要真实时间细分

不能时间过去了但粒子偷偷留在原地。失败区间分成真实的左右两半，右半从已完成左半的状态继续；失败本身不推进时间。
薄墙试验中大步的两个端点都可能分离，仍必须检查中间是否穿墙。三种请求步长均在 0.375 s 接触，随后继续沿墙运动，完整覆盖到 1 s。
所有合成、真实 patch 和 FEM 案例最大时间覆盖误差为 {v['time_coverage_error']:.6e} s。
最大细分深度 48 和时间 ULP 下限都是实现保护，不是生产步长；保护触发时保存明确错误，不接受失败状态。
胶囊轴和半径变化使用重新计算的连续支撑路径。
在最大可行半径随位置变化时，可用同体积、同面积预算的辅助连续胶囊证明每一时刻存在可行形状；接受状态仍取全区间搜索的最大可行半径，辅助形状不替代 RBC，也不缩小体积。

## 7. 为什么窄管 RBC 可以变形代理

这里不用完整膜求解，而使用低成本胶囊几何。先检查实际原椭球 WALL 冲突，再沿局部自由流方向寻找胶囊，不能仅凭一个手写管径触发。
选择不添加迟滞：每个接受状态先尝试原椭球，允许时立即使用 FREE_OBLATE；否则尝试 CAPILLARY_DEFORMED。
在 128 个固定分布样本 × 9 人工管径的 {v['deformation_validation_count']} 个组合中，FREE / DEFORMED / INFEASIBLE 为
{v['free_oblate_count']} / {v['capillary_deformed_count']} / {v['deformation_infeasible_count']}。
另外保留五个分位样本加两个合法极端样本的 63 个例子。管径、初始间隙比例与指定控制速度均为 VALIDATION_ONLY，不代表真实小鼠血管分布。
胶囊仍不可行时明确停止，不能把代理不可行直接叫作生理堵塞。

## 8. capsule 怎样保证不过度变形

原 RBC volume 必须保持，所需光滑面积不能超过由原 oblate 推导的模型预算；没有修改 D/V 分布、缩小体积或增加面积。
L 是两半球球心间的轴段长度，满足 L=V0/(πR²)−4R/3；半径范围由定体积和面积约束共同求得。
在整个允许区间寻找最大的可行半径，不能假设缩细总有利，因为胶囊同时变长。搜索使用包含关系与保守距离变化界排除区间。
最大相对体积误差 {v['rbc_volume_conservation_error']:.6e}；最大面积超预算 {v['max_area_budget_violation']:.6e} m²。
A_budget 是 MODEL_DERIVED_FROM_PARTICLE2_OBLATE_GEOMETRY，不是直接测到的真实 membrane area。
小于预算的剩余面积只代表未展开部分，本模型不计算膜褶皱。

## 9. 自动测试

| 检查 | 结果 | 最大误差 / 条件 | 测试文件或完整日志 |
|---|---|---|---|
{table}
| 球有限面/边/顶点 | {status} | {v['max_synthetic_sphere_gap_error_m']:.3e} m | test_convex_triangle_geometry.py |
| 椭球姿态、穿透、有限边/顶点 | {status} | 平面解析 {v['max_synthetic_ellipsoid_gap_error_m']:.3e} m | test_support_feature_consistency.py、test_convex_triangle_geometry.py |
| 胶囊解析几何 | {status} | {v['max_synthetic_capsule_gap_error_m']:.3e} m | test_capsule_area_budget.py、test_convex_triangle_geometry.py |
| 接触法向 / 切向 | {status} | {v['normal_velocity_constraint_error']:.3e} / {v['tangential_velocity_error']:.3e} m/s | test_hard_contact.py、test_real_wall_controlled.py |
| 无穿墙 / 真实时间覆盖 | {status} | {v['time_coverage_error']:.3e} s | test_physical_time.py、test_real_fem_evidence.py |
| 体积 / 面积 / 最小变形 / 明确不可行 | {status} | 相对体积 {v['rbc_volume_conservation_error']:.3e}，面积 {v['max_area_budget_violation']:.3e} m² | test_deformation_selection.py、test_capsule_area_budget.py |
| 输入、方向、范围、原分布不变 | {status} | 18 个所属四面体方向检查及 P2 全部历史 SHA | test_input_scope_and_normals.py |
| 13 图、源数据、重作图、报告 | {status} | SHA 与重新生成 PNG 一致 | test_particle3_artifacts.py |

开发中修复了三个实现问题：面法向混入数值最近点偏差、细长三角形的固定重心坐标阈值误判、半径搜索在全局误差尚未闭合时过早丢弃小区间。
固定阈值没有为通过测试扩大，失败日志和永久回归测试均保留。
独立数值参考求解器也增加了原始/对偶夹界认证，避免只相信优化器状态字符串；没有改变生产几何公式。

## 10. 人工审核图

{figures}

## 11. 当前限制

- WALL 是刚性静止壁；没有 glycocalyx、物理涂层、wall lubrication 或 near-field resistance。
- hard contact 是 kinematic constraint，不输出真实 contact force，没有反弹或摩擦。
- deformation surrogate 不是真实膜力学，不计算膜应力、FSI、真实压降或皱褶；A_budget 是模型推导值，不是直接膜面积测量。
- CAPILLARY_DEFORMED 中保留 q，但胶囊只用局部自由流向；JEFFERY_ORIENTATION_NOT_INTERPRETED_WHILE_DEFORMED。
- 当前真实 RBC 在入口即进入变形或不可行，没有可比较的上游 FREE_OBLATE 段；不能把冻结 q 的零差叫作姿态收敛。P2 原有明显 dt 敏感性仍然保留。
- 没有 transit-time penalty 或固定 60%–80% 减速；这些与 squeeze resistance 留给 Particle-5。
- 没有 RBC→MB margination、particle-particle contact、连续注入、hematocrit 或 LAMMPS。
- 没有 production timestep；三个 dt 只做验证比较，不是生产步长选择。
- Frozen FEM 仍然 one-way。RBC 强烈堵塞时真实流场会变化，本模型不会重新求 FEM。
- 真实 WALL patch 控制试验只验证局部有限特征接触；全局放不下的摆放另列负 gap。真实 FEM 使用整幅 WALL。
- 准静态模式切换和最大半径路径不描述连续膜变形动力学。代理不可行是模型结果，不能替代生理堵塞结论。

## 12. 人工审核状态

AUTOMATED_CHECKS = {status}

MANUAL_VISUAL_REVIEW = PENDING_USER_REVIEW

PRODUCTION_WALL_MODEL_STATUS = V0_VALIDATION_ONLY_PENDING_USER_REVIEW

PRODUCTION_PARTICLE_TIMESTEP_FROZEN = false

NO_CFD_EXECUTED = true

PARTICLE4_STARTED = false

Particle-3 在此停止。仅允许完成本地提交，不 push，不 merge main，不进入 Particle-4；人工 PASS 留给用户。
'''
    (REPORT/'PARTICLE3_REVIEW.md').write_text(text)


def persist(v):
    write_json(REPORT/'PARTICLE3_VALIDATION.json',v);write_review(v)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');p.add_argument('--bind-commit',action='store_true');p.add_argument('--precheck',action='store_true');args=p.parse_args()
    LOGS.mkdir(parents=True,exist_ok=True)
    if args.precheck:
        # Independent suites can finish while the long P3 trajectory sweep runs.
        # Reuse is permitted only when the ENTIRE source snapshot still matches.
        snapshot=source_hashes();commands=[];tests={}
        assert git('rev-parse','HEAD',cwd=HANDOFF)==FROZEN and not git('status','--porcelain',cwd=HANDOFF)
        for name,script in [('frozen_validate','validate_frozen.py'),('frozen_manifest','verify_manifest.py')]:
            commands.append(run(name,[sys.executable,'-B','scripts/fem_freeze_sync/'+script],FEM))
        paths=['tests/test_fem_sync_manifest.py']+[str(p.relative_to(HANDOFF)) for p in sorted((HANDOFF/'tests').glob('test_frozen_*.py'))]+['tests/test_particle_handoff_contract.py','tests/test_no_particle_implementation_yet.py','tests/test_main_history_preserved.py']
        suites=[('frozen',HANDOFF,paths)]+[(f'particle{i}',REPO,[f'particle_3d/tests/particle{i}']) for i in range(3)]
        for name,cwd,paths in suites:
            xml=LOGS/(name+'_junit.xml');commands.append(run(name+'_pytest',[sys.executable,'-B','-m','pytest','-q','-p','no:cacheprovider',f'--junitxml={xml}',*paths],cwd));tests[name]=junit(xml);print(name,tests[name],flush=True)
        assert source_hashes()==snapshot
        write_json(LOGS/'independent_regressions.json',dict(source_sha256=snapshot,commands=commands,tests=tests))
        if any(c['returncode'] for c in commands):raise SystemExit(1)
        return
    if args.bind_commit:
        v=json.loads((REPORT/'PARTICLE3_VALIDATION.json').read_text())
        assert v['automated_checks']=='PASS' and source_hashes()==v['source_sha256']
        assert not git('status','--porcelain','--','particle_3d/src','particle_3d/tests','particle_3d/scripts','particle_3d/PARTICLE3_README.md')
        v.update(git_commit=git('rev-parse','HEAD'),source_commit_status='EXACT_TESTED_SOURCE_COMMIT; source hashes unchanged since full regression')
        persist(v);print('Bound tested source commit',v['git_commit']);return
    v=collect();persist(v)
    if not args.check:return
    assert v['git_branch']==P3_BRANCH
    assert git('rev-parse','HEAD',cwd=HANDOFF)==FROZEN and not git('status','--porcelain',cwd=HANDOFF)
    assert git('branch','--show-current',cwd=HANDOFF)=='sync/fem-simvascular-stage-q-particle-handoff-20260920'
    commands=[];tests={};precheck=LOGS/'independent_regressions.json'
    prior=json.loads(precheck.read_text()) if precheck.exists() else None
    reusable=prior is not None and prior['source_sha256']==source_hashes() and all(c['returncode']==0 for c in prior['commands'])
    if reusable:commands=prior['commands'];tests=prior['tests']
    else:
        for name,script in [('frozen_validate','validate_frozen.py'),('frozen_manifest','verify_manifest.py')]:
            commands.append(run(name,[sys.executable,'-B','scripts/fem_freeze_sync/'+script],FEM))
    frozen_paths=['tests/test_fem_sync_manifest.py']+[str(p.relative_to(HANDOFF)) for p in sorted((HANDOFF/'tests').glob('test_frozen_*.py'))]+['tests/test_particle_handoff_contract.py','tests/test_no_particle_implementation_yet.py','tests/test_main_history_preserved.py']
    suites=[('frozen',HANDOFF,frozen_paths)]+[(f'particle{i}',REPO,[f'particle_3d/tests/particle{i}']) for i in range(4)]
    for name,cwd,paths in suites:
        if reusable and name in tests:continue
        xml=LOGS/(name+'_junit.xml');commands.append(run(name+'_pytest',[sys.executable,'-B','-m','pytest','-q','-p','no:cacheprovider',f'--junitxml={xml}',*paths],cwd));tests[name]=junit(xml)
        print(name,tests[name],flush=True)
    write_json(LOGS/'final_check_commands.json',commands)
    # All data checks are exercised in permanent tests. Retain explicit status
    # on unexpected real-trajectory stops instead of laundering finite states.
    allowed={'OUTLET_01','OUTLET_02','OUTLET_03','DEFORMATION_SURROGATE_INFEASIBLE','VALIDATION_HORIZON_REACHED'}
    passed=all(c['returncode']==0 for c in commands) and all(s['status'] in allowed for s in v['real_fem_mb_results']+v['real_fem_rbc_results'])
    status='PASS' if passed else 'FAIL'
    assert source_hashes()==v['source_sha256'],'Source changed during regression'
    for key in ['automated_checks','sphere_gap_tests','ellipsoid_gap_tests','capsule_gap_tests','no_tunneling_test']:v[key]=status
    v['real_wall_tests']['status']=status
    v.update(test_results=tests,p0_regression=tests['particle0'],p1_regression=tests['particle1'],p2_regression=tests['particle2'],p3_regression=tests['particle3'],
        frozen_integrity=dict(status='PASS' if all(c['returncode']==0 for c in commands[:2]) else 'FAIL',payload_files=4602,scientific_manifest_files=31,handoff_tests=tests['frozen'],original_checkout_unchanged=not git('status','--porcelain',cwd=HANDOFF)),
        command_evidence='logs/final_check_commands.json',log_sha256={str(p.relative_to(REPO)):sha256(p) for p in sorted(LOGS.iterdir()) if p.is_file()})
    persist(v)
    (REPORT/'FILES_CHANGED.txt').write_text(git('diff','--name-only','705dd0c')+'\n'+git('ls-files','--others','--exclude-standard')+'\n')
    print('AUTOMATED_CHECKS =',status)
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
