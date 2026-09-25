#!/usr/bin/env python3
"""Final P0–P4/Frozen regression, machine-readable audit and Chinese review."""
from pathlib import Path
import argparse,concurrent.futures,json,subprocess,sys,xml.etree.ElementTree as ET
import numpy as np
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parent
sys.path.insert(0,str(PACKAGE/'src'))
from particle_3d.audit import sha256
from particle_3d.particle3_cases import write_json
from particle_3d.particle4_cases import P4_BRANCH,P3_COMMIT,P3_ACCEPTANCE_COMMIT
REPORT=PACKAGE/'reports/particle4';DATA=REPORT/'data';LOGS=REPORT/'logs'
HANDOFF=Path('/home/lzy/projects/ulm-3d-vascular-model-generation/formal_3D_flow_solver/FEM_SimVascular')
FEM=REPO/'formal_3D_flow_solver/FEM_SimVascular';FROZEN='c83dccea0f1fe5cd9fd3f5939e2eba883aeba9c2'


def git(*args,cwd=REPO):return subprocess.check_output(['git',*args],cwd=cwd,text=True).strip()
def read(name):return json.loads((DATA/(name+'.json')).read_text())
def source_hashes():
    paths=[p for name in ['src','scripts','tests','contracts'] for p in (PACKAGE/name).rglob('*') if p.is_file() and p.suffix in ['.py','.json','.md','.toml']]
    paths += list(PACKAGE.glob('PARTICLE*_README.md'))+[PACKAGE/'pyproject.toml']
    return {str(p.relative_to(REPO)):sha256(p) for p in sorted(paths)}


def collect():
    scope=read('00_scope');matrix=read('01_geometry_matrix');reference=read('01_independent_references')
    projections=[read('02_head_on_projection'),read('03_glancing_projection'),read('06_offcenter_projection')]+read('07_capsule_projection')+read('09_simultaneous_projections')+read('11_wall_pair_projection')
    audits=[p['audit'] for p in projections]
    for path in DATA.glob('*_states.json'):
        audits.extend(r['projection'] for r in json.loads(path.read_text()) if r.get('projection'))
    for a in audits:
        assert a['max_primal_violation_m_s']<=a['velocity_budget_m_s']
        assert a['max_dual_violation']<=a['multiplier_budget']
        assert a['max_complementarity_residual']<=a['complementarity_budget']
    balance=max(a['pair_translation_balance_error_m_s'] for a in audits)
    glancing=read('03_glancing_projection');g=glancing['gaps'][0];n=np.array(g['normal_j_to_i']);tangent=0.
    for dv in glancing['translation_corrections'].values():
        d=np.array(dv);tangent=max(tangent,float(np.linalg.norm(d-(d@n)*n)))
    head=read('02_head_on_projection');head_error=float(np.linalg.norm(sum(np.array(v) for v in head['translation_corrections'].values())))
    mixed=read('15_timestep_comparison');real=read('14_real_two_mb_summary');no_tunnel=read('08_no_tunnel_summary');broad=read('12_broadphase_summary');orders=read('10_order_invariance')
    result=dict(stage='Particle-4',automated_checks='PENDING_FINAL_CHECKS',git_branch=git('branch','--show-current'),git_commit=git('rev-parse','HEAD'),source_commit_status='PENDING_SOURCE_COMMIT',
        particle0_dependency_commit='7cfe5141382600e28582f05bff712a6f09c38a39',particle1_dependency_commit='6e59c605cb8251b9fa2e80d2dbed0cfa6daaf03c',
        particle2_dependency_commit='c704e08d39f6134da300fc91cc93c8561314a17b',particle3_dependency_commit=P3_COMMIT,particle3_acceptance_evidence_commit=P3_ACCEPTANCE_COMMIT,
        pair_gap_engine=scope['pair_gap_engine'],pair_geometry_tests=dict(status='PENDING_FINAL_CHECKS',shape_pairs=6,matrix_cases=len(matrix),independent_reference_cases=len(reference)),
        max_sphere_sphere_gap_error_m=max(r['error_m'] for r in reference if r['pair_type']=='Sphere-Sphere'),
        max_sphere_capsule_gap_error_m=max(r['error_m'] for r in reference if r['pair_type']=='Sphere-Capsule'),
        max_capsule_capsule_gap_error_m=max(r['error_m'] for r in reference if r['pair_type']=='Capsule-Capsule'),
        max_ellipsoid_pair_reference_error_m=max(r['error_m'] for r in reference if r['pair_type'].startswith('Ellipsoid')),
        ellipsoid_reference_error_role='MAX_OF_PRIMAL_DUAL_BRACKET_ENDPOINT_DIFFERENCE_AND_EIGHT_SECTOR_PENETRATION_REFERENCE',
        reciprocity_max_gap_error_m=max(r['reciprocity_gap_error_m'] for r in matrix),reciprocity_max_normal_error=max(r['reciprocity_normal_error'] for r in matrix),
        reciprocity_max_point_error_m=max(r['reciprocity_point_error_m'] for r in matrix),contact_model=scope['contact_model'],contact_metric=scope['contact_metric'],physical_force_claimed=False,
        max_pair_translation_balance_error_m_s=balance,head_on_symmetry_error_m_s=head_error,
        max_normal_constraint_violation_m_s=max(a['max_primal_violation_m_s'] for a in audits),max_tangential_change_error_m_s=tangent,tangential_error_scope='ISOLATED_SPHERE_CONTACT_TRANSLATIONAL_TANGENT',
        max_kkt_primal_residual=max(a['max_primal_violation_m_s'] for a in audits),max_kkt_dual_residual=max(a['max_dual_violation'] for a in audits),
        max_kkt_complementarity_residual=max(a['max_complementarity_residual'] for a in audits),contact_solve_audit_count=len(audits),
        angular_contact_correction='FREE_OBLATE_ONLY',capillary_deformed_rotation_correction=False,
        multi_contact_solver='DETERMINISTIC_DUAL_ACTIVE_SET_V0',multi_contact_order_invariant=all(r['velocity_error_m_s']==0 and r['angular_error_s_inv']==0 and r['contact_set_identical'] for r in orders),
        input_permutation_cases=len(orders),no_tunneling_test='PENDING_FINAL_CHECKS',physical_time_coverage_error_s=max(r['time_coverage_error_s'] for r in no_tunnel+mixed['cases']+[real]),
        continuous_glancing_integration='EXACT_SAME_KINEMATIC_ODE; INDEPENDENT_SOLVE_IVP_REFERENCE; NO_POSITION_PROJECTION',
        broadphase=dict(method='SUPPORT_AABB_SWEEP_WITH_SWEPT_MOTION_BOUND',**broad),all_pairs_equivalence=broad['false_negatives']==0 and all(r['gap_error_m'] in [None,0.] and r['normal_error'] in [None,0.] for r in read('12_broadphase_comparison')),
        synthetic_mixed_case=mixed,real_fem_two_mb_case=real,real_rbc_pair_validation='DEFERRED_DUE_TO_P3_REAL_RBC_PASSAGE_LIMITATION',particle3_real_rbc_passage='NOT_ESTABLISHED',
        particle_particle_lubrication='DEFERRED_PARTICLE5',collision_induced_deformation=False,production_particle_timestep_frozen=False,
        manual_visual_review='PENDING_USER_REVIEW',particle5_started=False,no_cfd_executed=True,
        wall_sha256=scope['wall_sha256'],wall_tolerance_expanded=False,distribution_modified=False,source_sha256=source_hashes(),
        data_sha256={str(p.relative_to(REPO)):sha256(p) for p in sorted(DATA.iterdir()) if p.is_file()},
        figure_sha256={str(p.relative_to(REPO)):sha256(p) for p in sorted((REPORT/'figures').glob('*.png'))})
    return result


def write_review(v):
    tests=v.get('test_results',{});table='\n'.join(f"| {k} | {r['passed']} passed | 0 failures/errors | logs/{k}_pytest.log |" for k,r in tests.items())
    figures='\n\n'.join((REPORT/f'{i:02d}_step_notes.md').read_text() for i in range(16))
    real=v['real_fem_two_mb_case'];mixed=v['synthetic_mixed_case'];status=v['automated_checks']
    text=f'''# Particle-4 人工审核报告

分支：`{v['git_branch']}`。被测源码提交：`{v['git_commit']}`（{v['source_commit_status']}）。
机器记录：[PARTICLE4_VALIDATION.json](PARTICLE4_VALIDATION.json)。算法和复现：[PARTICLE4_README.md](../../PARTICLE4_README.md)。

## 1. 这一阶段做了什么

现在粒子之间也不能互相穿过去。复用 P3 的真实有限尺寸形状，增加成对间隙、对称硬接触，以及所有 WALL 和粒子接触的同时求解。
P0–P3 原源码和测试没有修改，Frozen FEM 和正式 WALL 保持原样，没有执行 CFD。
P3 用户验收单独保存于 PARTICLE3_MANUAL_ACCEPTANCE.json；原 P3 验证 JSON 是被测时刻的历史快照。
真实 RBC 通行仍为 NOT_ESTABLISHED，不能把 DEFORMATION_SURROGATE_INFEASIBLE 当作生理堵塞。

## 2. 什么叫 pair gap

它表示两个完整粒子表面还剩多少距离。正值表示分开，零表示接触，负值表示几何相交深度。
法向从第二只指向第一只；交换输入后间隙不变，两个接触点交换，法向反向。
计算使用原 shape support、Minkowski 差、GJK 和支撑函数平稳点细化。穿透方法经过独立数值验证，但不声称是任意凸体的形式化全局证明。
不使用 RBC 中心距离减有效半径的替代模型；不能可靠闭合接触点时明确报错。
几何舍入界沿用 P3 的 512×float64 eps×实际尺度，没有固定比例 overlap 或物理 padding。

## 3. 哪些 shape 可以互相碰

球–球、球–椭球、球–胶囊、椭球–椭球、椭球–胶囊、胶囊–胶囊，共六类。
RBC–MB 使用原 P2 的五个分位样本及不同姿态。胶囊取自原 P3 定体积、面积预算允许的几何。
已经 DEFORMATION_SURROGATE_INFEASIBLE 的 RBC 不是可继续模拟的状态。
粒子接触不修改 RBC 尺寸、分布、面积预算或胶囊轴；没有碰撞诱导变形。

## 4. 接触以后怎么处理

不弹、不摩擦，只去掉继续互相钻入的运动。每条粒子接触约束给两边成对相反的中心平移修正。
全部约束一起选择最小的几何加权速度修正，球的自旋和胶囊保存的角速度不受法向接触改变。
接触乘子是 KINEMATIC_CONTACT_MULTIPLIER，不是 force、impulse 或真实接触力，没有 N 或 pN 单位。
最大法向违反量为 {v['max_normal_constraint_violation_m_s']:.6e} m/s；孤立球接触最大切向变化误差为 {v['max_tangential_change_error_m_s']:.6e} m/s。

## 5. 为什么 RBC 会有角速度修正

偏心接触时，转动和中心平移可以共同去掉入侵速度。FREE_OBLATE 的角速度修正由接触点位置、法向和原尺寸自动算出，不手写转动角度。
比较平移与转动时使用 ell=max(a,b,c)，目标包含 ell²×角速度修正的平方；这是几何尺度，不是质量、惯性或物理能量。
CAPILLARY_DEFORMED 的轴由 P3 局部流向定义，所以不通过粒子接触转动它。

## 6. 三个粒子一起碰怎么办

把所有当前接触放进同一个约束问题，不能先左后右地逐对推开。
直线、三角形和四粒子案例均检查非负乘子、非负法向速度及互补条件；WALL 加粒子接触也一起解。
三粒子与四粒子的 {v['input_permutation_cases']} 个排列中，稳定 ID 保持原样，速度、角速度和接触集合完全一致。
求解后独立复算实际接触点速度；退化不相容或保护条件触发时报告 MULTI_CONTACT_INFEASIBLE。

## 7. 为什么还要真实时间细分

不能一步从另一只粒子的身体里穿过去，即使大步的两个端点看起来分离。
先用运动上界和连续分离平面证明区间安全，不能证明时真实分成左右两个时间区间；右半从左半接受状态继续。
最大时间覆盖误差为 {v['physical_time_coverage_error_s']:.6e} s。没有中心位置投影、推开重叠或失败后偷偷推进时间。
孤立双球持续擦边时，使用同一接触速度方程的解析积分，避免弦推进反复离开/返回接触面的数值抖动；独立 ODE 积分已验证。
改变 P3 形状模式或胶囊轴/尺寸时，旧 pair 运动证书作废，必须重新计算路径；本轮不绕过真实 RBC 通行限制去运行真实 RBC 对。

## 8. 自动测试

| 检查 | 结果 | 最大误差 / 条件 | 测试文件或完整日志 |
|---|---|---|---|
{table}
| 六类有限尺寸几何 | {status} | 球球 {v['max_sphere_sphere_gap_error_m']:.3e} m；球胶囊 {v['max_sphere_capsule_gap_error_m']:.3e} m；胶囊胶囊 {v['max_capsule_capsule_gap_error_m']:.3e} m | test_pair_geometry.py |
| 椭球独立原始/对偶参考 | {status} | 到参考界端点最大 {v['max_ellipsoid_pair_reference_error_m']:.3e} m | test_pair_geometry.py |
| reciprocity | {status} | gap {v['reciprocity_max_gap_error_m']:.3e} m；normal {v['reciprocity_max_normal_error']:.3e} | test_pair_geometry.py |
| 成对对称、切向、自旋、偏心转动 | {status} | 平移平衡 {v['max_pair_translation_balance_error_m_s']:.3e} m/s | test_contact_projection.py |
| 同时 WALL / pair；排列不变性 | {status} | primal {v['max_kkt_primal_residual']:.3e}；dual {v['max_kkt_dual_residual']:.3e}；complementarity {v['max_kkt_complementarity_residual']:.3e} | test_contact_projection.py |
| 无穿透、真实时间、初始重叠拒绝 | {status} | 时间 {v['physical_time_coverage_error_s']:.3e} s | test_pair_time_and_broadphase.py |
| broadphase 与穷举 | {status} | 漏检 {v['broadphase']['false_negatives']} | test_pair_time_and_broadphase.py |
| 真实双球、旧输入、图与报告 | {status} | 源文件 SHA、图重生成及有限状态检查 | test_scope_and_artifacts.py |

人工混合轨迹使用 dt / dt/2 / dt/4 = 0.2 / 0.1 / 0.05 s，只作 VALIDATION_ONLY。
三类接触事件均发生；最小接受 pair gap 为 {min(r['minimum_pair_gap_m'] for r in mixed['cases']):.6e} m，处于原有坐标舍入界内。
它使用中心线接触和零自由角速度，短轴差为零不代表真实 FEM 姿态收敛。
真实双球烟雾结果：{real['pair_contact_status']}；共 {real['rows']} 个接受状态，至 {real['last_accepted_time_s']:.9g} s。
真实 pair 最小 gap 为 {real['minimum_pair_gap_m']:.6e} m；WALL 最小 gap 为 {real['minimum_wall_gap_m']:.6e} m；出口事件：{real['outlet_events'] or '验证窗口内没有出口事件'}。
使用原 SonoVue 两个固定种子 20260920 / 20260921，没有缩小第二只 MB 或添加人工力。
初值取原 P3 路径中后续 256 个参考点均可容纳原大球的第一个宽段；最初仅检查瞬时可行的候选发生后续大球碰壁和大量 P3 子步，该开发尝试中止并保留日志。
正式烟雾检查不宣称验证了大球持续沿曲壁滑动。
开发中修复了 WALL 返回元组的适配、初始化搜索窗口过短和双球擦边积分抖动，失败日志保留；没有改变接触模型或科学阈值。
P4 平移更新精确保留原胶囊轴，避免重复构造时的单位化改变末位浮点数；独立接触步回归验证轴、自旋和尺寸不变。
补充深重叠检查时，独立参考优化器曾陷入较深的另一方向极值；参考改为覆盖八个角度区域的独立搜索，未使用正式解作为初值，也未放宽阈值。

## 9. 人工审核图

{figures}

## 10. 当前限制

- no contact force；no mass/inertia collision；乘子只是运动学求解量。
- no friction；no restitution；no lubrication；no particle-particle hydrodynamics。
- no collision-induced deformation；不修改 P2 distribution 或 P3 deformation。
- CAPILLARY_DEFORMED 不接受 contact angular rotation；模式或轴变化必须重新认证路径。
- no adhesion；no LAMMPS；no continuous injection 或 hematocrit。
- no production timestep；全部步长仅用于验证。
- real RBC passage from P3 is still NOT_ESTABLISHED；不能作生理堵塞结论。
- 真实 FEM two-MB smoke 不等于正式 suspension，也不表示这两个尺寸已完成整条血管通行。
- penetration 使用经过独立参考验证的固定种子 support-stationary 方法；没有任意凸体的形式化全局最优证明。
- 数值保护触发时明确停止，不允许重叠或静止状态消耗失败时间。

## 11. 人工审核状态

AUTOMATED_CHECKS = {status}

MANUAL_VISUAL_REVIEW = PENDING_USER_REVIEW

PRODUCTION_PARTICLE_TIMESTEP_FROZEN = false

NO_CFD_EXECUTED = true

PARTICLE5_STARTED = false

Particle-4 在此停止，仅本地提交，不 push，不 merge main，不进入 Particle-5。
'''
    (REPORT/'PARTICLE4_REVIEW.md').write_text(text)


def persist(v):write_json(REPORT/'PARTICLE4_VALIDATION.json',v);write_review(v)


def run(job):
    name,cwd,command=job;path=LOGS/(name+'.log')
    with path.open('w') as stream:
        stream.write('$ '+' '.join(map(str,command))+'\n');stream.flush();result=subprocess.run(command,cwd=cwd,stdout=stream,stderr=subprocess.STDOUT)
    print(name,result.returncode,flush=True)
    return dict(name=name,cwd=str(cwd),command=list(map(str,command)),returncode=result.returncode,log=str(path.relative_to(REPO)))


def junit(path):
    root=ET.parse(path).getroot();suites=[root] if root.tag=='testsuite' else list(root.findall('testsuite'))
    d={k:sum(int(s.attrib.get(k,0)) for s in suites) for k in ['tests','errors','failures','skipped']};d['passed']=d['tests']-d['errors']-d['failures']-d['skipped'];d['status']='PASS' if d['errors']==d['failures']==0 else 'FAIL';return d


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');p.add_argument('--bind-commit',action='store_true');args=p.parse_args();LOGS.mkdir(parents=True,exist_ok=True)
    if args.bind_commit:
        v=json.loads((REPORT/'PARTICLE4_VALIDATION.json').read_text());assert v['automated_checks']=='PASS' and source_hashes()==v['source_sha256']
        assert not git('status','--porcelain','--','particle_3d/src','particle_3d/scripts','particle_3d/tests','particle_3d/PARTICLE4_README.md')
        v.update(git_commit=git('rev-parse','HEAD'),source_commit_status='EXACT_TESTED_SOURCE_COMMIT; source hashes unchanged since full regression');persist(v);print(v['git_commit']);return
    v=collect();persist(v)
    if not args.check:return
    assert v['git_branch']==P4_BRANCH
    assert git('rev-parse','HEAD',cwd=HANDOFF)==FROZEN and not git('status','--porcelain',cwd=HANDOFF)
    jobs=[]
    for name in ['validate_frozen','verify_manifest']:jobs.append((name,FEM,[sys.executable,'-B','scripts/fem_freeze_sync/'+name+'.py']))
    paths=['tests/test_fem_sync_manifest.py']+[str(p.relative_to(HANDOFF)) for p in sorted((HANDOFF/'tests').glob('test_frozen_*.py'))]+['tests/test_particle_handoff_contract.py','tests/test_no_particle_implementation_yet.py','tests/test_main_history_preserved.py']
    suites=[('frozen',HANDOFF,paths)]+[(f'particle{i}',REPO,[f'particle_3d/tests/particle{i}']) for i in range(5)]
    for name,cwd,paths in suites:jobs.append((name+'_pytest',cwd,[sys.executable,'-B','-m','pytest','-q','-p','no:cacheprovider',f'--junitxml={LOGS}/{name}_junit.xml',*paths]))
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:commands=list(pool.map(run,jobs))
    write_json(LOGS/'final_check_commands.json',commands);tests={name:junit(LOGS/(name+'_junit.xml')) for name,_,_ in suites}
    assert source_hashes()==v['source_sha256'],'Source changed during regression'
    status='PASS' if all(r['returncode']==0 for r in commands) else 'FAIL';v['automated_checks']=status;v['pair_geometry_tests']['status']=status;v['no_tunneling_test']=status
    v.update(test_results=tests,**{f'p{i}_regression':tests[f'particle{i}'] for i in range(5)},frozen_integrity=dict(status='PASS' if all(c['returncode']==0 for c in commands[:2]) else 'FAIL',payload_files=4602,scientific_manifest_files=31,handoff_tests=tests['frozen'],original_checkout_unchanged=not git('status','--porcelain',cwd=HANDOFF)))
    v['log_sha256']={str(p.relative_to(REPO)):sha256(p) for p in sorted(LOGS.iterdir()) if p.is_file()};persist(v)
    (REPORT/'FILES_CHANGED.txt').write_text(git('diff','--name-only',P3_ACCEPTANCE_COMMIT)+'\n'+git('ls-files','--others','--exclude-standard')+'\n')
    print('AUTOMATED_CHECKS =',status)
    if status!='PASS':raise SystemExit(1)


if __name__=='__main__':main()
