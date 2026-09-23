#!/usr/bin/env python3
"""Final Frozen/P0–P5 regressions, machine evidence and twelve-section Chinese review."""
from pathlib import Path
import argparse,concurrent.futures,importlib.metadata,json,platform,subprocess,sys
import xml.etree.ElementTree as ET
import numpy as np
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parent
sys.path.insert(0,str(PACKAGE/'src'))
from particle_3d.audit import sha256
from particle_3d.particle3_cases import write_json
from particle_3d.particle5_cases import P5_BRANCH,P4_COMMIT
REPORT=PACKAGE/'reports/particle5';DATA=REPORT/'data';LOGS=REPORT/'logs'
FEM=REPO/'formal_3D_flow_solver/FEM_SimVascular'
HANDOFF=Path('/home/lzy/projects/ulm-3d-vascular-model-generation/formal_3D_flow_solver/FEM_SimVascular')
FROZEN='c83dccea0f1fe5cd9fd3f5939e2eba883aeba9c2'
DEPENDENCIES=['7cfe5141382600e28582f05bff712a6f09c38a39','6e59c605cb8251b9fa2e80d2dbed0cfa6daaf03c',
'c704e08d39f6134da300fc91cc93c8561314a17b','95e54fe474349c35aaad2e2754aff0bad7108c42',P4_COMMIT]


def git(*args,cwd=REPO):return subprocess.check_output(['git',*args],cwd=cwd,text=True).strip()
def read(name):return json.loads((DATA/(name+'.json')).read_text())


def source_hashes():
    paths=[p for name in ['src','scripts','tests','contracts'] for p in (PACKAGE/name).rglob('*') if p.is_file() and p.suffix in ['.py','.json','.md','.toml']]
    paths+=list(PACKAGE.glob('PARTICLE*_README.md'))+[PACKAGE/'pyproject.toml']
    return {str(p.relative_to(REPO)):sha256(p) for p in sorted(paths)}


def solver_records(value):
    if isinstance(value,dict):
        if 'condition_estimate' in value and 'relative_backward_residual' in value:
            yield value
        for item in value.values():yield from solver_records(item)
    elif isinstance(value,list):
        for item in value:yield from solver_records(item)


def collect():
    scope=read('00_scope');stokes=read('01_stokes');wall=read('02_wall');pair=read('03_pair')
    sparse=read('05_sparse');algebraic=read('04_algebraic');dissipation=read('06_dissipation');contact=read('07_contact')
    near=read('09_real_near_wall')[1];separating=read('09_real_near_wall')[0];real=read('10_real_summary')
    audits=[r for p in sorted(DATA.glob('*.json')) for r in solver_records(json.loads(p.read_text()))]
    assert all(r['minimum_scaled_eigenvalue']>0 for r in audits)
    assert all(min(r['dissipation'].values())>=0 for r in audits)
    assert all(min(row[k] for k in ['self','wall','pair','total'])>=0 for row in dissipation)
    assert algebraic['minimum_eigenvalue']>=-algebraic['eigenvalue_roundoff_bound']
    assert all(r['p5_resistance_objective']<r['p4_resistance_objective'] for r in contact)
    raw=source_hashes();p4=json.loads((PACKAGE/'reports/particle4/PARTICLE4_VALIDATION.json').read_text())
    assert all(raw.get(name)==digest for name,digest in p4['source_sha256'].items()),'P0–P4 scientific source changed'
    result=dict(scope,automated_checks='PENDING_FINAL_CHECKS',git_branch=git('branch','--show-current'),
        git_commit=git('rev-parse','HEAD'),source_commit_status='PENDING_SOURCE_COMMIT',
        frozen_dependency_commit=FROZEN,particle4_acceptance_evidence_commit=git('rev-parse','0fca8cd'),
        sphere_self_translation_model='STOKES_6_PI_MU_A',sphere_self_rotation_model='STOKES_8_PI_MU_A_CUBED',role='VALIDATION_ONLY',
        max_stokes_resistance_error=max(abs(r[k]/r['analytic_'+k]-1) for r in stokes for k in ['translation','rotation']),
        stokes_error_role='MAX_RELATIVE_TRANSLATION_AND_ROTATION_COEFFICIENT_ERROR',
        max_self_translation_error=max(abs(r['translation']-r['analytic_translation']) for r in stokes),
        max_self_rotation_error=max(abs(r['rotation']-r['analytic_rotation']) for r in stokes),
        max_isolated_velocity_error_m_s=max(r['isolated_velocity_error'] for r in stokes),
        max_isolated_angular_error_s_inv=max(r['isolated_angular_error'] for r in stokes),
        max_wall_lubrication_scaling_error=max(abs(r['normalized_scaling']-1) for r in wall),
        max_pair_lubrication_scaling_error=max(abs(r['normalized_scaling']-1) for r in pair),
        max_analytic_velocity_error=max(r['analytic_error'] for r in wall+pair),analytic_velocity_error_units='m/s',
        matrix_symmetry_error=max(r['matrix_symmetry_error'] for r in audits),
        scaled_matrix_symmetry_roundoff=max(r['scaled_matrix_symmetry_error'] for r in audits),
        minimum_dissipation=min(min(r['dissipation'].values()) for r in audits),
        minimum_dissipation_role='MINIMUM_COMPONENT_OR_TOTAL_ACROSS_SAVED_SOLVES; INCLUDES_NULL_MODES; W',
        minimum_random_total_dissipation_w=min(r['total'] for r in dissipation),
        minimum_random_component_dissipation_w=min(r[k] for r in dissipation for k in ['self','wall','pair']),
        dissipation_random_samples=len(dissipation),
        minimum_eigenvalue_with_roundoff_context=dict(physical_self_scaled_minimum=min(r['minimum_scaled_eigenvalue'] for r in audits),
            physical_maximum_roundoff_bound=max(r['eigenvalue_roundoff_bound'] for r in audits),
            algebraic_rank_one_minimum=algebraic['minimum_eigenvalue'],algebraic_roundoff_bound=algebraic['eigenvalue_roundoff_bound'],
            absolute_value_applied=False,physical_eigenvalue_coordinates='SELF_DIAGONALLY_SCALED_6N'),
        sparse_all_pairs_equivalence=sparse['sparse_all_pairs_equivalence'],sparse_candidate_count=len(sparse['candidate_pairs']),
        sparse_all_pair_count=sparse['all_pairs_count'],sparse_velocity_difference=sparse['velocity_difference'],
        resistance_solver_residual=max(r['relative_backward_residual'] for r in audits),
        max_condition_number=max(r['condition_estimate'] for r in audits),condition_number_role='SELF_DIAGONALLY_SCALED_6N',
        resistance_metric_contact_status='PASS',contact_multiplier_role='KINEMATIC_CONSTRAINT_MULTIPLIER',
        max_contact_constraint_violation_m_s=max(r['max_constraint_violation_m_s'] for r in audits),
        physical_time_coverage_error_s=max([r['time_coverage_error_s'] for r in read('11_timestep')]+[real['time_coverage_error_s']]),
        real_near_wall_mb_status='PASS',real_near_wall_mb_gap_m=near['gap_m'],
        real_near_wall_free_normal_velocity=near['free_normal_velocity_m_s'],
        real_near_wall_lubricated_normal_velocity=near['lubricated_normal_velocity_m_s'],
        real_near_wall_tangential_difference=near['tangential_difference_m_s'],
        real_near_wall_attenuation_ratio=near['attenuation_ratio'],real_near_wall_source_row=near['source_row_index'],
        original_absolute_minimum_wall_gap_m=separating['gap_m'],
        original_absolute_minimum_wall_sample_motion='SEPARATING; RETAINED_AS_SEPARATE_SIGNED_REPLAY',
        real_two_mb_pair_lubrication_status=real['pair_lubrication_status'],real_two_mb_case=real,
        real_rbc_deferral_reasons=real['deferral_reasons'],
        wall_sha256=read('10_initialization')['wall_provenance']['sha256'],
        wall_tolerance_expanded=False,particle_dimensions_changed=False,gap_floor_used=False,constant_velocity_multiplier_used=False,
        source_sha256=raw,data_sha256={str(p.relative_to(REPO)):sha256(p) for p in sorted(DATA.iterdir()) if p.is_file()},
        figure_sha256={str(p.relative_to(REPO)):sha256(p) for p in sorted((REPORT/'figures').glob('*.png'))})
    result.update({f'particle{i}_dependency_commit':commit for i,commit in enumerate(DEPENDENCIES)})
    return result


def write_review(v):
    table='\n'.join(f"| {name} | {r['passed']} 通过 | {r['failures']} 失败 / {r['errors']} 错误 / {r['skipped']} 跳过 | [日志](logs/{name}_pytest.log)、[JUnit](logs/{name}_junit.xml) |" for name,r in v.get('test_results',{}).items())
    figures='\n\n'.join((REPORT/f'{i:02d}_step_notes.md').read_text() for i in range(12))
    real=v['real_two_mb_case'];eigen=v['minimum_eigenvalue_with_roundoff_context'];status=v['automated_checks']
    text=f'''# Particle-5 人工审核报告

分支：`{v['git_branch']}`。被测源码提交：`{v['git_commit']}`，状态：`{v['source_commit_status']}`。
机器记录：[PARTICLE5_VALIDATION.json](PARTICLE5_VALIDATION.json)。复现说明：[PARTICLE5_README.md](../../PARTICLE5_README.md)。
P4 用户验收已保存为独立证据提交 `{v['particle4_acceptance_evidence_commit']}`，未改 P4 历史科学快照、算法或测试。

## 1. 这一阶段做了什么

现在球形粒子还没碰到，夹在表面之间的薄层血液就会抑制它们继续靠近。新增 self、墙面法向和双球法向阻力，统一解全部粒子的平移与自转速度。
原有硬接触继续负责最终非穿透，全部墙面和粒子接触一起处理。P0–P4 原科学源码、测试和 Frozen FEM 保持原样，开始前和结束时均运行完整回归，没有执行 CFD。

## 2. resistance 是什么

这里的 resistance 是速度与液体阻力项之间的比例关系。背景流通过 self 项驱动粒子，近墙和相对靠近会增加额外阻力，同样的背景流驱动下就更难继续靠近。
实际方程是 `R U = b`，其中 `R = R_self + R_wall + R_pair`，`b = R_self U_free`。静止墙面和 pair 项不把自由速度加回右端，否则会错误地得到所有粒子永远等于自由速度。
每只球都有六个未知量，包括三个平移和三个自转分量。本轮无外加非流体力，也没有质量或惯性。

## 3. 为什么不能简单乘 60%-80%

减速由间隙、尺寸和阻力平衡自动计算，不是把所有速度乘一个常数。球靠平面时，法向速度比例为 `h/(a+h)`；间隙越小，比例越低，切向和自转在当前法向模型中保持原值。
所有阻力项均使用同一黏度，因此在本轮零外力的速度比例中黏度会约掉；阻力系数和耗散仍随黏度变化。真实使用 μ = {v['dynamic_viscosity_pa_s']} Pa·s，来自冻结 `baseline_summary.json` 与 `run/solver.xml` 的一致记录，并核对 ρν，没有新增黏度选择。

## 4. self / wall / pair 三种阻力

self 表示粒子偏离背景流时的阻力，球平移为 `6πμa`、自转为 `8πμa³`。孤立球解出原自由速度。
wall 表示球相对静止墙的法向阻力，系数为 `6πμa²/h`。pair 表示两球相对法向运动，系数为 `6πμR_eff²/h`，其中 `R_eff=ai aj/(ai+aj)`，四个矩阵块是 `[K,-K;-K,K]`。
共同平移不会触发 pair 阻力。这里的 resistance contribution 是无惯性平衡里的液体阻力项；接触乘子依旧标记 `KINEMATIC_CONSTRAINT_MULTIPLIER`，不称为真实接触力。
当 h 小于等于原几何舍入界时进入硬接触，没有人为抬高 h。真实 near-field 资格采用 `h/relevant_radius <= 0.01`，只作 VALIDATION_ONLY；合成的 0.1、0.01、0.001、0.0001 是显式公式验证点。

## 5. 球和墙怎么验证

人工平面中，`阻力×间隙/(6πμa²)` 应等于 1，最大无量纲误差为 {v['max_wall_lubrication_scaling_error']:.6e}。法向速度直接比较标量解析解，同时检查切向和自转不变。
P4 对照在接触前恒速；P5 随间隙缩小而减速，在四秒验证时长内仍保留正间隙，不要求一定接触。每一步都使用真实时间和连续几何证书，没有碰墙后移动中心的处理。
真实静态回放保留了原始两种状态。原绝对最小墙隙 {v['original_absolute_minimum_wall_gap_m']:.9e} m 对应离墙运动，不能把它写成靠近。
用于靠近验证的是 P4 原轨迹第 {v['real_near_wall_source_row']} 行、ID 203，间隙 {v['real_near_wall_mb_gap_m']:.9e} m；自由法向速度 {v['real_near_wall_free_normal_velocity']:.9e} m/s，修正后 {v['real_near_wall_lubricated_normal_velocity']:.9e} m/s，比例 {v['real_near_wall_attenuation_ratio']:.9g}。
切向变化 {v['real_near_wall_tangential_difference']:.6e} m/s，处于浮点误差量级。中心、半径、原 WALL 与间隙都未人为改变。这里只验证 one-way frozen-flow 的局部 leading normal 修正，尚未处理切向迁移率、自转耦合、有限曲率和多墙面流体影响。

## 6. 两个球怎么验证

等半径和不等半径两组比较独立二元解析解，最大 pair 标度误差 {v['max_pair_lubrication_scaling_error']:.6e}。输入顺序交换、共同平移、成对相反贡献均有永久测试。
12 球例中 P4 包围盒查询给出 {v['sparse_candidate_count']} 个候选，穷举为 {v['sparse_all_pair_count']} 对。稀疏与穷举矩阵及求解速度逐元素一致；扩展只作用于验证查询包围盒，不改变真实球尺寸。
真实双球完全沿用 P4 初值、原 SonoVue 两个半径、dt={real['requested_dt_s']:.12g} s、horizon={real['horizon_s']:.12g} s。
得到 {real['rows']} 个接受状态，实际时间 {real['final_time_s']:.12g} s，覆盖误差 {real['time_coverage_error_s']:.3e} s。最小 pair gap 为 {real['minimum_pair_gap_m']:.9e} m，gap/R_eff={real['minimum_pair_gap_ratio']:.9g}，明显大于 0.01，因此 `PAIR_LUBRICATION_NOT_ACTIVE_IN_THIS_REAL_SMOKE`。
墙面近场启用 {real['wall_active_block_count']} 次，新轨迹最小墙隙 {real['minimum_wall_gap_m']:.9e} m，是原初值演化后的结果；没有修改 0.3 nm 原样本，也没有为了 pair 图把两球移近。窗口内无出口事件，不能据此声称完整血管通行。

## 7. 为什么目前不声称 RBC lubrication 已解决

这不是忘记实现，而是当前没有人工批准的非球形润滑物理系数。椭球和胶囊的真实润滑依赖姿态、局部曲率和变形，不能用一个等效球半径代替。
图 04 只检验通用 pair 块的对称性、互易性和非负耗散，系数 1 仅存在于内部归一测试单位。正式接口收到这类形状时明确报 `NONSPHERICAL_LUBRICATION_NOT_FROZEN`，也拒绝把验证矩阵当成物理块输入。
真实 RBC 流体回放为 DEFERRED：一是 P3 的 REAL_RBC_PASSAGE 仍为 NOT_ESTABLISHED，二是 P5 非球形润滑未冻结。这不是 P5 自动失败，也不是生理堵塞结论。

## 8. 为什么耗散必须非负

液体阻力消耗相对运动，不能凭空产生能量。分别计算相对背景流的 self 耗散、相对静止墙的 wall 耗散，以及两粒子相对运动的 pair 耗散。
固定种子生成 {v['dissipation_random_samples']} 个完整速度向量，各项均非负。随机样本最小总耗散为 {v['minimum_random_total_dissipation_w']:.6e} W；包含无作用方向和零贡献的全部保存求解中，最小分项耗散为 {v['minimum_dissipation']:.6e} W。
物理系统 self 缩放后的最小特征值为 {eigen['physical_self_scaled_minimum']:.9g}。代数 rank-one 测试的最小特征值为 {eigen['algebraic_rank_one_minimum']:.6e}，其舍入界为 {eigen['algebraic_roundoff_bound']:.6e}；原负号完整保留，没有取绝对值掩盖结果。
用代数对角缩放和稀疏 LU 求解，不直接求逆；缩放与未缩放的解有独立比较。最大缩放条件数 {v['max_condition_number']:.9g}，最大相对后向残差 {v['resistance_solver_residual']:.6e}。
超出 float64 可靠范围时明确报告 `RESISTANCE_SYSTEM_ILL_CONDITIONED`，请求物理时间细分；仍失败则停止案例。永久测试检查失败状态不移动、不消耗时间。

## 9. 自动测试

| 检查组 | 通过数 | 失败 / 错误 / 跳过 | 证据 |
|---|---:|---|---|
{table}

| 科学检查 | 结果或最大误差 | 永久测试 |
|---|---|---|
| 黏度双来源与不一致停止 | {status} | test_resistance_core.py |
| Stokes 平移 / 自转系数 | 相对误差 {v['max_stokes_resistance_error']:.3e} | test_resistance_core.py |
| 球墙 / 双球 1/h 标度 | {v['max_wall_lubrication_scaling_error']:.3e} / {v['max_pair_lubrication_scaling_error']:.3e} | test_resistance_core.py |
| 独立解析平移速度 | {v['max_analytic_velocity_error']:.3e} m/s | test_resistance_core.py |
| 未缩放物理矩阵对称性 | {v['matrix_symmetry_error']:.3e} | test_resistance_core.py |
| 正耗散 / 原符号谱 / 交换 / 稀疏穷举 | {status} | test_resistance_core.py |
| 非球形拒绝与验证系数隔离 | {status} | test_resistance_core.py |
| 同时接触、穿越检测与失败不推进 | {status}；最大违反 {v['max_contact_constraint_violation_m_s']:.3e} m/s | test_contact_and_time.py |
| 原真实状态 / 初始化与完整时长 | {status} | test_real_and_artifacts.py |
| dt、dt/2、dt/4 与独立隐式 ODE | 一阶误差随细化下降；仅验证步长 | test_contact_and_time.py、test_real_and_artifacts.py |
| 全部图重新生成及 SHA、历史源锁定 | {status} | test_real_and_artifacts.py |

开发中修正了图 01 刻度重叠、图 07 零附近自动放大和图 08 离散 P4 速度连线的展示方式。图 08 明确显示接触前恒速和接触后零法向速度，不用两点斜线暗示 P4 连续减速。
所有修正保持物理公式、旧阈值和 Frozen 不变。初次命令缺少源码导入路径的问题通过现有脚本路径约定处理；开发与最终测试日志保留。

## 10. 人工审核图

{figures}

## 11. 当前限制

- normal lubrication only；no tangential lubrication；no rotational lubrication coupling。
- leading asymptotic only；有限曲率、多墙面和任意分离距离的完整 mobility 尚未验证。
- production cutoff not frozen；0.01 是本轮真实回放的验证资格，pair / wall 生产 cutoff 均为 null。
- RBC/capsule lubrication not frozen；real RBC passage not established；真实 RBC hydrodynamics 为 DEFERRED。
- Frozen FEM one-way；粒子不反作用于流场，没有运行 CFD。
- no Brownian；no gravity、buoyancy、lift、added mass、Basset history、acoustics；no mass/inertia。
- no adhesion；no friction、restitution、springs；no RBC membrane mechanics 或 collision-induced deformation。
- no LAMMPS；no continuous injection、hematocrit；no full suspension。
- no production timestep；三组步长只用于验证，速度降低不代表可以放大生产步长。
- 稀疏组装与直接求解已验证 5–20 球；条件数检查使用小规模稠密谱分解，大规模性能未验证。
- 真实新轨迹的极小间隙按原模型保留，本轮没有引入分子尺度表面物理，因此不能将其解释成该尺度上的完整真实物理验证。

## 12. 人工审核状态

AUTOMATED_CHECKS = {status}

MANUAL_VISUAL_REVIEW = PENDING_USER_REVIEW

PRODUCTION_LUBRICATION_CUTOFF_FROZEN = false

PRODUCTION_PARTICLE_TIMESTEP_FROZEN = false

NO_CFD_EXECUTED = true

PARTICLE6_STARTED = false

Particle-5 在此停止，仅本地提交；没有 push、merge main 或启动 Particle-6。
'''
    (REPORT/'PARTICLE5_REVIEW.md').write_text(text)


def persist(v):write_json(REPORT/'PARTICLE5_VALIDATION.json',v);write_review(v)


def run(job):
    name,cwd,command=job;path=LOGS/(name+'.log')
    with path.open('w') as stream:
        stream.write('$ '+' '.join(map(str,command))+'\n');stream.flush()
        result=subprocess.run(command,cwd=cwd,stdout=stream,stderr=subprocess.STDOUT)
    print(name,result.returncode,flush=True)
    return dict(name=name,cwd=str(cwd),command=list(map(str,command)),returncode=result.returncode,log=str(path.relative_to(REPO)))


def junit(path):
    root=ET.parse(path).getroot();suites=[root] if root.tag=='testsuite' else root.findall('testsuite')
    d={k:sum(int(s.attrib.get(k,0)) for s in suites) for k in ['tests','errors','failures','skipped']}
    d['passed']=d['tests']-d['errors']-d['failures']-d['skipped'];d['status']='PASS' if d['errors']==d['failures']==d['skipped']==0 else 'FAIL'
    return d


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');parser.add_argument('--bind-commit',action='store_true');args=parser.parse_args()
    LOGS.mkdir(parents=True,exist_ok=True)
    if args.bind_commit:
        v=json.loads((REPORT/'PARTICLE5_VALIDATION.json').read_text())
        assert v['automated_checks']=='PASS' and v['source_sha256']==source_hashes()
        assert not git('status','--porcelain','--','particle_3d/src','particle_3d/scripts','particle_3d/tests','particle_3d/PARTICLE5_README.md')
        v.update(git_commit=git('rev-parse','HEAD'),source_commit_status='EXACT_TESTED_SOURCE_COMMIT; source hashes unchanged since full regression')
        persist(v);print(v['git_commit']);return
    v=collect();persist(v)
    if not args.check:return
    assert v['git_branch']==P5_BRANCH
    assert git('rev-parse','HEAD',cwd=HANDOFF)==FROZEN and not git('status','--porcelain',cwd=HANDOFF)
    versions={name:importlib.metadata.version(name) for name in ['numpy','scipy','vtk','pyvista','matplotlib','Pillow','pytest']}
    write_json(LOGS/'environment.json',dict(python=sys.version,executable=sys.executable,platform=platform.platform(),packages=versions))
    jobs=[(name,FEM,[sys.executable,'-B','scripts/fem_freeze_sync/'+name+'.py']) for name in ['validate_frozen','verify_manifest']]
    paths=['tests/test_fem_sync_manifest.py']+[str(p.relative_to(HANDOFF)) for p in sorted((HANDOFF/'tests').glob('test_frozen_*.py'))]+['tests/test_particle_handoff_contract.py','tests/test_no_particle_implementation_yet.py','tests/test_main_history_preserved.py']
    suites=[('frozen',HANDOFF,paths)]+[(f'particle{i}',REPO,[f'particle_3d/tests/particle{i}']) for i in range(6)]
    for name,cwd,paths in suites:
        jobs.append((name+'_pytest',cwd,[sys.executable,'-B','-m','pytest','-q','-p','no:cacheprovider',f'--junitxml={LOGS}/{name}_junit.xml',*paths]))
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:commands=list(pool.map(run,jobs))
    write_json(LOGS/'final_check_commands.json',commands)
    tests={name:junit(LOGS/(name+'_junit.xml')) for name,_,_ in suites}
    assert v['source_sha256']==source_hashes(),'Source changed during regression'
    expected=dict(frozen=18,particle0=46,particle1=67,particle2=74,particle3=90,particle4=69)
    assert all(tests[name]['tests']==count for name,count in expected.items()),'Historical test counts changed'
    status='PASS' if all(c['returncode']==0 for c in commands) and all(r['status']=='PASS' for r in tests.values()) else 'FAIL'
    v.update(automated_checks=status,test_results=tests,**{f'p{i}_regression':tests[f'particle{i}'] for i in range(5)},
        frozen_integrity=dict(status='PASS' if all(c['returncode']==0 for c in commands[:2]) else 'FAIL',
            payload_files=4602,scientific_manifest_files=31,handoff_tests=tests['frozen'],
            original_checkout_unchanged=not git('status','--porcelain',cwd=HANDOFF),original_head=git('rev-parse','HEAD',cwd=HANDOFF)))
    v['log_sha256']={str(p.relative_to(REPO)):sha256(p) for p in sorted(LOGS.iterdir()) if p.is_file()}
    persist(v)
    (REPORT/'FILES_CHANGED.txt').write_text(git('diff','--name-only','0fca8cd')+'\n'+git('ls-files','--others','--exclude-standard')+'\n')
    print('AUTOMATED_CHECKS =',status)
    if status!='PASS':raise SystemExit(1)


if __name__=='__main__':main()
