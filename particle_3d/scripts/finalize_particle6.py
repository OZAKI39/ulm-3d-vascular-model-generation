#!/usr/bin/env python3
"""Verify saved evidence, bind exact tested sources and write the Chinese review."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys,xml.etree.ElementTree as ET
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parent
sys.path.insert(0,str(PACKAGE/'src'))
from particle_3d.audit import sha256
from particle_3d.particle3_cases import write_json
from particle_3d.particle6_checkpoint import frozen_provenance,source_identity,DEPENDENCIES,MODEL_FLAGS
from particle_3d.lammps_state import schema_contract
from particle_3d.lammps_bridge import validate_command
REPORT=PACKAGE/'reports/particle6';DATA=REPORT/'data';LOGS=REPORT/'logs'
def read(name):return json.loads((DATA/(name+'.json')).read_text())
def git(*args):return subprocess.check_output(['git',*args],cwd=REPO,text=True).strip()
def maximum(rows,key):return max(r[key] for r in rows)
def hashes(paths):return {str(p.relative_to(REPO)):sha256(p) for p in sorted(paths) if p.is_file()}

def test_results(prefix):
    result={}
    for name in ['frozen']+[f'particle{i}' for i in range(7)]:
        path=LOGS/(prefix+'_'+name+'.xml')
        if not path.exists():result[name]=dict(status='PENDING');continue
        root=ET.parse(path).getroot();suites=list(root.iter('testsuite'))
        counts={k:sum(int(s.get(k,0)) for s in suites) for k in ['tests','failures','errors','skipped']}
        result[name]=dict(counts,status='PASS' if counts['tests']>0 and counts['failures']==counts['errors']==counts['skipped']==0 else 'FAIL',junit=str(path.relative_to(REPO)),sha256=sha256(path))
    return result

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--require-pass',action='store_true');parser.add_argument('--bind-source',action='store_true');args=parser.parse_args()
    env=read('00_lammps_environment');state=read('01_state_roundtrip');neighbor=read('02_neighbor_equivalence');one=read('03_one_step')
    multi=[read('04_sphere_multistep'),read('04_mixed_multistep')];resistance=read('05_resistance');force=read('06_force_audit');restart=[read('07_restart_sphere'),read('07_restart_mixed')];real=read('09_real_two_mb');stress=read('10_neighbor_rebuild')
    provenance=frozen_provenance(REPO/'formal_3D_flow_solver/FEM_SimVascular');identity=source_identity(REPO)
    baseline=test_results('baseline');final=test_results('final')
    for name,count in zip(['frozen']+[f'particle{i}' for i in range(6)],[18,46,67,74,90,69,63]):
        assert baseline[name]['status']=='PASS' and baseline[name]['tests']==count
        if args.require_pass:assert final[name]['status']=='PASS' and final[name]['tests']==count
    for command in read('06_command_audit'):validate_command(command['command'])
    assert state['errors']['exact_equal'] and neighbor['false_negatives']==neighbor['mismatch_count']==neighbor['filtered_mismatch']==0
    assert all(one[k]['exact_equal'] for k in ['sphere','mixed'])
    assert all(e['exact_equal'] and e['neighbor_mismatch_count']==0 for case in multi+restart+[real] for e in case['errors'])
    assert all(x==0 for x in resistance['max_errors'].values()) and resistance['active_constraints_equal'] and resistance['sparsity_equal']
    assert all(row['max_lammps_force']==row['max_lammps_torque']==row['lammps_step']==0 and not row['used_as_physics_input'] for row in force)
    assert all(d['physical_time_equal'] and d['step_index_equal'] and d['metadata']['errors']['exact_equal'] for d in restart)
    assert real['original_p5_state_max_error']==real['wall_gap_max_error_m']==real['pair_gap_max_error_m']==real['residual_max_error']==0
    assert real['nearfield_activation_equal'] and real['outlet_events_equal']
    assert all(r['mismatch_count']==0 for r in stress['rows'])
    figures=json.loads((REPORT/'FIGURE_MANIFEST.json').read_text());assert len(figures)==11
    for row in figures:
        assert sha256(REPORT/'figures'/row['figure'])==row['sha256']
        for path,digest in row['sources'].items():assert sha256(REPO/path)==digest
    # Assert the original P5 tested-source snapshot remains unchanged.
    old=json.loads((PACKAGE/'reports/particle5/PARTICLE5_VALIDATION.json').read_text())
    for path,digest in old['source_sha256'].items():assert sha256(REPO/path)==digest,path
    if args.require_pass:
        assert all(v['status']=='PASS' for v in final.values())
        runs=json.loads((LOGS/'final_commands.json').read_text());assert all(r['returncode']==0 for r in runs)
        tested=json.loads((LOGS/'final_tested_sources.json').read_text());assert tested['source_sha256']==identity['source_sha256'],'Source changed since final tests'
    if args.bind_source:
        assert args.require_pass,'Bind source only after all final tests'
        for path,digest in identity['source_sha256'].items():
            blob=subprocess.check_output(['git','show',identity['source_commit']+':'+path],cwd=REPO)
            assert hashlib.sha256(blob).hexdigest()==digest,path
        identity['source_commit_status']='EXACT_TESTED_SOURCE_COMMIT'
    for folder in [REPORT/'checkpoints/sphere',REPORT/'checkpoints/mixed']:
        manifest=json.loads((folder/'checkpoint_manifest.json').read_text())
        for path,digest in manifest['files'].items():assert sha256(folder/path)==digest
        assert manifest['source_sha256']==identity['source_sha256'],'Checkpoint source snapshot differs; regenerate using final source'
        if args.bind_source:
            manifest.update(source_commit=identity['source_commit'],source_commit_status=identity['source_commit_status']);write_json(folder/'checkpoint_manifest.json',manifest)
    errors=[e for case in multi for e in case['errors']];re=[e for case in restart for e in case['errors']]
    status='PASS' if args.require_pass else 'PENDING_FINAL_REGRESSION'
    v=dict(stage='Particle-6',automated_checks=status,stage_result='AUTOMATED_PASS_PENDING_USER_REVIEW' if args.require_pass else status,
        git_branch=git('branch','--show-current'),git_commit=identity['source_commit'],**{f'particle{i}_dependency_commit':c for i,c in enumerate(DEPENDENCIES)},
        source_commit_status=identity['source_commit_status'],source_sha256=identity['source_sha256'],
        frozen_integrity=dict(status='PASS',**provenance),historical_source_integrity='PASS; P5 tested source snapshot unchanged',
        **{k:env[k] for k in ['lammps_version','lammps_python_module','lammps_library_path','lammps_mpi_available']},
        **MODEL_FLAGS,pair_style='zero',custom_property_schema=schema_contract(),state_roundtrip_status='PASS',state_roundtrip_max_errors=state['errors'],
        neighbor_raw_pair_count=len(neighbor['lammps_candidates']),neighbor_standalone_pair_count=len(neighbor['standalone_candidates']),
        neighbor_filtered_pair_count=len(neighbor['lammps_exact_pairs']),neighbor_false_negatives=neighbor['false_negatives'],neighbor_mismatch_count=neighbor['mismatch_count'],neighbor_extra_candidates=neighbor['extra_candidates'],
        one_step_max_position_error_m=one['sphere']['position_m'],one_step_max_velocity_error_m_s=one['sphere']['velocity_m_s'],one_step_max_omega_error_s_inv=one['sphere']['omega_s_inv'],one_step_max_orientation_error_rad=one['sphere']['orientation_rad'],
        multistep_steps=100,multistep_case_count=2,multistep_max_position_error_m=maximum(errors,'position_m'),multistep_max_orientation_error_rad=maximum(errors,'orientation_rad'),multistep_max_velocity_error_m_s=maximum(errors,'velocity_m_s'),multistep_max_omega_error_s_inv=maximum(errors,'omega_s_inv'),multistep_neighbor_mismatch_count=sum(e['neighbor_mismatch_count'] for e in errors),
        resistance_matrix_max_error=resistance['max_errors']['R'],rhs_max_error=resistance['max_errors']['b'],solved_velocity_max_error=resistance['max_errors']['U'],contact_jacobian_max_error=resistance['max_errors']['J'],contact_active_constraints_equal=resistance['active_constraints_equal'],
        max_lammps_force=maximum(force,'max_lammps_force'),max_lammps_torque=maximum(force,'max_lammps_torque'),force_audit_count=len(force),
        torque_interpretation='UNDEFINED_FOR_ATOM_STYLE_ATOMIC; DIAGNOSTIC_ZERO_NOT_PHYSICS_INPUT',
        checkpoint_restart_status='PASS',checkpoint_step=40,restart_final_position_error_m=max(d['final_errors']['position_m'] for d in restart),restart_final_orientation_error_rad=max(d['final_errors']['orientation_rad'] for d in restart),restart_neighbor_mismatch_count=sum(e['neighbor_mismatch_count'] for e in re),restart_physical_time_preserved=True,restart_step_index_preserved=True,
        real_two_mb_bridge_status=real['status'],real_two_mb_max_position_error_m=maximum(real['errors'],'position_m'),real_two_mb_max_velocity_error_m_s=maximum(real['errors'],'velocity_m_s'),real_two_mb_max_omega_error_s_inv=maximum(real['errors'],'omega_s_inv'),real_two_mb_pair_gap_error_m=real['pair_gap_max_error_m'],real_two_mb_wall_gap_error_m=real['wall_gap_max_error_m'],real_two_mb_resistance_residual_error=real['residual_max_error'],real_two_mb_nearfield_activation_equal=real['nearfield_activation_equal'],real_two_mb_outlet_events_equal=real['outlet_events_equal'],
        parallel_lammps_validated=False,manual_visual_review='PENDING_USER_REVIEW',blocking_software_issue='NONE' if args.require_pass else 'FINAL_TESTS_PENDING',
        tests=dict(baseline=baseline,final=final),data_sha256=hashes(DATA.iterdir()),figure_sha256=hashes((REPORT/'figures').glob('*.png')),
        checkpoint_sha256=hashes((REPORT/'checkpoints').rglob('*')),log_sha256=hashes(LOGS.iterdir()),
        roundoff_comparison_policy='EXACT_EQUALITY_ACHIEVED; optional budget = 512*eps*steps*field_scale; orientation 512*eps*steps; no relaxed rtol',
        real_rbc_passage='NOT_ESTABLISHED',full_suspension=False)
    write_json(REPORT/'PARTICLE6_VALIDATION.json',v)
    sections=['# Particle-6 人工审核报告',
        '## 1. 这一阶段做了什么',
        f"把已经验证的 particle solver 接进实际 LAMMPS。没有修改 P0–P5 数值源码。P5 人工审核已按本轮授权单独记录为 PASS，并保留全部上游限制。\n\n分支：`{v['git_branch']}`。本轮源码提交：`{v['git_commit']}`；绑定状态：`{v['source_commit_status']}`。\n\n结果来自本轮实际运行；源码、数据、图、checkpoint 和日志的 SHA256 全部保存在 [机器验证记录](PARTICLE6_VALIDATION.json)。",
        '## 2. LAMMPS 到底负责什么',
        f"保存稳定 ID、类型、位置与自定义状态，生成邻居候选，输出 dump 和二进制 restart。版本 {env['lammps_version']}，Python {env['python_version'].split()[0]}，WSL CPU 单 MPI rank。\n\nPython 模块：`{env['lammps_python_module']}`。共享库：`{env['lammps_library_path']}`。LAMMPS 和 MPICH 仅安装在项目 venv；已有 FEM 数值依赖只读复用，未修改 Frozen FEM venv。安装审计见 [环境记录](data/00_lammps_environment.json)。",
        '## 3. LAMMPS 不负责什么',
        f"不算 F=ma，不做加速度积分，不计算 lubrication 或 contact physics，不重新决定 V / Omega。只用 pair_style zero、property/atom 和零步维护；所有 LAMMPS 时间步均为 0。\n\n{len(force)} 条重建前后审计中，force 最大值为 0。atomic 模式没有 torque 数组，诊断记为 0 并标注未定义、未使用；不能解读成物理力矩测量。mass = 1 只是 LAMMPS_BRIDGE_PLACEHOLDER_ONLY。永久测试还主动污染 force，确认在 run 0 清零之前被拒绝。",
        '## 4. particle state 怎样映射',
        'ID 和位置使用 LAMMPS 原生存储。类型编码 MB=1、RBC=2；形状编码 SPHERE_MB=1、FREE_OBLATE=2、CAPILLARY_DEFORMED=3。不可行变形状态拒绝插入。其余 q、V、Omega、半径、半轴、胶囊轴/半径/柱段长度和包围半径使用 custom properties。完整字段见 [V1 契约](../../contracts/particle6_state_v1.json)。\n\n额外保存原 P4 椭球旋转矩阵，避免从 q 重建产生舍入变化。胶囊几何由轴、半径、长度决定：RETAINED_QUATERNION_NOT_CAPSULE_ORIENTATION。六个固定粒子全部状态往返最大差均为 0。',
        '## 5. neighbor list 怎样比较',
        f"先使用相同 ValidationNeighborPolicy，再将 LAMMPS half list 的本地索引转换成稳定 ID，去重排序；原 P4/P5 精确几何与近场资格决定哪些候选进入物理计算。13 个粒子的 78 个可能配对中，standalone 和 LAMMPS 各产生 {v['neighbor_raw_pair_count']} 对候选，精确筛选后 {v['neighbor_filtered_pair_count']} 对；额外候选 {v['neighbor_extra_candidates']}，漏检 0，不匹配 0。六种形状组合均覆盖，ID 非连续且输入打乱。\n\n每次接受物理子步后写回并强制重建。跨范围压力测试证明旧候选会失效，而重建后逐步一致。所有 center cutoff 和 skin 都是 VALIDATION_NEIGHBOR_QUERY_ONLY，不是生产物理截断。",
        '## 6. 一个 timestep 是否完全一致',
        'P5 球体靠墙与成对近场场景在 dt = 0.001 s 后，x、q、V、Omega、类型与邻居完全一致，最大误差为 0。相同状态组装的 R、b、U、J 最大差均为 0，稀疏结构、接触 ID 和活动集相同。\n\n桥接只向原 trial 的私有依赖副本提供候选；使用完全相同 code object、原 P3 时间细分及原 P2 四元数积分。没有模块全局 monkeypatch，没有宽松 rtol。每步保存 float64 舍入预算，但本次实际达到逐位一致。',
        '## 7. 多步是否漂移',
        'P5 球体与 P4 混合形状各独立推进 100 个验证步，未用 reference 覆盖桥接状态。位置、q 轴、V、Omega、几何字段最大差均为 0，邻居不匹配总数为 0。混合形状仅复用 P4/P2 已有运动学；它不代表 P5 已具备 RBC 流体动力学。',
        '## 8. checkpoint / restart 是否连续',
        '球体和混合形状均实际运行 40 步，写 state.restart 与全局 sidecar，销毁旧 LAMMPS 实例，再在新实例中 read_restart，并继续 60 步。最终所有状态、邻居、physical_time_s = 0.1 s 和 step_index = 100 与连续路径完全一致。\n\n三种形状的全部 custom properties 从二进制恢复；sidecar 不保存逐粒子数组，也不重新采样。本例无运行时随机过程，rng_states = null。清单记录 binary/sidecar SHA、schema、LAMMPS 版本和源码提交；sidecar 记录 31 个 Frozen 文件哈希、依赖提交、黏度与模型/查询策略。两个交付 bundle 见 [混合形状](checkpoints/mixed/checkpoint_manifest.json) 和 [球体](checkpoints/sphere/checkpoint_manifest.json)。',
        '## 9. 真实 two-MB 是否一致',
        f"原两个 SonoVue 半径为 {real['bridge'][0]['particles'][0]['radius']:.17g} m、{real['bridge'][0]['particles'][1]['radius']:.17g} m。原初始位置未动，dt = {real['dt_s']:.17g} s，horizon = {real['horizon_s']:.17g} s，{real['requested_steps']} 个请求步、{real['accepted_substeps']} 个接受子步。\n\n本轮重新运行完整原 P5 simulate_resistance，与桥接逐接受子步比较：x、V、Omega、墙隙、pair gap、阻力残差最大差均为 0，近场启用和出口事件一致。此窗口墙面近场启用，pair 润滑未启用，没有出口事件。原计算出现的亚纳米墙隙没有被人为截断；其连续介质有效性仍未建立。真实 RBC 动力学记录为 DEFERRED_DUE_TO_UPSTREAM_PHYSICS。",
        '## 10. 自动测试',
        '| 测试组 | 开发前 | 最终回归 |\n|---|---:|---:|\n'+'\n'.join(f"| {name} | {baseline[name].get('tests','—')} / {baseline[name]['status']} | {final[name].get('tests','—')} / {final[name]['status']} |" for name in final),
        '原历史测试保持原内容。P6 永久测试实际调用 LAMMPS，包括 100 步、真实 FEM 和销毁重启；图测试重新生成 PNG 并比较 SHA256。初始/最终完整命令、返回码和 JUnit 都在 [logs](logs)。复现接口与环境说明见 [P6 README](../../PARTICLE6_README.md)。',
        '## 11. 人工审核图']
    for i,row in enumerate(figures):
        sections.extend([f"**{i:02d} — {row['figure']}**",f"![{i:02d}](figures/{row['figure']})",(REPORT/f'{i:02d}_step_notes.md').read_text().strip()])
    sections.extend(['## 12. 当前限制',
        '- no LAMMPS force physics；\n- no LAMMPS dynamics integrator；\n- no production neighbor cutoff；\n- no production neighbor skin；\n- no production lubrication cutoff；\n- no production timestep；\n- non-spherical lubrication not frozen；\n- real RBC passage not established；\n- real RBC dynamics bridge deferred；\n- serial correctness only，本轮没有 MPI 多 rank smoke；\n- no full suspension。\n\n亚纳米间隙连续介质有效性 NOT_ESTABLISHED。固定验证 box 越界会报错；私有依赖适配只支持锁定 P0–P5 版本。没有执行 CFD，没有开始 Particle-7，也没有 push 或 merge main。',
        '## 13. 人工审核状态',f"AUTOMATED_CHECKS = {status}\n\nMANUAL_VISUAL_REVIEW = PENDING_USER_REVIEW\n\n本报告提供自动证据和可审核图；没有代替用户签署 Particle-6 人工审核。"])
    (REPORT/'PARTICLE6_REVIEW.md').write_text('\n\n'.join(sections)+'\n')
    print(json.dumps(dict(automated_checks=status,git_commit=v['git_commit'],tests=final),indent=2))

if __name__=='__main__':main()
