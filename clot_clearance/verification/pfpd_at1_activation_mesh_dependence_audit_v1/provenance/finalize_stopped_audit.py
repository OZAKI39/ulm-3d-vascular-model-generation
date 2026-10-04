"""Finalize the paper-gate stop. Standard library only; no scientific solver."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

OUT = Path(__file__).resolve().parents[1]
ROOT = OUT.parents[1]
RESULTS = ROOT / 'results/pfpd_at1_activation_mesh_dependence_audit_v1'


def read(path):
    return json.loads(path.read_text())


def save(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def digest(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def tests(name):
    tree = ET.parse(OUT / f'tests/{name}.xml')
    cases = tree.getroot().findall('.//testcase')
    failures = [dict(id=c.get('classname') + '::' + c.get('name'),
                     message=c.find('failure').get('message'),
                     traceback=c.find('failure').text)
                for c in cases if c.find('failure') is not None]
    return dict(passed=sum(all(c.find(t) is None for t in ['failure', 'error', 'skipped']) for c in cases),
                failed=len(failures), errors=sum(c.find('error') is not None for c in cases),
                skipped=sum(c.find('skipped') is not None for c in cases), failures=failures)


def preserve():
    pre = read(OUT / 'PRE_RUN_MANIFEST.json')
    changed = []; missing = []; links = []
    for name, value in pre['files'].items():
        path = ROOT / name
        if not path.is_file():
            missing.append(name)
        elif digest(path) != value:
            changed.append(name)
    for name, value in pre['symlinks'].items():
        path = ROOT / name
        if not path.is_symlink() or os.readlink(path) != value:
            links.append(name)
    result = dict(files=len(pre['files']), symlinks=len(pre['symlinks']), modified=changed,
                  missing=missing, symlink_changes=links,
                  production_modified=[n for n in changed if n.startswith('pd_clot/')],
                  existing_tests_modified=[n for n in changed if n.startswith('tests/')],
                  passed=not (changed or missing or links))
    save(OUT / 'PRESERVATION_CHECK.json', result)
    assert result['passed'], result
    return result


def review():
    baseline = tests('BASELINE'); final = tests('FINAL')
    assert baseline['failures'] == final['failures']
    assert [baseline[k] for k in ['passed','failed','errors','skipped']] == [412,3,0,0]
    assert [final[k] for k in ['passed','failed','errors','skipped']] == [412,3,0,0]
    preservation = preserve()
    history = read(OUT / 'HISTORICAL_ACTIVATION_FACTS.json')
    resources = read(OUT / 'tests/FINAL_RESOURCES.json')
    plan = read(OUT / 'RESOURCE_PLAN.json')
    head = subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip()
    assert head == (OUT/'provenance/HEAD.txt').read_text().strip()
    status = dict(PFPD_PAPER_PHASE_MAPPING='AMBIGUOUS',DIRECT_PHASE_REGULARIZATION='ABSENT_BY_THEORY',
                  AFFINE_ORACLE='INCONCLUSIVE',EQUILIBRIUM_DRIVER_HOTSPOT='INCONCLUSIVE',
                  ACTIVATION_LAYER='FIRST_GRID_LAYER',PROFILE_SCALING='INCONCLUSIVE',
                  SUPPORT_TRUNCATION_ASSOCIATION='INCONCLUSIVE',CLAMP_KINEMATIC_ASSOCIATION='NOT_TESTED',
                  PHYSICAL_DAMAGE_SUPPORT_C5_C7='INCONCLUSIVE',ACTIVATION_THRESHOLD_OBJECTIVITY='NOT_ESTABLISHED',
                  ROOT_CAUSE_PRIMARY='INCONCLUSIVE',C7_STABILITY='NOT_RUN',C9='NOT_RUN',PRODUCTION='UNCHANGED')
    result = dict(task='PFPD_AT1_ACTIVATION_MESH_DEPENDENCE_AUDIT_V1',audit_status='STOPPED_AT_PAPER_MAPPING_GATE',
                  audit_complete=False,status=status,hard_stop=read(OUT/'HARD_STOP.json'),
                  baseline_tests=baseline,final_tests=final,known_failures_unchanged=True,
                  new_tests=0,new_test_reason='Stopped before numerical diagnosis; no scientific implementation changed.',
                  preservation=preservation,historical_facts=history,resources=plan,final_test_resources=resources,
                  observed_wall_seconds_since_pre_manifest=time.time()-read(OUT/'PRE_RUN_MANIFEST.json')['timestamp'],
                  wall_definition='Elapsed wall time from PRE_RUN_MANIFEST timestamp to report finalization, including reading, tools and two full test suites; excludes initial manifest hashing.',
                  peak_RAM_whole_task_bytes=None,peak_RAM_note='No whole-task peak monitor; final-suite observed peak is separately reported.',
                  new_equilibrium_solves=0,new_affine_oracle_evaluations=0,new_clamp_control_solves=0,
                  numerical_diagnosis='NOT_RUN_HARD_STOP',figures_generated=0,HEAD=head)
    save(OUT/'SUMMARY.json',result)
    lines = [
        '# PFPD_AT1_ACTIVATION_MESH_DEPENDENCE_AUDIT_V1',
        '**状态：STOPPED_AT_PAPER_MAPPING_GATE；未达到 AUDIT COMPLETE。** 本轮完成原文—代码映射并触发停止门槛，未开展新的网格、affine 或 clamp 数值实验。未发现应有而缺失的直接相场正则项，根因仍为 INCONCLUSIVE。',
        '## 1. 为什么现在不做 C7 stability',
        '当前问题是激活阈值的空间采样依赖，而非 C7 稳定性。本轮还必须先通过论文等价性门槛；该门槛未通过，所以也未进入共同应变力学诊断。C7 Hessian、Schur、critical cone 和 C9 均未运行。',
        '## 2. “2/5392 与 2/46059”为什么不能直接叫开裂比例',
        '它们是键计数比例，分母随空间采样改变，不是物理裂纹面积或体积分数。历史数值保留；本轮没有将它们转成开裂比例，也没有用两个分离热点的整体 RMS 解释单个损伤区宽度。',
        '## 3. q=2(Yc−Psi+) 告诉了我们什么',
        '在当前候选、s=0、h=1 下，q_b=2(Yc−Ψ⁺_b)。正的键积分权重及正的残差缩放都不改变它的符号。阈值变化需要检查真实机械平衡下的局部 Ψ⁺，不能单由 PETSc 选择解释。本轮没有继续计算共同应变场，因此尚无新的集中因素归因。',
        '## 4. 原 PFPD 论文的 phase spatial regularization 到底是什么',
        '本轮重新读取本地 PDF 第 6–9、14 页。式 (3.2) 是经典相场，含显式梯度项；式 (3.3) 后明确将该项从 PFPD 形式移除。式 (3.4)–(3.8) 的裂纹能对每根键局部可分。论文将正则化归于 horizon 内的非局部能量结构，而不是显式邻键相场梯度块。DIRECT_PHASE_REGULARIZATION=ABSENT_BY_THEORY（限定 h=1 的当前域）。',
        '固定 x 时当前 Hss=diag(2W_bΨ⁺_b)，这本身不代表漏项。核、c0 和 Gc/(c0δ) 均可找到对应实现。无限域完全断裂平面的归一化不能自动证明当前有限小试件中的局部损伤带宽度客观。详见 [逐项映射](references/PFPD_PHASE_FORMULATION_MAPPING.md) 和 PHASE_COUPLING_STRUCTURE.json；本轮没有生成 C7 Hessian。',
        '## 5. 当前代码与论文是否一致',
        '裂纹能和低损伤核归一化可以对应，但整体映射为 AMBIGUOUS。打印式 (3.4)/(3.10) 为 gΨ0 与 gP0；当前代码使用 gΨ⁺+Ψ⁻ 及对应应力。在 s>0、Ψ⁻>0 时，能量差为 [1−g(s)]Ψ⁻。原文式 (3.19) 指定拉伸驱动，但没有在上述总能量和应力公式中明确当前采用的替换。',
        '此外，式 (3.17) 的 AT1 约束最小值是 max(0,1−Yc/Y)，式 (3.25) 却给出 Yhistory/(Yhistory+Yc)。Y=Yc 时分别为 0 和 1/2；后者不是该 AT1 能量的内部驻值。这是此前已经记录并由用户选择 AT1 候选绕开的原文歧义，本轮重新读原文确认。历史候选验证仍然有效；不能由此声称其为全文精确复现，也不能把此歧义直接当成 C5/C7 偏移的原因。',
        '用户第 8 节要求 VERIFIED/EQUIVALENT 后才继续，第 46A 节要求能量不一致时停止。故本轮暂停后续数值诊断，没有自行放宽门槛或替换用户已选模型。没有找到 MISSING_REQUIRED_PHASE_TERM，未输出 PHASE_FORMULATION_INCOMPLETE。',
        '## 6. exact affine field 是否有 hotspot',
        'NOT_RUN_HARD_STOP。没有执行新的 affine oracle，因此 C3/C5/C7 的此次 Ψ⁺ CV 均为 N/A，不能填入历史不同任务的结果冒充本次共同应变审计。',
        '## 7. true equilibrium 为什么产生/不产生 hotspot',
        'NOT_RUN_HARD_STOP。ε=0.0020、0.0025、0.0030 的跨网格平衡比较及 A_psi/A_E 本轮未执行；未将非平衡仿射场当成机械平衡。',
        '## 8. hotspot 是否跟 d/h 缩放',
        '共同应变 profile 和 pairwise RMSE 未计算。单独两个激活键的 d/h≈0.5 不能替代用户要求的至少两个应变下 RMSE(d/h)<0.5*RMSE(d/δ) 检查。PROFILE_SCALING=INCONCLUSIVE。',
        '## 9. hotspot 是否跟 d/delta 缩放',
        '未执行 profile collapse，因此不支持或否定固定 horizon 尺度的层厚。保持 δ=0.25125 mm 不变。',
        '## 10. first activation 是否一直位于第一网格层',
        '只读取上一轮已保存的键中点与 clamp 距离并换算单位。C5：ε=0.004497390736212，d_act=31.25 μm，h=62.5 μm，d_act/h=0.5，d_act/δ=0.124378109453。C7：ε≈0.003137607455254，d_act=20.833333333 μm，h=41.666666667 μm，d_act/h≈0.5，d_act/δ=0.0829187396352。C7 的激活定位仍以历史 bracket [0.0031375732421875,0.00313763427734375] 为准。',
        'ACTIVATION_TRACKS_FIRST_GRID_LAYER 仅描述这两个历史网格首次激活键的几何位置；没有据此宣布 artifact 或完成根因归因。380 个 C5/C7 历史 NPZ 与上轮封存清单哈希一致。',
        '## 11. family/support truncation 的关联',
        'NOT_RUN_HARD_STOP。未执行 support-driver 散点、相关系数或候选端点与同层粒子比较。δ≈L，任意粒子至所有表面的最小距离都不可能达到 δ，所以此立方体没有真正完整球支持的 bulk；这是几何事实，不是截断导致集中之因果证据。',
        '## 12. clamp DOF 是否贡献集中',
        '源码确认：上下最外层粒子的 x/y/z 均固定；底层位移为零，顶层仅被指定 z=εL 的位移，顶层 x/y 仍固定。自由粒子保留三个平移分量。详情见 CLAMP_DOF_AUDIT.json。',
        'CLAMP_KINEMATIC_ASSOCIATION=NOT_TESTED。条件式控制的前置 affine/equilibrium 检查未运行，故没有修改 clamp，也没有求 axial-only 状态。',
        '## 13. C5/C7 真实物理 damage measure 是否比 bond fraction 更一致',
        'NOT_RUN_HARD_STOP。未挑选或计算本轮 matched max(s) 下的 D_weighted、Ecrack/Gc、上下热点单独尺寸及体积代理；这些项均为 N/A。历史分支文件保持不变。',
        '## 14. 当前 root cause',
        'ROOT_CAUSE_PRIMARY=INCONCLUSIVE。HYPOTHESIS D 中“漏了论文明确要求的邻键相场梯度项”没有得到支持；这不等于证明整体模型正确或完整。A/B/C/E 未完成本轮控制与尺度诊断。已知方程差异不是网格依赖的因果证明，两级旧结果使用的是同一候选。',
        '## 15. 对最终 clot simulation 的影响',
        '现有 C5/C7 结果仍是指定 AT1+SVK 谱分裂候选的离散驻定态证据。激活阈值客观性和空间收敛尚未建立；不得把两个活动键直接解释为真实血栓裂纹或据此进行材料标定。当前任务没有推进到 CFD 或侵蚀耦合。',
        '## 16. tests',
        f"起始全量：{baseline['passed']} passed / {baseline['failed']} known failed / {baseline['skipped']} skipped；最终全量：{final['passed']} passed / {final['failed']} known failed / {final['skipped']} skipped。原 3 项的 ID、message、traceback 逐字相同；没有 xfail、skip 或弱化断言。",
        '本轮没有修改或新增科学实现，且在论文硬门停止；因此后续 profile/affine/support/damage 的新增科学测试未创建、未假称通过。只按任务要求各执行一次起始和结束全量套件。',
        f"保护校验：{preservation['files']} 个原文件和 {preservation['symlinks']} 个符号链接保持一致，production/history/existing-tests modified=0/0/0。PRE/POST 清单另存；当前 HEAD={head}。",
        '## 17. CPU/GPU 使用',
        f"重新检测 {plan['physical_cores']} physical / {plan['logical_cores']} logical cores，RAM={plan['RAM_total_bytes']/1024**3:.3f} GiB。论文/源码审阅单进程，新增力学 worker=0，新 equilibrium solve=0，所有测试 BLAS/OMP/MKL/Numba 单线程。GPU_NOT_USED_BY_DESIGN，未重试安装。",
        f"从 PRE manifest 时间戳到报告生成的经过时间 {result['observed_wall_seconds_since_pre_manifest']:.3f} s（含阅读、工具往返和两次测试；不含最初哈希）。没有全任务 RAM 峰值监测，不能报全程峰值。最终测试独立监测到进程树峰值 {resources['peak_observed_tree_RSS_bytes']/1024**3:.3f} GiB，耗时 {resources['wall_seconds']:.3f} s，监测范围仅此测试过程。",
        '## 18. 唯一下一步',
        '等待人工审核论文映射与当前候选的规范范围，特别是是否在后续审计中认可既有 AT1+谱分裂选择，并将论文核对限定为裂纹能、核归一化与耦合结构。审核前不继续数值归因，不改模型，不运行 C7 stability 或 C9。',
        '本轮图目录已创建，但没有生成要求的后续数值图；FIGURE_MANIFEST.json 明确记为 NOT_GENERATED_HARD_STOP。',
        '```text\n' + '\n'.join(f'{k}: {v}' for k,v in status.items()) + '\n```',
    ]
    (OUT/'REVIEW_ACTIVATION_MESH_DEPENDENCE_V1.md').write_text('\n\n'.join(lines)+'\n')
    (RESULTS/'figures/NOT_GENERATED.txt').write_text('NOT_GENERATED_HARD_STOP\nThe paper/code mapping gate stopped numerical diagnosis. No scientific figure is claimed.\n')
    source_names = ['pd_clot/pfpd_candidate/kernel.py','pd_clot/pfpd_candidate/normalization.py',
                    'pd_clot/pfpd_candidate/at1_evolution.py','pd_clot/pfpd_candidate/baqp.py',
                    'pd_clot/pfpd_candidate/degradation.py','pd_clot/pfpd_candidate/svk_spectral_split.py',
                    'pd_clot/pfpd_candidate/coupling_stability_audit.py','pd_clot/pfpd_candidate/adr_adapter_at1.py',
                    'scripts/c5_activation_audit.py','scripts/c7_early_branch.py',
                    'scripts/spectral_split_audit.py','scripts/c5_stationarity_diagnostics.py']
    save(OUT/'provenance/AUDITED_SOURCE_HASHES.json',{n:digest(ROOT/n) for n in source_names})
    for name,args in [('HEAD_AFTER',['git','rev-parse','HEAD']),('STATUS_AFTER',['git','status','--short']),('DIFF_AFTER',['git','diff','--stat'])]:
        (OUT/f'provenance/{name}.txt').write_text(subprocess.check_output(args,cwd=ROOT,text=True))
    print(json.dumps(dict(audit_status=result['audit_status'],tests=[final['passed'],final['failed'],final['skipped']],
                         historical_preservation=preservation['passed'],known_failures_unchanged=True,
                         elapsed_seconds=result['observed_wall_seconds_since_pre_manifest'],
                         final_test_peak_RSS_bytes=resources['peak_observed_tree_RSS_bytes']),indent=2))


def post():
    files = {}; links = {}
    for directory, dirs, names in os.walk(ROOT,followlinks=False):
        dirs[:]=[name for name in dirs if name!='.git']
        for name in list(dirs)+names:
            path=Path(directory)/name
            key=str(path.relative_to(ROOT))
            if path==OUT/'POST_RUN_MANIFEST.json':
                continue
            if path.is_symlink():
                links[key]=os.readlink(path)
            elif path.is_file():
                files[key]=digest(path)
    save(OUT/'POST_RUN_MANIFEST.json',dict(root=str(ROOT),files=files,symlinks=links,self_excluded=True,timestamp=time.time()))
    pre=read(OUT/'PRE_RUN_MANIFEST.json')
    assert all(files.get(k)==v for k,v in pre['files'].items())
    assert all(links.get(k)==v for k,v in pre['symlinks'].items())
    print('POST_VERIFIED',len(files),'files;',len(links),'symlinks; historical modifications=0')


if __name__=='__main__':
    {'review':review,'post':post}[sys.argv[1]]()
