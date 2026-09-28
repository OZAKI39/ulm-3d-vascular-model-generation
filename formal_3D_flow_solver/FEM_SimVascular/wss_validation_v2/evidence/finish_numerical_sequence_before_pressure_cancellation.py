"""Local ordered gates: finish two-mesh stage 3, then execute medium-grid stage 4.

Never substitutes a queued job, interpolation, or solver exit for CFD evidence.
Run alongside the single stage-three result watcher. Resume from recorded state.
"""
from pathlib import Path
import csv, hashlib, json, os, subprocess, sys, time, traceback

V = Path(__file__).resolve().parents[1]
R = '/workspace/wss_validation_v2_20260927T1230Z'
env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
state = V/'logs/ordered_sequence_progress.json'

def save(**kw):
    state.write_text(json.dumps(dict(observed_unix=time.time(), **kw), indent=2)+'\n')

def read(p):
    return json.loads(Path(p).read_text())

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def run(args, log):
    with (V/'logs'/log).open('w') as f:
        subprocess.run(args, cwd=V, env=env, stdout=f, stderr=subprocess.STDOUT, check=True)

def script(name, extra=(), log=None):
    run([sys.executable, '-B', str(V/'scripts'/name), *map(str, extra)], log or name+'.log')

def main():
    if not (V/'stage3/completion.json').exists():
        save(stage=3, status='WAITING_INDEPENDENT_BASELINE_AND_MEDIUM_RESULTS')
        while True:
            p = V/'logs/stage3_analysis_watch.json'
            d = read(p) if p.exists() else {}
            if 'FAILURE' in d.get('status', ''):
                raise RuntimeError('Stage-three watcher requires review: '+str(d))
            if d.get('status') == 'COMPLETE':
                break
            time.sleep(30)
        for name in ['vessel_baseline','vessel_medium']:
            case = V/'stage3'/name
            assert read(case/'reports/flow_quality.json')['accepted_final_and_log_checks']
            if name != 'vessel_baseline':
                ex = read(case/'reports/execution.json')
                assert ex['status']=='PASS' and ex['input_hashes_unchanged']
                assert read(case/'reports/linear_attempt_acceptance.json')['accepted']
        script('collect_results.py', log='ordered_collect_stage3.log')
        script('collect_geometry.py', log='ordered_geometry_stage3.log')
        script('reproduce_local_wall_evidence.py', log='ordered_local_subset_reproduction.log')
        comparisons = list(csv.DictReader((V/'data/vessel_mesh_comparison.csv').open()))
        targets = [r for r in comparisons if r['metric'] in ['mean_Pa','p05_Pa','p95_Pa']]
        assert len(targets)==33
        result = dict(status='ACTUAL_TWO_MESH_SENSITIVITY_ANALYSIS_COMPLETE',
                      completed_unix=time.time(), cases=['vessel_baseline','vessel_medium'],
                      reuse='Original H0 identity-checked; medium is actual new CFD; fine SKIPPED_BY_USER',
                      regional_metrics_compared=len(targets),
                      fine_status='SKIPPED_BY_USER',
                      reason='User cancelled fine vascular CFD due to time cost; three-grid criteria cannot be evaluated',
                      not_a_claim_of_mesh_independence=True,
                      comparison_csv_sha256=sha(V/'data/vessel_mesh_comparison.csv'))
        (V/'stage3/completion.json').write_text(json.dumps(result, indent=2)+'\n')
        script('make_figures.py', ['--part','vessel'], 'ordered_vessel_figures.log')

    # User amended scope: only the actual accepted medium mesh is the pressure baseline.
    if not (V/'stage4/boundary_group.json').exists():
        save(stage=4, status='REGISTERING_FIXED_GRID_PRESSURE_GROUP')
        baseline = V/'stage3/vessel_medium'
        assert read(baseline/'reports/execution.json')['status']=='PASS'
        assert read(baseline/'reports/flow_quality.json')['accepted_final_and_log_checks']
        s3 = read(V/'stage3/completion.json')
        prereg = f'''# 阶段四：中档固定网格压力实验运行前登记

登记时间（Unix）：{time.time()}。中档实际CFD和独立检查已完成；阶段三记录SHA256 `{sha(V/'stage3/completion.json')}`。
用户已取消细档真实血管CFD，阶段三仅为两档网格敏感性分析，不能证明网格无关或确定离散误差。

依据用户最新指令，固定中档网格（159595节点、884975四面体），复用`vessel_medium`实际合格基线。
三例均为dt=1.5040600916882215e-7s、同一CPU8/LU二进制/PETSc/ASM-overlap2、P1/P1 VMS、零初值、相同入口/材料/网格及原收敛标准。
不采用此前原网格半步长候选基线；此前候选没有生成或运行压力算例。

仅将O2相对O3压差增加或减少1%，delta=29.320152710782013Pa；O1和O3不变。
`prepare_bc_cases.py`检查XML仅O2 Value一处语义变化、所有网格及PETSC_OPTIONS哈希一致。
不新增其他出口/幅度/网格、执行方式或时间步实验。

输出各出口通量/分流，以及固定区域WSS均值/P5/P50/P95/极值、<1/<2/>30Pa面积和位置。
对同一面片分别比较剪切向量和模长的正负响应及偶对称余量。
组内末三个完整保存间隔及终态的实际区域变化作为同网格迭代尺度；原网格dt/2结果仅作背景参照，不能替代中档时间步验证或作为其误差上界。
两档空间变化也不是严格误差上界。若信号与已有不确定性无法区分，如实报告，不通过新增大量计算得出确定结论。
'''
        (V/'stage4/PREREGISTRATION.md').write_text(prereg)
        script('prepare_bc_cases.py', ['--baseline',baseline], 'ordered_prepare_boundary_cases.log')

    names = ['O2_minus1pct','O2_plus1pct']
    if not (V/'stage4/remote_registration.json').exists():
        run(['rsync','-az',str(V/'stage4')+'/', 'vast4090:'+R+'/stage4/'], 'ordered_sync_stage4.log')
        run(['rsync','-az',str(V/'scripts/run_solver.py'),str(V/'scripts/log_acceptance.py'),str(V/'scripts/register_remote_program.py'),'vast4090:'+R+'/scripts/'], 'ordered_sync_runner.log')
        remote = '''from pathlib import Path
import subprocess
r=Path(%r)
for name,wait in [('O2_minus1pct',None),('O2_plus1pct','stage4/O2_minus1pct')]:
 program='boundary_'+name
 if '[program:'+program+']' in (r/'supervisord.conf').read_text(): continue
 cmd=['/root/particle8_2_runs/env/bin/python','-B',str(r/'scripts/register_remote_program.py'),'--name',program,'--case','stage4/'+name]
 if wait:cmd+=['--wait-for',wait]
 subprocess.run(cmd,check=True)
''' % R
        with (V/'logs/ordered_register_stage4.log').open('w') as f:
            subprocess.run(['ssh','vast4090','python3 -'],input=remote,text=True,stdout=f,stderr=subprocess.STDOUT,check=True)
        (V/'stage4/remote_registration.json').write_text(json.dumps(dict(unix=time.time(),cases=names,remote_root=R),indent=2)+'\n')
    for name in names:
        case=V/'stage4'/name
        if (case/'reports/independent_acceptance.json').exists():continue
        save(stage=4, status='WAITING_ACTUAL_CFD', case=name)
        code='from pathlib import Path; import json; p=Path(%r); print(p.read_text() if p.exists() else "{}")' % (R+'/stage4/'+name+'/reports/execution.json')
        while True:
            r=subprocess.run(['ssh','vast4090','python3 -'],input=code,text=True,capture_output=True)
            if r.returncode:time.sleep(30);continue
            ex=json.loads(r.stdout)
            if ex:
                assert ex['status']=='PASS', 'Actual stage-four solver failed: '+name
                break
            time.sleep(30)
        save(stage=4, status='RETRIEVING_AND_ANALYZING', case=name)
        run(['rsync','-az','vast4090:'+R+'/stage4/'+name+'/',str(case)+'/'],'ordered_retrieve_'+name+'.log')
        for script_name,extra in [('analyze_vessel.py',['--case',case]),('mesh_metrics.py',[case]),('diagnose_wall_stencil.py',['--case',case]),('control_volume_checks.py',['--case',case]),('export_local_wall_evidence.py',['--case',case])]:
            script(script_name,extra,'ordered_'+script_name+'_'+name+'.log')
        assert read(case/'reports/flow_quality.json')['accepted_final_and_log_checks']
        (case/'reports/independent_acceptance.json').write_text(json.dumps(dict(accepted=True,unix=time.time(),source_flow_sha256=sha(case/'frozen_flow/flow_arrays_si.npz')),indent=2)+'\n')
    script('collect_results.py', log='ordered_collect_all.log')
    script('analyze_bc_fields.py', log='ordered_boundary_field_analysis.log')
    script('analyze_boundary_resolution.py', log='ordered_boundary_resolution.log')
    script('make_figures.py',['--part','boundary'],'ordered_boundary_figure.log')
    (V/'stage4/completion.json').write_text(json.dumps(dict(status='ACTUAL_MEDIUM_FIXED_GRID_O2_PERTURBATION_CFD_COMPLETE',completed_unix=time.time(),baseline_case=read(V/'stage4/boundary_group.json')['baseline_case'],cases=names,pressure_origin_audit='separate_required_work_not_yet_claimed_complete'),indent=2)+'\n')
    save(stage=4,status='NUMERICAL_SEQUENCE_COMPLETE_PRESSURE_ORIGIN_AND_FINAL_REPORT_REMAIN')

if __name__=='__main__':
    try:main()
    except BaseException as e:
        save(status='REQUIRES_REVIEW',error=repr(e),traceback=traceback.format_exc())
        raise
