#!/usr/bin/python3
"""Independent offline replay and source audit. Never launches solver."""
import sys,os,csv,json,subprocess
from pathlib import Path
import numpy as np
sys.dont_write_bytecode=True
from remote_common import *
from convergence import *
from verify_run import verify,equivalent
from gpu_integrity import inspect_snapshot
R=Path(__file__).resolve().parents[1];S=R.parent/'20260914_stage4_batch_dispatch'
def audit():
    verify_bundle(R);verify_run_inputs(R,R);term=read(R/'RUN_TERMINAL.json')
    if term['returncode']!=0:raise ValueError('Solver error terminal: '+term['status'])
    summary=verify(R,R);n=summary['actual_steps']
    safety=np.genfromtxt(R/'runtime_safety.csv',names=True,delimiter=',',dtype=np.int64)
    assert len(safety)==n+1 and np.array_equal(safety['iteration'],np.arange(n+1))
    for col in safety.dtype.names[1:]:assert not np.any(safety[col]),'Nonfinite safety count'
    short=read(R/'diagnostics/SHORT_CASE_IDENTITY.json');assert short['status']=='PASS'
    checks=[]
    for step in range(5000,n+1,5000):
        result=inspect_snapshot(R,R,step)
        assert equivalent(result,read(R/f'diagnostics/gpu_integrity_{step}.json')),'Device monitor replay mismatch'
        checks.append(result)
    zero=inspect_snapshot(R,R,0);f=fields(field_path(R,0));assert np.all(f['rho']==1) and np.all(f['u']==0),'From-zero identity'
    with (R/'final_audit/STAGE4_PROTECTED_AUDIT.log').open('w') as log:
        ret=subprocess.run(['/usr/bin/python3','-B',str(S/'scripts/verify_stage4_archive.py')],stdout=log,stderr=subprocess.STDOUT)
    assert ret.returncode==0,'Protected Stage1-4 audit failed'
    result=dict(status='PASS',SOLVER_FINALIZER_IDENTITY='PASS',runtime_safety='PASS',gpu_numerical_integrity='PASS',source_integrity='PASS',
        full_convergence_replay=summary,every_evaluation_gpu_integrity_checks=checks,initial_snapshot=zero,short_case_identity=short,all_step_safety_rows=len(safety),
        physics_changed=False,numerical_contract_changed=False,geometry_changed=False,Qtarget_changed=False,RBC_included=False,IBM_included=False,human_review='PENDING')
    write(R/'FINAL_AUDIT.json',result);return result
def reports(result=None,error=None):
    term=read(R/'RUN_TERMINAL.json');timing=read(R/'diagnostics/solver_timing.json') if (R/'diagnostics/solver_timing.json').exists() else {}
    freeze=read(R/'FORMAL_STEP3C_GPU_FREEZE.json');c=contract(R);actual=term['actual_steps']
    summary=result['full_convergence_replay'] if result else {};row=summary.get('final_row',{});detail=summary.get('final_evaluation',{})
    state=read(R/'diagnostics/evaluator_state.json') if (R/'diagnostics/evaluator_state.json').exists() else {}
    windows=detail.get('windows',{});long=windows.get('long',{});qin=long.get('Qin')
    confirmed=bool(result and summary.get('auto_converged'));sps=timing.get('full_iteration_steps_per_second')
    reg=False;wp=R/'diagnostics/iteration_windows.csv'
    if wp.exists():
        with wp.open() as f:chunks=list(csv.DictReader(f))
        slow=0;maxslow=0
        for item in chunks:
            slow=slow+1 if float(item['steps_per_second'])<.75*125.8808630750868 else 0;maxslow=max(maxslow,slow)
        reg=maxslow>=10
        write(R/'performance/PERFORMANCE_REGRESSION.json',dict(PERFORMANCE_REGRESSION=reg,definition='10 consecutive 1000-step intervals below 75% of Stage4 baseline; includes output/audit costs; observation only',maximum_consecutive_slow_intervals=maxslow))
    vals=dict(STEP3C_FORMAL_GPU_STATUS=term['status'],PRODUCTION_GPU_CANDIDATE='STAGE4',GPU_MODEL=freeze['gpu_model'],GPU_COMPILER=freeze['compiler_version'].strip().replace('\n','; '),PALABOS_GPU_COMMIT=freeze['palabos_gpu_commit'],
      FORMAL_RUN_FROM_ZERO='YES',CALIBRATED_MULTIPLIER=1.1197286861799598,Q_TARGET=c['Qtarget_m3_s'],
      FORMAL_RUN_REPRODUCES_VALIDATED_SHORT_CASE='PASS' if (R/'diagnostics/SHORT_CASE_IDENTITY.json').exists() else 'FAIL_OR_NOT_REACHED',
      ACTUAL_STEPS=actual,FINAL_PHYSICAL_TIME_S=actual*c['dt_s'],ACTUAL_WALL_TIME_SECONDS=term['end_to_end_seconds'],AVERAGE_GPU_STEPS_PER_SEC=sps,
      FIRST_CONVERGED_ITERATION=state.get('first_converged_iteration'),CONFIRMATION_ITERATION=state.get('confirmation_iteration'),STOP_ITERATION=actual,
      FINAL_MEAN_QIN=qin,FINAL_QIN_OVER_QTARGET=qin/c['Qtarget_m3_s'] if qin is not None else None,FINAL_R_INLET=row.get('R_inlet'))
    for i,v in enumerate(long.get('Qout',[None]*3),1):vals[f'FINAL_MEAN_QOUTLET_0{i}']=v
    for name,key in [('R_FLOW','R_flow'),('R_CV','R_CV'),('R_VELOCITY','R_velocity'),('R_PRESSURE','R_pressure'),('FLOW_FRACTION_DRIFT','flow_fraction_drift'),('MULTIPLANE_SPREAD','cross_plane_flux_spread'),('QUADRATURE_ERROR','quadrature_error')]:vals['FINAL_'+name]=row.get(key)
    vals.update(ALL_OUTLET_MEANS_POSITIVE=all(v>0 for w in windows.values() for v in w['Qout']) if windows else 'UNVERIFIED',
      RUNTIME_SAFETY='PASS' if result else 'FAIL_OR_UNVERIFIED',GPU_NUMERICAL_INTEGRITY='PASS' if result else 'FAIL_OR_UNVERIFIED',SOLVER_FINALIZER_IDENTITY='PASS' if result else 'FAIL',
      STEP3B_REFERENCE_R_INLET=0.1069,STEP3C_FINAL_R_INLET=row.get('R_inlet'),INLET_CALIBRATION_EFFECTIVE=('YES' if row['R_inlet']<=c['gates']['R_inlet'] else 'NO') if row.get('R_inlet') is not None else 'UNVERIFIED',
      STEP3C_AUTO_CHECK='PASS' if confirmed else 'FAIL',STEP3C_STATUS=summary.get('status',term['status']),FAILED_CONVERGENCE_GATES=summary.get('failed_convergence_gates',[str(error)] if error else ['NOT_EVALUATED']),
      PHYSICS_CHANGED='NO',NUMERICAL_CONTRACT_CHANGED='NO',GEOMETRY_CHANGED='NO',Q_TARGET_CHANGED='NO',OUTLET_PRESSURES_CHANGED='NO',VISCOSITY_CHANGED='NO',FURTHER_PURE_FLUID_GPU_OPTIMIZATION='STOP',
      RBC_INCLUDED='NO',IBM_INCLUDED='NO',PLASMA_VISCOSITY_CONTRACT='PENDING',READY_TO_ADD_RBC='NO',HUMAN_PARAVIEW_REVIEW='PENDING',PERFORMANCE_REGRESSION='YES' if reg else 'NO',
      SHA256_STATUS='PENDING_FINAL_SEAL',REPORT_DIR=str(R),STEP3C_FORMAL_GPU_REPORT=str(R/'STEP3C_FORMAL_GPU_REPORT.md'))
    write(R/'FINAL_TERMINAL_SUMMARY.json',vals)
    (R/'FINAL_TERMINAL_SUMMARY.txt').write_text('\n'.join(f'{k} = {json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v if v is not None else "NOT_REACHED"}' for k,v in vals.items())+'\n')
    if (R/'diagnostics/flow_history.csv').exists():
        data=history(R);names=['iteration','time_s','total_mass','control_volume_mass','relative_mass_drift','control_volume_mass_balance']+[f'mass_outward_g{g}' for g in [2,8,14,20]]
        with (R/'mass_history.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(names)
            for item in data:w.writerow([item[n] for n in names])
    lines=['# Step3C 正式 GPU 长时间验证','',f'正式 solver 终态：**{term["status"]}**。独立终审：**{"PASS" if result else "FAIL"}**。',
      '',f'从零执行 {actual} 步，物理时间 {actual*c["dt_s"]:.17g} s；进程端到端实际耗时 {term["end_to_end_seconds"]:.6f} s。完整迭代区间实际平均 {sps} steps/s（含输出、同步身份核查与评价）。',
      '', '**Qtarget 没有增加 11.97%。** u_command = 1.1197286861799598 × u_nominal 只补偿离散入口实现偏差；实测 Qin 继续与原始 2.7369132390905703e-15 m3/s 比较。',
      '', 'Step3B 未校准参照：Qin/Qtarget ≈ 0.89307，R_inlet ≈ 0.1069。',
      f'本次长期 Qin/Qtarget = {vals["FINAL_QIN_OVER_QTARGET"]}，冻结双窗口最大 R_inlet = {vals["FINAL_R_INLET"]}；入口校准有效：{vals["INLET_CALIBRATION_EFFECTIVE"]}。',
      '', '完整复用原物理窗口、240000 步首评、5000 步评价频率及至少 120000 步连续通过确认。任一中间失败重置候选。未收敛到 700000 步正常报告失败，不调参重跑。',
      '', 'Stage4 原二进制只支持短基准。正式二进制仅接入长程执行控制、原冻结评价器、原稀疏输出和同步 5000 步身份核查；GPU adapter、LBM batch、Guo、halo、设备监测计算、原生头文件保持逐字一致。',
      '', '每步保留原设备密度/速度/质量/24 组流量计算与安全门；标量文件每 100 步。保留全部评价及滞后全场供离线重算。独立终审复核所有评价与完整确认状态，并重建流量、质量和设备统计。',
      '', '2dx/4dx/6dx 平面完整保留，4dx 仍是 reference plane；所有出口同时通过短、长窗口平均正向门。启动期回流不单独触发停止。',
      '', '| 指标 | 实测 | 冻结上限 |','|---|---:|---:|']
    for key,gate in [('R_inlet','R_inlet'),('R_flow','R_flow'),('R_CV','R_CV'),('R_velocity','R_velocity'),('R_pressure','R_pressure'),('flow_fraction_drift','flow_fraction_drift'),('cross_plane_flux_spread','cross_plane_spread'),('quadrature_error','quadrature_dx_vs_dx2')]:lines.append(f'| {key} | {row.get(key)} | {c["gates"][gate]} |')
    lines+=['',f'首次候选：{state.get("first_converged_iteration")}；接受候选：{state.get("candidate_iteration")}；确认：{state.get("confirmation_iteration")}；停止：{actual}。',f'失败门：{vals["FAILED_CONVERGENCE_GATES"]}。',
      '',f'性能退化记录：{reg}；不参与物理停止或调参。',
      '', '来源、清理清单、删除前哈希和回收量见 provenance。仅清理 7 个未链接的 CUDA 数学静态库，原有科研源码、二进制及原始结果完整保留。最终可用 sha256sum -c SHA256SUMS 核查。',
      '', '纯流体 GPU 优化停止。RBC/IBM=NO；PLASMA_VISCOSITY_CONTRACT=PENDING；READY_TO_ADD_RBC=NO；HUMAN_PARAVIEW_REVIEW=PENDING。']
    if error:lines+=['',f'终审错误：{error}。']
    (R/'STEP3C_FORMAL_GPU_REPORT.md').write_text('\n'.join(lines)+'\n')
def main():
    try:
        result=audit()
        env=dict(os.environ,TMPDIR=str(R/'tmp'),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
        with (R/'logs/PARAVIEW_EXPORT.log').open('w') as log:
            ret=subprocess.run(['/usr/bin/python3','-B',str(R/'scripts/export_remote_fields.py'),str(R),str(R)],env=env,stdout=log,stderr=subprocess.STDOUT)
        assert ret.returncode==0,'ParaView export failure; numerical audit preserved'
        reports(result);print('FINAL_AUDIT PASS',flush=True)
    except Exception as e:
        write(R/'final_audit/FINALIZER_FAILURE.json',{'status':'FAIL','error':str(e)})
        reports(error=e);raise
if __name__=='__main__':main()
