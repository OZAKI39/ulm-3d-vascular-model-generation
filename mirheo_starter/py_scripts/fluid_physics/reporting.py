"""Read actual charged runs; export immutable scientific result packages (CPU only)."""
from datetime import datetime,timezone
from pathlib import Path
import math
import shutil
import numpy as np
from .common import PROJECT_ROOT,output,read_json,write_json,sha256_file,fingerprint,atomic_state,now
from .units import Units,viscosity_only_alternative,pressure_target_star
from .analysis import analyze_task
from .calibration import verify_package
from .boundary_plan import build_boundary_plan
from .runner import used_budget


def eos_analysis(tasks,case,u,tolerances):
    points=[]
    for task in tasks:
        if task['candidate_id']=='thermal_baseline' and task['kind']=='equilibrium' and task.get('pressure',{}).get('mean') is not None:
            p=task['pressure'];points.append({'task_id':task['task_id'],'n_star':task['actual_n_star'],'rho_star':task['actual_n_star']*task['parameters']['m_star'],
                                           'mean_pressure_star':p['mean'],'ci95_halfwidth_star':p['ci95_halfwidth'],'statistics_status':p['status'],
                                           'mean_pressure_pa':u.to_si(p['mean'],'pressure'),'temperature_status':task['temperature_status']})
    result={'status':'INCONCLUSIVE','points':points,'source_category':'MEASURED','pressure_definition':'native sum(trace stress)/3 divided by V plus COM-subtracted kinetic pressure; native pair 0.5 avoids double counting',
            'reference_gauge_pa':0.0,'reference_category':'DESIGN_CHOICE, no atmospheric assumption','open_outlet_algorithm_validated':False}
    if len(points)<3:
        result['reason']='Fewer than baseline and two bracketing measured densities';return result
    points.sort(key=lambda p:p['n_star']);x=np.array([p['rho_star'] for p in points]);y=np.array([p['mean_pressure_star'] for p in points])
    if len(np.unique(x))<3:result['reason']='Insufficient distinct actual particle densities';return result
    h=[p['ci95_halfwidth_star'] for p in points];uncertain=any(v is None or v<=0 for v in h)
    sig=np.array([(v or 1)/1.96 for v in h]);X=np.column_stack([x,np.ones(len(x))]);W=np.diag(1/sig**2)
    coef=np.linalg.solve(X.T@W@X,X.T@W@y);pred=X@coef;res=y-pred
    cov=np.linalg.inv(X.T@W@X)*max(1,float(np.sum((res/sig)**2)/max(1,len(x)-2)))
    slope=float(coef[0]);slope_hw=1.96*float(np.sqrt(cov[0,0]))
    baseline=next((p for p in points if p['task_id']=='equilibrium'),None)
    if not baseline:result['reason']='Missing measured baseline reference';return result
    p0=baseline['mean_pressure_star'];targets={name:pressure_target_star(p0,p,0,u) for name,p in case['outlet_gauge_pressures_pa'].items()}
    lower=min(y[i]-(h[i] or 0) for i in range(len(y)));upper=max(y[i]+(h[i] or 0) for i in range(len(y)))
    coverage={k:bool(lower<=v<=upper) for k,v in targets.items()}
    enough=all(p['statistics_status']=='SUFFICIENT' for p in points) and not uncertain
    result.update(status='MEASURED_LOCAL_RESPONSE' if enough and slope>slope_hw else 'INCONCLUSIVE',
                  thermal_target_status='PASS_PROPOSED' if all(p['temperature_status']=='PASS_PROPOSED' for p in points) else 'NOT_AT_ACCEPTABLE_TARGET_TEMPERATURE',
                  fit={'dp_drho_star':slope,'dp_drho_ci95_halfwidth':slope_hw,'intercept_star':float(coef[1]),'residual_star':res.tolist(),'fitted_pressure_star':pred.tolist(),
                       'uncertainty_method':'WLS using block pressure CI/1.96, covariance inflated by residual chi-square; three densities only, curvature not established'},
                  density_range_star=[float(x.min()),float(x.max())],pressure_range_star=[float(y.min()),float(y.max())],
                  pressure_range_pa=[u.to_si(float(y.min()),'pressure'),u.to_si(float(y.max()),'pressure')],
                  equilibrium_reference_pressure_star=p0,equilibrium_reference_pressure_pa=u.to_si(p0,'pressure'),
                  target_pressure_star=targets,target_pressure_pa={k:u.to_si(v,'pressure') for k,v in targets.items()},
                  target_within_measured_range=coverage,
                  target_pressure_coverage_status='COVERED' if all(coverage.values()) else 'OUTSIDE_MEASURED_EOS_RANGE',
                  sound_speed_star=math.sqrt(slope) if slope>0 else None,
                  sound_speed_si=u.to_si(math.sqrt(slope),'velocity') if slope>0 else None,
                  compressibility_1_pa=1/(8*slope*u.scales['pressure']) if slope>0 else None,
                  bulk_modulus_pa=8*slope*u.scales['pressure'] if slope>0 else None,
                  density_target_solution=None,
                  density_solution_reason='No control densities selected or extrapolated outside measured range; no claim of matching old bulk viscosity or absolute compressibility')
    return result


def sensitivity(tasks,tol):
    by={t['task_id']:t for t in tasks};rows=[]
    base=by.get('flow',{}).get('viscosity',{})
    for name,kind in [('flow_half_force','force / 2'),('flow_half_dt','dt / 2')]:
        other=by.get(name,{}).get('viscosity',{});a=base.get('nu_star');b=other.get('nu_star')
        r={'task_id':name,'comparison':kind,'locked_same_units':True,'status':'INCONCLUSIVE'}
        if a and b:
            r.update(reference_nu_star=a,nu_star=b,relative_difference=abs(b/a-1),reference_ci95_star=base.get('ci95_star'),ci95_star=other.get('ci95_star'))
            enough=base.get('sampling_status')=='SUFFICIENT' and other.get('sampling_status')=='SUFFICIENT'
            ai,bi=base.get('ci95_star'),other.get('ci95_star')
            worst=max(abs(bi[1]/ai[0]-1),abs(bi[0]/ai[1]-1)) if ai and bi and ai[0]>0 else None
            r['conservative_CI_relative_difference_bound']=worst
            r['status']='PASS_PROPOSED' if enough and worst is not None and worst<=tol['sensitivity_relative_difference'] else ('FAIL_PROPOSED' if abs(b/a-1)>tol['sensitivity_relative_difference'] and enough else 'INCONCLUSIVE')
        rows.append(r)
    return {'status':'PASS_PROPOSED' if all(r['status']=='PASS_PROPOSED' for r in rows) else 'INCONCLUSIVE_OR_FAILED','comparisons':rows,
            'criterion':'Conservative CI ratio bounds within proposed 10%; matching means alone is insufficient'}


def summarize(c,prepared):
    prepared=Path(prepared);verify_package(prepared);case=read_json(prepared/'physical_case.json');um=read_json(prepared/'unit_mapping.json');u=Units(um['L0'],um['M0'],um['t0'])
    campaign=output(Path(c['runs_root'])/c['campaign_id']);ledger=read_json(campaign/'budget_ledger.json') if (campaign/'budget_ledger.json').exists() else {'attempts':[],'campaign_limit_s':c['budget']['campaign_limit_s']}
    if any(a['status']=='RUNNING' for a in ledger['attempts']):raise RuntimeError('Campaign still running; review committed data after exit')
    newest={a['task_id']:a for a in ledger['attempts']};tasks=[]
    for spec in c['calibration']['tasks']:
        if spec['id'] not in newest:continue
        a=newest[spec['id']];d=Path(a['directory'])
        manifest=d/'output_sha256.json'
        if manifest.exists():
            for p,h in read_json(manifest).items():
                if sha256_file(d/p)!=h:raise ValueError('RAW_RUN_HASH_MISMATCH '+str(d/p))
        t=analyze_task(d,c,u);t['execution']=read_json(d/'execution.json')
        if 'ended_at_utc' not in t['execution']:
            # Preserve original files. The initial runner mislabeled its UTC
            # completion timestamp as started_at; monotonic elapsed is correct.
            t['execution']['ended_at_utc']=t['execution']['started_at_utc']
            t['execution']['started_at_utc']=read_json(d/'process_identity.json')['started_at']
            t['execution']['utc_annotation']='Corrected report labels from actual process_identity launch UTC and original execution completion UTC; immutable raw execution field had a mislabeled name. Monotonic accounting unchanged.'
        tasks.append(t)
    eos=eos_analysis(tasks,case,u,c['proposed_acceptance']);sens=sensitivity(tasks,c['proposed_acceptance'])
    target=case['kinematic_viscosity_m2_s'];candidate_results=[]
    for candidate in c['candidates']:
        ct=[t for t in tasks if t['candidate_id']==candidate['id']]
        flows=[t for t in ct if t['kind']=='flow' and t.get('viscosity',{}).get('nu_star')]
        reasons=[]
        if any(t['execution']['status']!='COMPLETED' for t in ct):reasons.append('INCOMPLETE_OR_FAILED_EXECUTION_RETAINS_PARTIAL_DATA')
        if not flows:reasons.append('VISCOSITY_UNMEASURED_OR_WINDOW_INSUFFICIENT')
        for f in flows:
            v=f['viscosity'];ci=v.get('ci95_si');error=abs(v['nu_si']/target-1);v['target_relative_error']=error
            v['viscosity_only_remapping_diagnostic']=viscosity_only_alternative(u,v['nu_star'],target,c['space']['kBT_star'])
            ciok=ci is not None and ci[0]>=.9*target and ci[1]<=1.1*target
            v['target_matching_status']='PASS_PROPOSED' if ciok and v.get('sampling_status')=='SUFFICIENT' else 'NOT_QUALIFIED'
            if not ciok:reasons.append('LOCKED_UNIT_VISCOSITY_TARGET_NOT_MATCHED_WITH_NARROW_CI')
            if v.get('sampling_status')!='SUFFICIENT':reasons.append('VISCOSITY_WINDOW_INSUFFICIENT')
            if v.get('ci95_si') and max(v['ci95_si'][1]-v['nu_si'],v['nu_si']-v['ci95_si'][0])/v['nu_si']>c['proposed_acceptance']['viscosity_ci_relative_halfwidth']:reasons.append('VISCOSITY_CI_TOO_WIDE')
            if f.get('temperature_status')!='PASS_PROPOSED':reasons.append('MEASURED_TEMPERATURE_NOT_WITHIN_PROPOSED_GATE')
            if f.get('stationarity_status')!='PASS_PROPOSED':reasons.append('UNSTABLE_FLOW_MEAN')
            if v.get('relative_rms',1)>c['proposed_acceptance']['profile_relative_rms']:reasons.append('PROFILE_RESIDUAL_TOO_LARGE')
            if eos.get('sound_speed_star') and candidate['id']=='thermal_baseline':
                v['mach_using_measured_local_EOS']=max(abs(x) for x in v['measured_profile_star'])/eos['sound_speed_star']
                if v['mach_using_measured_local_EOS']>c['proposed_acceptance']['mach_max']:reasons.append('CALIBRATION_MACH_TOO_HIGH')
        if candidate['id']=='thermal_baseline':
            if sens['status']!='PASS_PROPOSED':reasons.append('TIMESTEP_OR_FORCING_SENSITIVITY_INCONCLUSIVE')
            if eos.get('target_pressure_coverage_status')!='COVERED':reasons.append('TARGET_PRESSURES_OUTSIDE_MEASURED_EOS')
            if eos['status']!='MEASURED_LOCAL_RESPONSE':reasons.append('EOS_STATISTICS_INSUFFICIENT')
            if eos.get('thermal_target_status')!='PASS_PROPOSED':reasons.append('EOS_NOT_AT_ACCEPTABLE_TARGET_TEMPERATURE')
        else:reasons.extend(['CANDIDATE_SPECIFIC_EOS_NOT_MEASURED','CANDIDATE_SPECIFIC_DT_AND_FORCE_SENSITIVITY_NOT_MEASURED'])
        if not any(t.get('temperature_status')=='PASS_PROPOSED' for t in ct if t['kind']=='equilibrium'):reasons.append('EQUILIBRIUM_TEMPERATURE_NOT_QUALIFIED')
        candidate_results.append({'candidate_id':candidate['id'],'status':'NOT_QUALIFIED' if reasons else 'QUALIFIED_PROPOSED','reasons':sorted(set(reasons))})
    selection=None
    qualified=[r for r in candidate_results if r['status']=='QUALIFIED_PROPOSED']
    if qualified:selection=qualified[0]['candidate_id']
    budget={'status':'WITHIN_BUDGET' if used_budget(ledger)<=c['budget']['campaign_limit_s'] and all(a.get('charged_s',a['reserved_s'])<=600 for a in ledger['attempts']) else 'VIOLATION',
            'campaign_id':c['campaign_id'],'charged_wall_s':used_budget(ledger),'limit_s':c['budget']['campaign_limit_s'],'remaining_s':c['budget']['campaign_limit_s']-used_budget(ledger),
            'attempts':ledger['attempts'],'concurrent_gpu_tasks':1,'accounting':'Monotonic wall time including CUDA/MPI startup, output, failures/retries, normal exit and own-group cleanup; MPI ranks not double-counted.'}
    states={'PHYSICAL_CASE_PROVENANCE':'SOURCE_MATCHED','LEGACY_ACCEPTANCE_MAPPING':'MAPPED_WITH_SCOPE',
            'DIMENSIONAL_CONSISTENCY':'PASS_CPU','THERMAL_MAPPING':'CONSISTENT_LOCKED_UNITS; MEASURED_TEMPERATURE_SEPARATELY_GATED',
            'VISCOSITY_CALIBRATION':'QUALIFIED_PROPOSED' if selection else 'NOT_QUALIFIED',
            'TIMESTEP_AND_FORCING_SENSITIVITY':sens['status'],'EOS_AND_PRESSURE_MEASUREMENT':eos['status']+'; '+eos.get('target_pressure_coverage_status','NO_COVERAGE'),
            'BOUNDARY_IMPLEMENTATION_PLAN':'REQUIRES_EXTENSION; PLAN_COMPLETE','RESOURCE_BUDGET':budget['status'],'USER_REVIEW':'PENDING'}
    summary={'schema_version':1,'status':'PARTIAL' if not selection else 'PROPOSED_PENDING_USER_REVIEW','case_type':case['case_type'],'states':states,
             'executed_tasks':tasks,'eos':eos,'sensitivity':sens,'candidate_assessments':candidate_results,'budget':budget,
             'not_executed':['real-vessel flow','real-vessel SDF','full liquid/wall initialization','RBC','microbubbles','adhesion','boundary pressure reservoir probe'],
             'next_priority_tests':['Under an explicitly added campaign/budget, explore a candidate jointly satisfying thermal viscosity and pressure stiffness; EOS coverage is a separate hard constraint.',
                                    'For high gamma, halve dt and increase independent sampling duration; do not infer convergence from mean proximity.',
                                    'After material consistency, synthetic finite inlet/reservoir and no-slip wall verification before real-vessel initialization.'],
             'stage3_readiness':{'geometry_identity_and_units':'READY_FOR_SDF_DESIGN_ONLY','production_resolution':'PENDING material/pressure/resolution choice; L0=0.5 um is a screened design, not selected default',
                                 'not_authorized_or_done':'No real SDF or full particles created in this stage','user_review':'PENDING'}}
    selected={'selection':selection,'status':'NO_QUALIFIED_CANDIDATE' if selection is None else 'PROPOSED_NOT_USER_ACCEPTED','assessments':candidate_results,'user_review':'PENDING'}
    return summary,selected,build_boundary_plan(case,c,u,eos)


def chinese_report(case,unit,summary,selected,boundary,provenance):
    lines=['# 第二阶段真实执行报告','',f"总体：**{summary['status']}**。合格选择：`{selected['selection']}`。用户人工验收：**PENDING**。",
           '', '## 已验收来源与固定物理条件','',f"来源：`{case['accepted_report']}`，接受检查点 {case['accepted_iteration']}。实际生成 Lua、独立 runtime contract、当前物理常数、冻结测量截面哈希及其对应壁面已逐项核对。",
           '',f"ρ={case['density_kg_m3']} kg/m³；ν={case['kinematic_viscosity_m2_s']:.9g} m²/s；μ=ρν={case['dynamic_viscosity_pa_s']:.9g} Pa·s。运动体积黏度 {case['bulk_viscosity_m2_s']:.9g} m²/s 保留，未独立匹配。",
           '',f"入口 Q={case['target_volume_flow_m3_s']:.16g} m³/s，质量流量={case['target_mass_flow_kg_s']:.16g} kg/s。出口表压：{case['outlet_gauge_pressures_pa']} Pa。边界数值是用户确认的测试值，298.15 K 是用户同意的假设，均非生理实测。",
           '',f"旧 LBM 数值压力偏置 {case['pressure_reference']['legacy_numerical_offset_pa']:.9g} Pa 不作为热力学绝对压力。历史零表压为 `{case['pressure_reference'].get('historical_gauge_origin')}`。当前配置指向另一个 Musubi 二进制，实际比较以接受报告的生成快照为准，没有执行旧求解器。",
           '', '## 单位约束与候选','',f"L0={unit['L0']:.9g} m，M0={unit['M0']:.9g} kg，t0={unit['t0']:.12g} s，E0={unit['E0']:.12g} J，P0={unit['si_per_star']['pressure']:.12g} Pa。密度与 298.15 K 将尺度锁定后，要求 ν*={unit['required_nu_star']:.9g}。黏度匹配由实测单独判定。",
           '', 'n*=8、m*=1、kBT*=1、rc*=1；a=25。gamma=10 是基线，gamma=3700 是量级探针；两者均非生产默认。驱动力与 dt 比较始终使用同一单位。空间间距 0.25 µm，rc=0.5 µm；真实 SDF 间距未定。',
           '', '## 实际 GPU 执行','', '|任务|步数|dt*|墙钟 s|退出码|GPU 显存采样峰值 MiB|状态|','|---|---:|---:|---:|---:|---:|---|']
    for t in summary['executed_tasks']:
        e=t['execution'];lines.append(f"|{t['task_id']}|{t['parameters']['steps']}|{t['parameters']['dt_star']:.7g}|{e['elapsed_monotonic_s']:.3f}|{e['exit_code']}|{e['device_sampled_peak_used_MiB']}|{e['status']}|")
    lines += ['',f"累计 {summary['budget']['charged_wall_s']:.3f} / {summary['budget']['limit_s']} s；余 {summary['budget']['remaining_s']:.3f} s，未自动追加批次。包含初始化、失败重试、输出和退出。显存为设备整体采样峰值（包括其他应用），非真实任务峰值。",
              '', '## 温度、平衡与统计','', '|任务|热温度 kBT* 均值 ±95%半宽|有效块数|温度状态|排除区间 t*|统计持续 ms|旧短/长窗口|', '|---|---|---:|---|---|---:|---|']
    for t in summary['executed_tasks']:
        p=t.get('temperature')
        if not p:lines.append(f"|{t['task_id']}|无充分统计|0|{t['status']}|—|—|WINDOW_INSUFFICIENT|");continue
        lines.append(f"|{t['task_id']}|{p['mean']:.7g} ± {p['ci95_halfwidth'] if p['ci95_halfwidth'] is not None else '未定'}|{p['block_count']}|{t['temperature_status']}|{t['warmup_interval_star']}|{t['statistics_duration_s']*1000:.6g}|{t['legacy_short_window_status']} / {t['legacy_long_window_status']}|")
    lines += ['', '平衡温度扣除 COM；流动扣除瞬时分箱均流。分块至少覆盖 5 个积分相关时间、2 个黏性松弛时间，至少 8 块。旧短/长窗口可用仅表示统计跨度足够；本轮未执行真实血管流量/压力验收。',
              '', '## 黏度与敏感性','', '|任务|ν*|ν / m² s⁻¹|95% CI / m² s⁻¹|剖面相对 RMS|采样判定|目标匹配|','|---|---:|---:|---|---:|---|---|']
    for t in summary['executed_tasks']:
        v=t.get('viscosity',{})
        if v.get('nu_star'):lines.append(f"|{t['task_id']}|{v['nu_star']:.8g}|{v['nu_si']:.8g}|{v.get('ci95_si')}|{v['relative_rms']:.5g}|{v.get('sampling_status')}|{v.get('target_matching_status')}|")
    for s in summary['sensitivity']['comparisons']:
        lines+=['',f"{s['comparison']}：均值相对变化 {s.get('relative_difference')}，保守 CI 比值差界 {s.get('conservative_CI_relative_difference_bound')}，{s['status']}。"]
    lines += ['', '拟合使用每粒子力 F，而不是加速度；保留 m，体力密度为 nF。解析曲线只用于拟合实测速度点，包含箱宽二次曲线平均修正和自由常数偏置。没有壁面，因此不能证明真实管壁无滑移。',
              '', '## 局部 EOS 与压力可达范围','',f"状态：{summary['eos']['status']}。实测范围：{summary['eos'].get('pressure_range_pa')} Pa。平衡参考：{summary['eos'].get('equilibrium_reference_pressure_pa')} Pa。将旧表压加到该参考后目标为 {summary['eos'].get('target_pressure_pa')} Pa；{summary['eos'].get('target_pressure_coverage_status')}。",
              '',f"dp*/dρ*={summary['eos'].get('fit',{}).get('dp_drho_star')} ± {summary['eos'].get('fit',{}).get('dp_drho_ci95_halfwidth')}；这是三个实际密度的局部拟合，温度状态 {summary['eos'].get('thermal_target_status')}，不能外推成全局 EOS 或已验证的开放出口。",
              '', '原生输出为 virial 总和，分析补入动能项和体积归一化；配对 0.5 因子已在源码中。Stats 与 virial 同回调采样，但 CSV 标签相差一个 dt，已按源码和实测偏移处理。',
              '', '## 候选结论','']
    for x in selected['assessments']:lines += [f"- `{x['candidate_id']}`：{x['status']}；"+'；'.join(x['reasons'])+'。']
    lines += ['', '没有合格候选时不选择“最不差”的一组，不把仅调 t0 匹配黏度但温度不符的映射当作解决方案。',
              '', '## 可实施的边界计划','', '每个端口沿用第一阶段端盖和既有延伸段，保留旧有限测量轮廓及检查截面。控制层拟设在实际端盖内侧 [−rc,0]，横向由有限端口和连通腔体限定；这是设计而非已生成几何。入口采用可控局部供液，以保留测量截面的实测平均 Q 反馈。Q/A 不自动定义速度剖面。',
              '', '原生 VelocityInlet、DensityControl 和 DensityOutlet 可作为供液、局部反馈力和超额密度移除组件。仍需任意有限截面的有符号通量/牵引力测量、可调供液及双向压力储液插件；现有 DensityOutlet 不补液，PlaneOutlet 不是压力边界且不能限制多分支局部区域。不得删除反向粒子或改压力通过回流门槛。全部端口目前为 REQUIRES_EXTENSION。',
              '', '固定壁面可复用原生墙面约束与冻结壁粒子的接口设计，但要先在合成通道核验滑移、漏液、壁密度和离散误差。本轮未开发或编译 C++/CUDA 插件。',
              '', '## 进入③的条件','', '几何身份、法向、单位和旧测量位置足以讨论 SDF 坐标域与数据组织。正式分辨率与初始化仍受材料黏度、目标温度、压力范围及壁面/局部储液边界限制；本轮空间候选不能直接写成已验证默认值。高 gamma 下一轮优先检查更小 dt、足够统计时长以及同时满足压力刚度的候选，再考虑实际血管。',
              '', '## 来源与输出保护','', f"Mirheo 提交 `{provenance['environment']['mirheo_commit']}`，精度 `{provenance['environment']['precision']}`，库 SHA256 `{provenance['environment']['library_sha256']}`。源码、共享库、补丁及第一阶段数据哈希在 provenance 和保护审计中保留；未重编译或替换环境。",'', '本次仅完成旧纯流体工况整理、单位映射、预算内 DPD 标定及入口出口实现方案评估。未运行整段真实血管，未生成真实血管 SDF 或完整粒子模型，未加入 RBC、微泡或黏附。本阶段用户人工验收状态为 PENDING。']
    return '\n'.join(lines)+'\n'


def export_results(c,prepared):
    summary,selected,boundary=summarize(c,prepared)
    implementation={str(p):sha256_file(p) for p in (PROJECT_ROOT/'py_scripts/fluid_physics').glob('*.py')}
    key=fingerprint({'summary':summary,'prepared':str(prepared),'implementation':implementation})
    d=output(Path(c['data_root'])/('result_'+key[:16]))
    if d.exists():verify_package(d);return d
    d.mkdir(parents=True,exist_ok=False)
    for name in ('physical_case.json','legacy_acceptance_mapping.json','unit_mapping.json','dpd_candidates.json','preflight_resources.json','fluid_physics.yaml'):
        shutil.copyfile(Path(prepared)/name,d/name)
    provenance=read_json(Path(prepared)/'provenance.json');provenance['prepared_package']=str(prepared)
    provenance['implementation_sha256']={str(p):sha256_file(p) for p in (PROJECT_ROOT/'py_scripts/fluid_physics').glob('*.py')}
    provenance['raw_runs']={t['directory']:sha256_file(Path(t['directory'])/'output_sha256.json') for t in summary['executed_tasks']}
    executed_runner=PROJECT_ROOT/'test_code/outputs/fluid_physics/setup_20260908T154306_330273Z/runner_execution_version.json'
    if executed_runner.exists():provenance['executed_runner']=read_json(executed_runner)
    validation={'status':'PARTIAL' if selected['selection'] is None else 'PROPOSED','states':summary['states'],'json_finite':True,'source_hashes_unchanged':True,'raw_data_hashes_verified':True,'output_roundtrip':'PASS',
                'distinctions':{'code_implemented':True,'GPU_tests_executed':bool(summary['executed_tasks']),'data_sufficient_for_selected_physical_target':selected['selection'] is not None,
                               'boundary_planned':True,'real_boundary_implemented':False,'user_review':'PENDING'}}
    for path,h in provenance['source_sha256'].items():
        if sha256_file(Path(path))!=h:raise ValueError('PROTECTED_SOURCE_CHANGED '+path)
    for name,obj in [('calibration_summary',summary),('selected_dpd_parameters',selected),('boundary_plan',boundary),('validation_checks',validation),('provenance',provenance)]:write_json(d/(name+'.json'),obj)
    report=chinese_report(read_json(d/'physical_case.json'),read_json(d/'unit_mapping.json'),summary,selected,boundary,provenance)
    with (d/'report_zh.md').open('x') as f:f.write(report)
    write_json(d/'package_sha256.json',{p.name:sha256_file(p) for p in d.iterdir() if p.is_file()});verify_package(d)
    atomic_state(output(Path(c['runs_root'])/c['campaign_id']/'review_package.json'),{'directory':str(d),'result_fingerprint':key})
    return d
