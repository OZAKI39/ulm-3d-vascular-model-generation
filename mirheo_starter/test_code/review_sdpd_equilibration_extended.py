"""Offline extended view, sharing the existing Plotly/table/browser stack."""
import csv
import html
import json
import os
from pathlib import Path
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.offline import get_plotlyjs
from py_scripts.fluid_physics.common import read_json,write_json,sha256_file,fingerprint,output
from py_scripts.fluid_physics.calibration import verify_package
from test_code.review_sdpd_diagnostics import table


def write_page(c,package):
    package=Path(package);verify_package(package)
    identity={'manifest':sha256_file(package/'package_sha256.json'),'viewer':sha256_file(Path(__file__))}
    d=output(Path(c['review_root'])/package.name/fingerprint(identity)[:12]);page=d/'sdpd_equilibration_extended_review.html'
    if page.exists():
        if read_json(d/'html_record.json')['sha256']!=sha256_file(page):raise ValueError('HTML_CHANGED')
        return page
    d.mkdir(parents=True)
    plan=read_json(package/'run_plan.json');summary=read_json(package/'equilibration_summary.json')
    budget=read_json(package/'budget_request.json');supp=read_json(package/'historical_supplement.json')
    old=read_json(package/'historical_observable_status.json');new=read_json(package/'observable_status.json')
    execution=read_json(package/'restart_execution.json') if (package/'restart_execution.json').exists() else None
    checkpoint=read_json(package/'checkpoint_failure_audit.json') if execution else None
    with (package/'historical_raw_statistics.csv').open() as f:rows=list(csv.DictReader(f))
    data={k:np.array([float(r[k]) if r[k] else np.nan for r in rows]) for k in rows[0]}
    scale=plan['parameters']['locked_units']['t0']*1e6;t=data['time_star'];physical=t*scale;mask=(t>.3+1e-12)&(t<=.4+1e-12)
    figures={};full=go.Figure()
    full.add_trace(go.Scatter(x=physical.tolist(),y=data['kBT_COM_star'].tolist(),name='历史参考：旧 400000 步完整轨迹',mode='lines',line={'color':'#19758a'}))
    full.add_hrect(y0=.98,y1=1.02,fillcolor='#31a66a',opacity=.18,line_width=0)
    full.add_vrect(x0=.6*scale,x1=.8*scale,fillcolor='#467bd3',opacity=.15,line_width=0,
                  annotation_text='新正式窗口 (0.60, 0.80] · 尚未运行',annotation_position='top left')
    full.add_vrect(x0=.4*scale,x1=.6*scale,fillcolor='#b7bfc5',opacity=.12,line_width=0)
    full.update_layout(title='完整启动过程与冻结后续窗口（独立轨迹不拼接）',xaxis={'title':'绝对物理时间 / µs','range':[0,.8*scale]},yaxis_title='T / 298.15 K',showlegend=True)
    figures['temperature']=full
    late=go.Figure(go.Scatter(x=physical[mask].tolist(),y=data['kBT_COM_star'][mask].tolist(),name='历史 (0.30,0.40]：500 点',mode='lines',line={'color':'#19758a'}))
    late.add_hrect(y0=.98,y1=1.02,fillcolor='#31a66a',opacity=.18,line_width=0)
    late.update_layout(title='历史固定末段：温度仍有下降证据；未给稳态 CI',xaxis_title='旧实验绝对物理时间 / µs',yaxis_title='T / 298.15 K')
    figures['late']=late
    pressure=go.Figure(go.Scatter(x=physical[mask].tolist(),y=data['pressure_pa'][mask].tolist(),name='历史瞬时机械压力（每 200 步一次）',mode='lines',line={'color':'#875aba'}))
    mean=supp['details']['pressure_star']['fixed_formal_description']['mean']*plan['parameters']['locked_units']['si_per_star']['pressure']
    pressure.add_hline(y=mean,line_dash='dash',annotation_text='固定窗口描述均值；稳态 CI 未建立')
    pressure.update_layout(title='压力背景与强波动：不将随机瞬时值要求为常数',xaxis_title='旧实验绝对物理时间 / µs',yaxis_title='完整机械压力 / Pa')
    figures['blocks']=pressure
    evolution=make_subplots(rows=3,cols=1,vertical_spacing=.12,subplot_titles=['核数密度均值','核数密度方差','保存帧最近邻分位数（分布，非粒子静止）'])
    for index,key in enumerate(['kernel_number_mean','kernel_number_variance'],1):
        evolution.add_trace(go.Scatter(x=physical[1:].tolist(),y=data[key][1:].tolist(),name=key,mode='lines'),row=index,col=1)
    frames=old['STRUCTURE_STATIONARITY']['description']['snapshots']
    for j,q in enumerate(['q10','q50','q90'],1):
        evolution.add_trace(go.Scatter(x=[x['step']*plan['dt_star']*scale for x in frames],
            y=[x['nearest_distance_quantiles_star'][j] for x in frames],name='最近邻 '+q,mode='lines+markers'),row=3,col=1)
    evolution.update_layout(height=850,title='历史结构：持续变化与仅三张末段快照的限制')
    evolution.update_xaxes(title_text='旧实验绝对物理时间 / µs',row=3,col=1);figures['evolution']=evolution
    if execution:
        diagnostic=make_subplots(rows=4,cols=1,vertical_spacing=.09,
            subplot_titles=['冷启动温度（仅恢复诊断）','完整瞬时机械压力 / Pa','核数密度均值','COM 速度模长*'])
        for job,color in zip(execution['jobs'],['#126b8c','#bc5228','#8b6aae']):
            samples=job['samples']
            if not samples:continue
            x=[r['time_star']*scale for r in samples]
            values=[[r['kBT_COM_star'] for r in samples],
                    [r['pressure_pa'] for r in samples], [r['kernel_number_mean'] for r in samples],
                    [float(np.linalg.norm([r['mean_v'+axis] for axis in 'xyz'])) for r in samples]]
            for i,y in enumerate(values,1):
                diagnostic.add_trace(go.Scatter(x=x,y=y,name=job['task_id']+' · 独立冷启动',
                    legendgroup=job['task_id'],showlegend=i==1,mode='lines+markers',line={'color':color},
                    marker={'size':4}),row=i,col=1)
        diagnostic.add_vline(x=.002*scale,line_dash='dot',line_color='#875aba',row='all',col=1)
        diagnostic.update_layout(height=1080,title='真实 A/B 短诊断：各从 t*=0 开始；竖线为 B 保存边界，恢复后无数据')
        diagnostic.update_xaxes(title_text='各独立实验的绝对物理时间 / µs',row=4,col=1)
        figures['restart-diagnostic']=diagnostic
    for f in figures.values():f.update_layout(template='plotly_white',margin={'l':75,'r':30,'t':80,'b':65},font={'size':13},legend={'orientation':'h','y':-0.22})
    labels={'TEMPERATURE_STATIONARITY':'温度稳定性','TEMPERATURE_TARGET_MATCH':'温度与目标匹配',
            'PRESSURE_STATIONARITY':'压力稳定性','STRUCTURE_STATIONARITY':'结构稳定性',
            'SAMPLING_SUFFICIENCY':'独立采样充分性','RESTART_VALIDITY':'保存／恢复有效性','OVERALL_EQUILIBRIUM':'整体液体平衡'}
    reasons={'TEMPERATURE_STATIONARITY':'两半变化约 0.790%，超过 0.5% 筛选；补充趋势分析支持末段继续冷却。',
             'TEMPERATURE_TARGET_MATCH':'尚未建立可信的稳定温度及其置信区间，不能认定稳定偏高。',
             'PRESSURE_STATIONARITY':'瞬时噪声强、采样精度不足；趋势区间跨零，未证明总压力持续漂移。',
             'STRUCTURE_STATIONARITY':'核密度及最近邻分布仍变化；仅三张末段快照不足以证明分布稳定。',
             'SAMPLING_SUFFICIENCY':'500 点最多容纳 12 个最短块；实际独立信息未验证。',
             'RESTART_VALIDITY':'真实恢复对照未运行；当前原生 SDPD 缺少相互作用 RNG 持久化。',
             'OVERALL_EQUILIBRIUM':'各项证据分别保留，尚不能宣布液体整体平衡或其他物性通过。'}
    if execution:reasons['RESTART_VALIDITY']='真实恢复进程在读取 saved_forces / Other 时失败；恢复前、首步及短段比较未测到。SDPD RNG 持久化仍缺失。'
    status_labels={'PASS':'通过筛选','FAIL':'明确未通过','INCONCLUSIVE':'证据不足','NOT_TESTED':'未测试'}
    statuses=table(['指标','历史固定窗口','本轮新测量','依据'],[[labels[k],status_labels.get(old[k]['status'],old[k]['status']),
        status_labels.get(new[k]['status'],new[k]['status']),reasons[k]] for k in old],id='outcomes')
    budget_table=table(['预算项','秒'],[
        ['全局当前总授权',budget['budget_snapshot']['limit_s']],['全局已计费',budget['budget_snapshot']['total_charged_or_reserved_s']],
        ['本范围当前可用余额',budget['current_scope_remaining_s']],['本范围已登记追加',budget['extra_authorized_gpu_seconds']],
        ['本轮真实计费',summary.get('new_GPU_elapsed_s',0)],['真实恢复对照上限（3×30）',90],
        ['仅恢复通过后可用的冷启动上限（4×340）',1360],['一次请求总上限',1450],
        ['需追加的整数额度',budget['additional_request_whole_seconds']]],id='budget-table')
    def link(name,label):return '<a href="'+html.escape(os.path.relpath(package/name,d))+'">'+html.escape(label)+'</a>'
    diagnostic_section=''
    notice='本轮 GPU：NOT_RUN。新绝对时间：未测量。最近实际进度来自旧实验：t*=0.40，12.662309354 µs。'
    measurement_note='新测量全部为空；CPU 合成测试不作为实验结果。'
    event_note='当前实际保存/恢复事件：无。'
    approval_note='一次申请追加 '+str(budget['additional_request_whole_seconds'])+' 秒；当前本范围已登记追加 '+str(budget['extra_authorized_gpu_seconds'])+' 秒。'
    if execution:
        notice=(f"本轮已执行 {execution['GPU_attempt_count']} 个 GPU 进程，计费 {execution['charged_s']:.6f} 秒。恢复读取失败；正式长实验 NOT_RUN。"
                '最近液体长观察仍为旧 t*=0.40（12.662309354 µs），不能把短恢复对照算成后段平衡证据。')
        measurement_note='下方 A/B 短诊断为本轮真实数据；正式窗口没有新样本。两条独立轨迹分别展示，不拼接旧实验，不用于物性验收。'
        event_note='B 在第 2000 步、t*=0.002 保存，额外第 2001 步只触发写盘；新进程在读取时失败，尚未取得恢复后状态。'
        approval_note=(f"已登记追加 {budget['extra_authorized_gpu_seconds']:g} 秒；当前范围余额 {budget['current_scope_remaining_s']:.6f} 秒。"
                       '本次不再申请额度；停止后不自动重试，余额仍受原用途和恢复条件限制。')
        job_rows=[]
        for j in execution['jobs']:
            time_us='未测到' if j['actual_physical_time_s'] is None else f"{j['actual_physical_time_s']*1e6:.12f}"
            job_rows.append([j['task_id'],j['status'],j['actual_steps'] if j['actual_steps'] is not None else '未到达',time_us,f"{j['charged_s']:.6f}",j['raw_rows']])
        ab=checkpoint['comparisons']['independent_A_vs_B_before_any_restore']
        diagnostic_section=('<h2>本轮真实恢复测试</h2>'+table(['进程','结果','实际步数','实际时间 / µs','计费 / 秒','采样行数'],job_rows,id='restart-jobs')+
            '<p>直接错误：<code>saved_forces / Other / Force</code>；HDF 只有 [4096,1]，原生读取器报 <code>Unrecognised form Other</code>。'
            '<code>saved_stresses</code> 为受支持的 Tensor6，[4096,6]。CPU 已加入导入 CUDA 前的格式拒绝检查；未修复原生随机状态，也未重试 GPU。</p>'+
            '<p>归档位置/速度与 B 自身保存边界按粒子 ID 比较，最大差均为 0。独立 A/B 保存前最大速度差为 '+f"{ab['max_velocity_error_star']:.9g}"+
            '*，已超过原 2e-6* 比较门槛；它不是恢复造成的跳变，原门槛未放宽。恢复后的温度、压力、密度和动量跳变仍未测量。</p>'+
            '<p>'+link('restart_execution.json','三进程原始来源与计费')+' · '+link('checkpoint_failure_audit.json','checkpoint 格式、保存状态与失败证据')+'</p>'+
            '<div class="plot" id="restart-diagnostic" style="min-height:1080px"></div>')
    text=('<h1>固定参数 SDPD：恢复准备与后段观察</h1>'
        '<div class="notice"><strong>'+html.escape(summary['status'])+'</strong><p>'+notice+'</p>'
        '<p>计划为独立冷启动；旧运行无完整 checkpoint。原生 SDPD 未保存交互 RNG，状态 RESTART_NOT_VALIDATED；当前正式恢复链被阻止，预算不能单独解除这一限制。</p></div>'
        '<p>'+link('report_zh.md','中文报告')+' · '+link('run_plan.json','冻结计划')+' · '+link('restart_contract.json','原生恢复证据')+' · '+link('budget_request.json','完整预算请求')+' · <button id="reset-plots">恢复图表范围</button></p>'
        '<h2>各项结论与证据边界</h2>'+statuses+
        '<p>温度：相关噪声修正的斜率仍为负，支持末段冷却。压力：斜率区间跨零，尚未证明持续漂移。结构：最近邻分布和核密度仍变化；三张末段快照不足以给出分布稳态置信度。</p>'
        '<p>旧实验结论与门槛保留。'+measurement_note+'生产 selection=null；人工验收 PENDING。</p>'+diagnostic_section+
        '<h2>完整温度与新旧实验边界</h2><p>蓝色带是尚未开始的新正式窗口；绿色带为目标 ±2%。没有把旧末态快照改成连续恢复。</p>'
        '<div class="plot" id="temperature"></div><div class="plot" id="late"></div>'
        '<h2>压力与采样成本</h2><p>每 200 步一个瞬时值，包含保守、耗散及随机应力，并非 200 步平均。正式窗口标准差约 '+f"{supp['pressure_sampling']['formal_sigma_pa']:.3f}"+' Pa；稳定性门槛约 0.26345 Pa。压力时间平均尚无合格 CI，点间波动不能当作平均值误差。</p>'
        '<div class="plot" id="blocks"></div><p>PROPOSED（未采用）：评估逐步累计完整机械压力的成本，保留随机/耗散贡献。不能靠删随机应力、用输入 EOS 替代实测或放宽门槛得到通过。</p>'
        '<h2>结构与独立信息</h2><div class="plot tall" id="evolution"></div>'
        '<p>500 点按物理最小块长最多容纳 12 块；独立性未验证，不声称已有 12 个独立块。1000 点后续窗口也不保证稳定。结构分布稳定不要求粒子位置和邻居停止变化。</p>'
        '<h2>保存位置、恢复状态与预算</h2><p>'+event_note+'计划的正式保存点为 0.20、0.40、0.60、0.80*，均尚未执行；以 checkpoint 中的全局步数/时间为准。缺 RNG、文件、哈希、ID、时间一致性或真实验证，均停止正式链。</p>'+budget_table+
        '<p>并发 1，单任务 ≤600 秒；恢复对照每进程 30 秒，条件长段各 340 秒。'+approval_note+'恢复未通过，1360 秒长段预算不支出，不自动从零重跑或突破 600 秒改成长进程。</p>'
        '<p>最长 t*=0.80，正式窗口 (0.60,0.80]，只观察窗口 (0.40,0.60]。到上限即停止，不自动移动窗口或延长。</p>')
    payload={'kind':'extended','planned_steps':800000,'formal_window_star':[.6,.8],
             'actual_steps':None,'new_temperature_mean':None,'stationary':None,'temperature_match':None,
             'raw_csv_rows':0,'historical_rows':len(rows),'historical_last_time_star':float(t[-1]),
             'restart_validity':'RESTART_NOT_VALIDATED','selection':None,'human_review':'PENDING',
             'extra_authorized_gpu_seconds':budget['extra_authorized_gpu_seconds'],'requested_additional_s':budget['additional_request_whole_seconds'],
             'execution_status':summary['execution_status'],'GPU_attempt_count':execution['GPU_attempt_count'] if execution else 0,
             'new_GPU_elapsed_s':execution['charged_s'] if execution else 0,
             'diagnostic_raw_rows':[j['raw_rows'] for j in execution['jobs']] if execution else [],
             'checkpoint_step':checkpoint['checkpoint_absolute_step'] if checkpoint else None,
             'restore_preadvance':checkpoint['restored_preadvance_state'] if checkpoint else 'NOT_RUN',
             'automatic_retry':False,'current_scope_remaining_s':budget['current_scope_remaining_s']}
    plots={k:json.loads(f.to_json()) for k,f in figures.items()}
    javascript='const figures='+json.dumps(plots,ensure_ascii=False).replace('</',r'<\/')+';window.equilibrationReady=false;Promise.all(Object.entries(figures).map(([id,f])=>Plotly.newPlot(id,f.data,f.layout,{responsive:true,scrollZoom:true,displaylogo:false}))).then(()=>window.equilibrationReady=true);document.getElementById("reset-plots").onclick=()=>Promise.all(Object.entries(figures).map(([id,f])=>Plotly.relayout(id,{"xaxis.autorange":true,"yaxis.autorange":true})));'
    page.write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>SDPD 恢复与后段观察准备</title>'
        '<style>body{font-family:Arial,"Microsoft YaHei",sans-serif;background:#f3f6f7;color:#18323d;margin:0}main{max-width:1260px;margin:auto;padding:30px}p{line-height:1.8}h2{margin-top:32px}.notice{background:#fff1d9;border-left:5px solid #b77b2c;padding:18px}.plot{background:white;min-height:480px;margin:20px 0}.tall{min-height:850px}.table-wrap{overflow:auto}table{border-collapse:collapse;width:100%;background:white}th,td{padding:10px;border:1px solid #d6e1e5;text-align:left;overflow-wrap:anywhere}th{background:#e5edf0}button{padding:8px;cursor:pointer}</style>'
        '<main>'+text+'</main><script id="equilibration-audit-data" type="application/json">'+json.dumps(payload,ensure_ascii=False)+'</script><script>'+get_plotlyjs()+'</script><script>'+javascript+'</script></html>',encoding='utf-8')
    write_json(d/'html_record.json',{'sha256':sha256_file(page),'package':str(package),'identity':identity,'human_review':'PENDING'})
    return page
