"""Offline recorded-evidence review. Never launches a Mirheo/GPU experiment."""
import argparse
import csv
import html
import json
import os
from pathlib import Path
import subprocess
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.offline import get_plotlyjs
from py_scripts.fluid_physics.common import PROJECT_ROOT,read_json,write_json,sha256_file,fingerprint,output
from py_scripts.fluid_physics.calibration import verify_package
from py_scripts.sdpd_diagnostics.workflow import load_config,export
from test_code.review_vessel_geometry import open_windows_file


def fmt(v):
    if v is None:return '未测/无有效 CI'
    if isinstance(v,(int,float)):return f'{v:.6g}'
    return html.escape(str(v))


def table(headers,rows,id=None):
    return '<div class="table-wrap"><table'+(' id="'+id+'"' if id else '')+'><thead><tr>'+''.join('<th>'+html.escape(x)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+fmt(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table></div>'


def make_figures(s,a):
    figs={};scale=s['units']['t0']*1e6;temps=a['temperature_audit'];by={t['task_id']:t for t in temps}
    temp=go.Figure()
    for r in temps:
        q=r['series'];temp.add_trace(go.Scatter(x=np.array(q['time_star'])*scale,y=q['thermal'],mode='lines',name=r['task_id'],visible=True if r['task_id'].startswith('sdpd_') else 'legendonly'))
    temp.add_hrect(y0=.98,y1=1.02,opacity=.15,line_width=0,fillcolor='green')
    temp.update_layout(title='完整温度轨迹：初始化峰值与末端均保留',xaxis_title='物理时间 / µs',yaxis_title='T / 298.15 K');figs['temperature']=temp
    windows=make_subplots(rows=1,cols=2,subplot_titles=['所有固定截段的均温','窗口两半的温度漂移'])
    for name in ['sdpd_equilibrium','sdpd_flow','sdpd_half_dt','sdpd_half_force']:
        rows=by[name]['windows'];x=[r['cut_fraction'] for r in rows]
        for col,key in [(1,'temperature_mean_all'),(2,'thermal_drift')]:
            windows.add_trace(go.Scatter(x=x,y=[r[key] for r in rows],mode='lines+markers',name=name,legendgroup=name,showlegend=col==1),row=1,col=col)
    windows.add_hline(y=.02,line_dash='dash',row=1,col=2)
    windows.update_xaxes(title_text='截去前段占总时长的比例');windows.update_layout(title='所有固定窗口');figs['windows']=windows
    diff=go.Figure()
    for r in temps:
        if not r['task_id'].startswith('sdpd_'):continue
        q=r['series'];diff.add_trace(go.Scatter(x=np.array(q['paired_time_star'])*scale,y=q['adjacent_step_native_minus_worker'],mode='lines',name=r['task_id']))
    diff.update_layout(title='原生 Stats − worker：源码确认相差一个积分步',xaxis_title='worker 物理时间 / µs',yaxis_title='温度比之差（并非同时刻误差）');figs['difference']=diff
    evo=make_subplots(rows=3,cols=1,shared_xaxes=True,vertical_spacing=.09,subplot_titles=['无驱动温度','机械压力：积分前应力 + 积分后动量','积分前局部核密度'])
    q=by['sdpd_equilibrium']['series']
    for row,x,y,label in [(1,'time_star','thermal','T/T目标'),(2,'pressure_time_star','mechanical_pressure_pa','Pa'),(3,'density_time_star','kernel_mean','核数密度*')]:
        evo.add_trace(go.Scatter(x=np.array(q[x])*scale,y=q[y],mode='lines',name=label),row=row,col=1)
        evo.update_yaxes(title_text=label,row=row,col=1)
    evo.update_xaxes(title_text='物理时间 / µs',row=3,col=1);evo.update_layout(title='sdpd_equilibrium：温度、压力和结构的对应关系',height=680,showlegend=False);figs['evolution']=evo
    profile=make_subplots(rows=2,cols=1,shared_xaxes=True,vertical_spacing=.17,subplot_titles=['完整块平均剖面与目标、拟合','相同完整块数据的残差'])
    options=[]
    for r in a['profile_audit']:
        indices=[];default=r['task_id']=='sdpd_flow'
        for y,label,mode in [(r['corrected_measured_star'],'实测点','markers'),(r['target_star'],'目标黏度曲线','lines'),(r['corrected_fit_star'],'拟合曲线','lines')]:
            indices.append(len(profile.data));profile.add_trace(go.Scatter(x=r['y_star'],y=y,name=label,mode=mode,visible=default),row=1,col=1)
        indices.append(len(profile.data));profile.add_trace(go.Scatter(x=r['y_star'],y=(np.array(r['corrected_measured_star'])-r['corrected_fit_star']),name='实测 − 拟合',mode='lines+markers',visible=default),row=2,col=1)
        options.append({'task_id':r['task_id'],'label':r['task_id']+' / '+str(r['valid_temporal_blocks'])+' 块，无合格 CI','indices':indices})
    profile.update_yaxes(title_text='uₓ*',row=1,col=1);profile.update_yaxes(title_text='残差*',row=2,col=1)
    profile.update_xaxes(title_text='y*（一个周期盒 0–8；变号面 0、4、8）',row=2,col=1)
    profile.update_layout(title='黏度为拟合估计：保留所有空间位置',height=620);figs['profile']=profile
    eos=go.Figure();margins=a['pressure_audit']['endpoint_margins'];x=[r['endpoint'] for r in margins]
    eos.add_trace(go.Scatter(x=x,y=[r['pressure_difference_pa'] for r in margins],mode='markers',name='机械压差与原保守半宽',error_y={'type':'data','array':[r['conservative_endpoint_plus_reference_halfwidth_pa'] for r in margins]},marker={'size':12}))
    eos.add_trace(go.Scatter(x=x,y=[r['required_gauge_pa'] for r in margins],mode='markers',name='旧目标表压差',marker={'size':13,'symbol':'x'}))
    eos.update_layout(title='均值名义覆盖，但误差条内缘不能覆盖目标',yaxis_title='相对同一参考态的压力 / Pa');figs['eos']=eos
    probe=go.Figure();p=s['new_probe']
    if p and p.get('series'):
        q=p['series']
        for prefix,label in [('baseline','旧无驱动 dt'),('half_dt','新无驱动 dt/2')]:
            probe.add_trace(go.Scatter(x=np.array(q[prefix+'_time_star'])*scale,y=q[prefix+'_temperature'],name=label,mode='lines'))
        for lo,hi in a['proposed_experiments']['selected_probe']['comparison_windows_star']:probe.add_vline(x=hi*scale,line_dash='dot',opacity=.4)
    else:probe.add_annotation(text='未执行/无完整可读取数据',showarrow=False)
    probe.update_layout(title='唯一新增探针：同物理时段的启动趋势（不作稳态验收）',xaxis_title='物理时间 / µs',yaxis_title='COM 温度比');figs['probe']=probe
    budget=go.Figure();b=s['shared_budget']
    for r in b['members']:budget.add_trace(go.Bar(y=['原授权 3600 秒'],x=[r['charged_or_reserved_s']],orientation='h',name=r['campaign_id'],text=[f"{r['charged_or_reserved_s']:.2f} s"]))
    budget.add_trace(go.Bar(y=['原授权 3600 秒'],x=[b['remaining_s']],orientation='h',name='剩余',text=[f"{b['remaining_s']:.2f} s"]))
    budget.update_layout(title='跨所有 campaign 累计的同一预算',barmode='stack',height=290,xaxis_title='含初始化、采样、失败和退出的墙钟 / s');figs['budget']=budget
    for f in figs.values():f.update_layout(template='plotly_white',font={'family':'Arial, Microsoft YaHei, sans-serif'},margin={'l':75,'r':25,'t':80,'b':70},legend={'orientation':'h','y':-.22})
    budget.update_layout(legend={'orientation':'h','y':1.08,'yanchor':'bottom','x':0})
    return figs,options


def write_page(c,package):
    verify_package(package);s=read_json(package/'diagnosis_summary.json')
    names=['evidence_matrix','temperature_audit','profile_audit','pressure_audit','proposed_experiments','repair_differences','matched_physical_time_audit']
    a={n:read_json(package/(n+'.json')) for n in names}
    with (package/'sampling_audit.csv').open() as f:sampling=list(csv.DictReader(f))
    checker=Path(__file__).with_name('check_sdpd_diagnostics_browser.cjs')
    key=fingerprint({'package':sha256_file(package/'package_sha256.json'),'viewer':sha256_file(Path(__file__)),'checker':sha256_file(checker)})[:12]
    d=output(Path(c['review_root'])/package.name/key);page=d/'sdpd_diagnostics_review.html'
    if page.exists():
        if read_json(d/'html_record.json')['sha256']!=sha256_file(page):raise ValueError('HTML_CHANGED')
        return page
    figures,options=make_figures(s,a)
    evidence=table(['异常','证据状态','一句话解释'],[[r['phenomenon'],r['status'],r['one_sentence']] for r in a['evidence_matrix']],id='evidence-table')
    block=table(['任务/指标','计划/分配/实际步数','截段开始*','剩余样本','间隔*','ACF 样本','物理下限','5τ 取整','块长','尾部','块数'],[[r['task_id']+'/'+r['observable'],'/'.join(r[k] for k in ['desired_steps','allocated_steps','actual_steps']),*[r[k] for k in ['selected_start_star','selected_samples','sample_interval_star','tau_samples','minimum_block_samples','ACF_block_samples','final_block_samples','discarded_tail_samples','block_count']]] for r in sampling if r['task_id'].startswith('sdpd_')],id='sampling-table')
    repair=table(['任务','温度：原 → 修复','ν SI：原 → 修复','剖面 RMS：原 → 修复'],[[r['task_id'],*[fmt(r['before'][k])+' → '+fmt(r['after'][k]) for k in ['temperature_mean','nu_si','profile_relative_RMS']]] for r in a['repair_differences']],id='repair-table')
    windows=table(['任务','截段比例','均温','温度漂移','信号漂移','原规则选用'],[[r['task_id'],w['cut_fraction'],w['temperature_mean_all'],w['thermal_drift'],w['signal_drift'],w['selected_by_original_rule']] for r in a['temperature_audit'] if r['task_id'].startswith('sdpd_') for w in r['windows']])
    plan=a['proposed_experiments'];probe=s['new_probe']
    experiments=table(['实验','状态','判断范围/未执行理由'],[[plan['selected_probe']['id'],plan['selected_probe']['status'],plan['selected_probe']['interpretation_rule']]]+[[r['id'],r['status'],r['reason']] for r in plan['other_probes']])
    costs=table(['情景（均需追加授权）','步数','物理时长 µs','预计墙钟 s','所缺预算 s','假设'],[[r['scenario'],r['steps'],r['physical_duration_s']*1e6,r['predicted_wall_s'],r['additional_budget_needed_s'],r['assumptions']] for r in plan['eight_block_cost_scenarios']])
    matched=table(['对照','共同物理窗口*','任务','均温','ν 拟合 SI','相对 RMS'],[[r['task_id'],w['interval_star'],t['task_id'],t['temperature_mean'],t['nu_fit_si'],t['relative_RMS']] for r in a['matched_physical_time_audit'] for w in r['windows'] for t in w['records']])
    intro='<h1>SDPD 原因排查</h1><p><strong>已证实：启动过渡进入统计段，统计均值与 CI 的样本不一致。</strong></p><p>温升更支持未稳定；原生积分的稳态温偏与有限尺度剖面效应仍未定论。selection=null · 人工验收 PENDING。</p>'
    text=intro+'<p><a href="'+html.escape(os.path.relpath(package/'report_zh.md',d))+'">完整中文报告</a> · <button id="reset-plots">恢复图表范围</button></p>'+evidence
    plot=lambda name:'<section><div class="plot" id="'+name+'"></div></section>'
    text+='<h2>温度与所有诊断窗口</h2>'+plot('temperature')+plot('windows')+'<details><summary>展开全部窗口数值与原选择规则</summary>'+windows+'</details>'
    text+='<p>温度目标为 1±2%。最后 10% 窗口只用于诊断，未替换正式验收窗口。旧逐粒子速度分量与异常粒子能量分布：未测，无法从现有总矩恢复。</p>'+plot('difference')+plot('evolution')
    text+='<p>Stats 与 worker 读取同一底层速度，CPU 分箱矩复核也非独立仪器；相邻步差异不能当成同时刻误差。密度与配对应力是积分前状态，动量为积分后状态。</p><h2>一块统计的计算账单</h2>'+block
    text+='<h2>速度剖面与物理时长匹配</h2><select id="profile-select"></select><p id="profile-caption"></p>'+plot('profile')+matched+'<p>三条流动轨迹仍在过渡中；这些相同物理窗口的拟合不构成稳态敏感性验证。没有删除变号平面附近的任何数据。</p>'
    text+='<h2>压力覆盖与分析修复</h2>'+plot('eos')+'<p>两端压差共享基准，误差相关。图示原规则的端点与基准半宽之和；独立运行 Welch 只是条件分析，不能视为已测运行间协方差。</p>'+repair
    text+='<h2>最少对照及预算</h2>'+experiments+plot('probe')
    if probe:
        explanation='两个相同物理时段的无驱动任务都保留强启动温升与衰减，减半时间步未消除这一过程，更支持初始化与结构松弛。短探针尚未检验稳定后的温度偏差。' if probe['summary']['status']=='COMPLETED_FIXED_PHYSICAL_WINDOW' else probe['summary']['interpretation']
        text+='<p>'+html.escape(explanation)+'</p>'
    for failure in s.get('preserved_failed_diagnostic_attempts',[]):
        text+='<p>保留首次观测代码失败：请求 forces，实际原生通道为 __forces；初始化退出计费 '+fmt(failure['execution']['elapsed_monotonic_s'])+' 秒。回归测试复现并修正后执行上面的独立探针。此失败不构成物理证据。</p>'
    priority=plan['next_priority_experiment']
    text+='<p><strong>唯一优先补测：</strong>原 dt、无驱动，固定 '+str(priority['steps'])+' 步 / '+fmt(priority['physical_duration_s']*1e6)+' µs；预计 '+fmt(priority['predicted_total_wall_s'])+' 秒，建议预约 600 秒，需追加授权 '+fmt(priority['additional_authorization_needed_for_allocation_s'])+' 秒。预先固定最后 t*=0.30–0.40 观测窗；是否能得到 8 个稳定块尚未保证。REQUIRES_ADDITIONAL_AUTHORIZATION。</p>'
    text+=plot('budget')+'<p>下一步只优先建立无驱动目标液体的稳定热平台和独立统计块。下表是条件预算，尚未测得可信稳态相关时间，不自动执行。</p>'+costs
    audit={'package':str(package),'manifest_sha256':sha256_file(package/'package_sha256.json'),'figure_ids':list(figures),
        'profile_options':options,'evidence_count':len(a['evidence_matrix']),'selection':None,'human_review':'PENDING',
        'probe_status':probe['summary']['status'] if probe else 'NOT_TESTED','shared_budget':s['shared_budget']}
    payload=json.dumps({name:json.loads(f.to_json()) for name,f in figures.items()},ensure_ascii=False).replace('</','<\\/')
    js="""const figures=FIGURES,options=OPTIONS;window.diagnosticsReady=false;
Promise.all(Object.entries(figures).map(([id,f])=>Plotly.newPlot(id,f.data,f.layout,{responsive:true,scrollZoom:true,displaylogo:false}))).then(()=>window.diagnosticsReady=true);
const select=document.getElementById('profile-select');options.forEach((o,i)=>{let e=document.createElement('option');e.value=i;e.textContent=o.label;select.appendChild(e)});
select.value=String(options.findIndex(o=>o.task_id==='sdpd_flow'));document.getElementById('profile-caption').textContent=options[+select.value].label;
select.addEventListener('change',()=>{Plotly.restyle('profile',{visible:figures.profile.data.map((_,i)=>options[+select.value].indices.includes(i))});document.getElementById('profile-caption').textContent=options[+select.value].label});
document.getElementById('reset-plots').addEventListener('click',()=>Object.entries(figures).forEach(([id,f])=>{let ranges={};Object.keys(f.layout).filter(k=>/^([xy]axis)(\\d*)$/.test(k)).forEach(k=>ranges[k+'.autorange']=true);Plotly.relayout(id,ranges)}));
""".replace('FIGURES',payload).replace('OPTIONS',json.dumps(options))
    page_text='<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>SDPD 原因排查</title><style>body{font:16px/1.6 Arial,"Microsoft YaHei",sans-serif;color:#17263c;background:#f3f5f8;margin:0}main{max-width:1370px;margin:auto;padding:28px}section,.table-wrap{background:white;border-radius:7px;margin:20px 0;padding:10px;overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:9px;border:1px solid #dbe2eb;vertical-align:top}th{background:#eef3fa}.plot{height:520px}#evolution{height:760px}#profile{height:700px}#budget{height:330px}select,button{padding:8px;font-size:15px}code{overflow-wrap:anywhere}summary{cursor:pointer}a{color:#1d4ed8}</style></head><body><main>'+text+'<p>完整来源、原始数据哈希及算法版本：<code>'+html.escape(str(package))+'</code>。此页面离线，不启动 GPU 实验。浏览器实测记录另存 browser_final/browser_checks.json；人工验收始终 PENDING。</p></main><script id="diagnostics-audit-data" type="application/json">'+json.dumps(audit,ensure_ascii=False).replace('</','<\\/')+'</script><script>'+get_plotlyjs()+'</script><script>'+js+'</script></body></html>'
    d.mkdir(parents=True,exist_ok=False);page.write_text(page_text)
    write_json(d/'html_record.json',{'path':str(page),'sha256':sha256_file(page),'package':str(package),'browser_status':'NOT_YET_CHECKED','human_review':'PENDING'})
    return page


def browser_check(page):
    out=page.parent/'browser_final';record=out/'browser_checks.json'
    if record.exists():
        value=read_json(record)
        if value['status']!='PASS' or value['html_sha256']!=sha256_file(page):raise ValueError('BROWSER_CACHE_MISMATCH')
        return value
    script=Path(__file__).with_name('check_sdpd_diagnostics_browser.cjs')
    paths=[subprocess.check_output(['wslpath','-w',str(p)],text=True).strip() for p in (script,page,out)]
    node=next(p for p in ['/mnt/d/Program Files/nodejs/node.exe','/mnt/c/Program Files/nodejs/node.exe'] if Path(p).exists())
    subprocess.run([node,*paths],check=True,timeout=90)
    value=read_json(record)
    if value['status']!='PASS' or value['html_sha256']!=sha256_file(page):raise ValueError('BROWSER_FAILED')
    return value


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',default='py_scripts/sdpd_diagnostics.yaml')
    parser.add_argument('--open',action='store_true');parser.add_argument('--browser-test',action='store_true');args=parser.parse_args()
    c=load_config(args.config);package=export(c);page=write_page(c,package)
    print('REPORT '+str(package/'report_zh.md'));print('HTML '+str(page),flush=True)
    if args.browser_test:browser_check(page)
    if args.open:open_windows_file(page)


if __name__=='__main__':main()
