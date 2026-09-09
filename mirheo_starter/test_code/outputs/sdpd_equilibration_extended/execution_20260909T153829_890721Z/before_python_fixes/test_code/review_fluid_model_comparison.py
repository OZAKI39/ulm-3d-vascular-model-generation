"""Offline review of recorded runs. This entry point never starts Mirheo or GPU jobs."""
import argparse
import html
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs
from plotly.subplots import make_subplots
from py_scripts.fluid_physics.common import PROJECT_ROOT,output,read_json,write_json,sha256_file,fingerprint
from py_scripts.fluid_physics.calibration import verify_package
from py_scripts.fluid_physics.units import Units
from py_scripts.fluid_comparison.experiments import load_config
from py_scripts.fluid_comparison.reporting import export_results
from test_code.review_vessel_geometry import open_windows_file


def make_figures(summary,case,mapping):
    u=Units(mapping['L0'],mapping['M0'],mapping['t0']);tasks=summary['tasks'];figures={}
    profile=go.Figure();options=[]
    for t in tasks:
        v=t.get('viscosity',{})
        if not v.get('nu_star'):continue
        label=f"{t['method']} / {t['task_id']} ({t['source_category']})";indices=[]
        for values,name,mode,dash in [(v['measured_profile_star'],'模拟点','markers',None),(v['target_profile_star'],'目标解析线','lines','dash'),(v['fit_profile_star'],'拟合实测黏度','lines','solid')]:
            indices.append(len(profile.data));profile.add_trace(go.Scatter(x=np.array(v['y_star'])*u.L0*1e6,y=np.array(values)*u.scales['velocity'],mode=mode,
               name=name,line={'dash':dash} if dash else None,visible=len(options)==0,hovertemplate='y=%{x:.3f} µm<br>u=%{y:.6g} m/s<extra>'+name+'</extra>'))
        options.append({'label':label,'indices':indices,'task_id':t['task_id']})
    preferred=next((o for o in options if o['task_id']=='sdpd_flow'),options[0] if options else {'indices':[]})
    for i,trace in enumerate(profile.data):trace.visible=i in preferred['indices']
    profile.update_layout(xaxis_title='y / µm',yaxis_title='uₓ / m·s⁻¹',title='真实速度剖面：目标预测与拟合结果分开')
    figures['profile']=profile
    temp=go.Figure()
    for t in tasks:
        s=t.get('series',{})
        if not s.get('time_star'):continue
        temp.add_trace(go.Scatter(x=np.array(s['time_star'])*u.t0*1e6,y=s['kBT_thermal_star'],mode='lines',name=t['task_id'],legendgroup=t['candidate_id'],
          visible=True if t['method']=='SDPD' else 'legendonly',hovertemplate='t=%{x:.5g} µs<br>T/T目标=%{y:.6g}<extra>'+t['task_id']+'</extra>'))
    temp.add_hrect(y0=.98,y1=1.02,fillcolor='#15803d',opacity=.15,line_width=0);temp.add_hline(y=1,line_dash='dash')
    temp.update_layout(title='温度：保留初始化与未通过记录；图例可显示历史 DPD',xaxis_title='实际物理时间 / µs',yaxis_title='实测 T / 298.15 K')
    figures['temperature']=temp
    visc=go.Figure()
    for t in tasks:
        v=t.get('viscosity',{})
        if not v.get('nu_si'):continue
        ci=v.get('ci95_si');error={'type':'data','symmetric':False,'array':[ci[1]-v['nu_si']],'arrayminus':[v['nu_si']-ci[0]]} if ci else None
        visc.add_trace(go.Scatter(x=[t['task_id']],y=[v['nu_si']],error_y=error,mode='markers',name=t['task_id'],marker={'size':10},
           hovertemplate='实测 ν=%{y:.6g} m²/s<extra>'+t['task_id']+('，无有效 CI' if ci is None else '')+'</extra>'))
    target=case['kinematic_viscosity_m2_s'];visc.add_hline(y=target,line_dash='dash');visc.add_hrect(y0=.9*target,y1=1.1*target,opacity=.1,fillcolor='green',line_width=0)
    visc.update_layout(title='独立拟合黏度及 95% CI；阴影为冻结的 ±10% 初筛范围',yaxis={'title':'ν / m²·s⁻¹','type':'log'},showlegend=False)
    figures['viscosity']=visc
    eos=go.Figure()
    for a in summary['candidate_assessments']:
        pts=a['eos']['points']
        if not pts:continue
        eos.add_trace(go.Scatter(x=[u.to_si(p['rho_star'],'mass_density') for p in pts],y=[p['mean_pressure_pa'] for p in pts],mode='markers+lines',name=a['candidate_id']+'：机械压力 MEASURED',
           error_y={'type':'data','array':[u.to_si(p['ci95_halfwidth_star'],'pressure') if p['ci95_halfwidth_star'] is not None else None for p in pts]},
           customdata=[p['statistics_status']+(' / CI unavailable' if p['ci95_halfwidth_star'] is None else '') for p in pts],hovertemplate='全盒ρ=%{x:.6g}<br>机械p=%{y:.7g} Pa<br>%{customdata}<extra></extra>'))
    sdpd=[t for t in tasks if t['method']=='SDPD' and t['kind']=='equilibrium']
    if sdpd:
        a=sdpd[0]['parameters']['candidate'];x=np.linspace(7.5,9.5,80)
        eos.add_trace(go.Scatter(x=u.to_si(x,'mass_density'),y=u.to_si(a['sound_speed_star']**2*(x-a['rho_0_star']),'pressure'),name='输入 Linear EOS（INPUT，横轴为公式的质量密度参数）',mode='lines',line={'dash':'dash'}))
        for t in sdpd:
            k=t.get('kernel_density',{})
            if k.get('local_mass_mean') is None:continue
            eos.add_trace(go.Scatter(x=[u.to_si(k['local_mass_mean'],'mass_density')],y=[u.to_si(k['model_EOS_pressure_mean_star'],'pressure')],mode='markers',marker={'symbol':'diamond','size':10},name=t['task_id']+'：局部核密度预测 DERIVED'))
    eos.update_layout(title='EOS 与机械压力：保留离散密度和背景压力差异',xaxis_title='质量密度 / kg·m⁻³（圆点：全盒；菱形：局部核估计）',yaxis_title='压力 / Pa（模型参考，不是测量绝对压力）')
    figures['eos']=eos
    cost=go.Figure()
    for t in tasks:
        p=t.get('performance',{})
        if not p.get('actual_physical_duration_s'):continue
        cost.add_trace(go.Scatter(x=[p['actual_physical_duration_s']*1e6],y=[p['total_wall_s']],mode='markers',name=t['task_id'],marker={'size':11},
            customdata=[[p['total_wall_s_per_step'],p['statistics_physical_duration_s']]],hovertemplate='演化=%{x:.6g} µs<br>墙钟=%{y:.6g} s<br>每步总成本=%{customdata[0]:.6g} s<br>统计时长=%{customdata[1]} s<extra>'+t['source_category']+'</extra>'))
    cost.update_layout(title='实际演化时长与总成本；不同物性设置不作目标液体速度排名',xaxis={'title':'实际物理时长 / µs','type':'log'},yaxis_title='总墙钟 / s')
    figures['cost']=cost
    numerical=go.Figure()
    for t in tasks:
        if t['task_id'] not in ('sdpd_probe','dpd_cost','dpd_cost_density') or t['source_category']!='MEASURED':continue
        p=t['performance'];numerical.add_trace(go.Bar(x=[t['task_id']],y=[1000*p['synchronized_steady_compute_s_per_step']],name=t['task_id']))
    numerical.update_layout(title='匹配短任务的同步计算成本（物性不同，无速度置信区间）',yaxis_title='ms / step',showlegend=False)
    figures['step_cost']=numerical
    memory=go.Figure()
    for t in tasks:
        samples=[s for s in t['execution'].get('resource_samples',[]) if s.get('gpus')]
        if samples:memory.add_trace(go.Scatter(x=[s['elapsed_s'] for s in samples],y=[s['gpus'][0]['used_MiB'] for s in samples],mode='lines+markers',name=t['task_id']))
    memory.update_layout(title='设备总显存采样（包含其他应用，非本任务独占峰值）',xaxis_title='各任务启动后墙钟 / s',yaxis_title='设备已用 / MiB')
    figures['memory']=memory
    budget=go.Figure();b=summary['shared_budget']
    for name,value,color in [('旧任务',b['historical_used_s'],'#64748b'),('本轮新增',b['new_used_s'],'#2563eb'),('剩余',b['remaining_s'],'#d1d5db')]:
        budget.add_trace(go.Bar(x=[value],y=['共享 3600 秒'],name=name,orientation='h',marker_color=color,text=[f'{value:.2f} s'],textposition='auto'))
    budget.update_layout(barmode='stack',title='同一授权预算，跨 campaign 累计',xaxis_title='墙钟 / s',height=240)
    figures['budget']=budget
    for f in figures.values():
        f.update_layout(template='plotly_white',font={'family':'Arial, Microsoft YaHei, sans-serif'},margin={'l':70,'r':30,'b':90,'t':65})
    return figures,options


def write_page(c,package):
    verify_package(package);package=Path(package);summary=read_json(package/'comparison_summary.json')
    case=read_json(package/'physical_case.json');mapping=read_json(package/'unit_mapping.json')
    checker=PROJECT_ROOT/'test_code/check_fluid_model_comparison_browser.cjs'
    key=fingerprint({'package':sha256_file(package/'package_sha256.json'),'viewer':sha256_file(Path(__file__)),'checker':sha256_file(checker)})[:12]
    d=output(Path(c['review_root'])/package.name/key);page=d/'comparison_review.html'
    if page.exists():
        record=read_json(d/'html_record.json')
        if record['sha256']!=sha256_file(page):raise ValueError('EXISTING_HTML_MODIFIED')
        return page
    d.mkdir(parents=True,exist_ok=False);figs,options=make_figures(summary,case,mapping)
    real=[a for a in summary['candidate_assessments'] if a['status']!='COST_ONLY']
    def fmt(value):return '无有效结果' if value is None else f'{value:.6g}'
    status_text={'NOT_QUALIFIED':'未合格','COVERED':'覆盖通过初筛','NOMINAL_COVERAGE_UNCERTAIN':'均值名义覆盖，计入不确定性未通过',
                 'OUTSIDE_MEASURED_EOS_RANGE':'超出已测压力范围','UNVERIFIED':'未验证'}
    gate_text={'planned_execution_complete':'计划执行不完整','temperature':'温度','stationarity':'稳态',
               'viscosity':'黏度及统计','step_and_force_sensitivity':'时间步/驱动力敏感性','EOS_statistics':'EOS 统计',
               'EOS_temperature':'EOS 温度','pressure_range':'压力范围','density_adjustment':'所需密度变化',
               'measured_EOS_Mach':'基于实测 EOS 的马赫数'}
    def val(a,k):
        v=a.get('base_viscosity') or {}
        if k=='nu':
            ci=v.get('ci95_si');uncertainty='无有效 95% CI' if ci is None else '95% CI ['+', '.join(fmt(x) for x in ci)+']'
            return fmt(v.get('nu_si'))+'<br>'+uncertainty
        if k=='mu_input':return '不适用' if v.get('input_mu_star') is None else fmt(v['input_mu_star'])
        if k=='error':return f"{100*v['target_relative_error']:.2f}%" if v.get('target_relative_error') is not None else '未测'
        if k=='T':return ', '.join(fmt(x) for x in a['equilibrium_temperature_star'])+' / '+fmt(a['flow_temperature_star'])
        if k=='samples':return str(v.get('block_valid_count',0))+' 个有效块；'+('充分' if v.get('sampling_status')=='SUFFICIENT' else '统计不足')
        if k=='density':
            values=[t['density_si'] for t in summary['tasks'] if t['method']==a['method'] and t['candidate_id']==a['candidate_id'] and t['test_kind']=='equilibrium' and 'density_si' in t]
            return ', '.join(fmt(x) for x in values) or '未测'
        if k=='EOS':
            code=a['eos'].get('target_pressure_coverage_status','UNVERIFIED');return html.escape(status_text.get(code,code))
        if k=='status':return status_text.get(a['status'],a['status'])
        return html.escape('、'.join(gate_text.get(k,k) for k in a['failed_or_unverified']))
    rows=''.join('<tr><th>'+name+'</th><td>'+target+'</td>'+''.join('<td>'+val(a,k)+'</td>' for a in real)+'</tr>' for k,name,target in [
      ('density','全盒质量密度 / kg·m⁻³',fmt(case['density_kg_m3'])),('nu','实测 ν / m²·s⁻¹',fmt(case['kinematic_viscosity_m2_s'])),('mu_input','SDPD 输入 μ*','3312.4601269（精确值见 JSON）'),('error','黏度绝对相对误差','±10% PROPOSED，另需 CI'),('samples','黏度统计','至少 8 个块及相关时间要求'),
      ('T','平衡 / 流动 T*','1 ±2%，含不确定性'),('EOS','机械压力覆盖','保留旧三个表压差'),('status','材料状态','独立通过全部初筛'),('failed','未通过或未验证','不得借用其他候选证据')])
    table='<table id="summary-table"><thead><tr><th>指标</th><th>目标/定义</th>'+''.join('<th>'+html.escape(a['method']+' / '+a['candidate_id'])+'</th>' for a in real)+'</tr></thead><tbody>'+rows+'</tbody></table>'
    graphs=''.join('<section id="section-'+name+'"><div class="plot" id="'+name+'"></div></section>' for name in figs)
    payload=json.dumps({name:json.loads(f.to_json()) for name,f in figs.items()},ensure_ascii=False).replace('</','<\\/')
    opts=json.dumps(options,ensure_ascii=False).replace('</','<\\/')
    audit={'package':str(package),'package_manifest_sha256':sha256_file(package/'package_sha256.json'),'task_count':len(summary['tasks']),
           'source_counts':{s:sum(t['source_category']==s for t in summary['tasks']) for s in ['MEASURED','HISTORICAL_REFERENCE']},
           'selection':summary['selection'],'human_review':'PENDING','figure_ids':list(figs),'profile_options':options,'checker_sha256':sha256_file(checker)}
    js="""
const figures=FIGURES,options=OPTIONS;
window.comparisonReady=false;
Promise.all(Object.entries(figures).map(([id,f])=>Plotly.newPlot(id,f.data,f.layout,{responsive:true,displaylogo:false,scrollZoom:true}))).then(()=>{window.comparisonReady=true;});
const select=document.getElementById('profile-select');
options.forEach((o,i)=>{let e=document.createElement('option');e.value=i;e.textContent=o.label;select.appendChild(e);});
select.value=String(Math.max(0,options.findIndex(o=>o.task_id==='sdpd_flow')));
if(options.length)document.getElementById('profile-caption').textContent=options[+select.value].label;
select.addEventListener('change',()=>{let visible=Array(figures.profile.data.length).fill(false);options[+select.value].indices.forEach(i=>visible[i]=true);Plotly.restyle('profile',{visible});document.getElementById('profile-caption').textContent=options[+select.value].label;});
document.getElementById('reset-plots').addEventListener('click',()=>Object.keys(figures).forEach(id=>Plotly.relayout(id,{'xaxis.autorange':true,'yaxis.autorange':true})));
""".replace('FIGURES',payload).replace('OPTIONS',opts)
    text='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>DPD / SDPD 纯液体对照核查</title><style>
body{font:16px/1.65 Arial,"Microsoft YaHei",sans-serif;background:#f4f6f9;color:#15233b;margin:0}main{max-width:1380px;margin:auto;padding:28px}h1{margin-bottom:6px}h2{margin-top:30px}.note{background:#fff4dc;border-left:5px solid #a16207;padding:14px}.plot{height:470px}#budget{height:260px}section,.table-wrap{background:white;margin:22px 0;padding:12px;border-radius:8px;overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{border:1px solid #dce2ea;padding:9px;vertical-align:top}th{background:#eef3fa}select,button{font-size:15px;padding:8px;max-width:100%;margin:8px}code{overflow-wrap:anywhere}.meta{font-size:13px;color:#475569}a{color:#1d4ed8}</style></head><body><main>
<h1>DPD 与原生 SDPD：纯液体对照</h1><p>结果状态：<strong>STATUS</strong> · 人工验收：<strong>PENDING</strong> · selection：SELECTION</p>
<div class="note">输入 viscosity 与实测黏度分开。历史 DPD 与本轮 SDPD 使用原保存的同一单位；物性未同时合格时不作目标液体速度排名。初始化、失败与统计不足的记录均保留。</div>
<div class="table-wrap">TABLE</div><h2>真实数据与可交互图</h2><p>剖面选择：<select id="profile-select"></select><button id="reset-plots">恢复图表范围</button></p><p id="profile-caption">选择实验可切换模拟点、目标解析线和拟合线。</p>
GRAPHS
<h2>范围与核查记录</h2><p>圆形 EOS 点为独立机械压力；菱形点为观测局部核密度代入 EOS 的预测。图中缺少误差条表示该记录没有有效 CI，不能视为零误差。温度和黏度图保留不足统计长度的实际数据。温度图默认显示本轮 SDPD；历史 DPD 的物理时长较长，可通过图例单独查看。</p>
<p>每个任务独立初始化。物理演化时长和有效统计时长另列在报告，未将多条短轨迹拼接。设备显存包含其他应用。未进行真实血管、壁面、RBC、微泡、SDF 或开放边界仿真；本页面不启动 CUDA 或 GPU 实验。</p>
<p class="meta">报告目录：<code>PACKAGE</code><br>精确输入与原始运行哈希见 provenance.json、comparison_summary.json 和 package_sha256.json。浏览器通过状态只在实际运行检查后记录到 browser_final/browser_checks.json；人工验收始终为 PENDING。</p>
</main><script id="comparison-audit-data" type="application/json">AUDIT</script><script>PLOTLY</script><script>JS</script></body></html>'''
    for key,value in [('STATUS',summary['status']),('SELECTION',html.escape(json.dumps(summary['selection'],ensure_ascii=False))),('TABLE',table),('GRAPHS',graphs),('PACKAGE',html.escape(str(package))),('AUDIT',json.dumps(audit,ensure_ascii=False)),('<script>PLOTLY</script>','<script>'+get_plotlyjs()+'</script>'),('<script>JS</script>','<script>'+js+'</script>')]:
        text=text.replace(key,value)
    with page.open('x') as f:f.write(text)
    write_json(d/'html_record.json',{'path':str(page),'sha256':sha256_file(page),'package':str(package),'status':'GENERATED','browser_status':'NOT_YET_CHECKED','human_review':'PENDING'})
    return page


def browser_check(page):
    out=page.parent/'browser_final'
    if (out/'browser_checks.json').exists():
        record=read_json(out/'browser_checks.json')
        if record['html_sha256']!=sha256_file(page) or record['status']!='PASS':raise ValueError('EXISTING_BROWSER_RECORD_NOT_VALID')
        return record
    script=PROJECT_ROOT/'test_code/check_fluid_model_comparison_browser.cjs'
    paths=[subprocess.check_output(['wslpath','-w',str(p)],text=True).strip() for p in [script,page,out]]
    node=next((p for p in ['/mnt/d/Program Files/nodejs/node.exe','/mnt/c/Program Files/nodejs/node.exe'] if Path(p).is_file()),None)
    if node is None:raise FileNotFoundError('Existing Windows Node executable not found')
    subprocess.run([node,*paths],check=True,timeout=90)
    record=read_json(out/'browser_checks.json')
    if record['html_sha256']!=sha256_file(page) or record['status']!='PASS':raise ValueError('BROWSER_CHECK_FAILED')
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',default='py_scripts/fluid_model_comparison.yaml')
    parser.add_argument('--open',action='store_true');parser.add_argument('--browser-test',action='store_true');args=parser.parse_args()
    c=load_config(args.config);package=export_results(c);page=write_page(c,package)
    print('REPORT '+str(package/'report_zh.md'));print('HTML '+str(page))
    if args.browser_test:browser_check(page)
    if args.open:open_windows_file(page)


if __name__=='__main__':main()
