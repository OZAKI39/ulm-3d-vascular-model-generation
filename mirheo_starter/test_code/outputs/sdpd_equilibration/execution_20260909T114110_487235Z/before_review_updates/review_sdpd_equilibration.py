"""Offline view of one frozen equilibration package; never invokes execute."""
import argparse
import html
import json
import os
from pathlib import Path
import subprocess
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.offline import get_plotlyjs
from py_scripts.fluid_physics.common import read_json, write_json, sha256_file, fingerprint, output
from py_scripts.fluid_physics.calibration import verify_package
from py_scripts.fluid_physics.analysis import read_csv
from py_scripts.sdpd_diagnostics.equilibration import load_config, export, read_raw
from test_code.review_sdpd_diagnostics import table, fmt
from test_code.review_vessel_geometry import open_windows_file


def make_figures(package, plan, summary):
    scale = plan['parameters']['locked_units']['t0']*1e6
    target = plan['parameters']['kBT_star']
    lo, hi = [x*scale for x in plan['formal_window_star']]
    actual = summary['actual_steps'] is not None
    data = read_raw(package) if actual and summary['actual_steps'] else None
    historical = read_csv(Path(plan['cost_evidence']['directory'])/'moments.csv')
    full = go.Figure()
    if data:
        full.add_trace(go.Scatter(x=data['physical_time_s']*1e6,y=data['kBT_COM_star']/target,
                                 mode='lines',name='本次连续演化',line={'color':'#087e8b'}))
    full.add_trace(go.Scatter(x=historical['time_star']*scale,y=historical['kBT_COM_star']/target,
                             mode='lines',name='历史参考：旧无驱动任务（非本轮结果）',
                             visible='legendonly' if data else True,line={'color':'#8b6f47','dash':'dot'}))
    full.update_layout(title='完整温度轨迹：启动峰值保留'+('' if data else ' · 本轮尚未运行'),
                       xaxis_title='物理时间 / µs',yaxis_title='T / 298.15 K',xaxis={'range':[0,hi]})
    late=go.Figure()
    if data:
        late.add_trace(go.Scatter(x=data['physical_time_s']*1e6,y=data['kBT_COM_star']/target,name='本次温度',mode='lines'))
    else:
        late.add_annotation(text='本轮待授权：正式窗口无新测量数据',showarrow=False,x=.5,y=.5,xref='paper',yref='paper')
    late.update_layout(title='固定正式窗口放大（不会替换全程图）',xaxis={'range':[lo,hi]},
                       xaxis_title='物理时间 / µs',yaxis_title='T / 298.15 K')
    for f in [full,late]:
        f.add_hrect(y0=.98,y1=1.02,fillcolor='#6dbb75',opacity=.22,line_width=0)
        f.add_hline(y=1,line_dash='dash',line_color='#236b36')
        f.add_vrect(x0=lo,x1=hi,fillcolor='#92bfe8',opacity=.14,line_width=0)
    blocks=go.Figure()
    stats=read_json(package/'block_statistics.json')
    if stats:
        t=stats['kBT_COM_star']
        blocks.add_trace(go.Scatter(x=np.array(t['block_times_star'])*scale,y=t['block_means'],mode='markers',name='完整块均值'))
        if t['ci95_halfwidth'] is not None:
            blocks.add_hrect(y0=t['mean']-t['ci95_halfwidth'],y1=t['mean']+t['ci95_halfwidth'],fillcolor='#087e8b',opacity=.2,line_width=0)
            blocks.add_hline(y=t['mean'],line_color='#087e8b')
    else:blocks.add_annotation(text='无新完整块；不生成 CI',showarrow=False)
    blocks.update_layout(title='完整块均值与整体均值 CI（非每块独立 CI）',xaxis_title='物理时间 / µs',yaxis_title='kBT*',xaxis={'range':[lo,hi]})
    evolution=make_subplots(rows=4,cols=1,shared_xaxes=True,vertical_spacing=.07,
                            subplot_titles=['机械压力','局部核密度均值与分位数','局部核密度方差','平均速度 / COM'])
    if data:
        time=data['physical_time_s']*1e6
        for row,key,label in [(1,'pressure_pa','机械压力 / Pa'),(2,'kernel_number_mean','核密度均值'),
                               (2,'kernel_q10','核密度 q10'),(2,'kernel_q90','核密度 q90'),
                               (3,'kernel_number_variance','核密度方差'),(4,'mean_vx','平均 vx'),(4,'mean_vy','平均 vy'),(4,'mean_vz','平均 vz')]:
            evolution.add_trace(go.Scatter(x=time,y=data[key],name=label,mode='lines'),row=row,col=1)
    else:
        for row in range(1,5):evolution.add_annotation(text='本轮未测',showarrow=False,row=row,col=1)
    evolution.update_yaxes(title_text='Pa',row=1,col=1)
    evolution.update_yaxes(title_text='核数密度*',row=2,col=1)
    evolution.update_yaxes(title_text='方差*',row=3,col=1)
    evolution.update_yaxes(title_text='速度*',row=4,col=1)
    evolution.update_xaxes(title_text='物理时间 / µs',row=4,col=1)
    evolution.update_layout(title='压力、局部结构与动量分别核查',height=880)
    figures={'temperature':full,'late':late,'blocks':blocks,'evolution':evolution}
    for f in figures.values():
        f.update_layout(template='plotly_white',font={'family':'Arial, Microsoft YaHei, sans-serif'},
                        showlegend=True,margin={'l':85,'r':30,'t':85,'b':100},legend={'orientation':'h','y':-.25})
    return figures


def write_page(c, package):
    package=Path(package);verify_package(package)
    plan=read_json(package/'run_plan.json');summary=read_json(package/'equilibration_summary.json')
    request=read_json(package/'budget_request.json');audit=read_json(package/'sampling_audit.json')
    checker=Path(__file__).with_name('check_sdpd_equilibration_browser.cjs')
    harness=Path(__file__).with_name('check_sdpd_diagnostics_browser.cjs')
    identity={'package':sha256_file(package/'package_sha256.json'),'viewer':sha256_file(Path(__file__)),
              'checker':sha256_file(checker),'harness':sha256_file(harness)}
    d=output(Path(c['review_root'])/package.name/fingerprint(identity)[:12])
    page=d/'sdpd_equilibration_review.html'
    if page.exists():
        if sha256_file(page)!=read_json(d/'html_record.json')['sha256']:raise ValueError('HTML_CHANGED')
        return page
    d.mkdir(parents=True,exist_ok=False)
    figures=make_figures(package,plan,summary)
    cards=table(['判断项','本轮结论'],[
        ['程序完成',summary['program_completed']],['液体已稳定',summary['liquid_stationary']],
        ['稳定温度匹配',summary['stable_temperature_matches_target']],['全部物性合格','未验证；selection=null']],id='outcomes')
    budget=table(['预算项','秒'],[['已有累计用量',request['total_charged_or_reserved_s']],['原授权余额',request['original_remaining_s']],
        ['本次预计耗时',request['expected_wall_s']],['建议单任务上限',request['suggested_task_limit_s']],
        ['需追加精确值',request['additional_required_s']],['建议追加整数值',request['additional_request_whole_seconds']]],id='budget-table')
    sampling=table(['统计项','数值'],[[k,audit.get(k)] for k in ['total_samples','formal_samples','block_samples','block_count','discarded_tail_samples','CI_basis']],id='sampling-table')
    text=('<h1>无驱动 SDPD：持续演化与稳定性</h1><p class="status">'+html.escape(summary['status'])+'</p>'
          '<p id="measurement-state">'+html.escape(summary['execution_status'])+'</p>'
          '<p>固定液体参数 · 单次连续演化 · 人工验收 PENDING</p>'+cards+
          '<p><strong>未运行时所有新测量为空。历史参考曲线不是本次结果。</strong></p>'
          '<p><a href="'+html.escape(os.path.relpath(package/'report_zh.md',d))+'">中文报告</a> · '
          '<a href="'+html.escape(os.path.relpath(package/'run_plan.json',d))+'">冻结计划</a> · <button id="reset-plots">恢复图表范围</button></p>'
          '<h2>完整启动曲线与固定正式窗口</h2><p>计划 400000 步，dt*=1e-6；正式窗口 t*=0.30–0.40。绿色带为目标 ±2%，蓝色带为正式观察范围。</p>'
          '<div class="plot" id="temperature"></div><div class="plot" id="late"></div>'
          '<h2>稳定性与温度匹配</h2><p>'+html.escape(summary['reason'])+'</p>'+sampling+
          '<div class="plot" id="blocks"></div><div class="plot tall" id="evolution"></div>'
          '<p>压力采用原生应力迹与 COM 扣除后的动能项；密度/应力来自积分前，速度来自积分后，相位差 dt 保留。压力稳定不代表出口目标或 EOS 验证通过。</p>'
          '<h2>原计划、实际完成量与预算</h2>'+budget+
          table(['项目','值'],[['计划物理时长 / µs',plan['physical_duration_s']*1e6],['实际步数',summary['actual_steps']],
                              ['实际物理时长 / s',summary['actual_physical_time_s']],['完成/停止状态',summary['execution_status']],
                              ['预算授权','当前未追加' if request['extra_authorized_gpu_seconds']==0 else request['extra_authorized_gpu_seconds']]])+
          '<p>0.30* 不保证已经平衡；8 个有效块和温度通过都不保证获得。预算不足不启动缩短实验，运行偏慢则按上限安全停止并标 PARTIAL。</p>'
          '<p>本轮生产参数 selection=null；不自动启动黏度测试、EOS 扫描或其他求解器。</p>')
    payload={'status':summary['status'],'execution_status':summary['execution_status'],'human_review':'PENDING',
             'planned_steps':plan['steps'],'formal_window_star':plan['formal_window_star'],
             'actual_steps':summary['actual_steps'],'new_temperature_mean':summary['temperature_mean_star'],
             'stationary':summary['liquid_stationary'],'temperature_match':summary['stable_temperature_matches_target'],
             'selection':None,'extra_authorized_gpu_seconds':request['extra_authorized_gpu_seconds'],
             'raw_csv_rows':sum(1 for _ in (package/'raw_statistics.csv').open())-1}
    plots={name:json.loads(f.to_json()) for name,f in figures.items()}
    javascript=('const figures='+json.dumps(plots,ensure_ascii=False).replace('</',r'<\/')+';'
        'window.equilibrationReady=false;'
        'Promise.all(Object.entries(figures).map(([id,f])=>Plotly.newPlot(id,f.data,f.layout,{responsive:true,scrollZoom:true,displaylogo:false})))'
        '.then(()=>{window.equilibrationReady=true});'
        'document.getElementById("reset-plots").onclick=()=>Promise.all(Object.keys(figures).map(id=>Plotly.relayout(id,{"xaxis.autorange":true,"yaxis.autorange":true})));')
    page.write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
        '<title>SDPD 无驱动稳定性</title><style>body{font-family:Arial,"Microsoft YaHei",sans-serif;color:#18323d;background:#f3f6f7;margin:0}main{max-width:1260px;margin:auto;padding:32px}h1{font-size:29px}h2{margin-top:34px}.status{font-weight:bold;color:#885419}p{line-height:1.7}.plot{background:white;margin:20px 0;min-height:470px;border:1px solid #dce5e8;border-radius:8px}.tall{min-height:880px}table{border-collapse:collapse;width:100%;background:white}td,th{border:1px solid #dce5e8;padding:10px;text-align:left;overflow-wrap:anywhere}th{background:#e7eff2}.table-wrap{overflow-x:auto;margin:15px 0}button{padding:8px 14px;cursor:pointer}</style>'
        '<main>'+text+'</main><script id="equilibration-audit-data" type="application/json">'+json.dumps(payload,ensure_ascii=False).replace('</',r'<\/')+
        '</script><script>'+get_plotlyjs()+'</script><script>'+javascript+'</script></html>')
    write_json(d/'html_record.json',{'sha256':sha256_file(page),'package':str(package),'identity':identity,'human_review':'PENDING'})
    return page


def browser_check(page):
    d=page.parent/'browser_final'
    if (d/'browser_checks.json').exists():
        r=read_json(d/'browser_checks.json')
        if r['status']!='PASS' or r['html_sha256']!=sha256_file(page):raise ValueError('BROWSER_CACHE_MISMATCH')
        return r
    paths=[subprocess.check_output(['wslpath','-w',str(p)],text=True).strip() for p in
           [Path(__file__).with_name('check_sdpd_equilibration_browser.cjs'),page,d]]
    node=next(p for p in ['/mnt/d/Program Files/nodejs/node.exe','/mnt/c/Program Files/nodejs/node.exe'] if Path(p).exists())
    subprocess.run([node,*paths],check=True,timeout=90)
    r=read_json(d/'browser_checks.json')
    if r['status']!='PASS' or r['html_sha256']!=sha256_file(page):raise ValueError('BROWSER_FAILED')
    return r


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',default='py_scripts/sdpd_equilibration.yaml')
    p.add_argument('--open',action='store_true');p.add_argument('--browser-test',action='store_true')
    args=p.parse_args();c=load_config(args.config);package=export(c,analyze=True);page=write_page(c,package)
    print('REPORT '+str(package/'report_zh.md'));print('HTML '+str(page),flush=True)
    if args.browser_test:browser_check(page)
    if args.open:open_windows_file(page)


if __name__=='__main__':main()
