"""Offline actual-data review. This module never launches a GPU worker."""
import argparse
from html import escape
import json
from pathlib import Path
import numpy as np
from py_scripts.fluid_physics.common import PROJECT_ROOT,config_read,output,read_json,write_json,sha256_file
from py_scripts.fluid_physics.reporting import export_results
from py_scripts.fluid_physics.calibration import verify_package
from py_scripts.vessel_geometry.export import load_package
from test_code.review_vessel_geometry import make_figure,open_windows_file


def jsdata(obj):return json.dumps(obj,ensure_ascii=False,allow_nan=False).replace('<','\\u003c')


def table(headers,rows):
    return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+escape(str(x))+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+escape(str(x))+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table></div>'


def real_figures(summary,unit,case):
    import plotly.graph_objects as go
    figs={k:go.Figure() for k in ('temperature','profile','viscosity','eos','sensitivity')}
    tasks=summary['executed_tasks'];scale=unit['si_per_star'];eos=summary['eos']
    for task in tasks:
        if 'series' in task:
            s=task['series'];t=np.array(s['time_star'])*unit['t0']*1000
            figs['temperature'].add_trace(go.Scatter(x=t.tolist(),y=(np.array(s['kBT_thermal_star'])*case['temperature_K']).tolist(),mode='lines',name=task['task_id']+' 实测热速度',line={'width':1},legendgroup=task['task_id']))
        v=task.get('viscosity',{})
        if v.get('nu_star'):
            y=(np.array(v['y_star'])*unit['L0']*1e6).tolist()
            figs['profile'].add_trace(go.Scatter(x=y,y=(np.array(v['measured_profile_star'])*scale['velocity']*1000).tolist(),mode='markers',name=task['task_id']+' 实测均值'))
            figs['profile'].add_trace(go.Scatter(x=y,y=(np.array(v['fit_profile_star'])*scale['velocity']*1000).tolist(),mode='lines',line={'dash':'dash'},name=task['task_id']+' 拟合解析式'))
            vals=v.get('block_nu_star',[]);times=v.get('block_times_star',[])
            figs['viscosity'].add_trace(go.Scatter(x=(np.array(times)*unit['t0']*1000).tolist(),y=[x*scale['kinematic_viscosity'] if x is not None else None for x in vals],mode='markers',name=task['task_id']+' 时间块'))
            ci=v.get('ci95_si');nu=v['nu_si']
            figs['sensitivity'].add_trace(go.Scatter(x=[task['task_id']],y=[nu],mode='markers',marker={'size':10},name=task['task_id'],error_y={'type':'data','symmetric':False,'array':[ci[1]-nu] if ci else [0],'arrayminus':[nu-ci[0]] if ci else [0],'visible':bool(ci)}))
    figs['temperature'].add_hline(y=case['temperature_K'],line_dash='dash',line_color='black',annotation_text='298.15 K 用户假设')
    figs['temperature'].add_hrect(y0=case['temperature_K']*.98,y1=case['temperature_K']*1.02,line_width=0,fillcolor='#20a880',opacity=.08,annotation_text='±2% PROPOSED')
    # Annotation positions on a Plotly log axis use log coordinates. Put labels
    # in paper coordinates so they cannot stretch the numerical axis up to 1.
    for key in ('viscosity','sensitivity'):
        figs[key].add_hline(y=case['kinematic_viscosity_m2_s'],line_dash='dash')
        figs[key].add_annotation(x=.99,y=.99,xref='paper',yref='paper',text='旧目标 3.27e-6 m²/s',showarrow=False,xanchor='right')
    points=eos.get('points',[])
    if points:
        figs['eos'].add_trace(go.Scatter(x=[p['rho_star']*scale['mass_density'] for p in points],y=[p['mean_pressure_pa'] for p in points],mode='markers',marker={'size':12},name='原生 virial + 动能项 / V 实测',error_y={'type':'data','array':[(p['ci95_halfwidth_star'] or 0)*scale['pressure'] for p in points],'visible':True}))
        if eos.get('fit'):
            figs['eos'].add_trace(go.Scatter(x=[p['rho_star']*scale['mass_density'] for p in points],y=[p*scale['pressure'] for p in eos['fit']['fitted_pressure_star']],mode='lines',line={'dash':'dash'},name='仅测量区间内局部拟合'))
    labels={
        'temperature':('实际热速度温度（流动扣除瞬时 y 分箱均速）','物理时间 / ms','T / K'),
        'profile':('周期双向 Poiseuille：实测点与拟合曲线','y / µm','uₓ / mm s⁻¹'),
        'viscosity':('独立时间块的黏度估计（合成数据不进入本页）','物理时间 / ms','ν / m² s⁻¹'),
        'eos':('实际压力—密度响应（未外推至出口目标）','ρ / kg m⁻³','p / Pa，新 DPD 平衡参考'),
        'sensitivity':('同一锁定单位下的驱动力与时间步比较','实际任务','ν / m² s⁻¹；95% CI')}
    for key,fig in figs.items():
        title,x,y=labels[key];fig.update_layout(template='plotly_white',title=title,xaxis_title=x,yaxis_title=y,height=480,margin={'l':70,'r':30,'t':70,'b':90},legend={'orientation':'h','y':-0.23})
        if not fig.data:fig.add_annotation(text='未执行 / 无可用实测数据',x=.5,y=.5,xref='paper',yref='paper',showarrow=False)
    figs['viscosity'].update_yaxes(type='log');figs['sensitivity'].update_yaxes(type='log')
    return {k:json.loads(fig.to_json()) for k,fig in figs.items()}


def write_html(package,destination,browser_status='NOT_TESTED'):
    from plotly.offline import get_plotlyjs
    package=Path(package);verify_package(package)
    case=read_json(package/'physical_case.json');unit=read_json(package/'unit_mapping.json');summary=read_json(package/'calibration_summary.json');boundary=read_json(package/'boundary_plan.json');prov=read_json(package/'provenance.json')
    geometry=load_package(Path(case['geometry_package']));fig,controls=make_figure(geometry)
    # Add only the actual frozen finite measurement contours (no invented buffers).
    for port in case['ports']:
        p=port['measurement_planes']['central'];uv=np.array(p['physical_aperture_contour_uv_m']);xyz=np.array(p['origin_m'])+uv[:,0,None]*np.array(p['basis_u'])+uv[:,1,None]*np.array(p['basis_v']);xyz=np.vstack([xyz,xyz[0]])*1e6
        index=len(fig['data']);fig['data'].append({'type':'scatter3d','mode':'lines','x':xyz[:,0].tolist(),'y':xyz[:,1].tolist(),'z':xyz[:,2].tolist(),'line':{'color':'#e36522','width':5},'name':port['name']+' 旧实测截面','legendgroup':str(port['patch']['entity_id'])})
        controls['ports'][str(port['patch']['entity_id'])].append(index)
    charts=real_figures(summary,unit,case)
    rows=[]
    names={'density_kg_m3':'密度','kinematic_viscosity_m2_s':'运动黏度','dynamic_viscosity_pa_s':'动力黏度 μ=ρν','bulk_viscosity_m2_s':'旧运动体积黏度（未独立匹配）','target_volume_flow_m3_s':'入口总体积流量','target_mass_flow_kg_s':'入口总质量流量','outlet_gauge_pressures_pa':'出口表压（负值保留）','pressure_numerical_offset_pa':'旧 LBM 数值压力偏置','temperature_K':'温度（25°C 建模假设）'}
    for key,v in case['trace'].items():rows.append([names.get(key,key),v['value'],v['unit'],v['category'],v['file']+' # '+v['field']])
    physics=table(['物理量','值','单位','属性','本地来源'],rows)
    units=table(['映射','数值','说明'],[
        ['1 模拟长度',f"{unit['L0']:.9g} m",'DESIGN_CHOICE；不是旧 0.20 µm 网格'],['1 模拟质量',f"{unit['M0']:.9g} kg",'由旧密度与 n* m* 确定'],['1 模拟时间',f"{unit['t0']:.9g} s",'由 298.15 K 与 kBT*=1 确定'],['1 模拟压力',f"{unit['si_per_star']['pressure']:.9g} Pa",'E0 / L0³'],['1 模拟能量',f"{unit['E0']:.9g} J",'kB × 298.15 K'],['目标 ν*',unit['required_nu_star'],'实测必须满足；不得逐任务重新调 t0'],['粒子平均间距 / rc',f"{unit['particle_spacing_m']*1e6:g} / {unit['interaction_cutoff_m']*1e6:g} µm",'未来 SDF 间距仍未确定'],['2/4 µm 相对间距',unit['bubble_to_spacing_ratio'],'仅分辨率约束；本轮没有微泡或 RBC']])
    taskrows=[];fitrows=[]
    for task in summary['executed_tasks']:
        e=task['execution'];taskrows.append([task['task_id'],task['parameters']['steps'],f"{e['elapsed_monotonic_s']:.3f}",e['exit_code'],e['timeout'],e['device_sampled_peak_used_MiB'],task.get('actual_N','—'),task.get('temperature_status',task['status'])])
        v=task.get('viscosity',{})
        if v.get('nu_star'):fitrows.append([task['task_id'],v['nu_star'],v['nu_si'],v.get('ci95_si'),v['relative_rms'],v.get('sampling_status'),v.get('target_matching_status'),task.get('stationarity_status')])
    resource=table(['任务','安排步数','墙钟 s','退出码','超时','显存采样峰值 MiB','实际粒子数','温度/统计状态'],taskrows)
    fits=table(['任务','ν*','ν / m² s⁻¹','95% CI / m² s⁻¹','剖面相对 RMS','采样','目标匹配','稳态'],fitrows)
    statuses=table(['项目','本轮状态'],summary['states'].items())
    options='<option value="">全部端口</option>'+''.join(f'<option value="{p["entity_id"]}">{escape(p["name"])}</option>' for p in boundary['ports'])
    plans=''
    for p in boundary['ports']:
        plans+=f'<details data-entity="{p["entity_id"]}" open><summary>{escape(p["name"])} · {escape(p["state"])}</summary><p><code>{escape(p["port_id"])}</code></p>'
        plans+=table(['要求','值'],[['物理要求（测试值）',p['physical_requirement']],['外法向 / 流量符号',str(p['outward_normal'])+' / '+p['flux_sign']],['端盖位置 m',p['cap_center_m']],['保留测量截面位置 m',p['measurement_plane']['origin_m']],['测量截面距端盖 µm',p['measurement_to_cap_axial_distance_m']*1e6],['已存在延伸段长度 µm',p['existing_extension']['actual_axial_length_um']],['控制区域设计',p['finite_control_region']],['可复用原生接口',p['existing_apis']],['下一阶段扩展',p['extension']],['下一阶段验证',p['next_validation']]])+'</details>'
    payload={'summary_status':summary['status'],'states':summary['states'],'figure':fig,'controls':controls,'charts':charts,'geometry_triangle_count':len(geometry.triangles),'package':str(package),'case_type':case['case_type'],'actual_task_ids':[t['task_id'] for t in summary['executed_tasks']]}
    text="""<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>旧纯流体工况 · DPD 第二阶段核查</title>
<style>*{box-sizing:border-box}body{margin:0;font:15px/1.65 system-ui,'Microsoft YaHei',sans-serif;background:#eef3f7;color:#173142}header{padding:30px 4%;background:#123247;color:white}h1{font-size:27px;margin:0 0 9px}h2{font-size:21px}main{max-width:1540px;margin:auto;padding:22px}section{background:white;padding:24px;margin:0 0 22px;border-radius:10px;border:1px solid #d5e0e7}.banner{padding:14px;background:#fff2d8;border-left:5px solid #c98118;color:#593708}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{border:1px solid #d3dfe7;padding:9px;text-align:left;vertical-align:top;overflow-wrap:anywhere}th{background:#edf4f8}td{min-width:90px;max-width:570px}details{margin:12px 0;border:1px solid #dae3ea;padding:12px}summary{cursor:pointer;font-weight:bold}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px}code{overflow-wrap:anywhere}select,button{padding:7px;font:inherit}.controls{display:flex;gap:18px;align-items:center;flex-wrap:wrap}.plot{min-height:480px}.selected{background:#fff0bf}#vessel-view{height:740px}.tag{display:inline-block;background:#dce9f2;color:#153a53;margin:4px;padding:5px 10px;border-radius:5px}a{color:#16718e}</style></head><body>
<header><h1>旧纯流体工况 · DPD 第二阶段核查</h1><div>真实本地小盒数据 · 旧边界值为测试值 · 25°C 为用户同意的假设 · 不代表真实小鼠血浆</div><span class=tag>总体 __STATUS__</span><span class=tag>USER_REVIEW = PENDING</span><span class=tag>浏览器 __BROWSER__</span></header><main>
<section><div class=banner>本轮没有合格参数时 selection=null。代码完成、GPU 执行、统计充分、物性匹配和边界实现分别判定。页面只读取数据，不启动 GPU。</div>__STATES__</section>
<section><h2>1. 可追溯物理工况</h2>__PHYSICS__<p>旧工况为恒定边界、稳定后的平均流场。入口限制总流量，未指定完整速度剖面。旧数值压力基线不是真实热力学绝对压力；负表压不被截为零。</p></section>
<section><h2>2. 锁定单位与分辨率约束</h2>__UNITS__<p>单位满足设定密度和热能约束；材料黏度仍须测量。实测温度偏离设定值会单独失败，不能凭 kBT*=1 宣布已达到 298.15 K。</p><details><summary>分量资源估算（未初始化真实血管）</summary><pre>__ESTIMATES__</pre></details></section>
<section><h2>3. 实际标定与统计限制</h2><p>点和时间序列来自本轮原始 CSV；虚线为拟合或来源目标。置信区间使用相关时间约束下的时间块；不足时显示 WINDOW_INSUFFICIENT。低 γ 与高 γ 结果在对数黏度轴上比较。</p><div id=temperature class=plot></div><div id=profile class=plot></div><div id=viscosity class=plot></div>__FITS__<div id=sensitivity class=plot></div><details><summary>时间步与驱动力敏感性判定</summary><pre>__SENSITIVITY__</pre></details><div id=eos class=plot></div><details open><summary>局部 EOS 与旧目标压力覆盖</summary><pre>__EOS__</pre></details><p>周期固定体积密度 N m / V 是构造值。平衡 EOS 不等于压力出口已验证；没有开放端口，也没有验证真实管壁无滑移。</p></section>
<section><h2>4. 实际血管与冻结测量截面</h2><p>继承完整第一阶段三角表面。橙色线是旧合同中的真实有限测量轮廓，已逆变换回解剖坐标。已有延伸段属于当前表面；没有新增缓冲区几何。</p><div class=controls><label><input id=wall-visible type=checkbox checked> 显示管壁</label><label>透明度 <input id=wall-opacity type=range min=0 max=1 step=.05 value=.35></label><label>突出端口 <select id=port-select>__OPTIONS__</select></label><button id=reset-view>恢复视角</button></div><p id=render-status>等待本地 Plotly 加载</p><div id=vessel-view></div>__PLANS__<details><summary>统一接口缺口与壁面方案</summary><pre>__GAPS__</pre></details></section>
<section><h2>5. 运行资源、证据与阶段状态</h2><p>累计 GPU 子任务墙钟 __USED__ / 3600 秒；同一任务的两个 MPI 进程只计一次。下列显存是设备整体的采样峰值，包含其他程序，不是精确任务峰值。</p>__RESOURCE__<details><summary>候选排除原因</summary><pre>__CANDIDATES__</pre></details><details><summary>源码、编译库与保存补丁身份</summary><pre>__ENVIRONMENT__</pre></details><details><summary>③ SDF 设计准备条件</summary><pre>__STAGE3__</pre></details><p>正式数据：<code>__PACKAGE__</code></p><p>本次仅完成旧纯流体工况整理、单位映射、预算内 DPD 标定及入口出口实现方案评估。未运行整段真实血管，未生成真实血管 SDF 或完整粒子模型，未加入 RBC、微泡或黏附。本阶段用户人工验收状态为 PENDING。</p></section></main>
<script>__PLOTLY__</script><script id=physics-audit-data type=application/json>__PAYLOAD__</script><script>
const a=JSON.parse(document.getElementById('physics-audit-data').textContent),v=document.getElementById('vessel-view');
window.physicsReady=false;
const jobs=Object.entries(a.charts).map(([id,f])=>Plotly.newPlot(id,f.data,f.layout,{responsive:true,displaylogo:false}));
jobs.push(Plotly.newPlot(v,a.figure.data,a.figure.layout,{responsive:true,displaylogo:false,scrollZoom:false}));
Promise.all(jobs).then(()=>{window.physicsReady=true;document.getElementById('render-status').textContent='已加载完整三角表面与实际标定曲线；人工验收仍为 PENDING';}).catch(e=>{document.getElementById('render-status').textContent='渲染失败：'+e;throw e;});
document.getElementById('wall-visible').addEventListener('change',e=>Plotly.restyle(v,{visible:e.target.checked},a.controls.walls));
document.getElementById('wall-opacity').addEventListener('input',e=>Plotly.restyle(v,{opacity:Number(e.target.value)},a.controls.walls));
document.getElementById('port-select').addEventListener('change',e=>{const id=e.target.value;for(const [p,indices] of Object.entries(a.controls.ports))Plotly.restyle(v,{opacity:!id||p===id?1:.12},indices);document.querySelectorAll('[data-entity]').forEach(x=>x.classList.toggle('selected',x.dataset.entity===id));});
document.getElementById('reset-view').addEventListener('click',()=>Plotly.relayout(v,{'scene.camera':{eye:{x:1.4,y:1.5,z:1.1}}}));
// Keep wheel zoom local to the vessel in a long document containing many plots.
// Capture before Plotly's GL handler, then use the same public relayout API as
// the reset control. This also avoids simultaneous document scrolling.
window.vesselWheelEvents=0;
v.addEventListener('wheel',e=>{e.preventDefault();e.stopImmediatePropagation();window.vesselWheelEvents++;const eye=v.layout.scene.camera.eye,f=e.deltaY>0?1.12:1/1.12,n=Math.hypot(eye.x,eye.y,eye.z),scale=Math.max(.2,Math.min(8,n*f))/n;Plotly.relayout(v,{'scene.camera.eye':{x:eye.x*scale,y:eye.y*scale,z:eye.z*scale}});},{passive:false,capture:true});
</script></body></html>"""
    replacements={'STATUS':summary['status'],'BROWSER':browser_status,'STATES':statuses,'PHYSICS':physics,'UNITS':units,'FITS':fits,'RESOURCE':resource,'OPTIONS':options,'PLANS':plans,'USED':f"{summary['budget']['charged_wall_s']:.3f}",'PACKAGE':escape(str(package)),
                  'ESTIMATES':escape(json.dumps(unit['resource_estimates'],ensure_ascii=False,indent=2)),
                  'EOS':escape(json.dumps(summary['eos'],ensure_ascii=False,indent=2)),
                  'SENSITIVITY':escape(json.dumps(summary['sensitivity'],ensure_ascii=False,indent=2)),
                  'CANDIDATES':escape(json.dumps(summary['candidate_assessments'],ensure_ascii=False,indent=2)),
                  'ENVIRONMENT':escape(json.dumps(prov['environment'],ensure_ascii=False,indent=2)),
                  'GAPS':escape(json.dumps({k:v for k,v in boundary.items() if k!='ports'},ensure_ascii=False,indent=2)),
                  'STAGE3':escape(json.dumps(summary['stage3_readiness'],ensure_ascii=False,indent=2)),
                  'PLOTLY':get_plotlyjs(),'PAYLOAD':jsdata(payload)}
    for key,value in replacements.items():text=text.replace('__'+key+'__',value)
    with Path(destination).open('x',encoding='utf-8') as f:f.write(text)
    if '<script src=' in text or 'NaN' in jsdata(payload):raise ValueError('HTML must be local and finite')
    return payload


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',default='py_scripts/fluid_physics.yaml');p.add_argument('--open',action='store_true');p.add_argument('--package',help='explicit immutable result package')
    a=p.parse_args();c=config_read(a.config);campaign=output(Path(c['runs_root'])/c['campaign_id'])
    if a.package:package=Path(a.package).resolve()
    else:
        if not (campaign/'prepared_package.json').exists():p.error('Run CPU prepare first; viewer does not initiate preparation or GPU execution')
        package=export_results(c,Path(read_json(campaign/'prepared_package.json')['directory']))
    directory=output(Path(c['review_root'])/package.name);directory.mkdir(parents=True,exist_ok=True);html=directory/'physics_review.html'
    if not html.exists():
        write_html(package,html);write_json(directory/'html_record.json',{'html':str(html),'sha256':sha256_file(html),'browser_status':'NOT_TESTED','package':str(package),'user_review':'PENDING'})
    else:
        record=read_json(directory/'html_record.json')
        if record['sha256']!=sha256_file(html):raise ValueError('EXISTING_HTML_HASH_MISMATCH')
    print('DATA',package);print('HTML',html);print('GPU_NOT_STARTED')
    if a.open:print(open_windows_file(html))


if __name__=='__main__':main()
