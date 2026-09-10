"""Standalone Chinese review from stored evidence, with no solver side effects."""
import html
import json
import time
from pathlib import Path
from plotly.offline import get_plotlyjs
from .analysis import collect
from .physics import paths
from py_scripts.fluid_physics.common import atomic_state,sha256_file,read_json,fingerprint


def fmt(x,n=3):
    return '未测' if x is None else f'{float(x):.{n}f}'


def table_row(values,status=None):
    return '<tr>'+''.join('<td'+(' class="'+html.escape(status)+'"' if status and i==1 else '')+'>'+html.escape(str(x))+'</td>' for i,x in enumerate(values))+'</tr>'


def build_html(data):
    here=Path(__file__).parent;lookup={r['task']:r for r in data['runs']}
    summary=[]
    for key,label in [('hemocell','HemoCell · CPU/MPI'),('mirheo','Mirheo DPD · GPU')]:
        r=data['repeat_statistics'][key]
        summary.append(table_row([label,r['count'],fmt(r['mean_s']),fmt(r['minimum_s'])+' – '+fmt(r['maximum_s']),fmt(r['sample_std_s']),fmt(r['mean_s']/4 if r['mean_s'] is not None else None)]))
    runrows=[]
    for r in data['runs']:
        cp=r['completion'];endpoint=fmt(cp.get('strain_end'),4) if cp.get('completed') else (fmt(r.get('observed_last_saved_strain'),4)+'（最后存帧，未完成）' if r.get('observed_last_saved_strain') is not None else '未完成 / 不适用')
        steps=f"{cp.get('actual_steps','未知')} / {cp.get('prep_steps','未知')}"
        runrows.append(table_row([r['task'],r['execution']['status'],fmt(r['execution']['elapsed_monotonic_s']),fmt(r.get('analysis_s')),cp.get('cell_count','1（已存膜）' if r['frames'] else '—'),endpoint,steps],r['execution']['status']))
    m=data.get('material') or {};mp=data.get('membrane_mapping') or {};dg=data.get('dimensionless') or {}
    cg=dg.get('HemoCell_effective_proxy',{});mg=dg.get('Mirheo_effective_proxy',{});g=data['geometry'] or {}
    model=[
        ['参考网格、中心、取向','共同导出的 HemoCell 原生初态','同一个 OFF，单位四元数',f"输入一致；A0={fmt(g.get('area'),6)}，V0={fmt(g.get('volume'),6)}，减缩体积={fmt(g.get('reduced_volume'),6)}。准备后形状另检。"],
        ['液体密度、黏度、λ','ρ*=8；ν*='+fmt(m.get('primary',{}).get('nu'),6),'n=8、m=1；独立 DPD 实测 ν*='+fmt(m.get('primary',{}).get('nu'),6),'λ=1：内内、外外和内外采用同一液体参数。黏度半步相差 '+fmt(100*m.get('strict_difference',0),3)+'%，5% 门槛未被放宽。'],
        ['原生小变形有效响应',fmt(mp.get('HemoCell',{}).get('Gs_affine_effective'),6),fmt(mp.get('Mirheo',{}).get('Gs_affine_effective'),6),'总增量力功定义；剪切约差 6%，伸长最大约差 18%。不是完整非线性本构等价。'],
        ['Re；Ca（有效响应代理）',fmt(cg.get('Re'),6)+'；'+fmt(cg.get('Ca'),6),fmt(mg.get('Re'),6)+'；'+fmt(mg.get('Ca'),6),'公共 Re；Ca 依赖各自测得的仿射有效响应，不把它称为已验证的精确连续膜 Gs。'],
        ['弯曲与 B','HO 局部参考曲率；B 未核实','Kantor kb=8，θ0=6.97°；B 代理='+fmt(mg.get('B'),6),'两种弯曲模型不同；同名数值 8 不构成等价证明。'],
        ['面积、体积约束','kArea=5、kVolume=0.2','ka=279.69313、kv_tot=1188.746；ka_tot=0','局部面积方向不同，体积只作源码一阶映射。非线性范围独立响应未通过。'],
        ['膜惯性与表面耗散','IBM 速度插值；eta_m=0','节点质量 1，共 642；gammaC=0','Mirheo 膜质量/排开液体质量≈0.9893，不能假定可忽略。没有独立证明释放耗散等效。'],
        ['热涨落','确定性 LBM，无分子热噪声','流体和膜耦合 DPD 保留 kBT=1','整体 kBT/(Gs_eff a²)≈'+fmt(mg.get('thermal_to_shear_energy'),5)+'，局部噪声影响和角度偏移仍进入可比性限制。'],
        ['边界及双向作用','X/Y 周期，Z 移动速度边界；实际 hc.iterate','MovingPlane、同速壁粒子；两侧 DPD + bounce_back','流体只生成一次再分类；实际网格、流体变化和探针均已保存。跨膜迹象不能由重新分类掩盖。'],
        ['空间分辨率与精度','48×48×49；dx*=0.5；double','110,592 个液体粒子；rc=1；single','CPU/MPI 与单 GPU 的实际选型比较。膜都是 642 顶点，流体分辨率不等价；空间收敛未验证。'],
        ['步长与逻辑更新','dt*=0.037037；半步 0.018519','dt*=0.001；半步 0.0005','每个流体步均更新膜及耦合，未用降低调用频率换速度。主 HemoCell 5,535 次 iterate（准备 135）；主 Mirheo 205,000 演化步另含 2,000 壁准备步；失败实际步数仅知保存下界。'],
        ['保护与范围','链接原有 HemoCell 库，新增独立案例','既有 Mirheo 库和 C++/CUDA 保留','固定提交、库哈希及每次新 worker/案例源码快照可核查；未提交或推送。']]
    labels={'empty_hemocell':'HemoCell 空通道','empty_dpd':'DPD 空通道','strict_hemocell':'HemoCell 半步窗口','strict_mirheo':'Mirheo 半步窗口','viscosity_half_dt':'DPD 黏度半步','native_membrane_response':'原生膜响应','prepared_shapes':'零剪切准备后的形状'}
    screens=[]
    for key,r in data['screen_checks'].items():
        if key.startswith('empty_') and r.get('metrics'):
            q=r['metrics'];note='速度 L2 '+fmt(q['empty_velocity_relative_l2']*100)+'%；剪切率 '+fmt(q['effective_shear_relative_error']*100)+'%；近壁密度 '+fmt(q['near_wall_density_relative_error']*100)+'%；质量漂移 '+fmt(q['mass_relative_drift']*100,6)+'%'
            if 'temperature_relative_error' in q:note+='；平均温度偏差 '+fmt(q['temperature_relative_error']*100)+'%'
        elif r.get('status')=='INCOMPLETE_WINDOW':note='计划窗口 Γ=0.5～2 不完整；原配置第 1 次最后存帧 Γ=0.6，半步 Γ=0.7。不能用不完整窗口的均值差作收敛结论。'
        elif 'D_absolute_difference' in r:note='Γ=0.5～2，平均 D 差 '+fmt(r['D_absolute_difference'],6)+' / 阈值 '+fmt(r['frozen_absolute_tolerance'],6)+'；平均倾角差 '+fmt(r['theta_mean_difference_deg'])+'°。非空间收敛证明。'
        elif 'relative_difference' in r:note='差 '+fmt(r['relative_difference']*100)+'%，门槛 '+fmt(r.get('tolerance',0)*100)+'%'
        elif key.endswith('_geometry'):
            q=r.get('details',{});note='面积最大漂移 '+fmt(q.get('area_drift',0)*100)+'%；体积 '+fmt(q.get('volume_drift',0)*100)+'%；存帧非相邻自交 '+str(q.get('nonadjacent_intersections','未测'))+'；穿墙 '+str(not q.get('no_wall_crossing',False))
        elif key.endswith('_fluid_membership'):
            q=r['details'];note=str(q['mismatches'])+' 次探针位置/成员不符 / '+str(q['tested_point_frames'])+' 点帧；不是严格渗透率测定。'
        elif key=='prepared_shapes':note='；'.join('重复 '+str(q['repeat'])+'：形状 RMS/a='+fmt(q['relative_shape_rms']*100)+'%，门槛 2%' for q in r.get('pairs',[])) or '双方对应的准备后形状尚不完整'
        else:note='未测或窗口不完整；详见嵌入 JSON'
        screens.append(table_row([labels.get(key,key),r['status'],note],r['status']))
    encoded=json.dumps(data,ensure_ascii=False,allow_nan=False).replace('</','<\\/')
    page=(here/'review_template.html').read_text()
    for key,value in {'SUMMARY_ROWS':''.join(summary),'RUN_ROWS':''.join(runrows),'MODEL_ROWS':''.join(table_row(r) for r in model),'SCREEN_ROWS':''.join(screens),'PLOTLY':get_plotlyjs(),'DATA':encoded,'REVIEW_JS':(here/'review.js').read_text()}.items():page=page.replace('@@'+key+'@@',value)
    return page


def export(c):
    start=time.perf_counter();base,runs,out=paths(c);out.mkdir(parents=True,exist_ok=True)
    inputs={str(p):sha256_file(p) for p in runs.rglob('output_sha256.json')}
    inputs.update({str(p):sha256_file(p) for p in Path(__file__).parent.iterdir() if p.suffix in ('.py','.html','.js')})
    evidence_names=('geometry.json','fluid_measured.json','membrane_measured.json','formal_frozen.json','dimensionless_audit.json','dimensionless_measured.json','material_diagnostics.json','runtime_candidate_frozen.json','cuda_queue_diagnostic_result.json','native_failure_diagnostic_plan.json','native_failure_final.json','screen_frozen.json','authorization.json','environment_authorized_start.json','protection_current.json','cpu_tests_final.json')
    inputs.update({str(base/name):sha256_file(base/name) for name in evidence_names if (base/name).exists()})
    if (base/'environment_authorized_start_corrected.json').exists():inputs['corrected_environment']=sha256_file(base/'environment_authorized_start_corrected.json')
    inputs['config']=c['_config_sha256'];key=fingerprint(inputs);latest=base/'review_latest.json'
    if latest.exists():
        r=read_json(latest)
        if r.get('cache_key')==key and Path(r['html']).is_file() and sha256_file(r['html'])==r['html_sha256']:
            return {**r,'cached':True}
    analysis_start=time.perf_counter();data=collect(c);collection_s=time.perf_counter()-analysis_start
    atomic_state(out/'results.json',data);render_start=time.perf_counter();page=out/'single_rbc_review.html';page.write_text(build_html(data),encoding='utf8')
    record=dict(html=str(page),results=str(out/'results.json'),html_sha256=sha256_file(page),html_render_s=time.perf_counter()-render_start,collection_wall_s=collection_s,total_review_build_s=time.perf_counter()-start,human_review='PENDING',browser_check='NOT_RUN',cache_key=key)
    atomic_state(base/'review_latest.json',record);return record
