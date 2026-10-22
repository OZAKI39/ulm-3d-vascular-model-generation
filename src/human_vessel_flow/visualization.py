"""Reproducible read-only figures. Display never changes CFD or particle arrays."""
from pathlib import Path
import html, os
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .paths import ROOT, read_json, fresh_directory, write_json, sha256
from .particles import historical_summary

def build_review(output,root=ROOT):
    output=fresh_directory(output,root)
    data=historical_summary(root)
    flow=read_json(root/'results/reports/original_cfd/ACTIVE_FLOW.json')
    wall=pv.read(root/'cases/balanced_pressure_final/postprocess/wall_wss_si.vtp')
    values=np.asarray(wall.cell_data['WSS_raw_Pa']);area=np.asarray(wall.cell_data['Area_m2'])
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    figure,axes=plt.subplots(1,3,figsize=(15,4.8),layout='constrained')
    labels=['O1 exit','O2 exit','O3 exit','Censored','Numerical\nfailure']
    counts=[458,527,437,141,2];colors=['#397fbd']*3+['#e2a538','#ba4146']
    bars=axes[0].bar(labels,counts,color=colors);axes[0].bar_label(bars,padding=3)
    axes[0].set(title='Historical population: 1,565 accepted',ylabel='Number of bubbles',ylim=(0,720))
    axes[0].tick_params(axis='x',labelsize=9);axes[0].grid(axis='y',alpha=.2)
    axes[0].text(.02,.96,'1,422 exits; 141 censored; 2 failures\nNot a fixed-window production run',transform=axes[0].transAxes,va='top',fontsize=10)
    fractions=flow['measurements']['outlet_fractions']
    q=[fractions[k]*100 for k in sorted(fractions)]
    bars=axes[1].bar(['O1','O2','O3'],q,color='#397fbd');axes[1].bar_label(bars,fmt='%.4f%%',padding=4)
    axes[1].set(title='Historical pressure-outlet CFD',ylabel='Share of total outlet flow (%)',ylim=(0,50))
    axes[1].axhline(100/3,color='gray',ls='--',lw=1);axes[1].grid(axis='y',alpha=.2)
    axes[1].text(.03,.95,'Inlet: 18 mL/min\nCalibrated pressures; local accuracy not certified',transform=axes[1].transAxes,va='top',fontsize=9)
    axes[2].hist(values,bins=50,weights=area/area.sum()*100,color='#397fbd')
    axes[2].set(title='Original wall-facet WSS',xlabel='WSS (Pa)',ylabel='Wall area (%)',ylim=(0,31))
    axes[2].text(.03,.96,f'{len(values):,} original facets\n{int((values==0).sum())} exactly zero; no smoothing',transform=axes[2].transAxes,va='top',fontsize=10)
    figure.suptitle('BraVa frozen evidence: migration changes file organization only',fontsize=15)
    figure.savefig(output/'frozen_evidence.png',dpi=300);figure.savefig(output/'frozen_evidence.pdf');plt.close(figure)
    # Direct links to the four unchanged media files; no old workspace URLs.
    views=root/'results/visualization/population_four_views'
    videos=sorted(views.rglob('*.mp4'))
    links=[]
    for p in videos:
        rel=os.path.relpath(p,output)
        links.append(f'<section><h3>{html.escape(p.parent.name)}</h3><video controls preload="metadata" src="{html.escape(rel)}"></video></section>')
    contact=views/'CONTACT_SHEET.png'
    contact_html=f'<img src="{os.path.relpath(contact,output)}" alt="Original four-view contact sheet">' if contact.exists() else ''
    page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>BraVa 工程迁移核查</title>
<style>body{font:17px/1.7 sans-serif;max-width:1250px;margin:32px auto;padding:0 18px;color:#202735}img{max-width:100%}video{width:100%;max-width:1000px;background:black}section{margin:24px 0}.note{padding:18px;background:#f0f3f7;border-radius:8px}a{color:#135a9b}</style>
<h1>BraVa 人体血管：迁移后结果核查</h1>
<p>本页读取原有结果。没有重新计算流场，没有新增正式微泡批次，也没有改动速度、壁面数据或终态。</p>
<img src="frozen_evidence.png" alt="Frozen flow, population and raw WSS"><p><a href="frozen_evidence.pdf">下载矢量 PDF</a></p>
<div class="note"><b>怎样理解这些结果：</b>1565 个微泡中，1422 个已从出口离开，141 个在观察结束时仍未出流，2 个发生数值失败。左右两个打印姿态使用同一物理解。
新协议的灌注时间是 1.240815660 秒，随后观察 10 秒，目前没有新协议的正式大批量结果。
全局流量平衡、CPU/GPU 一致或动画流畅，都不能证明每个近壁位置的流场准确。</div>
<h2>原有四视角动画</h2><p>黑底、蓝泡、英文、原始坐标网格和字号保持。逐帧低速筛选仅影响显示；科学账本仍包含全部 1565 个微泡。</p>'''+contact_html+''.join(links)+'''
<h2>仍未通过的科学检查</h2><ul><li>Gmsh V2 参考管局部场精度：FAIL。</li><li>P1 场误差及轨迹敏感性：原诊断结论保留。</li><li>P2/P1：资源限制停止，曲面压力回填仍有失败；没有正式高阶流场。</li><li>V3 网格：仅为提案。</li></ul></html>'''
    (output/'OPEN_RESULTS.html').write_text(page)
    manifest=dict(status='PASS_GENERATED_FROM_FROZEN_DATA',historical_population=data,
        files={p.name:sha256(p) for p in output.iterdir() if p.is_file()},
        unchanged_videos=[dict(path=str(p.relative_to(root)),sha256=sha256(p)) for p in videos],
        inputs={str(p.relative_to(root)):sha256(p) for p in [root/'results/population/realization_0000_contactfix/SUMMARY.json',root/'results/reports/original_cfd/ACTIVE_FLOW.json',root/'cases/balanced_pressure_final/postprocess/wall_wss_si.vtp']})
    write_json(output/'FIGURE_IDENTITY.json',manifest)
    return dict(status='PASS',page=str(output/'OPEN_RESULTS.html'),videos=len(videos))
