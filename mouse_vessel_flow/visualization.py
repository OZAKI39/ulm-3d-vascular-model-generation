"""Present existing scientific results; chart rendering never changes trajectories or flow."""
from pathlib import Path
import html,json,os
import numpy as np
from .paths import PROJECT_ROOT,new_output,project_path,case_config,CASE_NAMES
from .particles.analysis import summarize

def report_figures(output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=new_output(output);names=CASE_NAMES;rows=[summarize(n) for n in names]
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    colors={'O1':'#2c7fb8','O2':'#fdae61','O3':'#66bd63','SUPPORTED_STATIONARY':'#777777'}
    bottom=np.zeros(len(names))
    for label in colors:
        value=np.array([r['outlets'][label] if label.startswith('O') else r['statuses'][label] for r in rows])
        axes[0].bar(names,value,bottom=bottom,label=label,color=colors[label]);bottom+=value
    axes[0].set(ylabel='Number of original trajectories',title='Frozen campaign outcomes')
    axes[0].legend(fontsize=8)
    for i,r in enumerate(rows):axes[0].text(i,bottom[i]+15,str(r['count']),ha='center')
    for i,n in enumerate(names):
        cfg=case_config(n);events=json.loads(project_path(cfg['particles']['cohort']).read_text())['events']
        d=np.sort([e['diameter_um'] for e in events]);axes[1].plot(d,np.arange(1,len(d)+1)/len(d),label=n)
    axes[1].set(xlabel='Accepted bubble diameter (µm)',ylabel='Empirical CDF',title='Original accepted cohorts (no resampling)')
    axes[1].legend(fontsize=9);axes[1].grid(alpha=.25)
    fig.savefig(out/'frozen_campaign_comparison.png',dpi=300);fig.savefig(out/'frozen_campaign_comparison.pdf');plt.close(fig)
    links=[]
    for n in names:
        cfg=case_config(n);report=project_path(cfg['particles']['report'])
        for f in sorted(report.rglob('*.mp4')):
            links.append(f'<li>{html.escape(n)}: <a href="{html.escape(os.path.relpath(f,out))}">{html.escape(f.name)}</a></li>')
    page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>小鼠血管冻结结果</title>
<style>body{font:17px sans-serif;max-width:1200px;margin:2em auto;padding:0 1em}img{width:100%}li{margin:.6em}</style>
<h1>小鼠血管流场与微泡：冻结结果核验</h1><p>图表只读取原始事件与指标。接触支持静止不等于出流；动画按轨迹年龄对齐，不代表同时注入。</p>
<img src="frozen_campaign_comparison.png" alt="Frozen campaign comparison"><p>CPU 执行粒子动力学；原 GPU 证据涉及诊断与渲染。图表不构成 CFD 局部精度或时间步收敛证明。</p>
<h2>原始动画</h2><ul>'''+''.join(links)+'</ul></html>'
    (out/'index.html').write_text(page)
    (out/'source_summary.json').write_text(json.dumps(rows,indent=2)+'\n')
    return dict(output=str(out),figure_png='frozen_campaign_comparison.png',figure_pdf='frozen_campaign_comparison.pdf',html='index.html',new_simulation=False)
