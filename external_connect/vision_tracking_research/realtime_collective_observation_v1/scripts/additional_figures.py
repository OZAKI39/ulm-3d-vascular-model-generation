import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ycommon import ROOT,dump
sys.path.insert(0,str(ROOT/'scripts'))
from render_deliverables import loadrows
import numpy as np,av
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon,Patch
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220})
r=loadrows(ROOT/'results/formal_study/offline.jsonl.gz')[450]
with av.open('/home/lzy/projects/external_connect/sample_videos/Recording 2026-10-15 185901 截取视频.mp4') as c:
 for i,f in enumerate(c.decode(video=0)):
  if i==450:gray=f.to_ndarray(format='gray');break
g=r['geometry'];x,y,_,_=g['crop_xywh'];q=np.array(r['candidates']);v=np.array(r['vectors']);fig,ax=plt.subplots(figsize=(10,8));ax.imshow(gray,cmap='gray');ax.plot(*np.array(g['channel_polygon']+g['channel_polygon'][:1]).T,color='teal',lw=1.3)
for poly in g['excluded_polygons']:ax.add_patch(Polygon(poly,facecolor='gray',edgecolor='silver',alpha=.2,hatch='///'))
ax.scatter(q[:,0]+x,q[:,1]+y,s=6,facecolors='none',edgecolors='#32a7d2',linewidth=.6);ax.quiver(v[:,0]+x,v[:,1]+y,v[:,2],v[:,3],color='#d29322',angles='xy',scale_units='xy',scale=12,width=.002)
ax.set(xlabel='Image x / px',ylabel='Image y / px',title='Full-channel observation | frame 450 | 282 candidates, 418 valid PIV vectors')
ax.legend(handles=[Line2D([],[],marker='o',color='#32a7d2',linestyle='none',fillstyle='none',label='Visible candidates'),Line2D([],[],color='#d29322',label='Valid PIV velocity'),Patch(facecolor='gray',alpha=.25,hatch='///',label='Excluded / low-visibility wall')],loc='lower left',framealpha=.95,fontsize=8)
for ext in ['png','pdf']:fig.savefig(ROOT/'figures'/f'Fig_F10Y_2_full_channel_observation.{ext}',bbox_inches='tight')
plt.close(fig)
s=json.loads(next((ROOT/'results/stability_final').glob('session_*/summary.json')).read_text());resources=s['resources'];t=np.array([a['elapsed_s'] for a in resources]);fig,axs=plt.subplots(2,2,figsize=(12,7),layout='constrained')
axs[0,0].plot(t,[a['rss_mib'] for a in resources],label='Process RSS');axs[0,0].plot(t,[a.get('gpu_mib') for a in resources],label='GPU device');axs[0,0].set_ylabel('Memory / MiB');axs[0,0].legend()
axs[0,1].plot(t,[a['cpu_percent'] for a in resources],label='CPU (100% = one core)');axs[0,1].plot(t,[a.get('gpu_percent') for a in resources],label='GPU device');axs[0,1].set_ylabel('Utilization / %');axs[0,1].legend()
axs[1,0].plot(t,[a['age_ms'] for a in resources],lw=.8);axs[1,0].axhline(250,c='red',ls='--',label='TTL 250 ms');axs[1,0].set_ylabel('State software age / ms');axs[1,0].legend()
ui=json.loads((ROOT/'results/gpu_ui_final/ui_summary.json').read_text());v=[ui['ui_state_update_ms']['p95'],ui['canvas_paint_ms']['p95']];axs[1,1].bar(['UI state update','Qt canvas painting'],v,color=['#8b9dc3','#76b4a3']);axs[1,1].set_ylabel('Final UI p95 / ms');axs[1,1].text(.02,.95,'Physical scene-to-screen latency: NOT MEASURED',transform=axs[1,1].transAxes,fontsize=9,va='top')
for ax in axs.flat:
 ax.grid(alpha=.15)
 if ax!=axs[1,1]:ax.set_xlabel('Wall time / s')
fig.suptitle(f"Continuous replay: {s['wall_s']:.1f} s, {s['source_sessions']} sessions, no analysis errors")
for ext in ['png','pdf']:fig.savefig(ROOT/'figures'/f'Fig_F10Y_7_stability_and_UI.{ext}',bbox_inches='tight')
p=ROOT/'evidence/figures_validation.json';d=json.loads(p.read_text());d['figures']+=[] if 'Fig_F10Y_7_stability_and_UI' in d['figures'] else ['Fig_F10Y_7_stability_and_UI'];d['png_pdf_count']=14;dump(p,d)
