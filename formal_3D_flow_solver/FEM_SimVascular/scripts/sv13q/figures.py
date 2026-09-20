"""Decision figures from actual frozen observations; failed timings are never ranked."""
import json,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'));os.environ['MPLCONFIGDIR']=str(ROOT/'outputs/sv1_3q/plot_cache')
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from sv_validation.provenance import sha256
from sv_validation.sv13q import attempts
R=ROOT/'reports/sv1_3q';font=FontProperties(fname='/mnt/c/Windows/Fonts/msyh.ttc');figures=[]
def read(n):return json.loads((R/(n+'.json')).read_text())
def save(fig,name,title,note):
 fig.suptitle(title,fontproperties=font,fontsize=20,y=.98);fig.text(.5,.025,note,ha='center',fontproperties=font,fontsize=10,color='#334155')
 fig.savefig(R/name,dpi=150,facecolor='white');plt.close(fig);figures.append(dict(path='reports/sv1_3q/'+name,sha256=sha256(R/name),title=title))
def frame(ylabel):
 fig,ax=plt.subplots(figsize=(14,8));fig.subplots_adjust(top=.84,bottom=.18,left=.12,right=.96);ax.set_ylabel(ylabel,fontproperties=font);ax.grid(axis='y',alpha=.2);return fig,ax
base=read('baseline_r1');winner=read('winner');keys=['R1','R2','R3','R5','RA'];good={'R1':base['late']};allcases={}
for c in keys[1:]:
 p=R/(c+'_WINDOW_acceptance.json')
 if p.exists():
  d=read(c+'_WINDOW_acceptance');allcases[c]=d
  if d['status']=='PASS':good[c]=d
for metric,name,title,ylabel in [('wall','ilu_rebuild_candidates.png','ILU 多久重建一次最快？','10 步实际总耗时（秒）'),('builds','ilu_rebuild_count.png','减少了多少次 ILU 重建？','实际 ILU 数值分解次数'),('setup','ilu_setup_time.png','不同策略在准备 ILU 上花了多少时间？','准备 ILU 的时间（秒，含子事件）')]:
 fig,ax=frame(ylabel);vals=[]
 for i,k in enumerate(keys):
  d=good.get(k)
  if d:
   v=d['wall_time_s'] if metric=='wall' else d['profile']['events']['MatLUFactorNum']['count'] if metric=='builds' else d['profile']['events']['PCSetUpOnBlocks']['time_s'];vals.append(v)
   ax.bar(i,v,color='#137c66' if k==winner.get('candidate') else '#6685ad');ax.text(i,v+max(v*.025,.2),f'{v:.2f}' if metric!='builds' else str(v),ha='center')
  else:ax.text(i,2,'FAIL\n无健康完整窗口',ha='center',color='#b42334',fontproperties=font)
 ax.set_xticks(range(5),keys);ax.set_xlim(-.6,4.6);ax.set_ylim(0,max(vals+[1])*1.23)
 save(fig,name,title,'同一 step60 检查点的 10 步窗口。只排名健康完整结果；事件时间有嵌套，不相加。')
fig,axes=plt.subplots(1,2,figsize=(14,8));fig.subplots_adjust(top=.84,bottom=.18,wspace=.30)
for ax,metric in zip(axes,['iterations','delta']):
 names=list(good)
 vals=[(good[k].get('total_attempt_iterations',good[k]['statistics']['total_iterations']) if metric=='iterations' else base['late']['wall_time_s']-good[k]['wall_time_s']) for k in names]
 ax.bar(names,vals,color=['#137c66' if k==winner.get('candidate') else '#6685ad' for k in names]);ax.grid(axis='y',alpha=.2)
 for i,v in enumerate(vals):ax.text(i,v+(max(vals) or 1)*.02,f'{v:.1f}' if metric=='delta' else str(v),ha='center')
 ax.margins(y=.22)
axes[0].set_ylabel('GMRES 迭代总数（所有尝试）',fontproperties=font);axes[1].set_ylabel('相对 R1 节省的时间（秒）',fontproperties=font)
save(fig,'ilu_iterations_tradeoff.png','少重建 ILU 后，求解次数增加了多少？','左图包含所有尝试的 GMRES 迭代；右图显示实际净节省时间，不能只凭重建次数判断。')
fig,ax=frame('观察到的 10 步耗时（秒）');top=read('late_ranking')['top_two'];names=['R1']+top
for i,k in enumerate(names):
 late=base['late'] if k=='R1' else read(k+'_WINDOW_acceptance');early=base['early'] if k=='R1' else read(k+'_EARLY_acceptance')
 for offset,d,label,color in [(-.18,early,'早期：step 10→20','#f59e0b'),(.18,late,'晚期：step 60→70','#6685ad')]:
  if k=='R1' or d['status']=='PASS':
   v=d['wall_time_s'];ax.bar(i+offset,v,width=.34,color=color,label=label if i==0 else None);ax.text(i+offset,v+4,f'{v:.1f}',ha='center')
  else:ax.text(i+offset,10,'FAIL',ha='center',color='#b42334',fontproperties=font)
ax.set_xticks(range(len(names)),names);ax.margins(y=.2);ax.legend(prop=font)
save(fig,'early_vs_late_window.png','同一种策略在早期和接近稳态时表现一样吗？','R1 早期为历史连续运行区间；新策略为独立重启窗口，包含启动开销。只测晚期排名前两名。')
ra_case='REAL_VASCULAR_GPU_ILU_REUSE_WINNER' if winner.get('candidate')=='RA' and (R/'REAL_VASCULAR_GPU_ILU_REUSE_WINNER_acceptance.json').exists() else 'RA_WINDOW' if (R/'RA_WINDOW_acceptance.json').exists() else 'RA_SMOKE'
a=read(ra_case+'_acceptance');rows=attempts((ROOT/'logs/sv1_3q/remote'/(ra_case+'.log')).read_text());rows=[r for r in rows if r['outcome']]
fig,ax=frame('GMRES 迭代次数');x=[];per_step={}
for r in rows:
 j=per_step.get(r['step'],0);x.append(r['step']+.06*j);per_step[r['step']]=j+1
ax.plot(x,[r['outcome']['iterations'] for r in rows],'o-',label='实际迭代数',color='#6685ad');ax.step(x,[1.5*r['outcome']['ref'] for r in rows],where='post',ls='--',color='#b42334',label='1.5 × 参考迭代数')
seen_reasons=set()
for xx,r in zip(x,rows):
 if not r['reuse']:
  reason=r['rebuild_reason'];label={'FIRST_SOLVE':'首次求解','MAX_AGE':'年龄到限','ITERATION_GROWTH':'迭代增长','STALE_FAILURE':'失败后重建'}[reason]
  ax.scatter([xx],[r['outcome']['iterations']],marker='D',s=75,color={'FIRST_SOLVE':'#137c66','MAX_AGE':'#137c66','ITERATION_GROWTH':'#e07b00','STALE_FAILURE':'#b42334'}[reason],zorder=4,label=label if reason not in seen_reasons else None)
  if a['mode']!='full' or reason in ('ITERATION_GROWTH','STALE_FAILURE'):
   ax.annotate(label,(xx,r['outcome']['iterations']),xytext=(3,16),textcoords='offset points',fontsize=10,rotation=18,fontproperties=font)
  seen_reasons.add(reason)
ax.set_xlabel('时间步（小幅错开表示同一步内的不同求解）',fontproperties=font);ax.margins(y=.3);ax.legend(prop=font)
save(fig,'adaptive_rebuild_timeline.png','ILU 是什么时候重新构建的，为什么？',('完整运行；' if a['mode']=='full' else '晚期窗口；')+'菱形表示实际重建及原因；参考值为新 ILU 第一次健康求解，年龄达到 5 时强制重建。')
if winner['full_required'] and (R/'winner_steady_candidate.json').exists():
 h=read('steady_history');c=read('winner_steady_candidate');p=json.loads((ROOT/'configs/sv1_3/policy.json').read_text())
 fig,axes=plt.subplots(1,3,figsize=(16,8));fig.subplots_adjust(top=.84,bottom=.18,wspace=.33)
 for ax,key,limit in zip(axes,['E_u','E_Q','epsilon_mass'],['velocity_change_limit','flow_change_limit','mass_limit']):
  data=h['states'] if key=='epsilon_mass' else h['intervals'];ax.plot([r['step'] for r in data],[r[key] for r in data],'o-');ax.set_yscale('symlog',linthresh=1e-15);ax.axhline(p[limit],ls='--',color='#b42334');ax.axvline(c['first_full_steady_step'],ls=':',color='#137c66');ax.set_xlabel('正式监测时间步',fontproperties=font);ax.set_ylabel({'E_u':'速度变化指标','E_Q':'流量变化指标','epsilon_mass':'质量误差'}[key],fontproperties=font);ax.grid(alpha=.2)
 save(fig,'ilu_winner_steady.png','最快的 ILU 复用策略能否正常跑到稳态？',f"原阈值、连续 5 区间；首次通过 step{c['first_full_steady_step']}，实际停止 step{c['stop_step']}，只运行一次。")
 fig,ax=frame('从初始状态到安全停止的实际耗时（秒）');b=base['full'];vals=[b['wall_time_s'],c['wall_time_s']];ax.bar(['Stage P R1','Stage Q '+winner['candidate']],vals,color=['#6685ad','#137c66'])
 for i,v in enumerate(vals):ax.text(i,v*1.02,f'{v:.2f} s\n'+str(b['stop_step'] if i==0 else c['stop_step'])+' steps',ha='center')
 ax.set_ylim(0,max(vals)*1.2)
 save(fig,'ilu_final_speedup.png','新的 ILU 重建策略实际节省了多少时间？',f"OBSERVATIONAL DEVELOPMENT SPEEDUP = {b['wall_time_s']/c['wall_time_s']:.3f}×；各一次实测，包含实际停止步数差异。")
(R/'visuals.json').write_text(json.dumps(dict(status='PASS' if not winner['full_required'] or len(figures)==8 else 'WINDOW_FIGURES_READY_FULL_PENDING',files=[Path(f['path']).name for f in figures],figures=figures),indent=2,ensure_ascii=False)+'\n');print('Figures:',len(figures))
