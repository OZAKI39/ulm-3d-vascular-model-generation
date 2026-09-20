"""Decision plots from closed actual artifacts; unavailable times stay unavailable."""
import json,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
os.environ['MPLCONFIGDIR']=str(ROOT/'outputs/sv1_3p/plot_cache')
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from sv_validation.provenance import sha256
R=ROOT/'reports/sv1_3p';font=FontProperties(fname='/mnt/c/Windows/Fonts/msyh.ttc');figures=[]
def read(n):return json.loads((R/(n+'.json')).read_text())
def save(fig,name,title,note):
 fig.suptitle(title,fontproperties=font,fontsize=20,y=.98);fig.text(.5,.025,note,ha='center',fontproperties=font,fontsize=10,color='#334155')
 fig.savefig(R/name,dpi=150,facecolor='white');plt.close(fig);figures.append(dict(path='reports/sv1_3p/'+name,sha256=sha256(R/name),title=title,data_source='candidate_summary.json; original Stage O profile; winner artifacts'))
def frame(ylabel):
 fig,ax=plt.subplots(figsize=(14,8));fig.subplots_adjust(top=.84,bottom=.18,left=.12,right=.96);ax.set_ylabel(ylabel);ax.grid(axis='y',alpha=.2);return fig,ax
s=read('candidate_summary');w=read('winner');rows=s['candidates'];labels=['O Baseline','P1 Reuse','P2 hypre ILU0','P3 BoomerAMG','P4 GAMG','P5 AmgX'];keys=['O',*rows]
good={'O':s['baseline'],**{k:d['window'] for k,d in rows.items() if d.get('window') and d['window']['status']=='PASS'}}
fig,ax=frame('Complete 10-step process wall time (s)')
for i,k in enumerate(keys):
 if k in good:
  v=good[k]['wall_time_s'];ax.bar(i,v,color='#137c66' if k==w['fastest'] else '#6685ad');ax.text(i,v+4,f'{v:.2f} s',ha='center')
 else:ax.text(i,12,'No healthy\n10-step result',ha='center',color='#b42334')
ax.set_xticks(range(6),labels);ax.set_xlim(-.6,5.6);ax.set_ylim(0,max(d['wall_time_s'] for d in good.values())*1.2)
save(fig,'gpu_pc_candidates.png','哪一种 GPU 辅助求解方法最快？','只排名健康完成的 10 步窗口；P 候选包含 GPU 事件计时开销，O 基线直接复用历史结果。')
old=json.loads((ROOT/'reports/sv1_3o/profile_summary.json').read_text());op=old['profiles'][old['runs']['baseline']]['events']
for event,name,title in [('setup','gpu_pc_setup_cost.png','哪种方法准备辅助求解器最省时间？'),('PCApply','gpu_pc_apply_cost.png','真正求解时，哪种辅助方法最省时间？')]:
 fig,ax=frame('Inclusive PETSc event time (s)')
 vals=[]
 for i,k in enumerate(keys):
  e=op if k=='O' else good[k]['profile']['events'] if k in good else {}
  key='PCSetUpOnBlocks' if event=='setup' and k in ('O','P1') else 'PCSetUp' if event=='setup' else event
  v=e.get(key,{}).get('time_s');vals.append(v or 0)
  if v is not None:ax.bar(i,v,color='#137c66' if k==w['fastest'] else '#6685ad');ax.text(i,v+1,f'{v:.2f}\n{key}',ha='center',fontsize=9)
  else:ax.text(i,2,'N/A',ha='center',color='#b42334')
 ax.set_xticks(range(6),labels);ax.set_xlim(-.6,5.6);ax.set_ylim(0,max(vals+[1])*1.24)
 save(fig,name,title,'O 使用已保存的诊断窗口；失败或未完成的候选不填补时间。事件有嵌套，不能相加。')
fig,axes=plt.subplots(1,2,figsize=(14,8));fig.subplots_adjust(top=.84,bottom=.18,wspace=.30)
for ax,key in zip(axes,['iterations','cost']):
 names=list(good);vals=[good[k]['statistics']['total_iterations'] if key=='iterations' else 1000*good[k]['wall_time_s']/good[k]['statistics']['total_iterations'] for k in names]
 ax.bar(names,vals,color=['#137c66' if k==w['fastest'] else '#6685ad' for k in names]);ax.grid(axis='y',alpha=.2);ax.set_ylim(0,max(vals)*1.2)
 for i,v in enumerate(vals):ax.text(i,v*1.02,f'{v:.2f}' if key=='cost' else str(v),ha='center')
axes[0].set_ylabel('Total KSP iterations');axes[1].set_ylabel('Process wall / total iterations (ms)')
save(fig,'gpu_pc_iterations.png','更快的方法是因为迭代更少，还是每次更快？','只比较健康完成窗口；总时间除以迭代数包含组装、启动和输出，并非纯 GPU 单步成本。')
fig,ax=frame('Observed device-wide VRAM (MiB)');vals=[]
for i,k in enumerate(keys):
 d=s['baseline'] if k=='O' else rows[k]['window'] or rows[k]['smoke'];v=d.get('sampled_device_memory_max_MiB') if d else None;vals.append(v or 0)
 if v is not None:ax.bar(i,v,color='#6685ad');ax.text(i,v+50,f'{v:.0f}\n'+('window' if k=='O' else rows[k]['timing_scope']),ha='center')
 else:ax.text(i,80,'Not observed',ha='center',color='#b42334')
ax.set_xticks(range(6),labels);ax.set_xlim(-.6,5.6);ax.set_ylim(0,max(vals+[1])*1.23)
save(fig,'gpu_pc_memory.png','不同方法实际用了多少显存？','设备级采样最大值：smoke 每 2 秒、窗口每 15 秒；不是进程精确峰值，各候选观察时长不同。')
if w['full_required']:
 h=read('steady_history');c=read('winner_steady_candidate');policy=json.loads((ROOT/'configs/sv1_3/policy.json').read_text())
 fig,axes=plt.subplots(1,3,figsize=(16,8));fig.subplots_adjust(top=.84,bottom=.18,wspace=.33)
 for ax,key,limit in zip(axes,['E_u','E_Q','epsilon_mass'],['velocity_change_limit','flow_change_limit','mass_limit']):
  data=h['states'] if key=='epsilon_mass' else h['intervals'];ax.plot([r['step'] for r in data],[r[key] for r in data],'o-');ax.set_yscale('symlog',linthresh=1e-15);ax.axhline(policy[limit],ls='--',color='#b42334');ax.axvline(c['first_full_steady_step'],ls=':',color='#137c66');ax.set_xlabel('Formal monitor timestep');ax.set_ylabel(key);ax.grid(alpha=.2)
 save(fig,'gpu_pc_winner_steady.png','最快的方法能否正常跑到稳态？',f"原稳态门槛、连续 5 区间。首次通过 step{c['first_full_steady_step']}，安全停止 step{c['stop_step']}；仅一次完整运行。")
 fig,ax=frame('t=0 to safe steady stop process wall (s)');base=read('reference_freeze')['GPU_baseline'];values=[base['wall_time_s'],c['wall_time_s']];ax.bar(['Stage O','Stage P '+w['candidate']],values,color=['#6685ad','#137c66'])
 for i,v in enumerate(values):ax.text(i,v*1.02,f'{v:.2f} s\n'+str(base['stop_step'] if i==0 else c['stop_step'])+' steps',ha='center')
 ax.set_ylim(0,max(values)*1.2)
 save(fig,'gpu_pc_final_speedup.png','新的 GPU 辅助求解器实际节省了多少时间？',f"包含实际安全停止步数差异；各一次开发运行，完整用时比 {base['wall_time_s']/c['wall_time_s']:.3f}×；CPU/GPU 科学等价仍推迟。")
(R/'visuals.json').write_text(json.dumps(dict(status='PASS',figures=figures),indent=2,ensure_ascii=False)+'\n');print('Plots:',len(figures))
