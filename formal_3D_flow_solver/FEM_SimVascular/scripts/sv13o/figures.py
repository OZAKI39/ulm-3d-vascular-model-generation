"""Five development-performance plots from closed actual run artifacts."""
import json,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
os.environ['MPLCONFIGDIR']=str(ROOT/'outputs/sv1_3o/plot_cache')
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
from sv_validation.provenance import sha256
R=ROOT/'reports/sv1_3o';font=FontProperties(fname='/mnt/c/Windows/Fonts/msyh.ttc');figures=[]
def read(n):return json.loads((R/(n+'.json')).read_text())
def save(fig,name,title,note,source):
 fig.suptitle(title,fontproperties=font,fontsize=20,y=.98);fig.text(.5,.028,note,ha='center',fontproperties=font,fontsize=10,color='#334155')
 fig.savefig(R/name,dpi=150,facecolor='white');plt.close(fig);figures.append(dict(path='reports/sv1_3o/'+name,sha256=sha256(R/name),title=title,data_source=source))
w=read('winner');runs={k:read('PERF_'+k+'_acceptance') for k in 'ABCD'}
labels={'A':'A: ASM2 + ILU2 / R100','B':'B: direct ILU2 / R100','C':'C: direct ILU1 / R100','D':'D: winner PC / R200'}
fig,ax=plt.subplots(figsize=(14,8));fig.subplots_adjust(top=.86,bottom=.18,left=.12,right=.95)
for i,k in enumerate('ABCD'):
 d=runs[k]
 if d['status']=='PASS':
  t=d['wall_time_s'];ax.bar(i,t,color='#137c66' if k==w['candidate'] else '#6685ad',width=.6);ax.text(i,t+5,f'{t:.2f} s',ha='center')
 else:ax.text(i,10,'FAIL — incomplete' if d['status']=='FAIL' else 'NOT RUN',ha='center',color='#b42334')
ax.set_xticks(range(4),[labels[k] for k in 'ABCD'],fontsize=10);ax.set_xlim(-.6,3.6);ax.set_ylabel('Fixed-window process wall time (s)');ax.set_ylim(0,max(d.get('wall_time_s',0) for d in runs.values())*1.2);ax.grid(axis='y',alpha=.2)
save(fig,'gpu_tuning_candidates.png','哪些 GPU 求解设置真正节省了时间？',f"同一 step60→70 窗口；失败窗口不参与计时排名。Winner：{w['candidate']}；相对 A 减时 {w['fixed_window_reduction']*100:.2f}%。",'PERF_A/B/C/D_acceptance.json; winner.json')
good=[k for k,d in runs.items() if d['status']=='PASS'];fig,axes=plt.subplots(1,2,figsize=(14,8));fig.subplots_adjust(top=.85,bottom=.19,wspace=.30)
for a,key in zip(axes,['iterations','cost']):
 values=[runs[k]['statistics']['total_iterations'] if key=='iterations' else 1e3*runs[k]['wall_time_s']/runs[k]['statistics']['total_iterations'] for k in good]
 a.bar(good,values,color=['#137c66' if k==w['candidate'] else '#6685ad' for k in good]);a.set_xlabel('Healthy completed candidate');a.grid(axis='y',alpha=.2)
 for i,v in enumerate(values):a.text(i,v*1.015,f'{v:.2f}' if key=='cost' else str(v),ha='center')
 a.set_ylim(0,max(values)*1.18)
axes[0].set_ylabel('Total KSP iterations');axes[1].set_ylabel('Process wall / total KSP iterations (ms)')
save(fig,'ksp_cost_comparison.png','时间减少来自更少迭代还是每次迭代更快？','右图是总 wall / 迭代数，包含启动、组装和 I/O；不能视为纯 GPU 单次迭代耗时。','Completed candidate acceptance artifacts')
profiles=read('profile_summary');names=[profiles['runs']['baseline'],profiles['runs']['winner']]
events=['PCSetUpOnBlocks','MatLUFactorNum','KSPSolve','PCApply','MatILUFactorSym','MatMult','PCSetUp','MatAssemblyEnd'];fig,axes=plt.subplots(1,2,figsize=(16,8));fig.subplots_adjust(top=.84,bottom=.21,left=.13,right=.97,wspace=.48)
for ax,name,label in zip(axes,names,['Baseline A','Selected winner']):
 p=profiles['profiles'][name];values=[p['events'].get(e,{}).get('time_s',0) for e in events]
 ax.barh(events[::-1],values[::-1],color='#6685ad' if label=='Baseline A' else '#137c66');ax.set_xlabel('Inclusive PETSc event time (s)');ax.set_title(label+' — '+name+(' (same profile)' if names[0]==names[1] and label=='Selected winner' else ''));ax.grid(axis='x',alpha=.2)
 for i,v in enumerate(values[::-1]):ax.text(v,i,' '+(f'{v:.2g}' if 0<v<.01 else f'{v:.2f}'),va='center',fontsize=9)
 ax.set_xlim(0,max(values)*1.22)
save(fig,'gpu_time_breakdown.png','GPU 血流计算主要把时间花在哪里？','仅 baseline / winner 的轻量 profile；事件存在嵌套，不能相加。FEM 组装与 VTU/restart I/O 未单独计时。','profile_summary.json; native PETSc log_view')
h=read('steady_history');c=read('optimized_steady_candidate');policy=json.loads((ROOT/'configs/sv1_3/policy.json').read_text())
fig,axes=plt.subplots(1,3,figsize=(16,8));fig.subplots_adjust(top=.84,bottom=.18,wspace=.32)
for ax,key,limit in zip(axes,['E_u','E_Q','epsilon_mass'],['velocity_change_limit','flow_change_limit','mass_limit']):
 rows=h['states'] if key=='epsilon_mass' else h['intervals'];ax.plot([r['step'] for r in rows],[r[key] for r in rows],'o-');ax.set_yscale('symlog',linthresh=1e-15)
 ax.axhline(policy[limit],ls='--',color='#b42334');ax.axvline(c['first_full_steady_step'],ls=':',color='#137c66');ax.set_xlabel('Formal monitor timestep');ax.set_ylabel(key);ax.grid(alpha=.2)
save(fig,'optimized_steady_convergence.png','优化以后流场是否仍然正常达到稳态？',f"原 production 联合门槛：连续 5 个区间、每 10 步检查。首次完整通过 {c['first_full_steady_step']}；安全停止 {c['stop_step']}。",'steady_history.json; optimized_steady_candidate.json; frozen policy')
b=read('baseline_runtime');fig,ax=plt.subplots(figsize=(13,8));fig.subplots_adjust(top=.84,bottom=.20)
values=[b['wall_time_s'],c['wall_time_s']];ax.bar(['Stage N observed','Stage O observed'],values,color=['#6685ad','#137c66'],width=.52)
for i,v in enumerate(values):ax.text(i,v*1.025,f'{v:.2f} s\n'+str(b['stop_step'] if i==0 else c['stop_step'])+' steps',ha='center')
ax.set_ylabel('t=0 to safe steady stop: process wall time (s)');ax.set_ylim(0,max(values)*1.22);ax.grid(axis='y',alpha=.2)
save(fig,'gpu_runtime_before_after.png','GPU 优化前后实际 wall-clock time 相差多少？',f"OBSERVATIONAL DEVELOPMENT SPEEDUP：{c['development_speedup']:.3f}×。各一次运行；不是正式 benchmark，也不是科学等价验证。",'baseline_runtime.json; optimized_steady_candidate.json')
(R/'visuals.json').write_text(json.dumps(dict(status='PASS',figures=figures),indent=2,ensure_ascii=False)+'\n');print('Five Stage O development figures written.')
