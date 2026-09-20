#!/usr/bin/env python3
"""Evidence-driven review figures; unavailable GPU data never becomes a zero."""
import os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
os.environ['MPLCONFIGDIR']=str(ROOT/'outputs/sv1_3/plot_cache')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from sv_validation.sv13 import *
from sv_validation.provenance import write_json,sha256

FONT=FontProperties(fname='/mnt/c/Windows/Fonts/msyh.ttc')
plt.rcParams.update({'font.size':11,'figure.facecolor':'white','axes.unicode_minus':False})
ref=load('reference_freeze');cpu=load('cpu_execution');perf=load('cpu_performance')
eq=load('cpu_reference_equivalence');qc=load('cpu_saved_states');gpu=load('gpu_decision')
selection=load('production_selection');figures=[]
blue='#2878a4';green='#278873';gray='#a3adb8';red='#c44547'
def save(fig,name,title,caption,kind='measured CPU evidence'):
    fig.suptitle(title,fontproperties=FONT,fontsize=18,y=.97)
    fig.text(.5,.015,caption,fontproperties=FONT,fontsize=10,ha='center',va='bottom')
    fig.tight_layout(rect=(.02,.12,.98,.90));path=REPORT/name
    fig.savefig(path,dpi=160);plt.close(fig)
    figures.append({'path':str(path.relative_to(ROOT)),'sha256':sha256(path),'kind':kind,'caption':caption})
def unavailable(ax,message,sub=''):
    ax.axis('off');ax.text(.5,.6,message,fontproperties=FONT,fontsize=18,ha='center',va='center',wrap=True)
    if sub:ax.text(.5,.3,sub,fontproperties=FONT,fontsize=11,ha='center',va='center',wrap=True)

fig,axes=plt.subplots(1,2,figsize=(13,6))
pcts=np.array([r['linear_time_percent'] for r in cpu['linear'][1:]])
axes[0].hist(pcts,bins=np.arange(-.5,101.5,1),color=blue)
axes[0].set(xlabel='Native solve time / interval time (%)',ylabel='Linear solves',xlim=(max(0,pcts.min()-3),100.5))
axes[0].axvline(np.median(pcts),color=red,ls='--',label=f'Median: {np.median(pcts):g}%');axes[0].legend()
axes[1].axis('off')
rows=[['FEM assembly','Not separately measured'],['KSP + preconditioner',f'Native share median {np.median(pcts):g}%'],['Preconditioner alone','Included; no separate timing'],['I/O','Not separately measured'],['Other','Not separately measured']]
table=axes[1].table(cellText=rows,colLabels=['Category','Available evidence'],loc='center',cellLoc='left',colWidths=[.43,.57]);table.auto_set_font_size(False);table.set_fontsize(10);table.scale(1,2.1)
save(fig,'where_time_goes.png','现在一次血流计算，时间主要花在哪里？',
     '排除继承计时器影响的首条记录；原生占比含预条件、KSP 设置与同步；未测分项保持未知')

fig,axes=plt.subplots(1,2,figsize=(12,6));names=['SV1.2 reference','CPU automatic stop']
for ax,vals,label in [(axes[0],[perf['reference_endpoint'],perf['candidate_endpoint']],'Final timestep'),
                      (axes[1],[perf['reference_measured_wall_s']/3600,perf['candidate_measured_wall_s']/3600],'Measured continuation hours')]:
    bars=ax.bar(names,vals,color=[gray,blue]);ax.bar_label(bars,fmt='%.3g');ax.set_ylabel(label);ax.margins(y=.2)
save(fig,'early_stop_savings.png','达到稳态后自动停止，可以少算多少？',
     f"相同 step10 原生起点；首次合格 step{load('cpu_validation')['first_qualifying_step']}，实际安全停止 step{cpu['last_step']}；时间均为实测")

fig,ax=plt.subplots(figsize=(12,6))
runtime_path=REPORT/'benchmark_repeatability.json'
runtime=load('benchmark_repeatability') if runtime_path.exists() else {'status':'NOT_RUN'}
if runtime.get('status')=='PASS':
    values=[runtime[k]['reported_wall_seconds'] for k in ('CPU_1R','CPU_4R','GPU_1R')]
    bars=ax.bar(['CPU 1 rank','CPU 4 ranks','GPU 1 rank'],values,color=[gray,blue,green]);ax.bar_label(bars,fmt='%.1f');ax.set_ylabel('Seconds / fixed 20 steps')
else:
    unavailable(ax,'同机 20 步 CPU / GPU 比较未完成',gpu['reason']+'\n没有可报告的 GPU 加速比')
save(fig,'cpu_gpu_runtime.png','同一个计算，CPU 和 RTX 4090 谁更快？',
     '只接受同服务器、同原生 checkpoint、同 20 步窗口；未运行的数据不画成零',kind='benchmark evidence or explicit unavailability')

fig,ax=plt.subplots(figsize=(13,6));ax.axis('off')
labels=['CPU assembly','H2D','GPU KSP','D2H']
for i,label in enumerate(labels):
    x=.04+i*.245
    ax.text(x+.09,.60,label,ha='center',va='center',fontsize=15,bbox=dict(boxstyle='round,pad=.8',fc='#eef2f5',ec=gray))
    if i<3:ax.annotate('',xy=(x+.235,.60),xytext=(x+.19,.60),arrowprops={'arrowstyle':'->','color':gray})
ax.text(.5,.25,'目标数据流示意；实际 GPU 数据驻留尚未获得验证' if gpu['status']!='ADOPTED' else '详见实际 GPU 数据驻留审核',fontproperties=FONT,ha='center',fontsize=15)
save(fig,'gpu_residency.png','计算过程中数据是否一直留在显卡里？',gpu['reason'],kind='unverified target flow explicitly labeled')

fig,ax=plt.subplots(figsize=(12,6))
transfer_path=REPORT/'gpu_transfer_audit.json';transfer=load('gpu_transfer_audit') if transfer_path.exists() else {'status':'NOT_RUN'}
if transfer.get('status')=='PASS':
    vals=[transfer[k] for k in ('transfer_time_s','KSP_compute_time_s','PC_time_s')]
    bars=ax.bar(['Transfer','KSP compute','Preconditioner'],vals,color=[gray,blue,green]);ax.bar_label(bars,fmt='%.3g');ax.set_ylabel('Seconds; separate profiling run')
else:unavailable(ax,'搬运时间、GPU 计算时间和预条件时间：未测量',gpu['reason'])
save(fig,'gpu_transfer_cost.png','CPU 和 GPU 之间搬数据花了多少时间？',
     '未执行有效 GPU KSP 时不能推断搬运开销；GPU 利用率不能替代驻留证据',kind='profiling evidence or explicit unavailability')

baseline=load('solver_history','sv1_2')
fig,ax=plt.subplots(figsize=(12,6));old=[r['linear_iterations'] for r in baseline['linear'] if r['step']>cpu['initial_step']]
matched=[r['linear_iterations'] for r in baseline['linear'] if cpu['initial_step']<r['step']<=cpu['last_step']]
new=[r['linear_iterations'] for r in cpu['linear']]
ax.boxplot([old,matched,new],tick_labels=['SV1.2 steps 11-400','SV1.2 matched steps 11-80','CPU stop steps 11-80'],showfliers=False)
ax.set_ylabel('KSP iterations per linear solve');ax.grid(axis='y',alpha=.2)
save(fig,'ksp_iterations_comparison.png','优化后，每个线性问题需要多少次迭代？',
     '同一 PETSc 配置；中间列限定相同时间步，便于公平比较；只缩短终止时间，不声称预条件器变快')

fig,axes=plt.subplots(1,3,figsize=(15,6))
errors=eq['errors'];vals=[errors['velocity_relative_volume_L2'],errors['pressure_relative_volume_L2']]
axes[0].bar(['Velocity L2','Pressure L2'],vals,color=blue);axes[0].axhline(1e-5,color=red,ls='--',label='Limit 1e-5')
axes[0].set_yscale('symlog',linthresh=1e-16);axes[0].set_ylim(bottom=0);axes[0].legend();axes[0].set_ylabel('Relative difference')
roles=['OUTLET_01','OUTLET_02','OUTLET_03'];x=np.arange(3)
axes[1].bar(x-.18,[eq['reference']['outlet_fractions'][r] for r in roles],width=.36,color=gray,label='SV1.2')
axes[1].bar(x+.18,[eq['candidate']['outlet_fractions'][r] for r in roles],width=.36,color=blue,label='CPU stop')
axes[1].set_xticks(x,['Outlet 1','Outlet 2','Outlet 3']);axes[1].set_ylabel('Flow fraction');axes[1].legend()
axes[2].bar(['SV1.2','CPU stop'],[eq['reference']['epsilon_mass'],eq['candidate']['epsilon_mass']],color=[gray,blue])
axes[2].axhline(1e-6,color=red,ls='--',label='Limit 1e-6');axes[2].set_yscale('symlog',linthresh=1e-16);axes[2].set_ylim(0,1e-5);axes[2].set_ylabel('Mass error');axes[2].legend()
axes[2].scatter([0,1],[eq['reference']['epsilon_mass'],eq['candidate']['epsilon_mass']],color=[gray,blue],zorder=3,clip_on=False)
for i,key in enumerate(['reference','candidate']):axes[2].annotate(f"{eq[key]['epsilon_mass']:.2g}",(i,eq[key]['epsilon_mass']),xytext=(0,9),textcoords='offset points',ha='center')
save(fig,'solution_equivalence.png','加速之后，得到的流场还是同一个吗？',
     '真实原生 VTU 独立积分；压力不平移、流量不重标定；CPU 等价门槛全部通过；GPU 结论以实际证据为准')

fig,axes=plt.subplots(1,2,figsize=(13,6));orig=load('steady_stop_replay')['rows']
for ax,k,threshold in zip(axes,['E_u','E_Q'],[1e-5,1e-6]):
    r=[r['interval'] for r in orig if r['interval']];ax.plot([i['step'] for i in r],[i[k] for i in r],color=gray,label='SV1.2')
    ax.plot([i['step'] for i in qc['intervals']],[i[k] for i in qc['intervals']],'o-',color=blue,label='CPU automatic stop')
    ax.axhline(threshold,color=red,ls='--',label=f'Limit {threshold:g}');ax.set_yscale('symlog',linthresh=1e-16);ax.set_ylim(bottom=0)
    ax.set(xlabel='Saved timestep',ylabel=k);ax.grid(alpha=.2);ax.legend(fontsize=9)
save(fig,'steady_equivalence.png','提前停止和 GPU 计算仍然达到同样的稳态吗？',
     'CPU：连续五个保存区间及最终质量门槛均通过；GPU：没有有效结果时保持未验证')

fig,axes=plt.subplots(1,2,figsize=(12,6));rss=[ref['baseline_resources']['peak_single_process_RSS_KiB']/1024,cpu['peak_single_process_RSS_KiB']/1024]
bars=axes[0].bar(['SV1.2','CPU stop'],rss,color=[gray,blue]);axes[0].bar_label(bars,fmt='%.1f');axes[0].set_ylabel('Peak individual process RSS (MiB)');axes[0].margins(y=.2)
mem_path=REPORT/'gpu_memory.json';mem=load('gpu_memory') if mem_path.exists() else {'status':'NOT_RUN'}
if mem.get('status')=='PASS':
    axes[1].bar(['GPU peak'],[mem['peak_bytes']/2**30],color=green);axes[1].axhline(.9*mem['physical_bytes']/2**30,color=red,ls='--');axes[1].set_ylabel('GiB')
else:unavailable(axes[1],'GPU 峰值显存：未测量','设备总显存不等于求解峰值显存')
save(fig,'memory_usage.png','CPU 内存和 GPU 显存是否足够？','CPU RSS 来自 GNU time，为最大单进程值，不能当作所有 MPI ranks 内存之和')

fig,axes=plt.subplots(1,2,figsize=(13,6))
bars=axes[0].bar(['SV1.2 reference','CPU automatic stop'],[perf['reference_measured_wall_s']/60,perf['candidate_measured_wall_s']/60],color=[gray,blue]);axes[0].bar_label(bars,fmt='%.1f');axes[0].set_ylabel('Measured continuation minutes');axes[0].margins(y=.2)
axes[1].axis('off');axes[1].text(.5,.65,f"{perf['measured_continuation_speedup']:.2f}×",ha='center',fontsize=52,color=blue)
axes[1].text(.5,.32,f"含历史前 10 步成本：{perf['accounted_total_speedup']:.2f}×\n最终选择：{selection['selected']}",fontproperties=FONT,ha='center',fontsize=12)
save(fig,'final_speedup.png','最终一次正式计算快了多少？','主加速比：相同 step10 起点的两次实际续算；含起步成本的总计另列，不冒充新测 t=0 完整运行')
write_json(REPORT/'visuals.json',{'figures':figures,'unavailable_measurements_not_plotted_as_zero':True,
    'reference_sha256':ref['accepted_solution']['sha256'],'candidate_sha256':eq['candidate_sha256'],
    'rendering':'Matplotlib Agg CPU; no flow field edits','human_review_pending':True})
print('SV1.3 review figures: '+str(len(figures)),flush=True)
