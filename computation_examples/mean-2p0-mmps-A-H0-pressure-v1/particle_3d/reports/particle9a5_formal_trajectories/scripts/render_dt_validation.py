"""Human-reviewable dt revision: fixed first two births, actual saved states only."""
from pathlib import Path
import json,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import render_formal as v
import matplotlib.pyplot as plt
from particle_3d.formal_cohort_p9a5 import digest,canonical
R=v.R;root=v.ROOT
fig,axes=plt.subplots(1,3,figsize=(14,4.6));records=[];old_times=[];new_times=[]
for i,color in [(1,'#1767a6'),(2,'#d56821')]:
 old=root/'particle_3d/reports/particle9a4_population_inlet/outputs/smoke30/trajectories'/f'mb_{i:06d}.npz'
 folder=R/'validation/local_dt1ms_worker_parity/1'/str(i);new=folder/'trajectory.npz';a=np.load(old)['samples'];b=np.load(new)['samples'];m=json.loads((folder/'metrics.json').read_text())
 t=np.linspace(0,min(a[-1,0],b[-1,0]),501);x=np.column_stack([np.interp(t,a[:,0],a[:,j]) for j in [1,2,3]]);y=np.column_stack([np.interp(t,b[:,0],b[:,j]) for j in [1,2,3]])
 error=np.linalg.norm(x-y,axis=1)*1e6
 axes[0].plot(a[:,1]*1e6,a[:,2]*1e6,ls='--',color=color,alpha=.5,label=f'ID {i} · 0.25 ms');axes[0].plot(b[:,1]*1e6,b[:,2]*1e6,color=color,label=f'ID {i} · 1.0 ms')
 axes[1].plot(t*1000,error,color=color,label=f'ID {i}')
 old_times.append(a[-1,0]*1000);new_times.append(b[-1,0]*1000)
 records.append(dict(particle_id=i,old_samples_sha256=digest(old),new_samples_sha256=digest(new),new_outlet=m['outlet'],old_age_s=float(a[-1,0]),new_age_s=float(b[-1,0]),maximum_interpolated_position_difference_um=float(error.max()),safety={k:m[k] for k in ['penetration_count','handoff_violation_count','inlet_escape_count','nan_inf_count','unclassified_corruption_count']}))
axes[0].set(xlabel='x (µm)',ylabel='y (µm)',title='相同出生样本的实际路径（XY投影）');axes[0].legend(frameon=False,fontsize=8)
axes[1].set(xlabel='共同轨迹年龄 (ms)',ylabel='Interpolated position difference (µm)',title='相同物理时刻的位置差异');axes[1].legend(frameon=False)
x=np.arange(2);axes[2].bar(x-.18,old_times,.36,color='#91b0cb',label='0.25 ms');axes[2].bar(x+.18,new_times,.36,color='#d68a50',label='1.0 ms');axes[2].set(xticks=x,xticklabels=['ID 1','ID 2'],ylabel='Transit time (ms)',title='实际离开时间');axes[2].legend(frameon=False)
fig.suptitle('步长调整验证 · dt = 1.0 ms · 固定前两条样本的实轨迹检查',fontsize=16)
fig.text(.5,.01,'本地验证样本，非正式总体统计；多进程及同 dt 原积分器逐点一致。与旧 dt 的差异不构成时间收敛证明。',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.05,1,.93));v.save(fig,'Stage_B_dt1ms_validation')
(R/'data/dt_validation_figure.json').write_bytes(canonical(dict(dt_s=.001,particle_ids=[1,2],selection='FIRST_TWO_FROZEN_CORE500_MEMBERS',rows=records)))
p=R/'OPEN_RESULTS.html';s=p.read_text();section='''<section><h2>1.0 ms 本地验证通过</h2><p>新增永久测试23项通过；相同样本的1与2 workers科学结果一致；新调度器与原积分器在同一1.0 ms步长下逐点一致。服务器固定24条样本的worker scaling仍在运行。</p><img src="figures/Stage_B_dt1ms_validation.png"><p><a href="figures/Stage_B_dt1ms_validation.pdf">验证图 PDF</a> · <a href="data/local_dt1ms_validation.json">多进程一致性证据</a> · <a href="logs/new_tests.txt">23项测试日志</a> · <a href="data/dt_revision.json">步长修订记录</a></p></section>'''
s=s.replace('</html>',section+'</html>');p.write_text(s)
