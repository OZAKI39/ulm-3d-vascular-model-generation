#!/usr/bin/env python3
"""Observed MPI/build outcomes and explicit NOT RUN panels for unmet GPU gates."""
import json,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
os.environ['MPLCONFIGDIR']=str(ROOT/'outputs/sv1_3g/plot_cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from sv_validation.provenance import sha256,write_json
R=ROOT/'reports/sv1_3g';FONT=FontProperties(fname='/mnt/c/Windows/Fonts/msyh.ttc')
plt.rcParams.update({'font.size':12,'figure.facecolor':'#f8fafc','axes.unicode_minus':False})
green='#137c66';red='#b23c40';gray='#64748b';blue='#235b8b';figures=[]
def cn(fig,x,y,text,size=12,color='#243244',ha='left'):
    fig.text(x,y,text,fontproperties=FONT,fontsize=size,color=color,ha=ha,va='center')
def save(fig,name,title,measured,label,caption,bars=0):
    cn(fig,.06,.93,title,20)
    cn(fig,.06,.07,caption,11,gray)
    p=R/name;fig.savefig(p,dpi=150);plt.close(fig)
    figures.append({'path':str(p.relative_to(ROOT)),'sha256':sha256(p),'title':title,'has_measurements':measured,'label':label,'numeric_bars':bars,'caption':caption})
system=json.loads((R/'system_mpi_tests.json').read_text());gate=json.loads((R/'mpi_hard_gate.json').read_text())
fig=plt.figure(figsize=(12,7));ax=fig.add_axes([.06,.20,.88,.62]);ax.axis('off')
labels={'direct_true':'普通 /bin/true','direct_hostname':'普通 hostname','direct_singleton':'MPI singleton','default_true':'默认 launcher + true','default_hostname':'默认 launcher + hostname','default_hello':'默认 launcher + hello','root_true':'显式 root + true','root_hello':'显式 root + hello','no_binding':'关闭 binding','minimal_transport':'ob1 / self,tcp','tmpdir':'独立 TMPDIR'}
rows=[]
for t in system['tests']:
    status='TIMEOUT' if t['timeout'] else 'PASS' if t['exit_code']==0 else 'ROOT REFUSED'
    rows.append([labels[t['name']],status,f"{t['wall_time_s']:.3f} s"])
rows.extend([['独立 Open MPI：singleton','PASS',f"{gate['singleton']['wall_time_s']:.3f} s"],['独立 Open MPI：1 rank × 5','5 / 5 PASS',f"{min(x['wall_time_s'] for x in gate['rank1']):.3f}–{max(x['wall_time_s'] for x in gate['rank1']):.3f} s"],['独立 Open MPI：2 ranks','PASS',f"{gate['rank2']['wall_time_s']:.3f} s"]])
tab=ax.table(cellText=rows,colLabels=['实测测试','结果','墙钟时间'],cellLoc='left',colWidths=[.48,.28,.24],loc='center');tab.auto_set_font_size(False);tab.set_fontsize(11);tab.scale(1,1.6)
for (i,j),cell in tab.get_celld().items():
    cell.get_text().set_fontproperties(FONT);cell.set_edgecolor('#dce3eb')
    if i==0:cell.set_facecolor('#dfe9f3');cell.get_text().set_color(blue)
    elif j==1:cell.get_text().set_color(green if 'PASS' in rows[i-1][1] else red)
save(fig,'mpi_unblock.png','GPU 路线最开始卡在哪里，现在修好了吗？',True,'MEASURED','原系统 MPI 仍不可用；项目独立 Open MPI 4.1.6 通过验收；每个启动测试设置 10 秒超时。')

fig=plt.figure(figsize=(12,6));ax=fig.add_axes([0,0,1,1]);ax.axis('off')
stages=[('MPI runtime','PASS','1 rank 5/5；2 ranks PASS',green),('PETSc CUDA configure','PASS','PETSc 3.19.6 + CUDA 13.2',green),('PETSc CUDA compile','BLOCKED','CUDA / Thrust API 不兼容',red),('PETSc GPU smoke','NOT RUN','依赖 CUDA PETSc 库',gray),('svMultiPhysics GPU','NOT RUN','构建与官方 smoke 未启动',gray),('vascular proof','NOT RUN','GPU / CPU 20 步均未启动',gray)]
for i,(name,status,sub,color) in enumerate(stages):
    x=.06+(i%3)*.315;y=.67-(i//3)*.35
    ax.add_patch(plt.Rectangle((x,y-.12),.285,.24,facecolor='white',edgecolor=color,lw=1.6))
    cn(fig,x+.018,y+.055,name,13,color);cn(fig,x+.018,y-.007,status,18,color);cn(fig,x+.018,y-.072,sub,10,gray)
    if i%3<2:ax.annotate('',xy=(x+.31,y),xytext=(x+.287,y),arrowprops={'arrowstyle':'->','color':gray})
save(fig,'gpu_build_pipeline.png','GPU 求解链路已经走到哪一步？',True,'OBSERVED STATUS','实际停止点：PETSc make；未生成可用 CUDA PETSc 库，后续关卡未运行。')

panels=[
 ('cuda_types.png','PETSc 的矩阵和向量真的在显卡上吗？','Mat / Vec / KSP 类型：尚无 runtime evidence','configure 通过不等于 CUDA Mat/Vec 已运行。'),
 ('gpu_residency.png','求解过程中数据是否一直留在显卡里？','Matrix / vectors / PC apply 位置：未测量','GPU KSP 尚未运行；不能宣称整个求解过程驻留显卡。'),
 ('gpu_transfer_cost.png','CPU 和 GPU 之间搬数据是否成为瓶颈？','H2D / D2H 次数、字节数、耗时：未测量','没有传输轨迹；不能判为 transfer-bound，也不能判为无传输瓶颈。'),
 ('gpu_cpu_solution_difference.png','GPU 加速后，答案有没有改变？','GPU_PROOF_20 与 CPU_PROOF_20：均未运行','没有速度、压力或边界流量的 CPU/GPU 差异数据。'),
 ('cpu_gpu_runtime.png','同一台服务器上，CPU 和 RTX 4090 谁更快？','CPU 1R / CPU 4R / GPU 1R：无本阶段 benchmark','没有实际加速比；MPI 启动改善不能用来推断血流求解加速。'),
 ('gpu_memory.png','24 GB 显存够不够？','血流求解 peak GPU memory：未测量','设备容量可见；求解器未运行，无法判断血流计算的显存需求。')]
for name,title,line,caption in panels:
    fig=plt.figure(figsize=(12,5.5))
    cn(fig,.5,.60,'NOT RUN',36,gray,'center')
    cn(fig,.5,.41,line,16,blue,'center')
    cn(fig,.5,.27,'前置阻塞：CUDA_TOOLKIT_COMPATIBILITY',12,red,'center')
    save(fig,name,title,False,'NOT RUN',caption)
write_json(R/'visuals.json',{'figures':figures,'unknown_measurements_are_null':True,'human_review_status':'PENDING','stage_status':'BLOCKED'})
print('8 figures: 2 observed outcome panels, 6 explicit NOT RUN panels')
