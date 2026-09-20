"""Measured CUDA/MPI/build outcomes and explicit NOT RUN downstream panels."""
import json, os, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
os.environ['MPLCONFIGDIR']=str(ROOT/'outputs/sv1_3h/plot_cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from sv_validation.provenance import sha256,write_json
R=ROOT/'reports/sv1_3h'
FONT=FontProperties(fname='/mnt/c/Windows/Fonts/msyh.ttc')
plt.rcParams.update({'font.size':12,'figure.facecolor':'#f8fafc','axes.unicode_minus':False})
green='#137c66';red='#b23c40';gray='#64748b';blue='#235b8b';figures=[]
def cn(fig,x,y,text,size=12,color='#243244',ha='left'):
    fig.text(x,y,text,fontproperties=FONT,fontsize=size,color=color,ha=ha,va='center')
def save(fig,name,title,measured,label,caption,bars=0):
    cn(fig,.06,.93,title,21)
    cn(fig,.06,.07,caption,11,gray)
    p=R/name;fig.savefig(p,dpi=150);plt.close(fig)
    figures.append({'path':str(p.relative_to(ROOT)),'sha256':sha256(p),'title':title,'has_measurements':measured,'label':label,'numeric_bars':bars,'caption':caption})
def load(name):return json.loads((R/(name+'.json')).read_text())
assert load('stage_result')['status']=='FAIL'
fig=plt.figure(figsize=(12,7));ax=fig.add_axes([0,0,1,1]);ax.axis('off')
cards=[('保持科学基线','PETSc 3.19.6 原版','11,035 个源文件 SHA 一致\nsvMultiPhysics commit 未变\n生产路线仍为 CPU',blue),
       ('独立安装 CUDA 12.6.3','CUDA kernel：3 / 3 PASS','项目 prefix；仅 toolkit\n驱动与默认 CUDA 13.2 保持原样\nMPI：1 rank 5 / 5；2 ranks PASS',green),
       ('实际兼容性结果','PETSc make：FAIL','原 CUDA 13 三项错误未重现\n新增 tuple.get 接口不兼容\nC++17 clean build 仍然失败',red)]
for i,(title,status,body,color) in enumerate(cards):
    x=.06+i*.31
    ax.add_patch(plt.Rectangle((x,.26),.29,.49,facecolor='white',edgecolor=color,lw=1.5))
    cn(fig,x+.015,.69,title,15,color)
    cn(fig,x+.015,.59,status,13,color)
    cn(fig,x+.015,.42,body,11,gray)
save(fig,'cuda_version_strategy.png','为什么改用 CUDA 12.6，而不是修改 PETSc？',True,'OBSERVED STATUS','更换 toolkit 是本次受控验证变量；CUDA runtime 通过，尚未打通 PETSc GPU 求解。')
fig=plt.figure(figsize=(12,7));ax=fig.add_axes([0,0,1,1]);ax.axis('off')
stages=[('CUDA 12.6.3 toolkit','PASS','驱动 / CUDA13 无变化',green),
 ('MPI + CUDA kernel','PASS','MPI 5+1；kernel 3/3',green),
 ('PETSc configure','PASS','CUDA + sm_89 + C++17',green),
 ('PETSc make','FAIL','Thrust tuple.get 不兼容',red),
 ('PETSc / svMP GPU','NOT RUN','self-test / types / smoke',gray),
 ('血管 proof / benchmark','NOT RUN','无科学或性能结论',gray)]
for i,(name,status,sub,color) in enumerate(stages):
    x=.06+(i%3)*.31;y=.68-(i//3)*.34
    ax.add_patch(plt.Rectangle((x,y-.12),.29,.24,facecolor='white',edgecolor=color,lw=1.6))
    cn(fig,x+.015,y+.055,name,14,color);cn(fig,x+.015,y-.01,status,19,color)
    cn(fig,x+.015,y-.078,sub,10,gray)
save(fig,'gpu_build_pipeline.png','GPU 求解链路现在走通到哪一步？',True,'OBSERVED STATUS','停止点：PETSc make。所有后续阶段保持 NOT RUN；未补丁修复、未升级 PETSc。')
fig=plt.figure(figsize=(12,7));ax=fig.add_axes([.11,.23,.54,.56])
runs=load('cuda12_runtime')['runs'];times=[r['wall_time_s'] for r in runs]
ax.bar(['Run 1','Run 2','Run 3'],times,color=green,width=.54)
ax.set_ylim(0,max(times)*1.25);ax.set_ylabel('Process wall time (s)');ax.spines[['top','right']].set_visible(False)
for i,t in enumerate(times):ax.text(i,t+.007,f'{t:.3f} s',ha='center',fontsize=12)
cn(fig,.71,.69,'3 / 3 PASS',25,green)
cn(fig,.71,.51,'RTX 4090 / sm_89\nCUDA runtime 12.6\n1,024 个元素逐一校验\nlast_error = 0',13,blue)
save(fig,'cuda_runtime_smoke.png','CUDA 12.6 能正常驱动 RTX 4090 吗？',True,'MEASURED','实测包含进程启动、CUDA 初始化、H2D、kernel、D2H 与校验；不是血流 benchmark。',3)
panels=[
 ('petsc_gpu_backend.png','PETSc 的矩阵和向量真的在显卡里吗？','Mat / Vec runtime type：未获得','PETSc library 未构建成功；CUDA kernel 通过不能替代 PETSc Mat/Vec 证据。'),
 ('svmultiphysics_gpu_smoke.png','SimVascular 是否真的调用了 GPU PETSc？','svMultiPhysics GPU build / linkage / official smoke：未运行','没有 GPU PETSc 库；未启动新的 svMultiPhysics build，原 CPU solver 保持原样。'),
 ('gpu_cpu_solution_difference.png','GPU 和 CPU 算出的血流是否相同？','GPU_PROOF_20 与 CPU_PROOF_20_1R：未运行','速度、压力、入口 / 出口流量及质量守恒差异均未测量。'),
 ('gpu_residency.png','求解过程中数据是否一直留在显卡里？','matrix / vectors / PC execution：未测量','没有 PETSc GPU runtime；不能判断矩阵、向量或预条件器的执行位置。'),
 ('gpu_transfer_cost.png','CPU/GPU 搬数据有没有拖慢计算？','H2D / D2H 次数、字节数与耗时：未测量','没有实际求解传输轨迹；不能判断 transfer-bound。'),
 ('gpu_memory.png','RTX 4090 的 24 GB 显存够不够？','真实血管求解 peak VRAM：未测量','设备报告容量 24,564 MiB；此容量不代表求解峰值，也不能证明显存足够。'),
 ('cpu_gpu_runtime.png','同一台服务器上，CPU 和 RTX 4090 谁更快？','CPU 1R / CPU 4R / GPU 1R：无正式 benchmark','未测量中位时间和 speedup；不以 kernel 启动时间预测血流求解性能。')]
for name,title,line,caption in panels:
    fig=plt.figure(figsize=(12,6))
    cn(fig,.5,.62,'NOT RUN',38,gray,'center')
    cn(fig,.5,.43,line,16,blue,'center')
    cn(fig,.5,.28,'前置失败：PETSC_CUDA12_BUILD_FAIL',13,red,'center')
    save(fig,name,title,False,'NOT RUN',caption)
write_json(R/'visuals.json',{'figures':figures,'unknown_measurements_are_null':True,'human_review_status':'PENDING','stage_status':'FAIL'})
print('10 figures: 3 observed panels; 7 explicit NOT RUN panels')
