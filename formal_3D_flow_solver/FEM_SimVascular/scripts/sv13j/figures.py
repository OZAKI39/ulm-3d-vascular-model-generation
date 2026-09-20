"""Actual compatibility and GPU runtime evidence; no invented vascular measurements."""
import json,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
os.environ['MPLCONFIGDIR']=str(ROOT/'outputs/sv1_3j/plot_cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from sv_validation.provenance import sha256,write_json
R=ROOT/'reports/sv1_3j';FONT=FontProperties(fname='/mnt/c/Windows/Fonts/msyh.ttc')
plt.rcParams.update({'figure.facecolor':'#f8fafc','font.size':12,'axes.unicode_minus':False})
green='#137c66';red='#b23c40';gray='#64748b';blue='#235b8b';figures=[]
def load(n):return json.loads((R/(n+'.json')).read_text())
def cn(fig,x,y,text,size=12,color='#243244',ha='left'):
    fig.text(x,y,text,fontproperties=FONT,fontsize=size,color=color,ha=ha,va='center')
def save(fig,name,title,measured,label,caption,bars=0):
    cn(fig,.055,.93,title,21);cn(fig,.055,.065,caption,11,gray)
    p=R/name;fig.savefig(p,dpi=150);plt.close(fig)
    figures.append({'path':str(p.relative_to(ROOT)),'sha256':sha256(p),'title':title,'has_measurements':measured,'label':label,'numeric_bars':bars,'caption':caption})
def table(rows,cols,widths):
    fig=plt.figure(figsize=(12,7));ax=fig.add_axes([.055,.22,.90,.58]);ax.axis('off')
    tab=ax.table(cellText=rows,colLabels=cols,colWidths=widths,loc='center',cellLoc='left');tab.auto_set_font_size(False);tab.set_fontsize(12);tab.scale(1,2.4)
    for (i,j),cell in tab.get_celld().items():
        cell.get_text().set_fontproperties(FONT);cell.set_edgecolor('#dce3eb')
        if i==0:cell.set_facecolor('#dfe9f3');cell.get_text().set_color(blue)
        elif 'FAIL' in cell.get_text().get_text():cell.get_text().set_color(red)
        elif 'PASS' in cell.get_text().get_text():cell.get_text().set_color(green)
    return fig
a=load('compatibility_matrix')['candidates'][0]
rows=[['12.6.3（历史 H）','2.5.0','13.3.0','3/3 PASS','PASS','FAIL'],
      ['12.3.2（本阶段）','2.2.0',a['host_compiler'],'3/3 PASS','PASS','PASS'],
      ['12.2.2','NOT RUN','NOT RUN','NOT RUN','NOT RUN','NOT RUN'],
      ['12.1.1','NOT RUN','NOT RUN','NOT RUN','NOT RUN','NOT RUN']]
fig=table(rows,['CUDA','Thrust','GCC','kernel','configure','make'],[.23,.14,.14,.18,.17,.14])
cn(fig,.055,.17,'Winner：CUDA 12.3.2；后续两个候选 NOT_REQUIRED，未下载、未安装、未构建。',13,green)
save(fig,'compatibility_matrix.png','哪一版 CUDA 能和 PETSc 3.19.6 正常一起编译？',True,'OBSERVED STATUS','历史 H 作为对照，未重复构建。当前 PETSc 原版源文件修改数 = 0。')
rows=[['13.2（历史 G）','PETSc make FAIL','clockRate / memoryClockRate / unary_function'],
      ['12.6.3（历史 H）','PETSc make FAIL','thrust::tuple 没有成员 get'],
      ['12.3.2（本阶段）','PETSc GPU PASS','已越过两类旧编译错误；GPU smoke 3/3'],
      ['当前 svMultiPhysics','官方 smoke FAIL','MPI_Bcast: MPI_ERR_TYPE；未进入求解'],
      ['12.2.2 / 12.1.1','NOT RUN','已找到首个 winner，不再搜索']]
fig=table(rows,['版本 / 层次','实测状态','主要证据'],[.23,.25,.52])
save(fig,'compatibility_failure_signatures.png','不同 CUDA 版本分别卡在哪里？',True,'OBSERVED STATUS','当前新增失败位于应用所需 MPI datatype 支持；不能归因于未通过验证的 CUDA/PETSc 组合。')
fig=plt.figure(figsize=(12,7));ax=fig.add_axes([0,0,1,1]);ax.axis('off')
stages=[('GCC12 + CUDA12.3.2','PASS','host smoke / kernel / MPI',green),('PETSc build + self-test','PASS','原版 3.19.6，C++17，sm_89',green),
        ('PETSc GPU sparse solve','PASS','CUDA Mat / Vec；3/3',green),('svMultiPhysics build / link','PASS','固定 commit 与 winner PETSc',green),
        ('官方 fluid GPU smoke','FAIL','MPI Fortran 类型不可用',red),('血管 proof / 性能','NOT RUN','遵守前置 smoke hard gate',gray)]
for i,(name,status,sub,color) in enumerate(stages):
    x=.055+(i%3)*.312;y=.68-(i//3)*.34
    ax.add_patch(plt.Rectangle((x,y-.12),.29,.24,facecolor='white',edgecolor=color,lw=1.5))
    cn(fig,x+.012,y+.058,name,13,color);cn(fig,x+.012,y-.011,status,19,color);cn(fig,x+.012,y-.078,sub,10,gray)
save(fig,'gpu_build_pipeline.png','GPU 求解链路目前已经走到哪一步？',True,'OBSERVED STATUS','原 MPI C rank 测试通过，应用所需 Fortran 预定义类型仍缺失；两项结论分别保留。')
runs=load('petsc_gpu_smoke')['runs']
rows=[[str(i),r['mat_type'],r['vec_type'],str(r['iterations']),f"{r['true_relative_residual']:.3e}",f"{r['solution_error_inf']:.3e}"] for i,r in enumerate(runs,1)]
fig=table(rows,['Run','实际 Mat','实际 Vec','GMRES 次数','真实相对残差','解误差 inf'],[.07,.23,.14,.15,.20,.21])
cn(fig,.055,.19,'3 / 3 PASS；KSP reason = 2（CONVERGED_RTOL）；PC = Jacobi。',14,green)
save(fig,'petsc_gpu_backend.png','PETSc 的矩阵和向量真的进入 RTX 4090 了吗？',True,'MEASURED','这是 128 维稀疏系统的实际 GPU 证据；不能替代血管 ASM/ILU(2) 的驻留或性能测量。')
fig=table([['同 commit 原版构建','PASS','未改 solver source'],['链接 winner GPU PETSc','PASS','ldd / readelf / CMakeCache 一致'],
    ['官方 fluid/newtonian，1 rank','FAIL','exit 3；MPI_ERR_TYPE；未输出 VTU'],['C 类型 Bcast，1 / 2 ranks','PASS','MPI_INT / DOUBLE / CHAR / CXX_BOOL'],
    ['Fortran 类型 Bcast，1 / 2 ranks','FAIL','INTEGER / DOUBLE_PRECISION / CHARACTER 等']],
    ['核验内容','状态','证据'],[.36,.14,.50])
save(fig,'svmultiphysics_gpu_smoke.png','SimVascular 已经真正调用 GPU PETSc 了吗？',True,'OBSERVED STATUS','链接已通过；应用在初始化广播失败，尚无 svMultiPhysics runtime CUDA Mat/Vec 证据。')
panels=[
 ('gpu_cpu_solution_difference.png','CPU 和 GPU 算出的答案是否一致？','GPU / CPU t=0 → step20：均未运行','速度、压力、Qin、各 Qout、mass error 与 max velocity 差异均未测量。'),
 ('gpu_residency.png','计算数据有没有一直留在显卡里？','血管 matrix / vectors / MatMult / ASM / ILU(2)：未测量','PETSc 小系统使用 Jacobi；不能据此推断血管 ASM/ILU(2) 的执行位置。'),
 ('gpu_transfer_cost.png','CPU/GPU 搬数据是不是瓶颈？','H2D / D2H 数量、字节、耗时与迭代增长：未测量','未做血管 profiling；不能判为 transfer-bound 或无传输瓶颈。'),
 ('gpu_memory.png','24 GB 显存实际用了多少？','真实血管求解 peak VRAM：未测量','24,564 MiB 是设备容量；未使用容量、初始化占用或小系统占用替代峰值。'),
 ('cpu_gpu_runtime.png','同一台服务器上 GPU 到底快多少？','CPU1 / CPU4 / GPU1 的正式中位时间：未测量','没有 speedup；不使用小系统 smoke 时间估计真实血管速度。')]
for name,title,line,caption in panels:
    fig=plt.figure(figsize=(12,6));cn(fig,.5,.62,'NOT RUN',38,gray,'center');cn(fig,.5,.43,line,16,blue,'center')
    cn(fig,.5,.28,'前置失败：SVMULTIPHYSICS_GPU_FAIL',13,red,'center')
    save(fig,name,title,False,'NOT RUN',caption)
write_json(R/'visuals.json',{'figures':figures,'stage_status':'FAIL','human_review_status':'PENDING','unknown_measurements_are_null':True})
print('10 figures: 5 actual evidence panels; 5 explicit NOT RUN panels')
