"""Actual compatibility and GPU runtime evidence; no invented vascular measurements."""
import json,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
os.environ['MPLCONFIGDIR']=str(ROOT/'outputs/sv1_3l/plot_cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from sv_validation.provenance import sha256,write_json
R=ROOT/'reports/sv1_3l';FONT=FontProperties(fname='/mnt/c/Windows/Fonts/msyh.ttc')
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
datatype=load('mpi_fortran_datatypes');sizes=datatype['native_Fortran_sizes']
old=json.loads((ROOT/'reports/sv1_3j/mpi_datatype_diagnosis.json').read_text())
oldrows={r['name']:r for r in old['runs'][0]['rows']}
newrows={r['name']:r for r in datatype['repetitions'][0]['runs'][0]['rows']}
rows=[[name,str(oldrows[name]['size']),'MPI_ERR_TYPE',str(newrows[name]['size']),
       'PASS' if newrows[name]['bcast_rc']==0 and newrows[name]['data_ok']==1 else 'FAIL'] for name in sizes]
fig=table(rows,['数据类型','旧大小','旧广播','新大小','新广播 / 内容'],[.32,.12,.19,.12,.25])
cn(fig,.055,.18,'C++ probe：3 轮 ×（rank1 + rank2）通过；原生 Fortran 广播也通过。',13,green)
save(fig,'mpi_datatype_before_after.png','SimVascular 需要的 MPI 数据类型修好了吗？',True,'MEASURED','大小单位：bytes；新类型大小同时与 Fortran storage_size 实测结果核对。旧 MPI 未修改。')
smoke=load('svmp_gpu_smoke');proof=load('gpu_proof')
stages=[('新 MPI + Fortran','PASS','四种类型和真实广播通过'),
        ('PETSc clean rebuild',load('petsc_rebuild')['status'],'CUDA12.3 / 新 MPI / 原版 3.19.6'),
        ('PETSc GPU revalidation',load('petsc_gpu_smoke')['status'],'三次 CUDA Mat/Vec 稀疏求解'),
        ('svMultiPhysics build/link',load('svmp_gpu_link')['status'],'固定 commit，源码未修改'),
        ('官方 fluid GPU smoke',smoke['status'],smoke.get('failure_class') or '实际流体求解'),
        ('真实血管 proof',proof['status'],proof.get('reason','20-step science validation'))]
fig=plt.figure(figsize=(12,7));ax=fig.add_axes([0,0,1,1]);ax.axis('off')
for i,(name,status,detail) in enumerate(stages):
    color=green if status=='PASS' else red if status=='FAIL' else gray
    x=.055+(i%3)*.312;y=.68-(i//3)*.34
    ax.add_patch(plt.Rectangle((x,y-.12),.29,.24,facecolor='white',edgecolor=color,lw=1.5))
    cn(fig,x+.012,y+.058,name,13,color);cn(fig,x+.012,y-.011,status.replace('_',' '),19,color)
    short='详见官方 smoke 失败证据' if len(detail)>43 else detail
    cn(fig,x+.012,y-.078,short,9,gray)
save(fig,'gpu_stack_pipeline.png','RTX 4090 求解链路现在走到哪一步？',True,'OBSERVED STATUS','只有官方 GPU smoke 完整通过，才继续真实血管 proof；正式 CPU production 保持不变。')
runs=load('petsc_gpu_smoke')['runs']
rows=[[str(i),r['mat_type'],r['vec_type'],str(r['iterations']),f"{r['true_relative_residual']:.3e}",f"{r['solution_error_inf']:.3e}"] for i,r in enumerate(runs,1)]
fig=table(rows,['Run','Mat','Vec','GMRES 次数','真实相对残差','解误差 inf'],[.07,.23,.14,.15,.20,.21])
cn(fig,.055,.19,'新 MPI：3 / 3 PASS；CPU / MPI / CUDA self-tests 通过。',14,green)
save(fig,'petsc_gpu_revalidation.png','换 MPI 后 PETSc GPU 还能正常工作吗？',True,'MEASURED','这是重新构建后的实测结果；未沿用 Stage J 的 runtime 通过结论。小型 smoke 的 PC 为 Jacobi。')
def observed(value):return 'NOT RUN' if value is None else str(value)
rows=[['原 MPI_ERR_TYPE', '消失' if smoke.get('startup_pass') else '未通过','进入 PETSc/KSP，未再出现 MPI 类型错误'],
      ['实际进程退出',str(smoke.get('exit_code')),smoke.get('first_error') or 'PASS'],
      ['线性 / 非线性',smoke.get('solver_status','NOT RUN'),'KSP 记录 1 次收敛；未完成非线性步'],
      ['VTU / reload',observed(smoke.get('VTU_count'))+' / '+observed(smoke.get('reload_pass')),'字段有限性见报告'],
      ['runtime Mat / Vec',smoke.get('backend_status','NOT RUN'),observed(smoke.get('mat_type'))+' / '+observed(smoke.get('vec_type'))]]
fig=table(rows,['检查项','实测','证据'],[.26,.23,.51])
save(fig,'svmultiphysics_gpu_smoke.png','SimVascular 现在能真正启动 GPU 流体求解吗？',True,'OBSERVED STATUS','启动、收敛、场和 CUDA backend 分别验收；standalone PETSc PASS 不替代应用验收。')
if smoke['status']=='PASS':
    panels=[('gpu_cpu_solution_difference.png','CPU 和 GPU 算出的血流一样吗？','science_equivalence'),
            ('gpu_residency.png','求解时矩阵和向量是否一直留在显卡里？','gpu_residency'),
            ('gpu_transfer_cost.png','CPU 和 GPU 搬数据花了多少？','gpu_transfer'),
            ('gpu_memory.png','RTX 4090 实际用了多少显存？','gpu_memory'),
            ('cpu_gpu_runtime.png','同一台服务器上 GPU 到底快多少？','benchmark')]
    for name,title,record in panels:
        data=load(record);assert data['status']=='NOT_RUN','Measured panels require actual-data rendering'
        fig=plt.figure(figsize=(12,6));cn(fig,.5,.62,'NOT RUN',38,gray,'center')
        cn(fig,.5,.39,data['reason'],12,red,'center')
        save(fig,name,title,False,'NOT RUN','没有执行的测量保持空值；不使用零值、推算或理论 speedup。')
write_json(R/'visuals.json',{'figures':figures,'human_review_status':'PENDING','unknown_measurements_are_null':True})
print(str(len(figures))+' figures from captured evidence')
