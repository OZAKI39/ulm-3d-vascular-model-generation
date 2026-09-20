"""Five measured Stage M figures; downstream unexecuted stages stay NOT RUN."""
import json,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'));os.environ['MPLCONFIGDIR']=str(ROOT/'outputs/sv1_3m/plot_cache')
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from sv_validation.provenance import sha256
R=ROOT/'reports/sv1_3m';FONT=FontProperties(fname='/mnt/c/Windows/Fonts/msyh.ttc');plt.rcParams.update({'figure.facecolor':'#f8fafc','font.size':12,'axes.unicode_minus':False})
load=lambda n:json.loads((R/(n+'.json')).read_text());figures=[]
def table(name,title,columns,rows,caption,source,widths=None):
 fig=plt.figure(figsize=(13,7));ax=fig.add_axes([.045,.23,.91,.57]);ax.axis('off')
 tab=ax.table(cellText=rows,colLabels=columns,colWidths=widths,loc='center',cellLoc='left');tab.auto_set_font_size(False);tab.set_fontsize(11);tab.scale(1,2.5)
 for (i,j),cell in tab.get_celld().items():
  cell.get_text().set_fontproperties(FONT);cell.set_edgecolor('#dce3eb')
  if i==0:cell.set_facecolor('#dfe9f3')
  if 'FAIL' in cell.get_text().get_text():cell.get_text().set_color('#b23c40')
  elif 'PASS' in cell.get_text().get_text():cell.get_text().set_color('#137c66')
 fig.text(.045,.92,title,fontproperties=FONT,fontsize=21,color='#243244');fig.text(.045,.11,caption,fontproperties=FONT,fontsize=11,color='#475569')
 p=R/name;fig.savefig(p,dpi=150);plt.close(fig);figures.append(dict(path=str(p.relative_to(ROOT)),sha256=sha256(p),title=title,data_source=source,caption=caption))
rows=[]
for variant in ('CPU_seq','CUDA_seq','CPU_MPI','CUDA_MPI'):
 vals=[]
 for stage in ('before','repair_01','repair_02','repair_03'):
  runs=[r for r in load('ghost_probe_'+stage)['runs'] if r['variant']==variant];vals.append(('PASS' if all(r['accepted'] for r in runs) else 'FAIL')+f" ({sum(r['accepted'] for r in runs)}/{len(runs)})")
 rows.append([variant,*vals])
table('ghost_failure_before_after.png','CUDA ghost vector 的报错修好了吗？',['测试','原版 3.19.6','repair_01','repair_02','repair_03'],rows,'最终 ghost 回归全部通过；repair_02 虽已消除接口错误，反向数据仍错误，因此保留 FAIL。','ghost_probe_before/repair_01/repair_02/repair_03.json')
rows=[]
for variant in ('CPU_seq','CUDA_seq','CPU_MPI','CUDA_MPI'):
 runs=[r for r in load('ghost_probe_after')['runs'] if r['variant']==variant];checks=[c for r in runs for c in r['checks']];rows.append([variant,'3/3 PASS',str(max(float(c[1]) for c in checks)),str(max(float(c[2]) for c in checks)),f"{checks[0][3]} / {checks[0][4]}"])
table('ghost_data_correctness.png','Ghost 数据同步以后真的正确吗？',['测试','重复结果','正向最大误差','反向最大误差','实际 / 预期 local size'],rows,'单进程没有远端 ghost；双进程每 rank 有 2 个 ghost 标量。另有 GPU 运算、写回和复制检查全部通过。','ghost_probe_after.json; ghost_coherence.json')
rows=[['commonmpvec.c','4','类型识别 + local form 主机数据同步'],['vecreg.c','1','MPI → CUDA 保留原有数据结构'],['vecmpicupm.cu','1','复用原有 CUDA 初始化接口'],['vecimpl.h','声明','内部转换函数声明']]
table('patch_scope.png','这次到底改了 PETSc 的哪些地方？',['文件','函数数','作用'],rows,'累计 4 个文件、6 个函数、+54 / −5 行。两个上游来源 + 明确标注的本地同步适配；solver source diff = 0。','repair_03.json; patch_integrity.json',[.28,.13,.59])
rows=[['旧错误复现','PASS','实际 seqcuda + ghost error'],['CPU / GPU ghost','PASS','4 类 × 3 次 + 数据一致性'],['PETSc GPU / svMP build','PASS','3 次稀疏求解；正确动态链接'],['official GPU smoke','FAIL','4 次线性收敛、1 步完成；退出清理失败'],['真实血管 / 科学对比','NOT RUN','official smoke 门槛未通过'],['驻留 / 传输 / 性能','NOT RUN','不填估算值']]
table('gpu_stack_pipeline.png','RTX 4090 求解链路现在走到哪一步？',['步骤','状态','实际证据'],rows,'Stage M: FAIL — SVMULTIPHYSICS_GPU_SMOKE_FAIL。CPU production 保持不变。','stage_result.json; svmp_gpu_smoke.json',[.32,.17,.51])
s=load('svmp_gpu_smoke');rows=[['Mat / Vec','seqaijcusparse / seqcuda','PASS'],['KSP / PC','GMRES / ASM overlap 2 / ILU(2)','PASS'],['线性 / 非线性','4 次收敛 / 1 个时间步完成','PASS'],['VTU / finite / reload','1 个 VTU / 全部有限 / fresh reload','PASS'],['进程正常退出','MPI_Finalize 后仍访问 MPI；exit=1','FAIL']]
table('svmultiphysics_gpu_smoke.png','SimVascular 现在能完整跑完 GPU 流体测试吗？',['检查项','实际结果','状态'],rows,'生成正确格式的场不能抵消退出错误。独立小型 KSP 也复现：缺失清理失败，正确清理后 exit=0。','svmp_gpu_smoke.json; finalize_ksp_probe.json',[.22,.63,.15])
(R/'visuals.json').write_text(json.dumps(dict(status='PASS',figures=figures,downstream_figures='NOT_RUN: official smoke did not pass'),indent=2,ensure_ascii=False)+'\n')
