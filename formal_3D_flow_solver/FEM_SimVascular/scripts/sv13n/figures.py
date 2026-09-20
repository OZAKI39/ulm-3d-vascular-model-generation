"""Five development figures from actual GPU evidence; no performance comparison."""
import json,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
os.environ['MPLCONFIGDIR']=str(ROOT/'outputs/sv1_3n/plot_cache');os.environ['LIBGL_ALWAYS_SOFTWARE']='1'
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
from sv_validation.provenance import sha256
R=ROOT/'reports/sv1_3n';FONT=FontProperties(fname='/mnt/c/Windows/Fonts/msyh.ttc');figures=[]
def read(n):
 p=R/(n+'.json');return json.loads(p.read_text()) if p.exists() else dict(status='NOT RUN')
def save(fig,name,title,note,source):
 fig.suptitle(title,fontproperties=FONT,fontsize=20,y=.98)
 fig.text(.5,.025,note,ha='center',fontproperties=FONT,fontsize=10,color='#334155')
 fig.savefig(R/name,dpi=150,facecolor='white');plt.close(fig)
 figures.append(dict(path='reports/sv1_3n/'+name,sha256=sha256(R/name),title=title,data_source=source))
def table(name,title,rows,note,source):
 fig,ax=plt.subplots(figsize=(14,8));ax.axis('off')
 t=ax.table(cellText=rows,colLabels=['检查项目','实际结果'],colWidths=[.38,.62],loc='center',cellLoc='left');t.auto_set_font_size(False);t.set_fontsize(12);t.scale(1,2.1)
 for (i,j),cell in t.get_celld().items():
  cell.get_text().set_fontproperties(FONT);cell.set_edgecolor('#dbe2eb')
  if i==0:cell.set_facecolor('#dfe9f3')
  if 'FAIL' in cell.get_text().get_text():cell.get_text().set_color('#b42334')
 save(fig,name,title,note,source)
candidate=read('gpu_steady_candidate');smoke=read('svmp_gpu_smoke')
table('gpu_stack_pipeline.png','SimVascular 的 RTX 4090 求解链路已经走到哪一步？',[
 ['PETSc / CUDA',read('petsc_gpu13_build')['status']+' — v3.25.5 / CUDA 13.2'],
 ['svMultiPhysics 构建 / 链接',read('svmp_gpu_build')['status']+' / '+read('svmp_gpu_link')['status']],
 ['单 rank ghost / 生命周期',read('gpu_ghost_target')['status']+' / '+read('solver_finalize_lifecycle')['status']],
 ['官方 GPU 流体测试',smoke['status']],['真实血管 GPU',read('REAL_VASCULAR_GPU_acceptance')['status']],
 ['稳态 / 安全停止 / reload',candidate['status']],['结果用途','GPU_STEADY_CANDIDATE；科学等价验证 DEFERRED']],
 '限定 1 MPI rank、1 RTX 4090；CPU production 保留。','petsc_gpu13_build; svmp_gpu_build; svmp_gpu_link; gpu_steady_candidate')
ghost=read('ghost_probe_new');rows=[]
for r in ghost.get('runs',[]):
 forward=max((float(x[1]) for x in r['checks']),default=float('nan'));reverse=max((float(x[2]) for x in r['checks']),default=float('nan'))
 rows.append([r['variant'],f"{'PASS' if r['accepted'] else 'FAIL'}; forward error={forward:g}, reverse error={reverse:g}"])
table('cuda_ghost_status.png','新版 PETSc 能正常处理 CUDA ghost vector 吗？',rows or [['检查','NOT RUN']],
 '单 rank 求解路径通过；额外 2-rank CUDA reverse/local-form 检查失败，未宣称 MPI CUDA 支持。','ghost_probe_new.json; gpu_ghost_target.json')
run=read('REAL_VASCULAR_GPU_acceptance');linear=run.get('history',{}).get('linear_solves',[])
fig,axes=plt.subplots(1,2,figsize=(14,8));fig.subplots_adjust(top=.85,bottom=.17,wspace=.28)
if linear:
 for nl in sorted({r['nonlinear_iteration'] for r in linear}):
  rr=[r for r in linear if r['nonlinear_iteration']==nl]
  axes[0].plot([r['step'] for r in rr],[r['linear_iterations'] for r in rr],'.-',ms=3,label=f'NL {nl}')
 final={r['step']:r for r in linear};axes[1].semilogy(list(final),[max(min(r['nonlinear_Ri_over_R0'],r['nonlinear_Ri_over_R1']),1e-30) for r in final.values()],'.-')
 axes[1].axhline(1e-10,color='#b42334',ls='--',label='Frozen nonlinear tolerance')
 axes[0].legend(fontsize=9);axes[1].legend(fontsize=9)
else:
 for a in axes:a.text(.5,.5,'NOT RUN',ha='center',transform=a.transAxes)
for a,y in zip(axes,['KSP iterations / solve','Final nonlinear residual ratio']):a.set_xlabel('Timestep');a.set_ylabel(y);a.grid(alpha=.2)
save(fig,'gpu_solver_progress.png','真实血管 GPU 求解是否持续稳定运行？',f"线性失败 {run.get('linear_failures','NOT RUN')}；非线性失败 {run.get('nonlinear_failures','NOT RUN')}。KSP 收敛不等于流场稳态。",'REAL_VASCULAR_GPU_acceptance.json')
h=read('gpu_steady_history');p=json.loads((ROOT/'configs/sv1_3/policy.json').read_text());fig,axes=plt.subplots(1,3,figsize=(16,8));fig.subplots_adjust(top=.83,bottom=.17,wspace=.32)
for a,key,limit in zip(axes,['E_u','E_Q','epsilon_mass'],['velocity_change_limit','flow_change_limit','mass_limit']):
 rows=h.get('states',[]) if key=='epsilon_mass' else h.get('intervals',[])
 if rows:a.plot([r['step'] for r in rows],[r[key] for r in rows],'o-')
 else:a.text(.5,.5,'NOT RUN',ha='center',transform=a.transAxes)
 a.set_yscale('symlog',linthresh=1e-15);a.axhline(p[limit],color='#b42334',ls='--');a.set_xlabel('Formal monitor timestep');a.set_ylabel(key);a.grid(alpha=.2)
 if candidate['status']=='PASS':a.axvline(candidate['first_full_steady_step'],color='#137c66',ls=':')
save(fig,'gpu_steady_convergence.png','GPU 流场什么时候达到稳态？',f"production 原门槛；每 10 步检查，连续 5 个合格区间。首次完整通过：{candidate.get('first_full_steady_step','NOT RUN')}；安全停止：{candidate.get('stop_step','NOT RUN')}。",'gpu_steady_history.json; gpu_steady_candidate.json')
fig,ax=plt.subplots(figsize=(15,8));ax.axis('off');fig.subplots_adjust(top=.88,bottom=.12,left=.02,right=.98)
if candidate['status']=='PASS':
 import pyvista as pv
 pv.OFF_SCREEN=True;grid=pv.read(ROOT/candidate['VTU']);speed=np.linalg.norm(np.asarray(grid['Velocity']),axis=1)
 scene=pv.Plotter(shape=(1,2),off_screen=True,window_size=(1800,760));scene.set_background('white',all_renderers=True)
 scene.subplot(0,0);surface=grid.extract_surface(algorithm='dataset_surface');scene.add_mesh(surface,color='#9fb4c7',opacity=.06)
 cloud=pv.PolyData(grid.points.copy());cloud['Speed']=speed
 scene.add_mesh(cloud,scalars='Speed',cmap='viridis',point_size=3,opacity=.12,scalar_bar_args={'title':'Speed (m/s)','fmt':'%.2e','n_labels':4,'width':.60,'position_x':.20,'title_font_size':18,'label_font_size':16});scene.view_isometric();scene.reset_camera();scene.add_axes()
 scene.subplot(0,1);scene.add_mesh(surface,scalars='Pressure',cmap='coolwarm',scalar_bar_args={'title':'Native pressure (Pa)','fmt':'%.2e','n_labels':4,'width':.60,'position_x':.20,'title_font_size':18,'label_font_size':16});scene.view_isometric();scene.reset_camera();scene.add_axes()
 frame=scene.screenshot(return_img=True);scene.close();ax.imshow(frame)
else:ax.text(.5,.5,'NOT RUN — no complete GPU steady candidate',ha='center')
save(fig,'gpu_final_field.png','GPU 收敛后的速度和压力场是否完整？','原生节点速度模长（透明显示全部节点）与外表面原生压力；未平移压力。GPU candidate 尚未通过 CPU/GPU 科学等价验证。','gpu_steady_candidate.json; native VTU')
(R/'visuals.json').write_text(json.dumps(dict(status='PASS',figures=figures),indent=2,ensure_ascii=False)+'\n')
print('Five development figures written; no performance statistics plotted.')
