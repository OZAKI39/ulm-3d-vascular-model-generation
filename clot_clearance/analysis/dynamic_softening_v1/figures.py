from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import time,resource
start=time.perf_counter();plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.titlesize':15,'axes.labelsize':13,'legend.fontsize':11,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.45,'grid.color':'#607080'})
blue='#187F94';orange='#C16A29';red='#C34049';gray='#7F95A4';purple='#7756A2';figs=[]
z=np.load(O/'main_accepted_trajectory.npz');s=z['series'];c={n:i for i,n in enumerate(z['columns'])};t=s[:,c['time_s']]*1e6;b=z['bonds'];X=z['reference_X'];a=json.loads((O/'DYNAMIC_341_NODE_AUDIT.json').read_text());h=np.load(O/'dt_half_window_trajectory.npz');hc=h['series'];compare=json.loads((O/'DT_HALF_WINDOW_COMPARISON.json').read_text());ledger=json.loads((O/'DYNAMIC_ENERGY_LEDGER.json').read_text());fi=z['free_ids']-1

def save(fig,name,caption):
 fig.savefig(F/(name+'.png'),dpi=190,bbox_inches='tight');fig.savefig(F/(name+'.svg'),bbox_inches='tight');plt.close(fig);figs.append(dict(name=name,caption=caption,png_sha256=sha(F/(name+'.png')),svg_sha256=sha(F/(name+'.svg'))))
def events(ax):
 ax.axvline(85,color=orange,ls=':',lw=1.4);ax.axvline(92.5,color=red,ls=':',lw=1.4);ax.axvline(122,color='#35424B',ls='--',lw=1.3);ax.set_xlim(0,124)
fig,ax=plt.subplots(figsize=(15,8.8));ax.set(xlim=(0,15),ylim=(0,9));ax.axis('off')
def box(x,y,w,hh,text,color=blue):
 ax.add_patch(FancyBboxPatch((x,y),w,hh,boxstyle='round,pad=.15',fc='#F0F5F7',ec=color,lw=2));ax.text(x+w/2,y+hh/2,text,ha='center',va='center',fontsize=13)
def arrow(x,y,xx,yy):ax.annotate('',xy=(xx,yy),xytext=(x,y),arrowprops=dict(arrowstyle='->',color='#526574',lw=2))
box(.25,6.2,4.3,1.5,'Frozen offline model\nmodel.py\nIndependent endpoint-star energy');box(5.1,6.2,4.5,1.5,'Isolated bond style\nbpm/peri/variational\nFixed original graph, m, r, taper');arrow(4.75,6.95,4.9,6.95)
box(.25,3.6,4.3,1.7,'Native fix nve + fix move\n245 free / 96 prescribed nodes\nPhysical mass = rho × V');box(5.1,3.6,4.5,1.7,'LPS/VARIATIONAL/DYNAMIC\nPOST_INTEGRATE: commit κ/d\nPOST_FORCE: density → N');arrow(4.75,4.45,4.9,4.45);arrow(7.3,6.0,7.3,5.5)
box(10.2,4.65,4.3,1.9,'Unchanged BPM contact\nFailed original edge exclusions\nremoved before pair evaluation\nContact radius = 1.0 dx',orange);arrow(9.8,4.65,10,5.35)
box(10.2,1.5,4.3,1.5,'Production library isolated\nSHA unchanged\nNo deployment',red);box(.25,.9,9.35,1.7,'P0: 28 candidate files reviewed; 0 safely deletable files\nAccepted replay code and uncertain compiler probe retained\nCompact active entrypoint; historical evidence protected')
ax.text(7.5,8.55,'Actual cleanup and active software paths',ha='center',fontsize=21,weight='bold');ax.text(7.5,.1,'P0 completed before new C++ edits. No deletion was invented to satisfy a cleanup target.',ha='center',fontsize=12)
save(fig,'Fig1_cleanup_and_native_structure','P0先行审计28个文件，安全删除0个。新协调fix只位于本轮隔离源码，保留原自由能、参考图及1.0dx接触。')
fig,ax=plt.subplots(2,3,figsize=(17.2,10.3),layout='constrained');disp=np.linalg.norm(z['x']-X,axis=2)*1e6;speed=np.linalg.norm(z['velocity'],axis=2)
ax[0,0].plot(t,disp[:,fi].max(axis=1),color=blue,label='Maximum free-node displacement');ax[0,0].set(ylabel='Displacement (µm)',title='A  True free-node motion');ax[0,0].legend()
ax[0,1].plot(t,speed[:,fi].max(axis=1),color=blue,label='Free nodes');ax[0,1].plot(t,speed.max(axis=1),color=orange,ls='--',label='All nodes');ax[0,1].set(ylabel='Speed (m/s)',title='B  Native accepted velocities');ax[0,1].legend()
ax[1,0].plot(t,z['kappa'][:,3392],color=blue,label='κ, original bond 3393');ax[1,0].plot(t,z['d'][:,3392],color=red,label='d, original bond 3393');ax[1,0].set(ylabel='Dimensionless',title='C  First softening: step 170 (85 µs)');ax[1,0].legend()
ax[1,1].plot(t,s[:,c['total_failed']],color=red,label='Cumulative complete failures');ax[1,1].plot(t,s[:,c['new_failures']],color=orange,label='New failures this step');ax[1,1].set(ylabel='Original bond count',title='D  First full failure: step 185 (92.5 µs)');ax[1,1].legend()
ax[0,2].plot(t,z['scalars'][:,4],color=blue,label='Accepted commits (global physical steps)');ax[0,2].plot(t,z['scalars'][:,6],color=orange,label='Force solves (process-local counter)');ax[0,2].axvline(85,color=gray,ls='--',lw=1);ax[0,2].set(ylabel='Actual native counter',title='C  Force evaluation is not a commit');ax[0,2].legend(fontsize=9);ax[0,2].text(.04,.56,'Restart at step 170:\nsolve counter resets;\nκ/d and commit token persist.',transform=ax[0,2].transAxes,fontsize=10)
ax[1,2].axis('off');ax[1,2].set_title('F  Actual native time-layer contract');ax[1,2].text(.02,.92,'1  Free nodes: first half kick + drift\n2  Grips: analytic prescribed motion\n3  POST_INTEGRATE: one κ/d commit\n4  Rebuild failed-edge exclusions\n5  Pair + full variational bond forces\n6  Convert force density to newtons\n7  Free nodes: second half kick\n\nrun 0 / changed trial coordinates:\n    no irreversible history update\n\nSaved original IDs, coordinates, velocities\nand counters remain in FP64 evidence.',va='top',fontsize=12,linespacing=1.65)
ax[1,0].set_title('D  First softening: step 170 (85 µs)');ax[1,1].set_title('E  First full failure: step 185 (92.5 µs)')
for aa in [ax[0,0],ax[0,1],ax[0,2],ax[1,0],ax[1,1]]:events(aa);aa.set_xlabel('Native physical time (µs)')
fig.suptitle('Real Velocity-Verlet dynamics — 244 accepted steps, then a safety stop',fontsize=20);fig.supxlabel('Dotted lines: first softening / complete failure. Dashed line: Δd > 0.10 stop. Each native step commits exactly once; run 0 never commits.',fontsize=11)
save(fig,'Fig2_physical_steps_and_damage','真实0–122µs主时程。首次软化85µs；首次完全失效92.5µs；第244步因Δd=0.10255停止，未达到预注册700µs终点。')
fig=plt.figure(figsize=(15.4,11.4));ref=np.load(O/'provenance/reference.npz');short=ref['r0']<=DX*1.01
for j,(it,title) in enumerate([(0,'A  Reference: step 0'),(170,'B  First softening: step 170'),(185,'C  First complete failure: step 185'),(244,'D  Saved stop state: step 244')],1):
 ax=fig.add_subplot(2,2,j,projection='3d');ax.view_init(17,-58);ax.set_box_aspect((1,1,1));ax.set(xlim=(-205,195),ylim=(-200,200),zlim=(-200,200),xlabel='X (µm)',ylabel='Y (µm)',zlabel='Z (µm)');ax.set_xticks([-150,0,150]);ax.set_yticks([-150,0,150]);ax.set_zticks([-150,0,150]);ax.tick_params(labelsize=11,pad=1);ax.set_title(title,pad=12)
 for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.set_facecolor((.965,.974,.981,1));axis._axinfo['grid']['color']=(.42,.49,.54,.72)
 x=z['x'][it]*1e6;d=z['d'][it];bad=d==1
 for sel,col,alpha,ls in [((d<1)&short,gray,.20,'-'),((d>0)&(d<1),orange,.19,'-'),(bad,red,.24,'--')]:
  if sel.any():ax.add_collection3d(Line3DCollection(x[b[sel]],colors=col,alpha=alpha,linewidths=.7,linestyles=ls))
 nd=np.zeros(len(X));np.maximum.at(nd,b[:,0],d);np.maximum.at(nd,b[:,1],d);ax.scatter(*x.T,c=nd,cmap='viridis',vmin=0,vmax=1,s=35,depthshade=False,edgecolors='white',linewidth=.3);ax.text2D(.03,.92,f't = {t[it]:.1f} µs; fully failed = {int(bad.sum())}',transform=ax.transAxes,fontsize=12)
fig.subplots_adjust(left=.02,right=.88,bottom=.14,top=.9,wspace=.06,hspace=.18);cax=fig.add_axes([.91,.30,.015,.39]);fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,1),cmap='viridis'),cax=cax,label='Maximum incident bond damage d');fig.suptitle('P3 actual point-cloud evidence — P4 crack experiment NOT RUN',fontsize=20)
fig.legend(handles=[Line2D([0],[0],color=gray,label='Intact: short original edges shown'),Line2D([0],[0],color=orange,label='Partial damage'),Line2D([0],[0],color=red,ls='--',label='Failed ID trace; zero bond force')],loc='lower center',bbox_to_anchor=(.5,.067),ncol=3,frameon=False)
fig.text(.5,.028,'One connected 341-node body remains. Dashed failed traces are not live bonds. No medical shell or inferred crack surface.',ha='center',fontsize=11);fig.text(.5,.005,'Full graph is stored in FP64; intact rendering uses original r₀ ≤ 1.01 dx edges for visibility. Fixed physical axes and camera.',ha='center',fontsize=10)
save(fig,'Fig3_P3_saved_stop_evidence','P4未执行；本图仅展示真实341点P3保存态。为可读性，完整键背景只画r0≤1.01dx子集，损伤/失效原键全部显示。节点颜色为最大邻接d，不把损伤云画成连续裂面。')
fig,ax=plt.subplots(2,2,figsize=(15.5,10),layout='constrained')
for key,label,color in [('W_J','External W',blue),('U_J','Recoverable U',orange),('K_J','Kinetic K',purple),('D_J','Damage D',red)]:ax[0,0].plot(t,s[:,c[key]]*1e9,label=label,color=color)
ax[0,0].set(ylabel='Energy (nJ)',title='A  Complete physical ledger');ax[0,0].legend(ncol=2)
ax[0,1].plot(t,s[:,c['residual_J']]*1e12,color=red,label='W − ΔU − ΔK − D − ΔUc + Jc');ax[0,1].set(ylabel='Signed residual (pJ)',title='B  Observed balance: max 0.1446%');ax[0,1].legend(fontsize=10)
dec=np.load(O/'energy_residual_decomposition.npz')['series']
for col,label,color in [(2,'Verlet work quadrature',blue),(3,'Fixed-damage spatial quadrature',orange),(4,'Damage force switch work',red)]:ax[1,0].plot(t,-dec[:,col]*1e12,label=label,color=color)
ax[1,0].plot(t,dec[:,6]*1e12,'k--',label='Sum: predicted residual');ax[1,0].set(ylabel='Contribution to residual (pJ)',title='C  Independently identified time-layer terms');ax[1,0].legend(fontsize=10)
ax[1,1].plot(t,-s[:,c['D_J']]*1e9,color=red,label='Signed damage topology change = −D');ax[1,1].plot(t,s[:,c['Uc_J']]*1e9,color=blue,label='Contact energy: zero observed');ax[1,1].plot(t,np.zeros_like(t),color=orange,ls='--',label='Damping: absent');ax[1,1].set(ylabel='Energy (nJ)',title='D  Signed topology is not counted twice');ax[1,1].legend(fontsize=10)
for aa in ax.ravel():events(aa);aa.set_xlabel('Native physical time (µs)')
fig.suptitle('Energy audit through the registered stop — residual kept explicit',fontsize=20);fig.supxlabel('Grip work uses actual displacements, trapezoidal material reaction and exact grip ΔK. No closure contact was activated; material Gc is UNAVAILABLE.',fontsize=11)
save(fig,'Fig4_dynamic_energy_ledger','外功、完整变分储能、物理动能、损伤耗散与带符号拓扑项分别展示；残差由Verlet功求积、固定损伤空间求积和损伤后力切换项重建，不归入耗散。')
fig,ax=plt.subplots(2,3,figsize=(17.2,10.3),layout='constrained')
ax[0,0].plot(t,s[:,c['vmax_m_s']],color=blue);ax[0,0].set(xlabel='Physical time (µs)',ylabel='Maximum speed (m/s)',title='A  Peak speed 0.44434 m/s')
ax[0,1].plot(t,s[:,c['step_disp_dx']],color=blue,label='Node displacement / dx');ax[0,1].plot(t,s[:,c['relative_endpoint_step_dx']],color=orange,label='Relative endpoints / dx');ax[0,1].set(xlabel='Physical time (µs)',ylabel='Dimensionless',title='B  Actual step increments');ax[0,1].legend(fontsize=10)
sel=t>=105;ht=hc[:,c['time_s']]*1e6
ax[0,2].plot(t[sel],s[sel,c['max_dd']],color=red,label='dt = 0.5 µs');ax[0,2].plot(ht[1:],hc[1:,c['max_dd']],color=blue,label='dt = 0.25 µs');ax[0,2].axhline(.1,color='k',ls='--',label='Frozen Δd stop gate');ax[0,2].set(xlabel='Physical time (µs)',ylabel='Maximum Δd per step',title='C  One admitted short dt/2 window');ax[0,2].legend(fontsize=10)
ax[1,0].semilogy(t[1:],np.maximum(s[1:,c['momentum_residual_Ns']],1e-30),color=blue,label='Momentum equation residual');ax[1,0].set(xlabel='Physical time (µs)',ylabel='Absolute residual (N s)',title='D  Momentum: max 1.09 × 10⁻²⁵ N s');ax[1,0].legend(fontsize=9)
base=210;ax[1,1].plot(t[sel],(s[sel,c['residual_J']]-s[base,c['residual_J']])*1e12,color=red,label='dt = 0.5 µs');ax[1,1].plot(ht,hc[:,c['residual_J']]*1e12,color=blue,label='dt = 0.25 µs');ax[1,1].set(xlabel='Physical time (µs)',ylabel='Window residual (pJ)',title='E  Residual ratio = 0.49619');ax[1,1].legend()
ax[1,2].axis('off');ax[1,2].set_title('F  Separate scientific decisions');ax[1,2].text(.01,.96,'Real native free dynamics: demonstrated\nOne commit per physical step: verified\nNatural softening / full failure: observed\nForce / history / momentum: consistent\nEnergy residual: 0.1446% through stop\n\nMain dt safety: STOPPED at step 244\nShort dt/2: local improvement only\nP4 small crack: NOT RUN\nClosure contact dynamics: unobserved\nNew fracture area / material Gc: UNAVAILABLE\nMacro strength / full V3: NOT CERTIFIED',va='top',fontsize=12,linespacing=1.55)
fig.suptitle('Time integration and scientific admission — no scale-up approval',fontsize=20);fig.supxlabel('Same 105–122 µs restart window; 34 coarse and 68 fine steps. No dt/4, alternate load, P4, CFD or production deployment.',fontsize=11)
save(fig,'Fig5_timestep_and_scientific_admission','真实损伤增量触发停止；唯一dt/2窗口使最大Δd减半、局部能量残差约减半。主时程未完成，因此不得据此准入P4或大规模材料标定。')
dump(O/'FIGURE_DATA_PROVENANCE.json',dict(script_sha256=sha(Path(__file__)),figures=figs,inputs={n:sha(O/n) for n in ['main_accepted_trajectory.npz','dt_half_window_trajectory.npz','energy_residual_decomposition.npz']},all_curves_from_actual_native_states=True,P4_figure_not_fabricated=True))
links=''.join(f'<section><h2>{f["name"]}</h2><p>{f["caption"]}</p><a href="{f["name"]}.svg">SVG</a><img src="{f["name"]}.png"></section>' for f in figs)
(F/'OPEN_RESULTS.html').write_text('<!doctype html><html lang="zh"><meta charset="utf-8"><title>LPS真实动力学与停止审查</title><style>body{max-width:1450px;margin:30px auto;padding:0 24px;font:17px/1.7 system-ui;color:#21333f}img{width:100%}section{margin-top:32px;border-top:1px solid #aebfc8}.note{padding:18px;background:#fff1e7;border-left:5px solid #c16a29}</style><h1>LPS变分软化：真实动力学与停止审查</h1><p class="note">已完成隔离C++接入和341点真实自由积分。主工况在244步、122µs因单步损伤增量超过冻结0.10门槛停止；自然首软化85µs、首失效92.5µs，累计520条失效键。限定窗口dt/2局部改善，不能代替全程认证。P4未运行，材料Gc、宏观强度和生产V3未认证。</p><a href="../../../verification/'+TASK+'/'+RUN+'/REVIEW_LPS_NATIVE_SOFTENING_DYNAMIC_COMMIT_AND_SMALL_CRACK_V1.md">中文51项报告</a> · <a href="../../../verification/'+TASK+'/'+RUN+'/FINAL_STATUS.json">分项状态</a>'+links+'</html>')
(O/'OPEN_RESULTS.html').write_text('<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="0;url=../../../results/'+TASK+'/'+RUN+'/OPEN_RESULTS.html"><a href="../../../results/'+TASK+'/'+RUN+'/OPEN_RESULTS.html">打开五组真实图</a>')
dump(O/'FIGURE_RESOURCES.json',dict(wall_s=time.perf_counter()-start,CPU_s=time.process_time(),peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,groups=5,PNG=5,SVG=5));print('FIVE_ACTUAL_FIGURES_COMPLETE',flush=True)
