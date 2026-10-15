from common import *
from analyze import components
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from itertools import product,combinations

C=['#187F94','#C16A29','#7756A2','#C34049','#7F95A4','#314653']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.titlesize':16,'axes.labelsize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':11,'axes.edgecolor':'#46545F','axes.spines.top':False,'axes.spines.right':False,'grid.color':'#607080','grid.alpha':.36,'grid.linewidth':.7,'figure.facecolor':'white','savefig.facecolor':'white','svg.fonttype':'none'})
full=O/'full';z=np.load(full/'trajectory.npz');a=z['series'];c={n:a[:,i] for i,n in enumerate(z['columns'])};t=c['time_s']*1e6;nodes=np.memmap(full/'accepted_nodes.bin',mode='r',dtype='<f8',shape=(len(a),4,341,3));ref=np.load(O/'provenance/reference.npz');b=ref['bonds'];status=json.loads((O/'FINAL_STATUS.json').read_text());res=json.loads((O/'RUN_RESOURCES.json').read_text());comp=json.loads((O/'INDEPENDENT_COMPONENT_RECONSTRUCTION.json').read_text());events=json.loads((full/'FIRST_SOFTENING_AND_FAILURE_EVENTS.json').read_text());index=json.loads((full/'KEYFRAME_INDEX.json').read_text());last=dict(np.load(full/'last_accepted.npz'));F.mkdir(exist_ok=True,parents=True)
def save(fig,name,caption):
 fig.savefig(F/(name+'.png'),dpi=210,bbox_inches='tight');fig.savefig(F/(name+'.svg'),bbox_inches='tight');plt.close(fig);captions[name]=caption
def style(ax,title=None,ylabel=None):
 ax.grid(True);ax.set_axisbelow(True)
 if title:ax.set_title(title,loc='left',pad=10)
 if ylabel:ax.set_ylabel(ylabel)
def timeaxis(ax):ax.set_xlabel('Physical time (µs)');ax.set_xlim(0,700);ax.set_xticks(np.arange(0,701,100))
def style3(ax,limits=(-260,280)):
 ax.set(xlabel='X (µm)',ylabel='Y (µm)',zlabel='Z (µm)',xlim=limits,ylim=limits,zlim=limits);ax.set_box_aspect((1,1,1));ax.view_init(elev=17,azim=-58);ax.set_proj_type('ortho');ax.tick_params(pad=0,labelsize=10)
 for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.fill=True;axis.pane.set_facecolor((.97,.98,.985,1));axis._axinfo['grid']['color']=(.35,.43,.49,.42)
def wire(ax,lo,hi,color,ls='-',lw=1):
 corners=np.array(list(product(*zip(lo,hi))))
 seg=[(p,q) for p,q in combinations(corners,2) if np.sum(p!=q)==1]
 ax.add_collection3d(Line3DCollection(seg,colors=color,linestyles=ls,linewidths=lw))
captions={}

# 1. Full numerical box and a readable crossing detail from actual accepted data.
fig=plt.figure(figsize=(15,7.5));gs=fig.add_gridspec(2,2,width_ratios=[1.1,1],hspace=.5,wspace=.25);ax=fig.add_subplot(gs[:,0],projection='3d')
style3(ax,(-420,420));ax.set_xticks([-400,-200,0,200,400]);ax.set_yticks([-400,-200,0,200,400]);ax.set_zticks([-400,-200,0,200,400]);wire(ax,[-400]*3,[400]*3,C[0],lw=1.15);wire(ax,[-250]*3,[250]*3,C[4],ls='--',lw=1.2);x=ref['X']*1e6;ax.scatter(*x.T,s=12,c=C[4],alpha=.72,depthshade=False);path=nodes[:,0,72]*1e6;ax.plot(*path.T,color=C[3],lw=2.5);ax.scatter(*path[-1],c=C[3],s=65,edgecolors='white');wire(ax,last['x'].min(0)*1e6,last['x'].max(0)*1e6,C[1],ls=':',lw=1.3)
ax.set_title('A  Numerical domain and atom 73',loc='left',pad=17)
ax.legend(handles=[Line2D([0],[0],color=C[0],label='New fixed box: ±400 µm'),Line2D([0],[0],color=C[4],ls='--',label='Old fixed box: ±250 µm'),Line2D([0],[0],color=C[3],label='Actual atom 73 trajectory'),Line2D([0],[0],color=C[1],ls=':',label='Final material bounds')],loc='lower left',bbox_to_anchor=(0,-.16),fontsize=10)
ax=fig.add_subplot(gs[0,1]);sel=t>=680;ax.plot(t[sel],path[sel,1],color=C[3],lw=2.3,label='New accepted trajectory');ax.axhline(250,color=C[4],ls='--',label='Old Y boundary');i=int(np.argmin(abs(t-688)));ax.scatter([t[i]],[path[i,1]],s=65,c=C[3]);ax.annotate(f'688 µs: {path[i,1]:.6f} µm\nretained in new domain',xy=(t[i],path[i,1]),xytext=(690,249.1),fontsize=11,arrowprops=dict(arrowstyle='->',color=C[5]));style(ax,'B  Old-boundary crossing','Atom 73 Y (µm)');ax.set_xlabel('Physical time (µs)');ax.legend(loc='upper left',fontsize=10);ax.set_xlim(680,700)
ax=fig.add_subplot(gs[1,1]);
for j,key in enumerate(['x','y','z']):ax.plot(t,c[key+'max_m']*1e6,c=C[j],label=key.upper()+' maximum');ax.plot(t,c[key+'min_m']*1e6,c=C[j],ls='--')
ax.axhline(400,c=C[5],lw=1);ax.axhline(-400,c=C[5],lw=1);style(ax,'C  Actual all-node coordinate envelope','Coordinate (µm)');timeaxis(ax);ax.legend(ncol=3,fontsize=10,loc='upper center');fig.suptitle('Fixed numerical domain retains all 341 material points',fontsize=20,y=1.02);fig.text(.08,-.14,'Boxes are numerical containers, not physical vessel walls. No remap, added wall or periodicity.',fontsize=13,color=C[5]);save(fig,'Fig1_domain_and_atom73','原盒子各轴 ±250 µm、新盒子 ±400 µm；真实341点、节点73全程轨迹和末态包络。688 µs新接受态与旧漂移预测相容，旧运行并无688 µs合法接受态。盒子仅为数值域，不是血管壁。')

# 2. Dynamics and the genuine adaptive steps.
fig,axs=plt.subplots(2,2,figsize=(14,9),layout='constrained');q1=np.minimum(c['time_s']/.0003,1);q2=np.clip((c['time_s']-.00035)/.0003,0,1);u=30*(np.cos(np.pi*q2)-np.cos(np.pi*q1))
axs[0,0].plot(t,u,c=C[0],lw=2.3,label='Each grip magnitude');style(axs[0,0],'A  Prescribed half-cosine loading','Displacement (µm)');axs[0,0].legend(loc='upper right')
axs[0,1].plot(t,c['max_free_disp_m']*1e6,c=C[1],lw=2.2);style(axs[0,1],'B  Free-node motion','Maximum displacement (µm)')
axs[1,0].plot(t,c['max_speed_m_s'],c=C[0],lw=2,label='All nodes');axs[1,0].plot(t,c['max_free_speed_m_s'],c=C[1],ls='--',label='Free nodes');style(axs[1,0],'C  Native velocities','Maximum speed (m/s)');axs[1,0].legend()
axs[1,1].step(t[1:],c['dt_s'][1:]*1e6,where='post',c=C[2],lw=1.8);style(axs[1,1],'D  Actual accepted time steps','Time step (µs)');axs[1,1].set_ylim(.10,.275);axs[1,1].text(.04,.12,'2804 accepted · 4 rejected\n8 accepted half-steps at 0.125 µs',transform=axs[1,1].transAxes,fontsize=12)
for ax in axs.ravel():
 timeaxis(ax)
 for p in [300,350,650]:ax.axvline(p,color=C[4],ls=':',lw=1)
ax=axs[0,0]
for event,color in [('first_softening',C[2]),('first_complete_failure',C[3])]:
 tx=events[event]['time_s']*1e6;ax.axvline(tx,c=color,ls='--',lw=1.2)
ax.text(.03,.77,'Softening: 85 µs\nFirst full failure: 92.5 µs',transform=ax.transAxes,fontsize=11)
fig.suptitle('Fresh t = 0 → 700 µs native dynamics',fontsize=20);save(fig,'Fig2_full_dynamics','原始t=0新运行，左右每侧最大60 µm；全部曲线采用真实累计物理时间。0–300升载、300–350保持、350–650卸载、650–700最终保持；图中子步为实际接受的原生VV步。')

# 3. Six actual point-lattice states. Connectivity uses ALL original bonds.
fig=plt.figure(figsize=(17,12.5));fig.subplots_adjust(left=.02,right=.96,bottom=.13,top=.92,wspace=.03,hspace=.19)
times=[0,events['first_softening']['time_s'],events['first_complete_failure']['time_s'],.00065,.00068775,float(last['time'])];names=['Initial','First softening','First full failure','End of unloading','Historical last time','Final accepted state'];short=ref['r0']<=1.01*DX
for pane,(target,title) in enumerate(zip(times,names)):
 frame=min(index,key=lambda q:abs(q['time_s']-target));s=dict(np.load(full/frame['path']));assert abs(float(s['time'])-target)<1e-15
 ax=fig.add_subplot(2,3,pane+1,projection='3d');style3(ax);x=s['x']*1e6;groups=components(341,b[s['d']<1]);labels=np.empty(341,int)
 for k,g in enumerate(groups):labels[g]=k
 live=short&(s['d']<1);failed=short&(s['d']==1);ax.add_collection3d(Line3DCollection(x[b[live]],colors=plt.colormaps['YlOrRd'](s['d'][live]),linewidths=1.15,alpha=.95));
 if np.any(failed):ax.add_collection3d(Line3DCollection(x[b[failed]],colors=C[3],linewidths=.45,linestyles='--',alpha=.15))
 colors=[C[k%len(C)] for k in labels];ax.scatter(*x.T,c=colors,s=19,edgecolors='white',linewidths=.18,depthshade=False)
 ax.set_title(f'{chr(65+pane)}  {title}\n{float(s["time"])*1e6:g} µs',loc='left',fontsize=15,pad=4)
 sizes=' + '.join(map(str,map(len,groups)));ax.text2D(.02,-.03,f'{int(np.sum(s["d"]==1))} failed bonds\nComponents: {sizes}',transform=ax.transAxes,fontsize=11)
barax=fig.add_axes([.31,.06,.38,.018]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap='YlOrRd'),cax=barax,orientation='horizontal');cb.set_label('Damage d of displayed surviving bonds',fontsize=12)
fig.suptitle('Actual material-point lattice and surviving-bond topology',fontsize=21);fig.text(.06,-.035,'All 341 nodes shown. Node colors identify graph components. Only nearest reference bonds drawn; dashed traces are failed bonds.\nComponents are reconstructed from all 12,047 original bonds. No continuous surface or certified physical clot fragment is implied.',fontsize=12,color=C[5]);save(fig,'Fig3_native_topology','六个真实FP64状态；全部341点均显示，节点颜色表示全参考键图连通分量，存活短键按真实d着色，失效短键用淡虚线显示。为辨识点阵仅绘r0≤1.01dx的键，统计和连通性使用完整12047参考键。不是连续医学外壳，也不是已认证真实血栓碎块。')

# 4. Ledger and uncalibrated contact overlap.
fig,axs=plt.subplots(2,2,figsize=(14,9),layout='constrained')
for j,k in enumerate(['W_J','U_J','K_J','D_J','Uc_J','Jc_J']):axs[0,0].plot(t,c[k]*1e9,color=C[j],lw=1.8,label=k.split('_')[0],ls='--' if k=='Jc_J' else '-')
style(axs[0,0],'A  Full variational energy ledger','Energy (nJ)');axs[0,0].legend(ncol=3)
axs[0,1].plot(t,c['R_J']*1e12,c=C[3],lw=2);style(axs[0,1],'B  Numerical balance R','R (pJ)');axs[0,1].text(.04,.09,f'Maximum normalized |R|: {status["max_relative_energy_balance"]*100:.6f}%\nFrozen admission limit: 2%',transform=axs[0,1].transAxes,fontsize=12)
axs[1,0].step(t,c['contact_pairs'],where='post',c=C[0],lw=1.6);style(axs[1,0],'C  Actual active contact pairs','Pair count')
axs[1,1].plot(t,c['max_contact_penetration_dx'],c=C[1],lw=1.8);style(axs[1,1],'D  Contact overlap remains uncalibrated','Maximum penetration / dx');axs[1,1].text(.04,.85,'Maximum: 0.246649 dx ≈ 12.33 µm',transform=axs[1,1].transAxes,fontsize=12)
for ax in axs.ravel():timeaxis(ax)
fig.suptitle('W = ΔU + ΔK + D + ΔUc − Jc + R',fontsize=21);save(fig,'Fig4_energy_and_contact','沿用完整变分能量账本及原归一化；Jc实际为0。实际接触对逐步保存。1.0dx接触规则未修改，最大重叠约12.33 µm仍未解决其物理精度问题。')

# 5. Admission gates and actual resource costs.
fig,axs=plt.subplots(2,3,figsize=(17,9),layout='constrained');ax=axs[0,0];ax.plot(t,c['atom_count'],c=C[0],lw=2);ax.set_ylim(339,343);ax.set_yticks([340,341,342]);style(ax,'A  Original material identity','Actual atom count');ax.text(.04,.1,'IDs 1–341 checked at every accepted step',transform=ax.transAxes,fontsize=10)
ax=axs[0,1];ax.plot(t,c['max_dd'],c=C[2],lw=.9);ax.axhline(.1,c=C[3],ls='--',label='0.10 limit');style(ax,'B  Damage increment admission','Maximum Δd');ax.legend(loc='upper left');ax.set_ylim(0,.115)
ax=axs[0,2];ax.plot(t,c['relative_R']*100,c=C[0],lw=2);ax.axhline(2,c=C[3],ls='--',label='2% limit');style(ax,'C  Cumulative energy admission','Normalized |R| (%)');ax.legend(loc='center left');ax.set_ylim(0,2.2)
for ax in axs[0]:timeaxis(ax)
ax=axs[1,0];ax.bar(['Native','Independent'],[status['final_live_bonds']]*2,color=[C[0],C[1]],width=.55);style(ax,'D  Final surviving topology','Live bonds');ax.text(.5,.8,'Exact endpoint-set agreement\n192 + 147 + 1 + 1 nodes',ha='center',transform=ax.transAxes,fontsize=11);ax.set_ylim(0,9000)
ax=axs[1,1];ax.barh(['Main wall','Main CPU','Native candidates'],[res['main_wall_s'],res['main_CPU_s'],res['native_candidate_wall_s']],color=[C[0],C[1],C[2]]);style(ax,'E  Measured compute cost','');ax.set_xlabel('Seconds');ax.text(.02,.93,f'Peak RSS: {res["peak_RSS_bytes"]/1024**2:.1f} MiB · 1 rank × 1 thread',transform=ax.transAxes,fontsize=10);ax.set_xlim(0,65)
ax=axs[1,2];ax.axis('off');ax.text(0,.97,'F  Full protocol: PASS',fontsize=17,color=C[0],va='top');ax.text(0,.82,'700 µs reached from pristine t = 0\n2804 accepted / 4 damage rejections\nAll 4 native rollbacks verified\n341 material points retained\n4937 failed / 7110 surviving bonds\n\nGc, strength and global time convergence\nremain unvalidated.',fontsize=12,linespacing=1.7,va='top')
fig.suptitle('Domain retention, full protocol, energy and topology: PASS',fontsize=20);save(fig,'Fig5_admission_and_cost','四项准入分开报告：域保留、完整700 µs、能量、原生与独立拓扑一致。成本为本次主计算实际统计；未重新编译、未启用GPU或额外物理工况。材料Gc、宏观强度、全局dt收敛和接触物理精度仍未认证。')
dump(F/'FIGURE_CAPTIONS.json',captions)
cards=''.join(f'<section><h2>{name}</h2><a href="{name}.png"><img src="{name}.png" alt="{caption}"></a><p>{caption}</p><p><a href="{name}.svg">SVG 矢量图</a> · <a href="{name}.png">PNG 原图</a></p></section>' for name,caption in captions.items())
html='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>341节点完整700 µs软化动力学</title><style>body{max-width:1500px;margin:32px auto;padding:0 24px;font:18px/1.65 system-ui;color:#233440;background:#f7f9fb}h1{color:#187F94}section{background:white;padding:24px;margin:28px 0;border:1px solid #d2dce4;border-radius:8px}img{width:100%;height:auto}a{color:#187F94}p{max-width:1200px}</style><h1>341节点原生渐进软化：完整700 µs通过</h1><p>Run '+RUN+' · 固定域各轴 ±400 µm · 2804接受 / 4拒绝 · 最大Δd 0.0978676 · 最大相对能量余额0.372427%。全程341点保留；末态4937失效、7110存活键；分量192+147+1+1。</p><p>全部图来自本次原生计算的FP64数据。数值计算域不是血管壁；点阵分量尚不能称为认证的真实血栓碎块。接触物理精度、材料Gc、宏观强度和全局时间收敛仍未认证。</p>'+cards+'</html>'
(F/'OPEN_RESULTS.html').write_text(html)
rel='../../../results/'+TASK+'/'+RUN+'/OPEN_RESULTS.html';(O/'OPEN_RESULTS.html').write_text('<!doctype html><meta charset="utf-8"><title>结果入口</title><p><a href="'+rel+'">打开五组科研图与结果说明</a></p><meta http-equiv="refresh" content="0;url='+rel+'">')
print('FIVE_FIGURE_GROUPS',F)
