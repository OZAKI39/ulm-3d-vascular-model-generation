"""Six honest figure groups from existing accepted FP64 data only."""
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection,Poly3DCollection

C=['#187F94','#C16A29','#7756A2','#C34049','#7F95A4','#314653']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.titlesize':16,'axes.labelsize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':11,'axes.edgecolor':'#46545F','axes.spines.top':False,'axes.spines.right':False,'grid.color':'#607080','grid.alpha':.4,'figure.facecolor':'white','savefig.facecolor':'white','svg.fonttype':'none'})
ref=dict(np.load(O/'provenance/reference.npz'));s=dict(np.load(O/'full/last_accepted.npz'));h=np.genfromtxt(O/'NATIVE_DYNAMIC_HISTORY.csv',delimiter=',',names=True);t=h['time_s']*1e6
fi=np.genfromtxt(O/'full/native_frames.csv',delimiter=',',names=True);nodes=np.memmap(O/'full/native_frames.bin',mode='r',dtype='<f8',shape=(len(fi),4,len(ref['X']),3))
X=ref['X'];b=ref['bonds'];cut=ref['initial_cut'];near=ref['r0']<=DX*1.01;grip=~ref['free_dof'][:,2]
events=json.loads((O/'FIRST_SOFTENING_AND_FAILURE.json').read_text());loc=json.loads((O/'CRACK_LOCALIZATION_AUDIT.json').read_text());plane=dict(np.load(O/'PLANE_CELL_AUDIT.npz'));energy=json.loads((O/'DYNAMIC_ENERGY_LEDGER.json').read_text());resource=json.loads((O/'RUN_RESOURCES.json').read_text());comp=json.loads((O/'NATIVE_COMPONENT_TOPOLOGY.json').read_text())
F.mkdir(exist_ok=True,parents=True);captions={};provenance={}
def save(fig,name,caption):
 fig.savefig(F/(name+'.png'),dpi=230,bbox_inches='tight');fig.savefig(F/(name+'.svg'),bbox_inches='tight');plt.close(fig);captions[name]=caption
def style(ax,title,ylabel=None):
 ax.set_title(title,loc='left',pad=12);ax.grid(True);ax.set_axisbelow(True)
 if ylabel:ax.set_ylabel(ylabel)
def timeaxis(ax):ax.set_xlabel('Physical time (µs)');ax.set_xlim(0,50);ax.set_xticks(np.arange(0,51,10))
def style3(ax):
 ax.set(xlabel='X (µm)',ylabel='Y (µm)',zlabel='Z (µm)',xlim=(-380,380),ylim=(-290,290),zlim=(-670,670));ax.set_box_aspect((760,580,1340));ax.view_init(elev=15,azim=-58);ax.set_proj_type('ortho');ax.set_xticks([-300,0,300]);ax.set_yticks([-200,0,200]);ax.set_zticks([-500,0,500]);ax.tick_params(labelsize=11,pad=0)
 for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.fill=True;axis.pane.set_facecolor((.97,.98,.985,1));axis._axinfo['grid']['color']=(.35,.43,.49,.48)
def front(ax):ax.plot([-50,-50],[-250,250],[0,0],c=C[3],lw=2.8)
def lattice(ax,x,point_colors=None,speed=False):
 x=x*1e6;ax.add_collection3d(Line3DCollection(x[b[near&~cut]],colors=C[4],linewidths=.4,alpha=.35));ax.add_collection3d(Line3DCollection(x[b[near&cut]],colors=C[1],linewidths=.9,linestyles='--',alpha=.85))
 if point_colors is None:point_colors=np.where(grip,C[2],C[0])
 kw=dict(s=6,depthshade=False,linewidths=0,alpha=.9)
 if speed:kw.update(cmap='viridis',vmin=0,vmax=float(h['max_speed_m_s'].max()))
 sc=ax.scatter(*x.T,c=point_colors,**kw);front(ax);style3(ax);return sc

# 1. Actual full specimen and every original reference segment in a thin reference layer.
fig=plt.figure(figsize=(15,9));gs=fig.add_gridspec(2,2,width_ratios=[1,1.05],height_ratios=[1,.45],hspace=.35,wspace=.3)
ax=fig.add_subplot(gs[:,0],projection='3d');lattice(ax,X);ax.set_title('A  All 3,080 material points',loc='left',pad=12)
poly=np.array([[-350,-250,0],[-50,-250,0],[-50,250,0],[-350,250,0]]);ax.add_collection3d(Poly3DCollection([poly],facecolors=C[1],alpha=.15,edgecolors='none'))
for z,sgn in [(545,1),(-545,-1)]:ax.quiver(0,0,z,0,0,sgn*105,color=C[3],arrow_length_ratio=.3,linewidth=2)
fig.legend(handles=[Line2D([],[],marker='o',ls='',color=C[0],label='All-DOF free internal nodes'),Line2D([],[],marker='o',ls='',color=C[2],label='z-controlled grips; x,y free'),Line2D([],[],color=C[4],label='Live nearest reference bonds'),Line2D([],[],ls='--',color=C[1],label='Initial precut records'),Line2D([],[],color=C[3],label='Initial front / load directions')],loc='upper left',bbox_to_anchor=(.04,-.06),ncol=2,fontsize=11)
ax=fig.add_subplot(gs[0,1]);layer=np.isclose(X[:,1],-DX/2,rtol=0,atol=1e-15);sel=layer[b[:,0]]&layer[b[:,1]];xy=X[:,[0,2]]*1e6
ax.add_collection(LineCollection(xy[b[sel&~cut]],colors=C[4],linewidths=.35,alpha=.45));ax.add_collection(LineCollection(xy[b[sel&cut]],colors=C[1],linewidths=.7,linestyles='--',alpha=.75));ax.scatter(*xy[layer].T,c=np.where(grip[layer],C[2],C[0]),s=12,zorder=4);ax.plot([-350,-50],[0,0],c=C[3],lw=2.5);ax.scatter([-50],[0],s=65,c=C[3],zorder=6)
ax.add_patch(plt.Circle((-50,0),HORIZON*1e6,fill=False,color=C[2],ls=':',lw=1.8));ax.text(115,65,'δ = 150.5 µm',fontsize=11,color=C[2]);ax.set(xlim=(-380,380),ylim=(-570,570),xlabel='X (µm)',ylabel='Z (µm)');ax.set_aspect('equal');style(ax,'B  Actual bonds in Y = −25 µm layer')
ax=fig.add_subplot(gs[1,1]);ax.axis('off');ax.text(0,1,'700 × 500 × 1,100 µm cell-volume dimensions\n144,448 original bonds · 3,524 initial cuts\n4,576 live cross-plane bonds after precutting\nFront: x = −50 µm, z = 0; propagation toward +X\n2,800 fully free nodes · 280 normal-only grips\nLigament = 2.66δ; tip-to-grip = 3.49δ',va='top',fontsize=13,linespacing=1.7)
fig.suptitle('One frozen three-dimensional precracked specimen',fontsize=21,y=1.01)
save(fig,'Fig1_actual_specimen','全部3080点，左图仅显示原始最近邻键以便辨认，橙色虚线是已预切的参考记录，不是仍传力的连接。右图显示Y=−25 µm层内全部原始键（包括长键）；统计使用完整144448原始键。半透明矩形只是初始预切几何定义，不是生成的物理裂面。厚度边界和有限horizon效应仍存在。')
provenance['Fig1_actual_specimen']=dict(all_nodes=True,left_bond_ids=(np.flatnonzero(near)+1).tolist(),right_layer_bond_ids=(np.flatnonzero(sel)+1).tolist())

# 2. No first event exists; show actual subthreshold evidence explicitly.
fig=plt.figure(figsize=(15,8.5));gs=fig.add_gridspec(2,2,width_ratios=[1,1.2],height_ratios=[.8,1],hspace=.3,wspace=.3)
ax=fig.add_subplot(gs[:,0],projection='3d');lattice(ax,s['x']);r=events['highest_unfailed_kappa_bonds_at_last_state'][0];bid=r['bond_id']-1;q=s['x'][b[bid]]*1e6;ax.plot(*q.T,c=C[3],lw=4);ax.scatter(*q.T,c=C[3],s=35,depthshade=False);ax.set_title('A  Highest unfailed κ at 50 µs',loc='left');ax.text2D(.02,-.30,f'Bond {bid+1}; atom IDs {r["atom_ids"]}\nκ = {r["kappa"]:.7f}; d = 0\nRegion: {r["region"]}',transform=ax.transAxes,fontsize=12)
ax=fig.add_subplot(gs[0,1]);ax.axis('off');ax.text(0,.97,'First natural softening: NOT OBSERVED\nFirst complete failure: NOT OBSERVED\nFirst localized damage band: NOT OBSERVED',va='top',fontsize=16,linespacing=1.9,color=C[3]);ax.text(0,.23,'Observed window: 0–50 µs of planned 2,500 µs.\nInitial cuts (κ = sc, d = 1) are excluded.\nHighlighted bond remains undamaged; it is not a failure event.',va='top',fontsize=12,linespacing=1.5)
ax=fig.add_subplot(gs[1,1]);labels=['Tip interior','Tip side surface','Crack flank','Grip','Other free surface','Remote interior'];vals=[z['max_natural_kappa'] or 0 for z in loc['regions']];ax.barh(labels,vals,color=[C[0],C[0],C[1],C[2],C[4],C[5]]);ax.invert_yaxis();ax.set_xlabel('Maximum natural κ at 50 µs');ax.set_xlim(0,max(vals)*1.26);style(ax,'B  Subthreshold history by frozen region');ax.text(.5,.07,'Softening threshold s₀ = 0.15\n(outside this magnified scale)',transform=ax.transAxes,fontsize=11)
fig.suptitle('No natural damage within the accepted early-time window',fontsize=21,y=1.02);save(fig,'Fig2_event_not_observed','本次实际数据中不存在首次软化、首次完全失效或集中损伤事件。显示的是50 µs尚未失效键的最大κ与真实位置、ID、d=0；其位置不能替代首次破坏位置。所有初始预切键均从自然事件统计中排除。')

# 3. Six existing physical snapshots, not invented fracture milestones.
fig=plt.figure(figsize=(17,13));fig.subplots_adjust(left=.025,right=.97,top=.91,bottom=.15,wspace=.04,hspace=.17)
for k,idx in enumerate([0,2,4,6,8,10]):
 ax=fig.add_subplot(2,3,k+1,projection='3d');speed=np.linalg.norm(nodes[idx,1],axis=1);sc=lattice(ax,nodes[idx,0],speed,True);ax.set_title(f'{chr(65+k)}  {fi["time_s"][idx]*1e6:g} µs',loc='left',fontsize=17)
cax=fig.add_axes([.3,.065,.4,.019]);fig.colorbar(sc,cax=cax,orientation='horizontal',label='Actual node speed (m/s), common scale')
fig.suptitle('Actual early-time motion; no fracture milestones reached',fontsize=22);fig.text(.07,-.025,'All frames: 0 natural failures, one component. All 3,080 points shown; common camera, bounds and color scale.\nOnly nearest reference bonds drawn. Orange dashed traces are initial cuts. Numerical stop: 50 µs / planned 2,500 µs.',fontsize=13,color=C[5]);save(fig,'Fig3_accepted_early_motion','使用真实0、10、20、30、40、50 µs六个已保存态；相机、坐标范围和速度色标完全相同，位移无放大。显示全部节点和最近邻键；完整图拓扑使用全部参考键。实际尚未到达首次软化、断裂或传播阶段，不能把这些时间帧称为断裂阶段。')

# 4. Plane grids and actual opening. No fabricated crack surface.
fig,axs=plt.subplots(2,2,figsize=(15,9.8),layout='constrained');extent=(-350,350,-250,250)
for ax,key,title,cmap in [(axs[0,0],'initial_cut_count','A  Initial precut crossing records','Oranges'),(axs[0,1],'live_count','B  Surviving crossings at 50 µs','Blues')]:
 im=ax.imshow(plane[key].reshape(14,10).T,origin='lower',extent=extent,cmap=cmap,aspect='equal',interpolation='nearest');fig.colorbar(im,ax=ax,label='Bond intersections / 50 µm cell',shrink=.8);style(ax,title);ax.set(xlabel='X (µm)',ylabel='Y (µm)');ax.axvline(-50,c=C[3],ls='--',lw=2)
ax=axs[1,0];opening=plane['relative_opening_m'].reshape(14,10).T*1e9;lim=max(abs(opening).max(),1e-12);im=ax.imshow(opening,origin='lower',extent=extent,cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim),interpolation='nearest',aspect='equal');fig.colorbar(im,ax=ax,label='Actual relative opening (nm)',shrink=.8);style(ax,'C  Paired-face Δz − initial gap');ax.set(xlabel='X (µm)',ylabel='Y (µm)');ax.axvline(-50,c=C[3],ls='--',lw=2)
ax=axs[1,1];ax.axis('off');ax.text(0,.98,'CRACK AREA UNRESOLVED',color=C[3],fontsize=19,weight='bold',va='top');ax.text(0,.84,'8,100 original cross-plane bonds\n3,524 initial cuts (excluded from new area)\n4,576 surviving bridges; 0 new failures\n\n0 qualifying new separation cells\nNo observed front advance\nOpening criterion: 5 µm (5,000 nm)\nInitial projected cut area: 0.15 mm²\n\nNew physical crack area: UNAVAILABLE\nNo bond-count × dx² area conversion',va='top',fontsize=13,linespacing=1.5)
fig.suptitle('All original cross-plane bonds included, including long-range bridges',fontsize=20);save(fig,'Fig4_plane_and_bridging','按原始线段与z=0的真实交点统计全部8100条跨面键，并非按键中点猜测。3524初始切断、4576仍存活；新增自然失效0。开口为对应上下节点的真实z差减初始50 µm间距，单位nm，未放大成裂缝。网格边界浮点交点采用冻结的floor归属，单元计数纹理不是不对称损伤。没有满足冻结开口与无桥接条件的新生区域；投影新面积不可用，初始0.15 mm²不能用于D/A。')

# 5. All requested ledger channels and real physical time.
fig,axs=plt.subplots(3,2,figsize=(15,13),layout='constrained')
ax=axs[0,0];ax.plot(t,h['reaction_top_N']*1e6,c=C[0],lw=2,label='Material reaction');ax.plot(t,h['actuator_top_N']*1e6,c=C[1],lw=1.8,ls='--',label='Actuator incl. grip inertia');style(ax,'A  Top normal reaction','Force (µN)');ax.legend()
ax=axs[0,1];ax.plot(t,h['max_speed_m_s'],c=C[0],lw=2);style(ax,'B  Maximum actual node speed','Speed (m/s)')
ax=axs[1,0];ax.step(t,h['total_natural_failed'],where='post',c=C[3],lw=2,label='Cumulative natural failure');ax.set_ylim(-.1,1.1);style(ax,'C  New natural failure and damage','Natural failed bonds');ax.text(.05,.75,'Every step: new softening = 0\nEvery step: new complete failure = 0\nCumulative D = 0 J',transform=ax.transAxes,fontsize=12)
ax=axs[1,1]
for j,key in enumerate(['W_J','U_J','K_J','D_J','Uc_J','Jc_J']):ax.plot(t,h[key]*1e12,c=C[j],lw=2 if j<3 else 1.2,ls='-' if j<3 else ['--',':','-.'][j-3],label=key.split('_')[0])
style(ax,'D  Complete energy ledger','Energy (pJ)');ax.legend(ncol=3,fontsize=11)
ax=axs[2,0];ax.plot(t,h['R_J']*1e18,c=C[3],lw=2);style(ax,'E  Retained numerical balance R','R (aJ)');ax.text(.05,.12,'W = ΔU + ΔK + D + ΔUc − Jc + R',transform=ax.transAxes,fontsize=12)
ax=axs[2,1];ax.plot(t,h['relative_R']*100,c=C[0],lw=2);style(ax,'F  Cumulative normalized balance','Normalized |R| (%)');ax.text(.42,.08,f'Maximum: {energy["max_relative_R_percent"]:.6f}%\nFrozen gate: 2% (outside plotted scale)\nK / U at 50 µs = {energy["last_K_over_U"]:.3f}\nEarly dynamic transient, not quasistatic',transform=ax.transAxes,fontsize=12)
for ax in axs.flat:timeaxis(ax)
fig.suptitle('Accepted 0–50 µs mechanics only; full fracture protocol incomplete',fontsize=21);save(fig,'Fig5_reaction_motion_energy','全部横轴为真实µs，不是步号；上端材料反力与含夹持规定质量加速度的执行器反力分开。能量单位pJ，余额单位aJ；D、Uc、Jc实际均为0，重合曲线已说明。50 µs时K/U约2.61，不能称准静态。能量通过仅限已接受前缀，不能替代资源或完整协议通过。')

# 6. Actual last accepted topology and the failed admissions.
fig=plt.figure(figsize=(16,10));gs=fig.add_gridspec(2,2,width_ratios=[1,1.2],height_ratios=[.8,1],wspace=.3,hspace=.32)
ax=fig.add_subplot(gs[:,0],projection='3d');lattice(ax,s['x'],np.full(len(X),C[0]));ax.set_title('A  Last accepted topology: one component',loc='left');span=np.array(comp['components'][0]['actual_center_span_m'])*1e6;ax.text2D(.02,-.30,f'3,080 nodes; 140,924 live bonds\nActual center spans: {span[0]:.3f} × {span[1]:.3f} × {span[2]:.3f} µm\nAll original IDs retained; no new detached fragment',transform=ax.transAxes,fontsize=12)
ax=fig.add_subplot(gs[0,1]);mem=[4,resource['peak_RSS_bytes']/1024**3,resource['preflight']['peak_RSS_bytes']/1024**3];bars=ax.barh(['Frozen budget','Formal peak RSS','Preflight peak RSS'],mem,color=[C[4],C[3],C[1]]);ax.invert_yaxis();ax.set_xlim(0,13);ax.set_xlabel('GiB');style(ax,'B  Memory admission failed');ax.bar_label(bars,fmt='%.2f',padding=5);ax.text(.02,.02,'Measured high-water RSS',transform=ax.transAxes,fontsize=11,bbox=dict(facecolor='white',alpha=.95,edgecolor='none'))
ax=fig.add_subplot(gs[1,1]);ax.axis('off');rows=[('Accepted 0–50 µs mechanics','PASS (limited window)'),('Planned 2,500 µs protocol','NOT COMPLETED'),('Natural failure: tip / grip / other','0 / 0 / 0; order unknown'),('Spatial crack propagation','NOT OBSERVED'),('New projected crack area','UNRESOLVED'),('Protocol separation energy','UNAVAILABLE'),('Material Gc / macro strength','NOT VALIDATED'),('Final native endpoint gather','NOT PERFORMED after stop')];table=ax.table(cellText=rows,colWidths=[.54,.46],loc='center',cellLoc='left');table.auto_set_font_size(False);table.set_fontsize(11);table.scale(1,2.15)
for (i,j),cell in table.get_celld().items():cell.set_edgecolor('#D1DCE2');cell.set_facecolor('#F0F6F8' if i%2==0 else 'white')
ax.set_title('C  Scientific and engineering status',loc='left',fontsize=16)
fig.suptitle('Resource stop; natural crack / area–energy validation not completed',fontsize=21,y=1.01);save(fig,'Fig6_topology_and_admission','末个已接受状态仍为3080点单分量。独立全原始图重建与逐步原生存活键数一致，但异常后没有再调用原生末态端点gather；不把未做检查写成端点失配。主程序峰值约11.24 GiB超过4 GiB预算，P1早已约11.26 GiB而资源准入漏检。本轮未完成自然裂纹/面积/耗能验证，不改变前轮341点完整700 µs已通过的基线。')
dump(F/'FIGURE_CAPTIONS.json',captions);dump(F/'FIGURE_DATA_PROVENANCE.json',provenance)
cards=''.join(f'<section><h2>{name}</h2><a href="{name}.png"><img src="{name}.png" alt="{cap}"></a><p>{cap}</p><a href="{name}.svg">SVG</a> · <a href="{name}.png">PNG 原图</a></section>' for name,cap in captions.items())
report='../../../verification/'+TASK+'/'+RUN+'/REVIEW_'+TASK.upper()+'.md'
html='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>小型非均匀裂纹：资源停止与真实前缀</title><style>body{max-width:1550px;margin:30px auto;padding:0 25px;font:18px/1.65 system-ui;color:#233440;background:#f7f9fb}h1{color:#a73941}section{background:white;padding:22px;margin:28px 0;border:1px solid #d2dce4;border-radius:8px}img{width:100%;height:auto}a{color:#187F94}</style><h1>自然裂纹验证未完成：50 µs 资源预算停止</h1><p>3080真实材料点；144448原始键；3524初始预切。仅一次正式路径，200接受步，0拒绝；尚无自然软化或断裂，4576条跨面桥接存活。峰值RSS 11.24 GiB超过冻结4 GiB预算。此结论不能解释为材料不会断裂。</p><p>最大归一化能量余额0.002013%，仅认证前50 µs机械计算；新生裂面面积、协议分离耗能、材料Gc均不可用。所有图来自现有FP64数据，无人工裂面，无断裂阶段补画。</p><p><a href="'+report+'">完整34项科研审查报告</a></p>'+cards+'</html>'
(F/'OPEN_RESULTS.html').write_text(html);rel='../../../results/'+TASK+'/'+RUN+'/OPEN_RESULTS.html';(O/'OPEN_RESULTS.html').write_text('<!doctype html><meta charset="utf-8"><title>结果入口</title><p><a href="'+rel+'">打开六组真实科研图与说明</a></p><meta http-equiv="refresh" content="0;url='+rel+'">')
print('SIX_FIGURE_GROUPS',F)
