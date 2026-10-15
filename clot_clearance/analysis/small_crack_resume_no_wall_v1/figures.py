from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
import shutil
C=['#187F94','#C16A29','#7756A2','#C34049','#7F95A4','#314653']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.titlesize':15,'axes.labelsize':13,'xtick.labelsize':11,'ytick.labelsize':11,'legend.fontsize':11,'axes.edgecolor':'#46545F','axes.spines.top':False,'axes.spines.right':False,'grid.color':'#607080','grid.alpha':.4,'figure.facecolor':'white','savefig.facecolor':'white','svg.fonttype':'none'})
ref=dict(np.load(O/'provenance/reference.npz'));last=dict(np.load(O/'full/last_accepted.npz'));h=np.genfromtxt(O/'NATIVE_DYNAMIC_HISTORY.csv',delimiter=',',names=True);t=h['time_s']*1e3;status=json.loads((O/'FINAL_STATUS.json').read_text());events=json.loads((O/'FIRST_SOFTENING_AND_FAILURE.json').read_text());loc=json.loads((O/'CRACK_LOCALIZATION_AUDIT.json').read_text());energy=json.loads((O/'DYNAMIC_ENERGY_LEDGER.json').read_text());resource=json.loads((O/'RUN_RESOURCES.json').read_text());topology=json.loads((O/'NATIVE_COMPONENT_TOPOLOGY.json').read_text());contact=json.loads((O/'CONTACT_STREAM_AUDIT.json').read_text());ph=np.load(O/'PLANE_HISTORY.npz');X=ref['X'];b=ref['bonds'];cut=ref['initial_cut'];near=ref['r0']<=1.01*DX;N=len(X);captions={};F.mkdir(exist_ok=True,parents=True)
def save(fig,name,caption):
 fig.savefig(F/(name+'.png'),dpi=230,bbox_inches='tight');fig.savefig(F/(name+'.svg'),bbox_inches='tight');plt.close(fig);captions[name]=caption
def style(ax,title,ylabel):ax.set_title(title,loc='left',pad=12);ax.set_ylabel(ylabel);ax.grid(True);ax.set_axisbelow(True)
def timeaxis(ax):ax.set_xlabel('Physical time (ms)');ax.set_xlim(0,float(last['time'])*1e3)
# Common bounds include every actual saved node, not just reference geometry.
fi=np.genfromtxt(O/'full/native_frames.csv',delimiter=',',names=True);nodes=np.memmap(O/'full/native_frames.bin',dtype='<f8',mode='r',shape=(len(fi),4,N,3));lo=nodes[:,0].min(axis=(0,1))*1e6;hi=nodes[:,0].max(axis=(0,1))*1e6;lo=np.minimum(lo,[-350,-250,-550])-40;hi=np.maximum(hi,[350,250,550])+40
def style3(ax):
 ax.set(xlabel='X (µm)',ylabel='Y (µm)',zlabel='Z (µm)',xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]),zlim=(lo[2],hi[2]));ax.set_box_aspect(hi-lo);ax.view_init(15,-58);ax.set_proj_type('ortho');ax.tick_params(labelsize=10,pad=0)
 ax.set_xticks([-300,0,300]);ax.set_yticks([-200,0,200]);ax.set_zticks([-600,-300,0,300,600])
 for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.fill=True;axis.pane.set_facecolor((.97,.98,.985,1));axis._axinfo['grid']['color']=(.35,.43,.49,.45)
def lattice(ax,s,highlight=None):
 q=s['x']*1e6;d=s['d'];alive=near&(d<1);fail=near&(d==1)&~cut
 ax.add_collection3d(Line3DCollection(q[b[alive]],colors=plt.colormaps['YlOrRd'](d[alive]),linewidths=.8,alpha=.65))
 if fail.any():ax.add_collection3d(Line3DCollection(q[b[fail]],colors=C[3],linestyles='--',linewidths=.5,alpha=.28))
 ax.add_collection3d(Line3DCollection(q[b[near&cut]],colors=C[1],linestyles=':',linewidths=.55,alpha=.4))
 bb=b[d<1];g=coo_matrix((np.ones(len(bb)),(bb[:,0],bb[:,1])),shape=(N,N)).tocsr();nc,lab=connected_components(g,directed=False);sizes=np.bincount(lab);rank=np.empty(nc,int);rank[np.argsort(sizes)[::-1]]=np.arange(nc);colors=[C[int(rank[j])%len(C)] for j in lab];ax.scatter(*q.T,s=7,c=colors,depthshade=False,linewidths=0,alpha=.9)
 ax.plot([-50,-50],[-250,250],[0,0],c=C[3],lw=2.5,ls=':')
 if highlight is not None and len(highlight):
  hb=b[np.asarray(highlight)-1];ax.add_collection3d(Line3DCollection(q[hb],colors=C[3],linewidths=3.2,alpha=1));ax.scatter(*q[np.unique(hb)].T,s=75,c=C[3],edgecolors='white',linewidths=.7,depthshade=False)
 style3(ax);return nc

oldfig=R/'results'/TASK/'20261014T210608Z'
for ext in ['png','svg']:shutil.copy2(oldfig/('Fig1_actual_specimen.'+ext),F/('Fig1_actual_specimen.'+ext))
captions['Fig1_actual_specimen']='几何、原始键与初始缺口完全不变，复用上一归档的结构图；全部3080点，三维图只绘最近邻参考键，右侧层内绘全部参考键。橙色为初始预切记录，不是新自然断裂。'

fig=plt.figure(figsize=(15,10));fig.subplots_adjust(left=.02,right=.98,bottom=.28,top=.89,wspace=.16)
for k,(field,title) in enumerate([('first_softening','First natural softening'),('first_complete_failure','First complete failure')]):
 ev=events[field];s=last if ev is None else dict(np.load(O/'full/keyframes'/f'step_{ev["step"]:06d}.npz'));ax=fig.add_subplot(1,2,k+1,projection='3d');lattice(ax,s,None if ev is None else ev['bond_ids']);ax.set_title(title+('\nNOT OBSERVED' if ev is None else f'\n{ev["time_s"]*1e3:.6f} ms'),loc='left',pad=15)
 if ev is None:txt='No event within the accepted observation window.'
 else:
  ids=ev['bond_ids'];reg=', '.join(sorted(set(ev['regions'])));txt=f'Original bond IDs: '+', '.join(map(str,ids[:8]))+(' …' if len(ids)>8 else '')+f'\n{len(ids)} simultaneous bonds; region: {reg}\nκ: {min(ev["kappa"]):.6f}–{max(ev["kappa"]):.6f}; d: {min(ev["d"]):.6f}–{max(ev["d"]):.6f}\nActual stretch: {min(ev["stretch"]):.6f}–{max(ev["stretch"]):.6f}'
 fig.text(.055+.48*k,.12,txt,fontsize=12,linespacing=1.7,va='top')
fig.suptitle('Actual event locations and immutable original bond identities',fontsize=21);save(fig,'Fig2_first_events','首次自然软化与完全失效的实际原生状态、原始bond ID、参考区域和κ/d。全部初始3524预切被排除；图中事件位置不移动，若未发生则明确标记。所有同时发生事件的完整ID与伸长在FIRST_SOFTENING_AND_FAILURE.json。')

keys=json.loads((O/'full/KEYFRAME_INDEX.json').read_text());choices=[dict(time_s=0,step=0,title='Initial precut')]
for field,title in [('first_softening','First softening'),('first_complete_failure','First failure'),('first_localized_cloud','Localized cloud'),('first_candidate_front_advance','Front advance candidate')]:
 ev=events[field]
 if ev is not None:
  if field=='first_candidate_front_advance':choices.append(dict(time_s=ev['time_s'],step=ev['step'],title=title,state_path=str(O/'front_candidate_visual_state.npz')))
  else:
   key=min(keys,key=lambda z:abs(z['time_s']-ev['time_s']));choices.append(dict(time_s=key['time_s'],step=key['step'],title=title if abs(key['time_s']-ev['time_s'])<1e-15 else 'Saved state near candidate'))
for q in keys:
 if q['reason']=='first_open_no_bridge_cell' and len(choices)<5:choices.append(dict(time_s=q['time_s'],step=q['step'],title='First open/no-bridge cell'))
for q in keys:
 if q['reason']=='fixed_observation_time' and len(choices)<5:choices.append(dict(time_s=q['time_s'],step=q['step'],title='Fixed observation'))
for q in sorted(keys,key=lambda z:z['time_s']):
 if len(choices)>=5:break
 if q['step'] not in [z['step'] for z in choices]:choices.append(dict(time_s=q['time_s'],step=q['step'],title='Actual saved state'))
choices=sorted(choices[:5],key=lambda z:z['time_s'])+[dict(time_s=float(last['time']),step=int(last['step']),title='Last accepted')]
fig=plt.figure(figsize=(17,13));fig.subplots_adjust(left=.02,right=.98,bottom=.15,top=.92,wspace=.04,hspace=.28)
for k,q in enumerate(choices):
 p=Path(q['state_path']) if 'state_path' in q else O/'full/keyframes'/f'step_{q["step"]:06d}.npz';s=dict(np.load(p));ax=fig.add_subplot(2,3,k+1,projection='3d');nc=lattice(ax,s);ax.set_title(f'{chr(65+k)}  {q["title"]}\n{q["time_s"]*1e3:.6g} ms · {nc} components',loc='left',fontsize=14,pad=8)
cax=fig.add_axes([.32,.05,.36,.018]);fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap='YlOrRd'),cax=cax,orientation='horizontal',label='Surviving bond damage d (common scale)');fig.suptitle('Actual lattice evolution; no displacement magnification',fontsize=21);fig.text(.05,-.015,'All 3,080 points shown. Nearest bonds drawn; connectivity uses all original bonds. Dashed red: new full failure.\nDotted orange: initial precut. Dotted front line is a reference marker. Missing fracture stages are not fabricated.',fontsize=12);save(fig,'Fig3_actual_evolution','从实际完整态选取初态、真实首次事件、已保存时刻及末接受态，全部采用相同物理范围、相机和d色标；位移不放大。图中只绘最近邻，分量用全部原始存活图重建。未发生的集中或传播阶段不补画；分量颜色按当前大小排序，不能据此跨帧追踪同一碎块身份。')
dump(F/'EVOLUTION_FRAME_SELECTION.json',choices)

fig,axs=plt.subplots(2,2,figsize=(15,10),layout='constrained');extent=(-350,350,-250,250)
for ax,data,title,label,cmap in [(axs[0,0],ph['live_counts'][-1].T,'A  All surviving cross-plane bridges','Live crossings / cell','Blues'),(axs[0,1],(ph['initial_cross_counts']-ph['initial_cut_counts']-ph['live_counts'][-1]).T,'B  New fully failed cross-plane bonds','Natural failed crossings / cell','Reds'),(axs[1,0],ph['opening_m'][-1].T*1e6,'C  Actual paired-face relative opening','Opening (µm)','RdBu_r')]:
  im=ax.imshow(data,origin='lower',extent=extent,aspect='equal',interpolation='nearest',cmap=cmap);fig.colorbar(im,ax=ax,shrink=.8,label=label);ax.axvline(-50,c=C[3],ls='--',lw=2,label='Initial front');ys=ph['cell_centers_xy_m'].reshape(14,10,2)[0,:,1]*1e6;ax.plot(ph['front_x_m'][-1]*1e6,ys,c='#191919',lw=1.7,marker='.',label='No-bridge/opening prefix');style(ax,title,'Y (µm)');ax.set_xlabel('X (µm)')
axs[0,0].legend(fontsize=8,loc='upper left')
ax=axs[1,1];ax.plot(ph['time_s']*1e3,ph['notch_connected'].sum(axis=(1,2))*DX**2*1e6,c=C[0],label='No-bridge/opening prefix proxy');ax.plot(ph['time_s']*1e3,ph['conservative'].sum(axis=(1,2))*DX**2*1e6,c=C[1],ls='--',label='Conservative horizon-clear subset');style(ax,'D  CRACK AREA UNRESOLVED','Grid proxy area (mm²)');ax.set_xlabel('Physical time (ms)');ax.legend(fontsize=10);ax.text(.03,.75,'Initial cut area excluded\nGrid proxies are not a unique physical surface',transform=ax.transAxes,fontsize=11)
fig.suptitle(f'Last accepted: {status["live_cross_bonds"]} surviving bridges; all original lengths counted',fontsize=20);save(fig,'Fig4_plane_and_bridging','全部跨面原始线段按真实交点统计，初始预切与自然完全失效分开；所有长程桥接纳入。原冻结floor网格边界归属保留，计数纹理不等于不对称损伤。真实开口减去初始点距；新单元代理要求无桥接、开口和从缺口连续，保守子集再检查horizon邻域。代理面积尚非唯一物理裂面，不计算任意D/截面Gc。')

fig,axs=plt.subplots(3,2,figsize=(15,13),layout='constrained')
ax=axs[0,0];ax.plot(t,h['reaction_top_N']*1e6,c=C[0],label='Material reaction');ax.plot(t,h['actuator_top_N']*1e6,c=C[1],ls='--',label='Actuator incl. grip inertia');style(ax,'A  Actual top normal reaction','Force (µN)');ax.legend()
ax=axs[0,1];ax.plot(t,h['max_speed_m_s'],c=C[0]);style(ax,'B  Maximum node speed','Speed (m/s)')
ax=axs[1,0];ax.step(t,h['total_natural_failed'],where='post',c=C[3],label='Cumulative natural failure');style(ax,'C  Natural failure and dissipation','New failed bonds (cumulative)');ax2=ax.twinx();ax2.plot(t,h['D_J']*1e9,c=C[2],label='D');ax2.set_ylabel('Cumulative D (nJ)',color=C[2]);ax.legend(loc='upper left')
ax=axs[1,1]
for j,key in enumerate(['W_J','U_J','K_J','D_J','Uc_J','Jc_J']):ax.plot(t,h[key]*1e9,c=C[j],lw=1.8,ls='-' if j<4 else '--',label=key.split('_')[0])
style(ax,'D  Full dynamic energy ledger','Energy (nJ)');ax.legend(ncol=3)
ax=axs[2,0];ax.plot(t,h['R_J']*1e12,c=C[3]);style(ax,'E  Retained numerical balance R','R (pJ)')
ax=axs[2,1];ax.plot(t,h['relative_R']*100,c=C[0]);ax.axhline(2,c=C[3],ls='--',label='Unchanged 2% gate');style(ax,'F  Cumulative normalized balance','Normalized |R| (%)');ax.legend();ax.text(.04,.12,f'Max accepted: {energy["max_relative_R_percent"]:.6g}%\nW = ΔU + ΔK + D + ΔUc − Jc + R',transform=ax.transAxes,fontsize=11)
for ax in axs.flat:timeaxis(ax)
fig.suptitle('Measured dynamics and energy; no residual relabeled as damage',fontsize=21);save(fig,'Fig5_reaction_motion_energy','全程采用真实物理ms；材料反力与含夹持惯性的执行器反力分开。D来自真实变分更新，R独立保留；初始预切不计D。图只包含已接受态，若发生门槛停止，不将未接受候选绘入合法曲线；动态惯性仍需单独讨论。')

fig=plt.figure(figsize=(16,10));gs=fig.add_gridspec(2,2,width_ratios=[1,1.2],hspace=.4,wspace=.35);ax=fig.add_subplot(gs[:,0],projection='3d');lattice(ax,last);ax.set_title(f'A  Last accepted: {len(topology["components"])} components',loc='left');sizes=topology['component_sizes'];ax.text2D(.02,-.25,'Sizes: '+', '.join(map(str,sizes[:10]))+(' …' if len(sizes)>10 else '')+'\nAll 3,080 original material IDs retained\nMain component span: '+' × '.join(f'{v*1e3:.3f}' for v in topology['components'][0]['center_span_m'])+' mm',transform=ax.transAxes,fontsize=11)
ax=fig.add_subplot(gs[0,1]);labels=['Tip interior','Tip side','Flank','Grip','Other surface','Interior'];ax.barh(labels,[z['naturally_failed'] for z in loc['regions']],color=[C[0],C[0],C[1],C[2],C[4],C[5]]);ax.invert_yaxis();style(ax,'B  Actual competing failure regions','');ax.set_xlabel('Natural fully failed bonds')
ax=fig.add_subplot(gs[1,1]);ax.axis('off');rows=[('Full 2.5 ms protocol','PASS' if status['full_protocol_completed'] else 'NOT COMPLETED'),('Observed classification',loc['classification'].replace('_',' ')),('Stop reason','END OF FROZEN PROTOCOL' if status['full_protocol_completed'] else status['stop_category']),('RSS / wall-clock ceiling','NONE / NONE'),('Protocol peak RSS',f'{resource["protocol_peak_RSS_bytes"]/1024**3:.3f} GiB'),('Projected new surface','UNRESOLVED'),('G separation / material Gc','UNAVAILABLE'),('Max contact overlap',f'{contact["max_overlap_dx"]:.6g} dx'),('Native final endpoint check','PASS' if topology['native_final_endpoint_match'] else 'NOT AVAILABLE')];table=ax.table(cellText=rows,colWidths=[.49,.51],loc='center',cellLoc='left');table.auto_set_font_size(False);table.set_fontsize(10.5);table.scale(1,2)
for (i,j),cell in table.get_celld().items():cell.set_edgecolor('#D1DCE2');cell.set_facecolor('#F0F6F8' if i%2==0 else 'white')
ax.set_title('C  Numerical and scientific status',loc='left');fig.suptitle('Native restart continuation of the original 2.5 ms protocol',fontsize=21);save(fig,'Fig6_topology_and_admission','末接受状态全原始图连通分量、各预注册区域的真实自然失效数及独立准入判断。取消RSS和墙钟上限不改变能量/损伤/运动学门槛；数值协议通过与自然裂纹面积/耗能认证分开。接触重叠来自实际数据，不因规则未改就宣称物理校准。')
dump(F/'FIGURE_CAPTIONS.json',captions)
cards=''.join(f'<section><h2>{name}</h2><a href="{name}.png"><img src="{name}.png" alt="{cap}"></a><p>{cap}</p><a href="{name}.svg">SVG</a> · <a href="{name}.png">PNG</a></section>' for name,cap in captions.items())
html='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>取消墙钟上限后的原生续算结果</title><style>body{max-width:1550px;margin:32px auto;padding:0 24px;font:18px/1.65 system-ui;color:#233440;background:#f7f9fb}section{background:white;padding:22px;margin:28px 0;border:1px solid #d2dce4}img{width:100%;height:auto}a{color:#187F94}</style><h1>同一预裂纹试样：从2.4 ms原生检查点续算</h1><p>Run '+RUN+' · '+status['status']+' · '+loc['classification']+'</p><p>实际到 '+f'{float(last["time"])*1e3:.6g}'+' ms；'+str(status['natural_failed'])+' 个新增完全失效键；'+str(status['live_cross_bonds'])+' 条跨面桥接存活。新裂面和协议分离耗能仍须独立验收，当前Gc不可用。</p>'+cards+'</html>'
html=html.replace('</p>','</p><p><a href="../../../verification/'+TASK+'/'+RUN+'/REVIEW_'+TASK.upper()+'.md">完整科研报告（34项回答）</a></p>',1)
(F/'OPEN_RESULTS.html').write_text(html);rel='../../../results/'+TASK+'/'+RUN+'/OPEN_RESULTS.html';(O/'OPEN_RESULTS.html').write_text('<!doctype html><meta charset="utf-8"><p><a href="'+rel+'">六组科研图及说明</a></p><meta http-equiv="refresh" content="0;url='+rel+'">');print('SIX_FIGURES_COMPLETE')
