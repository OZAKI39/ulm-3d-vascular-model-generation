from common import *
from external import graph_components,NativeYWallQuery
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection,Poly3DCollection
from matplotlib.ticker import MaxNLocator
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.titlesize':15,'axes.labelsize':12,'xtick.labelsize':11,'ytick.labelsize':11,'legend.fontsize':11,'axes.edgecolor':'#46545F','axes.spines.top':False,'axes.spines.right':False,'grid.color':'#607080','grid.alpha':.4,'figure.facecolor':'white','savefig.facecolor':'white','svg.fonttype':'none'})
C=['#187F94','#D38635','#7756A2','#C34049','#91A0AB','#314653'];caps={}
def read(name):return json.loads((O/name).read_text())
ref=dict(np.load(O/'provenance/reference.npz'));last=dict(np.load(O/'full/last_accepted.npz'));h=np.atleast_1d(np.genfromtxt(O/'V3_NATIVE_PILOT_HISTORY.csv',delimiter=',',names=True));t=h['time_s']*1e3;st=read('FINAL_STATUS.json');ev=read('FIRST_EVENTS.json');en=read('V3_ENERGY_LEDGER.json');top=read('V3_FRAGMENT_TOPOLOGY.json');rs=read('RUN_RESOURCES.json');ha=read('V3_BACKGROUND_DRAG_HANDOFF.json');co=read('HISTORICAL_COMPARISON.json');au=read('NUMERICAL_AUDIT.json');X=ref['X'];b=ref['bonds'];bt=ref['bond_types'];n=int(ref['clot_count']);N=len(X);F.mkdir(exist_ok=True,parents=True)
keys=read('full/KEYFRAME_INDEX.json');fi=np.atleast_1d(np.genfromtxt(O/'full/native_frames.csv',delimiter=',',names=True));frames=np.memmap(O/'full/native_frames.bin',dtype='<f8',mode='r',shape=(len(fi),4,N,3));lo=frames[:,0].min(axis=(0,1))*1e3-.06;hi=frames[:,0].max(axis=(0,1))*1e3+.06
# Physical coordinates, same orthographic camera and aspect for all evolution panels.
def axes3(ax,bounds=None):
 ll,hh=(lo,hi) if bounds is None else bounds
 ax.set(xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)',xlim=(ll[0],hh[0]),ylim=(ll[1],hh[1]),zlim=(ll[2],hh[2]));ax.set_box_aspect(hh-ll);ax.view_init(13,-67);ax.set_proj_type('ortho');ax.tick_params(labelsize=10,pad=0)
 for a in [ax.xaxis,ax.yaxis,ax.zaxis]:a.set_major_locator(MaxNLocator(4));a.pane.fill=True;a.pane.set_facecolor((.97,.98,.985,1));a._axinfo['grid']['color']=(.32,.4,.46,.65)
def state(step):
 p=O/'full/keyframes'/f'step_{step:06d}.npz'
 return last if step==int(last['step']) else dict(np.load(p))
def cloud(ax,s,highlight=None,material=True,bounds=None):
 q=s['x']*1e3;lab,sz,at=graph_components(n,b,bt,s['d']<1);det=~at[lab];tiny=det&(sz[lab]<=2)
 ax.scatter(*q[n:].T,s=4,c=C[4],alpha=.32,depthshade=False,linewidths=0)
 for mask,color,ss in [(~det,C[0],9),(det&~tiny,C[2],21),(tiny,C[3],20)]:
  if mask.any():ax.scatter(*q[:n][mask].T,s=ss,c=color,depthshade=False,linewidths=.15,edgecolors='white')
 if highlight is not None and len(highlight):
  ids=np.asarray(highlight,int);ax.scatter(*q[ids-1].T,s=80,c=C[1],depthshade=False,edgecolors='black',linewidths=.6)
 axes3(ax,bounds);return int(det.sum()),len(sz[~at])
def save(fig,name,caption):
 fig.savefig(F/(name+'.png'),dpi=240,bbox_inches='tight');fig.savefig(F/(name+'.svg'),bbox_inches='tight');plt.close(fig);caps[name]=caption;print('FIGURE',name,flush=True)
def legend(fig,y=.02):
 fig.legend([Line2D([0],[0],marker='o',ls='',c=c,markersize=8) for c in [C[0],C[4],C[2],C[3]]],['Attached material','Fixed anchors','Detached component (≥3 nodes)','Unresolved debris (1–2 nodes)'],loc='lower center',bbox_to_anchor=(.5,y),ncol=4,frameon=False,fontsize=12)
def style(ax,title,label):ax.set_title(title,loc='left',pad=10);ax.set_ylabel(label);ax.grid(True);ax.set_axisbelow(True)
def times(ax):ax.set_xlabel('Physical time (ms)');ax.set_xlim(0,t[-1])

fig=plt.figure(figsize=(18,9));fig.subplots_adjust(left=.015,right=.98,bottom=.22,top=.87,wspace=.1)
ax=fig.add_subplot(1,2,1,projection='3d');initial=state(0);bounds=(X.min(0)*1e3-.04,X.max(0)*1e3+.04);cloud(ax,initial,bounds=bounds)
near=ref['r0']<=1.01*DX
for typ,color in [(1,C[0]),(2,C[1])]:
 keep=near&(bt==typ);ax.add_collection3d(Line3DCollection(X[b[keep]]*1e3,colors=color,linewidths=.55,alpha=.5))
ax.set_title('A  Actual V3 lattice and nearest typed bonds\n1,927 material points + 1,561 fixed anchors',loc='left')
ax=fig.add_subplot(1,2,2,projection='3d');query=NativeYWallQuery();tri=query.tri*1e3;ax.add_collection3d(Poly3DCollection(tri,facecolors='#65869A',edgecolors='none',alpha=.22,rasterized=True));ax.scatter(*X[:n].T*1e3,s=12,c=C[3],depthshade=False);wallpts=tri.reshape(-1,3);axes3(ax,(wallpts.min(0),wallpts.max(0)));ax.view_init(29,-75);ax.set_title('B  Exact location in the native Y vessel\nMother wall near X = −6 mm, +Z side',loc='left')
fig.legend([Line2D([0],[0],marker='o',ls='',c=C[0]),Line2D([0],[0],marker='o',ls='',c=C[4]),Line2D([0],[0],c=C[0],lw=2),Line2D([0],[0],c=C[1],lw=2)],['Material point','Fixed anchor','Bulk bond','Interface bond'],loc='center',bbox_to_anchor=(.27,.255),ncol=4,frameon=False,fontsize=11)
fig.suptitle('V3 variational softening: geometry and one frozen material card',fontsize=22)
fig.text(.055,.12,'Bulk: 82,270 bonds · s₀ = 0.15 · s_c = 0.18\nInterface: 22,040 bonds · s₀ = 0.140625 · s_c = 0.16875 (uncalibrated)\nFixed–fixed scaffold: 53,327 bonds · no integrated anchor motion',fontsize=13,linespacing=1.6)
fig.text(.56,.12,'G = 940 Pa · K = 33,333.33 Pa · ρ = 1,000 kg/m³\ndx = 50 µm · horizon = 150.5 µm\nNearest bonds drawn; all 157,637 reference bonds used in mechanics.',fontsize=13,linespacing=1.6)
save(fig,'Fig1_V3_material_and_bonds','真实1927材料点、1561锚点、157637原始键和原始Y形壁面三角网格。左图只显示最近邻本体/界面线段以便看清，力学与图连通性始终使用全部参考键；未改变几何。界面sc=0.16875为明确未标定工程假设。')

fig=plt.figure(figsize=(21,14));fig.subplots_adjust(left=.02,right=.98,bottom=.11,top=.91,wspace=.12,hspace=.52)
items=[('bulk_softening','Bulk softening'),('bulk_failure','Bulk complete failure'),('interface_softening','Interface softening'),('interface_failure','Interface complete failure'),('detachment','First detachment')]
for k,(key,title) in enumerate(items):
 ax=fig.add_subplot(2,3,k+1,projection='3d');event=ev[key]
 if event is None:
  cloud(ax,last);ax.set_title(title+'\nNOT OBSERVED',loc='left');txt=f'Accepted window: 0–{t[-1]:.6g} ms\nNo fabricated event or location.'
 else:
  ss=state(event['step']);ids=np.asarray(event['atom_ids']).ravel();cloud(ax,ss,ids);ax.set_title(title+f'\n{event["time_s"]*1e3:.6g} ms · step {event["step"]}',loc='left');bs=event.get('bond_ids',[])
  if bs:ax.add_collection3d(Line3DCollection(ss['x'][b[np.asarray(bs)-1]]*1e3,colors=C[1],linewidths=3))
  txt='Original atom IDs: '+', '.join(map(str,ids[:7]))+(' …' if len(ids)>7 else '')+'\nOriginal bond IDs: '+(', '.join(map(str,bs[:6]))+(' …' if len(bs)>6 else '') if bs else 'topology event; see JSON')
 ax.text2D(.03,-.15,txt,transform=ax.transAxes,fontsize=11,linespacing=1.6)
ax=fig.add_subplot(2,3,6);ax.axis('off');ax.text(.03,.8,'Every marker is from an accepted native state.\n\nStable original atom / bond IDs are retained.\n\nFull simultaneous-event IDs and κ / d:\nFIRST_EVENTS.json\n\nNo displacement magnification.\nNo artificial continuous clot shell.',fontsize=15,linespacing=1.7,va='top');legend(fig,.025);fig.suptitle('First typed damage and actual topological detachment',fontsize=22)
save(fig,'Fig2_first_damage_events','五类首次事件采用该接受步FP64节点及κ/d全态；标注原始atom/bond ID。不存在的事件明确写未观察到，拓扑脱离不是某一条任意键的失效替代。完整同时事件ID在FIRST_EVENTS.json。')

choices=[dict(step=0,time_s=0,title='Initial')]
for key,title in [('bulk_softening','First damage'),('detachment','First detachment')]:
 e=ev[key]
 if e:choices.append(dict(step=e['step'],time_s=e['time_s'],title=title))
 else:choices.append(dict(step=int(last['step']),time_s=float(last['time']),title=title+' · NOT OBSERVED'))
obs=[q for q in keys if q['reason']=='fixed_observation']
for target,title in [(float(last['time'])*.35,'Intermediate'),(min(.00125,float(last['time'])*.8),'Later accepted state')]:
 q=min(obs,key=lambda q:abs(q['time_s']-target)) if obs else dict(step=int(last['step']),time_s=float(last['time']));choices.append(dict(step=q['step'],time_s=q['time_s'],title=title))
choices.append(dict(step=int(last['step']),time_s=float(last['time']),title='Last accepted'))
fig=plt.figure(figsize=(19,13));fig.subplots_adjust(left=.015,right=.98,bottom=.12,top=.92,wspace=.06,hspace=.25)
for k,q in enumerate(choices):
 ax=fig.add_subplot(2,3,k+1,projection='3d');nd,nc=cloud(ax,state(q['step']));ax.set_title(f'{chr(65+k)}  {q["title"]}\n{q["time_s"]*1e3:.6g} ms · {nd} detached nodes / {nc} components',loc='left',fontsize=13)
legend(fig,.025);fig.suptitle('Actual V3 fracture point clouds · fixed physical bounds and camera',fontsize=22)
save(fig,'Fig3_actual_V3_evolution','各帧使用同一真实物理坐标、正交相机及轴比例，显示全部3488节点，不放大位移，不绘连续医学外壳。附着/游离分类用全部d<1原始键与固定锚点连接关系；1–2点为未解析数值碎屑，非临床微栓子。');dump(F/'EVOLUTION_FRAME_SELECTION.json',choices)

fig,axs=plt.subplots(2,3,figsize=(18,11));fig.subplots_adjust(left=.065,right=.98,bottom=.21,top=.88,wspace=.32,hspace=.43);labels=['Original\ninstant-break','Dual-fix\ninstant-break','New typed\nsoftening'];cases=co['cases']
for ax,(key,title,unit,scale) in zip(axs.flat,[('detached_mass_kg','Detached material mass','µg',1e9),('components','Free connected components','Count',1),('largest_nodes','Largest free component','Material nodes',1),('tiny_nodes','Unresolved 1–2-node debris','Material nodes',1),('N_ge200','Free components ≥200 µm','Count',1),('N_ge300','Free components ≥300 µm','Count',1)]):
 vals=[q[key]*scale for q in cases];bars=ax.bar(labels,vals,color=[C[4],C[1],C[0]],width=.64);style(ax,title,unit);ax.bar_label(bars,fmt='%.3g',padding=4);ax.set_ylim(0,max(max(vals)*1.25,1.2));ax.tick_params(axis='x',labelsize=11)
fig.suptitle('Historical comparison using attachment-aware, surviving-bond graphs',fontsize=21)
fig.text(.07,.065,'Same initial geometry and C1 artificial load. Original historical run: 217 nodes / 208 components; corrected native run: 199 / 190.\nNew bulk and interface softening laws, interface onset rounding (+1 ppm), event integration and MPI layout differ.\nExploratory comparison only; changes cannot be attributed solely to bulk softening. '+('All three windows end at 2.33 ms.' if st['full_protocol_completed'] else f'WARNING: new case stops at {t[-1]:.6g} ms; terminal counts are not a same-time comparison.'),fontsize=12,linespacing=1.65)
save(fig,'Fig4_historical_fragment_comparison','只读复用旧217/208及双Bug修正199/190两套结果，未重跑旧模型。比较脱离质量、图分量、最大节点数、未解析碎屑及200/300µm阈值。等效直径由N dx³定义；旧新损伤律、界面onset约1ppm、dt及MPI不同，因此不作单一机制定量归因。')

fig,axs=plt.subplots(4,2,figsize=(18,19));fig.subplots_adjust(left=.075,right=.93,bottom=.09,top=.94,hspace=.49,wspace=.32)
ax=axs[0,0];ax.plot(t,h['load_sum_norm_N']*1e6,c=C[0],label='Σ |actual nodal load|');style(ax,'A  Current-position artificial loading','Force (µN)');ax2=ax.twinx();ax2.plot(t,h['load_amp'],c=C[1],ls='--');ax2.set_ylabel('Frozen temporal envelope',color=C[1]);ax.legend(loc='lower right')
ax=axs[0,1]
for key,c,lab in [('bulk_failed',C[3],'Bulk'),('interface_failed',C[1],'Interface')]:ax.step(t,h[key],where='post',c=c,label=lab)
style(ax,'B  Fully failed original bonds','Cumulative count');ax.legend()
ax=axs[1,0];ax.plot(t,h['max_speed_m_s'],c=C[0]);style(ax,'C  Maximum material-point speed','Speed (m/s)')
ax=axs[1,1];ax.plot(t,h['max_dd'],c=C[2],lw=.9);ax.axhline(.1,ls='--',c=C[3],label='Unchanged Δd ≤ 0.10 gate');style(ax,'D  Accepted damage increments','Maximum Δd');ax.legend()
ax=axs[2,0]
for j,key in enumerate(['W_J','U_J','K_J','D_J','Uc_J','Uw_J']):ax.plot(t,h[key]*1e9,c=C[j],label=key.split('_')[0],ls='-' if j<4 else '--')
style(ax,'E  Full mechanical ledger','Energy (nJ)');ax.legend(ncol=3)
ax=axs[2,1]
for j,(key,lab) in enumerate([('W_load_J','Artificial load'),('W_drag_J','Background drag'),('W_numeric_J','Attached numerical damping'),('W_wall_damping_J','Wall damping')]):ax.plot(t,h[key]*1e9,c=C[j],label=lab)
style(ax,'F  External work kept separate from damage','Work (nJ)');ax.legend(fontsize=10)
ax=axs[3,0]
for j,(key,lab) in enumerate([('Uc_J','Particle Uc'),('Uw_J','Wall Uw'),('Jc_J','Contact activation Jc'),('R_J','Numerical residual R')]):ax.plot(t,h[key]*1e12,c=C[j],label=lab)
style(ax,'G  Contact and numerical contributions','Energy (pJ)');ax.legend(ncol=2,fontsize=10)
ax=axs[3,1];ax.plot(t,h['relative_R']*100,c=C[0]);ax.axhline(2,ls='--',c=C[3],label='Unchanged 2% gate');style(ax,'H  Cumulative normalized energy balance','|R| / scale (%)');ax.legend();ax.text(.04,.45,f'Max accepted: {en["max_relative_R_percent"]:.6g}%',transform=ax.transAxes)
for ax in axs.flat:times(ax)
fig.suptitle('Accepted native dynamics and independently checked energy ledger',fontsize=22)
fig.text(.08,.025,'Wload + Wdrag + Wnumeric + Wwall,damp = ΔU + ΔK + D + ΔUc + ΔUw − Jc + R\nD is the actual irreversible variational energy release. R is retained separately; rejected candidates are excluded.',fontsize=13,linespacing=1.6)
save(fig,'Fig5_load_motion_energy','实际外力按当前坐标计算且按历史时程卸载；真实物理ms。绘接受态的本体/界面完全失效、速度、Δd和完整W/U/K/D/接触/R；背景阻力、附着数值阻尼、壁面阻尼功单列，绝不混入材料D。R使用冻结动态能量尺度；拒绝候选不进入合法曲线。')

fig,axs=plt.subplots(1,2,figsize=(18,10),gridspec_kw={'width_ratios':[1,1.12]});fig.subplots_adjust(left=.045,right=.97,bottom=.08,top=.86,wspace=.12)
rows=[('Physical window',f'0–{t[-1]:.6g} ms'),('Total pilot wall time',f'{rs["main_wall_s"]:.1f} s / {rs["main_wall_s"]/60:.2f} min'),('Accepted physical steps',str(rs['accepted_steps'])),('Accepted event substeps',str(rs['substep_count'])),('Rejected damage trials',str(rs['rejected_trials'])),('Actual dt range',f'{rs["actual_dt_min_s"]*1e6:.5g}–{rs["actual_dt_max_s"]*1e6:.5g} µs'),('Peak process RSS',f'{rs["peak_RSS_bytes"]/1024**3:.3f} GiB'),('Wall cost / physical ms',f'{rs["wall_s_per_physical_ms"]:.1f} s/ms'),('MPI / threads / GPU','1 / 1 / no'),('Wall / RSS ceilings','None / none'),('Formal cases','Exactly one')]
rows2=[('Numerical protocol','PASS' if st['full_protocol_completed'] else 'STOPPED; accepted prefix only'),('Scientific label','EXPLORATORY_UNCALIBRATED'),('Exploratory sample schema','Saved with full provenance'),('Formal AI training admission','NOT GRANTED'),('Critical load','UNAVAILABLE'),('Particle contact max overlap',f'{ha["maximum_particle_contact_overlap_dx"]:.3f} dx; uncalibrated'),('Y-wall penalty contact','EXERCISED' if ha['actual_wall_contact_exercised'] else 'NOT EXERCISED'),('Material Gc','UNAVAILABLE'),('Bulk macro strength','UNVALIDATED'),('Interface softening tail','UNCALIBRATED'),('Fluid model','One-way background COM drag'),('True microbubble streaming','NOT CONNECTED'),('Artificial Pa → ultrasound power','NO SUCH MAPPING')]
for ax,rows,title in [(axs[0],rows,'A  Measured computational cost'),(axs[1],rows2,'B  Numerical and scientific admission')]:
 ax.axis('off');tab=ax.table(cellText=rows,colWidths=[.51,.49],cellLoc='left',loc='center');tab.auto_set_font_size(False);tab.set_fontsize(11);tab.scale(1,2.65)
 for (i,j),cell in tab.get_celld().items():cell.set_edgecolor('#D1DCE2');cell.set_facecolor('#F0F6F8' if i%2==0 else 'white')
 ax.set_title(title,loc='left',pad=10,fontsize=17)
fig.suptitle('One real V3 pilot: measured cost and a bounded exploratory AI record',fontsize=22)
save(fig,'Fig6_cost_and_AI_admission','统计包括原生计算、独立力/能量检查、回滚及检查点的正式pilot墙钟，单进程RSS实测。未启用固定RSS/墙钟上限；原输出和数值门槛保留。成功数值窗口不等于材料Gc、宏观强度、界面尾部或临界载荷认证。仅一个可追溯探索样本，不能直接形成正式AI训练集。')
dump(F/'FIGURE_CAPTIONS.json',caps)
cards=''.join(f'<section><h2>{name}</h2><a href="{name}.png"><img src="{name}.png" alt="{cap}"></a><p>{cap}</p><a href="{name}.svg">SVG</a> · <a href="{name}.png">PNG</a></section>' for name,cap in caps.items())
html='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>V3渐进软化单工况结果</title><style>body{max-width:1700px;margin:32px auto;padding:0 24px;font:18px/1.7 system-ui;color:#233440;background:#f7f9fb}section{background:white;padding:22px;margin:28px 0;border:1px solid #d2dce4}img{width:100%;height:auto}a{color:#187F94}</style><h1>真实V3血栓：变分软化集成与唯一正式pilot</h1><p>'+RUN+' · '+st['status']+' · EXPLORATORY_UNCALIBRATED</p><p>实际到 '+f'{t[-1]:.6g}'+' ms，脱离 '+str(top['detached_material_points'])+' 个材料点、'+str(top['detached_component_count'])+' 个游离分量。等效直径≥200µm分量 '+str(top['N_detached_ge200'])+'。1–2点碎屑不能称为临床微栓子。</p><p><a href="../../../verification/'+TASK+'/'+RUN+'/REVIEW_'+TASK.upper()+'.md">完整中文报告（28项回答）</a></p>'+cards+'</html>'
(F/'OPEN_RESULTS.html').write_text(html);rel='../../../results/'+TASK+'/'+RUN+'/OPEN_RESULTS.html';(O/'OPEN_RESULTS.html').write_text('<!doctype html><meta charset="utf-8"><p><a href="'+rel+'">六组科研图及说明</a></p><meta http-equiv="refresh" content="0;url='+rel+'">');print('SIX_FIGURES_COMPLETE',flush=True)
