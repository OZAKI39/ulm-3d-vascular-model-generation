"""All formal members, explicit CORE/FULL comparisons, 300 dpi PNG and PDF."""
from pathlib import Path
import argparse,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection,Line3DCollection
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.formal_cohort_p9a5 import *
R=ROOT/REL;D=R/'data';F=R/'figures'
font=R/'assets/NotoSansCJKsc-Regular.otf'
if font.exists():font_manager.fontManager.addfont(str(font));plt.rcParams['font.family']=font_manager.FontProperties(fname=font).get_name()
plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white','pdf.fonttype':42,'axes.unicode_minus':False})
COLORS={'O1':'#1767a6','O2':'#1c8066','O3':'#d56821','SUPPORTED_STATIONARY':'#b94148','LONG_RESIDENCE_CENSORED':'#8059a4','SOLVER_FAILURE':'#282c34'}
LABELS={'O1':'O1','O2':'O2','O3':'O3','SUPPORTED_STATIONARY':'有接触支持的静止','LONG_RESIDENCE_CENSORED':'长驻留 / 删失','SOLVER_FAILURE':'求解失败'}


def load():
 summary=json.loads((D/'final_summary.json').read_text());co=json.loads((D/'FINAL_FORMAL_COHORT.json').read_text())
 rows=merge_rows(json.loads((D/'analysis_metrics.json').read_text()),[e['particle_id'] for e in co['events']])
 assert len(rows)==summary['final_formal_N']==co['count']>=500
 assert digest(D/'FINAL_FORMAL_COHORT.json')==summary['final_formal_cohort_sha']
 return summary,co,rows


def category(r):return r['outlet'] if r['status']=='COMPLETED' else r['status']


def save(fig,name):
 for extension in ['png','pdf']:fig.savefig(F/(name+'.'+extension),dpi=300,bbox_inches='tight',facecolor='white')
 plt.close(fig);print('FIGURE_READY',name,flush=True)


def geometry():
 import pyvista as pv
 base=ROOT/'particle_3d/reports/network_derived_flow_mb_validation_v1/server_bundle/formal_3D_flow_solver/FEM_SimVascular'
 manifest=json.loads((base/'frozen_reference/boundary_manifest.json').read_text())
 surfaces={k:pv.read(base/v['path']) for k,v in manifest['boundaries'].items()}
 wall=surfaces['WALL'];triangles=wall.points[wall.faces.reshape(-1,4)[:,1:]]*1e6
 return wall,triangles,surfaces


def bounds(ax,points):
 lo=points.min(axis=0);hi=points.max(axis=0);margin=(hi-lo)*.07
 ax.set_xlim(lo[0]-margin[0],hi[0]+margin[0]);ax.set_ylim(lo[1]-margin[1],hi[1]+margin[1]);ax.set_zlim(lo[2]-margin[2],hi[2]+margin[2])
 ax.set_box_aspect(hi-lo);ax.view_init(elev=22,azim=-65)
 ax.set(xlabel='x (µm)',ylabel='y (µm)',zlabel='z (µm)')


def main():
 s,co,rows=load();n=len(rows);present=[k for k in COLORS if any(category(r)==k for r in rows)]
 reps={k:min(r['particle_id'] for r in rows if category(r)==k) for k in present}
 arrays={r['particle_id']:np.load(R/r['track_relative_path']/'trajectory.npz')['samples'] for r in rows}
 wall,triangles,surfaces=geometry()
 # 01: EVERY formal member; a second panel exposes deterministic representatives.
 fig=plt.figure(figsize=(14,6.5));axes=[fig.add_subplot(121,projection='3d'),fig.add_subplot(122,projection='3d')]
 for ax in axes:
  ax.add_collection3d(Poly3DCollection(triangles,facecolor='#bdc8d1',edgecolor='none',alpha=.025,rasterized=True));bounds(ax,wall.points*1e6)
 for k in present:
  paths=[]
  for r in rows:
   if category(r)==k:
    xyz=arrays[r['particle_id']][:,1:4]*1e6;indices=np.unique(np.linspace(0,len(xyz)-1,min(180,len(xyz)),dtype=int));paths.append(xyz[indices])
  axes[0].add_collection3d(Line3DCollection(paths,colors=COLORS[k],linewidths=.55,alpha=max(.04,min(.2,50/n))))
  r=next(r for r in rows if r['particle_id']==reps[k]);xyz=arrays[r['particle_id']][:,1:4]*1e6
  axes[1].plot(*xyz.T,color=COLORS[k],lw=2.2,label=f'{LABELS[k]} · ID {r["particle_id"]}')
 for name,surf in surfaces.items():
  if name.startswith('OUTLET_'):
   center=surf.points.mean(axis=0)*1e6;role='O'+str(int(name[-2:]))
   for ax in axes:ax.text(*center,role,color=COLORS[role],fontsize=11,fontweight='bold')
 axes[0].set_title(f'全体 {n:,} 条；低透明度显示，不删除任何成员');axes[1].set_title('各结果类别中 particle_id 最小的代表轨迹');axes[1].legend(loc='upper left',fontsize=8,frameon=False)
 fig.suptitle(f'Figure 01 · 正式 {n:,} 条微泡在三维血管中的运动路径',fontsize=16);fig.tight_layout();save(fig,'Figure_01_all_formal_trajectories')
 # 02: sample order, not execution finish order. Coverage is not a quota.
 fig,ax=plt.subplots(figsize=(11,5));ids=np.arange(1,n+1)
 for k in ['O1','O2','O3']:
  cumulative=np.cumsum([r['status']=='COMPLETED' and r['outlet']==k for r in rows]);ax.plot(ids,cumulative,color=COLORS[k],lw=2,label=f'{k}: {cumulative[-1]}')
  first=s['first_'+k+'_at_N']
  if first is not None:ax.scatter([first],[1],color=COLORS[k],s=35);ax.annotate(f'FIRST_{k}: N={first}',(first,1),xytext=(6,12+18*['O1','O2','O3'].index(k)),textcoords='offset points',color=COLORS[k],fontsize=9)
  else:ax.text(.02,.92-.07*['O1','O2','O3'].index(k),f'{k} 未观察到',color=COLORS[k],transform=ax.transAxes)
 ax.axvline(500,color='#6b737c',ls='--',label='500 minimum threshold')
 for checkpoint in s['population_checkpoints'][1:]:ax.axvline(checkpoint['N'],color='#bfc8d0',lw=.6,alpha=.5)
 ax.set(xlabel='Formal particle count N (accepted order)',ylabel='Cumulative completed count',title='Figure 02 · 随着无偏样本增加，三个出口何时自然出现？');ax.legend(frameon=False);ax.grid(alpha=.12);fig.tight_layout();save(fig,'Figure_02_natural_outlet_coverage')
 # 03: both denominators and both uncertainty conventions remain explicit.
 fig,axes=plt.subplots(1,2,figsize=(12,4.8));x=np.arange(3);width=.36
 for j,(label,g) in enumerate([('CORE500',s['core500']),(f'FULL {n}',s['full_formal'])]):
  vals=np.array([g['fractions'][o]['fraction'] for o in ['O1','O2','O3']]);ci=np.array([[g['fractions'][o]['lower'],g['fractions'][o]['upper']] for o in ['O1','O2','O3']])
  pos=x+(j-.5)*width;axes[0].bar(pos,[g[o] for o in ['O1','O2','O3']],width,label=label,color=['#6a91b5','#db935e'][j]);axes[1].bar(pos,vals*100,width,label=label,color=['#6a91b5','#db935e'][j]);axes[1].errorbar(pos,vals*100,yerr=np.maximum(0,np.array([vals-ci[:,0],ci[:,1]-vals])*100),fmt='none',ecolor='#34485d',capsize=3)
 for ax in axes:ax.set_xticks(x,['O1','O2','O3']);ax.legend(frameon=False)
 axes[0].set(ylabel='Completed trajectories',title='出口数量（CORE 与 FULL 分开）');axes[1].set(ylabel='Fraction of all formal members (%)',title='出口比例与不确定区间')
 fig.suptitle('Figure 03 · CORE500 与完整正式样本的出口结果');fig.text(.5,.01,'CORE：固定样本 Wilson 95%；FULL：覆盖6类别×13个预定检查点的保守同时区间（允许按覆盖情况停止）',ha='center',fontsize=8);fig.tight_layout(rect=(0,.05,1,.95));save(fig,'Figure_03_core_full_outlets')
 # 04: actual checkpoints, including the reason for the final N.
 checkpoints=s['population_checkpoints'];fig,ax=plt.subplots(figsize=(11,max(3.5,.65*len(checkpoints)+1.5)));ax.axis('off')
 for i,c in enumerate(checkpoints):
  y=1-(i+1)/(len(checkpoints)+1);why='三出口全部自然观察到，停止扩样' if c['all_outlets_naturally_observed'] else ('达到预先声明的资源上限' if i==len(checkpoints)-1 else '仍缺出口，继续后续连续 births')
  ax.text(.03,y,f'N = {c["N"]:,}\n本批 +{c["added"]}',va='center',fontsize=12,bbox=dict(boxstyle='round,pad=.5',facecolor='#eaf0f5',edgecolor='#b9c9d6'))
  ax.text(.24,y,f'O1={c["O1"]}   O2={c["O2"]}   O3={c["O3"]}   静止={c["stationary"]}   长驻留={c["censored"]}\n{why}',va='center',fontsize=11)
 ax.set_title(f'Figure 04 · 最后为什么停在 {n:,} 颗？',fontsize=15,pad=15)
 fig.text(.5,.02,'样本增加不是为了改变比例，只是为了提高看到稀有出口实例的机会。入口分布和已选成员保持不变。',ha='center',fontsize=10);save(fig,'Figure_04_sample_expansion_timeline')
 # 05: size versus outcome with all members, no outcome-selected diameter draw.
 fig,axes=plt.subplots(1,2,figsize=(12,4.5));diam=np.array([r['diameter_um'] for r in rows]);bins=np.linspace(diam.min()-.03,diam.max()+.03,24)
 axes[0].hist([[r['diameter_um'] for r in rows if category(r)==k] for k in present],bins=bins,stacked=True,color=[COLORS[k] for k in present],label=[LABELS[k] for k in present],edgecolor='white',lw=.35);axes[0].set(xlabel='Diameter (µm)',ylabel='Formal members',title='全部成员的粒径与结果')
 for j,k in enumerate(present):
  rr=[r for r in rows if category(r)==k];jitter=np.array([((r['particle_id']*37)%101-50)/300 for r in rr]);axes[1].scatter([r['diameter_um'] for r in rr],j+jitter,s=10,alpha=.35,color=COLORS[k])
 axes[1].set(yticks=range(len(present)),yticklabels=[LABELS[k] for k in present],xlabel='Diameter (µm)',title='每个点是一条正式轨迹');axes[0].legend(fontsize=8,frameon=False);fig.suptitle('Figure 05 · 微泡大小与最终结果');fig.tight_layout();save(fig,'Figure_05_diameter_vs_outcome')
 # 06: project original centers onto a deterministic best-fit inlet plane.
 xyz=np.array([e['birth_center_m'] for e in co['events']]);origin=xyz.mean(axis=0);_,_,basis=np.linalg.svd(xyz-origin,full_matrices=False);xy=(xyz-origin)@basis[:2].T*1e6
 fig,ax=plt.subplots(figsize=(8,6))
 for k in present:
  mask=np.array([category(r)==k for r in rows]);ax.scatter(xy[mask,0],xy[mask,1],s=13,alpha=.55,color=COLORS[k],label=f'{LABELS[k]} ({int(mask.sum())})')
 ax.set(xlabel='Inlet in-plane coordinate 1 (µm)',ylabel='Inlet in-plane coordinate 2 (µm)',title='Figure 06 · 原始入口位置与事后观察的结果',aspect='equal');ax.legend(frameon=False,bbox_to_anchor=(1,1));fig.tight_layout();save(fig,'Figure_06_inlet_position_vs_outlet')
 # 07: completed transit times; retain non-completions elsewhere, never assign them fake exit times.
 fig,axes=plt.subplots(1,2,figsize=(11,4.5));outlets=[o for o in ['O1','O2','O3'] if s['full_formal'][o]>0]
 for o in outlets:
  a=np.sort([r['transit_time_s'] for r in rows if r['outlet']==o and r['status']=='COMPLETED']);axes[0].step(a,np.arange(1,len(a)+1)/len(a),where='post',color=COLORS[o],label=f'{o} (n={len(a)})')
 axes[0].set(xlabel='Transit time (s)',ylabel='Empirical CDF',title='已完成轨迹的通行时间');axes[0].legend(frameon=False)
 if outlets:axes[1].boxplot([[r['transit_time_s'] for r in rows if r['outlet']==o and r['status']=='COMPLETED'] for o in outlets],tick_labels=outlets,showfliers=True)
 axes[1].set(ylabel='Transit time (s)',title='出口之间的分布');fig.suptitle('Figure 07 · 通行时间：只对真正穿过出口的轨迹定义');fig.tight_layout();save(fig,'Figure_07_transit_by_outlet')
 # 08: all wall gap and exposure metrics, in physical units.
 fig,axes=plt.subplots(1,2,figsize=(12,4.5))
 for k in present:
  rr=[r for r in rows if category(r)==k];axes[0].scatter([r['diameter_um'] for r in rr],[r['minimum_wall_gap_m']*1e9 for r in rr],s=12,alpha=.5,color=COLORS[k],label=LABELS[k]);axes[1].scatter([r['residence_time_s'] for r in rr],[r['nearwall_exposure_s'] for r in rr],s=12,alpha=.5,color=COLORS[k])
 axes[0].set(xlabel='Diameter (µm)',ylabel='Minimum real-wall gap (nm)',title='独立原始壁面距离检查');axes[0].set_yscale('symlog',linthresh=1.);axes[0].axhline(0,color='black',lw=.7);axes[0].legend(frameon=False,fontsize=8)
 axes[1].set(xlabel='Observed trajectory age (s)',ylabel='Near-wall exposure (s)',title='近壁定义：gap / radius ≤ 0.1');fig.suptitle('Figure 08 · 最小壁间隙与近壁停留');fig.tight_layout();save(fig,'Figure_08_wall_and_nearwall')
 # 09: same cohort inspected at historical/new age limits, not an old/new flow comparison.
 fig,ax=plt.subplots(figsize=(9,4.6));h=s['horizon_resolution'];bottom=np.zeros(4)
 for key,color,label in [('completed','#1c8066','已完成'),('stationary','#b94148','接触支持静止'),('failed','#282c34','失败'),('still_active_or_censored','#8059a4','仍运动 / 删失')]:
  vals=[v[key] for v in h];ax.bar(range(4),vals,bottom=bottom,color=color,label=label);bottom+=vals
 ax.set(xticks=range(4),xticklabels=['1.5 s（旧上限）','3 s','6 s','12 s'],ylabel='Formal members',title='Figure 09 · 同一正式队列在不同年龄上限下的结果可辨识性');ax.legend(frameon=False,bbox_to_anchor=(1,1));fig.text(.5,.01,f'由保存轨迹作截断审核；实际最大使用上限 {s["maximum_horizon_used_s"]:g} s。已终止状态不会为了填满时间而继续积分。',ha='center',fontsize=8);fig.tight_layout(rect=(0,.05,1,1));save(fig,'Figure_09_horizon_resolution')
 # 10: post-cohort point routing; the longer frozen point diagnostic budget is explicit.
 point_labels=['O1','O2','O3','NO_EXIT'];matrix=np.array([[sum(r['point_outlet']==p and category(r)==k for r in rows) for k in present] for p in point_labels]);fig,ax=plt.subplots(figsize=(9,5));im=ax.imshow(matrix,cmap='Blues',aspect='auto')
 for i in range(matrix.shape[0]):
  for j in range(matrix.shape[1]):ax.text(j,i,str(matrix[i,j]),ha='center',va='center',color='white' if matrix[i,j]>matrix.max()*.55 else '#203044')
 ax.set(xticks=range(len(present)),xticklabels=[LABELS[k] for k in present],yticks=range(4),yticklabels=point_labels,xlabel='Finite-size MB outcome',ylabel='Point-tracer outcome',title='Figure 10 · 同一出生中心：微泡与流体路径的事后比较');fig.colorbar(im,ax=ax,label='Count');fig.text(.5,.01,'Point 使用原有 30 s / 2 mm 诊断上限；在最终 MB cohort 冻结后计算，未用于选择成员。',ha='center',fontsize=8);fig.tight_layout(rect=(0,.05,1,1));save(fig,'Figure_10_mb_vs_point')
 # 11: measured CPU scaling (including oversubscription controls).
 benchmark_figure(s['worker_scaling_results'])
 # 12: only when an actual expansion occurred.
 if n>500:
  fig,axes=plt.subplots(1,2,figsize=(12,4.5));keys=['O1','O2','O3','stationary','censored','failed'];delta=[100*s['core_full_fraction_difference'][k] for k in keys]
  axes[0].barh(keys,delta,color=['#1767a6' if v>=0 else '#d56821' for v in delta]);axes[0].axvline(0,color='#65717d',lw=.7);axes[0].set(xlabel='FULL − CORE500 (percentage points)',title='结果比例变化；两组样本相互嵌套')
  for label,g,color in [('CORE500',rows[:500],'#1767a6'),(f'FULL {n}',rows,'#d56821')]:
   a=np.sort([r['diameter_um'] for r in g]);axes[1].step(a,np.arange(1,len(a)+1)/len(a),where='post',color=color,label=label)
  axes[1].set(xlabel='Diameter (µm)',ylabel='CDF',title='粒径分布稳定性');axes[1].legend(frameon=False);fig.suptitle('Figure 12 · 扩样后，关键统计发生了多大变化？');fig.tight_layout();save(fig,'Figure_12_core_full_stability')
 manifest=dict(formal_N=n,formal_cohort_sha256=digest(D/'FINAL_FORMAL_COHORT.json'),all_drawn_particle_ids=[r['particle_id'] for r in rows],
   representative_particle_ids=reps,representative_rule='SMALLEST_PARTICLE_ID_WITHIN_COMPLETED_OUTLET_OR_TERMINAL_CATEGORY',
   figure01_all_members=True,geometry_only_display_decimation='At most 180 vertices per path; includes endpoints; all particles retained; raw data unchanged',
   output_files={p.name:digest(p) for p in sorted(F.glob('Figure_*')) if p.suffix in ['.png','.pdf']})
 (D/'visualization_manifest.json').write_bytes(canonical(manifest))


def benchmark_figure(bench):
 a=bench['configurations'];workers=[r['workers'] for r in a];rates=[r['tracks_per_second'] for r in a]
 fig,axes=plt.subplots(1,2,figsize=(11,4.6));axes[0].plot(workers,rates,'o-',color='#1767a6');best=bench['workers_production'];idx=workers.index(best);axes[0].scatter([best],[rates[idx]],s=130,facecolor='none',edgecolor='#d56821',lw=2,label=f'Production: {best} workers');axes[0].legend(frameon=False)
 axes[0].axvline(bench['CPU_quota_cores'],color='#65717d',ls='--',label='CPU quota');axes[0].set(xlabel='CPU workers',ylabel='Tracks / second',xticks=workers,title='固定同一24条轨迹的吞吐量')
 axes[1].plot(workers,[r['peak_tree_pss_bytes']/2**30 for r in a],'o-',color='#1c8066');axes[1].set(xlabel='CPU workers',ylabel='Peak process-tree PSS (GiB)',xticks=workers,title='实测内存占用（共享页按比例计）')
 fig.suptitle(f'Figure 11 · CPU worker scaling；容器配额 {bench["CPU_quota_cores"]:.2f} 核；轨迹无 GPU kernel');fig.tight_layout();save(fig,'Figure_11_worker_scaling')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--benchmark-only',action='store_true');a=p.parse_args()
 if a.benchmark_only:benchmark_figure(json.loads((D/'worker_scaling_results.json').read_text()))
 else:main()
