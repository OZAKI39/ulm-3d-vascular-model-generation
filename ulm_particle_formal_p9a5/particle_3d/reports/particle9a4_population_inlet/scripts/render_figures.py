"""English scientific figures: white background, 300 dpi PNG and vector PDF."""
from pathlib import Path
import argparse,csv,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.injection_method_c import TruncatedSonoVue
from particle_3d.particle82a_geometry import maximum_handoff_radius
R=ROOT/'particle_3d/reports/particle9a4_population_inlet';D=R/'data';F=R/'figures'
plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white','pdf.fonttype':42})
BLUE='#1767a6';ORANGE='#d56821';GRAY='#65717d';GREEN='#1c8066';RED='#b94148'

def save(fig,name):
 for ext in ['png','pdf']:
  p=F/(name+'.'+ext)
  if p.exists():raise FileExistsError(p)
  fig.savefig(p,dpi=300,bbox_inches='tight',facecolor='white')
 plt.close(fig)

def read(p):return json.loads(Path(p).read_text())
def write_csv(path,rows):
 with Path(path).open('x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def inlet_figures():
 rows=[json.loads(line) for line in (D/'inlet100k/proposal_ledger.jsonl').read_text().splitlines()];s=read(D/'inlet100k/summary.json');c=read(ROOT/'particle_3d/contracts/P9A4_CONTINUOUS_INFUSION_V1.json');qacc=read(D/'qacc_diagnostic.json');dist=TruncatedSonoVue(ROOT/'sonovue_size_distribution_v0')
 diam=np.array([r['diameter_um'] for r in rows]);keep=np.array([r['particle_id'] is not None for r in rows]);xyz=np.array([r['position_m'] for r in rows]);t=np.array([r['proposal_time_s'] for r in rows]);tri=np.array([r['inlet_triangle_id'] for r in rows]);reasons=np.array([r['rejection_reason'] or 'ACCEPTED' for r in rows]);g=np.load(D/'inlet_geometry.npz');xy=(xyz-g['origin'])@g['basis'].T*1e6
 # 01: vector diagram, exact algorithm semantics.
 fig,ax=plt.subplots(figsize=(12,4.8));ax.set(xlim=(0,12),ylim=(0,4.7));ax.axis('off')
 def box(x,y,w,h,label,color):
  ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.12',facecolor='white',edgecolor=color,linewidth=1.7));ax.text(x+w/2,y+h/2,label,ha='center',va='center',color=color,fontsize=11)
 def arrow(a,b,color=GRAY):ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=15,color=color,lw=1.5))
 box(.2,2.6,2.4,1.2,'Poisson source proposals\n'+r'$\lambda_0=C_{\leq4}\,Q_{in}$'+'\nOne source_event_id',BLUE)
 box(3.05,2.6,2.2,1.2,'Draw diameter once\n'+r'$D\sim f(D\mid D\leq4\,\mu m)$',BLUE)
 box(5.7,2.6,2.2,1.2,'Draw position once\n'+r'$p(x)\propto(\mathbf{u}\cdot\mathbf{n})_+$',BLUE)
 box(8.35,2.6,3.1,1.2,'One finite-size check\nReal WALL + original handoff\nOpen INLET is not solid',GRAY)
 arrow((2.75,3.2),(2.95,3.2));arrow((5.4,3.2),(5.6,3.2));arrow((8.05,3.2),(8.25,3.2))
 box(5.2,.35,2.7,1.05,'REJECTED source event\nRetain time, size and position\nNo retry; no particle_id',RED)
 box(8.65,.35,2.7,1.05,'ACCEPTED birth\nAssign next particle_id\nUnchanged P9-A.1 dynamics',GREEN)
 arrow((9.3,2.4),(6.6,1.5),RED);arrow((10.4,2.4),(10.,1.5),GREEN)
 ax.text(.2,1.35,'Direct Poisson thinning\n\nAll proposals stay in the ledger.\nActive set = {} (independent tracks).',va='top',fontsize=11)
 ax.set_title('P9-A.4 | Steady continuous infusion with finite-size flux admission',fontsize=15,pad=12);save(fig,'Figure_01_p9a4_concept')
 # 02: distinguish source and conditional entering laws.
 fig,ax=plt.subplots(1,2,figsize=(11,4.2));bins=np.linspace(.5,4,36)
 ax[0].hist(diam,bins=bins,density=True,histtype='step',lw=2,color=BLUE,label='Source proposals');ax[0].hist(diam[keep],bins=bins,density=True,histtype='step',lw=2,color=ORANGE,label='Accepted births');ax[0].set(xlabel='Diameter (µm)',ylabel='Probability density (µm⁻¹)',title='Source and entering size PDFs');ax[0].legend(frameon=False)
 grid=np.linspace(.5,4,150);control=np.load(D/'qacc_positions.npz');dmax=2*maximum_handoff_radius(control['wall_distance_m']+qacc['wall_roundoff_m']);den=dist.cdf(dmax).mean()
 ax[1].plot(grid,dist.cdf(grid*1e-6),color=BLUE,lw=2,label='Conditional SonoVue source');ax[1].plot(np.sort(diam[keep]),np.arange(1,keep.sum()+1)/keep.sum(),color=ORANGE,lw=2,label='Accepted empirical CDF')
 expected=np.array([dist.cdf(np.minimum(d*1e-6,dmax)).mean()/den for d in grid]);ax[1].plot(grid,expected,'--',color=GRAY,lw=1.4,label='Independent flux prediction')
 ax[1].set(xlabel='Diameter (µm)',ylabel='Cumulative probability',title='Finite-size admission favors smaller diameters',ylim=(0,1.02));ax[1].legend(frameon=False,fontsize=8.5)
 fig.suptitle(f'100,000 source proposals → {keep.sum():,} accepted births');fig.tight_layout();save(fig,'Figure_02_source_entering_size')
 # 03: normalized 2D position laws; same color scale and physical coordinates.
 fig,axes=plt.subplots(1,2,figsize=(10,4.2),layout='constrained');edges=[np.linspace(xy[:,j].min(),xy[:,j].max(),41) for j in range(2)];dens=[]
 for mask in [np.ones(len(xy),bool),keep]:
  count,_,_=np.histogram2d(xy[mask,0],xy[mask,1],bins=edges);dens.append(count/(mask.sum()*np.diff(edges[0])[0]*np.diff(edges[1])[0]))
 vmax=max(x.max() for x in dens)
 for ax,density,title in zip(axes,dens,['Source: full positive-flux cap','Accepted: finite-size accessible population']):
  im=ax.pcolormesh(*edges,density.T,cmap='viridis',vmin=0,vmax=vmax,shading='flat');ax.set(xlabel='Inlet coordinate ξ (µm)',ylabel='Inlet coordinate η (µm)',title=title,aspect='equal')
 fig.colorbar(im,ax=list(axes),label='Normalized spatial density (µm⁻²)',shrink=.86);save(fig,'Figure_03_inlet_position_density')
 # 04: common-position diagnostic, pointwise Wilson 95% intervals.
 fig,ax=plt.subplots(figsize=(7.5,4.6));q=qacc['rows'];x=np.array([r['diameter_um'] for r in q]);y=np.array([r['fraction'] for r in q]);ci=np.array([r['ci95'] for r in q])
 ax.errorbar(x,y,yerr=np.maximum(0,np.array([y-ci[:,0],ci[:,1]-y])),fmt='o-',color=GREEN,capsize=3,lw=1.6);ax.set(xlabel='Bubble diameter (µm)',ylabel=r'$Q_{acc}(D)/Q_{in}$',title='Accessible inlet flux decreases with finite size',ylim=(-.015,1.02));ax.grid(alpha=.16);ax.text(.98,.95,'20,000 independent full-flux positions\nCommon positions across sizes\nPointwise Wilson 95% intervals',ha='right',va='top',transform=ax.transAxes,fontsize=9);fig.tight_layout();save(fig,'Figure_04_accessible_flux')
 # 05: first 120 physical proposals. Rejected time is never recycled.
 fig,axes=plt.subplots(2,1,figsize=(11,4.8),gridspec_kw={'height_ratios':[1.2,1]});n=120;chosen=keep[:n];axes[0].eventplot([t[:n],t[:n][chosen]],lineoffsets=[1,0],linelengths=.65,colors=[BLUE,ORANGE],linewidths=1.1);axes[0].set(yticks=[0,1],yticklabels=['Accepted births','Source proposals'],xlabel='Physical source time (s)',title='Poisson thinning: rejected proposals remain on the source clock');axes[0].set_xlim(0,t[n-1])
 dt=np.array([r['arrival_delta_t_s'] for r in rows]);x=np.linspace(0,np.quantile(dt,.995),160);axes[1].hist(dt,bins=np.linspace(0,x[-1],65),density=True,color=BLUE,alpha=.35,label='Source interarrival distribution');axes[1].plot(x,c['lambda_source_s_inv']*np.exp(-c['lambda_source_s_inv']*x),color=BLUE,lw=2,label='Exponential model');axes[1].set(xlabel='Interarrival time (s)',ylabel='Density (s⁻¹)');axes[1].legend(frameon=False);fig.tight_layout();save(fig,'Figure_05_poisson_thinning_times')
 # 06 and saved diameter/location rejection tables.
 b=np.linspace(.5,4,19);binrows=[]
 for lo,hi in zip(b[:-1],b[1:]):
  m=(diam>=lo)&(diam<hi);count=int(m.sum());nacc=int((m&keep).sum());wall=int((m&(reasons=='WALL_REJECTED')).sum());nf=int((m&(reasons=='WALL_NEARFIELD_REJECTED')).sum())
  binrows.append(dict(D_low_um=lo,D_high_um=hi,source=count,accepted=nacc,rejected=count-nacc,wall=wall,nearfield=nf,rejection_fraction=(count-nacc)/count if count else None))
 write_csv(D/'rejections_by_diameter.csv',binrows)
 loc=[]
 for i in range(len(g['triangles_m'])):
  mask=tri==i;count=int(mask.sum());loc.append(dict(triangle_id=i,source=count,accepted=int((mask&keep).sum()),wall=int((mask&(reasons=='WALL_REJECTED')).sum()),nearfield=int((mask&(reasons=='WALL_NEARFIELD_REJECTED')).sum()),rejection_fraction=float((mask&~keep).sum()/count) if count else None))
 write_csv(D/'rejections_by_inlet_triangle.csv',loc)
 fig,axes=plt.subplots(1,2,figsize=(10.5,4.2));mid=(b[:-1]+b[1:])/2
 axes[0].plot(mid,[r['rejection_fraction'] for r in binrows],'o-',color=RED);axes[0].set(xlabel='Diameter (µm)',ylabel='Rejected source fraction',ylim=(0,1.03),title='Finite-size rejection vs diameter');axes[0].grid(alpha=.16)
 counts=np.array([[r['wall'],r['nearfield']] for r in binrows]);axes[1].bar(mid,counts[:,0],width=np.diff(b)*.85,color=RED,label='Real WALL overlap');axes[1].bar(mid,counts[:,1],bottom=counts[:,0],width=np.diff(b)*.85,color=ORANGE,label='Original handoff rule');axes[1].set(xlabel='Diameter (µm)',ylabel='Rejected proposal count',title='Rejection reasons (not solver failures)');axes[1].legend(frameon=False,fontsize=9);fig.tight_layout();save(fig,'Figure_06_rejection_by_diameter')
 # 07: semantics, with separate empirical entry-size distributions.
 B=read(D/'legacy_b_inlet_500.json');C=read(D/'legacy_c_inlet_500.json');fig=plt.figure(figsize=(12,7));ax=fig.add_axes([.03,.42,.94,.49]);ax.axis('off')
 body=[['Source size','Original full SonoVue','SonoVue D ≤ 4 µm','SonoVue D ≤ 4 µm'],['Source time','Legacy entering clock','Deterministic entering clock','Poisson source clock'],['Size retry','Yes, fixed anchor','If globally infeasible','Never'],['Position retry','No','Yes, fixed feasible size','Never'],['Rejected physical event','Not represented','Not represented','Retained with source_event_id'],['Entering rate','Legacy interpretation','C_MB × Q by definition','C_≤4 × E[Q_acc(D)]'],['Trajectory comparison','Not used','Not used','No outlet quota']]
 table=ax.table(cellText=body,colLabels=['Property','Legacy Method B','Legacy Method C','P9-A.4'],cellLoc='left',loc='center',colWidths=[.21,.24,.27,.28]);table.auto_set_font_size(False);table.set_fontsize(10);table.scale(1,1.85)
 for (i,j),cell in table.get_celld().items():cell.set_edgecolor('#d9e0e6');cell.set_facecolor('#eaf0f5' if i==0 else 'white')
 ax2=fig.add_axes([.09,.10,.82,.25]);ax2.hist([r['result']['diameter_um'] for r in B if r['result']['accepted']],bins=bins,density=True,histtype='step',lw=1.8,label='B accepted (306/500 requested)',color=GRAY);ax2.hist([r['diameter_um'] for r in C],bins=bins,density=True,histtype='step',lw=1.8,label='C accepted (500 requested)',color=BLUE);ax2.hist(diam[keep],bins=bins,density=True,histtype='step',lw=1.8,label='P9-A.4 accepted (19,221/100,000 source)',color=ORANGE);ax2.set(xlabel='Accepted diameter (µm)',ylabel='Density (µm⁻¹)');ax2.legend(frameon=False,fontsize=9,ncol=3);fig.suptitle('Legacy conditioning and physical single-shot thinning have different semantics',fontsize=14,y=.98);save(fig,'Figure_07_legacy_semantic_comparison')
 print('FIGURES_01_TO_07_COMPLETE',flush=True)

def smoke_figure():
 import pyvista as pv
 from mpl_toolkits.mplot3d.art3d import Poly3DCollection
 from matplotlib.lines import Line2D
 reference=ROOT/'particle_3d/reports/network_derived_flow_mb_validation_v1/server_bundle/formal_3D_flow_solver/FEM_SimVascular'
 manifest=read(reference/'frozen_reference/boundary_manifest.json');spec=manifest['boundaries']['WALL'];wall=pv.read(reference/spec['path']);triangles=wall.points[wall.faces.reshape(-1,4)[:,1:]]*1e6
 rows=read(D/'smoke30_metrics.json');colors={'O1':BLUE,'O2':GREEN,'O3':ORANGE,'STATIONARY':RED,'SOLVER_FAIL':'black'}
 fig=plt.figure(figsize=(9,7));ax=fig.add_subplot(111,projection='3d');ax.add_collection3d(Poly3DCollection(triangles,facecolor='#cad3db',edgecolor='none',alpha=.035,rasterized=True))
 for r in rows:
  a=np.load(R/'outputs/smoke30/trajectories'/f"mb_{r['bubble_id']:06d}.npz")['samples'];xyz=a[:,1:4]*1e6;color=colors.get(r['outlet'] or r['status'],GRAY);ax.plot(*xyz.T,color=color,lw=1.4,alpha=.9);ax.scatter(*xyz[0],color=color,s=12)
 ax.set(xlabel='x (µm)',ylabel='y (µm)',zlabel='z (µm)',title='First 30 accepted P9-A.4 births | NEW Network-H0 flow');ax.view_init(elev=22,azim=-65);limits=np.ptp(wall.points,axis=0);ax.set_box_aspect(limits);ax.legend(handles=[Line2D([0],[0],color=v,lw=2,label=k) for k,v in colors.items() if any((r['outlet'] or r['status'])==k for r in rows)],frameon=False,loc='upper left');fig.tight_layout();save(fig,'Figure_08_smoke30_trajectories');print('FIGURE_08_COMPLETE')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['inlet','smoke']);a=p.parse_args();inlet_figures() if a.mode=='inlet' else smoke_figure()
