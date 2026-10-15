"""Four review figure groups from frozen full-precision outputs; no solver/integration."""
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
import time, resource, html

start=time.perf_counter(); F.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.titlesize':16,'axes.labelsize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':11,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.color':'#b4bdc5','grid.alpha':.65,'grid.linewidth':.6,'savefig.facecolor':'white','svg.fonttype':'none'})
blue='#2b7099'; orange='#d88927'; red='#b84c50'; gray='#9aa7b1'; purple='#9067a8'
npz=np.load(O/'LAMBDA_AND_ENERGY_DENSITY.npz'); raw=np.load(O/'source_inputs/reference_topology.npz')
X=npz['X_m']; b=raw['bonds']; roi=npz['roi']; active=npz['0_active_ROI']; cross=npz['0_cross_plane_bonds']; lam=npz['0_lambda']; n=len(X)
res=json.loads((O/'DISCRETE_SURFACE_ENERGY_RESULTS.json').read_text()); main,control=res['cases']; theory=json.loads((O/'THEORY_FORMULA_AUDIT.json').read_text()); topo=json.loads((O/'CROSS_PLANE_REFERENCE_TOPOLOGY_AUDIT.json').read_text())['cases'][0]
files=[]; plotdata={}
def save(fig,name):
    for ext in ['png','svg']:
        p=F/(name+'.'+ext); fig.savefig(p,dpi=210,bbox_inches='tight'); files.append(p.name)
    plt.close(fig)
def axes3(a,aspect=(1,1,1)):
    a.set(xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)');a.view_init(elev=17,azim=-64);a.set_proj_type('ortho');a.set_box_aspect(aspect)
    a.xaxis.labelpad=15; a.yaxis.labelpad=12; a.zaxis.labelpad=11
    for aa in [a.xaxis,a.yaxis,a.zaxis]:
        aa.pane.set_facecolor((.97,.98,.99,.6)); aa._axinfo['grid']['color']=(.55,.6,.65,.55)
def plane(a,xlim,ylim,z=0,alpha=.15):
    a.add_collection3d(Poly3DCollection([[(xlim[0],ylim[0],z),(xlim[1],ylim[0],z),(xlim[1],ylim[1],z),(xlim[0],ylim[1],z)]],facecolors=blue,edgecolors=blue,alpha=alpha))

# Fig1: true point set (all 16000 points) and one actual complete neighborhood.
fig=plt.figure(figsize=(15,7.2));ax=fig.add_subplot(121,projection='3d',computed_zorder=False);ay=fig.add_subplot(122,projection='3d')
for mask,color,size,alpha in [(~roi,gray,2,.09),(roi&~active,gray,3,.20),(active&(X[:,2]<0),blue,6,.85),(active&(X[:,2]>0),red,6,.85)]:
    ax.scatter(*(X[mask]*1e3).T,s=size,c=color,alpha=alpha,depthshade=False,rasterized=True)
plane(ax,[-2.5,2.5],[-.2,.2],alpha=.08)
roioutline=np.array([[-.75,-.2,0],[.75,-.2,0],[.75,.2,0],[-.75,.2,0],[-.75,-.2,0]])
ax.plot(*roioutline.T,c=orange,lw=2.8,zorder=10)
for yy in [-.2,.2]:ax.plot([-2.5,2.5],[yy,yy],[0,0],c=purple,ls='--',lw=1.2,zorder=9)
sample=np.where(cross&(abs(X[b[:,0],0])<.12e-3)&(abs(X[b[:,1],0])<.12e-3)&(abs(X[b[:,0],1])<.076e-3)&(abs(X[b[:,1],1])<.076e-3))[0][::4]
ax.add_collection3d(Line3DCollection(X[b[sample]]*1e3,colors=orange,linewidths=.5,alpha=.5))
axes3(ax,(5,1.3,1.5));ax.set(xlim=(-2.6,2.6),ylim=(-.23,.23),zlim=(-.53,.53),xticks=[-2.5,0,2.5],yticks=[-.2,0,.2],zticks=[-.5,0,.5]);ax.set_title('(a) Pristine specimen and frozen ROI',pad=15)
ax.legend(handles=[Line2D([],[],marker='o',ls='',color=blue,label='Lower active nodes'),Line2D([],[],marker='o',ls='',color=red,label='Upper active nodes'),Line2D([],[],color=orange,lw=2,label='ROI and displayed cross-plane bonds'),Line2D([],[],color=purple,ls='--',label='Periodic Y cell boundaries')],loc='lower left',bbox_to_anchor=(-.04,-.11),frameon=False,ncol=2)
center=np.array([.000025,.000175,.000025]); i=int(np.argmin(np.linalg.norm(X-center,axis=1))); bix=np.where((b==i).any(1))[0]; js=np.where(b[bix,0]==i,b[bix,1],b[bix,0]); neigh=X[js].copy(); shifted=-np.rint((neigh[:,1]-X[i,1])/PERIOD); neigh[:,1]+=shifted*PERIOD
segments=np.stack([np.broadcast_to(X[i],neigh.shape),neigh],axis=1)*1e3; localcross=cross[bix]
ay.add_collection3d(Line3DCollection(segments[~localcross],colors=gray,linewidths=.55,alpha=.5));ay.add_collection3d(Line3DCollection(segments[localcross],colors=orange,linewidths=1.15,alpha=.9))
ay.scatter(*(neigh*1e3).T,c=np.where(neigh[:,2]<0,blue,red),s=16,depthshade=False)
ghost=shifted!=0;ay.scatter(*(neigh[ghost]*1e3).T,s=42,facecolors='none',edgecolors=purple,linewidths=1)
ay.scatter(*(X[i]*1e3),s=130,c='black',marker='*',depthshade=False)
for zz in [-.125,0,.175]:
    rr=np.sqrt(max(0,DELTA**2-(zz*1e-3-X[i,2])**2))*1e3; tt=np.linspace(0,2*np.pi,90)
    ay.plot(X[i,0]*1e3+rr*np.cos(tt),X[i,1]*1e3+rr*np.sin(tt),np.full(len(tt),zz),c=gray,ls=':',lw=.7)
plane(ay,[-.14,.19],[.01,.34],alpha=.10)
ay.plot([-.14,.19],[.2,.2],[0,0],c=purple,ls='--',lw=1.4)
axes3(ay);ay.set_title(f'(b) Actual full neighborhood of atom {i+1}',pad=15);ay.set(xlim=(-.14,.19),ylim=(.01,.34),zlim=(-.14,.19));ay.locator_params(nbins=4)
fig.text(.58,.14,f'{len(js)} neighbors; {int(localcross.sum())} cross-plane\nOpen rings: periodic images of existing IDs',va='top',fontsize=12)
fig.suptitle('True 3D reference topology — no continuous clot shell',fontsize=20,y=.98)
fig.text(.5,.055,'16,000 nodes  |  909,392 unique bonds  |  55,408 cross-plane bonds  |  ROI = 1.5 × 0.4 mm²',ha='center',fontsize=13)
fig.subplots_adjust(left=.035,right=.965,bottom=.22,top=.88,wspace=.15)
save(fig,'Fig1_reference_plane_and_bonds');plotdata['Fig1']={'all_material_points_shown':n,'displayed_cross_bond_reference_IDs':(sample+1).tolist(),'local_center_atom_ID':i+1,'local_reference_bond_IDs':(bix+1).tolist(),'periodic_images_are_not_extra_physical_nodes':True,'visual_box_aspect_not_equal_physical_scale':True}

# Fig2: actual one-Y-layer maps and a full-ROI layer profile (all 8 Y layers).
fig,axs=plt.subplots(2,3,figsize=(16,8.4),layout='constrained')
ys=np.unique(X[:,1]); yi=ys[4]; keep=(X[:,1]==yi)&(abs(X[:,0])<.001)&(abs(X[:,2])<.00023)
spec=[(lam,r'(a) $\lambda=m_{cross}/m_{full}$','fraction','viridis'),(npz['0_m_cross_m4']*1e16,r'(b) Cross-plane weighted volume',r'$m_{cross}$ ($10^{-16}$ m$^4$)','viridis'),(npz['0_W_K_J_m3'],r'(c) Volumetric energy density $W_K$','J m$^{-3}$','magma'),(npz['0_W_G_J_m3'],r'(d) Deviatoric energy density $W_G$','J m$^{-3}$','magma'),(npz['0_W_ideal_J_m3'],r'(e) Total ideal energy density','J m$^{-3}$','magma')]
for a,(arr,title,unit,cmap) in zip(axs.ravel(),spec):
    sc=a.scatter(X[keep,0]*1e3,X[keep,2]*1e3,c=arr[keep],s=9,marker='o',cmap=cmap,vmin=0,vmax=float(arr[keep].max()),linewidths=0)
    for xx in [-.75,.75]:a.axvline(xx,c=orange,ls='--',lw=1.3)
    a.axhline(0,c='#67727b',ls=':',lw=.9);a.set(xlabel='X (mm)',ylabel='Z (mm)',title=title,xlim=(-1.02,1.02),ylim=(-.24,.24));cb=fig.colorbar(sc,ax=a,shrink=.82,pad=.02);cb.set_label(unit,fontsize=12)
a=axs[1,2];profiles=[p for p in json.loads((O/'LAYER_PROFILES.json').read_text()) if p['plane_index']==0 and p['z_m']>0 and p['distance_m']<DELTA]
dist=np.array([p['distance_m'] for p in profiles])*1e6; val=[p['lambda_mean'] for p in profiles]; qq=np.linspace(0,1,200)
a.plot(qq*DELTA*1e6,.5-2*qq/3+qq**4/6,c=gray,lw=2,label='Continuous cap')
a.plot(dist,val,'o-',c=blue,ms=8,lw=1.5,label='Actual discrete layers');a.set(xlabel='Distance from plane (µm)',ylabel=r'$\lambda$',title='(f) Full-ROI layer averages',xlim=(0,DELTA*1e6),ylim=(0,.52));a.legend(frameon=False)
fig.suptitle('Main plane z = 0: actual node values and native partial-volume weights',fontsize=19)
fig.text(.02,-.02,'Maps: Y = +25 µm layer. Dashed lines delimit ROI; all 8 periodic Y layers enter the energy sum. No smoothing of node values.',fontsize=12)
save(fig,'Fig2_lambda_and_local_energy')

# Fig3: zero-based absolute scale, contributions, direct quadrature discrepancy.
fig,axs=plt.subplots(1,3,figsize=(16,5.6),layout='constrained');a=axs[0]
kk=np.array([theory['K_contribution_J_m2'],main['K_contribution_J_m2'],control['K_contribution_J_m2']]);gg=np.array([theory['G_contribution_J_m2'],main['G_contribution_J_m2'],control['G_contribution_J_m2']]); xx=np.arange(3)
a.bar(xx,kk,color=blue,label='K contribution',width=.62);a.bar(xx,gg,bottom=kk,color=orange,label='G contribution',width=.62)
for j,v in enumerate(kk+gg):a.text(j,v+.0018,f'{v:.8f}',ha='center',fontsize=11)
a.set(xticks=xx,xticklabels=['Continuous\ntheory','Discrete\nz = 0','Discrete\nz = 50 µm'],ylabel=r'Nominal energy (J m$^{-2}$)',ylim=(0,.09),title='(a) Nominal surface energy');a.legend(loc='upper right',frameon=False)
a=axs[1];a.bar(xx,100*kk/(kk+gg),color=blue,width=.62);a.bar(xx,100*gg/(kk+gg),bottom=100*kk/(kk+gg),color=orange,width=.62)
for j in xx:
    a.text(j,45,f'{100*kk[j]/(kk[j]+gg[j]):.3f}% K',ha='center',color='white',fontsize=12);a.text(j,95.5,f'{100*gg[j]/(kk[j]+gg[j]):.3f}%',ha='center',fontsize=11)
a.set(xticks=xx,xticklabels=['Theory','Main','Control'],ylabel='Fraction of nominal energy (%)',ylim=(0,100),title='(b) K / G decomposition')
a=axs[2];names=['K term','G term','Total'];diff=np.array([main['K_term_difference_from_theory_J_m2'],main['G_term_difference_from_theory_J_m2'],main['difference_from_theory_J_m2']])*1e3
a.bar(names,diff,color=[blue,orange,gray],width=.6);a.axhline(0,c='#333333',lw=1)
for j,v in enumerate(diff):a.text(j,v-.075,f'{v:.5f}',ha='center',va='top',fontsize=11)
a.set(ylabel=r'Discrete − theory (mJ m$^{-2}$)',ylim=(-2.45,.8),title='(c) Actual finite-quadrature difference')
a.text(.04,.97,f'Main / theory = {main["G0_over_theory"]:.6f}\nRelative difference = {100*main["relative_difference_from_theory"]:.4f}%',transform=a.transAxes,va='top',fontsize=12)
fig.suptitle('Nominal energy identified; actual crack-propagation resistance remains unmeasured',fontsize=18)
fig.text(.02,-.03,'Two-form relative error: 3.16×10⁻¹⁶. One translated-gap control differs by 2.78×10⁻¹⁷ J/m². Neither result is a calibrated Gc.',fontsize=12)
save(fig,'Fig3_theory_and_discrete_energy')

# Fig4: traceable method, real taper shells, honest counting and limits.
fig,axs=plt.subplots(2,2,figsize=(15,9.5),layout='constrained')
a=axs[0,0];a.axis('off');a.set_title('(a) Paper → source → ideal discrete state',loc='left')
method=[('Paper, intact 3D family',r'$W=\frac{K}{2}\theta^2+\frac{15G}{2m}\sum_j\omega(e-\theta r/3)^2 V_j h$'),('Influence normalization',r'$\omega_{paper}=\delta/r,\quad \omega_{source}=1/r,\quad m_{paper}=\delta m_{source}$'),('Same reference family in both forms',r'$m_i=\sum_j r_0 V_j h,\quad e_{cross}=s_{00}r_0,\quad e_{other}=0$')]
for y,(title,formula) in zip([.85,.54,.23],method):
    a.text(.02,y,title,transform=a.transAxes,fontsize=13,fontweight='bold');a.text(.02,y-.115,formula,transform=a.transAxes,fontsize=14)
a=axs[0,1];shells=json.loads((O/'TAPER_SHELL_AUDIT.json').read_text())['shells']; rr=[v['r_over_dx_mean'] for v in shells]; hh=[v['taper_mean'] for v in shells]
xx=np.linspace(0,3.01,301);hcurve=np.where(xx<2.51,1,3.51-xx);a.plot(xx,hcurve,c=blue,lw=2,label='Native taper rule');a.scatter(rr,hh,c=orange,s=[v['neighbors']*6 for v in shells],edgecolors='white',zorder=4,label='Actual neighbor shells')
for v in shells:a.annotate(str(v['neighbors']),(v['r_over_dx_mean'],v['taper_mean']),xytext=(0,10),textcoords='offset points',ha='center',fontsize=10)
a.set(xlabel=r'Reference length $r_0/\Delta x$',ylabel='Partial-volume weight h',xlim=(0,3.16),ylim=(0,1.25),title='(b) Actual taper: 122 neighbors / interior node');a.axvline(3.01,c=gray,ls='--');a.legend(frameon=False,loc='lower left');a.text(.05,.25,'Labels = shell multiplicity; size ∝ multiplicity',transform=a.transAxes,fontsize=11)
a=axs[1,0];a.axis('off');a.set_title('(c) Periodic identity and area accounting',loc='left')
items=['909,392 unique reference interactions','131,196 cross the periodic Y seam','55,408 cross z = 0; 7,324 also cross Y seam','ROI: 720 lower + 720 upper active nodes','Full families retained at ROI edges',r'$A=(30\Delta x)\,B=0.6\ \mathrm{mm}^2$',r'$G_0=\sum_{i\in ROI}W_i V_i/A$  (no extra ×2)']
for y,t in zip(np.linspace(.88,.08,len(items)),items):a.text(.03,y,t,transform=a.transAxes,fontsize=13)
a=axs[1,1];a.axis('off');a.set_title('(d) Confirmed scope and remaining uncertainty',loc='left')
lines=[('CONFIRMED',blue),('Full pristine topology; two energy forms agree',None),('One interior translation gives the same result',None),('LIMITED',orange),('δ / Δx = 3.01; no grid or phase convergence',None),('Ideal extensions need not be globally realizable',None),('NOT MEASURED',red),('Actual Gc, natural crack area, dissipated work',None),('Next: controlled Mode-I separation (not run)',None)]
for y,(t,col) in zip(np.linspace(.90,.05,len(lines)),lines):a.text(.03,y,t,transform=a.transAxes,color=col or '#263746',fontweight='bold' if col else 'normal',fontsize=12.5)
fig.suptitle('Method audit — no new dynamics, equilibrium solve, parameter scan, or native modification',fontsize=18)
save(fig,'Fig4_method_and_uncertainty')
dump(O/'FIGURE_DATA_PROVENANCE.json',{'figures':files,'all_computed_inputs_from':str(O),'source_data_sha256':sha(O/'LAMBDA_AND_ENERGY_DENSITY.npz'),'details':plotdata,'new_energy_evaluations':0})
intro=f'<p>真实 50 µm 点阵、s00=0.15：<strong>G0_discrete = {main["G0_discrete_J_m2"]:.12f} J/m²</strong>；连续理论 {theory["G0_theory_J_m2"]:.12f} J/m²。仅为理想标称能量，GC_EFFECTIVE 尚未测量。</p>'
cards=''.join(f'<section><h2>{html.escape(name)}</h2><a href="{name}.svg">SVG</a> · <a href="{name}.png">PNG</a><img src="{name}.png" alt="{name}"></section>' for name in [Path(f).stem for f in files if f.endswith('.png')])
page='<!doctype html><html lang="zh"><meta charset="utf-8"><title>LPS 3D 裂面标称能量</title><style>body{font:17px system-ui;max-width:1480px;margin:32px auto;padding:0 24px;color:#263746}img{width:100%;height:auto}section{margin:32px 0;border-top:1px solid #cad3da}h1{font-size:28px}h2{font-size:20px}a{color:#22658a}</style><h1>LPS 三维裂面标称能量识别</h1>'+intro+f'<p><a href="../../../verification/{C["task"]}/{C["run"]}/REVIEW_LPS_3D_FRACTURE_SURFACE_ENERGY_IDENTIFICATION_V1.md">中文审查报告与28项回答</a> · <a href="../../../verification/{C["task"]}/{C["run"]}/DISCRETE_SURFACE_ENERGY_RESULTS.json">完整精度结果</a></p>'+cards+'</html>'
(F/'OPEN_RESULTS.html').write_text(page)
(O/'OPEN_RESULTS.html').write_text(f'<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="0;url=../../../results/{C["task"]}/{C["run"]}/OPEN_RESULTS.html"><a href="../../../results/{C["task"]}/{C["run"]}/OPEN_RESULTS.html">打开四组审查图</a>')
previous=json.loads((O/'VISUALIZATION_RESOURCES.json').read_text()) if (O/'VISUALIZATION_RESOURCES.json').exists() else None
dump(O/'VISUALIZATION_RESOURCES.json',{'wall_s':time.perf_counter()-start,'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'files':files,'GPU_used':False,'native_or_formal_integral_calls':0,'previous_render':previous,'revision_reason':'Readability of point markers, ROI outline and text placement; numerical data unchanged' if previous else None})
print(json.dumps({'figures':len(files),'wall_s':time.perf_counter()-start}))
