"""Four figures from actual native topology and accepted equilibrium data."""
from analysis.lps_pvce_v1.core import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
RED='#a93040';BLUE='#2375a3';TEAL='#208e84';GOLD='#c58a2c';GRAY='#777777'
COLORS=[TEAL,BLUE,GOLD]
OUT=ROOT/'results/lps_pvce_energy_release_feasibility_v1/20261011T132103Z'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':14,'axes.titlesize':16,'axes.labelsize':14,'legend.fontsize':12,'axes.grid':True,'grid.color':'#777777','grid.alpha':.55,'svg.fonttype':'none'})

def view(ax,title):
 ax.view_init(elev=12,azim=-58);ax.set(xlabel='x (µm)',ylabel='y (µm)',zlabel='z (µm)',title=title,xlim=[-330,330],ylim=[-330,330],zlim=[-520,520],xticks=[-300,0,300],yticks=[-300,0,300],zticks=[-500,0,500]);ax.set_box_aspect([.66,.66,1.04]);ax.tick_params(labelsize=11,pad=1)
 for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.set_facecolor((.97,.97,.97,1));axis._axinfo['grid'].update(color=(.42,.42,.42,.65),linewidth=.7)
def savefig(fig,name,svg=True):
 fig.savefig(OUT/(name+'.png'),dpi=240,bbox_inches='tight',facecolor='white')
 if svg:fig.savefig(OUT/(name+'.svg'),bbox_inches='tight',facecolor='white')
 plt.close(fig)
def main():
 OUT.mkdir(parents=True,exist_ok=True);top=json.loads((BASE/'PVCE_TOPOLOGY_AUDIT.json').read_text());R=json.loads((BASE/'PVCE_ENERGY_RELEASE_RESULTS.json').read_text());F=json.loads((BASE/'PVCE_DESIGN_FREEZE.json').read_text());refs=[dict(np.load(BASE/f'A{j}/reference.npz')) for j in range(3)];states=[dict(np.load(BASE/f'A{j}/display_fields.npz')) for j in range(3)]
 # True deleted bonds and full point sets. Separate close views reveal selection.
 fig=plt.figure(figsize=(18,11));fig.subplots_adjust(top=.89,bottom=.13,hspace=.25,wspace=.20)
 for j,z in enumerate(refs):
  X=z['X']*1e6;b=z['bonds'];deleted=z['alive']==0;inherited=refs[j-1]['alive']==0 if j else np.zeros(len(b),bool);new=deleted&~inherited;ax=fig.add_subplot(2,3,j+1,projection='3d');mask=z['mask'][:,2]
  ax.scatter(*X[~mask].T,c=RED,s=5,alpha=.30,depthshade=False);ax.scatter(*X[mask].T,c=BLUE,s=10,alpha=.75,depthshade=False)
  for m,col in [(inherited,GOLD),(new,TEAL if j==0 else BLUE)]:ax.add_collection3d(Line3DCollection(X[b[m]],colors=col,linewidths=.4,alpha=.22))
  tip=top['records'][j]['tip_m']*1e6;ax.plot([tip,tip],[-300,300],[0,0],c='black',lw=2);view(ax,f'A{j}: a = {200+50*j} µm\n{int(deleted.sum()):,} reference bonds cut')
  ax=fig.add_subplot(2,3,j+4);slab=(abs(X[b[:,0],1])<=26)&(abs(X[b[:,1],1])<=26)&(abs(X[b[:,0],2])<=150)&(abs(X[b[:,1],2])<=150)
  from matplotlib.collections import LineCollection
  for m,col in [(~deleted,'#c9c9c9'),(inherited,GOLD),(new,TEAL if j==0 else BLUE)]:ax.add_collection(LineCollection(X[b[m&slab]][:,:,[0,2]],colors=col,linewidths=.7,alpha=.70))
  sel=(abs(X[:,1])<=26)&(abs(X[:,2])<=150);ax.scatter(X[sel,0],X[sel,2],s=14,c=RED,zorder=4);ax.plot([-300,tip],[0,0],c='black',lw=2);ax.axvline(tip,c='black',ls=':',lw=1.1);ax.set(xlim=[-305,305],ylim=[-155,155],xlabel='Reference x (µm)',ylabel='Reference z (µm)',title=f'Central y slab | '+('initial cut' if j==0 else '+766 actual virtual cuts'));ax.set_aspect('equal')
 fig.legend(handles=[Line2D([],[],c=GOLD,label='Inherited cut bonds'),Line2D([],[],c=TEAL,label='A0 prescribed pre-cut'),Line2D([],[],c=BLUE,label='Current added cuts'),Line2D([],[],c='black',label='Declared crack / front')],loc='lower center',bbox_to_anchor=(.5,.01),ncol=4)
 fig.suptitle('1 | Actual reference lattice and nested virtual cuts\nB = 600 µm through y; Δa = 50 µm; ΔA = BΔa = 3.0 × 10⁻⁸ m² (one advance area)',fontsize=19);savefig(fig,'Fig1_virtual_crack_topologies')
 # Actual optimizer histories. No kinetic trace for static minimization.
 fig,axs=plt.subplots(2,2,figsize=(15,11));fig.subplots_adjust(top=.89,bottom=.12,hspace=.31,wspace=.27)
 for j in range(3):
  h=np.genfromtxt(BASE/f'A{j}/history.csv',delimiter=',',names=True);c=COLORS[j];it=h['iteration'];axs[0,0].semilogy(it,h['relative_force'],c=c,label=f'A{j}');axs[0,1].plot(it,h['U_J']*1e12,c=c,label=f'A{j}');axs[1,0].semilogy(it,h['max_force_N'],c=c);axs[1,1].semilogy(it[9:],h['energy_span_relative'][9:],c=c)
 axs[0,0].axhline(1e-6,c='black',ls='--',label='Frozen force gate');axs[1,0].axhline(1e-13,c='black',ls='--');axs[1,1].axhline(1e-8,c='black',ls='--');axs[0,0].legend();axs[0,1].legend()
 for ax in axs.ravel():ax.set_xlabel('Accepted L-BFGS iteration')
 axs[0,0].set(ylabel='||F free||₂ / ||F controlled||₂',title='Actual residual; all three accepted');axs[0,1].set(ylabel='Independent elastic energy U (pJ)',title='Fixed opening = 2 µm');axs[1,0].set(ylabel='Maximum free force component (N)',title='Absolute force gate');axs[1,1].set(ylabel='Last 10 states: energy span / U',title='Plateau is checked alongside force')
 fig.suptitle('2 | Three bounded static equilibria',fontsize=20);fig.text(.5,.025,'K/U: not applicable to a static optimizer. No artificial kinetic curve or pseudo physical time.\nNative checks use zero velocities; acceptance also requires stable coordinates and balanced reactions.',ha='center',fontsize=13);savefig(fig,'Fig2_equilibrium_convergence')
 # Errorbars display the empirical solver-resolution scale, not material uncertainty.
 fig,axs=plt.subplots(2,2,figsize=(15,10));fig.subplots_adjust(top=.88,bottom=.14,hspace=.38,wspace=.30);U=np.array(R['U_J']);err=np.array(R['empirical_U_resolution_J']);aa=np.array([200,250,300]);axs[0,0].errorbar(aa,U*1e12,yerr=err*1e12,fmt='o',c=BLUE,ms=9,capsize=7);axs[0,0].set(xlabel='Actual crack length a (µm)',ylabel='U (pJ)',title='Three actual equilibria; error bars to scale',xticks=aa)
 for j in range(3):axs[0,0].annotate(f'A{j}',(aa[j],U[j]*1e12),xytext=(10,4),textcoords='offset points')
 du=-np.diff(U);de=err[:-1]+err[1:];axs[0,1].errorbar([0,1],du*1e12,yerr=de*1e12,fmt='o',c=TEAL,ms=9,capsize=7);axs[0,1].set(xticks=[0,1],xticklabels=['U0 − U1','U1 − U2'],ylabel='Elastic energy decrease (pJ)',title='Both differences exceed solver noise')
 vals=np.array([R['G_two_J_m2'],R['G_three_J_m2']])*1e6;ee=np.array([R['G_two_empirical_resolution_J_m2'],R['G_three_empirical_resolution_J_m2']])*1e6;axs[1,0].errorbar([0,1],vals,yerr=ee,fmt='o',c=RED,ms=9,capsize=7);axs[1,0].set(xticks=[0,1],xticklabels=['Two point','Three point'],ylabel='G drive (µJ/m²)',title=f'One local estimate at a0 = 200 µm\nDifference = {R["relative_difference_to_G_three"]*100:.4f}%')
 axs[1,1].bar(['A0','A1','A2'],err*1e18,color=COLORS);axs[1,1].set(ylabel='Empirical energy resolution (aJ)',title='Residual-work + last-10-state drift scale')
 fig.suptitle('3 | Specimen-dependent PVCE energy release: G drive ≠ material Gc',fontsize=19);fig.text(.5,.025,'Error bars estimate solver resolution only; they are not rigorous bounds and exclude finite-size / horizon / crack-increment error.\nNo load, dx or horizon scan. These three states do not provide a complete G(a) curve.',ha='center',fontsize=12);savefig(fig,'Fig3_PVCE_energy_release')
 # Shared scales; accepted coordinates, deformation scale1; point-cloud only.
 fig=plt.figure(figsize=(16,12));fig.subplots_adjust(top=.90,bottom=.1,wspace=.2,hspace=.12);vmax=max(np.linalg.norm(states[j]['u_m'],axis=1).max()*1e6 for j in [0,2]);emax=max(states[j]['energy_density_J_m3'].max() for j in [0,2]);axes=[]
 for col,j in enumerate([0,2]):
  z=states[j];x=z['x_m']*1e6;m=z['mask'][:,2];tip=top['records'][j]['tip_m']*1e6
  for row in range(2):
   ax=fig.add_subplot(2,2,row*2+col+1,projection='3d');axes.append(ax);val=np.linalg.norm(z['u_m'],axis=1)*1e6 if row==0 else z['energy_density_J_m3'];sc=ax.scatter(*x.T,c=val,cmap='viridis' if row==0 else 'magma',vmin=0,vmax=vmax if row==0 else emax,s=16,alpha=.95,depthshade=False);ax.scatter(*x[m].T,s=23,facecolors='none',edgecolors=BLUE,linewidths=.7,depthshade=False)
   # Front is a declared reference line only, never rendered as a fitted surface.
   ax.plot([tip,tip],[-300,300],[0,0],c='limegreen',lw=2,ls='--');view(ax,f'A{j} | '+('Displacement magnitude' if row==0 else 'Independent LPS energy density'))
   if col==1:
    cb=fig.colorbar(sc,ax=[axes[row*2] if len(axes)>row*2 else ax,ax],fraction=.022,pad=.03,shrink=.65);cb.set_label('|u| (µm)' if row==0 else 'W (J/m³)')
 fig.suptitle('4 | Actual equilibria: deformation scale = 1\nPoint-cloud fields; shared color ranges for A0 and A2',fontsize=20);fig.text(.5,.025,'Blue outlines: end-plane z-controlled nodes (x/y free). Other nodes: all DOFs free.\nGreen dashed lines mark the declared reference crack fronts; no continuous crack shell is synthesized.',ha='center',fontsize=13);savefig(fig,'Fig4_equilibrium_fields')
 html='''<!doctype html><meta charset="utf-8"><title>LPS PVCE 实际结果</title><style>body{max-width:1400px;margin:auto;font:18px sans-serif;line-height:1.6;background:#fafafa;color:#222}img{width:100%;background:white}a{color:#16749a}section{margin:36px 0}</style><h1>LPS PVCE：三个平衡态的试样相关 G_drive</h1><p>固定位移总开度 2 µm；a=200/250/300 µm。三个状态均通过预注册平衡检查。G_two=6.18596×10⁻⁶ J/m²；G_three=6.17031×10⁻⁶ J/m²；差异0.2536%。结果不是材料固有 Gc，未进行离散收敛认证。</p>'''
 for p in sorted(OUT.glob('Fig*.png')):html+=f'<section><a href="{p.name}"><img src="{p.name}"></a><p><a href="{p.stem}.svg">SVG</a> · {p.stem}</p></section>'
 (OUT/'OPEN_RESULTS.html').write_text(html)
 # Verification index uses a stable relative link to actual results.
 rel='../../../results/lps_pvce_energy_release_feasibility_v1/20261011T132103Z/'
 (BASE/'OPEN_RESULTS.html').write_text(html.replace('src="Fig','src="'+rel+'Fig').replace('href="Fig','href="'+rel+'Fig'))
 print('generated',len(list(OUT.glob('Fig*.png'))),'PNG and SVG figures')
if __name__=='__main__':main()
