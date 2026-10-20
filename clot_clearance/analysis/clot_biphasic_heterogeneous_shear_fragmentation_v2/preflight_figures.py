from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import pandas as pd

plt.rcParams.update({'font.size':12,'axes.titlesize':14,'axes.labelsize':12,'pdf.fonttype':42})
colors=np.array(plt.get_cmap('tab10').colors)
def save(fig,name):
 fig.savefig(F/'figures'/f'{name}.png',dpi=220,bbox_inches='tight',facecolor='white')
 fig.savefig(F/'figures'/f'{name}.pdf',bbox_inches='tight',facecolor='white');plt.close(fig)
def main():
 ref=load_npz(O/'provenance/reference.npz');field=load_npz(O/'BIPHASIC_FIELD.npz');audit=json.loads((O/'PREFLIGHT_STRUCTURE_AUDIT.json').read_text())
 X=ref['X']*1e3;n=int(ref['clot_count']);dom=field['strong_domain_id'];df=pd.read_csv(O/'STRONG_DOMAINS.csv')
 fig=plt.figure(figsize=(17,8))
 for j,az in enumerate([-58,125]):
  ax=fig.add_subplot(1,2,j+1,projection='3d',computed_zorder=False);weak=np.flatnonzero(dom[:n]==0)
  ax.scatter(*X[weak].T,s=4,c='#adb5bd',alpha=.16,depthshade=False,rasterized=False)
  ax.scatter(*X[n:].T,s=3,c='#43474c',alpha=.22,depthshade=False)
  for c in range(1,11):
   p=X[dom==c];ax.scatter(*p.T,s=8,c=[colors[c-1]],alpha=.82,depthshade=False)
   center=p.mean(0);ax.text(*center,str(c),fontsize=14,fontweight='bold',zorder=100,bbox=dict(facecolor='white',alpha=.8,pad=1,edgecolor='none'))
  b=ref['bonds'];near=np.flatnonzero((ref['r0']<DX*1.01)&~ref['initial_cut']);chosen=near[::8]
  ax.add_collection3d(Line3DCollection(X[b[chosen]],colors='#6d747d',alpha=.10,linewidths=.35))
  ax.set(xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)',xlim=(-.67,.67),ylim=(-.67,.67),zlim=(-.5,.35),title=f'Common reference specimen · view {j+1}')
  ax.set_box_aspect((1.34,1.34,.85));ax.view_init(17,az);ax.tick_params(labelsize=11)
 fig.suptitle(f'10 separated strong domains in one connected weak matrix\n{n:,} material nodes + {len(X)-n:,} anchors | strong volume {audit["strong_volume_fraction"]:.3%} | dx = 50 µm, horizon = 150.5 µm',fontsize=17,y=.98)
 handles=[Line2D([],[],marker='o',ls='',color='#adb5bd',label='Weak matrix'),Line2D([],[],marker='o',ls='',color='#43474c',label='Fixed anchors')]+[Line2D([],[],marker='o',ls='',color=colors[i],label=f'Strong {i+1}') for i in range(10)]
 fig.legend(handles=handles,loc='lower center',ncol=6,frameon=False);fig.subplots_adjust(top=.82,bottom=.14,wspace=.08);save(fig,'FIG_01_BIPHASIC_STRUCTURE_3D')
 grid=field['grid_domain_id'];cmap=ListedColormap(['#e2e5e9']+list(colors));norm=BoundaryNorm(np.arange(-.5,11),11)
 fig,axs=plt.subplots(2,3,figsize=(16,10),constrained_layout=True)
 specs=[(2,3,'XY','X','Y'),(1,4,'XZ','X','Z'),(0,12,'YZ','Y','Z'),(2,9,'XY','X','Y'),(1,12,'XZ','X','Z'),(0,20,'YZ','Y','Z')]
 coords=[np.unique(X[:n,k]) for k in range(3)]
 for ax,(axis,index,plane,xa,ya) in zip(axs.flat,specs):
  data=np.take(grid,index,axis=axis);other=[q for q in range(3) if q!=axis];xx,yy=np.meshgrid(coords[other[0]],coords[other[1]],indexing='ij')
  ax.pcolormesh(xx,yy,data,cmap=cmap,norm=norm,shading='nearest',edgecolors='white',linewidth=.15)
  for c in np.unique(data):
   if c:
    ax.contour(xx,yy,(data==c).astype(float),levels=[.5],colors=[colors[c-1]],linewidths=1.4)
    q=np.argwhere(data==c).mean(0);ax.text(np.interp(q[0],np.arange(len(coords[other[0]])),coords[other[0]]),np.interp(q[1],np.arange(len(coords[other[1]])),coords[other[1]]),str(c),ha='center',va='center',fontweight='bold',bbox=dict(facecolor='white',alpha=.85,edgecolor='none',pad=1))
  ax.set(xlabel=f'{xa} (mm)',ylabel=f'{ya} (mm)',title=f'{plane}: {"XYZ"[axis]} = {coords[axis][index]:+.2f} mm',aspect='equal')
 fig.suptitle('Actual binary phase slices · numbered boundaries are material-property regions\nAll original cross-phase bonds retained except the common, phase-independent edge notch',fontsize=16);save(fig,'FIG_02_BIPHASIC_STRUCTURE_SLICES')
 fig,axs=plt.subplots(2,2,figsize=(14,9),constrained_layout=True);ids=df.domain_id
 axs[0,0].bar(ids,df.nodes,color=colors);axs[0,0].set(ylabel='Material nodes',title='Actual region sizes')
 axs[0,1].bar(ids,df.equivalent_volume_diameter_m*1e6,color=colors);axs[0,1].axhline(HORIZON*1e6,color='k',ls='--',label='Horizon 150.5 µm');axs[0,1].set(ylabel='Equivalent diameter (µm)',title='Equal reference-volume sphere');axs[0,1].legend()
 axs[1,0].bar(ids,df.reference_volume_m3*1e9,color=colors);axs[1,0].set(ylabel='Reference volume (mm³)',title='Nodal quadrature: N × dx³')
 for i,axis in enumerate('xyz'):axs[1,1].bar(ids+(i-1)*.24,df[f'{axis}_span_m']*1e6,width=.24,label=axis.upper())
 axs[1,1].set(ylabel='Node-center bounding span (µm)',title='Actual 3D domain extents');axs[1,1].legend(ncol=3)
 for ax in axs.flat:ax.set_xticks(ids);ax.set_xlabel('Strong domain ID');ax.grid(axis='y',alpha=.25);ax.set_axisbelow(True)
 fig.suptitle('Exploratory material scales · these regions are not predefined fragments\n6-face adjacency defines region count; the 150.5 µm material horizon can bridge weak separators',fontsize=15);save(fig,'FIG_03_DOMAIN_SIZE_DISTRIBUTION')
 dump(O/'P0_FIGURE_AUDIT.json',{'source_field_sha256':sha(O/'BIPHASIC_FIELD.npz'),'figures':{p.name:sha(p) for p in sorted((F/'figures').glob('*'))},'geometry':'Physical node-center coordinates, common equal 3D scale, no displacement or phase-boundary separation.'})
if __name__=='__main__':main()
