from common import *
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import MaxNLocator
plt.rcParams.update({'font.size':12,'axes.titlesize':13,'axes.labelsize':11,'savefig.dpi':220,'svg.fonttype':'none'})
z=np.load(O/'HISTORICAL_SAME_COORDINATE_DECOMPOSITION.npz');x=z['x']*1e6;X=z['X']*1e6;b=z['bonds'];br=z['broken_bond_id']-1
fig=plt.figure(figsize=(16,10),layout='constrained');axes=[fig.add_subplot(2,3,i+1,projection='3d') for i in range(6)]
def base(ax,coords,alive):
 seg=coords[b[alive]];ax.add_collection3d(Line3DCollection(seg,colors='#808890',linewidths=.23,alpha=.025,rasterized=True))
 ax.set(xlim=(-175,185),ylim=(-175,175),zlim=(-175,175),xlabel='X (um)',ylabel='Y (um)',zlabel='Z (um)');ax.view_init(18,-55);ax.set_box_aspect((1,1,1))
 for a in [ax.xaxis,ax.yaxis,ax.zaxis]:a.set_major_locator(MaxNLocator(3));a._axinfo['grid'].update(color='#a8adb3',linewidth=.6)
 return ax
for i,coord in enumerate([X,x]):
 ax=base(axes[i],coord,z['alive_pre'] if i==0 else z['alive_post']);ax.scatter(*coord.T,s=11,color='#889ca9',alpha=.65,depthshade=False)
 ax.add_collection3d(Line3DCollection(coord[b[br]],colors='#cb3036',linewidths=2.5))
 ids=np.unique(b[br]);ax.scatter(*coord[ids].T,s=35,color='#cb3036',depthshade=False)
 for j in ids:ax.text(*coord[j],str(j+1),fontsize=10,color='#75131a')
 ax.set_title(['(a) Reference locations; 5 true broken bonds','(b) Same actual event coordinates; post topology'][i])
axes[2].set_axis_off()
axes[2].text2D(.0,.9,'Historical static event\n341 material nodes; 12,047 initial bonds\nAtom IDs: 1-based original identities\nBond IDs: original input rows\n\nF_pre: B- theta, B- force\nF_mixed: B- theta, B+ force\nF_post: B+ theta, B+ force\n\nNo physical event time or integrated velocity\nNo displacement magnification\nForce arrows below: one common scale',transform=axes[2].transAxes,fontsize=12,va='top',linespacing=1.6)
vmin=min(z['theta_pre'].min(),z['theta_post'].min());vmax=max(z['theta_pre'].max(),z['theta_post'].max());norm=Normalize(vmin,vmax)
for ax,key,title in zip(axes[3:5],['theta_pre','theta_post'],['(c) Pre-break dilatation theta','(d) Post-break dilatation theta']):
 base(ax,x,z['alive_post']);ax.scatter(*x.T,c=z[key],cmap='viridis',norm=norm,s=16,depthshade=False);ax.set_title(title)
fig.colorbar(ScalarMappable(norm=norm,cmap='viridis'),ax=axes[3:5],orientation='horizontal',shrink=.72,pad=.03,label='Dilatation (dimensionless)')
ax=base(axes[5],x,z['alive_post']);ax.scatter(*x.T,s=6,color='#a4a9b2',alpha=.35,depthshade=False)
for key,col,label in [('DeltaF_topology_N','#257ba4','Topology force change'),('DeltaF_timelayer_N','#cf4e38','Time-layer force difference')]:
 f=z[key];mask=np.linalg.norm(f,axis=1)>1e-9;ax.quiver(*x[mask].T,*(f[mask]*8e6).T,color=col,arrow_length_ratio=.2,linewidth=1.2,label=label)
ax.set_title('(e) Two distinct force differences (N)');ax.legend(loc='upper left',fontsize=9);ax.text2D(.01,.02,'Arrow scale: 1 uN = 8 um on axes',transform=ax.transAxes,fontsize=10)
fig.suptitle('LPS fracture time-layer decomposition | saved native static event, same coordinates',fontsize=17)
fig.savefig(F/'Fig1_static_force_decomposition.png');fig.savefig(F/'Fig1_static_force_decomposition.svg');plt.close(fig)
print('Saved Fig1_static_force_decomposition')
