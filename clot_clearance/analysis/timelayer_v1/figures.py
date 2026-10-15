from common import *
import numpy as np,time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d.art3d import Line3DCollection
beg=time.perf_counter();plt.rcParams.update({'font.size':12,'axes.titlesize':13,'axes.labelsize':12,'legend.fontsize':11,'savefig.dpi':230,'svg.fonttype':'none','axes.grid':True,'grid.alpha':.45,'grid.color':'#9ca5ad','figure.facecolor':'white'})
names=['dt_0p5us','dt_0p25us'];colors=['#136a98','#ce5337'];labels=['dt = 0.50 us, every step','dt = 0.25 us, every step'];D={n:np.genfromtxt(O/n/'analysis_scalars.csv',delimiter=',',names=True) for n in names};S={n:np.load(O/n/'native_states.npz') for n in names};g=np.load(O/'historical_static/reference.npz');b=g['bonds'];X=g['X']*1e6
def save(fig,name):
 fig.savefig(F/(name+'.png'));fig.savefig(F/(name+'.svg'));plt.close(fig)
def line(ax,key,ylabel,scale=1.,post=False):
 for n,col,label in zip(names,colors,labels):
  d=D[n];fn=ax.step if post else ax.plot
  fn(d['time_s']*1e6,d[key]*scale,label=label,color=col,linewidth=1.8,**({'where':'post'} if post else {}))
 ax.set(xlabel='Physical time (us)',ylabel=ylabel)
def base(ax,x,alive):
 ax.add_collection3d(Line3DCollection(x[b[alive]],colors='#738693',linewidths=.3,alpha=.03,rasterized=True));ax.set(xlim=(-175,185),ylim=(-175,175),zlim=(-175,175),xlabel='X (um)',ylabel='Y (um)',zlabel='Z (um)');ax.set_box_aspect((1,1,1));ax.view_init(18,-55)
 for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.set_major_locator(MaxNLocator(3));axis._axinfo['grid'].update(color='#a1a9b1',linewidth=.6)
# Fig1: actual first dynamic event, same coordinates in all force panels.
z=np.load(O/names[0]/'decomposition_0010.npz');x=z['x']*1e6;changed=z['alive_pre']&~z['alive_post'];ends=np.unique(b[changed]);fig=plt.figure(figsize=(17,10),layout='constrained');axes=[fig.add_subplot(2,3,k+1,projection='3d') for k in range(6)]
for ax,coords,title,live in [(axes[0],X,'(a) First broken bond: ID 202; atoms 4-90',z['alive_pre']),(axes[1],x,'(b) Actual x*: pre/post topology',z['alive_post'])]:
 base(ax,coords,live);ax.scatter(*coords.T,s=12,c='#8da1ad',alpha=.65,depthshade=False);ax.add_collection3d(Line3DCollection(coords[b[changed]],colors='#c7333d',linewidths=3));ax.scatter(*coords[ends].T,s=50,color='#c7333d',depthshade=False)
 for j in ends:ax.text(*coords[j],str(j+1),fontsize=12,color='#7c1c22')
 ax.set_title(title)
norm=Normalize(min(z['theta_pre'].min(),z['theta_post'].min()),max(z['theta_pre'].max(),z['theta_post'].max()))
for ax,key,title in [(axes[2],'theta_pre','(c) theta from pre-break topology B-'),(axes[3],'theta_post','(d) theta from post-break topology B+')]:
 base(ax,x,z['alive_post']);ax.scatter(*x.T,c=z[key],cmap='viridis',norm=norm,s=17,depthshade=False);ax.set_title(title)
fig.colorbar(ScalarMappable(norm=norm,cmap='viridis'),ax=[axes[2],axes[3]],orientation='horizontal',shrink=.55,pad=.02,label='Dilatation (dimensionless), shared scale')
for ax,key,title,col in [(axes[4],'DeltaF_topology_N','(e) F_post - F_pre: topology change','#166b99'),(axes[5],'DeltaF_timelayer_N','(f) F_mixed - F_post: time-layer difference','#cd5138')]:
 base(ax,x,z['alive_post']);ax.scatter(*x.T,s=7,color='#a1aab1',alpha=.4,depthshade=False);f=z[key];sel=np.linalg.norm(f,axis=1)>1e-9
 ax.quiver(*x[sel].T,*(f[sel]*40e6).T,color=col,arrow_length_ratio=.22,linewidth=1.4);ax.set_title(title);ax.text2D(.03,.90,'Common arrow scale: 1 uN = 40 um',transform=ax.transAxes,fontsize=10,bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
fig.suptitle('Actual native first fracture: t = 5.00 us, dt = 0.50 us | same-coordinate force decomposition\n341 original material IDs; original bond IDs; no displacement magnification',fontsize=16)
save(fig,'Fig1_native_timelayer_force')
# Fig2: time alignment and actual event observations.
fig,ax=plt.subplots(2,3,figsize=(16,9),layout='constrained');aa=ax.ravel()
bond_rows={bid:k for k,bid in enumerate([202,5702,5703,5704,5715])}
for n,col,label in zip(names,colors,labels):
 e=np.genfromtxt(O/n/'native_break_events.csv',delimiter=',',names=True);aa[0].scatter(e['time_s']*1e6,[bond_rows[int(v)] for v in e['bond_id']],label=label,color=col,s=65,marker='o' if n==names[0] else 'x')
aa[0].set_yticks(list(bond_rows.values()),list(bond_rows.keys()))
aa[0].set(xlabel='Actual event time (us)',ylabel='Original bond ID',title='(a) Discrete native events; no interpolation',xlim=(4.3,9))
line(aa[1],'max_prebreak_stretch_utilization','max stretch / s00');aa[1].axhline(1,color='#444',ls=':');aa[1].set(xlim=(0,12),title='(b) Evaluated on topology before this step breaks')
line(aa[2],'cumulative_breaks','Broken bonds',post=True);aa[2].set(xlim=(0,12),title='(c) Same final 5 bonds; different grouping')
line(aa[3],'neighbor4_vx_m_s','Atom 4 vx (m/s)');aa[3].set_title('(d) Direct native full-step velocity: first endpoint')
line(aa[4],'vmax_free_m_s','Max free-node speed (m/s)');aa[4].set_title('(e) Free nodes; prescribed atom 90 excluded')
line(aa[5],'max_step_displacement_m','Max one-step displacement (um)',1e6);aa[5].set_title('(f) Actual accepted position increments')
aa[0].legend();fig.suptitle('Timestep comparison | same approved 15 um / 10 us loading, 50 us observation\nNative accepted states: 0.50 us and 0.25 us sampling; first breaks at 5.00 us and 4.75 us',fontsize=16)
save(fig,'Fig2_timestep_comparison')
# Fig3: impulses remain labelled as proxies, actual velocity not a proxy.
fig,ax=plt.subplots(2,3,figsize=(16,9),layout='constrained');aa=ax.ravel()
line(aa[0],'lag_force_max_free_N','Max free-node force difference (uN)',1e6);aa[0].set(xlim=(0,12),title='(a) F_mixed - F_post at actual states')
line(aa[1],'J_lag_proxy_max_free_Ns','Max nodal lag impulse proxy (pN s)',1e12);aa[1].set_title('(b) Actual VV endpoint weights; free DOFs')
line(aa[2],'neighbor4_vx_m_s','Atom 4 actual vx (m/s)');aa[2].set(xlim=(3,11),title='(c) Direct accepted velocity near events')
window_v=np.concatenate([d['neighbor4_vx_m_s'][(d['time_s']>=3e-6)&(d['time_s']<=11e-6)] for d in D.values()]);pad=.12*np.ptp(window_v);aa[2].set_ylim(window_v.min()-pad,window_v.max()+pad)
for n,col,label in zip(names,colors,labels):
 d=D[n];aa[3].plot(d['time_s']*1e6,d['DeltaP_node4_x_Ns']*1e12,c=col,label=label+' actual');aa[3].plot(d['time_s']*1e6,d['J_lag_node4_x_Ns']*1e12,c=col,ls='--',label=label+' proxy')
aa[3].set(xlabel='Physical time (us)',ylabel='Atom 4 x-momentum / proxy (pN s)',title='(d) Actual momentum includes all internal response');aa[3].legend(fontsize=8)
line(aa[4],'K_free_J','Free-node kinetic energy (pJ)',1e12);aa[4].set_title('(e) Native full-step velocity, physical mass')
line(aa[5],'damping_work_J','Numerical damping work (pJ)',1e12);aa[5].set_title('(f) Exactly zero: no damping in this specimen');aa[5].text(.06,.76,'Pair zero: no contact work\nNo background fluid force\nOriginal tiny-specimen settings retained',transform=aa[5].transAxes,fontsize=11)
aa[0].legend(fontsize=10);fig.suptitle('Force-difference impulse proxy ≠ measured velocity jump\nNative velocities after final_integrate; force differences evaluated on the actual trajectory',fontsize=16)
save(fig,'Fig3_impulse_and_dynamics')
# Fig4: signed energy ledger and actual unchanged final connectivity.
fig=plt.figure(figsize=(17,14),layout='constrained');grid=fig.add_gridspec(3,3);aa=[fig.add_subplot(grid[i//3,i%3]) for i in range(6)]
for ax,key,ylabel,title in [(aa[0],'U_LPS_J','Independent recoverable U (pJ)','(a) Full gradient-consistent LPS potential'),(aa[1],'W_boundary_J','Boundary work (pJ)','(b) Work into deformable system'),(aa[2],'K_free_J','Free-node kinetic energy (pJ)','(c) Driven-node inertia excluded'),(aa[3],'cumulative_DeltaU_topology_J','Cumulative signed DeltaU (pJ)','(d) Same coordinates: U_post - U_pre'),(aa[4],'energy_residual_J','Energy residual R (pJ)','(e) R includes initial prestrain energy')]:
 line(ax,key,ylabel,1e12,post=(key=='cumulative_DeltaU_topology_J'));ax.set_title(title)
for n,col,label in zip(names,colors,labels):
 d=D[n];aa[5].plot(d['time_s']*1e6,d['W_lag_on_actual_path_J']*1e12,c=col,label=label+' lag work');aa[5].plot(d['time_s']*1e6,d['W_topology_switch_J']*1e12,c=col,ls='--',label=label+' topology time placement')
aa[5].set(xlabel='Physical time (us)',ylabel='On-path residual terms (pJ)',title='(f) Distinct signed contributions');aa[5].legend(fontsize=8)
for k,n in enumerate(names):
 ax=fig.add_subplot(grid[2,k],projection='3d');s=S[n];xx=s['x'][-1]*1e6;alive=np.unpackbits(s['alive_packed'][-1])[:len(b)].astype(bool);base(ax,xx,alive);ax.scatter(*xx.T,s=12,color=colors[k],depthshade=False);ax.add_collection3d(Line3DCollection(xx[b[~alive]],colors='#c42d37',linewidths=2));ax.scatter(*xx[89],s=60,c='#222',marker='s');ax.set_title(f'({chr(103+k)}) t=50 us, dt={0.5 if k==0 else 0.25} us\nOne connected 341-node component')
ax=fig.add_subplot(grid[2,2]);ax.axis('off');ax.text(.02,.93,'Energy ledger (SI joules)\nR = K_free + U - (K0 + U0)\n      - W_boundary - sum(DeltaU_topology)\n\nDamping = 0; contact = 0\nIndependent U, not native PE x volume\n\nBoth runs: same five broken bond IDs\n202, 5702, 5703, 5704, 5715\nRed segments mark deleted bonds\nGray segments are true surviving bonds\n\nTopology release is not material Gc',va='top',fontsize=12,linespacing=1.5)
aa[0].legend(fontsize=10);fig.suptitle('Energy and actual fracture topology | approved small-specimen dynamics, 0-50 us\nSigned topology release, lag work and integration error are separate; no material Gc inferred',fontsize=16)
save(fig,'Fig4_energy_and_topology')
dump(O/'FIGURE_PERFORMANCE.json',{'wall_s':time.perf_counter()-beg,'groups':4,'formats':['PNG','SVG'],'source':'actual native states and same-coordinate independent diagnostics','no_interpolated_fracture_events':True});print('Four figure groups saved')
