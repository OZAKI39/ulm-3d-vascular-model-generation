"""Actual accepted coordinates x1; solved template fields; explicit absent events."""
from analyze import *
from acoustic_streaming_shear import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.ticker import MaxNLocator
import textwrap
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d.art3d import Poly3DCollection,Line3DCollection
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':14,'axes.titlesize':16,'axes.labelsize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':12,'svg.fonttype':'none','savefig.facecolor':'white'})
MANIFEST=[]
def export(fig,name,caption):
 F.mkdir(parents=True,exist_ok=True);fig.savefig(F/(name+'.png'),dpi=190,bbox_inches='tight');fig.savefig(F/(name+'.svg'),bbox_inches='tight');plt.close(fig);MANIFEST.append(dict(name=name,caption=caption))
def axes3(ax,title):
 ax.set(xlabel='X (µm)',ylabel='Y (µm)',zlabel='Z (µm)',title=title);ax.view_init(16,-63);ax.set_box_aspect([1,.65,1.2]);ax.grid(True)
 ax.xaxis.set_major_locator(MaxNLocator(3));ax.yaxis.set_major_locator(MaxNLocator(3));ax.zaxis.set_major_locator(MaxNLocator(3));ax.set_xlabel('X (µm)',labelpad=10);ax.set_ylabel('Y (µm)',labelpad=10);ax.set_zlabel('Z (µm)',labelpad=12)
 for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis._axinfo['grid'].update(color='#a3a9b1',linewidth=.65)
def cloud(ax,ref,s,title,broken=False):
 n=int(ref['clot_count']);x=s['x']*1e6;lab,sizes,att=graph.topology(ref,s['d']);attached=att[lab]
 ax.scatter(*x[:n][attached].T,c='#a7acb4',s=7,alpha=.5,depthshade=False,axlim_clip=True);ax.scatter(*x[n:].T,c='#3a4655',s=7,alpha=.6,depthshade=False,axlim_clip=True)
 free=np.flatnonzero(~att)
 for i,label in enumerate(free):
  ids=np.flatnonzero(lab==label);ax.scatter(*x[ids].T,c=[plt.cm.tab20(i%20)],s=20 if len(ids)>2 else 12,depthshade=False,axlim_clip=True)
 if broken:
  ids=np.flatnonzero((s['d']==1)&~ref['initial_cut']&np.isin(ref['bond_types'],[1,2]));ends=np.unique(ref['bonds'][ids]);ends=ends[ends<n];ends=ends[attached[ends]]
  if len(ends):ax.scatter(*x[ends].T,c='#be4b39',s=8,alpha=.75,depthshade=False,axlim_clip=True)
 axes3(ax,title);ax.set_xlim(-600,600);ax.set_ylim(-500,500);ax.set_zlim(-650,650)
 outside=np.any((x[:n]<[-600,-500,-650])|(x[:n]>[600,500,650]),axis=1).sum()
 ax.text2D(.01,.92,f't={float(s["time"])*1e3:.3f} ms | free nodes={int(sizes[~att].sum())}\n20–200: {int(np.sum((sizes[~att]>=20)&(sizes[~att]<=200)))} | outside common ROI: {outside}',transform=ax.transAxes,fontsize=11)
 return free
def main(inputs_only=False):
 ref=load_npz(O/'provenance/reference.npz');analysis=read(O/'ANALYSIS.json') if not inputs_only else {};cases={k:O/'cases'/k for k in ['stokes','inverse']};rows={k:histories(p) for k,p in cases.items()} if not inputs_only else {};load=ShearLoad(ref,'stokes');audit=read(O/'ACOUSTIC_AND_BUBBLE_INPUT_AUDIT.json')
 # 1. Data status is part of the diagram, not an implied completed coupling.
 fig,ax=plt.subplots(figsize=(20,5.5));ax.axis('off');boxes=[('Acoustics','ACTUAL_DATA','F1 wave solve / F10S proxy\n2.5 MHz; distinct geometries'),('Bubble response','ACTUAL_DATA','APECSS volume harmonics\nF10S center trajectories'),('Bubble streaming','NOT_AVAILABLE','Clot registration + near-wall\nfirst/second-order flow missing'),('Viscous shear','PHYSICS_APPROXIMATION','Prescribed Stokes counter-shear\n25 µm gap; U=4 mm/s'),('Surface → native LPS','COMPUTED','Original-node triangles → N\nIrreversible damage + live bonds')]
 for i,(title,status,body) in enumerate(boxes):
  x=.008+i*.2;color='#fbe9e5' if status=='NOT_AVAILABLE' else '#eaf0f5'
  ax.add_patch(FancyBboxPatch((x,.30),.177,.61,boxstyle='round,pad=.008',transform=ax.transAxes,facecolor=color,edgecolor='#546575'))
  ax.text(x+.007,.84,title,transform=ax.transAxes,va='top',fontsize=12)
  ax.text(x+.007,.68,status,transform=ax.transAxes,va='top',fontsize=10.5,weight='bold')
  body='\n'.join(textwrap.fill(line,26) for line in body.split('\n'));ax.text(x+.007,.54,body,transform=ax.transAxes,va='top',fontsize=11.5)
  if i<4:ax.annotate('',xy=(x+.196,.62),xytext=(x+.179,.62),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',color='#536878',lw=1.8,linestyle='--' if i in [1,2] else '-'))
 ax.text(.02,.08,'Missing streaming link: no acoustic pressure → clot shear conversion has been established.\nThe high-amplitude case is INVERSE_SHEAR_REQUIREMENT / MECHANICAL_REQUIREMENT_ONLY.',transform=ax.transAxes,fontsize=14)
 export(fig,'Fig1_pipeline_status','Actual acoustic/bubble inputs and missing streaming link are explicitly separated from the computed Stokes approximation.')
 # 2. Actual analytic-mode fluid solution, not reconstructed bubbles.
 fig,axs=plt.subplots(1,3,figsize=(18,5.2));y=np.linspace(-350e-6,350e-6,101);z=np.linspace(0,GAP,51);Y,Z=np.meshgrid(y,z,indexing='ij');field=load.flow.field(Y,Z,HOLD)
 im=axs[0].pcolormesh(y*1e6,z*1e6,field['velocity_m_s'][...,0].T*1e3,cmap='coolwarm',vmin=-4,vmax=4,shading='auto');fig.colorbar(im,ax=axs[0],label='uₓ (mm/s)');axs[0].set(xlabel='Transverse y (µm)',ylabel='Gap z (µm)',title='Prescribed Stokes counter-flow')
 for tt in [.0005,.001,.0024,.0029,.0034]:axs[1].plot(load.flow.modal_velocity(z,tt)*1e3,z*1e6,label=f'{tt*1e3:g} ms')
 axs[1].set(xlabel='uₓ at y=175 µm (mm/s)',ylabel='Gap z (µm)',title='No-slip startup / decay');axs[1].legend()
 evidence=load_npz(O/'MEVO_ACTUAL_EVIDENCE.npz');pos=evidence['positions_m'][0]*1e6;axs[2].scatter(pos[:,0],pos[:,1],s=12,c='#56819d');axs[2].set(xlabel='F10S x (µm)',ylabel='F10S y (µm)',title='Actual initial bubble centers');axs[2].text(.02,.02,'Own coordinate system\nClot gap/registration: UNAVAILABLE',transform=axs[2].transAxes,fontsize=12)
 fig.tight_layout();export(fig,'Fig2_stokes_flow_and_separate_bubbles','Solved prescribed flow and actual unregistered bubble positions. Bubble center velocity is not fluid velocity.')
 # 3. Actual CFD triangles and both mapping choices.
 z=load_npz(O/'V3_SURFACE_MAPPING_DATA.npz');origin=z['X'].mean(0);fig=plt.figure(figsize=(19,6));ax=fig.add_subplot(131,projection='3d');mag=np.linalg.norm(z['shear_traction_Pa'],axis=1);norm=Normalize(0,max(mag));polys=Poly3DCollection((z['triangles_xyz_m']-origin)*1e6,facecolors=plt.cm.viridis(norm(mag)),linewidth=0,alpha=.9);ax.add_collection3d(polys);x=(z['X']-origin)*1e6;ax.auto_scale_xyz(x[:,0],x[:,1],x[:,2]);axes3(ax,'Actual CFD exposed triangles\nlocal-patch Stokes approximation');fig.colorbar(ScalarMappable(norm=norm,cmap='viridis'),ax=ax,shrink=.6,label='Shear (Pa)')
 for j,(key,title) in enumerate([('signed_node_force_N','Signed affine map'),('nonnegative_node_force_N','Nonnegative local comparison')],start=2):
  ax=fig.add_subplot(1,3,j,projection='3d');f=z[key];mag=np.linalg.norm(f,axis=1);active=mag>1e-15;sc=ax.scatter(*x.T,c=mag*1e9,s=7,cmap='viridis');stride=np.flatnonzero(active)[::8];ax.quiver(*x[stride].T,*(f[stride]*2e10).T,color='#bf472f',arrow_length_ratio=.3,linewidth=.7);axes3(ax,title);fig.colorbar(sc,ax=ax,shrink=.6,label='Node force (nN)')
 allx=np.vstack([z['X'],z['triangles_xyz_m'].reshape(-1,3)])-origin;lo=allx.min(0)*1e6;hi=allx.max(0)*1e6
 for panel in fig.axes:
  if hasattr(panel,'view_init'):panel.set_xlim(lo[0]-30,hi[0]+30);panel.set_ylim(lo[1]-30,hi[1]+30);panel.set_zlim(lo[2]-30,hi[2]+30)
 fig.text(.5,.015,'Arrows: 1 nN → 20 µm | pressure increment = 0 | all coordinates relative to material centroid',ha='center',fontsize=12)
 fig.subplots_adjust(bottom=.23,left=.02,right=.96,wspace=.22);export(fig,'Fig3_actual_V3_surface_mapping','3611 actual CFD faces; no artificial outer shell. Existing signed and new positive maps are compared, neither is admitted as validated V3 streaming.')
 # 4. Common mode with actual triangle force budget.
 fig,axs=plt.subplots(1,4,figsize=(21,5.1),layout='constrained');cent=load.cent*1e6
 for ax,tt in zip(axs[:3],[.0005,.0015,.0024]):
  shear=load.faces(tt)[:,0]/load.area;im=ax.scatter(cent[:,0],cent[:,1],c=shear,cmap='coolwarm',vmin=-.16,vmax=.16,s=27);skip=np.arange(0,len(cent),8);ax.quiver(cent[skip,0],cent[skip,1],shear[skip],np.zeros(len(skip)),angles='xy',scale_units='xy',scale=.003,color='#323c48');ax.set(xlabel='X (µm)',ylabel='Y (µm)',title=f'{tt*1e3:g} ms, reference top');ax.set_aspect('equal')
 fig.colorbar(im,ax=list(axs[:3]),shrink=.7,label='Signed shear τₓ (Pa)');fig.suptitle('Reference shear map | arrows: 0.1 Pa → 33.3 µm',fontsize=15)
 for k in cases:
  tt=np.linspace(0,END,721);axs[3].plot(tt*1e3,load.amplitude(tt)*np.max(abs(load.mode))/load.physical_peak,label=k,ls='-' if k=='stokes' else '--')
 axs[3].set(xlabel='Time (ms)',ylabel='Peak / own maximum',title='Same frozen time protocol');axs[3].axvline(OFF*1e3,ls='--',c='.4');axs[3].legend()
 export(fig,'Fig4_spatiotemporal_shear','Continuous transverse counter-shear on the original top triangles. B has the same shape/time and a mechanical-only 1000 Pa peak.')
 if inputs_only:
  dump(O/'INPUT_FIGURE_MANIFEST.json',dict(figures=MANIFEST));print('INPUT_FIGURES_DONE',len(MANIFEST));return
 # 5. Only events that actually exist.
 for name,case in cases.items():
  keys=['initial','first_softening','first_failure','first_multinode','max_multifragment','final'];fig=plt.figure(figsize=(18,11))
  for j,key in enumerate(keys,1):
   ax=fig.add_subplot(2,3,j,projection='3d');p=case/'keys'/(key+'.npz')
   title='best objective snapshot' if key=='max_multifragment' else key.replace('_',' ')
   if p.exists():cloud(ax,ref,load_npz(p),title,broken=True)
   else:
    axes3(ax,title);ax.set_xlim(-600,600);ax.set_ylim(-500,500);ax.set_zlim(-650,650);ax.text2D(.15,.45,'NOT_OBSERVED',transform=ax.transAxes,fontsize=20)
  fig.suptitle(f'{name.upper()} | accepted native states, displacement ×1\nRed: attached nodes incident to natural bond failures; colors: free components',fontsize=17);fig.tight_layout();export(fig,'Fig5_'+name+'_fracture','Original live material topology and actual positions. Red attached nodes are endpoints of naturally failed original bonds; no nonexistent broken connections are drawn across space. Common-ROI clipping is counted; all final nodes appear in Fig7. Missing stages remain NOT_OBSERVED.')
 # 6. Full requested bins including zero counts.
 fig,axs=plt.subplots(1,2,figsize=(15,5));bins=['1–2','3–13','14–33','34–113','114–200','>200'];table=list(csv.DictReader((O/'FRAGMENT_SIZE_DISTRIBUTION.csv').open()))
 for j,name in enumerate(cases):
  q=[r for r in table if r['case']==name and r['state']=='final'];xx=np.arange(6)+(j-.5)*.35;axs[0].bar(xx,[float(r['component_count']) for r in q],.35,label=name);axs[1].bar(xx,[float(r['mass_fraction_initial'])*100 for r in q],.35,label=name)
 for ax in axs:ax.set_xticks(np.arange(6),bins);ax.set_xlabel('Original nodes per free component');ax.legend()
 axs[0].set(ylabel='Component count',title='Final actual free components');axs[1].set(ylabel='Initial material mass (%)',title='Mass in each size bin');fig.tight_layout();export(fig,'Fig6_fragment_size_distribution','Empty bins are measured zeros. Mass is original material-node volume, anchors excluded.')
 # 7. Common camera and full extent, true component IDs only.
 fig=plt.figure(figsize=(16,7))
 for j,(name,case) in enumerate(cases.items(),1):
  ax=fig.add_subplot(1,2,j,projection='3d');s=load_npz(case/'keys/final.npz');cloud(ax,ref,s,name+' final');x=s['x']*1e6;ax.set_xlim(min(-600,x[:,0].min()-30),max(600,x[:,0].max()+30));ax.set_ylim(min(-500,x[:,1].min()-30),max(500,x[:,1].max()+30));ax.set_zlim(min(-650,x[:,2].min()-30),max(650,x[:,2].max()+30))
 fig.suptitle('True material connected components | gray: attached | dark: anchors',fontsize=17);fig.tight_layout();export(fig,'Fig7_true_3D_components','Full final extents, actual coordinates ×1 with independently fitted axis ranges. Inverse-case free nodes travel tens of millimeters under dead reference traction without solved fluid drag; this is not vessel transport. Spatial neighbors/contact pairs never define fragments. Colors may cycle for numerous small components.')
 # 8. Capability not confused with inverse demand.
 fig,axs=plt.subplots(1,2,figsize=(16,5.5));axs[0].bar(['Stokes template','Mechanical requirement'],[load.physical_peak,1000],color=['#4b8197','#be6a40']);axs[0].set_yscale('log');axs[0].set(ylabel='Peak actual-face shear (Pa)',title='Same specimen / mode / duration');axs[0].text(.03,.88,'Bubble streaming capability: UNAVAILABLE',transform=axs[0].transAxes,fontsize=13)
 axs[1].axis('off');text='Observed native response\n\n'
 for name,q in analysis.items():text+=f'{name}: target 20–200 count max {q["history_max_N20_200"]}\nfinal free nodes {int(q["final"]["detached_nodes"])}, largest history {q["history_max_free_nodes"]}\n\n'
 text+=f'{1000/load.physical_peak:.2f}× compares imposed templates only.\nNo target-success threshold or acoustic shortfall is inferred.\nAPECSS cache μ=2 mPa·s; this baseline μ=1 mPa·s.';axs[1].text(0,.97,text,va='top',fontsize=14)
 fig.tight_layout();export(fig,'Fig8_capability_and_requirement','Physical acoustic capability remains undetermined; the large ratio is a controlled mechanical template scaling, not a validated bubble deficit.')
 # 9. Work and shear torque without net-force cancellation hiding loading.
 fig,axs=plt.subplots(2,3,figsize=(18,10))
 for row,(name,rr) in enumerate(rows.items()):
  t=np.array([q['time_s'] for q in rr])*1e3
  for key,label in [('W_J','External work'),('U_J','Stored LPS'),('K_J','Kinetic'),('D_J','Damage dissipation'),('Uc_J','Pair contact')]:axs[row,0].plot(t,np.array([q[key] for q in rr])*1e9,label=label)
  axs[row,0].set(title=name+' energy',xlabel='Time (ms)',ylabel='Energy (nJ)');axs[row,0].legend(fontsize=10)
  if name=='inverse':
   inset=axs[row,0].inset_axes([.55,.12,.42,.34])
   for key,color in [('U_J','tab:orange'),('D_J','tab:red'),('Uc_J','tab:purple')]:inset.plot(t,np.array([q[key] for q in rr])*1e9,label=key[:2],color=color)
   inset.set_title('Stored / damage / contact zoom',fontsize=8);inset.tick_params(labelsize=7);inset.set_ylabel('nJ',fontsize=8)
  axs[row,1].plot(t,[q['sum_force_norm_N']*1e6 for q in rr],label='Sum |node force|');axs[row,1].plot(t,[abs(q['net_force_x_N'])*1e6 for q in rr],label='|Net force|');axs[row,1].set(title='Opposing forces and cancellation',xlabel='Time (ms)',ylabel='Force (µN)');axs[row,1].legend()
  axs[row,2].plot(t,[q['relative_R']*100 for q in rr],label='Normalized numerical residual');axs[row,2].axhline(1,ls='--',c='#a6332a',label='1% gate');axs[row,2].set(title='Residual stays a residual',xlabel='Time (ms)',ylabel='Energy residual (%)');axs[row,2].legend(fontsize=10)
 fig.tight_layout();export(fig,'Fig9_energy_and_force','Actual material external work, storage, kinetic/damage/contact energy and independent numerical balance; zero net force does not mean no shear work.')
 fig=plt.figure(figsize=(18,6));bridge=read(O/'CRACK_BRIDGE_AUDIT.json');networks=read(O/'DAMAGE_NETWORK_HISTORY.json')
 for j,(name,case) in enumerate(cases.items(),1):
  ax=fig.add_subplot(1,3,j,projection='3d');s=load_npz(case/'keys/final.npz');region=bridge['cases'][name]['final']['regions'][0];ids=np.array(region['original_surviving_bond_ids'],int)-1;cloud(ax,ref,s,name+' surviving bridge bonds');ax.add_collection3d(Line3DCollection(s['x'][ref['bonds'][ids]]*1e6,colors='#2d858d',linewidths=.5,alpha=.23));ax.set_title(f'{name}: {len(ids)} live bridge bonds\nUpper diagnostic region: {region["candidate_nodes"]} nodes',fontsize=14)
 ax=fig.add_subplot(133)
 for name,q in networks.items():
  h=q['all_material_including_interface']['history'];ax.plot([r['time_s']*1e3 for r in h],[r['separate_failed_material_networks'] for r in h],label=name)
 ax.set(xlabel='Time (ms)',ylabel='Separate failed-edge networks',title='Natural damage-network connectivity');ax.legend();ax.text(.02,.03,'Failed-edge merging is not proof\nof a closed fracture surface.',transform=ax.transAxes,fontsize=12)
 fig.tight_layout();export(fig,'Fig10_3D_bridges_and_damage_networks','Live original bridging bonds and failed-edge networks; diagnostic volume is not a prescribed fragment or loading mask.')
 dump(O/'FIGURE_MANIFEST.json',dict(figures=MANIFEST,actual_positions_scale=1,implied_missing_events=False));print('FIGURES_DONE',len(MANIFEST))
if __name__=='__main__':
 import argparse
 parser=argparse.ArgumentParser();parser.add_argument('--inputs-only',action='store_true');args=parser.parse_args();main(args.inputs_only)
