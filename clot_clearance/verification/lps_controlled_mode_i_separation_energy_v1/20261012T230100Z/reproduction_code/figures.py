"""Static evidence only. No synthetic fracture frames or physical-time curves."""
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import time,resource,html
start=time.perf_counter();plt.rcParams.update({'font.family':'DejaVu Sans','font.size':14,'axes.titlesize':17,'axes.labelsize':15,'xtick.labelsize':13,'ytick.labelsize':13,'legend.fontsize':12,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.color':'#9eabb5','grid.alpha':.65,'grid.linewidth':.65,'savefig.facecolor':'white','svg.fonttype':'none'})
blue='#2b7099';red='#b84c50';orange='#d88927';gray='#9aa7b1';purple='#9067a8';dark='#263746'
rt=np.load(O/'reference_topology.npz');s=np.load(O/'elastic_precheck.npz');t=np.load(O/'topology_audit.npz');X=rt['X'];b=rt['bonds'];alive=t['alive'];cross=t['cross_z0'];N=len(X);x=s['x'];p=json.loads((O/'ELASTIC_PRECHECK.json').read_text());load=json.loads((O/'LOADING_AND_REACTION_AUDIT.json').read_text());regions=json.loads((O/'SPATIAL_REPORTING_AUDIT.json').read_text())['regions'];files=[];samples={}
def save(fig,name):
 for ext in ['png','svg']:
  fig.savefig(F/(name+'.'+ext),dpi=240,bbox_inches='tight');files.append(name+'.'+ext)
 plt.close(fig)
def spatial(ax):
 ax.set(xlabel='X (mm)',ylabel='Z (mm)',xlim=(-2.55,2.55),ylim=(-1.1,1.1));ax.set_aspect('equal');ax.axvline(-1,color=orange,lw=1.4,ls='--');ax.grid(alpha=.6)
sheet=np.isclose(X[:,1],-25e-6,atol=1e-12);sheetb=alive&sheet[b[:,0]]&sheet[b[:,1]];bid=np.flatnonzero(sheetb);samples['XZ_sheet_atom_ids']=np.flatnonzero(sheet)+1;samples['XZ_sheet_alive_bond_ids']=bid+1
fig=plt.figure(figsize=(17,8));fig.subplots_adjust(left=.03,right=.98,bottom=.25,top=.88,wspace=.50);ax=fig.add_subplot(121,projection='3d');ay=fig.add_subplot(122)
ax.scatter(*(X*1e3).T,s=2,c=gray,alpha=.17,depthshade=False,rasterized=True)
for mask,col in [(rt['high'],red),(rt['low'],blue)]:ax.scatter(*(X[mask]*1e3).T,s=8,c=col,alpha=.8,depthshade=False,rasterized=True)
sample=bid[::7];ax.add_collection3d(Line3DCollection(X[b[sample]]*1e3,colors=blue,alpha=.20,linewidths=.45));samples['Fig1_3D_displayed_alive_bond_ids']=sample+1
outline=np.array([[-2.5,-.2,0],[-1,-.2,0],[-1,.2,0],[-2.5,.2,0],[-2.5,-.2,0]])
ax.plot(*outline.T,c=red,lw=2.3);unit=np.array([[-1,-.2,0],[-.95,-.2,0],[-.95,.2,0],[-1,.2,0],[-1,-.2,0]]);ax.plot(*unit.T,c=orange,lw=3)
for yy in [-.2,.2]:ax.plot([-2.5,2.5],[yy,yy],[0,0],ls='--',c=purple,lw=1.2)
for xx in [-2,-.8,.5,1.8]:
 ax.quiver(xx,0,.5,0,0,.22,color=red,arrow_length_ratio=.3,lw=1.5);ax.quiver(xx,0,-.5,0,0,-.22,color=blue,arrow_length_ratio=.3,lw=1.5)
ax.set(xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)',xlim=(-2.6,2.6),ylim=(-.24,.24),zlim=(-.8,.8),xticks=[-2.5,0,2.5],yticks=[-.2,.2],zticks=[-.5,0,.5]);ax.view_init(elev=17,azim=-64);ax.set_proj_type('ortho');ax.set_box_aspect((5,1.1,1.6))
ax.xaxis.labelpad=18;ax.yaxis.labelpad=12;ax.zaxis.labelpad=12
for aa in [ax.xaxis,ax.yaxis,ax.zaxis]:aa.pane.set_facecolor((.97,.98,.99,.6));aa._axinfo['grid']['color']=(.45,.52,.58,.6)
ax.set_title('(a) Actual 16,000-point A0 specimen',pad=16)
ay.scatter(X[sheet,0]*1e3,X[sheet,2]*1e3,s=5,c=gray);ay.add_collection(LineCollection(X[b[bid]][:,:,[0,2]]*1e3,colors=blue,linewidths=.35,alpha=.22));ay.plot([-2.5,-1],[0,0],c=red,lw=3);ay.annotate('advance +X',xy=(-.6,0),xytext=(-1,-.28),arrowprops={'arrowstyle':'->','color':orange},color=orange);spatial(ay);ay.set_title('(b) Real slice, Y = −25 µm; no shell')
fig.suptitle('One candidate boundary: uniform normal traction; all lateral DOFs free',fontsize=20,y=.96)
fig.text(.06,.12,'20 Pa outward traction on each face; 40 µN per face (arrow lengths schematic).\n5 × 0.4 × 1 mm (X × Y × Z); Y periodic; dx = 50 µm; δ = 150.5 µm.\nInitial pre-cut: 1.50 mm. Orange column area: 50 µm × 0.4 mm = 2 × 10⁻⁸ m².',fontsize=14)
fig.legend(handles=[Line2D([],[],marker='o',ls='',color=red,label='Top load nodes / initial pre-cut'),Line2D([],[],color=blue,label='Surviving bonds / bottom load'),Line2D([],[],color=purple,ls='--',label='Periodic Y boundaries')],loc='lower center',bbox_to_anchor=(.5,.015),ncol=3,frameon=False)
save(fig,'Fig1_specimen_and_loading')

# Only the initial accepted native state exists. Do not invent event panels.
fig,axs=plt.subplots(2,1,figsize=(15,8),gridspec_kw={'height_ratios':[3,1.2]},layout='constrained');ax=axs[0]
native=np.load(O/'native_initial_A0.npz');xx=native['x'];ax.scatter(xx[sheet,0]*1e3,xx[sheet,2]*1e3,s=12,c=blue);crosssheet=np.flatnonzero(sheetb&cross);ax.add_collection(LineCollection(xx[b[crosssheet]][:,:,[0,2]]*1e3,colors=orange,linewidths=.6,alpha=.8));ax.plot([-2.5,-1],[0,0],lw=3,c=red);spatial(ax);ax.set_title('Only accepted state: native zero-load A0; coordinates unchanged (×1)')
for id_,offset in [(4710,(-.55,-.67)),(4871,(.2,.60))]:
 assert sheet[id_-1]
 pos=X[id_-1,[0,2]]*1e3;ax.annotate(f'atom {id_}\n(reference position)',xy=pos,xytext=pos+offset,arrowprops={'arrowstyle':'-','color':dark},fontsize=11,bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
axs[1].axis('off');axs[1].text(.03,.88,'NO NATIVE FRACTURE EVOLUTION RECORDED',fontsize=20,color=red,weight='bold',va='top');axs[1].text(.03,.56,'First break, crack advance and final fractured state: NOT MEASURED.\nDesign stopped at the single elastic precheck; no dynamic case was launched.',fontsize=16,va='top');fig.suptitle('Accepted native state and explicit limits of the record',fontsize=21)
samples['Fig2_displayed_surviving_cross_bond_ids']=crosssheet+1;save(fig,'Fig2_accepted_state_no_fracture_evolution')

fig,axs=plt.subplots(2,1,figsize=(15,9),gridspec_kw={'height_ratios':[1.3,1]},layout='constrained');xc=-2.5+(np.arange(100)+.5)*.05
axs[0].bar(xc,t['column_surviving_cross_bonds'],width=.043,color=blue,label='Surviving cross-plane interactions');axs[0].plot(xc,t['column_all_cross_bonds'],c=gray,lw=1.5,label='Original reference interactions');axs[0].axvspan(-2.5,-1,color=red,alpha=.12,label='Initial pre-cut (excluded from new area)');axs[0].axvline(-1,c=orange,ls='--');axs[0].set(xlim=(-2.5,2.5),ylabel='Unique bond count / X column',title='(a) Every original cross-plane interaction retained in accounting');axs[0].legend(frameon=True,facecolor='white',framealpha=1,fontsize=12,ncol=2,loc='lower left')
im=axs[1].imshow(t['column_Y_surviving_cross_bonds'].T,origin='lower',extent=[-2.5,2.5,-.2,.2],aspect='auto',cmap='Blues',vmin=0);axs[1].axvline(-1,c=orange,ls='--');axs[1].set(xlabel='Reference X (mm)',ylabel='Periodic Y (mm)',title='(b) Surviving interactions in all 8 periodic Y bins');fig.colorbar(im,ax=axs[1],label='Surviving bonds / bin')
fig.suptitle('55,408 cross-plane references = 16,416 initial cuts + 38,992 survivors',fontsize=20);fig.text(.06,-.045,'One connected 16,000-point material body. No new cuts in static checks. New dynamic area: NOT MEASURED.\nFull X-projection bridge counts are preserved separately; long bonds are never reduced to broken-count × dx².',fontsize=13);save(fig,'Fig3_crack_area_and_surviving_topology')

fig,axs=plt.subplots(2,2,figsize=(16,10),layout='constrained');ax=axs[0,0];ax.scatter([0],[0],c=blue,s=85,label='Accepted zero-load state');ax.scatter([load['mean_opening_m']*1e6],[load['force_resultant_top_N']*1e6],c=red,marker='x',s=110,lw=2.5,label='Unaccepted static endpoint');ax.set(xlabel='Mean end opening (µm)',ylabel='Applied force per face (µN)',xlim=(-5,155),ylim=(-2,46),title='(a) Static endpoints only — no connecting curve');ax.legend(frameon=False,loc='upper left')
ax=axs[0,1];q=np.genfromtxt(O/'precheck_iterations.csv',delimiter=',',names=True);ax.semilogy(q['evaluation'],q['residual_force_max_N'],c=blue,lw=1.5);ax.axhline(1e-13,c=red,ls='--',label='Frozen absolute gate');ax.scatter([q['evaluation'][-1]],[p['max_residual_force_N']],c=red,s=45,zorder=5);ax.set(xlabel='Optimizer evaluation (NOT physical time)',ylabel='Max residual force (N)',title='(b) Static numerical convergence');ax.legend(frameon=False)
ax=axs[1,0];ax.bar([0,1],[0,p['U_LPS_J']*1e9],color=[blue,red],width=.55);ax.scatter([0],[0],c=blue,s=50);ax.set(xticks=[0,1],xticklabels=['Accepted initial','Unaccepted endpoint'],ylabel='Independent LPS U (nJ)',ylim=(0,3.1),title='(c) Recoverable energy of saved coordinates');ax.text(1,2.65,f'{p["U_LPS_J"]*1e9:.6f} nJ',ha='center');ax.text(.04,.87,'Fixed pristine m; native force cross-check passed',transform=ax.transAxes,fontsize=12)
ax=axs[1,1];ax.axis('off');ax.set_title('(d) Physical trajectory ledger',loc='left');items=[('Boundary work','NOT MEASURED'),('Dynamic K / damping / contact','NOT MEASURED'),('Signed topology energy','NOT MEASURED'),('Energy balance / lag work','NOT MEASURED'),('Native static K and contact U','0 (static checks only)')]
for yy,(a,z) in zip(np.linspace(.82,.10,5),items):ax.text(.01,yy,a,fontsize=13);ax.text(.64,yy,z,color=red if 'NOT' in z else blue,fontsize=12)
fig.suptitle('No physical loading path: do not interpret static load potential as external work',fontsize=20);save(fig,'Fig4_static_loading_energy_and_missing_ledger')

fig,ax=plt.subplots(figsize=(13,7),layout='constrained');vals=[.07170070833333332,.06969623734233767];ax.barh([3,2],vals,color=[blue,orange],height=.6)
for y,v in zip([3,2],vals):ax.text(v+.002,y,f'{v:.8f}',va='center',fontsize=16)
for y in [1,0]:ax.add_patch(Rectangle((0,y-.3),.082,.6,fill=False,edgecolor=red,lw=1.5,ls='--'));ax.text(.041,y,'UNAVAILABLE',ha='center',va='center',color=red,fontsize=18,weight='bold')
ax.set(yticks=[3,2,1,0],yticklabels=['G0_THEORY\nideal continuum','G0_DISCRETE\nideal point-state','G_SEPARATION_PROTOCOL\nactual protocol','GC_EFFECTIVE\nnatural propagation'],xlabel='Energy per projected crack area (J m⁻²)',xlim=(0,.09),ylim=(-.65,3.7));ax.grid(axis='y',visible=False);ax.set_title('Nominal references retained; no measured actual separation resistance',pad=22);fig.text(.02,-.10,'Outlined boxes mean unavailable, not zero energy. No actual-to-nominal ratio or material Gc is inferred.\nEvidence concerns static admission; it does not establish insufficient material dissipation.',fontsize=14);save(fig,'Fig5_nominal_vs_unavailable_actual_energy')

fig,axs=plt.subplots(2,1,figsize=(15,9),layout='constrained');ax=axs[0];ax.scatter(x[sheet,0]*1e3,x[sheet,2]*1e3,s=10,c=blue);ax.add_collection(LineCollection(x[b[crosssheet]][:,:,[0,2]]*1e3,colors=orange,linewidths=.55,alpha=.7));spatial(ax);ax.set_title('(a) Unaccepted static endpoint — actual coordinates, deformation ×1')
ax=axs[1];rn=list(regions);v=[regions[k]['max_utilization'] for k in rn];ax.bar(rn,v,color=[blue,gray,red,gray,gray]);ax.axhline(.03/.15,c=orange,ls='--',label='Frozen precheck stretch ceiling / s00');ax.axhline(1,c=red,ls=':',label='Native breaking threshold');ax.set(xticks=np.arange(5),xticklabels=['Notch front','Notch mouth','Loading ends','Other X surface','Interior'],ylabel='Maximum stretch / s00',ylim=(0,1.12),title='(b) Saved elastic trend; no measured first-break location');ax.legend(frameon=False,fontsize=12,ncol=2);fig.suptitle('Supplement: one 20 Pa elastic precheck failed the fixed admission gates',fontsize=20);save(fig,'FigS1_unaccepted_precheck_and_spatial_competition')
np.savez_compressed(O/'FIGURE_DISPLAY_SUBSETS.npz',**samples)
dump(O/'FIGURE_DATA_PROVENANCE.json',{'figures':files,'data_sources':{str(z.relative_to(O)):sha(z) for z in [O/'reference_topology.npz',O/'elastic_precheck.npz',O/'topology_audit.npz',O/'native_initial_A0.npz',O/'precheck_iterations.csv']},'exact_display_IDs':'FIGURE_DISPLAY_SUBSETS.npz','no_physics_computed':True,'unaccepted_endpoint_appears_only_as_explicit_static_evidence':True,'no_fracture_event_images_invented':True,'PNG_DPI':240})
captions=['真实点阵、预切口和单次候选边界','唯一接受原生初态：没有断裂演化帧','真实存活跨面键与面积识别范围','静态端点、求解残差及未测的动力账本','保留标称值：实际分离能不可用','补充：未接受预检态与区域伸长竞争']
cards=''.join(f'<section><h2>{html.escape(cap)}</h2><a href="{name}.png">PNG</a> · <a href="{name}.svg">SVG</a><img src="{name}.png" alt="{html.escape(cap)}"></section>' for name,cap in zip([Path(f).stem for f in files if f.endswith('.png')],captions))
intro=f'<h1>Mode-I 分离能：单次预检停止</h1><p>Run {C["run"]}。20 Pa 分布牵引预检：最大伸长 0.138705（门槛 0.03）；最大残余力 1.954×10⁻¹² N（门槛 1×10⁻¹³ N）。正式原生断裂算例 <b>0</b>。<b>GC_EFFECTIVE = UNAVAILABLE</b>，不是零。</p><p>已核对原生零载平衡、固定参考m和独立力；未认证20 Pa端点为平衡态。优化迭代不是物理时间，静态势能不是路径外功。</p>'
page='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Controlled Mode-I review</title><style>body{font-family:system-ui,sans-serif;margin:30px auto;max-width:1400px;padding:0 24px;color:#263746;background:#f7f9fb}h1{color:#2b7099}section{background:white;padding:20px;margin:25px 0;border:1px solid #d6dde4;border-radius:8px}img{width:100%;height:auto}p{line-height:1.7}a{color:#2b7099}</style>'+intro+cards+'</html>'
(F/'OPEN_RESULTS.html').write_text(page); rel='../../../results/'+C['task']+'/'+C['run']+'/';(O/'OPEN_RESULTS.html').write_text(page.replace('href="Fig','href="'+rel+'Fig').replace('src="Fig','src="'+rel+'Fig'))
dump(O/'FIGURE_RESOURCES.json',{'wall_s':time.perf_counter()-start,'CPU_user_s':resource.getrusage(resource.RUSAGE_SELF).ru_utime,'CPU_system_s':resource.getrusage(resource.RUSAGE_SELF).ru_stime,'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'physics_runs':0});print('FIGURES',len(files),flush=True)
