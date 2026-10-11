from analysis.dual_track_v1.common import *
from analysis.dual_track_v1.measure_b import Measurements
from pd_clot.lammps_fragmentation.y_bifurcation_sources import YSources
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection,Line3DCollection
from matplotlib.lines import Line2D
import pyvista as pv,sys
RED='#a63345';BLUE='#287ca4';TEAL='#24877c';GOLD='#bf8a32';GRAY='#69717a'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':14,'axes.titlesize':16,'axes.labelsize':14,'legend.fontsize':12,'axes.grid':True,'grid.color':'#6f7479','grid.alpha':.50,'svg.fonttype':'none'})
def write(fig,name,svg=True):
 FB.mkdir(parents=True,exist_ok=True)
 for ext in (['png','svg'] if svg else ['png']):fig.savefig(FB/(name+'.'+ext),dpi=230,bbox_inches='tight',facecolor='white')
 plt.close(fig)
def grid3(ax,title,limits=None,aspect=None,elev=18,azim=-65):
 ax.set(xlabel='x (mm)',ylabel='y (mm)',zlabel='z (mm)',title=title);ax.view_init(elev,azim);ax.tick_params(labelsize=11,pad=2);ax.xaxis.labelpad=10;ax.yaxis.labelpad=10;ax.zaxis.labelpad=8
 if limits:ax.set(xlim=limits[0],ylim=limits[1],zlim=limits[2])
 if aspect:ax.set_box_aspect(aspect,zoom=1.03)
 for a in [ax.xaxis,ax.yaxis,ax.zaxis]:a.pane.set_facecolor((.96,.96,.96,1));a._axinfo['grid'].update(color=(.4,.4,.4,.7),linewidth=.8)
def line_segments(poly,dims=[0,2]):
 ls=np.asarray(poly.lines);segs=[];i=0
 while i<len(ls):
  n=ls[i];ids=ls[i+1:i+1+n];segs.extend(np.stack([poly.points[ids[:-1]][:,dims],poly.points[ids[1:]][:,dims]],axis=1));i+=n+1
 return np.array(segs)*1e3

def geometry():
 s=YSources();z=np.load(V3/'V3_MATERIAL_POINTS_AND_BONDS.npz');X=z['X'][:1927]*1e3;points=np.load(B/'geometry_arrays.npz');P=points['points_m']*1e3;T=points['triangles'];tags=points['tags'];audit=json.loads((B/'GEOMETRY_AUDIT.json').read_text());fig=plt.figure(figsize=(18,10));fig.subplots_adjust(top=.88,bottom=.12,wspace=.2);ax=fig.add_subplot(1,2,1,projection='3d');surface=Poly3DCollection(P[T[tags!=5]],facecolors=(.58,.66,.7,.075),edgecolors=(.35,.42,.48,.035),linewidths=.15);ax.add_collection3d(surface);ax.scatter(*X.T,c=RED,s=5,alpha=.95,depthshade=False)
 for xyz,label,col in [([-14,0,0],'INLET',BLUE),([8.6,7,1.5],'OUTLET 1',TEAL),([10.392,-6,0],'OUTLET 2',TEAL)]:ax.text(*xyz,label,color=col,fontsize=12)
 ax.text(-6,0,2,'V3 clot',color=RED,fontsize=13);grid3(ax,'Original native Y lumen and 1,927 clot points',([-15,12],[-8,8],[-2,3]),[26,15,4],elev=30,azim=-73)
 ax=fig.add_subplot(1,2,2,projection='3d');wall=s.faces['WALL'];wP=wall.points*1e3;wT=wall.faces.reshape(-1,4)[:,1:];cent=wP[wT].mean(1);sel=(abs(cent[:,0]+6)<1)&(cent[:,2]>.45)&(abs(cent[:,1])<.9);ax.add_collection3d(Poly3DCollection(wP[wT[sel]],facecolors=(.6,.65,.7,.09),edgecolors=(.4,.45,.5,.3),linewidths=.45));ax.scatter(*X.T,c=RED,s=11,alpha=.8,depthshade=False);rim=P[points['contact_line_nodes']];ef=np.sort(np.concatenate([T[tags==5][:,[0,1]],T[tags==5][:,[1,2]],T[tags==5][:,[2,0]]]),axis=1);ee,counts=np.unique(ef,axis=0,return_counts=True);ax.add_collection3d(Line3DCollection(P[ee[counts==1]],colors=TEAL,lw=1.2,label='Shared CFD contact line'));grid3(ax,'Original V3 point cloud on the +z wall',([-6.85,-5.10],[-.8,.8],[.35,1.45]),[1.8,1.6,1.1],elev=16,azim=-65);ax.legend(loc='lower left',fontsize=11)
 fig.suptitle('B1 | Actual native Y geometry and original V3 morphology',fontsize=21);fig.text(.5,.035,f'CFD-triangulated clot volume = {audit["CFD_triangulated_clot_volume_m3"]*1e9:.6f} mm³; discrete material volume = 0.240875 mm³.\nGrey surfaces are numerical fluid-domain boundaries. Clot body remains the actual material-point representation.',ha='center',fontsize=13);write(fig,'FigB1_native_Y_and_V3_points')
 mesh=pv.read(B/'baseline/SV_MESH/mesh-complete.mesh.vtu');cut=mesh.slice(normal=[0,1,0],origin=[0,0,0]).triangulate();tr=cut.faces.reshape(-1,4)[:,1:];c=cut.points[tr].mean(1);use=(c[:,0]>-.0075)&(c[:,0]<-.0045);tr=tr[use];segs=np.concatenate([cut.points[tr[:,[0,1]]][:,:,[0,2]],cut.points[tr[:,[1,2]]][:,:,[0,2]],cut.points[tr[:,[2,0]]][:,:,[0,2]]])*1e3
 fig,axs=plt.subplots(1,2,figsize=(18,8));fig.subplots_adjust(top=.85,bottom=.16,wspace=.25);ax=axs[0];ax.add_collection(LineCollection(segs,colors='#6a7c88',lw=.25,alpha=.75));sel=abs(X[:,1])<.04;ax.scatter(X[sel,0],X[sel,2],s=12,c=RED,zorder=5);ax.set(xlim=[-7.5,-4.5],ylim=[-1.35,1.35],xlabel='x (mm)',ylabel='z (mm)',title='Actual fluid tetrahedra sliced at y = 0');ax.set_aspect('equal')
 ax=axs[1];surf=pv.read(B/'clot_present_lumen.vtp');old=s.exterior.slice(normal=[1,0,0],origin=[-.006,0,0]);new=surf.slice(normal=[1,0,0],origin=[-.006,0,0]);ax.add_collection(LineCollection(line_segments(old,[1,2]),colors=GRAY,lw=2,label='Original lumen'));ax.add_collection(LineCollection(line_segments(new,[1,2]),colors=TEAL,lw=2,label='Actual clot-present fluid contour'));sel=abs(X[:,0]+6)<.03;ax.scatter(X[sel,1],X[sel,2],s=18,c=RED,label='Actual clot points near slice');ax.set(xlim=[-1.4,1.4],ylim=[-1.4,1.4],xlabel='y (mm)',ylabel='z (mm)',title='x = −6 mm: open bypass channel');ax.set_aspect('equal');ax.legend(loc='lower right',fontsize=10)
 meshq=json.loads((B/'baseline/MESH_QUALITY_AUDIT.json').read_text());fig.suptitle('B2 | Actual fluid mesh and clot bypass section',fontsize=21);fig.text(.5,.035,f'{meshq["points"]:,} nodes; {meshq["tetrahedra"]:,} tetrahedra; minimum signed quality = {meshq["quality"]["q_min"]:.3g}.\nNo inverted/degenerate cells. Contact line shared, internal attachment face excluded; surface geometry fixed for the two meshes.',ha='center',fontsize=13);write(fig,'FigB2_mesh_and_bypass')

def fields(level):
 out=B/level;comp=np.load(out/'COMMON_POSITION_COMPARISON.npz');X=np.load(V3/'V3_MATERIAL_POINTS_AND_BONDS.npz')['X'][:1927]*1e3;traction=json.loads((out/'CFD_CLOT_TRACTION_AUDIT.json').read_text());data=np.load(out/'CLOT_SURFACE_TRACTIONS.npz');result=json.loads((out/'CFD_RESULT.json').read_text())
 # Longitudinal sections at identical physical sampling points; strict masks.
 fig,axs=plt.subplots(3,2,figsize=(17,16));fig.subplots_adjust(top=.92,bottom=.11,wspace=.28,hspace=.40);coords=comp['longitudinal_points_m']*1e3;xx=coords[:,0].reshape(130,151);zz=coords[:,2].reshape(130,151);vels=[np.linalg.norm(comp['longitudinal_'+k+'_u_m_s'],axis=1)*1e3 for k in ['old','new']];prs=[comp['longitudinal_'+k+'_p_Pa'] for k in ['old','new']];vnorm=Normalize(0,max(np.nanmax(v) for v in vels));pnorm=Normalize(min(np.nanmin(p) for p in prs),max(np.nanmax(p) for p in prs))
 for row,values,norm,cmap,unit in [(0,vels,vnorm,'viridis','Speed (mm/s)'),(1,prs,pnorm,'cividis','Gauge pressure (Pa)')]:
  for j,k in enumerate(['old','new']):
   ax=axs[row,j];field=values[j].copy();valid=comp['longitudinal_'+k+'_valid'];field[~valid]=np.nan;plot=ax.pcolormesh(xx,zz,field.reshape(130,151),norm=norm,cmap=cmap,shading='auto',rasterized=True);ax.set(xlabel='x (mm)',ylabel='z (mm)',title=('Original no-clot CFD' if j==0 else 'New clot-present CFD')+' | y = 0');ax.set_aspect('equal');fig.colorbar(plot,ax=ax,pad=.025,label=unit,shrink=.82)
   if j==1:
    ids=abs(X[:,1])<.035;ax.scatter(X[ids,0],X[ids,2],c=RED,s=5,zorder=5)
 for key,col,label in [('old',GRAY,'Original no-clot'),('new',TEAL,'Actual clot-present')]:
  axs[2,0].plot(comp['centerline_m'][:,0]*1e3,np.linalg.norm(comp['centerline_'+key+'_u_m_s'],axis=1)*1e3,c=col,label=label);axs[2,1].plot(comp['centerline_m'][:,0]*1e3,comp['centerline_'+key+'_p_Pa'],c=col,label=label)
 axs[2,0].set(xlabel='x (mm)',ylabel='Speed (mm/s)',title='Same physical parent centerline');axs[2,1].set(xlabel='x (mm)',ylabel='Gauge pressure (Pa)',title='Same physical parent centerline')
 for ax in axs[2]:ax.axvline(-6,c=RED,ls=':',label='Clot center');ax.legend(fontsize=11)
 fig.suptitle('B3 | Native velocity and pressure compared at identical physical positions',fontsize=21);fig.text(.5,.03,f'New field: {level} accepted steady solution. Strict P1 interpolation is an analysis operation, with no extrapolation.\nWhite regions are outside that model’s fluid domain; red points belong to the original clot. Pressure is outlet-referenced gauge pressure.',ha='center',fontsize=12);write(fig,'FigB3_common_sections_velocity_pressure')
 # Real exposed computational boundary, distinct from the clot point body.
 fig=plt.figure(figsize=(18,10));fig.subplots_adjust(top=.87,bottom=.15,wspace=.25)
 for j,(key,cmap,title,units) in enumerate([('pressure_Pa','cividis','Pressure on exposed fluid boundary','Gauge pressure (Pa)'),('tangential_magnitude_Pa','magma','Tangential shear on exposed fluid boundary','Shear traction (Pa)')]):
  ax=fig.add_subplot(1,2,j+1,projection='3d');v=data[key];norm=Normalize(v.min(),v.max());colors=plt.get_cmap(cmap)(norm(v));poly=Poly3DCollection(data['triangles_xyz_m']*1e3,facecolors=colors,edgecolors=(.15,.15,.15,.1),linewidths=.10);ax.add_collection3d(poly);ax.scatter(*X.T,c=RED,s=4,alpha=.08,depthshade=False);grid3(ax,title,([-6.75,-5.18],[-.75,.75],[.50,1.36]),[1.6,1.5,.9],elev=16,azim=-68);fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=ax,shrink=.72,pad=.10,label=units,orientation='horizontal',fraction=.045)
 fig.suptitle('B4 | Actual CFD exposed boundary: pressure and shear are separated',fontsize=21);fig.text(.5,.03,f'Area-weighted mean / 95th percentile: pressure {traction["pressure_Pa"]["mean"]:.5g} / {traction["pressure_Pa"]["p95"]:.5g} Pa; shear {traction["shear_Pa"]["mean"]:.5g} / {traction["shear_Pa"]["p95"]:.5g} Pa.\nShown triangles are the numerical fluid boundary; internal attachment face is absent. Adjacent-tet P1 gradients, no smoothing.\nContact-line peak stresses and the one-refinement sensitivity are reported separately; this is not a medical surface reconstruction.',ha='center',fontsize=12);write(fig,'FigB4_exposed_pressure_and_shear')
 # True field streamlines, no animation and no dynamics claim.
 m=Measurements(out/'SV_MESH/mesh_arrays.npz',2.6354060534101725e-8,.005);u,p=m.read(ROOT/result['snapshot']);mesh=pv.read(out/'SV_MESH/mesh-complete.mesh.vtu');mesh.point_data['Velocity']=u;mesh.point_data['Pressure']=p;seeds=np.array([[-.008,y,z] for y in np.linspace(-.0006,.0006,5) for z in [-.0005,0,.0004,.0008] if y*y+z*z<.0011**2]);stream=mesh.streamlines_from_source(pv.PolyData(seeds),vectors='Velocity',integration_direction='forward',initial_step_length=.05,max_steps=800,max_length=.006,terminal_speed=1e-8,compute_vorticity=False);stream.save(out/'native_field_streamlines.vtp');lines=stream.lines;segments=[];speed=[];arrows=[];k=0
 while k<len(lines):
  n=lines[k];ids=lines[k+1:k+n+1];xyz=stream.points[ids];vel=stream['Velocity'][ids];segments.extend(np.stack([xyz[:-1],xyz[1:]],axis=1)*1e3);speed.extend(.5*(np.linalg.norm(vel[:-1],axis=1)+np.linalg.norm(vel[1:],axis=1))*1e3)
  if n>5:
   j=n//3;arrows.append((xyz[j]*1e3,vel[j]/np.linalg.norm(vel[j])))
  k+=n+1
 fig=plt.figure(figsize=(18,9));fig.subplots_adjust(top=.86,bottom=.22);ax=fig.add_subplot(111,projection='3d');norm=Normalize(0,max(speed));lc=Line3DCollection(segments,cmap='viridis',norm=norm,linewidths=1.1);lc.set_array(np.array(speed));ax.add_collection3d(lc);ax.scatter(*X.T,c=RED,s=11,alpha=.8,depthshade=False)
 for pt,direction in arrows:ax.quiver(*pt,*direction,length=.20,color=BLUE,linewidth=.9,arrow_length_ratio=.5)
 grid3(ax,'Forward streamlines of the accepted stationary CFD field',([-8.1,-2.0],[-1.05,1.05],[-1.15,1.45]),[6.1,2.1,2.6],elev=13,azim=-76);ax.set_box_aspect([6.1,2.1,2.6],zoom=1.15);ax.set_yticks([-1,0,1]);fig.colorbar(lc,ax=ax,shrink=.50,pad=.12,fraction=.035,orientation='horizontal',label='Speed (mm/s)');fig.suptitle('B6 | Actual clot points, local streamlines and flow direction',fontsize=21);fig.text(.5,.025,'Streamlines integrate the native P1 velocity field; arrow direction follows that field.\nRigid geometry and steady water flow: these are not microbubble streaming, detached-fragment trajectories or two-way FSI.',ha='center',fontsize=13);write(fig,'FigB6_local_streamlines_and_direction')
 # Transverse bypass profiles at the same coordinates, plus traction directions.
 fig,axs=plt.subplots(1,3,figsize=(19,7));fig.subplots_adjust(top=.84,bottom=.18,wspace=.30);coords=comp['section1_points_m']*1e3;yy=coords[:,1].reshape(111,111);zz=coords[:,2].reshape(111,111);norm=Normalize(0,max(np.nanmax(np.linalg.norm(comp['section1_'+key+'_u_m_s'],axis=1)*1e3) for key in ['old','new']))
 for j,key in enumerate(['old','new']):
  vals=np.linalg.norm(comp['section1_'+key+'_u_m_s'],axis=1)*1e3;vals[~comp['section1_'+key+'_valid']]=np.nan;im=axs[j].pcolormesh(yy,zz,vals.reshape(111,111),norm=norm,cmap='viridis',shading='auto',rasterized=True);axs[j].set(xlabel='y (mm)',ylabel='z (mm)',title=('Original no-clot' if j==0 else 'Actual clot-present')+' | x = −6 mm');axs[j].set_aspect('equal')
  if j==1:
   ids=abs(X[:,0]+6)<.03;axs[j].scatter(X[ids,1],X[ids,2],c=RED,s=10)
 fig.colorbar(im,ax=axs[:2].tolist(),shrink=.75,label='Speed (mm/s)',pad=.025)
 ax=axs[2];cent=data['surface_coordinates_m']*1e3;sv=data['tangential_traction_Pa'];ids=np.arange(0,len(cent),50);ax.scatter(cent[:,0],cent[:,1],c=data['tangential_magnitude_Pa'],cmap='magma',s=4);ar=ax.quiver(cent[ids,0],cent[ids,1],sv[ids,0],sv[ids,1],color=BLUE,angles='xy',scale_units='xy',scale=.20,width=.0045);ax.quiverkey(ar,.10,.08,.02,'0.02 Pa',labelpos='E',coordinates='axes',fontproperties={'size':10});ax.set(xlabel='x (mm)',ylabel='y (mm)',title='Actual shear direction, xy projection');ax.set_aspect('equal')
 fig.suptitle('B7 | Cross-section bypass and tangential surface-force direction',fontsize=21);fig.text(.5,.035,'Both transverse panels use the same physical sampling locations and color range. Main flow travels in +x (out of these transverse panels).\nRight: exposed-face centroids and unnormalized tangential traction vectors; pressure is not included in these arrows.',ha='center',fontsize=12);write(fig,'FigB7_transverse_flow_and_shear_direction')

def summary():
 base=json.loads((B/'baseline/CFD_RESULT.json').read_text());levels=['baseline']+(['refined'] if (B/'refined/CFD_RESULT.json').exists() and json.loads((B/'refined/CFD_RESULT.json').read_text())['status']=='PASS_WITH_SCOPE' else []);s=YSources();old=s.states['simvascular/steady']['metrics'];fig,ax=plt.subplots(2,3,figsize=(20,12));fig.subplots_adjust(top=.88,bottom=.14,hspace=.4,wspace=.55)
 for level,col in zip(levels,[BLUE,TEAL]):
  progress=B/level/('ADAPTER_PROGRESS.json' if level=='baseline' else 'PROGRESS.json');d=json.loads(progress.read_text());st=d['states'];it=d['intervals'];ax[0,0].semilogy([v['step'] for v in st],[max(v['epsilon_mass'],1e-16) for v in st],c=col,label=level);ax[0,1].semilogy([v['step'] for v in it],[v['E_u'] for v in it],c=col,label=level+' velocity');ax[0,1].semilogy([v['step'] for v in it],[max(v['E_Q'],1e-16) for v in it],c=col,ls='--',label=level+' flow')
 ax[0,0].axhline(1e-5,c='black',ls='--');ax[0,0].set(xlabel='Native time step',ylabel='|Q out − Q in| / Q target',title='Actual mass balance');ax[0,0].legend();ax[0,1].axhline(1e-5,c='gray',ls=':');ax[0,1].axhline(1e-6,c='black',ls=':');ax[0,1].set(xlabel='Native time step',ylabel='Relative change per 10 steps',title='Five-interval steady acceptance');ax[0,1].legend(fontsize=10)
 labels=['Old no-clot']+levels;metrics=[old]+[json.loads((B/l/'CFD_RESULT.json').read_text())['metrics'] for l in levels];xx=np.arange(len(labels))
 ax[0,2].bar(xx-.16,[m['outlet_flows_m3_s']['OUTLET_A']*6e7 for m in metrics],.30,color=BLUE,label='Outlet 1');ax[0,2].bar(xx+.16,[m['outlet_flows_m3_s']['OUTLET_B']*6e7 for m in metrics],.30,color=TEAL,label='Outlet 2');ax[0,2].set(xticks=xx,xticklabels=labels,ylabel='Outlet flow (mL/min)',title='Fixed inlet flow, actual outlet split');ax[0,2].legend()
 comps=[json.loads((B/l/'CFD_BASELINE_COMPARISON.json').read_text()) for l in levels];drop=[comps[0]['parent_pressure_drop'][0]['drop_Pa']]+[v['parent_pressure_drop'][1]['drop_Pa'] for v in comps];ax[1,0].bar(xx,drop,color=[GRAY,BLUE,TEAL][:len(labels)]);ax[1,0].set(xticks=xx,xticklabels=labels,ylabel='Area-mean pressure drop (Pa)',title='Parent: x = −9 to −3 mm')
 for j,l in enumerate(levels):
  tr=json.loads((B/l/'CFD_CLOT_TRACTION_AUDIT.json').read_text());ax[1,1].bar(np.arange(3)+(j-.5)*.3,np.array(tr['total_force_N'])*1e6,.28,color=[BLUE,TEAL][j],label=l)
 ax[1,1].set(xticks=range(3),xticklabels=['F x','F y','F z'],ylabel='Integrated fluid force (µN)',title='Pressure + viscous force on clot');ax[1,1].legend()
 sens=B/'MESH_SENSITIVITY.json'
 if sens.exists():
  v=json.loads(sens.read_text());names=['Pressure drop','Mean pressure','Mean shear','Shear P95','Peak speed','Force vector'];errors=list(v['relative_errors'].values());ax[1,2].barh(names,100*np.array(errors),color=TEAL);ax[1,2].axvline(15,c='black',ls='--');ax[1,2].tick_params(axis='y',labelsize=12);ax[1,2].set(xlabel='Refined vs baseline change (%)',title='One local volume refinement')
 else:ax[1,2].axis('off');ax[1,2].text(0,1,'LOCAL REFINEMENT\nNOT ACCEPTED / NOT RUN\n\nSee explicit status and stop evidence.',va='top')
 fig.suptitle('B5 | Native steady gates, pressure drop, integrated force and mesh comparison',fontsize=21);fig.text(.5,.035,'One baseline and at most one refinement; original unobstructed CFD reused. All tractions use only the exposed clot boundary.\nA two-mesh engineering screen is separate from contact-line peak convergence, clinical validity, streaming and two-way FSI.',ha='center',fontsize=12);write(fig,'FigB5_flow_force_and_mesh_audit')
if __name__=='__main__':
 if sys.argv[1]=='geometry':geometry()
 elif sys.argv[1]=='summary':summary()
 else:fields(sys.argv[1])
