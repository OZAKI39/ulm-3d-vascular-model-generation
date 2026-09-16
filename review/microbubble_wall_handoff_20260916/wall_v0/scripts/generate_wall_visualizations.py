"""Static QA figures and ParaView artifacts. Blocked cases never gain fake motion."""
from pathlib import Path
import sys,json,hashlib
import numpy as np,h5py
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import vtk
from vtk.util.numpy_support import numpy_to_vtk
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'src'));from reference_wall_hydrodynamics_v0 import Geometry
V=R/'visualization';V.mkdir(exist_ok=True);PV=V/'paraview';PV.mkdir(exist_ok=True);outputs=[];sources=set()
plt.rcParams.update({'font.size':9,'figure.dpi':120,'axes.grid':True,'grid.alpha':.2})
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def track(p):sources.add(str(p.relative_to(R)));return p
def arr(p):return np.atleast_1d(np.genfromtxt(track(p),delimiter=',',names=True,dtype=None,encoding='utf-8'))
def figsave(name,fig):
 fig.tight_layout();p=V/(name+'.png');fig.savefig(p,dpi=160);plt.close(fig);outputs.append(p)
def lineaxes(ax):
 for a in np.ravel(ax):a.set_xlabel('h/a');a.set_xscale('log')
def case(name):return arr(R/'cases'/name/'TRAJECTORIES.csv')
def xyz(t):return np.column_stack([t[k] for k in ['x_m','y_m','z_m']])
def field(poly,name,value):
 a=vtk.vtkStringArray();a.SetName(name);a.InsertNextValue(value);poly.GetFieldData().AddArray(a)
def write(poly,name,status):
 field(poly,'STATUS',status);field(poly,'HUMAN_VISUAL_REVIEW','PENDING');p=PV/name;w=vtk.vtkXMLPolyDataWriter();w.SetFileName(str(p));w.SetInputData(poly);assert w.Write()==1;outputs.append(p)
def pointpoly(points,data):
 poly=vtk.vtkPolyData();p=vtk.vtkPoints();p.SetData(numpy_to_vtk(np.ascontiguousarray(points,float),deep=True));poly.SetPoints(p);verts=vtk.vtkCellArray()
 for i in range(len(points)):verts.InsertNextCell(1);verts.InsertCellPoint(i)
 poly.SetVerts(verts)
 for k,v in data.items():a=numpy_to_vtk(np.ascontiguousarray(v),deep=True);a.SetName(k);poly.GetPointData().AddArray(a)
 return poly
with h5py.File(track(R/'tables/RMBW_WALL_RESISTANCE_TABLE_V0.h5')) as f:
 e=f['epsilon'][:];T=f['R_total_scaled'][:];eig=f['M_scaled_eigenvalues'][:,0];rec=f['reciprocity_error'][:];cond=f['condition_number'][:]
fig,ax=plt.subplots(2,1,figsize=(9,7));
for label,i,j in [('TT normal',2,2),('TT parallel',0,0),('RR parallel',3,3),('RR normal',5,5)]:ax[0].loglog(e,T[:,i,j],label=label)
ax[0].set_ylabel('Total resistance / bulk');ax[0].legend(ncol=2);ax[1].semilogx(e,T[:,0,4],label='TR / work-conjugate drag scale');ax[1].set_ylabel('Scaled cross coupling');ax[1].legend()
for a in ax:
 for lo,hi,color in [(.001,.2,'#3dba83'),(.2,5,'#e2c956'),(5,20,'#de806e')]:a.axvspan(lo,hi,alpha=.12,color=color)
 for value in [.001,.2,5,20]:a.axvline(value,color='gray',lw=.6)
 a.set_xlabel('h/a')
fig.suptitle('RMBW offline lookup | near-field / qualified / limited ranges\nRR/TR retain source accuracy limitations');figsave('VIS_WALL_LOOKUP_TABLE',fig)
errors=arr(R/'validation/TABLE_HOLDOUT_ERRORS.csv');fig,ax=plt.subplots(1,2,figsize=(11,4));sel=errors['size']=='d50'
for name in ['action_relative_error','mobility_action_relative_error']:ax[0].loglog(errors['epsilon'][sel],np.maximum(errors[name][sel],1e-17),'.',ms=1,label=name.replace('_relative_error',''))
ax[0].axhline(1e-3,color='r',ls='--',label='action gate');ax[0].set_ylabel('Relative action error');ax[0].legend()
for name in ['normal_TT','parallel_TT','RR_parallel','RR_normal','TR']:ax[1].loglog(errors['epsilon'][sel],np.maximum(errors[name+'_relative_error'][sel],1e-17),'.',ms=1,label=name)
ax[1].set_ylabel('Relative coefficient error (diagnostic)');ax[1].legend();lineaxes(ax);fig.suptitle('5,500 independent holdouts per size; all three sizes archived\nSmall TR coefficients can have larger relative error than action error');figsave('VIS_WALL_LOOKUP_ERROR',fig)
fig,ax=plt.subplots(1,3,figsize=(12,3.8));ax[0].loglog(e,eig);ax[0].set_ylabel('Minimum scaled mobility eigenvalue');ax[1].loglog(e,np.maximum(rec,1e-20));ax[1].axhline(1e-10,color='r',ls='--');ax[1].set_ylabel('Reciprocity error');ax[2].loglog(e,cond);ax[2].set_ylabel('Scaled condition number');lineaxes(ax);figsave('VIS_WALL_MATRIX_PROPERTIES',fig)
fig,ax=plt.subplots(2,1,figsize=(9,6))
for name,label in [('F_no_wall','No wall (intentional penetrating control)'),('F_resistance_only','Resistance only (intentional penetrating control)'),('F_hard_cw2','Resistance + hard wall; Cwall=0.2')]:
 t=case(name);ax[0].plot(t['time_s']*1e3,t['surface_wall_gap_m']*1e6,'o-',ms=2,label=label);ax[1].plot(t['time_s']*1e3,t['vz_m_s']*1e3,'o-',ms=2,label=label)
ax[0].axhline(0,color='k',lw=.7);ax[0].set_ylabel('Wall gap (um)');ax[0].legend(fontsize=8);ax[1].set_ylabel('Normal velocity (mm/s)');ax[1].set_xlabel('Time (ms)');figsave('VIS_FLAT_WALL_NORMAL_APPROACH',fig)
fig,ax=plt.subplots(1,3,figsize=(12,4));t=case('F_hard_cw2');ax[0].plot(t['time_s']*1e3,t['vx_m_s']*1e3);ax[0].set_ylabel('Tangential velocity (mm/s)');ax[1].plot(t['time_s']*1e3,t['omega_y_rad_s']);ax[1].set_ylabel('Omega_y (rad/s)');ax[2].semilogx(np.maximum(t['epsilon'],.001),t['wall_excess_TR_scaled'],'o-');ax[2].set_ylabel('Scaled TR excess');ax[2].set_xlabel('h/a (clamped at 0.001)');ax[0].set_xlabel('Time (ms)');ax[1].set_xlabel('Time (ms)');figsave('VIS_FLAT_WALL_TANGENTIAL_ROTATION',fig)
fig,ax=plt.subplots(1,2,figsize=(11,4))
for cw in [4,2,1]:
 t=case('F_hard_cw'+str(cw));ax[0].plot(t['time_s']*1e3,t['surface_wall_gap_m']*1e9,'o-',ms=2,label='Cwall='+str(cw/10));ax[1].plot(t['time_s']*1e3,t['surface_wall_gap_m']*1e12,'o-',ms=2)
for a in ax:a.axhline(0,color='k',lw=.7);a.set_xlabel('Time (ms)')
ax[0].set_ylabel('Gap (nm)');ax[0].legend();ax[1].set_ylabel('Gap (pm), contact zoom');ax[1].axhline(-1,color='red',ls='--',label='-1e-12 m tolerance');ax[1].set_ylim(-1.5,3);ax[1].legend();figsave('VIS_HARD_WALL_NONOVERLAP',fig)
real=np.load(track(R/'raw/REAL_GEOMETRY_FIRST.npz'));curve=np.load(track(R/'raw/CURVED_GEOMETRY_FIRST.npz'));geo=Geometry(track(R/'provenance/closed_geometry_m.stl'));tri=geo.tri*1e6
fig=plt.figure(figsize=(10,5));ax=fig.add_subplot(121,projection='3d');theta=np.linspace(0,2*np.pi,40);phi=np.linspace(0,np.pi,20);xx=20*np.outer(np.cos(theta),np.sin(phi));yy=20*np.outer(np.sin(theta),np.sin(phi));zz=20*np.outer(np.ones_like(theta),np.cos(phi));ax.plot_wireframe(xx,yy,zz,color='gray',alpha=.25,lw=.4);sel=np.arange(0,10000,150);p=curve['positions'][sel]*1e6;q=curve['cpp_query'][sel,:3]*1e6;n=curve['cpp_query'][sel,3:6];ax.scatter(*p.T,s=4,label='sphere centers');ax.scatter(*q.T,s=4,label='closest wall points');ax.quiver(*q.T,*n.T,length=2,normalize=True,color='red',lw=.5);ax.set(xlabel='x (um)',ylabel='y (um)',zlabel='z (um)');ax.legend(fontsize=7);b=fig.add_subplot(122);b.hist((curve['cpp_query'][:,12]-curve['analytic_distance'])*1e9,bins=40);b.set(xlabel='STL distance - analytic sphere distance (nm)',ylabel='Queries',title='Mesh discretization, separate from same-STL identity');figsave('VIS_CURVED_WALL_QUERY_QA',fig)
q=real['cpp_query'];categories=np.zeros(len(tri),int);seenvalid=np.zeros(len(tri),bool);seeninvalid=seenvalid.copy();seenvalid[q[q[:,20]==1,14].astype(int)]=True;seeninvalid[q[q[:,20]==0,14].astype(int)]=True;categories[seenvalid]=1;categories[seeninvalid]=2;categories[seenvalid&seeninvalid]=3
fig,ax=plt.subplots(1,3,figsize=(14,4.5));colors=np.array(['#d9dfe5','#168348','#d86558','#a25abd']);collection=PolyCollection(tri[:,:,:2],facecolors=colors[categories],edgecolors='none',alpha=.8);ax[0].add_collection(collection);ax[0].autoscale();ax[0].set_aspect('equal');ax[0].set(xlabel='x (um)',ylabel='y (um)',title='Queried surface triangles; gray = not evaluated');ax[1].scatter(q[:,15],q[:,16],c=np.where(q[:,20]>0,'#168348','#d86558'),s=3,alpha=.3);ax[1].axvline(.1,color='k',ls='--');ax[1].axhline(20,color='k',ls='--');ax[1].set(xlabel='RMS / a',ylabel='P95 normal spread (deg)',title='Frozen local-plane gates');ax[2].bar(['valid','invalid'],[int(q[:,20].sum()),int((q[:,20]==0).sum())],color=['#168348','#d86558']);ax[2].set_ylabel('Safe interior query count');ax[2].set_title('4 / 10,000 valid; all d10 at a flat cap');figsave('VIS_REAL_LOCAL_PLANE_VALIDITY',fig)
def context(a,i,j):
 points=tri.reshape(-1,3)[::8];a.scatter(points[:,i],points[:,j],s=.12,color='#a6aeb8',alpha=.2);a.set_aspect('equal');a.set_xlabel('xyz'[i]+' (um)');a.set_ylabel('xyz'[j]+' (um)')
fig,ax=plt.subplots(1,3,figsize=(13,4))
for a,(i,j) in zip(ax,[(0,1),(0,2),(1,2)]):context(a,i,j)
fig.suptitle('Case I BLOCKED: no valid d50 near-wall start found\nNo trajectory and no timestep were generated');figsave('VIS_REAL_SINGLE_BUBBLE_WALL',fig)
fig,ax=plt.subplots(figsize=(9,3));ax.axis('off');ax.text(.5,.6,'CASE I: BLOCKED_NO_VALID_NEAR_WALL_START',ha='center',fontsize=14,color='#ad3f31');ax.text(.5,.33,'Gap / h/a / wall activity time histories are unavailable.\nTarget duration remains 0.02 s; actual duration is 0 s.\nNo threshold, radius, or geometry was changed.',ha='center',fontsize=11);figsave('VIS_REAL_SINGLE_BUBBLE_GAP',fig)
old=arr(R/'provenance/PREVIOUS_REAL_8_TRAJECTORIES.csv');new=case('J_real_8');fig,ax=plt.subplots(1,3,figsize=(13,4.5));palette=plt.get_cmap('tab10')
for a,(i,j) in zip(ax,[(0,1),(0,2),(1,2)]):
 context(a,i,j)
 for k,id in enumerate(np.unique(old['particle_id'])):
  p=xyz(old[old['particle_id']==id])*1e6;a.plot(p[:,i],p[:,j],color=palette(k),lw=1.5);p=xyz(new[new['particle_id']==id])*1e6;a.scatter(p[:,i],p[:,j],c=[palette(k)],marker='x',s=25)
fig.suptitle('Case J: prior no-wall trajectories (lines), new initial positions only (crosses)\nNew wall run blocked at t=0; no 0.02 s comparison is available');figsave('VIS_REAL_8BUBBLE_WALL_COMPARISON',fig)
fig,ax=plt.subplots(2,2,figsize=(10,6));labels=['Wall-active bubbles','Minimum wall gap (um)','Maximum resistance correction','Wall constraint activations']
for a,label in zip(ax.ravel(),labels):a.axis('off');a.text(.5,.8,label,ha='center',fontsize=11);a.text(.5,.4,'NO ACCEPTED REAL TIMESTEPS\nBlocked before the first solve',ha='center',color='#ad3f31')
fig.suptitle('Real wall activity is unavailable; initial geometric eligibility is not runtime activity');figsave('VIS_REAL_WALL_ACTIVITY',fig)
# Mandatory point fields; static/nonexistent velocities explicitly NaN.
def metadata(n):return {'particle_id':np.arange(1,n+1,dtype=np.int32),'time_s':np.full(n,np.nan),'diameter_um':np.full(n,np.nan),'velocity_m_s':np.full((n,3),np.nan),'angular_velocity_rad_s':np.full((n,3),np.nan),'gap_m':np.full(n,np.nan),'epsilon':np.full(n,np.nan),'wall_normal':np.full((n,3),np.nan),'triangle_id':np.full(n,-1,np.int32),'rms_over_a':np.full(n,np.nan),'normal_spread_deg':np.full(n,np.nan),'wall_active':np.full(n,-1,np.int32),'constraint_active':np.full(n,-1,np.int32),'accepted_state':np.zeros(n,np.int32)}
write(pointpoly(np.empty((0,3)),metadata(0)),'CASE_I_TRAJECTORY.vtp','BLOCKED; EMPTY; NO_VALID_START_OR_TRAJECTORY')
m=metadata(len(new));m.update(particle_id=new['particle_id'].astype(np.int32),time_s=np.zeros(len(new)),diameter_um=new['diameter_um'],gap_m=new['surface_wall_gap_m'],epsilon=new['epsilon'],wall_normal=np.column_stack([new['normal_'+s] for s in 'xyz']),triangle_id=new['triangle_id'].astype(np.int32),rms_over_a=new['rms_over_a'],normal_spread_deg=new['normal_spread_deg']);write(pointpoly(xyz(new),m),'CASE_J_TRAJECTORIES.vtp','BLOCKED_AT_T0; INITIAL_POINTS_ONLY; NO_LINES; VELOCITY_NOT_SOLVED')
m=metadata(len(q));m.update(diameter_um=real['radii']*2e6,gap_m=q[:,13],epsilon=q[:,13]/real['radii'],wall_normal=q[:,3:6],triangle_id=q[:,14].astype(np.int32),rms_over_a=q[:,15],normal_spread_deg=q[:,16],local_plane_valid=q[:,20].astype(np.int32));write(pointpoly(real['positions'],m),'LOCAL_PLANE_VALIDITY.vtp','STATIC_QUERY_POINTS; VALIDITY_DEPENDS_ON_RADIUS; NO_TIME_EVOLUTION')
sel=np.arange(0,len(q),20);m2={k:v[sel] for k,v in m.items()};m2['nearest_wall_point_m']=q[sel,:3];write(pointpoly(real['positions'][sel],m2),'WALL_QUERY_NORMALS.vtp','STATIC; GLYPH wall_normal OR USE nearest_wall_point_m')
surface=vtk.vtkPolyData();surface.DeepCopy(geo.poly);a=numpy_to_vtk(categories,deep=True);a.SetName('query_validity_class_0_unseen_1_valid_2_invalid_3_mixed');surface.GetCellData().AddArray(a);write(surface,'REAL_SURFACE_VALIDITY_SAMPLES.vtp','STATIC_NEAREST_TRIANGLE_LABELS_ONLY; UNQUERIED_SURFACE_NOT_INFERRED')
provenance={'status':'GENERATED','human_visual_review':'PENDING','script_sha256':sha(Path(__file__)),'sources':{s:sha(R/s) for s in sorted(sources)},'outputs':{str(p.relative_to(R)):sha(p) for p in outputs},'Case_I':'empty trajectory artifact; no dynamics','Case_J':'initial points only; prior 0.02 s baseline separately overlaid; no new trajectory comparison','surface_validity':'radius-dependent sampled nearest-triangle classification; gray unqueried; no spatial extrapolation','NaN_in_VTP':'Intentional not-applicable physical values for static/blocked cases, not solver nonfinite values'};(R/'VISUALIZATION_PROVENANCE.json').write_text(json.dumps(provenance,indent=2)+'\n');print('VIS_GENERATED',len(outputs),'HUMAN_REVIEW_PENDING')
