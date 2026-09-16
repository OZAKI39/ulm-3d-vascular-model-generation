"""QA figures and ASCII VTP from actual CSV states, exact frozen STL, no smoothing."""
from pathlib import Path
import sys,json,hashlib,csv,xml.etree.ElementTree as ET
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from scipy.integrate import solve_ivp
import vtk
from vtk.util.numpy_support import vtk_to_numpy
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'scripts'));sys.path.insert(0,str(R/'src'))
from validate_cases import table,xyz,qvals
from reference_rigid_sphere_lubrication import coefficients,dense_solve
V=R/'visualization';V.mkdir(exist_ok=True);provenance={}
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.2,'savefig.dpi':180})
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def register(name,inputs):provenance[name]={'sha256':sha(V/name),'source_files':[{'path':str(p.relative_to(R)),'sha256':sha(p)} for p in inputs],'role':'NUMERICAL_QA_ONLY','HUMAN_VISUAL_REVIEW':'PENDING','smoothing':'NONE'}
def figsave(fig,name,inputs):fig.tight_layout();fig.savefig(V/name);plt.close(fig);register(name,inputs)
def paths(*names):return [R/'cases'/n/'TRAJECTORIES.csv' for n in names]
# Dimensional coefficients normalized for a readable comparison; actual contract cutoff/floor.
fig,axs=plt.subplots(1,3,figsize=(13,4));z=np.geomspace(.001,.199,200)
for a,b,label in [(.375e-6,.375e-6,'dmin/dmin'),(.375e-6,2.625e-6,'dmin/dmax'),(.9683592065545495e-6,1.5212443438914028e-6,'d50/d90')]:
 re=a*b/(a+b);c=coefficients(a,b,z*re)
 for ax,k,scale in zip(axs,[0,1,2],[6*np.pi*.001*re,6*np.pi*.001*re,8*np.pi*.001*re**3]):ax.loglog(z,c[:,k]/scale,label=label);ax.set_xlabel('gap / effective radius')
for ax,title in zip(axs,['Normal / (6 pi mu Reff)','Tangential / (6 pi mu Reff)','Pump / (8 pi mu Reff^3)']):ax.set_title(title);ax.legend(fontsize=8)
fig.suptitle('Frozen near-field resistance modes; axial twist remains pending');figsave(fig,'VIS_RESISTANCE_MODES.png',[R/'contracts/RIGID_SPHERE_LUBRICATION_REFERENCE_CONTRACT.json',R/'src/reference_rigid_sphere_lubrication.py'])
# Counterfactual passive and lubrication-only curves are clearly independent reference curves.
A=table(paths('A_normal_cg04')[0]);a=np.array(json.loads((R/'cases/A_normal_cg04/CASE_CONTRACT.json').read_text())['radii_m']);x0=xyz(A[A['step']==0]);end=float(A['time_s'].max());grad=np.diag([-100.,50.,50.]);pair=[(0,1)]
def fun(t,y):
 x=y.reshape(2,3);U=x@grad.T;return dense_solve(x,a,U,np.zeros((2,3)),np.repeat(grad[None],2,axis=0),pair)[:,:3].reshape(-1)
def collision(t,y):x=y.reshape(2,3);return np.linalg.norm(x[0]-x[1])-sum(a)
collision.terminal=True;collision.direction=-1
sol=solve_ivp(fun,[0,end],x0.reshape(-1),method='DOP853',rtol=1e-10,atol=1e-16,events=collision,max_step=1e-4)
xref=sol.y.T.reshape(-1,2,3);gref=np.linalg.norm(xref[:,0]-xref[:,1],axis=1)-sum(a);tpass=np.linspace(0,np.log(np.linalg.norm(x0[0]-x0[1])/sum(a))/100,120);gpass=np.linalg.norm(x0[0]-x0[1])*np.exp(-100*tpass)-sum(a);aa=A[A['particle_id']==1]
fig,axs=plt.subplots(1,2,figsize=(12,4.5));axs[0].plot(tpass*1e3,gpass*1e9,label='passive analytic: ends at contact');axs[0].plot(sol.t*1e3,gref*1e9,'--',label='lubrication only: Python DOP853');axs[0].plot(aa['time_s']*1e3,aa['nearest_gap_m']*1e9,label='actual lubrication + hard constraint');axs[0].set(xlabel='time (ms)',ylabel='surface gap (nm)');axs[0].legend(fontsize=8)
e=table(R/'cases/A_normal_cg04/PAIR_HYDRODYNAMIC_EVENTS.csv');axs[1].plot(e['time_s']*1e3,e['normal_relative_speed_m_s']*1e6);axs[1].set(xlabel='time (ms)',ylabel='relative normal speed (um/s)');fig.suptitle('A: unequal spheres approaching; geometry constraint continues after contact');figsave(fig,'VIS_NORMAL_APPROACH.png',paths('A_normal_cg04')+[R/'src/reference_rigid_sphere_lubrication.py'])
np.savez_compressed(V/'NORMAL_REFERENCE_CURVES.npz',passive_time=tpass,passive_gap=gpass,lubrication_only_time=sol.t,lubrication_only_gap=gref);register('NORMAL_REFERENCE_CURVES.npz',[R/'src/reference_rigid_sphere_lubrication.py'])
# Tangential slide, center relative position, velocity and active shear coefficient.
B=table(paths('B_tangential')[0]);b1=B[B['particle_id']==1];b2=B[B['particle_id']==2];ev=table(R/'cases/B_tangential/PAIR_HYDRODYNAMIC_EVENTS.csv');fig,axs=plt.subplots(1,3,figsize=(13,4));axs[0].plot(b1['time_s']*1e3,(b1['y_m']-b2['y_m'])*1e6);axs[0].set(ylabel='relative y (um)');axs[1].plot(ev['time_s']*1e3,ev['tangential_relative_speed_m_s']*1e6,label='tangential');axs[1].plot(ev['time_s']*1e3,ev['normal_relative_speed_m_s']*1e6,label='normal');axs[1].set(ylabel='relative center speed (um/s)');axs[1].legend();axs[2].plot(ev['time_s']*1e3,ev['shear_resistance_kg_s']);axs[2].set(ylabel='shear resistance (kg/s)')
for ax in axs:ax.set_xlabel('time (ms)')
fig.suptitle('B: initially pure tangential driving, with subsequent coupled motion');figsave(fig,'VIS_TANGENTIAL_SLIDE.png',paths('B_tangential')+[R/'cases/B_tangential/PAIR_HYDRODYNAMIC_EVENTS.csv'])
fig,axs=plt.subplots(2,2,figsize=(11,7))
for row,n in enumerate(['C1_equal_transverse','C2_opposite_transverse']):
 c=table(paths(n)[0]);ev=table(R/'cases'/n/'PAIR_HYDRODYNAMIC_EVENTS.csv')
 for tag in [1,2]:rr=c[c['particle_id']==tag];axs[row,0].plot(rr['time_s']*1e3,rr['omega_z_rad_s'],label='sphere '+str(tag))
 for key,label in [('torque_iz','sphere 1'),('torque_jz','sphere 2')]:axs[row,1].plot(ev['time_s']*1e3,ev[key],label=label)
 axs[row,0].set(ylabel='omega_z (rad/s)',title=n);axs[row,1].set(ylabel='pair torque_z (N m)')
for ax in axs.flat:ax.set_xlabel('time (ms)');ax.legend()
fig.suptitle('C1/C2: applied transverse torque drives; full 6-DOF response. C3 twist blocked');figsave(fig,'VIS_ROTATIONAL_LUBRICATION.png',paths('C1_equal_transverse','C2_opposite_transverse'))
D=table(paths('D_pass_by')[0]);fig,axs=plt.subplots(1,2,figsize=(11,4.5))
for tag in [1,2]:
 rr=D[D['particle_id']==tag];axs[0].plot(rr['x_m']*1e6,rr['y_m']*1e6,label=f'sphere {tag}');sc=axs[0].scatter(rr['x_m'][::10]*1e6,rr['y_m'][::10]*1e6,c=rr['angular_speed_rad_s'][::10],cmap='viridis',s=12);axs[1].plot(rr['time_s']*1e3,rr['angular_speed_rad_s'],label=f'sphere {tag}')
axs[0].set(xlabel='x (um)',ylabel='y (um)',aspect='equal');fig.colorbar(sc,ax=axs[0],label='angular speed (rad/s)');axs[1].set(xlabel='time (ms)',ylabel='angular speed (rad/s)')
for ax in axs:ax.legend()
fig.suptitle('D: tangential pass-by induces particle rotation');figsave(fig,'VIS_TRANSLATION_ROTATION_COUPLING.png',paths('D_pass_by'))
# Real unmodified lumen and both old/new actual trajectories.
reader=vtk.vtkSTLReader();reader.SetFileName(str(R/'provenance/closed_geometry_m.stl'));reader.Update();mesh=reader.GetOutput();coords=vtk_to_numpy(mesh.GetPoints().GetData()).astype(float);polys=vtk_to_numpy(mesh.GetPolys().GetData()).reshape(-1,4)[:,1:];tri=coords[polys]*1e6
K=table(paths('K_real_8')[0]);old=table(R/'provenance/old_case5/trajectory_rank0.csv');fig,axs=plt.subplots(1,3,figsize=(16,5.5));dims=[(0,1),(0,2),(1,2)]
for ax,(i,j) in zip(axs,dims):
 ax.add_collection(PolyCollection(tri[:,:,[i,j]],facecolors='#b7c1cc',edgecolors='none',alpha=.10,rasterized=True));ax.autoscale_view();ax.set(xlabel='xyz'[i]+' (um)',ylabel='xyz'[j]+' (um)',aspect='equal')
 for tag in range(1,9):
  rr=K[K['particle_id']==tag];oo=old[old['particle_id']==tag];x=xyz(rr)*1e6;xo=xyz(oo)*1e6;line=ax.plot(x[:,i],x[:,j],lw=1.7,label=f'ID {tag}')[0];ax.plot(xo[:,i],xo[:,j],'--',color=line.get_color(),lw=1.2);ax.scatter(xo[-1,i],xo[-1,j],marker='x',c='black',s=22,zorder=4)
 ax.legend(fontsize=6,ncol=2)
fig.suptitle('K: real lumen; solid=new 0.020 s, dashed=old passive, x=old stop at 0.005675914 s');figsave(fig,'VIS_REAL_8BUBBLE_RIGID_HYDRO.png',paths('K_real_8')+[R/'provenance/old_case5/trajectory_rank0.csv',R/'provenance/closed_geometry_m.stl'])
# Exact ASCII arrays, including true offline nearest gap for isolated neighbors beyond raw halo.
gaps=table(R/'validation/CASE_K_NEAREST_PAIR_GAPS.csv');assert np.array_equal(gaps['particle_id'],K['particle_id']) and np.array_equal(gaps['step'],K['step'])
data={'particle_id':K['particle_id'],'diameter_um':K['diameter_um'],'time_s':K['time_s'],'velocity':qvals(K)[:,:3],'angular_velocity':qvals(K)[:,3:],'angular_speed':K['angular_speed_rad_s'],'nearest_gap':gaps['nearest_gap_m'],'active_pair_count':gaps['active_pair_count']};lines=[np.flatnonzero(K['particle_id']==tag) for tag in range(1,9)];coordinates=xyz(K)
root=ET.Element('VTKFile',type='PolyData',version='0.1',byte_order='LittleEndian');poly=ET.SubElement(root,'PolyData');piece=ET.SubElement(poly,'Piece',NumberOfPoints=str(len(K)),NumberOfVerts='0',NumberOfLines='8',NumberOfStrips='0',NumberOfPolys='0');pd=ET.SubElement(piece,'PointData')
for name,values in data.items():
 arr=np.asarray(values);e=ET.SubElement(pd,'DataArray',type='Float64',Name=name,NumberOfComponents=str(arr.shape[1] if arr.ndim==2 else 1),format='ascii');e.text=' '.join(format(float(v),'.17g') for v in arr.ravel())
pts=ET.SubElement(piece,'Points');e=ET.SubElement(pts,'DataArray',type='Float64',NumberOfComponents='3',format='ascii');e.text=' '.join(format(float(v),'.17g') for v in coordinates.ravel());ll=ET.SubElement(piece,'Lines');e=ET.SubElement(ll,'DataArray',type='Int64',Name='connectivity',format='ascii');e.text=' '.join(str(i) for line in lines for i in line);e=ET.SubElement(ll,'DataArray',type='Int64',Name='offsets',format='ascii');e.text=' '.join(str(v) for v in np.cumsum([len(l) for l in lines]));ET.indent(root);path=V/'CASE_K_TRAJECTORIES.vtp';ET.ElementTree(root).write(path,encoding='utf-8',xml_declaration=True)
rr=vtk.vtkXMLPolyDataReader();rr.SetFileName(str(path));rr.Update();out=rr.GetOutput();assert out.GetNumberOfLines()==8 and out.GetNumberOfPoints()==len(K);assert np.array_equal(vtk_to_numpy(out.GetPoints().GetData()),coordinates)
for name,values in data.items():assert np.array_equal(vtk_to_numpy(out.GetPointData().GetArray(name)),values),name
register('CASE_K_TRAJECTORIES.vtp',paths('K_real_8')+[R/'validation/CASE_K_NEAREST_PAIR_GAPS.csv'])
# Single indexed event artifact, while preserving per-case CSV unchanged.
with open(R/'PAIR_HYDRODYNAMIC_EVENTS.csv','w') as f:
 w=csv.writer(f);header=False
 for case in json.loads((R/'contracts/CASE_INDEX.json').read_text()):
  with open(R/'cases'/case['name']/'PAIR_HYDRODYNAMIC_EVENTS.csv') as inp:
   reader=csv.reader(inp);names=next(reader)
   if not header:w.writerow(['case']+names);header=True
   for row in reader:w.writerow([case['name']]+row)
(R/'VISUALIZATION_PROVENANCE.json').write_text(json.dumps({'status':'PASS','HUMAN_VISUAL_REVIEW':'PENDING','VTP_READBACK':'PASS_EXACT_ARRAYS','artifacts':provenance},indent=2)+'\n');print('VISUALIZATION PASS',len(provenance))
