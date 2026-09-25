"""Generate all figures/VTP from immutable numerical CSV/HDF5/STL, no trajectory smoothing."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import vtk
from vtk.util.numpy_support import vtk_to_numpy,numpy_to_vtk
R=Path(sys.argv[1]);sys.path.insert(0,str(R/'src'))
from finalize_passive_transport_v0 import read,trajectory,xyz,vel,integrity
from reference_passive_transport import ReferenceField,dop853,velocity
assert integrity(R)
out=R/'visualization';out.mkdir(exist_ok=True);ct=read(R/'contracts/PASSIVE_MICROBUBBLE_TRANSPORT_V0_CONTRACT.json');audit=read(R/'LOCAL_NUMERICAL_FINALIZER.json');geo=read(R/'LOCAL_GEOMETRY_AND_REPLAY_AUDIT.json');assert audit['status']==geo['status']=='PASS';label=audit['selected_label'];g=ct['gates'];syn=read(R/'contracts/SYNTHETIC_FIELDS_CONTRACT.json');provenance={}
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.25,'savefig.dpi':180})
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def sources(names,extra=[]):
 files=[R/'contracts/PASSIVE_MICROBUBBLE_TRANSPORT_V0_CONTRACT.json',R/'contracts/PASSIVE_TRANSPORT_TIMESTEP_CONTRACT.json']+[R/'cases'/n/'CASE_CONTRACT.json' for n in names]+[f for n in names for f in sorted((R/'cases'/n).glob('trajectory_rank*.csv'))]+[R/e for e in extra]
 return [{'file':str(p.relative_to(R)),'sha256':sha(p)} for p in dict.fromkeys(files)]
def write_ascii_vtp(path,coordinates,point_data,polylines=None,polygons=None):
 import xml.etree.ElementTree as ET
 polylines=[] if polylines is None else polylines;polygons=[] if polygons is None else polygons
 root=ET.Element('VTKFile',type='PolyData',version='0.1',byte_order='LittleEndian');pd=ET.SubElement(root,'PolyData');piece=ET.SubElement(pd,'Piece',NumberOfPoints=str(len(coordinates)),NumberOfVerts='0',NumberOfLines=str(len(polylines)),NumberOfStrips='0',NumberOfPolys=str(len(polygons)))
 arrays=ET.SubElement(piece,'PointData')
 for name,values in point_data.items():
  arr=ET.SubElement(arrays,'DataArray',type='Float64',Name=name,format='ascii');arr.text=' '.join(format(float(x),'.17g') for x in values)
 ET.SubElement(piece,'CellData');points=ET.SubElement(piece,'Points');arr=ET.SubElement(points,'DataArray',type='Float64',NumberOfComponents='3',format='ascii');arr.text=' '.join(format(float(x),'.17g') for x in np.asarray(coordinates).ravel())
 for tag,cells in [('Lines',polylines),('Polys',polygons)]:
  group=ET.SubElement(piece,tag);connect=ET.SubElement(group,'DataArray',type='Int64',Name='connectivity',format='ascii');connect.text=' '.join(str(int(i)) for line in cells for i in line);offset=ET.SubElement(group,'DataArray',type='Int64',Name='offsets',format='ascii');offset.text=' '.join(str(int(i)) for i in np.cumsum([len(line) for line in cells]))
 ET.indent(root);ET.ElementTree(root).write(path,encoding='utf-8',xml_declaration=True)
 ET.parse(path)
 reader=vtk.vtkXMLPolyDataReader();reader.SetFileName(str(path));reader.Update();output=reader.GetOutput();assert output.GetNumberOfPoints()==len(coordinates);assert np.array_equal(vtk_to_numpy(output.GetPoints().GetData()),coordinates)
 for key,values in point_data.items():assert np.array_equal(vtk_to_numpy(output.GetPointData().GetArray(key)),values)
 assert output.GetNumberOfLines()==len(polylines) and output.GetNumberOfPolys()==len(polygons)
 return output
def register(name,inputs):
 p=out/name;provenance[name]={'output_sha256':sha(p),'sources':inputs,'plot_script':str(Path(__file__).relative_to(R)),'plot_script_sha256':sha(__file__),'trajectory_smoothing':'NONE','HUMAN_VISUAL_REVIEW':'PENDING'}
def savefig(fig,name,inputs):fig.tight_layout();fig.savefig(out/name);plt.close(fig);register(name,inputs)
# Uniform analytic trajectory: exact same frozen samples/timestamps.
n='case0_'+label;t=trajectory(R/'cases'/n);sp=read(R/'cases'/n/'CASE_CONTRACT.json');x=xyz(t);ref=np.array(sp['initial_positions_m'][0])+t['time_s'][:,None]*syn['uniform_velocity_m_s'];err=np.linalg.norm(x-ref,axis=1)
fig,axes=plt.subplots(2,2,figsize=(11,7));
for a,axis in enumerate(axes.flat[:3]):axis.plot(t['time_s']*1e3,x[:,a]*1e6,label='LAMMPS RK2',lw=2);axis.plot(t['time_s']*1e3,ref[:,a]*1e6,'--',label='analytic',lw=1.2);axis.set(xlabel='time (ms)',ylabel='XYZ'[a]+' (um)');axis.legend()
axes.flat[3].plot(t['time_s']*1e3,err);axes.flat[3].axhline(g['uniform_position_m'],ls='--',color='red',label='gate 1e-12 m');axes.flat[3].set(xlabel='time (ms)',ylabel='absolute position error (m)');axes.flat[3].legend();axes.flat[3].set_yscale('symlog',linthresh=1e-18);axes.flat[3].set_ylim(bottom=0);fig.suptitle('Case 0: passive uniform advection')
savefig(fig,'VIS_CASE0_UNIFORM_ANALYTIC.png',sources([n],['fields/UNIFORM_FIELD.h5','contracts/SYNTHETIC_FIELDS_CONTRACT.json']))
# Affine field independent high-accuracy reference, no C++ calls.
n='case1_'+label;t=trajectory(R/'cases'/n);sp=read(R/'cases'/n/'CASE_CONTRACT.json');field=ReferenceField(R/sp['field']);ref,meta=dop853(field,sp['initial_positions_m'][0],t['time_s'],affine=(syn['linear_intercept_m_s'],syn['linear_matrix_per_s']));x=xyz(t);fig,axes=plt.subplots(2,2,figsize=(11,7))
for axis,(a,b) in zip(axes.flat[:3],[(0,1),(0,2),(1,2)]):axis.plot(x[:,a]*1e6,x[:,b]*1e6,label='LAMMPS RK2',lw=2);axis.plot(ref[:,a]*1e6,ref[:,b]*1e6,'--',label='Python DOP853');axis.set(xlabel='XYZ'[a]+' (um)',ylabel='XYZ'[b]+' (um)');axis.legend()
axes.flat[3].plot(t['time_s']*1e6,np.linalg.norm(x-ref,axis=1));axes.flat[3].axhline(g['linear_max_position_m'],ls='--',c='red',label='max trajectory gate');axes.flat[3].set(xlabel='time (us)',ylabel='absolute position error (m)');axes.flat[3].legend();axes.flat[3].set_yscale('symlog',linthresh=1e-18);axes.flat[3].set_ylim(bottom=0);fig.suptitle('Case 1: affine flow, independent DOP853 reference')
savefig(fig,'VIS_CASE1_LINEAR_REFERENCE.png',sources([n],['fields/LINEAR_FIELD.h5','contracts/SYNTHETIC_FIELDS_CONTRACT.json']))
# Three step sizes at common physical times; .125 is a comparison reference, not zero true error.
conv=read(R/'contracts/PASSIVE_TRANSPORT_TIMESTEP_CONTRACT.json');fig,axes=plt.subplots(1,2,figsize=(12,4.8));cadv=np.array([.5,.25,.125]);metrics={}
for case,kind in [(1,'linear'),(2,'real')]:
 series={}
 for s in ['c050','c025','c0125']:
  a=trajectory(R/'cases'/f'case{case}_{s}');base=.5*ct['field_speed_statistics'][kind]['dt_base_s'];k=np.rint(a['time_s']/base).astype(int);valid=np.abs(a['time_s']/base-k)<1e-6;series[s]={int(i):p for i,p in zip(k[valid],xyz(a)[valid])}
 finals=[];maxes=[]
 for s in ['c050','c025','c0125']:
  keys=sorted(set(series[s])&set(series['c0125']));diff=np.array([series[s][i]-series['c0125'][i] for i in keys]);finals.append(float(np.linalg.norm(diff[-1])));maxes.append(float(np.max(np.linalg.norm(diff,axis=1))))
 axes[0].plot(cadv,np.array(finals)/ct['dx_m'],'o-',label=kind);axes[1].plot(cadv,np.array(maxes)/ct['dx_m'],'o-',label=kind);metrics[kind]={'C_adv':cadv.tolist(),'final_difference_vs_0125_m':finals,'max_difference_vs_0125_m':maxes}
for ax,gate,ylabel in zip(axes,[.05,.10],['final position difference / dx','max trajectory difference / dx']):ax.axhline(gate,color='red',ls='--',label='gate');ax.axvline(conv['selected_c_adv'],color='green',ls=':',label='selected production C_adv');ax.set(xlabel='C_adv',ylabel=ylabel,xticks=cadv);ax.legend();ax.set_yscale('symlog',linthresh=1e-9);ax.set_ylim(bottom=0)
fig.suptitle('Timestep convergence at matched times (.125 is reference; zero there is by definition)');savefig(fig,'VIS_TIMESTEP_CONVERGENCE.png',sources([f'case{i}_{s}' for i in [1,2] for s in ['c050','c025','c0125']]))
(out/'TIMESTEP_PLOT_DATA.json').write_text(json.dumps(metrics,indent=2)+'\n');register('TIMESTEP_PLOT_DATA.json',sources([f'case{i}_{s}' for i in [1,2] for s in ['c050','c025','c0125']]))
# Read the unchanged closed lumen. Orthographic projections show complete geometry context.
reader=vtk.vtkSTLReader();reader.SetFileName(str(R/'provenance/closed_geometry_m.stl'));reader.Update();mesh=reader.GetOutput();pts=vtk_to_numpy(mesh.GetPoints().GetData()).astype(float)*1e6;tri=vtk_to_numpy(mesh.GetPolys().GetData()).reshape(-1,4)[:,1:];triangles=pts[tri]
write_ascii_vtp(out/'FROZEN_LUMEN.vtp',vtk_to_numpy(mesh.GetPoints().GetData()).astype(float),{},polygons=tri);register('FROZEN_LUMEN.vtp',sources([],['provenance/closed_geometry_m.stl']))
def geometry_axes(title):
 fig,axes=plt.subplots(1,3,figsize=(15,5.5))
 for ax,(a,b) in zip(axes,[(0,1),(0,2),(1,2)]):
  ax.add_collection(PolyCollection(triangles[:,:,[a,b]],facecolors='#cad0d6',edgecolors='none',alpha=.10,rasterized=True));ax.autoscale_view();ax.set_aspect('equal');ax.set(xlabel='XYZ'[a]+' (um)',ylabel='XYZ'[b]+' (um)')
 fig.suptitle(title);return fig,axes
n='case2_'+label;t=trajectory(R/'cases'/n);x=xyz(t)*1e6;reason=read(R/'cases'/n/'trajectory_state_rank0.json')['termination_reason'];fig,axes=geometry_axes('Case 2: frozen transient field, passive bubble center trajectory')
for ax,(a,b) in zip(axes,[(0,1),(0,2),(1,2)]):ax.plot(x[:,a],x[:,b],color='#1266b1',lw=2.4,label='center trajectory');ax.scatter(x[0,a],x[0,b],c='green',s=35,label='start',zorder=5);ax.scatter(x[-1,a],x[-1,b],c='red',marker='x',s=50,label='end: '+reason,zorder=5);ax.legend(fontsize=7)
inputs=sources([n],['provenance/closed_geometry_m.stl','fields/FROZEN_FLOW_FIELD_V0.h5']);savefig(fig,'VIS_CASE2_REAL_TRAJECTORY.png',inputs)
fig,axes=plt.subplots(2,1,figsize=(9,7),sharex=True);axes[0].plot(t['time_s'],t['speed_m_s'],label='particle speed',lw=2);axes[0].plot(t['time_s'],t['fluid_speed_m_s'],'--',label='sampled fluid speed');axes[0].set(ylabel='speed (m/s)');axes[0].legend();axes[1].plot(t['time_s'],np.linalg.norm(vel(t)-np.column_stack([t[k] for k in ['fluid_ux_m_s','fluid_uy_m_s','fluid_uz_m_s']]),axis=1));axes[1].set(xlabel='time (s)',ylabel='|v_particle - u_fluid| (m/s)');fig.suptitle('Case 2: V = U(X) at every completed state');savefig(fig,'VIS_CASE2_SPEED_TIME.png',inputs)
fig,ax=plt.subplots(figsize=(9,5));ax.plot(t['time_s'],t['surface_wall_gap_um'],label='sphere surface-wall gap');ax.axhline(5*ct['dx_m']*1e6,c='green',ls='--',label='5dx initial target');ax.axhline(3*ct['dx_m']*1e6,c='red',ls='--',label='3dx safety boundary');ax.scatter(t['time_s'][-1],t['surface_wall_gap_um'][-1],marker='x',c='black',label=reason);ax.set(xlabel='time (s)',ylabel='surface-wall gap (um)',title='Case 2: exact STL distance, validated independently');ax.legend();savefig(fig,'VIS_CASE2_WALL_CLEARANCE.png',inputs)
def vtp_trajectories(name,rows,inputs):
 fields=['particle_id','time_s','diameter_um','speed_m_s','fluid_speed_m_s','center_wall_distance_um','surface_wall_gap_um']
 polylines=[np.flatnonzero(rows['particle_id']==pid) for pid in np.unique(rows['particle_id'])]
 write_ascii_vtp(out/name,xyz(rows),{key:rows[key] for key in fields},polylines=polylines)
 register(name,inputs)
vtp_trajectories('CASE2_TRAJECTORY.vtp',t,inputs)
fig,axes=plt.subplots(1,3,figsize=(14,4.5));names=[]
for q,style in [('d10','-'),('d50','--'),('d90',':')]:
 n='case3_'+q+'_'+label;names.append(n);a=trajectory(R/'cases'/n)
 for k,ax in enumerate(axes):ax.plot(a['time_s']*1e6,xyz(a)[:,k]*1e6,style,label=q+f' ({a["diameter_um"][0]:.4f} um)',lw=2);ax.set(xlabel='time (us)',ylabel='XYZ'[k]+' (um)');ax.legend(fontsize=8)
fig.suptitle('Expected in passive Stokes-only overdamped V0');savefig(fig,'VIS_CASE3_DIAMETER_INDEPENDENCE.png',sources(names,['fields/LINEAR_FIELD.h5']))
n='case5_'+label
if (R/'cases'/n/'RUN_METRICS.json').exists():
 t=trajectory(R/'cases'/n);reason=read(R/'cases'/n/'trajectory_state_rank0.json')['termination_reason'];fig,axes=geometry_axes('Case 5: independent passive paths; '+reason)
 for pid in np.unique(t['particle_id']):
  a=t[t['particle_id']==pid];x=xyz(a)*1e6
  for ax,(i,j) in zip(axes,[(0,1),(0,2),(1,2)]):line=ax.plot(x[:,i],x[:,j],lw=2,label=f'ID {int(pid)}')[0];ax.scatter(x[0,i],x[0,j],c=[line.get_color()],s=15);ax.scatter(x[-1,i],x[-1,j],c=[line.get_color()],marker='x',s=25)
 for ax in axes:ax.legend(fontsize=6,ncol=2)
 inputs=sources([n],['provenance/closed_geometry_m.stl','fields/FROZEN_FLOW_FIELD_V0.h5']);savefig(fig,'VIS_CASE5_MULTIBUBBLE_TRAJECTORIES.png',inputs);vtp_trajectories('CASE5_TRAJECTORIES.vtp',t,inputs)
(out/'FROZEN_LUMEN_REFERENCE.json').write_text(json.dumps({'file':'FROZEN_LUMEN.vtp','sha256':sha(out/'FROZEN_LUMEN.vtp'),'original_STL':'../provenance/closed_geometry_m.stl','original_STL_sha256':sha(R/'provenance/closed_geometry_m.stl'),'coordinate_units':'m','source_geometry_modified':False},indent=2)+'\n')
(R/'VISUALIZATION_PROVENANCE.json').write_text(json.dumps({'status':'PASS','HUMAN_VISUAL_REVIEW':'PENDING','artifacts':provenance},indent=2)+'\n');print('VISUALIZATION_GENERATION=PASS',len(provenance))
