from analysis.dual_track_v1.common import *
from analysis.dual_track_v1.measure_b import Measurements
from pd_clot.lammps_fragmentation.y_bifurcation_sources import YSources
import pyvista as pv,vtk,csv,sys
from scipy.spatial import cKDTree

def weighted(values,area):
 order=np.argsort(values);v=values[order];a=area[order];q=(np.cumsum(a)-a/2)/a.sum();return dict(mean=float(np.dot(values,area)/area.sum()),p05=float(np.interp(.05,q,v)),p50=float(np.interp(.5,q,v)),p95=float(np.interp(.95,q,v)),min=float(values.min()),max=float(values.max()))
class Sampler:
 def __init__(self,grid,u,p):
  self.grid=grid;self.u=u;self.p=p;self.tet=grid.cells_dict[pv.CellType.TETRA];xyz=grid.points[self.tet];self.origin=xyz[:,0];self.inverse=np.linalg.inv((xyz[:,1:]-xyz[:,:1]).transpose(0,2,1));self.loc=vtk.vtkStaticCellLocator();self.loc.SetDataSet(grid);self.loc.BuildLocator()
 def sample(self,pts):
  pts=np.asarray(pts);ids=np.array([self.loc.FindCell(x) for x in pts]);w=np.einsum('nij,nj->ni',self.inverse[np.maximum(ids,0)],pts-self.origin[np.maximum(ids,0)]);w=np.c_[1-w.sum(1),w];valid=(ids>=0)&(w.min(1)>=-1e-10)&(w.max(1)<=1+1e-10)
  for k in np.flatnonzero(~valid):
   cand=vtk.vtkIdList();self.loc.FindCellsWithinBounds(np.c_[pts[k]-1e-14,pts[k]+1e-14].ravel(),cand);ii=np.array([cand.GetId(j) for j in range(cand.GetNumberOfIds())],int)
   if len(ii):
    ww=np.einsum('nij,nj->ni',self.inverse[ii],pts[k]-self.origin[ii]);ww=np.c_[1-ww.sum(1),ww];good=(ww.min(1)>=-1e-10)&(ww.max(1)<=1+1e-10)
    if good.any():j=np.flatnonzero(good)[0];ids[k]=ii[j];w[k]=ww[j];valid[k]=True
  u=np.full((len(pts),3),np.nan);p=np.full(len(pts),np.nan);tet=self.tet[ids[valid]];u[valid]=np.einsum('ni,nij->nj',w[valid],self.u[tet]);p[valid]=np.einsum('ni,ni->n',w[valid],self.p[tet]);return u,p,valid

def process(level):
 out=B/level;r=json.loads((out/'CFD_RESULT.json').read_text());assert r['status']=='PASS_WITH_SCOPE';m=Measurements(out/'SV_MESH/mesh_arrays.npz',2.6354060534101725e-8,.005);u,p=m.read(ROOT/r['snapshot']);a=np.load(out/'SV_MESH/mesh_arrays.npz');P=a['points_m'];T=a['tetra'];faces=a['boundary_triangles'];tags=a['facet_tags'];adj=a['adjacent_tetra'];sel=tags==5;f=faces[sel];ad=adj[sel];xyz=P[f];centers=xyz.mean(1);av=.5*np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]);area=np.linalg.norm(av,axis=1);nf=av/area[:,None];nc=-nf;assert np.all(np.einsum('ij,ij->i',nf,centers-P[T[ad]].mean(1))>0)
 edge=P[T[ad,1:]]-P[T[ad,0]][:,None];du=u[T[ad,1:]]-u[T[ad,0]][:,None];grad=np.linalg.solve(edge,du).transpose(0,2,1);strain=grad+grad.transpose(0,2,1);tv=.001*np.einsum('nij,nj->ni',strain,nc);pm=p[f].mean(1);tp=-pm[:,None]*nc;tvn=np.einsum('ij,ij->i',tv,nc)[:,None]*nc;ts=tv-tvn;total=tp+tv;shear=np.linalg.norm(ts,axis=1)
 assert np.max(abs(np.einsum('ij,ij->i',ts,nc)))<1e-12;assert np.allclose(tp+tvn+ts,total,rtol=1e-13,atol=1e-15)
 geom=np.load(B/'geometry_arrays.npz');gt=geom['triangles'][geom['tags']==5];ee=np.sort(np.concatenate([gt[:,[0,1]],gt[:,[1,2]],gt[:,[2,0]]]),axis=1);edge,counts=np.unique(ee,axis=0,return_counts=True);line=geom['points_m'][edge[counts==1]];direction=line[:,1]-line[:,0];length2=np.sum(direction**2,axis=1);rimdist=np.empty(len(centers))
 for begin in range(0,len(centers),512):
  delta=centers[begin:begin+512,None]-line[None,:,0];tt=np.clip(np.einsum('fkj,kj->fk',delta,direction)/length2,0,1);rimdist[begin:begin+512]=np.sqrt(np.min(np.sum((delta-tt[:,:,None]*direction)**2,axis=2),axis=1))
 rim=rimdist<=100e-6;z=np.load(V3/'V3_MATERIAL_POINTS_AND_BONDS.npz');ref=z['X'][:1927].mean(0)
 bary=np.array([[2/3,1/6,1/6],[1/6,2/3,1/6],[1/6,1/6,2/3]]);qxyz=np.einsum('qk,fkj->fqj',bary,xyz);qp=np.einsum('qk,fk->fq',bary,p[f]);qtp=-qp[:,:,None]*nc[:,None];qt=qtp+tv[:,None];mom=np.sum(np.cross(qxyz-ref,qt)*area[:,None,None]/3,axis=(0,1));fp=np.einsum('i,ij->j',area,tp);fv=np.einsum('i,ij->j',area,tv);ft=np.einsum('i,ij->j',area,total);assert np.allclose(ft,fp+fv,rtol=1e-13,atol=1e-20)
 np.savez_compressed(out/'CLOT_SURFACE_TRACTIONS.npz',surface_coordinates_m=centers,triangles_xyz_m=xyz,triangle_mesh_node_ids=f+1,adjacent_mesh_tetra_ids=ad+1,triangle_area_m2=area,n_clot=nc,n_fluid=nf,pressure_Pa=pm,pressure_traction_Pa=tp,viscous_traction_Pa=tv,viscous_normal_traction_Pa=tvn,tangential_traction_Pa=ts,tangential_magnitude_Pa=shear,total_traction_Pa=total,velocity_gradient_s_inv=grad,quadrature_coordinates_m=qxyz,quadrature_pressure_Pa=qp,quadrature_total_traction_Pa=qt,quadrature_area_m2=np.repeat(area[:,None]/3,3,axis=1),integrated_pressure_force_N=fp,integrated_viscous_force_N=fv,integrated_total_force_N=ft,reference_point_m=ref,integrated_total_torque_Nm=mom,contact_rim_distance_m=rimdist,contact_rim_mask=rim)
 surface=pv.PolyData(xyz.reshape(-1,3),np.c_[np.full(len(f),3),np.arange(len(f)*3).reshape(-1,3)]);
 for key,value in dict(Pressure_Pa=pm,Shear_Pa=shear,Shear_vector_Pa=ts,Total_traction_Pa=total,Clot_normal=nc,Contact_rim=rim.astype(np.uint8)).items():surface.cell_data[key]=value
 surface.save(out/'clot_traction_surface.vtp')
 traction=dict(status='EXTRACTED_FROM_ACCEPTED_NATIVE_FLOW',level=level,faces=len(f),exposed_area_m2=float(area.sum()),pressure_Pa=weighted(pm,area),shear_Pa=weighted(shear,area),pressure_force_N=fp,viscous_force_N=fv,viscous_normal_force_N=np.einsum('i,ij->j',area,tvn),viscous_tangential_force_N=np.einsum('i,ij->j',area,ts),total_force_N=ft,total_force_norm_N=float(np.linalg.norm(ft)),reference_point_m=ref,torque_Nm=mom,normal_check='n_fluid outward by adjacent tet; n_clot=-n_fluid; applied fluid-on-clot traction is sigma@n_clot. Tangential dot normal checked.',pressure_scope='Gauge pressure with original zero-reference Neumann traction outlets, not absolute blood pressure',gradient_scope='P1 adjacent native tetrahedral gradient; no smoothing/recovered-gradient model. Mesh sensitivity checked separately.',contact_rim=dict(distance_method='Exact centroid-to-segment distance on actual CFD contact-line edges',width_m=.0001,faces=int(rim.sum()),area_fraction=float(area[rim].sum()/area.sum()),force_N=np.einsum('i,ij->j',area[rim],total[rim]),traction_magnitude_integral_fraction=float(np.sum(area[rim]*np.linalg.norm(total[rim],axis=1))/np.sum(area*np.linalg.norm(total,axis=1))),shear_mean_Pa=float(np.dot(area[rim],shear[rim])/area[rim].sum()),shear_max_Pa=float(shear[rim].max()),shear_nonrim_max_Pa=float(shear[~rim].max()),contains_global_shear_peak=bool(rim[shear.argmax()])),mesh_sha256=r['mesh_sha256'],snapshot_sha256=r['snapshot_sha256'],data_sha256=sha(out/'CLOT_SURFACE_TRACTIONS.npz'),native_field_precision='float64 velocity and pressure, preserved in full precision; vtkPoints float32 identity mapped to frozen input doubles',no_attachment_face_integrated=True,no_LAMMPS_production_load_written=True)
 save(out/'CFD_CLOT_TRACTION_AUDIT.json',traction)
 s=YSources();oldm=s.states['simvascular/steady']['metrics'];old=Sampler(s.flow,s.u,s.p);grid=pv.read(out/'SV_MESH/mesh-complete.mesh.vtu');new=Sampler(grid,u,p)
 # Shared physical parent centerline and three identical transverse planes.
 line=np.c_[np.linspace(-.0135,-.0005,261),np.zeros((261,2))];uo,po,mo=old.sample(line);un,pn,mn=new.sample(line);assert np.all(mo&mn);data=dict(centerline_m=line,centerline_old_u_m_s=uo,centerline_new_u_m_s=un,centerline_old_p_Pa=po,centerline_new_p_Pa=pn)
 axis=np.linspace(-.00129,.00129,111);yy,zz=np.meshgrid(axis,axis,indexing='xy')
 for i,xpos in enumerate([-.007,-.006,-.005]):
  pts=np.c_[np.full(yy.size,xpos),yy.ravel(),zz.ravel()];uo,po,mo=old.sample(pts);un,pn,mn=new.sample(pts);data.update({f'section{i}_points_m':pts,f'section{i}_old_u_m_s':uo,f'section{i}_new_u_m_s':un,f'section{i}_old_p_Pa':po,f'section{i}_new_p_Pa':pn,f'section{i}_old_valid':mo,f'section{i}_new_valid':mn,f'section{i}_common_valid':mo&mn})
 xx,zz=np.meshgrid(np.linspace(-.0075,-.0045,151),np.linspace(-.00129,.00129,130),indexing='xy');pts=np.c_[xx.ravel(),np.zeros(xx.size),zz.ravel()];uo,po,mo=old.sample(pts);un,pn,mn=new.sample(pts);data.update(longitudinal_points_m=pts,longitudinal_old_u_m_s=uo,longitudinal_new_u_m_s=un,longitudinal_old_p_Pa=po,longitudinal_new_p_Pa=pn,longitudinal_old_valid=mo,longitudinal_new_valid=mn)
 # Cross-section area-weighted pressure drop around parent clot: actual mesh slices.
 slice_data=[]
 for sampler,name in [(old,'old'),(new,'new')]:
  field=sampler.grid.copy();field.point_data['P']=sampler.p;field.point_data['U']=sampler.u;means=[]
  for xpos in [-.009,-.003]:
   cut=field.slice(normal=[1,0,0],origin=[xpos,0,0]).triangulate();tri=cut.faces.reshape(-1,4)[:,1:];q=cut.points[tri];aa=.5*np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1);means.append(float(np.dot(aa,cut['P'][tri].mean(1))/aa.sum()))
  slice_data.append(dict(model=name,x_sections_m=[-.009,-.003],mean_pressure_Pa=means,drop_Pa=means[0]-means[1]))
 np.savez_compressed(out/'COMMON_POSITION_COMPARISON.npz',**data)
 comparison=dict(scope='Read-only strict P1 interpolation at identical physical positions; outside each native domain remains NaN and is never extrapolated. Original unclotted CFD reused, never rerun.',old_field_sha=s.flow_id,new_field_sha=r['snapshot_sha256'],old_metrics=oldm,new_metrics=r['metrics'],parent_pressure_drop=slice_data,old_peak_location_m=s.flow.points[np.linalg.norm(s.u,axis=1).argmax()],new_peak_location_m=P[np.linalg.norm(u,axis=1).argmax()],source_meshes_distinct=True,source_geometry='Frozen native Y; clot-present domain removes original V3 analytic continuous shape')
 save(out/'CFD_BASELINE_COMPARISON.json',comparison);rows=[]
 for label,key in [('Q inlet (m3/s)','Q_in_m3_s'),('Q outlets total (m3/s)','Q_out_total_m3_s'),('Maximum speed (m/s)','velocity_max_m_s')]:rows.append(dict(metric=label,unobstructed=oldm[key],clot_present=r['metrics'][key]))
 rows.append(dict(metric='Parent pressure drop -9 to -3 mm (Pa)',unobstructed=slice_data[0]['drop_Pa'],clot_present=slice_data[1]['drop_Pa']))
 for name in ['OUTLET_A','OUTLET_B']:rows.append(dict(metric=name+' flow fraction',unobstructed=oldm['outlet_fractions'][name],clot_present=r['metrics']['outlet_fractions'][name]))
 with (out/'CFD_BASELINE_COMPARISON.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=['metric','unobstructed','clot_present']);w.writeheader();w.writerows(rows)
 save(out/'TRACTION_INTERFACE.json',dict(status='READ_ONLY_OFFLINE_DATA_READY',NPZ='CLOT_SURFACE_TRACTIONS.npz',data_sha256=sha(out/'CLOT_SURFACE_TRACTIONS.npz'),mesh_sha256=r['mesh_sha256'],snapshot_sha256=r['snapshot_sha256'],surface_tag=5,coordinate='global native xyz metres; identity transform',pressure_and_traction='Pa = N/m2',force='N',torque='N m',normal='solid into fluid',mapping_to_LAMMPS='Not implemented/validated; 1927 points and1561 wall anchors are not CFD surface quadrature nodes',streaming=False,two_way_FSI=False))
 print(level,'postprocessed',len(faces[sel]),'exposed faces',flush=True)

def compare():
 rows={k:json.loads((B/k/'CFD_CLOT_TRACTION_AUDIT.json').read_text()) for k in ['baseline','refined']};comp={k:json.loads((B/k/'CFD_BASELINE_COMPARISON.json').read_text()) for k in rows};metrics={k:json.loads((B/k/'CFD_RESULT.json').read_text())['metrics'] for k in rows};a,b=rows['baseline'],rows['refined'];aa,bb=comp['baseline'],comp['refined'];errors={}
 for name,va,vb in [('mother_pressure_drop',aa['parent_pressure_drop'][1]['drop_Pa'],bb['parent_pressure_drop'][1]['drop_Pa']),('surface_pressure_area_mean',a['pressure_Pa']['mean'],b['pressure_Pa']['mean']),('surface_shear_area_mean',a['shear_Pa']['mean'],b['shear_Pa']['mean']),('surface_shear_area_p95',a['shear_Pa']['p95'],b['shear_Pa']['p95']),('peak_velocity',metrics['baseline']['velocity_max_m_s'],metrics['refined']['velocity_max_m_s'])]:errors[name]=abs(vb-va)/abs(va)
 errors['integrated_force_vector_norm_difference']=float(np.linalg.norm(np.array(b['total_force_N'])-a['total_force_N'])/np.linalg.norm(a['total_force_N']));fractions={n:abs(metrics['refined']['outlet_fractions'][n]-metrics['baseline']['outlet_fractions'][n]) for n in ['OUTLET_A','OUTLET_B']};checks={k:v<=.15 for k,v in errors.items()};checks['outlet_split']=max(fractions.values())<=.01
 save(B/'MESH_SENSITIVITY.json',dict(status='PASS_ENGINEERING_SCREEN' if all(checks.values()) else 'MESH_SENSITIVITY_UNRESOLVED',relative_errors=errors,outlet_fraction_absolute_difference=fractions,checks=checks,thresholds=dict(relative=.15,outlet_fraction_absolute=.01),surface_fixed=True,scope='One local volume refinement, unchanged continuous geometry and pre-refined surface. Engineering two-mesh screen, not asymptotic convergence or contact-line peak convergence.',baseline_peak_shear_Pa=a['shear_Pa']['max'],refined_peak_shear_Pa=b['shear_Pa']['max'],baseline_rim=a['contact_rim'],refined_rim=b['contact_rim']))
if __name__=='__main__':compare() if sys.argv[1]=='compare' else process(sys.argv[1])
