"""CFD-only PLC: native Y triangles with the ORIGINAL V3 height law applied locally.
No convex hull, ideal Y reconstruction, smoothing, or production geometry edits.
"""
from analysis.dual_track_v1.common import *
import pyvista as pv
from pd_clot.lammps_fragmentation.y_bifurcation_sources import YSources
from pd_clot.lammps_fragmentation.y_clot_morphology_v3 import profile

def surface():
 s=YSources();c=json.loads((V3/'V3_GEOMETRY_CONFIG.json').read_text());surf=s.exterior;P=np.array(surf.points,float);T=surf.faces.reshape(-1,4)[:,1:].copy();tags=np.array(surf['ModelFaceID'],int);initialN=len(P);counts=[]
 # Conforming edge midpoint refinement in a fixed local wall patch; neighboring
 # triangles share each midpoint. Red/green subdivision remains on original planes.
 for it in range(8):
  center=P[T].mean(1);hot=(tags==1)&(abs(center[:,0]-c['x0_m'])<.0010)&(abs(center[:,1])<.0010)&(center[:,2]>.0004);edges=np.sort(np.concatenate([T[hot][:,[0,1]],T[hot][:,[1,2]],T[hot][:,[2,0]]]),axis=1);edges=np.unique(edges,axis=0);edges=edges[np.linalg.norm(P[edges[:,0]]-P[edges[:,1]],axis=1)>50e-6*(1+1e-10)]
  if not len(edges):break
  pts=P.tolist();emap={tuple(e):len(P)+i for i,e in enumerate(edges)};pts.extend(P[edges].mean(1).tolist());newT=[];newtag=[]
  for tri,tag in zip(T,tags):
   a,b,cc=map(int,tri);m=[emap.get(tuple(sorted(e))) for e in [(a,b),(b,cc),(cc,a)]];n=sum(v is not None for v in m)
   if n==0:children=[[a,b,cc]]
   elif n==3:
    ab,bc,ca=m;children=[[a,ab,ca],[ab,b,bc],[ca,bc,cc],[ab,bc,ca]]
   else:
    k=next(k for k in range(3) if m[k] is not None and (n==1 or m[(k+1)%3] is not None));a,b,cc=[int(tri[(k+j)%3]) for j in range(3)];ab=m[k]
    if n==1:children=[[a,ab,cc],[ab,b,cc]]
    else:
     bc=m[(k+1)%3];children=[[ab,b,bc],[a,ab,cc],[ab,bc,cc]]
   newT.extend(children);newtag.extend([tag]*len(children))
  P=np.array(pts);T=np.array(newT,int);tags=np.array(newtag,int);counts.append(dict(iteration=it,points=len(P),triangles=len(T),split_edges=len(edges)))
 else:raise RuntimeError('Local geometry triangulation refinement budget exceeded')
 q,h=profile(P,c);pts=P.tolist();boundary=set();intersections={};newT=[];newtags=[]
 def inter(a,b):
  key=tuple(sorted((a,b)))
  if key not in intersections:
   t=(1-q[a])/(q[b]-q[a]);pt=P[a]+t*(P[b]-P[a]);j=len(pts);pts.append(pt.tolist());intersections[key]=j;boundary.add(j)
  return intersections[key]
 def clip(tri,inside):
  out=[]
  for i in range(3):
   a=int(tri[i]);b=int(tri[(i+1)%3]);ina=q[a]<1;inb=q[b]<1
   if ina==inside:out.append(a)
   if ina!=inb:out.append(inter(a,b))
  return out
 for tri,tag in zip(T,tags):
  if tag!=1:newT.append(tri);newtags.append(tag);continue
  for inner in [False,True]:
   poly=clip(tri,inner)
   for j in range(1,len(poly)-1):newT.append([poly[0],poly[j],poly[j+1]]);newtags.append(5 if inner else 1)
 P=np.array(pts);T=np.array(newT,int);tags=np.array(newtags,int);original=P.copy();height=profile(P,c)[1];height[np.array(sorted(boundary),int)]=0.;affected=np.unique(T[tags==5]);notclot=np.ones(len(P),bool);notclot[affected]=False;height[notclot]=0.;radius=np.linalg.norm(P[:,1:],axis=1);P[:,1:]*=((radius-height)/radius)[:,None]
 out=pv.PolyData(P,np.c_[np.full(len(T),3),T]);out.cell_data['ModelFaceID']=tags.astype(np.int32);out=out.compute_normals(consistent_normals=True,auto_orient_normals=True,split_vertices=False,inplace=False);T=out.faces.reshape(-1,4)[:,1:];tags=np.array(out['ModelFaceID']);assert out.n_open_edges==0
 edges=np.sort(np.concatenate([T[:,[0,1]],T[:,[1,2]],T[:,[2,0]]]),axis=1);_,ec=np.unique(edges,axis=0,return_counts=True);assert np.all(ec==2);assert len(np.unique(np.sort(T,axis=1),axis=0))==len(T)
 areas=.5*np.linalg.norm(np.cross(P[T[:,1]]-P[T[:,0]],P[T[:,2]]-P[T[:,0]]),axis=1);assert areas.min()>0;assert set(tags)=={1,2,3,4,5}
 # Independent signed-volume difference on identical closed triangulation.
 vol=lambda X:abs(np.einsum('ij,ij->',X[T[:,0]],np.cross(X[T[:,1]],X[T[:,2]]))/6)
 originalV=vol(original);fluidV=vol(P);clotV=originalV-fluidV;oldV=json.loads((V3/'GEOMETRY_SUMMARY.json').read_text())['continuous_mesh_volume_m3'];assert abs(clotV/oldV-1)<.02
 out.save(B/'clot_present_lumen.vtp');np.savez_compressed(B/'geometry_arrays.npz',points_m=P,reference_wall_points_m=original,triangles=T,tags=tags,height_m=height,contact_line_nodes=np.array(sorted(boundary)),clot_surface_nodes=affected)
 top=out.extract_cells(tags==5).extract_surface().triangulate();top.save(B/'clot_fluid_surface.vtp');oldsolid=pv.read(ROOT/'results/y_clot_realistic_morphology_wall_contact_and_flow_transport_v1/paraview/clot_solid.vtp')
 # Actual native lumen membership checks, without a newly constructed cylinder.
 inside=pv.PolyData(P[height>1e-8]).select_enclosed_points(s.exterior,tolerance=1e-8,check_surface=True)['SelectedPoints'];assert np.all(inside)
 from scipy.spatial import cKDTree
 sections=[]
 for x in np.linspace(c['x0_m']-.0007,c['x0_m']+.0007,29):
  # Closed slices integrated after triangulation; original native external surface.
  oldslice=s.exterior.slice(normal=[1,0,0],origin=[x,0,0]).triangulate_contours();newslice=out.slice(normal=[1,0,0],origin=[x,0,0]).triangulate_contours()
  sections.append(dict(x_m=float(x),original_area_m2=float(oldslice.area),fluid_area_m2=float(newslice.area),occlusion_fraction=float(1-newslice.area/oldslice.area)))
 audit=dict(status='VALID_CFD_PLC_PENDING_VOLUME_MESH',source_definition='Original y_clot_morphology_v3.profile + V3_GEOMETRY_CONFIG + radial inward height from actual native WALL triangulation; no point-cloud hull or smoothing',source_sha={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'pd_clot/lammps_fragmentation/y_clot_morphology_v3.py',ROOT/'pd_clot/lammps_fragmentation/y_mural_clot_geometry.py',V3/'V3_GEOMETRY_CONFIG.json',ROOT/'verification/y_bifurcation_clot_placement_visual_integration_v1/SOURCE_ASSET_REGISTER.json']},original_native_exterior_sha=sha(s.path(str(Path(s.states['simvascular/mesh']['mesh']).parent/'mesh-complete.exterior.vtp'))),units='global xyz metres, identity transform',original_V3_continuous_volume_m3=oldV,V3_discrete_volume_m3=1927*DX**3,CFD_triangulated_clot_volume_m3=clotV,relative_to_original_V3_mesh=clotV/oldV-1,raw_analytic_solid_scope='Native-wall-fitted height law; no outside-wall solid component or new clipping into lumen. Historical and new finite surface triangulations have distinct documented volumes.',fluid_volume_m3=fluidV,original_lumen_volume_m3=originalV,removed_volume_identity=originalV-fluidV,clot_surface_area_m2=float(areas[tags==5].sum()),wall_attachment_surface_excluded=True,shared_contact_line_nodes=len(boundary),gap_at_contact_line_m=float(abs(P[list(boundary)]-original[list(boundary)]).max()),closed=True,edge_incidence_two=True,single_connected_surface=int(out.connectivity().cell_data['RegionId'].max()+1)==1,duplicate_triangles=0,minimum_triangle_area_m2=float(areas.min()),surface_points=len(P),surface_triangles=len(T),geometry_refinement_history=counts,max_inward_height_m=float(height.max()),minimum_clot_z_m=float(P[affected,2].min()),cross_sections=sections,tag_names={1:'VESSEL_WALL',2:'INLET',3:'OUTLET_1',4:'OUTLET_2',5:'CLOT_FLUID_SURFACE'},ports_preserved=True,geometry_change='CFD-only finite triangulation of the same frozen analytic law; original assets and point cloud untouched. This same PLC is shared by both volume meshes.')
 save(B/'GEOMETRY_AUDIT.json',audit);save(B/'CLOT_PRESENT_GEOMETRY_AUDIT.json',audit);print(json.dumps({k:v for k,v in audit.items() if k not in ['source_sha','cross_sections','geometry_refinement_history']},indent=2))
if __name__=='__main__':surface()
