"""One mesh-fitted mural solid, fresh global cubic material points and bonds."""
import itertools
import time
import numpy as np
import pyvista as pv
import vtk
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from .config import dump
from .y_bifurcation_sources import OUT, RESULTS

CONFIG = dict(geometry_type='Y_CYLINDRICAL_MURAL_CLOT', geometry_id='CURVED_MURAL_DOME_V1',
              x0_m=-.006, a_m=.0005625, b_m=.0005625, H_m=.0005625, Rv_m=.0013,
              dx_m=50e-6,delta_m=3.01*50e-6,K_Pa=33333.333333333336,G_Pa=940.,rho_kg_m3=1000.,
              s00_bulk=.15,s00_interface=.140624859375,alpha=0.,anchor_layers=3,
              contact_model='MODIFIED_1P0DX',contact_distance_factor=1.,
              coordinate_frame='Y_SOURCE_GLOBAL_XYZ_SI', global_lattice_origin_m=[0.,0.,0.],
              clot_axial_center_m=[-.006,0.,0.],nominal_wall_footpoint_m=[-.006,0.,.0013],
              nominal_flow_direction=[1.,0.,0.],nominal_wall_to_lumen_direction=[0.,0.,-1.],
              CURVED_INTERFACE_STRENGTH_VALIDATION='NOT_ESTABLISHED',CLOT_PRESENT_CFD_SOLVED='NOT_RUN')


def profile(x, c=CONFIG):
    q=np.atleast_2d(x);theta=np.arctan2(q[:,1],q[:,2])
    s2=((q[:,0]-c['x0_m'])/c['a_m'])**2+(c['Rv_m']*theta/c['b_m'])**2
    return s2,c['H_m']*np.maximum(1-s2,0)**1.5


def locator(surface):
    loc=vtk.vtkStaticCellLocator();loc.SetDataSet(surface);loc.BuildLocator();return loc


def hit(loc,a,b):
    t=vtk.mutable(0.);x=[0.,0.,0.];p=[0.,0.,0.];sub=vtk.mutable(0);cell=vtk.mutable(0)
    ok=loc.IntersectWithLine(a,b,1e-13,t,x,p,sub,cell)
    return bool(ok),np.array(x),int(cell)


def inside(surface, positions):
    return np.asarray(pv.PolyData(np.asarray(positions)).select_enclosed_points(surface,tolerance=1e-8,check_surface=True)['SelectedPoints'],bool)


_DISTANCE_CACHE={}


def distances(surface,positions):
    key=id(surface)
    if key not in _DISTANCE_CACHE:
        imp=vtk.vtkImplicitPolyDataDistance();imp.SetInput(surface)
        _DISTANCE_CACHE[key]=(surface,imp)
    imp=_DISTANCE_CACHE[key][1]
    return np.array([imp.EvaluateFunction(p) for p in positions])


def make_surface(sources,c,profile_fn=profile):
    profile=profile_fn
    wall=sources.faces['WALL'];xyz=wall.cell_centers().points
    selected=(abs(xyz[:,0]-c['x0_m'])<c['a_m']+.0004)&(abs(xyz[:,1])<c['b_m']+.0004)&(xyz[:,2]>.0006)
    patch=wall.extract_cells(selected).extract_surface().triangulate().clean().subdivide(3,subfilter='linear')
    patch['s2']=profile(patch.points,c)[0]
    patch=patch.clip_scalar(scalars='s2',value=1.,invert=True).triangulate().clean()
    edge=patch.extract_feature_edges(boundary_edges=True,feature_edges=False,manifold_edges=False,non_manifold_edges=False)
    ids=cKDTree(patch.points).query(edge.points)[1]
    h=profile(patch.points,c)[1];h[ids]=0.
    p=np.array(patch.points,dtype=float);radius=np.linalg.norm(p[:,1:],axis=1)
    exposed=p.copy();exposed[:,1:]*=((radius-h)/radius)[:,None]
    faces=patch.faces.reshape(-1,4)[:,1:]
    attach=pv.PolyData(p,np.c_[np.full(len(faces),3),faces]).triangulate()
    top=pv.PolyData(exposed,np.c_[np.full(len(faces),3),faces[:,::-1]]).triangulate()
    attach.cell_data['surface_tag']=np.full(len(faces),2,np.int32)
    top.cell_data['surface_tag']=np.full(len(faces),1,np.int32)
    solid=attach.merge(top,merge_points=True,tolerance=1e-12).clean(tolerance=1e-12,absolute=True).triangulate()
    area=solid.compute_cell_sizes()['Area'];solid=solid.extract_cells(area>1e-22).extract_surface().clean().triangulate()
    solid=solid.compute_normals(consistent_normals=True,auto_orient_normals=True,split_vertices=False)
    assert solid.n_open_edges==0,'Mural reference solid must be watertight'
    top=solid.extract_cells(solid['surface_tag']==1).extract_surface().triangulate()
    attach=solid.extract_cells(solid['surface_tag']==2).extract_surface().triangulate()
    loc=locator(wall)
    ok,foot,cell=hit(loc,[c['x0_m'],0.,0.],[c['x0_m'],0.,.002])
    assert ok
    assert np.linalg.norm(foot-np.array(c['nominal_wall_footpoint_m']))<=c['dx_m'],'Wall discrepancy exceeds dx: stop for explicit projection review'
    c.update(actual_wall_footpoint_m=foot.tolist(),footpoint_wall_triangle=int(cell),
             projection_rule='Attachment patch is clipped/subdivided actual WALL triangles; profile height measured radially inward from each actual wall point',
             nominal_footpoint_correction_m=float(np.linalg.norm(foot-np.array(c['nominal_wall_footpoint_m']))))
    return solid,top,attach


def build(sources,c=None,*,out=OUT,results=RESULTS,profile_fn=profile):
    OUT=out;RESULTS=results;profile=profile_fn
    started=time.perf_counter();c=CONFIG.copy() if c is None else dict(c)
    solid,top,attach=make_surface(sources,c,profile_fn=profile_fn)
    dx=c['dx_m'];lo=np.array([c['x0_m']-c['a_m']-dx,-c['b_m']-dx,c['Rv_m']-c['H_m']-2*dx])
    hi=np.array([c['x0_m']+c['a_m']+dx,c['b_m']+dx,c['Rv_m']+4*dx])
    axes=[np.arange(np.floor(lo[k]/dx),np.ceil(hi[k]/dx)+1,dtype=int) for k in range(3)]
    ijk=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1).reshape(-1,3);cand=ijk*dx
    body_mask=inside(solid,cand)
    X=cand[body_mask];initial=len(X)
    # Whole nominal voxel corners must remain in the original lumen. The smooth
    # boundary mesh still represents the actual wall attachment; dx^3 is PD quadrature.
    corners=np.array(list(itertools.product([-.5,0,.5],repeat=3)))*dx
    corner_ids=sources.mesh.find_containing_cell((X[:,None,:]+corners).reshape(-1,3)).reshape(-1,27)
    corner_signed=distances(sources.exterior,(X[:,None,:]+corners).reshape(-1,3)).reshape(-1,27)
    X=X[(corner_ids>=0).all(1)&(corner_signed.max(1)<-1e-9)];conservative_removed=initial-len(X)
    # Remove only disconnected / one-neighbour discrete whiskers, no new shapes.
    removed=0
    while True:
        index={tuple(v):i for i,v in enumerate(np.rint(X/dx).astype(int))}
        degrees=np.array([sum(tuple(q+e) in index for e in np.r_[np.eye(3,dtype=int),-np.eye(3,dtype=int)]) for q in np.rint(X/dx).astype(int)])
        keep=degrees>=2
        if keep.all():break
        removed+=int((~keep).sum());X=X[keep]
    n=len(X);assert n>100
    wall_dist=distances(sources.exterior,X);assert wall_dist.max()<0
    s2,h=profile(cand,c);sd=distances(sources.exterior,cand)
    # Three signed-normal-distance shells, selected on the SAME global lattice.
    amask=(s2<1)&(h>=dx/2)&(sd>1e-10)&(sd<=3*dx+1e-12)
    A=cand[amask];ad=sd[amask];layers=np.minimum(np.ceil(ad/dx).astype(int),3)
    assert len(A)>0 and set(layers)=={1,2,3}
    allx=np.vstack([X,A]);N=len(allx);assert len(np.unique(allx,axis=0))==N
    top_loc=locator(top);solid_loc=locator(solid);attach_loc=locator(attach);wall_loc=locator(sources.exterior)
    candidates=cKDTree(allx).query_pairs(c['delta_m']*(1+1e-12),output_type='ndarray')
    kept=[];types=[];rejected={1:0,2:0,3:0}
    for i,j in candidates:
        kind=1 if j<n else (2 if i<n else 3)
        a,b=allx[i],allx[j]
        if kind==1:
            valid=not hit(solid_loc,a,b)[0]
        elif kind==2:
            valid=hit(attach_loc,a,b)[0] and not hit(top_loc,a,b)[0]
        else:
            # Reject chords cutting through the lumen or leaving the local shell footprint.
            q=a[None,:]+np.linspace(0,1,7)[:,None]*(b-a)
            valid=bool((profile(q,c)[0]<1).all() and (distances(sources.exterior,q)>-1e-12).all() and not hit(wall_loc,a,b)[0])
        if valid:kept.append((i,j));types.append(kind)
        else:rejected[kind]+=1
    bonds=np.asarray(kept,int);types=np.asarray(types,int);order=np.lexsort((bonds[:,1],bonds[:,0]));bonds=bonds[order];types=types[order]
    assert set(types)=={1,2,3}
    bb=bonds[types==1];adj=coo_matrix((np.ones(len(bb)*2),(np.r_[bb[:,0],bb[:,1]],np.r_[bb[:,1],bb[:,0]])),shape=(n,n))
    assert connected_components(adj)[0]==1
    index={tuple(q):i for i,q in enumerate(np.rint(X/dx).astype(int))}
    boundary=np.array([any(tuple(q+e) not in index for e in np.r_[np.eye(3,dtype=int),-np.eye(3,dtype=int)]) for q in np.rint(X/dx).astype(int)])
    attachment_nodes=np.zeros(n,bool);attachment_nodes[bonds[types==2,0]]=True
    # Only actual free-boundary candidates get area. Attachment node group may
    # overlap topologically at edges, but geometric face areas never overlap.
    dtop=np.abs(distances(top,X));dattach=np.abs(wall_dist)
    surface=boundary&(dtop<dattach);assert surface.any()
    exposed_nodes=np.flatnonzero(surface)
    top=top.compute_normals(cell_normals=True,point_normals=False,consistent_normals=True,split_vertices=False)
    cent=top.cell_centers().points;areas=top.compute_cell_sizes()['Area']
    owner=exposed_nodes[cKDTree(X[exposed_nodes]).query(cent)[1]]
    nodal_areas=np.bincount(owner,weights=areas,minlength=n)
    normal=np.stack([np.bincount(owner,weights=areas*top.cell_normals[:,k],minlength=n) for k in range(3)],axis=1)
    nn=np.linalg.norm(normal,axis=1);normal[nn>0]/=nn[nn>0,None]
    surface=nodal_areas>0
    info=dict(geometry_type=c['geometry_type'],geometry_id=c['geometry_id'],node_count=N,clot_nodes=n,anchor_nodes=len(A),
              bond_count=len(bonds),bulk_bonds=int((types==1).sum()),interface_bonds=int((types==2).sum()),wall_bonds=int((types==3).sum()),
              average_bonds_per_node=2*len(bonds)/N,clot_volume_m3=n*dx**3,clot_mass_kg=n*dx**3*c['rho_kg_m3'],
              continuous_mesh_volume_m3=float(solid.volume),old_hemisphere_analytic_volume_m3=2*np.pi*c['H_m']**3/3,
              old_hemisphere_nominal_volume_m3=None,
              bbox_material_centers_m=[X.min(0),X.max(0)],bbox_continuous_m=list(solid.bounds),centroid_material_m=X.mean(0),
              exposed_area_m2=float(top.area),attachment_area_m2=float(attach.area),total_surface_area_m2=float(solid.area),
              footprint_axial_extent_m=float(np.ptp(attach.points[:,0])),footprint_circumferential_arc_extent_m=float(np.ptp(c['Rv_m']*np.arctan2(attach.points[:,1],attach.points[:,2]))),
              maximum_mesh_height_m=float(np.max(profile(attach.points,c)[1])),
              conservative_wall_voxel_removals=conservative_removed,removed_lattice_whiskers=removed,
              rejected_void_crossing_bonds=rejected,wall_layers_counts={str(i):int((layers==i).sum()) for i in [1,2,3]},
              discrete_body_components=1,anchor_mass_in_clot_totals=False,
              nodal_volume_rule='Full dx^3 for retained cubic centers; distinct from boundary-fitted mesh volume. Conservative lumen-corner mask; unresolved thin rim has no particles.',
              surface_area_rule='Actual exposed triangulation area assigned once to nearest free-boundary material node; attachment triangles excluded',
              CURVED_INTERFACE_STRENGTH_VALIDATION='NOT_ESTABLISHED',CLOT_PRESENT_CFD_SOLVED='NOT_RUN')
    c.update(source_geometry_id=sources.geometry_id,source_mesh_id=sources.mesh_id,source_flow_id=sources.flow_id)
    dump(OUT/'PLACEMENT_CONFIGURATION.json',c);dump(OUT/'GEOMETRY_SUMMARY.json',info)
    np.savez_compressed(OUT/'CURVED_MURAL_CLOT_GEOMETRY.npz',X=X,lattice_ijk=np.rint(X/dx).astype(int),candidate_lattice_ijk=ijk,
                        candidate_solid_mask=body_mask,node_id=np.arange(1,n+1),node_volume_m3=np.full(n,dx**3),node_mass_kg=np.full(n,c['rho_kg_m3']*dx**3),
                        exposed_surface=surface,attachment_surface=attachment_nodes,areas_m2=nodal_areas,normals=normal)
    np.savez_compressed(OUT/'CURVED_WALL_ANCHORS.npz',X=A,node_id=np.arange(n+1,N+1),layer=layers,wall_distance_m=ad,fixed=np.ones(len(A),bool),lattice_ijk=np.rint(A/dx).astype(int))
    np.savez_compressed(OUT/'BOND_GROUPS.npz',bonds_zero_based=bonds,bond_type=types,bond_id=np.arange(1,len(bonds)+1),reference_length_m=np.linalg.norm(allx[bonds[:,1]]-allx[bonds[:,0]],axis=1))
    out=RESULTS/'paraview';out.mkdir(exist_ok=True)
    solid.save(out/'clot_solid.vtp');top.save(out/'exposed_surface.vtp');attach.save(out/'attachment_surface.vtp')
    sources.exterior.save(out/'vessel_surface.vtp');sources.faces['WALL'].save(out/'vessel_wall.vtp')
    for name,points in [('clot_points',X),('anchors',A)]:
        grid=pv.PolyData(points).cast_to_unstructured_grid()
        if name=='clot_points':grid['node_id']=np.arange(1,n+1);grid['exposed_surface']=surface.astype(int);grid['attachment_surface']=attachment_nodes.astype(int)
        else:grid['node_id']=np.arange(n+1,N+1);grid['layer']=layers
        grid.save(out/(name+'.vtu'))
    lines=pv.PolyData();lines.points=allx;lines.lines=np.c_[np.full(len(bonds),2),bonds];lines.cell_data['bond_type']=types
    lines.save(out/'all_bonds.vtp')
    g=dict(X=allx,bonds=bonds,bond_types=types,clot_count=n,anchors=np.arange(N)>=n,
           volume=np.full(N,dx**3),mass=np.full(N,c['rho_kg_m3']*dx**3),
           surface=np.r_[surface,np.zeros(len(A),bool)],base=np.r_[attachment_nodes,np.zeros(len(A),bool)],
           areas=np.r_[nodal_areas,np.zeros(len(A))],normals=np.vstack([normal,np.zeros_like(A)]),info=info)
    print('GEOMETRY_BUILT',n,len(A),len(bonds),time.perf_counter()-started,flush=True)
    dump(OUT/'GEOMETRY_BUILD_PERFORMANCE.json',dict(wall_seconds=time.perf_counter()-started,work='One fixed geometry; new cubic mask and bonds; no time integration'))
    return g,c,solid,top,attach
