"""Closed geometry -> regular cubic volume points -> fresh horizon bonds."""
from pathlib import Path
import itertools
import numpy as np
from scipy.spatial import cKDTree

def _surface(path):
    import trimesh
    if Path(path).suffix.lower()=='.stl':
        surface=trimesh.load(path,force='mesh',process=True)
    else:
        import pyvista as pv
        if Path(path).suffix.lower()=='.npz':
            z=np.load(path);X=z['X'] if 'X' in z else z['mesh_X'];cells=z['cells'] if 'cells' in z else z['mesh_cells']
            if cells.shape[1] != 4: raise ValueError('NPZ FEM mesh needs tetrahedral cells')
            grid=pv.UnstructuredGrid(np.c_[np.full(len(cells),4),cells].ravel(),np.full(len(cells),10,np.uint8),X)
        else: grid=pv.read(path)
        surf=grid.extract_surface().triangulate()
        surface=trimesh.Trimesh(surf.points,surf.faces.reshape(-1,4)[:,1:],process=True)
    if not surface.is_watertight: raise ValueError('Provide a closed watertight volume boundary')
    return surface

def prepare_geometry(request):
    r=request;dx=float(r['dx_m']);kind=r['geometry_type'];origin=np.asarray(r.get('center_m',[0.,0.,0.]),float)
    if kind=='HEMISPHERE':
        R=r['radius_m'];lo=origin+[-R,-R,0.];hi=origin+[R,R,R];surf=None
    else:
        surf=_surface(r['geometry_path']);lo,hi=surf.bounds;R=None
    # Half-cell offset normal to base; x/y stay on a global dx lattice.
    axes=[np.arange(np.ceil(lo[k]/dx),np.floor(hi[k]/dx)+1)*dx for k in [0,1]]
    axes.append(np.arange(lo[2]+dx/2,hi[2]+1e-14,dx))
    if np.prod([len(a) for a in axes])>2_000_000: raise ValueError('Geometry exceeds the declared V1 preprocessing size')
    X=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1).reshape(-1,3)
    inside=np.linalg.norm(X-origin,axis=1)<=R+1e-14 if R else surf.contains(X)
    X=X[inside];N=len(X)
    if N<2: raise ValueError('Geometry has fewer than two material points')
    lattice=np.rint((X-X.min(0))/dx).astype(int);index={tuple(q):i for i,q in enumerate(lattice)}
    base=X[:,2]<lo[2]+dx*(1+1e-9);surface=np.zeros(N,bool)
    for i,q in enumerate(lattice):
        surface[i]=any(tuple(q+np.asarray(e)) not in index for e in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)])
    surface &= ~base
    normals=np.zeros_like(X)
    if R:
        normals[surface]=(X[surface]-origin)/np.linalg.norm(X[surface]-origin,axis=1)[:,None]
        area=2*np.pi*R*R;analytic_volume=2*np.pi*R**3/3
    else:
        _,nearest=cKDTree(surf.triangles_center).query(X[surface]);normals[surface]=surf.face_normals[nearest]
        above=surf.triangles_center[:,2]>lo[2]+.25*dx;area=float(surf.area_faces[above].sum());analytic_volume=float(abs(surf.volume))
    if not base.any() or not surface.any(): raise ValueError('Empty base or free surface group')
    areas=np.zeros(N);areas[surface]=area/surface.sum()
    bonds=cKDTree(X).query_pairs(r['delta_m']*(1+1e-12),output_type='ndarray');bonds=bonds[np.lexsort((bonds[:,1],bonds[:,0]))]
    volume=np.full(N,dx**3);mass=volume*r['rho_kg_m3']
    info=dict(geometry_type=kind,geometry_id=r['geometry_id'],node_count=N,bond_count=len(bonds),average_bonds_per_node=2*len(bonds)/N,
        total_volume_m3=float(volume.sum()),analytic_volume_m3=analytic_volume,volume_ratio=float(volume.sum()/analytic_volume),
        total_mass_kg=float(mass.sum()),surface_nodes=int(surface.sum()),base_nodes=int(base.sum()),dx_m=dx,delta_m=r['delta_m'],
        equal_cubic_spacing=True,nodal_volume_rule='Full dx^3 per retained lattice point; no partial boundary volume fitting',
        surface_area_rule='Equal weights on free surface nodes summing to curved analytic area (hemisphere) or mapped mesh area',
        wall_interface_model='FIXED_BASE_DEMO',base_definition='Lowest one-dx layer above minimum-z plane',
        validation='PASS' if .5<volume.sum()/analytic_volume<1.5 and np.isclose(mass.sum(),r['rho_kg_m3']*volume.sum(),rtol=1e-12) else 'FAIL')
    return dict(X=X,bonds=bonds,volume=volume,mass=mass,base=base,surface=surface,normals=normals,areas=areas,info=info)
