from common import *
from p2 import *
import math,xml.etree.ElementTree as ET

def test_native_vtk_global_ids_int32():
    import pyvista as pv
    files=[CASE/'SV_MESH/mesh-complete.mesh.vtu',*sorted((CASE/'SV_MESH/mesh-surfaces').glob('*.vtp'))]
    for path in files:
        grid=pv.read(path)
        assert grid.point_data['GlobalNodeID'].dtype==np.dtype('int32')
        assert grid.cell_data['GlobalElementID'].dtype==np.dtype('int32')

def test_tet4_to_tet10_unique_edge_mapping(elevated):
    t,e=elevated['tetra10'],elevated['edges'];n=len(elevated['points_m'])-len(e)
    assert np.array_equal(np.sort(t[:,:4][:,EDGES],axis=-1),e[t[:,4:]-n])
    assert len(np.unique(e,axis=0))==len(e)

def test_tet10_midpoints_exact(elevated):
    p,e=elevated['points_m'],elevated['edges'];n=len(p)-len(e)
    assert np.array_equal(p[n:],p[e].mean(axis=1))

def test_tri6_uses_volume_edge_nodes(elevated):
    n=len(elevated['points_m'])-len(elevated['edges'])
    for name in ['WALL','INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
        f=elevated[name]
        assert np.array_equal(elevated['edges'][f[:,3:]-n],np.sort(f[:,:3][:,TRI_EDGES],axis=-1))

def test_tet10_local_node_order(synthetic):
    corners,nodes=synthetic
    # Source FE lattice origin-first, then exact pinned nn.cpp permutation.
    base=np.array([[0.,0,0],[1.,0,0],[0,1.,0],[0,0,1.]])
    lattice=np.vstack([base,base[np.array([[0,1],[1,2],[2,0],[0,3],[1,3],[2,3]])].mean(axis=1)])
    assert np.array_equal(nodes,lattice[[1,2,3,0,5,9,8,4,6,7]])
    lam=barycentric(nodes,np.broadcast_to(corners,(10,4,3)))
    assert np.allclose(basis(lam),np.eye(10),atol=1e-15)

def test_tet10_positive_jacobian(elevated):
    p,t=elevated['points_m'],elevated['tetra10'];q,w=tetra_rule15()
    minimum=float('inf')
    for i in range(0,len(t),10000):
        x=p[t[i:i+10000]];J=np.einsum('tai,qaj->tqij',x-x[:,:1],reference_gradient(q))
        minimum=min(minimum,float(np.linalg.det(J).min()))
    assert minimum>0

def test_tet4_tet10_geometry_volume_equivalence(elevated):
    import pyvista as pv
    original=pv.read(BASE/'SV_MESH/mesh-complete.mesh.vtu')
    assert np.array_equal(original.points,elevated['points_m'][:original.n_points])
    assert np.array_equal(np.sort(original.cells.reshape(-1,5)[:,1:],axis=1),np.sort(elevated['tetra4'],axis=1))
    a=baseline_arrays();x=a['points_m'][a['tetra']]
    volume=np.abs(np.linalg.det(x[:,1:]-x[:,:1])).sum()/6
    r=json.loads((REPORT/'data/mesh_elevation_audit.json').read_text())
    assert abs(volume/r['total_volume_P2_m3']-1)<1e-12
    assert np.array_equal(np.sort(a['tetra'],axis=1),np.sort(elevated['tetra4'],axis=1))

def test_boundary_area_equivalence(elevated):
    p=elevated['points_m'];lam=TRI_Q
    g=np.array([[1.,0],[0,1.],[-1.,-1.]])
    dn=np.concatenate([(4*lam-1)[...,None]*g,4*(lam[:,TRI_EDGES[:,0],None]*g[TRI_EDGES[:,1]]+lam[:,TRI_EDGES[:,1],None]*g[TRI_EDGES[:,0]])],axis=1)
    for name in ['WALL','INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
        x=p[elevated[name]];J=np.einsum('tai,qaj->tqij',x-x[:,:1],dn)
        area2=np.linalg.norm(np.cross(J[:,:,:,0],J[:,:,:,1]),axis=-1)@TRI_W/2
        area1=np.linalg.norm(np.cross(x[:,1]-x[:,0],x[:,2]-x[:,0]),axis=1)/2
        assert np.allclose(area2,area1,rtol=2e-12,atol=1e-30)

def test_wall_midnodes_zero_bc(elevated):
    tree=ET.parse(CASE/'run/solver.xml');bc=tree.find(".//Add_BC[@name='WALL']")
    assert bc.find('Type').text=='Dir' and float(bc.find('Value').text)==0
    n=len(elevated['points_m'])-len(elevated['edges'])
    assert np.all(elevated['WALL'][:,3:]>=n)
    profile=np.load(REPORT/'data/tet10_inlet_profile.npz')
    wall=np.isin(profile['node_ids'],np.unique(elevated['WALL']))
    assert np.all(profile['velocity_m_s'][wall]==0)

def test_p2_inlet_flux_matches_target(elevated):
    from p2_measure import P2Measurements
    from native_inlet import native_profile
    m=P2Measurements(CASE);a=native_profile(m);f,av=m.boundary['INLET']
    u=np.zeros_like(m.points);u[a['node_ids']]=a['velocity_m_s']
    uq=np.einsum('qa,tai->tqi',tri_basis(TRI_Q),u[f]);q=-math.fsum(np.einsum('tqi,ti->tq',uq,av)@TRI_W)
    assert abs(q/1.551359160440232e-14-1)<1e-14

def test_solver_config_only_changes_expected_fields():
    old=ET.parse(BASE/'run/solver.xml');new=ET.parse(CASE/'run/solver.xml')
    eq=new.find('Add_equation');th=eq.find('Use_taylor_hood_type_basis');assert th.text=='true';eq.remove(th)
    a=old.find(".//Add_BC[@name='INLET']/Value");b=new.find(".//Add_BC[@name='INLET']/Value")
    assert float(b.text)==-1.551359160440232e-14;b.text=a.text
    def flat(t):return [(x.tag,tuple(sorted(x.attrib.items())),(x.text or '').strip()) for x in t.iter()]
    assert flat(old)==flat(new)
    assert sha(BASE/'run/PETSC_OPTIONS.txt')==sha(CASE/'run/PETSC_OPTIONS.txt')

def test_p2_divergence_free_polynomial(synthetic):
    c,n=synthetic;x,y,z=n.T
    u=np.column_stack([x*x,-2*x*y,np.zeros_like(x)])
    q,w=tetra_rule15();g=gradient_at(q,corner_gradients(c));div=np.einsum('qai,ai->q',g,u)
    assert abs(div).max()<2e-14
    assert not np.allclose(u[4:],u[EDGES].mean(axis=1))

def test_p2_known_divergence_polynomial(synthetic):
    c,n=synthetic;u=n*n;q,w=tetra_rule15();g=gradient_at(q,corner_gradients(c))
    div=np.einsum('qai,ai->q',g,u);xyz=q@c
    assert np.allclose(div,2*xyz.sum(axis=1),atol=3e-15)
    assert abs(w@(div**2)/6-.4)<1e-13 # exact integral of [2(x+y+z)]^2 on unit simplex

def test_p2_gauss_theorem(synthetic):
    c,n=synthetic;u=n*n;q,w=tetra_rule15();D=w@np.einsum('qai,ai->q',gradient_at(q,corner_gradients(c)),u)/6
    flux=0.
    for ids in [[0,1,2],[0,1,3],[0,2,3],[1,2,3]]:
        v=c[ids];av=np.cross(v[1]-v[0],v[2]-v[0])/2
        if av@(v.mean(axis=0)-c.mean(axis=0))<0:av=-av
        xyz=TRI_Q@v;velocity=xyz*xyz;flux+=TRI_W@(velocity@av)
    assert abs(flux-D)<1e-14

def test_p2_section_quadrature_degree2():
    # Unit right triangle z=0: exact x^a y^b moment = a! b!/(a+b+2)!.
    v=np.array([[0.,0,0],[1.,0,0],[0,1.,0]]);xyz=TRI_Q@v
    for a,b in [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]:
        integral=.5*(TRI_W@(xyz[:,0]**a*xyz[:,1]**b))
        assert abs(integral-math.factorial(a)*math.factorial(b)/math.factorial(a+b+2))<1e-15

def test_clipped_p2_gauss_and_first_moment(synthetic):
    from p2_audit import clipped_centroid_volume,polygon_p2_flux
    from particle_3d.flowfield_conservation_diagnosis import clip_polygon,tetra_plane_polygon
    c,n=synthetic;a=.2;center=np.array([a,0.,0.]);normal=np.array([1.,0.,0.])
    u=np.column_stack([n[:,0]**2,np.zeros((10,2))])
    volume,centroid=clipped_centroid_volume(c,center,normal)
    exactD=a*a/2-2*a**3/3+a**4/4
    assert abs(volume-(a-a*a+a**3/3)/2)<1e-15
    assert abs(2*volume*centroid[0]-exactD)<1e-15
    polygon,_=tetra_plane_polygon(c,np.zeros((4,3)),center,normal)
    flux,_=polygon_p2_flux(polygon,normal,c,u)
    for ids in [[0,1,2],[0,1,3],[0,2,3],[1,2,3]]:
        v=c[ids];av=np.cross(v[1]-v[0],v[2]-v[0]);nv=av/np.linalg.norm(av)
        if nv@(v.mean(axis=0)-c.mean(axis=0))<0:nv=-nv
        polygon,_=clip_polygon(v,np.zeros((3,3)),center,normal)
        f,_=polygon_p2_flux(polygon,nv,c,u);flux+=f
    assert abs(flux-exactD)<2e-15

def test_baseline_protection_manifest_unchanged():
    manifest=json.loads((REPORT/'baseline_protection_manifest.json').read_text())
    for rel,item in manifest['baseline_case_relative_files'].items():
        assert sha(BASE/rel)==item['sha256'],rel
