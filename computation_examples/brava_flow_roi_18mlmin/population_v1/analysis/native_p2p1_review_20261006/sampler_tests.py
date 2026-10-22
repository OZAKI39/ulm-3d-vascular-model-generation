"""Phase C tests independent analytic fields, native basis, topology and containment."""
from common import *
from p2_sampler import *

def main():
    assert read(HERE/'P2_REPRESENTATION_FIXTURES.json')['status']=='P2_REPRESENTATION_PASS'
    tol=read(HERE/'PROTOCOL.json')['FP64_roundoff_tolerance']
    native=np.load(HERE/'NATIVE_BASIS_SAMPLES.npz');b=np.load(HERE/'P2_REPRESENTATION_FIXTURES.npz')
    tests=[]
    def check(name,passed,**data):
        tests.append(dict(name=name,status='PASS' if passed else 'FAIL',**data))
    def err(name,x,y,scale=1.):
        value=float(np.max(np.abs(x-y))/scale);check(name,value<tol,scaled_max_error=value)
    N,D,L=basis(native['xi'])
    err('native_TET10_node_order_and_values',N,native['N'])
    err('native_TET10_derivatives',D,native['dN'])
    err('native_TET4_pressure_basis',L,native['P'])
    err('partition_of_unity',N.sum(1),1.)
    err('Kronecker_all_10_nodes',N[:10],np.eye(10))
    err('gradient_sum_zero',D.sum(1),0.)
    err('affine_coordinate_reproduction',N@native['nodes'],native['xi'])
    X=b['B1_nodes'];U=b['B1_exact_u'][:10];p=1.+X@np.array([.3,-.7,.8])
    sampler=P2Sampler(X,np.arange(10)[None],U,p)
    result=sampler.sample(b['B1_points'])
    check('all_nodes_and_4096_internal_inside',bool(result['inside_lumen'].all()))
    for k in ['u','G','curl','strain']:err('quadratic_'+k,result[k],b['B1_exact_'+k],max(1.,np.max(abs(b['B1_exact_'+k]))))
    err('linear_pressure',result['pressure'],1+b['B1_points']@np.array([.3,-.7,.8]))
    # Rigid motion: exact antisymmetric gradient, zero strain, finite curl.
    W=np.array([[0,-.6,.4],[.6,0,-.2],[-.4,.2,0]])
    rigid=P2Sampler(X,np.arange(10)[None],X@W.T+np.array([.1,.3,.5]),p).sample(b['B1_points'])
    err('rigid_velocity',rigid['u'],b['B1_points']@W.T+np.array([.1,.3,.5]))
    err('rigid_gradient',rigid['G'],W)
    err('rigid_strain_zero',rigid['strain'],0.)
    # Two positively oriented tetra share the z=0 face, with one off-face apex each.
    points=np.array([[1.,0,0],[0,1,0],[0,0,1],[0,0,0],[0,0,-1]])
    cells=np.array([[0,1,2,3],[1,0,4,3]])
    XP,CP,unique=elevate_affine(points,cells)
    check('unique_shared_edge_midpoints',len(XP)==5+9 and len(set(CP[0,4:])&set(CP[1,4:]))==3)
    test=P2Sampler(XP,CP,XP**2,XP@np.array([1.,2.,3.]))
    probes=np.array([[.2,.3,0],[.5,.5,0],[0,0,0],[.2,.3,.1],[.2,.3,-.1],[.2,.3,-2.],[np.nan,0,0]])
    sampled=test.sample(probes)
    check('shared_face_edge_vertex_smallest_cell',sampled['element_id'].tolist()==[0,0,0,0,1,-1,-1],ids=sampled['element_id'].tolist())
    err('face_edge_internal_u',sampled['u'][:5],probes[:5]**2)
    check('outside_and_nonfinite_not_extrapolated',bool(np.isnan(sampled['u'][5:]).all()) and not sampled['inside_lumen'][5:].any())
    check('valid_samples_finite',all(np.isfinite(result[k]).all() for k in ['u','G','curl','strain','pressure']))
    # Piecewise gradient jump is retained, not recovered/averaged across shared face.
    velocities=np.zeros_like(XP);velocities[CP[0,2],0]=1.
    discontinuous=P2Sampler(XP,CP,velocities,np.zeros(len(XP)))
    sides=discontinuous.sample([[.2,.3,1e-5],[.2,.3,-1e-5],[.2,.3,0.]])
    check('no_gradient_averaging',abs(sides['G'][0,0,2]-sides['G'][1,0,2])>.9 and sides['element_id'][2]==0)
    # Positive mildly curved geometry fixture; physical coordinate field reproduced exactly.
    XC=native['nodes'].copy();XC[4,2]+=.025
    NN,DD,LL=basis(native['xi']);world=NN@XC
    curved=P2Sampler(XC,np.arange(10)[None],XC.copy(),np.arange(10,dtype=float))
    curved_out=curved.sample(world)
    check('curved_fixture_inside',bool(curved_out['inside_lumen'].all()))
    err('curved_coordinate_identity',curved_out['u'],world)
    err('curved_identity_gradient',curved_out['G'],np.eye(3))
    put('P2_FIELD_SAMPLER_CONTRACT.json',dict(status='VALIDATION_ONLY',implementation='p2_sampler.py',
        dtype='float64',velocity='native TET10, ten nodal DOFs, G_ij=du_i/dx_j',pressure='TET4 corner values at native shared reference point; edge pressure serialization values ignored',
        geometry='Isoparametric TET10; inverse Newton mapping, analytic derivatives; affine midpoints reproduce affine mapping',
        containment='Closed simplex; inverse residual <=128 eps coordinate scale; barycentric tolerance derived from inverse Jacobian and capped at1e-10. No clipping, relocation or extrapolation. Boundary points accepted only to documented FP64 roundoff.',
        tie_break='smallest input element index among all admissible candidates; Bernstein-control-point bounds preserve curved candidates',
        invalid='inside_lumen=false, element_id=-1, field quantities NaN; retained not dropped',
        limits='Validation CPU only. Does not certify curved solver mixed geometry, global mesh injectivity or production particle integration.',
        smoothing=False,recovery=False,cell_to_point_average=False,GPU_kernel=False,tolerance=tol))
    put('P2_FIELD_SAMPLER_TESTS.json',dict(status='PASS' if all(t['status']=='PASS' for t in tests) else 'FAIL',tests=tests,
        native_reference_sha256=sha(HERE/'NATIVE_BASIS_SAMPLES.npz')))
    print([(t['name'],t['status']) for t in tests])
    if any(t['status']=='FAIL' for t in tests):raise SystemExit(1)

if __name__=='__main__':main()
