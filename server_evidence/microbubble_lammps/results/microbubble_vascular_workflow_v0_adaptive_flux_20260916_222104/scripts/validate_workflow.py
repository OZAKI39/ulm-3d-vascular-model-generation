"""Independent audit validator only. Runtime validation remains NOT_RUN after the required STOP."""
from pathlib import Path
import os,json,hashlib,subprocess
import numpy as np,vtk
from flow_geometry import FrozenSampler,Geometry
S=Path(__file__).resolve().parents[1]
a=json.loads((S/'validation/INLET_FLUX_AUDIT.json').read_text())
f=FrozenSampler(S/'fields/FROZEN_FLOW_FIELD_V0.h5');g=Geometry(S/'geometry/GEOMETRY_ARRAYS.npz');sec=np.load(S/'geometry/INJECTION_SECTION.npz')
# Independently intersect original section triangles with Cartesian interpolation cells.
# The actual plane is slightly tilted, so verify it lies in one z grid interval.
zgrid=(sec['points'][:,2]-f.origin[2])/f.dx
assert np.floor(zgrid).min()==np.floor(zgrid).max()
xy=(sec['triangles'][:,:,:2]-f.origin[:2])/f.dx
lo=np.floor(xy.min(axis=(0,1))).astype(int);hi=np.floor(xy.max(axis=(0,1))).astype(int)

def clip(poly,axis,bound,lower):
    if len(poly)==0:return []
    result=[]
    for A,B in zip(poly,np.roll(poly,-1,axis=0)):
        ina=A[axis]>=bound if lower else A[axis]<=bound
        inb=B[axis]>=bound if lower else B[axis]<=bound
        if ina:result.append(A)
        if ina!=inb:
            t=(bound-A[axis])/(B[axis]-A[axis]);result.append(A+t*(B-A))
    return result

def polygon_area(p):
    if len(p)<3:return 0.
    p=np.asarray(p);p=p-p[0]
    return abs(np.cross(p,np.roll(p,-1,axis=0)).sum())*.5

cell_areas=np.zeros(5)
for x in range(lo[0],hi[0]+1):
    for y in range(lo[1],hi[1]+1):
        area=0.
        for t in xy:
            p=t
            for axis,bound,lower in [(0,x,True),(0,x+1,False),(1,y,True),(1,y+1,False)]:p=clip(p,axis,bound,lower)
            area+=polygon_area(p)
        if not area:continue
        # Orthogonal projected area to true planar area Jacobian.
        area*=f.dx*f.dx/abs(sec['normal'][2])
        point=f.origin+f.dx*np.array([x+.5,y+.5,zgrid.mean()]);point[2]=sec['center'][2]-np.dot(sec['normal'][:2],point[:2]-sec['center'][:2])/sec['normal'][2]
        status,_=f.query([point]);cell_areas[int(status[0])]+=area
assert abs(cell_areas.sum()-float(sec['area']))<1e-8*float(sec['area'])
exact_fraction=float(cell_areas[0]/cell_areas.sum());quad_fraction=a['quadrature'][-1]['valid_area_fraction']
assert abs(exact_fraction-quad_fraction)<.001
assert exact_fraction<.9
# Verify invalid support includes points strictly inside the continuous closed lumen.
raw=np.load(S/'raw/INLET_QUADRATURE_FINE.npz');invalid=np.flatnonzero(raw['status']==2)
selected=invalid[np.linspace(0,len(invalid)-1,128,dtype=int)];implicit=vtk.vtkImplicitPolyDataDistance();implicit.SetInput(g.closed)
dist=np.array([implicit.EvaluateFunction(p) for p in raw['points'][selected]])
assert (dist<0).all()
# Area integrals recomputed from raw point values, excluding invalid statuses explicitly.
valid=raw['status']==0;un=raw['velocity']@(-sec['normal']);qpositive=float(raw['weights'][valid]@np.maximum(un[valid],0));qback=float(raw['weights'][valid]@np.maximum(-un[valid],0))
assert abs(qpositive-a['quadrature'][-1]['partial_Q_positive_m3_s'])<1e-25
assert abs(qback-a['quadrature'][-1]['partial_Q_backflow_m3_s'])<1e-25
assert a['Q_positive_m3_s'] is None and a['INLET_FLUX_INTEGRATION_STATUS']=='BLOCKED'
native=np.loadtxt(S/'raw/NATIVE_SAMPLER_OUTPUT.txt');ref=np.load(S/'raw/SAMPLER_PYTHON_REFERENCE.npz')
assert np.array_equal(native[:,0].astype(int),ref['status']);assert np.max(abs(native[:,1:]-ref['velocity']))<1e-15
# Pure quadrature checks using independent scalar polynomials on a reference triangle.
from flow_geometry import quadrature,refine
tri=np.array([[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]]]);analytic=[]
for level in range(3):
    p,w=quadrature(tri)
    vals=[w.sum(),float(w@p[:,0]),float(w@(p[:,0]**2)),float(w@(p[:,0]*p[:,1]))]
    assert np.allclose(vals,[.5,1/6,1/12,1/24],rtol=1e-13,atol=1e-15)
    analytic.append({'level':level,'integrals_1_x_x2_xy':vals});tri=refine(tri)
for kind in ['DONOR_HASHES_BEFORE','ORIGINAL_HASHES_BEFORE']:
    for path,data in json.loads((S/f'provenance/{kind}.json').read_text()).items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==data['sha256'],path
before=json.loads((S/'provenance/GIT_BEFORE.json').read_text());after={}
for repo in before:
    after[repo]={name:subprocess.check_output(['git','-C',repo]+cmd,env={**os.environ,'GIT_OPTIONAL_LOCKS':'0'},text=True) for name,cmd in [('branch',['branch','--show-current']),('commit',['rev-parse','HEAD']),('status',['status','--porcelain=v1','--untracked-files=all'])]}
assert before==after
(S/'provenance/GIT_AFTER.json').write_text(json.dumps(after,indent=2)+'\n')
for p in (S/'scripts').glob('*.py'):compile(p.read_text(),str(p),'exec')
notrun=['A injector arithmetic','B source/admitted semantics','D inlet size admissibility','E transient congestion','F multi-injection','G low-flux workflow','H longer workflow','I high-flux capacity','J outlet lifecycle','K inlet backflow','L MPI1/MPI4','flux-weighted position histogram','near-wall temporal event recomputation']
res={'PREFLIGHT_VALIDATION_STATUS':'PASS','WORKFLOW_TECHNICAL_STATUS':'BLOCKED','STOP_REASON':a['STOP_REASON'],'independent_cell_clipping_area_m2':cell_areas.tolist(),'independent_valid_area_fraction':exact_fraction,'independent_invalid_area_fraction':1-exact_fraction,'fine_quadrature_valid_fraction':quad_fraction,'absolute_fraction_difference':abs(exact_fraction-quad_fraction),'strictly_interior_but_invalid_points_verified':len(dist),'interior_distance_to_closed_surface_range_m':[float(-dist.max()),float(-dist.min())],'native_python_query_comparison':'PASS','partial_Q_positive_recomputed_m3_s':qpositive,'partial_Q_backflow_recomputed_m3_s':qback,'quadrature_polynomial_integrals':analytic,'prior_stage_96_files_unchanged':True,'original_40_inputs_unchanged':True,'both_git_repositories_unchanged':True,'NOT_RUN':notrun,'no_runtime_PASS_inferred':True}
(S/'validation/PREFLIGHT_VALIDATION.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
