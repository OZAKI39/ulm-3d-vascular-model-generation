"""Geometry-only coverage/force checks on a known detached fixture; no dynamics."""
from common import *
from current_surface import CurrentSurface
from shear_field_library import ShearField
import time,resource
start=time.perf_counter();cpu=time.process_time()
ref=load_npz(O/'provenance/reference.npz');X=ref['X'];n=int(ref['clot_count']);b=ref['bonds']
graph=module('surface_semantics_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')
upper=X[:,2]>0;d=ref['initial_cut'].astype(float);d[upper[b[:,0]]!=upper[b[:,1]]]=1
x=X.copy();x[upper,2]+=20e-6
lab,sizes,attached=graph.topology(ref,d);free=~attached[lab];assert free.sum()==896
results={}
for mode in CASES:
    surf=CurrentSurface(ref,mode);v=np.zeros_like(X);v[upper]=[.1,.03,0]
    z=surf.evaluate(x,v,d,.002)
    free_faces=z['active']&np.r_[free,np.zeros(len(X)-n,bool)][surf.tri[:,0]]
    assert free_faces.sum()>0 and z['accessible'].sum()>0 and np.linalg.norm(z['force'][:n][free])>0
    assert np.all(lab[surf.tri[free_faces]]==lab[surf.tri[free_faces,0]][:,None])
    off=surf.evaluate(x,v,d,END);assert not off['drive_force'].any() and np.sum(off['force']*v)<0
    far=surf.evaluate(x+[.003,0,0],v,d,.002);assert not far['drive_force'].any() and np.sum(far['force']*v)<0
    # One/two-node fine cannot be assigned any fully supported triangle/tetra.
    pair=np.array([0,1]);df=ref['initial_cut'].astype(float);sel=np.isin(b,pair)
    df[sel.any(1)&~sel.all(1)]=1
    fine=surf.evaluate(X,np.zeros_like(X),df,.002)
    assert not fine['force'][pair].any()
    # Field continuity and analytic gradient, independent central differences.
    pts=np.array([[-.2e-3,-.08e-3,375e-6],[.02e-3,.06e-3,375e-6],[.16e-3,-.03e-3,375e-6]])
    normals=np.tile([0.,0,1.],(3,1));eps=1e-9
    def psi(p):
        a=p[:,0];c=p[:,1]
        if mode=='counter_rotating':
            sigma=125e-6;sep=300e-6
            return sigma*np.exp(.5)/2*(np.exp(-((a+sep/2)**2+c*c)/(2*sigma*sigma))-np.exp(-((a-sep/2)**2+c*c)/(2*sigma*sigma)))
        k=np.pi/250e-6
        return np.sin(k*(a+375e-6))*np.sin(k*(c+125e-6))/k
    ga=(psi(pts+[eps,0,0])-psi(pts-[eps,0,0]))/(2*eps)
    gb=(psi(pts+[0,eps,0])-psi(pts-[0,eps,0]))/(2*eps)
    exact=np.c_[-gb,ga,np.zeros(3)];error=float(abs(exact-surf.field.spatial(pts,normals)).max());assert error<1e-7
    results[mode]=dict(free_fixture_nodes=int(free.sum()),free_supported_faces=int(free_faces.sum()),new_accessible_faces=int(z['accessible'].sum()),off_drive_power_W=off['drive_power'],passive_power_W=off['drag_power'],outside_zone_passive_power_W=far['drag_power'],fine_force_status='UNREPRESENTED_NOT_PHYSICAL_ZERO',fine_geometry_no_fake_shell=True,streamfunction_gradient_error=error,native_steps=0)
dump(O/'SURFACE_SEMANTICS_AUDIT.json',dict(status='PASS',scope='Known detached geometry fixture only; no material dynamics or tuning',cases=results,wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))
print('SURFACE_SEMANTICS_PASS',results)
