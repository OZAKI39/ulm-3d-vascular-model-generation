"""Geometry-only fixtures and historical replay; no main dynamics calibration."""
from common import *
from current_surface import CurrentSurface
import time,resource
start=time.perf_counter();cpu=time.process_time()
ref=load_npz(O/'provenance/reference.npz');X=ref['X'];n=int(ref['clot_count']);b=ref['bonds']
surf=CurrentSurface(ref,'counter_rotating');d=ref['initial_cut'].astype(float);v=np.zeros_like(X)
results={};z=surf.evaluate(X,v,d,.002)
assert len(np.unique(np.sort(surf.tri,axis=1),axis=0))==len(surf.tri)
assert not z['accessible'].any() and not z['open_gap_tetra'].any()
assert not z['force'][n:].any()
results['initial']=surf.metrics(z,X)
# Known opening crack: all original crossing bonds cut in fixture only, no native call.
upper=X[:,2]>0;cut=upper[b[:,0]]!=upper[b[:,1]];df=d.copy();df[cut]=1
xp=X.copy();xp[upper,2]+=20e-6
opened=surf.evaluate(xp,v,df,.002)
assert opened['accessible'].sum()>0 and opened['fluid_tetra'].sum()>0
assert not np.any(opened['original']&opened['accessible'])
results['known_open_fixture']=surf.metrics(opened,xp)
closed=surf.evaluate(X,v,df,.002)
assert not closed['accessible'].any()
results['closed_fixture_new_loaded']=int(closed['accessible'].sum())
# Isolated internal pocket: no gap-cell route to the exterior => no new wet face.
inside=(abs(X[:,0])<100e-6)&(abs(X[:,1])<80e-6)&(abs(X[:,2]-150e-6)<80e-6)
di=d.copy();di[inside[b[:,0]]!=inside[b[:,1]]]=1
xi=X.copy();xi[inside]*=.8
zi=surf.geometry(xi,di)
assert not zi['accessible'].any()
results['sealed_pocket_new_loaded']=int(zi['accessible'].sum())
# Exact mapping virtual-work, resultant, torque and rigid-translation covariance.
rng=np.random.default_rng(631);vv=rng.normal(size=X.shape);zz=surf.evaluate(xp,vv,df,.002)
ferr=np.max(abs(zz['force'].sum(0)-zz['face_force'].sum(0)))
merr=np.max(abs(np.cross(xp,zz['force']).sum(0)-np.cross(zz['centroid'],zz['face_force']).sum(0)))
perr=abs(float(np.sum(zz['force']*vv))-float(np.sum(zz['face_force']*zz['face_velocity'])))
assert ferr<1e-16 and merr<1e-19 and perr<1e-16 and zz['drag_power']<=0
results['mapping']=dict(force_error_N=ferr,moment_error_Nm=merr,virtual_power_error_W=perr,anchor_force_N=float(abs(zz['force'][n:]).max()))
# Historical accepted states: read-only surface reconstruction, no mechanics rerun.
replays=[]
for mode in ['uniform','counter_rotating']:
 for name in ['first_failure','peak_load','final']:
  p=OLD/'cases'/mode/'keys'/f'{name}.npz';s=load_npz(p)
  zz=surf.evaluate(s['x'],s['velocity'],s['d'],float(s['time']))
  replays.append(dict(mode=mode,key=name,**surf.metrics(zz,s['x'])))
results['historical_replays']=replays
dump(O/'GEOMETRY_CHECKS.json',dict(status='PASS',checks=results,scope='Geometry fixtures and immutable old states only; zero new native dynamics.',wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))
dump(O/'SURFACE_MAPPING_AUDIT.json',results['mapping'])
print(json.dumps(results,default=cv),flush=True)
