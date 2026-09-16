"""Bounded static experiments; upstream hydrodynamic operators only."""
from pathlib import Path
import os,sys,time,json,signal,traceback,resource
import numpy as np,h5py,psutil
ROOT=Path(os.environ['AUDIT_ROOT']);OUT=Path(os.environ['AUDIT_RESULT']);sys.path.insert(0,str(ROOT/'src'))
from pecnut_fixed_wall_adapter import FixedWall,lattice,phases
C=json.loads((ROOT/'contracts/FIXED_MULTIBLOB_WALL_FEASIBILITY_CONTRACT.json').read_text())
beta=float(sys.argv[1]);mode=sys.argv[2] if len(sys.argv)>2 else 'planar';start=time.monotonic()
OUT.mkdir(parents=True,exist_ok=True)
fn=OUT/f'RAW_{mode}_beta{beta:g}.h5'
if fn.exists():raise RuntimeError('NO_AUTOMATIC_OVERWRITE_OR_RETRY')
def timeout(signum,frame):raise TimeoutError('FROZEN_PER_QUERY_TIME_LIMIT')
signal.signal(signal.SIGALRM,timeout)
def guard(N):
 estimated=48*(11*N)**2;avail=psutil.virtual_memory().available
 if estimated>min(20*2**30,.6*avail):raise MemoryError(f'PREFLIGHT N={N} estimated={estimated} available={avail}')
 if time.monotonic()-start>3600:raise TimeoutError('FROZEN_TOTAL_LIMIT')
def write(g,d):
 for k,v in d.items():
  if isinstance(v,np.ndarray):g.create_dataset(k,data=v,compression='gzip' if v.size>100 else None)
  elif v is None:g.attrs[k]='NONE'
  else:g.attrs[k]=v
 g.attrs['rss_bytes']=psutil.Process().memory_info().rss;g.attrs['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024

def geometry(f,layout,layers,extent,spacing,category='PLANAR',count=None):
 name=f'{category}_{layout}_L{layers}_R{extent:g}_S{spacing:g}' + (f'_N{count}' if count else '')
 x,a,ids,meta=lattice(layout,beta,layers,extent,spacing,count)
 g=f.create_group(name);write(g,dict(coords=x,radii=a,layer=ids,beta=beta,layout=layout,layers=layers,extent=extent,spacing_ratio=spacing,category=category,N_wall=len(a),lattice_metadata=json.dumps(meta)))
 try:
  guard(len(a));signal.alarm(180);t=time.perf_counter();wall=FixedWall(x,a);signal.alarm(0)
  g.attrs['wall_total_seconds']=time.perf_counter()-t;g.attrs['wall_assembly_seconds']=wall.assembly_seconds;g.attrs['wall_factor_seconds']=wall.factor_seconds;g.attrs['operator_bytes']=wall.peak_operator_bytes
  from scipy.spatial import cKDTree
  g.attrs['minimum_wall_wall_gap']=float((cKDTree(x).query(x,k=2)[0][:,1]-2*beta).min())
  g.attrs['status']='ASSEMBLED';f.flush();print('ASSEMBLED',name,len(a),g.attrs['wall_total_seconds'],flush=True)
  return g,wall
 except Exception as e:
  signal.alarm(0);g.attrs['status']='NOT_COMPLETED';g.attrs['reason']=str(e);f.flush();print('BLOCKED',name,repr(e),flush=True);return g,None

def query(f,g,wall,name,x,a,attrs):
 q=g.create_group(name);write(q,dict(targets=np.atleast_2d(x),target_radii=np.atleast_1d(a),**attrs))
 try:
  guard(wall.N);signal.alarm(180);d=wall.query(x,a,store_forces=False);signal.alarm(0);write(q,d);q.attrs['status']='COMPLETED'
 except Exception as e:signal.alarm(0);q.attrs['status']='FAILED';q.attrs['reason']=repr(e);print('QUERY_FAIL',name,repr(e),flush=True)
 f.flush()

with h5py.File(fn,'x') as f:
 f.attrs['contract_sha256']=__import__('hashlib').sha256((ROOT/'contracts/FIXED_MULTIBLOB_WALL_FEASIBILITY_CONTRACT.json').read_bytes()).hexdigest();f.attrs['source']='UNMODIFIED_PECNUT_OPERATORS';f.attrs['mode']=mode;f.attrs['beta']=beta
 if mode=='planar':
  for layout in C['layouts']:
   for layers in C['layers']:
    g,w=geometry(f,layout,layers,4,1.01)
    if w:
     for ei,eps in enumerate(C['epsilon']):
      for pi,(label,xy) in enumerate(phases(layout,beta)):
       query(f,g,w,f'Q{ei:02}_{pi:02}',[*xy,1+eps],[1.],dict(epsilon=eps,phase=label,phase_index=pi))
      print('EPS_DONE',g.name,eps,flush=True)
    del w
 elif mode=='patch':
  # Both layouts, all layer counts; special positions. Shared r=4 queries are already in planar data.
  for extent in [6,8,12]:
   for layout in C['layouts']:
    for layers in C['layers']:
     g,w=geometry(f,layout,layers,extent,1.01,'PATCH')
     if w:
      for ei,eps in enumerate(C['epsilon']):
       for pi,(label,xy) in enumerate(phases(layout,beta)[:3]):query(f,g,w,f'Q{ei:02}_{pi:02}',[*xy,1+eps],[1.],dict(epsilon=eps,phase=label,phase_index=pi))
     del w
 elif mode=='diagnostic':
  for layout in C['layouts']:
   for layers in C['layers']:
    for spacing in C['spacing_ratio_sweep']:
     g,w=geometry(f,layout,layers,4,spacing,'LEAKAGE_SPACING')
     if w:
      lower=float(w.x[:,2].min()-beta-2);targets=np.array([[0,0,2],[0,0,lower]])
      query(f,g,w,'LEAKAGE',targets,[1.,1.],dict(kind='TWO_PROBE_TRANSMISSION'))
      free=FixedWall(np.empty((0,3)),np.empty(0));query(f,g,free,'LEAKAGE_FREE',targets,[1.,1.],dict(kind='NO_WALL_NORMALIZER'))
      for ei,eps in enumerate([.01,.05,.2]):
       for pi,(label,xy) in enumerate(phases(layout,beta,spacing)[:3]):query(f,g,w,f'Q{ei:02}_{pi:02}',[*xy,1+eps],[1.],dict(epsilon=eps,phase=label,phase_index=pi))
     del w
 elif mode=='performance':
  for N in C['performance_wall_counts']+C['performance_optional_counts']:
   extent=float(np.ceil(np.sqrt(N)*2*beta));g,w=geometry(f,'HEX',1,extent,1.01,'PERFORMANCE',N)
   if w:
    rng=np.random.default_rng(C['seeds']['repeated_queries']);xy=rng.uniform(-.5,.5,(100 if N==100 else 3,2))
    for i,p in enumerate(xy):query(f,g,w,f'Q{i:03}',[*p,1.1],[1.],dict(epsilon=.1,query_index=i))
   del w
 f.attrs['terminal_state']='COMPLETED_WITH_RECORDED_GATES';f.attrs['wall_seconds']=time.monotonic()-start
print('TERMINAL',fn,flush=True)
