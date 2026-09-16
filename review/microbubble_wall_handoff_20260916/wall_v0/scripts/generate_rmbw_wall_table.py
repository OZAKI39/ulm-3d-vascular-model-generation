"""OFFLINE ONLY: unmodified frozen RMBW binding; outputs never imported at runtime.

Writes first raw generator outputs as NumPy arrays. HDF5 packaging/validation is
performed separately so the frozen Python 3.8 reference environment stays intact.
"""
from pathlib import Path
import json,sys,hashlib,time
import numpy as np
R=Path(__file__).resolve().parents[1];P=json.loads((R/'provenance/TASK_PATHS.json').read_text())
audit=Path(P['reference_audit']);sys.path.insert(0,str(audit/'src'))
from rmbw_reference_adapter import get_wall_mobility
C=json.loads((R/'provenance/reference_audit/WALL_REFERENCE_CONVENTION.json').read_text())
seed=2026091601; n_holdout=5500
# Dense log grid plus original nodes. Near-neighbors resolve existing branch
# boundaries without smoothing or modifying reference coefficients.
extra=[.01,.1,.2,5.,9.018296,10.,20.]+C['gaps_h_over_a']
for e in [.01,.1,9.018296]:extra += [e*(1-1e-10),e*(1+1e-10)]
grid=np.unique(np.r_[np.geomspace(.001,20.,8192),extra]);grid.sort()
rng=np.random.default_rng(seed);hold=np.exp(rng.uniform(np.log(.001),np.log(20),n_holdout))
assert not np.any(np.isin(hold,grid))
plan={'status':'FROZEN_BEFORE_GENERATION','seed':seed,'grid_base_points':8192,'grid_points':len(grid),
    'grid':grid.tolist(),'holdout_epsilon':hold.tolist(),'holdout_count':n_holdout,
    'radii_m':C['radii_m'],'mu_Pa_s':.001,'action_vectors_per_epsilon':128,
    'action_error_gate':1e-3,'qualified_max_epsilon':5.,'reciprocity_gate':1e-10,
    'positive_scaled_eigenvalue_gate':0.,'interpolation':'piecewise linear in log(epsilon)',
    'reference_commit':'8d41e464d7de9b6514a85a741dd227f1219e3a49',
    'reference_adapter_sha256':hashlib.sha256((audit/'src/rmbw_reference_adapter.py').read_bytes()).hexdigest(),
    'convention_sha256':hashlib.sha256((R/'provenance/reference_audit/WALL_REFERENCE_CONVENTION.json').read_bytes()).hexdigest(),
    'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'no_reference_modification':True,'runtime_dependency':'NONE'}
planpath=R/'contracts/TABLE_GENERATION_PLAN.json'
assert not planpath.exists(),'Preserve first generation plan/output; no silent overwrite'
planpath.write_text(json.dumps(plan,indent=2)+'\n')
def evaluate(epsilon,a):
    M,meta=get_wall_mobility(a,a*epsilon,.001,implementation='lubrication')
    scale=np.sqrt([6*np.pi*.001*a]*3+[8*np.pi*.001*a**3]*3)
    B=M*np.outer(scale,scale);T=np.linalg.solve(B,np.eye(6))
    return M,T*np.outer(scale,scale),B,T
start=time.time()
for label,ee,a in [('table',grid,C['radii_m']['d50'])]+[(name,hold,a) for name,a in C['radii_m'].items()]:
    dest=R/'raw'/('RMBW_'+label+'_FIRST.npz');assert not dest.exists()
    matrices=[[],[],[],[]]
    for k,e in enumerate(ee):
        result=evaluate(float(e),float(a))
        for out,value in zip(matrices,result):out.append(value)
        if (k+1)%1000==0:print(label,k+1,'wall_s',time.time()-start,flush=True)
    np.savez_compressed(dest,epsilon=ee,radius_m=a,mu_Pa_s=.001,
        M_SI=matrices[0],R_total_SI=matrices[1],M_scaled=matrices[2],R_total_scaled=matrices[3])
    print(label,'DONE',len(ee),flush=True)
print('OFFLINE_REFERENCE_COMPLETE',time.time()-start,flush=True)
