"""Independent, readonly geometric diagnostics, without mesher substitutes."""
import importlib.util
import sys
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from .validation import require
from .provenance import sha256


def unique_edges(tetra):
    t=np.asarray(tetra,int)
    return np.unique(np.sort(np.concatenate([t[:,[i,j]] for i,j in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]]),axis=1),axis=0)


def volume_validity(points,tetra):
    p=np.asarray(points,float);t=np.asarray(tetra,int)
    require(p.ndim==2 and p.shape[1]==3 and t.ndim==2 and t.shape[1]==4 and len(t)>0,'Tetrahedral volume required')
    require(np.isfinite(p).all() and t.min()>=0 and t.max()<len(p),'Invalid vertex coordinates or connectivity')
    x=p[t];volume=np.linalg.det(x[:,1:]-x[:,0,None,:])/6
    require(np.isfinite(volume).all() and np.all(volume>0),'Inverted, degenerate or non-finite tetrahedra')
    faces=np.sort(np.concatenate([t[:,ids] for ids in [[0,1,2],[0,1,3],[0,2,3],[1,2,3]]]),axis=1)
    unique,inverse,counts=np.unique(faces,axis=0,return_inverse=True,return_counts=True)
    require(np.all(counts<=2),'Non-manifold volume facets')
    order=np.argsort(inverse);starts=np.r_[0,np.cumsum(counts)[:-1]];interior=np.flatnonzero(counts==2)
    first=order[starts[interior]]%len(t);second=order[starts[interior]+1]%len(t)
    graph=coo_matrix((np.ones(2*len(first)),(np.r_[first,second],np.r_[second,first])),shape=(len(t),len(t)))
    ncomponents=connected_components(graph,directed=False,return_labels=False);require(ncomponents==1,'Disconnected tetra volume')
    boundary=unique[counts==1];edges=np.sort(np.concatenate([boundary[:,[0,1]],boundary[:,[0,2]],boundary[:,[1,2]]]),axis=1)
    _,edgecounts=np.unique(edges,axis=0,return_counts=True);require(np.all(edgecounts==2),'Exterior is not closed')
    return {'tetra':len(t),'vertices':len(p),'edges':len(unique_edges(t)),'connected_volumes':ncomponents,
            'inverted':0,'degenerate':0,'nonfinite':0,'exterior_closed':True,'volume_m3':float(volume.sum())}


def quality_summary(q):
    q=np.asarray(q,float);require(q.ndim==1 and len(q)>0 and np.isfinite(q).all(),'Invalid independent quality values')
    vals=np.quantile(q,[0,.01,.05,.5,.95])
    return {**dict(zip(('q_min','P1','P5','median','P95'),map(float,vals))),'N_low':int(np.count_nonzero(q<.1)),'low_fraction':float(np.mean(q<.1))}


def quality_policy(baseline,candidate):
    limits={name:factor*baseline[name] for name,factor in [('q_min',.80),('P1',.90),('P5',.90),('median',.95)]}
    limits['low_fraction']=max(5*baseline['low_fraction'],.0005)
    passed=all(candidate[name]>=limit for name,limit in limits.items() if name!='low_fraction') and candidate['low_fraction']<=limits['low_fraction']
    return {'status':'PASS' if passed else 'MESH_QUALITY_REVIEW_REQUIRED','limits':limits}


def readonly_pressure_support(old_fem_root,data):
    path=Path(old_fem_root)/'src/fem3d/vascular_rank_diagnostics.py'
    before=sha256(path);saved=sys.dont_write_bytecode;sys.dont_write_bytecode=True
    try:
        spec=importlib.util.spec_from_file_location('sv0_readonly_formal_pressure_audit',path)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        result=module.constrained_pressure_support(data)
    finally:sys.dont_write_bytecode=saved
    require(sha256(path)==before,'Old audit source changed')
    return {'source_path':str(path),'source_sha256':before,'pde_solve_called':False,'copied_solver_code':False,
            'unsupported_pressure_dofs':len(result['uncoupled_pressure_vertices']),'details':result}
