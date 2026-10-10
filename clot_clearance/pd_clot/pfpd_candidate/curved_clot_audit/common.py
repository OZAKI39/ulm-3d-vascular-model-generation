from pathlib import Path
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'verification/curved_clot_wall_application_discrepancy_audit_v1'
OLD=ROOT/'verification/curved_clot_wall_intact_mechanics_v1'
RESULTS=ROOT/'results/curved_clot_wall_application_discrepancy_audit_v1'
def plain(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,Path):return str(v)
    raise TypeError(type(v).__name__)
def dump(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(data,default=plain,indent=2,allow_nan=False)+'\n')
def read(path):return json.loads(Path(path).read_text())
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def config():return read(OUT/'AUDIT_CONFIGURATION.json')
def ring_distance(X,R):return np.hypot(np.linalg.norm(X[...,:2],axis=-1)-R,X[...,2])
def positive_energy(F,xp=np):
    E=.5*(xp.swapaxes(F,-1,-2)@F-xp.eye(3,dtype=xp.float64));e=xp.linalg.eigvalsh(E)
    return 1000*xp.sum(xp.maximum(e,0)**2,axis=-1)
def gradient_norm(F,xp=np):return xp.sqrt(xp.sum((F-xp.eye(3,dtype=xp.float64))**2,axis=(-2,-1)))
def weighted_stats(values,weights,points=None,mask=None,ids=None,save_order=None):
    values=np.asarray(values,dtype=np.float64);weights=np.broadcast_to(np.asarray(weights,dtype=np.float64),values.shape)
    assert values.ndim==1 and np.isfinite(values).all() and np.isfinite(weights).all() and np.all(weights>0)
    mask=np.ones(len(values),bool) if mask is None else np.asarray(mask,bool);ix=np.flatnonzero(mask)
    if ids is None:ids=np.arange(len(values))
    if not len(ix):raise ValueError('Empty measure')
    order=ix[np.argsort(values[ix],kind='stable')];cumulative=np.cumsum(weights[order]);total=float(cumulative[-1]);out=dict(count=len(ix),weight_sum=total,weighted_RMS=float(np.sqrt(np.average(values[ix]**2,weights=weights[ix]))),mean=float(np.average(values[ix],weights=weights[ix])))
    def entry(k):
        i=order[k];return dict(id=int(ids[i]),array_index=int(i),value=float(values[i]),position_m=None if points is None else points[i],cumulative_weight=float(cumulative[k]),cumulative_fraction=float(cumulative[k]/total))
    for name,q in [('P95',.95),('P99',.99)]:
        k=min(int(np.searchsorted(cumulative,q*total,side='left')),len(order)-1);s=entry(k);value=s['value'];ties=order[values[order]==value]
        out[name]=value;out[name+'_crossing']=s;out[name+'_bracket']=dict(previous=None if k==0 else entry(k-1),crossing=s,target_fraction=q,tied_sample_count=len(ties),tie_ids=ids[ties],location_kind='actual first-crossing representative; no percentile interpolation')
    # V1 max uses first original occurrence, independent of sorted tie order.
    i=ix[np.argmax(values[ix])];out['max']=float(values[i]);out['max_id']=int(ids[i]);out['max_position_m']=None if points is None else points[i]
    if save_order is not None:np.savez_compressed(save_order,sorted_array_indices=order,sorted_sample_ids=ids[order],sorted_values=values[order],sorted_weights=weights[order],cumulative_weights=cumulative,cumulative_fractions=cumulative/total)
    return out
def bijective_map(source_points,target_points,tolerance):
    from scipy.spatial import cKDTree
    d,ids=cKDTree(source_points).query(target_points)
    if len(source_points)!=len(target_points) or len(np.unique(ids))!=len(ids) or d.max()>tolerance:raise ValueError('Non-bijective/ambiguous native coordinate mapping')
    if len(source_points)>1 and cKDTree(source_points).query(source_points,k=2)[0][:,1].min()<=2*tolerance:raise ValueError('Coordinate identity not uniquely separable')
    return ids,float(d.max())
def owned_rows(values,global_ids,owned):
    values=np.asarray(values);ids=np.asarray(global_ids);owned=np.asarray(owned,bool)
    if len(np.unique(ids[owned]))!=owned.sum():raise ValueError('Duplicate owned global DOFs/cells: MPI double counting')
    return values[owned],ids[owned]
def recover_traces(trace,cell_vertex_ids,cell_volume,n_nodes):
    weight=np.zeros(n_nodes);F=np.zeros((n_nodes,3,3));B=np.zeros(n_nodes)
    ids=cell_vertex_ids.ravel();w=np.repeat(cell_volume,4)
    np.add.at(weight,ids,w);np.add.at(F,ids,(trace*cell_volume[:,None,None,None]).reshape(-1,3,3));np.add.at(B,ids,(positive_energy(trace)*cell_volume[:,None]).ravel())
    if np.any(weight<=0):raise ValueError('Orphan recovery node')
    F/=weight[:,None,None];B/=weight;A=positive_energy(F);dispersion=np.zeros(n_nodes)
    np.add.at(dispersion,ids,(np.sum((trace-F[cell_vertex_ids])**2,axis=(2,3))*cell_volume[:,None]).ravel());dispersion=np.sqrt(dispersion/weight)
    return F,A,B,dispersion,weight
