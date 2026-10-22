"""Common adapters and frozen scoring; all geometric estimates reuse F10U."""
import time,copy
import numpy as np
from wcommon import old,NAMES
from motion_models import DenseOfficial,CoTrackerOfficial,optimal_transport,point_array
LIMITS=np.array([5,5,.2,.2,.1,.1])
def construct(name):
    return CoTrackerOfficial() if name=='CoTracker3' else (DenseOfficial(name) if name in ['LRAFT_ST','UnLiteFlowNet'] else None)
def predict(name,model,images,objects,times,config):
    start=time.perf_counter()
    if model is not None:r=model.pair(images,objects,times,**config)
    elif name.startswith('POT_'):r=optimal_transport(images,objects,times,method=name[4:],**config)
    elif name=='PIV':r=old.piv(images[-3],images[-1],times[-1]-times[-3],mode='candidate',objects_a=objects[-3],objects_b=objects[-1])
    elif name=='CPD':r=old.cpd(point_array(objects[-3]),point_array(objects[-1]),times[-1]-times[-3],images[-1].shape)
    else:raise ValueError(name)
    e=r['estimate'];r.setdefault('quality_score',1/(1+e['residual']) if e.get('residual') is not None else 0.)
    r.setdefault('timing',{})['adapter_total_ms']=(time.perf_counter()-start)*1000
    r.setdefault('vectors',[])
    return r
def scored(r,truth,label):
    e=r['estimate'];raw=e.get('raw_estimate');finite=raw is not None and np.isfinite(raw).all()
    err=np.abs(np.asarray(raw)-truth) if finite and label==0 else None
    valid=all(e['valid'].values());bad=label!=0 or not finite or bool(np.any(err>LIMITS))
    return dict(raw=raw,absolute_error=err,base_valid=valid,bad_if_accepted=bad,joint_normalized_error=float(np.max(err/LIMITS)) if err is not None else None)
def apply_gate(estimate,quality,threshold):
    e=copy.deepcopy(estimate)
    if threshold is None or quality<threshold:
        e['values']={k:None for k in NAMES};e['valid']={k:False for k in NAMES};e['reasons']={k:'VALIDATION_RISK_GATE_NOT_ADMITTED' for k in NAMES};e['observability_status']='NOT_OBSERVABLE'
    return e
