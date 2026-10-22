"""Independent target labels, fixed IDs, and separate scalar/structural recovery."""
import numpy as np
from scipy.spatial.distance import pdist
from scripts.straight_tube_exp1.f10ab.metrics import measure,geometry
from .common import R,read,CONFIG
def target_metrics(x0,x,numerical=True):
    a=measure(x0,x0);b=measure(x0,x);cfg=read(CONFIG)
    gain=b['aspect_3d']/a['aspect_3d']-1;com=float(np.linalg.norm(x.mean(0)-x0.mean(0)))
    delta=np.asarray(b['principal_rms_lengths_3d_m'])-a['principal_rms_lengths_3d_m']
    out=dict(aspect_gain=gain,Rg_gain=b['Rg_3d_m']/a['Rg_3d_m']-1,COM_m=com,axis_delta_m=delta,axis_relative_change=np.asarray(b['principal_rms_lengths_3d_m'])/a['principal_rms_lengths_3d_m']-1,NUMERICAL_PASS=bool(numerical),COM_TARGET_PASS=com<=cfg['T1']['COM_limit_m'],endpoint_connected=b['proximity_components']==1,components=b['proximity_components'],kabsch_m=b['kabsch_rmsd_3d_m'],orientation_valid=b['orientation_valid_3d'],raw_extension_above_floor=bool(np.max(delta)>cfg['T3']['axis_extension_floor_m']),axis_class='ANISOTROPIC_CONTRACTION' if np.all(delta<0) else 'EXTENSION_CANDIDATE_REQUIRES_CONVERGENCE')
    out['SHAPE_TARGET_PASS']=bool(abs(gain-cfg['T1']['aspect_gain'])<=cfg['T1']['tolerance'])
    out['T1_DISCOVERY_PASS']=out['SHAPE_TARGET_PASS'] and out['COM_TARGET_PASS'] and numerical
    out['T2_SHAPE_PASS']=bool(abs(gain-cfg['T2']['aspect_gain'])<=cfg['T2']['tolerance'])
    out['T2_DISCOVERY_PASS']=out['T2_SHAPE_PASS'] and out['COM_TARGET_PASS'] and out['endpoint_connected'] and numerical
    out['T1_loss']=float(((gain-cfg['T1']['aspect_gain'])/cfg['T1']['tolerance'])**2)
    out['T2_loss']=float(((gain-cfg['T2']['aspect_gain'])/cfg['T2']['tolerance'])**2)
    return out
def recovery(x0,xa,xb):
    a=measure(x0,xa);b=measure(x0,xb);g=geometry(x0)
    qa=g['covariance_3d_m2'];scale=np.linalg.norm(qa)
    edges0=pdist(x0)<=2.1*R
    def errors(y,m):
        edges=pdist(y)<=2.1*R
        jac=float(np.sum(edges & edges0)/max(np.sum(edges | edges0),1))
        # Covariance in initial frame after proper Kabsch alignment.
        aligned=m['aligned_positions_m'];q=aligned.T@aligned/len(y)
        return dict(Rg=abs(m['Rg_3d_m']/g['Rg_3d_m']-1),axes=float(np.max(abs(m['principal_rms_lengths_3d_m']/g['principal_rms_lengths_3d_m']-1))),covariance=float(np.linalg.norm(q-qa)/scale),Kabsch_m=m['kabsch_rmsd_3d_m'],pair_distance_m=m['pair_distance_rms_change_3d_m'],neighbor_Jaccard=jac,neighbor_error=1-jac)
    e1=errors(xa,a);e2=errors(xb,b)
    gain={k:(1-e2[k]/v if v>1e-30 else None) for k,v in e1.items() if k!='neighbor_Jaccard'}
    cfg=read(CONFIG)['T4'];tol=dict(Rg=cfg['Rg_relative_tolerance'],axes=cfg['axis_relative_tolerance'],covariance=cfg['covariance_relative_tolerance'],Kabsch_m=cfg['fixed_ID_Kabsch_tolerance_m'],neighbor_error=1-cfg['neighbor_Jaccard_min'])
    passes={k:bool(e2[k]<=v and gain[k] is not None and gain[k]>=cfg['minimum_recovery_fraction']) for k,v in tol.items()}
    return dict(first_errors=e1,final_errors=e2,recovery_fraction=gain,raw_recovery_pass=passes,scope='Requires dt/cache and independent-state admission; Rg/axes do not substitute for fixed-ID or topology recovery')
