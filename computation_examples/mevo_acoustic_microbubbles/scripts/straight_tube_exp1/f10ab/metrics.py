"""Fixed-ID Euclidean shape decomposition, independent of the particle solver.

Column-vector convention y(t)=F y(0). No particle rematching. Principal axes are
invalid near eigenvalue degeneracy; a useful residual is still returned when a
single global affine map is inadequate. Uniform-scale and rigid residuals differ.
"""
import numpy as np
from scipy.linalg import polar
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import connected_components

RADIUS = 1.78e-6

def rms(a):
    return float(np.sqrt(np.mean(np.sum(np.asarray(a)**2,axis=-1))))

def geometry(x, radius=RADIUS):
    x=np.asarray(x,float)
    if x.ndim!=2 or x.shape[1]!=3 or len(x)<4 or not np.isfinite(x).all():
        raise ValueError('Finite N by 3 fixed-ID positions required')
    c=x.mean(0);y=x-c;result={'center_3d_m':c,'center_xy_m':c[:2]}
    for dim,key in [(3,'3d'),(2,'xy')]:
        yy=y[:,:dim];q=yy.T@yy/len(x);lam,vec=np.linalg.eigh(q);lam=lam[::-1];vec=vec[:,::-1]
        top=max(float(lam[0]),1e-40);valid=np.array([min([abs(lam[k]-lam[j])/top for j in range(dim) if j!=k])>.03 for k in range(dim)])
        result.update({f'covariance_{key}_m2':q,f'eigenvalues_{key}_m2':lam,
            f'principal_axes_{key}':np.where(valid[None,:],vec,np.nan),f'orientation_valid_{key}':valid,
            f'principal_rms_lengths_{key}_m':np.sqrt(np.maximum(lam,0)),
            f'Rg_{key}_m':float(np.sqrt(max(0,np.trace(q)))),
            f'aspect_{key}':float(np.sqrt(lam[0]/lam[-1])) if lam[-1]>top*1e-12 else np.nan})
    distances=pdist(x);d=squareform(distances);np.fill_diagonal(d,np.inf)
    edges=d<=2.1*radius;nc,labels=connected_components(edges,directed=False)
    result.update(pair_distance_statistics_m=np.quantile(distances,[0,.1,.5,.9,1]),nearest_median_m=float(np.median(d.min(1))),
       minimum_center_distance_m=float(distances.min()),proximity_components=int(nc),largest_component=int(np.bincount(labels).max()),
       proximity_edges=int(edges.sum()//2),hard_contact_pairs=int(np.count_nonzero(distances<=2*radius+1e-9)),
       component_labels=labels,extent_3d_m=np.ptp(x,axis=0))
    return result

def deformation(x0,x):
    x0=np.asarray(x0,float);x=np.asarray(x,float)
    if x.shape!=x0.shape or not np.isfinite(x).all() or not np.isfinite(x0).all():raise ValueError('Same finite fixed-ID arrays required')
    a=x0-x0.mean(0);b=x-x.mean(0);out={}
    for dim,key in [(3,'3d'),(2,'xy')]:
        aa=a[:,:dim];bb=b[:,:dim];q,s,vt=np.linalg.svd(aa.T@bb)
        det=np.ones(dim);det[-1]=np.linalg.det(vt.T@q.T);rot=(vt.T*det)@q.T
        sv=np.linalg.svd(aa,compute_uv=False);rank=int(np.sum(sv>max(sv[0],1e-40)*1e-8));cond=sv[0]/sv[-1] if sv[-1]>0 else np.inf
        rig=aa@rot.T;scale=float(np.sum(rig*bb)/np.sum(aa*aa)) if np.sum(aa*aa)>0 else np.nan
        aligned=bb@rot
        out.update({f'rigid_rotation_{key}':rot if rank>=dim-1 else np.full((dim,dim),np.nan),
            f'kabsch_rmsd_{key}_m':rms(bb-rig),f'similarity_scale_{key}':scale,
            f'similarity_rmsd_{key}_m':rms(bb-scale*rig),f'translation_removed_rms_{key}_m':rms(bb-aa),
            f'unregistered_rms_{key}_m':rms(x[:,:dim]-x0[:,:dim]),
            f'pair_distance_rms_change_{key}_m':float(np.sqrt(np.mean((pdist(bb)-pdist(aa))**2))),
            f'initial_condition_{key}':float(cond),f'initial_rank_{key}':rank})
        fit=np.linalg.lstsq(aa,bb,rcond=1e-8)[0].T
        residual=rms(bb-aa@fit.T);rel=residual/max(rms(bb),1e-40)
        valid=rank==dim and cond<100 and np.linalg.det(fit)>0 and rel<=.10
        out[f'affine_valid_{key}']=bool(valid);out[f'nonaffine_residual_{key}_m']=residual
        out[f'nonaffine_relative_{key}']=rel
        out[f'affine_deformation_{key}']=fit if valid else np.full((dim,dim),np.nan)
        if valid:
            rr,stretch=polar(fit);eig=np.linalg.eigvalsh(stretch);iso=float(np.mean(eig));dev=stretch-iso*np.eye(dim)
            shear=dev-np.diag(np.diag(dev))
            out.update({f'polar_rotation_{key}':rr,f'stretch_tensor_{key}':stretch,f'principal_stretches_{key}':eig[::-1],
                f'isotropic_stretch_{key}':iso,f'deviatoric_stretch_{key}':float(np.linalg.norm(dev)),
                f'shear_norm_{key}':float(np.linalg.norm(shear)),f'directional_xy_stretch_{key}':float(stretch[0,0]-stretch[1,1])})
        else:
            for k in ['polar_rotation','stretch_tensor']:out[f'{k}_{key}']=np.full((dim,dim),np.nan)
            out[f'principal_stretches_{key}']=np.full(dim,np.nan)
            for k in ['isotropic_stretch','deviatoric_stretch','shear_norm','directional_xy_stretch']:out[f'{k}_{key}']=np.nan
    out['aligned_positions_m']=b@out['rigid_rotation_3d']
    return out

def measure(x0,x):
    return {**geometry(x),**deformation(x0,x)}

def series(times,positions,stride=1):
    ix=np.unique(np.r_[np.arange(0,len(times),stride),len(times)-1]).astype(int)
    rows=[measure(positions[0],positions[k]) for k in ix]
    return {'time_s':np.asarray(times)[ix], 'sample_index':ix,**{key:np.asarray([row[key] for row in rows]) for key in rows[0]}}
