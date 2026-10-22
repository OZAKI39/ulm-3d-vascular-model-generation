"""Conservative Eulerian 2-D trends; no particle identities or physical fusion claims."""
from collections import deque
import numpy as np
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN

def distribution(q):
 q=np.asarray(q,float).reshape(-1,2);n=len(q)
 if n<2:return dict(count=n,center=None,rg_px=None,nn_px=None,knn_density=None,shape_tensor=None,aspect=None,angle_deg=None)
 c=q.mean(0);d=q-c;cov=d.T@d/n;e,v=np.linalg.eigh(cov);nn=cKDTree(q).query(q,k=min(7,n))[0]
 return dict(count=n,center=c.tolist(),rg_px=float(np.sqrt(np.trace(cov))),nn_px=float(np.median(nn[:,1])),knn_density=float(np.median((nn.shape[1]-1)/(np.pi*np.maximum(nn[:,-1],1)**2))),shape_tensor=cov.tolist(),aspect=float(np.sqrt(e[-1]/max(e[0],1e-6))),angle_deg=float(np.degrees(np.arctan2(v[1,-1],v[0,-1]))))
class Aggregation:
 def __init__(self,window_s=.8,min_span=.5):self.history={};self.window=window_s;self.min_span=min_span;self.prev_regions=[]
 def update(self,pts,q,xy,velocity,shape,focus=1.,illum=0.,boundary_fraction=0.,cells=True,support_mask=None):
  q=np.asarray(q,float).reshape(-1,2);xy=np.asarray(xy,float).reshape(-1,2);velocity=np.asarray(velocity,float).reshape(-1,2)
  global_d=distribution(q);regions=[];keys=[('all',0,0,shape[1],shape[0])]
  if cells:
   keys += [(f'{x}:{y}',x,y,min(300,shape[1]-x),min(300,shape[0]-y)) for y in range(0,shape[0],300) for x in range(0,shape[1],300)]
  for key,x,y,w,h in keys:
   sel=(q[:,0]>=x)&(q[:,0]<x+w)&(q[:,1]>=y)&(q[:,1]<y+h);p=q[sel];d=distribution(p)
   area=int(np.sum(support_mask[y:y+h,x:x+w])) if support_mask is not None else w*h
   d.update(key=key,bounds=[x,y,w,h],observed_area_px=area,density_per_10k_px=len(p)/max(area,1)*1e4,focus=focus,illum=illum,boundary_fraction=boundary_fraction,pts=pts)
   choose=(xy[:,0]>=x)&(xy[:,0]<x+w)&(xy[:,1]>=y)&(xy[:,1]<y+h);pos=xy[choose];vel=velocity[choose];support=len(pos)
   d.update(motion_support=support,radial_px_s=None,divergence_s=None,mean_velocity=None,coherence=None)
   if support>=8 and d['center'] is not None:
    rad=pos-np.asarray(d['center']);den=np.linalg.norm(rad,axis=1);vm=vel.mean(0);relative=vel-vm
    speed=np.linalg.norm(vel,axis=1);use=den>3
    X=np.c_[np.ones(support),pos-pos.mean(0)];coef=np.linalg.lstsq(X,vel,rcond=None)[0]
    d.update(mean_velocity=vm.tolist(),coherence=float(np.linalg.norm(vm)/max(np.mean(speed),1e-9)),radial_px_s=float(np.median(np.sum(relative[use]*rad[use],axis=1)/den[use])) if use.any() else None,divergence_s=float(coef[1,0]+coef[2,1]))
   hist=self.history.setdefault(key,deque(maxlen=90));hist.append(d.copy())
   while hist and pts-hist[0]['pts']>self.window:hist.popleft()
   reasons=[];slope_nn=slope_rg=None;status='NOT_OBSERVABLE';change='INSUFFICIENT_HISTORY';converge=False
   if len(hist)>=8 and pts-hist[0]['pts']>=self.min_span and d['count']>=8:
    counts=np.array([a['count'] for a in hist]);f=np.array([a['focus'] for a in hist]);il=np.array([a['illum'] for a in hist]);times=np.array([a['pts'] for a in hist]);good=[a for a in hist if a['rg_px'] and a['nn_px']]
    if np.ptp(counts)/max(np.mean(counts),1)>.15:reasons.append('VISIBLE_COUNT_CHANGE')
    if np.ptp(f)/max(np.mean(f),1e-9)>.25:reasons.append('FOCUS_CHANGE')
    if np.ptp(il)>12:reasons.append('ILLUMINATION_CHANGE')
    if boundary_fraction>.15:reasons.append('BOUNDARY_FLUX_RISK')
    if support<8:reasons.append('LOW_MOTION_SUPPORT')
    if len(good)!=len(hist):reasons.append('LOW_CANDIDATE_SUPPORT')
    else:
     slope_nn=float(np.polyfit(times-times[0],np.log(np.maximum([a['nn_px'] for a in hist],1e-6)),1)[0]);slope_rg=float(np.polyfit(times-times[0],np.log(np.maximum([a['rg_px'] for a in hist],1e-6)),1)[0])
     change='CONTRACTING_DISTRIBUTION' if slope_nn<-.025 and slope_rg<-.025 else 'EXPANDING_DISTRIBUTION' if slope_rg>.025 else 'NO_CONSISTENT_CONTRACTION'
     validmotion=[a for a in hist if a['divergence_s'] is not None]
     converge=len(validmotion)>=8 and np.median([a['divergence_s'] for a in validmotion])<-.04 and np.median([a['radial_px_s'] or 0 for a in validmotion])<-1.
     if not reasons:status='VISIBLE_CONVERGENCE' if change=='CONTRACTING_DISTRIBUTION' and converge else 'VISIBLE_DISPERSION' if change=='EXPANDING_DISTRIBUTION' and d['divergence_s'] and d['divergence_s']>.04 else 'NO_CONSISTENT_CONVERGENCE'
   else:reasons.append('WARMUP_OR_LOW_SUPPORT')
   d.update(visible_distribution_change=change,local_convergence_evidence=bool(converge),aggregation_trend_status=status,observation_quality='OBSERVABLE_2D_ONLY' if not reasons else 'NOT_OBSERVABLE',reasons=reasons,nn_log_slope=slope_nn,rg_log_slope=slope_rg)
   regions.append(d)
  clusters=[]
  if len(q)>=4:
   labels=DBSCAN(eps=55,min_samples=4,n_jobs=1).fit_predict(q)
   for label in set(labels)-{-1}:
    a=q[labels==label];d=distribution(a);d['extent']=[*a.min(0),*a.max(0)];d['physical_id']=None
    d['spatial_persistence']=bool(any(np.linalg.norm(np.asarray(d['center'])-p)<55 for p in self.prev_regions));clusters.append(d)
  self.prev_regions=[np.asarray(a['center']) for a in clusters]
  return dict(global_distribution=global_d,regions=regions,clusters=clusters,status='VISIBLE_CONVERGENCE' if any(r['aggregation_trend_status']=='VISIBLE_CONVERGENCE' for r in regions) else regions[0]['aggregation_trend_status'],scientific_accuracy=False,scope='2D_VISIBLE_DISTRIBUTION_ONLY')
