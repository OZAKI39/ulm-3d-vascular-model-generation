"""SI action-state adapter; no images, inferred material IDs or time filling."""
import numpy as np,h5py,json
from scipy.spatial.distance import pdist,squareform
from scipy.sparse.csgraph import connected_components
from scripts.straight_tube_exp1.f10o.analyze import metrics
from .common import R,write,sha

def validate_trajectory(t,x,ids):
 t=np.asarray(t);x=np.asarray(x);ids=np.asarray(ids)
 if x.shape!=(len(t),32,3) or not len(t) or not np.isfinite(x).all() or not np.isfinite(t).all():raise ValueError('Finite 32-particle three-dimensional trajectory required')
 if t[0]!=0 or np.any(np.diff(t)<=0) or not np.array_equal(ids,np.arange(32)):raise ValueError('No duplicate times or changed fixed IDs')
 return t,x

def series(t,x):
 t,x=validate_trajectory(t,x,np.arange(32));m=metrics(t,x)
 centered=m['centered_positions_m'];qxy=m['tensor_m2'][:,:2,:2];e,vec=np.linalg.eigh(qxy)
 eigen_valid=(e[:,-1]-e[:,0])>1e-3*np.maximum(e[:,-1],1e-30)
 # Interval velocities are defined on actual physical time differences, no resampling.
 v=np.diff(x,axis=0)/np.diff(t)[:,None,None] if len(t)>1 else np.empty((0,32,3))
 vr=v-v.mean(1,keepdims=True);rad=centered[:-1];den=np.linalg.norm(rad,axis=-1)
 radial=np.sum(vr*rad,axis=-1)/np.maximum(den,1e-30);radial=np.where(den>1e-15,radial,np.nan)
 pairs=np.array([pdist(a) for a in x]);pair_xy=np.array([pdist(a[:,:2]) for a in x]);edge=pairs<=2.1*R
 nn=[];nnxy=[];components=[];contact_count=[];largest=[];labels=[];local_density=[]
 for a,aa,ed in zip(pairs,pair_xy,edge):
  d=squareform(a);np.fill_diagonal(d,np.inf);nn.append(np.median(d.min(1)))
  d=squareform(aa);np.fill_diagonal(d,np.inf);nnxy.append(np.median(d.min(1)));local_density.append((np.sum(d<=10e-6,axis=1)+1)/(np.pi*(10e-6)**2))
  ncomp,lab=connected_components(squareform(ed),directed=False);components.append(ncomp);labels.append(lab);largest.append(int(np.bincount(lab).max()));contact_count.append(int(np.sum(a<=2*R+1e-9)))
 return dict(**m,xy_tensor_m2=qxy,xy_rg_m=np.sqrt(np.trace(qxy,axis1=1,axis2=2)),xy_aspect=np.sqrt(e[:,-1]/np.maximum(e[:,0],1e-30)),xy_axis=vec[:,:,-1],xy_axis_identifiable=eigen_valid,pair_distance_m=pairs,xy_pair_distance_m=pair_xy,nn_m=np.array(nn),xy_nn_m=np.array(nnxy),pair_quantiles_m=np.quantile(pairs,[.1,.5,.9],axis=1).T,edges=edge,edge_count=edge.sum(1),components=np.array(components),largest_component=np.array(largest),component_labels=np.array(labels),contact_count=np.array(contact_count),velocity_m_s=v,relative_velocity_m_s=vr,interval_time_s=t[:-1],interval_dt_s=np.diff(t),radial_mean_m_s=np.nanmean(radial,axis=1),radial_median_m_s=np.nanmedian(radial,axis=1),pair_speed_m_s=np.diff(pairs,axis=0)/np.diff(t)[:,None],xy_density_m2=np.full(len(t),32/(120e-6)**2),xy_local_density_10um_m2=np.array(local_density))

def export(path,out):
 with h5py.File(path/'trajectory.h5') as h:
  t=h['time_s'][:];x=h['positions_m'][:];ids=h['particle_id'][:];attrs=dict(h.attrs)
 validate_trajectory(t,x,ids);m=series(t,x)
 with h5py.File(out,'w') as h:
  h.create_dataset('particle_id',data=ids)
  h.create_dataset('radius_m',data=np.full(32,R))
  for k,v in m.items():h.create_dataset(k,data=v,compression='gzip')
  h.create_dataset('positions_m',data=x,compression='gzip')
  h.create_dataset('accepted_state_valid',data=np.ones(len(t),bool))
  h.attrs['interval_velocity_definition']='Actual accepted position difference / actual dt; interval average, not instantaneous camera velocity'
  h.attrs['camera_calibrated']=False
  h.attrs['state_sources']='SIMULATION_3D_TRUTH;SIMULATION_XY_PROJECTION'
  h.attrs['imaging_status']='NOT_GENERATED; no calibrated PSF, depth, visibility or Chronos observation'
  h.attrs['xy_local_density_definition']='True projected centres within a 10 um disk, including self / disk area; not calibrated visible detections'
  h.attrs['flow']='STRICT_ZERO';h.attrs['raw_trajectory_sha256']=sha(path/'trajectory.h5')
  h.attrs['action_schedule_json']=attrs['schedule_json'];h.attrs['scope']='fixed 32 simulation IDs; projected density uses fixed 120 um square; visual candidate counts differ'
  h.attrs['requested_s']=json.loads((path/'RESULT.json').read_text())['requested_s']
  h.attrs['completed_s']=float(t[-1]);h.attrs['no_padding']=True
 return m
