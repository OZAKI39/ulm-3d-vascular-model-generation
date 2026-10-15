"""Post-run geometric audit; never changes the frozen binning or native state."""
from common import *
from itertools import product

r=dict(np.load(O/'provenance/reference.npz'));s=dict(np.load(O/'full/last_accepted.npz'));ph=np.load(O/'PLANE_HISTORY.npz')
ci=r['cross_ids'];alive=s['d'][ci]<1;xy=r['cross_intersections'][:,:2]
scaled=(xy+SIZE[:2]/2)/DX;nearest=np.rint(scaled);ties=np.abs(scaled-nearest)<1e-10
# Count each alive exact-boundary crossing in every incident cell. This is a
# conservative diagnostic, not a new assignment of the original interactions.
incident=np.zeros(COUNTS[:2],int)
for q,tie in zip(scaled[alive],ties[alive]):
 choices=[]
 for axis in range(2):
  n=int(round(q[axis]));bins=[n-1,n] if tie[axis] else [int(np.floor(q[axis]))]
  choices.append([v for v in bins if 0<=v<COUNTS[axis]])
 for ix,iy in product(*choices):incident[ix,iy]+=1
centers=ph['cell_centers_xy_m'].reshape(*COUNTS[:2],2);new=centers[:,:,0]>TIP
eligible=new&(incident==0)&(ph['opening_m'][-1]>=.1*DX);prefix=np.zeros_like(eligible)
tipcol=int(np.flatnonzero(new[:,0])[0])
for iy in range(COUNTS[1]):
 for ix in range(tipcol,COUNTS[0]):
  if not eligible[ix,iy]:break
  prefix[ix,iy]=True
assert np.all(~ph['conservative'][-1]|(incident==0)), 'Horizon-clear cells must also be free of boundary-tied bridges'
nat=~r['initial_cut'];failed=nat&(s['d']==1);mid=.5*(r['X'][r['bonds'][:,0]]+r['X'][r['bonds'][:,1]])
audit=dict(status='POSTPROCESS_GEOMETRIC_CHECK_COMPLETE',frozen_floor_assignment_unchanged=True,boundary_tolerance_cell_units=1e-10,boundary_tolerance_m=DX*1e-10,live_cross=int(alive.sum()),live_cross_on_grid_boundaries=int(np.sum(alive&np.any(ties,axis=1))),all_incident_cell_counts=incident,original_floor_empty_new_cells=int(np.sum(new&(ph['live_counts'][-1]==0))),all_incident_empty_new_cells=int(np.sum(new&(incident==0))),frozen_prefix_cells=int(ph['notch_connected'][-1].sum()),tie_robust_prefix_cells=int(prefix.sum()),tie_robust_prefix_grid_proxy_m2=float(prefix.sum()*DX**2),horizon_clear_subset_cells=int(ph['conservative'][-1].sum()),naturally_failed_off_plane=int(np.sum(failed&~r['cross'])),failed_midpoints_outside_plane_horizon=int(np.sum(failed&(abs(mid[:,2])>HORIZON))),interpretation='All-incident binning removes spurious empty cells from floating-point boundary ownership. It cannot certify a unique 3D crack surface or attribute coupled dissipation to a surface.')
dump(O/'GEOMETRIC_BRIDGING_REVIEW.json',audit)
area=json.loads((O/'CRACK_AREA_IDENTIFIABILITY.json').read_text());area['boundary_tie_review']='GEOMETRIC_BRIDGING_REVIEW.json';area['tie_robust_prefix_grid_proxy_m2']=audit['tie_robust_prefix_grid_proxy_m2'];area['frozen_prefix_lost_under_tie_review']=int(ph['notch_connected'][-1].sum()-prefix.sum());dump(O/'CRACK_AREA_IDENTIFIABILITY.json',area)
print(json.dumps({k:v for k,v in audit.items() if k!='all_incident_cell_counts'},default=cv,indent=2))
