"""Render-only QA using available accepted frames, never a new dynamics run."""
from render_videos import *
ref=load_npz(O/'provenance/reference.npz');field=load_npz(O/'HETEROGENEITY_FIELD.npz');states={}
for m in CASES:
 paths=sorted((O/'cases'/m/'keys').glob('frame_*.npz'));p=paths[-1];states[m]=dict(path=p,state=load_npz(p))
row=dict(requested_time_s=float(states['uniform']['state']['time']),**states);pts=np.concatenate([z['state']['x'] for z in states.values()])*1e6;lo=pts.min(0)-50;hi=pts.max(0)+50;bounds=np.array([lo[0],hi[0],lo[1],hi[1],lo[2],hi[2]]);center=(lo+hi)/2;span=np.max(hi-lo);camera=[center+np.array([1.4,-2.9,1.0])*span,center,[0,0,1]]
registries={}
for m,z in states.items():
 lab,sz,att=graph.topology(ref,z['state']['d']);registries[m]={tuple((np.flatnonzero(lab==label)+1).tolist()):int(label+1) for label in np.flatnonzero(~att)}
render(2,[row],ref,field,bounds,camera,registries,only_frame=True)
