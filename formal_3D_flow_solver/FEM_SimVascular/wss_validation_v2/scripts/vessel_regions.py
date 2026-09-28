"""Frozen physical region definitions; no cross-grid node-ID matching."""
import json
import numpy as np
from pathlib import Path
V=Path(__file__).resolve().parents[1];A=V.parent/'wss_audit'
PORTS=json.loads((A/'inputs/network/roi_ports_in_a.json').read_text())['ports']
JUMP=np.array([90.9159025468,48.5593069849,112.109730253])*1e-6
CENTERS={'J1_r5um':[92,49,111],'J2_r5um':[130.04,82.04,87.18],'bend_r5um':[154.34,94.58,83.0],'worst_tet_r5um':[98.658078,42.075917,131.953726]}

def masks(centers):
 out={'all_wall':np.ones(len(centers),bool)}
 for name,xyz in CENTERS.items():out[name]=np.linalg.norm(centers-np.array(xyz)*1e-6,axis=1)<5e-6
 extension=np.zeros(len(centers),bool)
 for p in PORTS:
  origin=np.array(p['real_cut_xyz_um'])*1e-6;n=np.array(p['outward_normal']);s=(centers-origin)@n;r=np.linalg.norm(centers-origin-s[:,None]*n,axis=1)
  m=(s>0)&(s<p['fem_cap']['real_cut_to_cap_axial_um']*1e-6)&(r<3*p['radius_um']*1e-6)
  out[p['name']+'_extension']=m;extension|=m
 out['real_ROI']=~extension;out['artificial_extensions']=extension
 return out

def path_coordinates(centers):
 path=json.loads((A/'data/outlet2_path.json').read_text());xyz=np.array(path['path_xyz_um'])*1e-6
 g=np.load(A/'inputs/network/analysis_A_H0_graph_si.npz');ids=g['ids'];radii=[float(g['radius_m'][np.flatnonzero(ids==j)[0]]) for j in path['node_ids']];radii.append(next(p['fem_cap']['equivalent_radius_um']*1e-6 for p in PORTS if p['name']=='O2'))
 radii=np.array(radii);vec=xyz[1:]-xyz[:-1];length=np.linalg.norm(vec,axis=1);arc=np.r_[0,np.cumsum(length)];sbest=np.zeros(len(centers));dbest=np.full(len(centers),np.inf);rbest=np.zeros(len(centers))
 for i,(start,v,L) in enumerate(zip(xyz[:-1],vec,length)):
  f=np.clip((centers-start)@v/L**2,0,1);dist=np.linalg.norm(centers-start-f[:,None]*v,axis=1);sel=dist<dbest
  sbest[sel]=arc[i]+f[sel]*L;dbest[sel]=dist[sel];rbest[sel]=(1-f[sel])*radii[i]+f[sel]*radii[i+1]
 return sbest,dbest<2.5*rbest,arc[-1]

def section_definitions():
 rows=[]
 for p in PORTS:
  n=np.array(p['outward_normal']);origin=np.array(p['real_cut_xyz_um'])*1e-6;L=p['fem_cap']['real_cut_to_cap_axial_um']*1e-6
  for f in [.25,.5]:rows.append(dict(name=p['name']+'_extension_'+str(f),origin_m=(origin+f*L*n).tolist(),normal=n.tolist(),radius_limit_m=3*p['radius_um']*1e-6,compare_boundary=p['name'],expected_sign=1))
 # True ROI trunk between the inlet cut and first junction, oriented inlet -> J1.
 g=np.load(A/'inputs/network/analysis_A_H0_graph_si.npz');ids=g['ids'];edges=g['edges'][g['roi_internal_edge_mask']];xyz=g['xyz_m'];r=g['radius_m'];adj={}
 for i,j in edges:adj.setdefault(int(i),[]).append(int(j));adj.setdefault(int(j),[]).append(int(i))
 root=int(np.flatnonzero(ids==-10001)[0]);parent={root:-1};queue=[root]
 for i in queue:
  for j in adj[i]:
   if j not in parent:parent[j]=i;queue.append(j)
 def path(a,b):
  cur=int(np.flatnonzero(ids==b)[0]);stop=int(np.flatnonzero(ids==a)[0]);out=[cur]
  while out[-1]!=stop:out.append(parent[out[-1]])
  return out[::-1]
 for label,a,b,bound,sgn in [('inlet_trunk',-10001,3238,'INLET',-1),('J1_to_J2',3238,3274,'O1_plus_O3',1)]:
  ii=path(a,b);xx=xyz[ii];length=np.linalg.norm(np.diff(xx,axis=0),axis=1);arc=np.r_[0,np.cumsum(length)];k=min(np.searchsorted(arc,.5*arc[-1])-1,len(length)-1);f=(.5*arc[-1]-arc[k])/length[k];origin=xx[k]+f*(xx[k+1]-xx[k]);n=(xx[k+1]-xx[k])/length[k]
  rows.append(dict(name=label,origin_m=origin.tolist(),normal=n.tolist(),radius_limit_m=float(3*((1-f)*r[ii[k]]+f*r[ii[k+1]])),compare_boundary=bound,expected_sign=sgn))
 return rows
