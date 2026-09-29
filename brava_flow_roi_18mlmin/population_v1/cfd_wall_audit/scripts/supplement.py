"""Additional explicit diagnostics on audit arrays; no new CFD/trajectory run."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import pyvista as pv
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr
from audit_field import js
ap=argparse.ArgumentParser();ap.add_argument('--population-root',type=Path,required=True);root=ap.parse_args().population_root;out=root/'cfd_wall_audit'
z=np.load(out/'recomputed_arrays.npz');f=pd.read_csv(out/'all_cells.csv');wf=pd.read_csv(out/'wall_face_structure.csv');ident=json.loads((out/'AUDIT_IDENTITY.json').read_text());c=z['connectivity'];xyz=z['xyz'];u=z['velocity'];p=z['pressure'];adj=z['adjacency'];N=len(c);owner=z['owner'];area=z['wall_area'];layer=z['layer'];eta=z['eta'];vol=z['volume'];G=z['G'];X=xyz[c];cent=z['centroid']
ports=[];all_boundary_keys=[]
for name,path in ident['input_paths'].items():
 if name in ['mesh','flow']:continue
 b=pv.read(path);ids=np.asarray(b['GlobalNodeID'],dtype=np.int64)-1;faces=ids[np.asarray(b.faces).reshape(-1,4)[:,1:]];tri=xyz[faces];raw_n=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])/2
 # Surface GlobalElementID is verified against the containing canonical tetra.
 parent=np.asarray(b.cell_data['GlobalElementID'],dtype=np.int64)-1
 assert np.all([len(set(face)&set(c[t]))==3 for face,t in zip(faces,parent)])
 raw_n*=np.where(np.einsum('ij,ij->i',raw_n,cent[parent]-tri.mean(1))>0,-1,1)[:,None]
 q=np.einsum('ij,ij->i',u[faces].mean(1),raw_n).sum()
 ports.append(dict(boundary=name,triangle_count=len(faces),area_m2=np.linalg.norm(raw_n,axis=1).sum(),flow_outward_m3_s=q,mean_node_pressure_Pa=float(p[ids].mean()),center_x_m=tri.mean((0,1))[0],center_y_m=tri.mean((0,1))[1],center_z_m=tri.mean((0,1))[2]))
 all_boundary_keys.extend([tuple(sorted(face)) for face in faces])
pd.DataFrame(ports).to_csv(out/'boundary_identity_flux.csv',index=False)
faces4=np.array([[1,2,3],[0,2,3],[0,1,3],[0,1,2]]);bi,bf=np.where(adj<0);expected=set(map(tuple,np.sort(c[bi[:,None],faces4[bf]],axis=1)))
assert len(all_boundary_keys)==len(set(all_boundary_keys)) and set(all_boundary_keys)==expected
# Moderate screening components eta >= .1 across W0-W2.
flag=(layer<=2)&(eta>=.1);ids=np.flatnonzero(flag);idx=np.full(N,-1);idx[ids]=np.arange(len(ids));aa=np.repeat(ids,4);bb=adj[ids].ravel();keep=(bb>=0);aa=aa[keep];bb=bb[keep];keep=flag[bb];aa=aa[keep];bb=bb[keep]
graph=coo_matrix((np.ones(len(aa)),(idx[aa],idx[bb])),shape=(len(ids),len(ids))).tocsr();n,lab=connected_components(graph,directed=False);sizes=np.bincount(lab);order=np.argsort(-sizes);rank=np.empty(n,dtype=int);rank[order]=np.arange(n);labels=np.full(N,-1);labels[ids]=rank[lab];rows=[]
for k in range(n):
 ci=ids[labels[ids]==k];wi=np.flatnonzero(labels[owner]==k);bounds=np.stack([X[ci].min((0,1)),X[ci].max((0,1))])
 rows.append(dict(cluster_id=k,cell_count=len(ci),volume_m3=vol[ci].sum(),wall_area_m2=area[wi].sum(),eta_median=float(np.median(eta[ci])),eta_max=float(eta[ci].max()),bbox_min_m=json.dumps(bounds[0].tolist()),bbox_max_m=json.dumps(bounds[1].tolist())))
pd.DataFrame(rows).to_csv(out/'moderate_eta_clusters.csv',index=False);pd.DataFrame(dict(tetra_id=ids,cluster_id=labels[ids])).to_csv(out/'moderate_eta_cluster_cells.csv',index=False)
wi=np.flatnonzero(labels[owner]>=0);pd.DataFrame(dict(wall_triangle_id=wi,cluster_id=labels[owner[wi]])).to_csv(out/'moderate_eta_cluster_wall_triangles.csv',index=False)
# Spatial octants, not anatomy labels: count-weighted flagged fractions in each equal bounding-box half.
mid=(cent.min(0)+cent.max(0))/2;octant=np.sum((cent>=mid)*np.array([1,2,4]),axis=1);sp=[]
for l in range(3):
 for k in range(8):
  mask=(layer==l)&(octant==k)
  sp.append(dict(layer=f'W{l}',bbox_octant=k,cell_count=int(mask.sum()),eta_ge_0p1=int((mask&(eta>=.1)).sum()),eta_ge_0p5=int((mask&(eta>=.5)).sum())))
pd.DataFrame(sp).to_csv(out/'spatial_octants.csv',index=False)
curv=wf.max_adjacent_wall_normal_angle_deg.to_numpy();phi=wf.interior_vertex_normal_m_s.to_numpy();h=wf.interior_vertex_height_m.to_numpy();wall_eta=eta[owner]
rows=[]
for name,mask in [('all',np.ones(len(wf),bool)),('eta_ge_0p1',wall_eta>=.1),('eta_ge_0p5',wall_eta>=.5)]:
 for key in ['interior_vertex_speed_m_s','interior_vertex_height_m','wall_face_node_speed_max_m_s','max_adjacent_wall_normal_angle_deg','hmax_m']:
  rows.append(dict(region=name,metric=key,**dict(zip(['min','median','P95','max'],np.quantile(wf.loc[mask,key],[0,.5,.95,1])))))
pd.DataFrame(rows).to_csv(out/'wall_structure_distributions.csv',index=False)
q=[r['flow_outward_m3_s'] for r in ports];qi=-next(r['flow_outward_m3_s'] for r in ports if r['boundary']=='INLET')
summary=dict(port_coverage_all_exterior_exact=True,external_face_count=len(expected),boundary_sum_m3_s=sum(q),boundary_relative_imbalance=abs(sum(q))/qi,
 pressure_range_Pa=[float(p.min()),float(p.max())],velocity_node_speed_range_m_s=[float(np.linalg.norm(u,axis=1).min()),float(np.linalg.norm(u,axis=1).max())],
 canonical_unique_node_coordinates=len(np.unique(xyz,axis=0)),canonical_node_count=len(xyz),
 P1_shared_face_pressure_and_velocity_continuity='Exact by shared three canonical nodal values; no discontinuous duplicate flow field at the same face',
 W0_identity_eta_equals_normal_fraction_max_abs=float(np.max(abs(f.loc[layer==0,'eta_div']-f.loc[layer==0,'normal_fraction']))),
 W0_identity_div_equals_inward_u_over_height_max_abs_s_inv=float(np.max(abs(z['div'][owner]-phi/h))),
 eta_vs_local_normal_angle_spearman=float(spearmanr(wall_eta,curv).statistic),moderate_eta_cluster_count=n,
 max_moderate_eta_cluster=rows[0] if False else pd.read_csv(out/'moderate_eta_clusters.csv').iloc[0].to_dict(),
 bbox_octant_split_m=mid,warning='Centroid wall-normal metric for W1/W2 uses a nearest triangle; normal is not unique at shared edges. Primary W0 own-face normal is unambiguous.')
js(out/'SUPPLEMENT_SUMMARY.json',summary);print(json.dumps(summary,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item()),flush=True)
