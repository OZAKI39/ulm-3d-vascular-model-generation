from pathlib import Path
import numpy as np,json,collections
S=Path(__file__).resolve().parents[1]
g=dict(np.load(S/'geometry/GEOMETRY_ARRAYS.npz'));faces=g['faces'][g['classes']<=2];regions=g['classes'][g['classes']<=2];coords=g['points'][faces]
vertices,inv=np.unique(coords.reshape(-1,3),axis=0,return_inverse=True);f=inv.reshape(-1,3);cross=np.cross(coords[:,1]-coords[:,0],coords[:,2]-coords[:,0]);area=np.linalg.norm(cross,axis=1)/2;n=cross/(2*area[:,None]);edges=collections.defaultdict(list)
for i,t in enumerate(f):
 for k in range(3):
  a,b=int(t[k]),int(t[(k+1)%3]);edges[tuple(sorted((a,b)))].append((i,1 if a<b else -1))
adj=[[] for _ in f];angles=[];seams=collections.Counter();flips=0;bound=0;nonman=[];elist=[]
for (a,b),incident in edges.items():
 if len(incident)==1:bound+=1
 elif len(incident)>2:nonman.append(dict(vertices=[a,b],faces=[x[0] for x in incident]))
 else:
  (i,s),(j,t)=incident;adj[i].append(j);adj[j].append(i);angle=float(np.degrees(np.arccos(np.clip(n[i]@n[j],-1,1))));angles.append(angle);flips+=s==t;elist.append([a,b,i,j,angle])
  if regions[i]!=regions[j]:seams[tuple(sorted((int(regions[i]),int(regions[j]))))]+=1
seen=set();components=[]
for i in range(len(f)):
 if i in seen:continue
 todo=[i];seen.add(i);count=0
 while todo:
  j=todo.pop();count+=1
  for k in adj[j]:
   if k not in seen:seen.add(k);todo.append(k)
 components.append(count)
angles=np.array(angles);np.savez_compressed(S/'raw/TOPOLOGY_ARRAYS.npz',vertices=vertices,faces=f,normals=n,area=area,regions=regions,manifold_edges=np.array(elist))
result=dict(status='PASS' if not nonman and not flips and np.all(area>0) else 'REQUIRES_REVIEW',triangle_count=len(f),vertex_count=len(vertices),edge_count=len(edges),manifold_edge_count=len(angles),boundary_edge_count=bound,nonmanifold_edge_count=len(nonman),nonmanifold_edges=nonman,connected_components=len(components),component_triangles=sorted(components,reverse=True),orientation_inconsistent_edges=int(flips),degenerate_triangles=int((area<=0).sum()),dihedral_deg_percentiles=dict(zip(['min','p25','p50','p75','p90','p95','p99','max'],np.percentile(angles,[0,25,50,75,90,95,99,100]).tolist())),region_seams={str(k):v for k,v in seams.items()},region_counts={str(i):int((regions==i).sum()) for i in range(3)},ports_excluded=True,production_geometry_modified=False)
(S/'validation/WALL_SURFACE_TOPOLOGY_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
# Engine-side region annotations index exact original WALL_ONLY triangles.
(S/'geometry/WALL_TRIANGLE_REGIONS.txt').write_text('\n'.join(map(str,regions.tolist()))+'\n')
