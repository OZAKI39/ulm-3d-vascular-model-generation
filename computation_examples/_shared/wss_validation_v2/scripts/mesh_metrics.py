"""SI mesh/geometry summaries, with surface area weights and explicit region units."""
from pathlib import Path
import sys,json,csv
import numpy as np
from case_common import *

def summarize(case):
 case=Path(case);a=np.load(case/'SV_MESH/mesh_arrays.npz');x=a['points_m'];t=a['tetra'];b=a['boundary_triangles'];tags=a['facet_tags'];own=a['adjacent_tetra'];q=a['min_sicn'];vol=np.linalg.det(x[t[:,1:]]-x[t[:,:1]])/6;xyz=x[b];area=.5*np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1);h=3*vol[own]/area
 rows=[]
 for tag in [None,*np.unique(tags)]:
  sel=np.ones(len(b),bool) if tag is None else tags==tag;aa=area[sel];hh=h[sel];qq=q[own[sel]];i=np.argsort(hh);hs=np.interp([.05,.5,.95],(np.cumsum(aa[i])-.5*aa[i])/aa.sum(),hh[i])
  row=dict(case=case.name,mesh=case.name,region='all_exterior' if tag is None else ROLES[int(tag)],statistical_weight='area_for_surface_height;cell_count_for_quality',nodes=len(x),tetra=len(t),facets=int(sel.sum()),volume_um3=float(vol.sum()*1e18),area_um2=float(aa.sum()*1e12),equivalent_radius_um=float(np.sqrt(aa.sum()/np.pi)*1e6) if tag is not None and tag!=1 else '',nearwall_height_mean_um=float(np.average(hh,weights=aa)*1e6),height_p05_um=float(hs[0]*1e6),height_p50_um=float(hs[1]*1e6),height_p95_um=float(hs[2]*1e6),owner_minSICN_min=float(qq.min()),owner_minSICN_p05=float(np.quantile(qq,.05)),owner_minSICN_median=float(np.median(qq)),global_minSICN_min=float(q.min()),global_minSICN_p05=float(np.quantile(q,.05)),global_minSICN_median=float(np.median(q)),nonpositive_tetra=int(np.sum(vol<=0)),geometry_reference='explicit case mesh; compare identical physical source geometry')
  rows.append(row)
 (case/'reports').mkdir(exist_ok=True)
 with (case/'reports/mesh_geometry.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 return rows
if __name__=='__main__':
 rows=[]
 for arg in sys.argv[1:]:rows.extend(summarize(arg))
 rows=[]
 for p in sorted(V.glob('stage[234]/*/reports/mesh_geometry.csv')):
  rows.extend(list(csv.DictReader(p.open())))
 out=V/'data/mesh_geometry.csv'
 with out.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 print(json.dumps({'mesh_cases':len(sys.argv)-1,'region_rows':len(rows)},indent=2))
