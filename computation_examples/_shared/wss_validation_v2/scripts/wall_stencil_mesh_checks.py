"""Detect P1 tetrahedra constrained to zero velocity at all four vertices."""
from pathlib import Path
import csv
import numpy as np
V=Path(__file__).resolve().parents[1];rows=[]
for name in ['vessel_baseline','vessel_medium','vessel_fine']:
 m=np.load(V/'stage3'/name/'SV_MESH/mesh_arrays.npz');x=m['points_m'];t=m['tetra'];b=m['boundary_triangles'];wall=b[m['facet_tags']==1];flag=np.zeros(len(x),bool);flag[np.unique(wall)]=True;count=flag[t].sum(axis=1)
 for k in range(5):rows.append(dict(case=name,region='whole_volume',wall_constrained_vertices_per_tetra=k,tetra_count=int(np.count_nonzero(count==k)),fraction_by_tetra_count=float(np.mean(count==k)),unit='count_or_fraction',reference='actual tag1 no-slip wall; four constrained vertices imply identically zero P1 velocity',statistical_weight='tetra_count'))
with (V/'data/wall_constrained_tetra_counts.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print([r for r in rows if r['wall_constrained_vertices_per_tetra']==4])
