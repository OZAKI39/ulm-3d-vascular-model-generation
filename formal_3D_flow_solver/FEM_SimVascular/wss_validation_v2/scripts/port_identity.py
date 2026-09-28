"""Area-weighted actual cap coordinates/normals and XML pressure labels."""
import csv,json,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
V=Path(__file__).resolve().parents[1];roles={2:'OUTLET_03',3:'OUTLET_01',4:'INLET',5:'OUTLET_02'};rows=[]
for name in ['vessel_baseline','vessel_medium','vessel_fine']:
 c=V/'stage3'/name;m=np.load(c/'SV_MESH/mesh_arrays.npz');x=m['points_m'];b=m['boundary_triangles'];tags=m['facet_tags'];root=ET.parse(c/'run/solver.xml')
 for tag,role in roles.items():
  tri=x[b[tags==tag]];av=.5*np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);area=np.linalg.norm(av,axis=1);n=av.sum(axis=0);n/=np.linalg.norm(n);center=np.average(tri.mean(axis=1),axis=0,weights=area)*1e6;bc=root.find('.//Add_BC[@name="'+role+'"]')
  rows.append(dict(case=name,facet_tag=tag,boundary=role,center_x_um=center[0],center_y_um=center[1],center_z_um=center[2],outward_normal_x=n[0],outward_normal_y=n[1],outward_normal_z=n[2],area_um2=area.sum()*1e12,BC_type=bc.find('Type').text,XML_value=float(bc.find('Value').text),XML_value_unit='m3/s' if role=='INLET' else 'Pa',statistical_weight='actual_cap_triangle_area',reference='same_case_mesh_and_XML'))
with (V/'data/port_identity.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps(rows[:4],indent=2))
