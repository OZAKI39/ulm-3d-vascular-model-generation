from pathlib import Path
import numpy as np,json,trimesh,csv
ROOT=Path(__file__).resolve().parents[1]
a=json.loads((ROOT/'inputs/port_attachment_summary.json').read_text())
g=json.loads((ROOT/'inputs/geometry.json').read_text())
swc=np.loadtxt(ROOT/'inputs/fitted_centerline.swc')
T=np.array(json.loads((ROOT/'inputs/print_transform.json').read_text())['transform_4x4'])
xyz=swc[:,2:5]@T[:3,:3].T+T[:3,3]
index={int(row[0]):i for i,row in enumerate(swc)}
parents={i:index[int(row[6])] for i,row in enumerate(swc) if row[6]>=0}
children={i:[] for i in range(len(swc))}
for i,j in parents.items():children[j].append(i)
root=next(i for i in children if i not in parents);leaves=[i for i in children if not children[i]]
assert len(leaves)==3
mu=g['material']['dynamic_viscosity_Pa_s'];flows={};paths={};tips={};tubes={}
for name,v in a['per_port'].items():
 ring=np.array(v['cap']['aligned_ring_print_mm']);c=ring.mean(0)
 pool=[root] if name=='I1' else leaves
 tip=min(pool,key=lambda i:np.linalg.norm(xyz[i]-c));err=np.linalg.norm(xyz[tip]-c)
 assert err<.2,(name,err)
 tips[name]=dict(swc_node=int(swc[tip,0]),position_source_mm=xyz[tip].tolist(),cap_distance_mm=float(err))
 if name!='I1':
  path=[];i=tip
  while i!=root:path.append(i);i=parents[i]
  paths[name]=path
  for i in path:flows[i]=flows.get(i,0)+1e-7
 src=Path(v['port_stl'].replace('/home/lzy/projects/ulm_3D_vascular/','/home/lzy/projects/vascular_printing/'))
 tube=trimesh.load_mesh(src,process=True)
 area=g['ports'][name]['area_m2'];length=tube.volume*1e-9/area
 radius=np.sqrt(area/np.pi);R=8*mu*length/(np.pi*radius**4)
 tubes[name]=dict(length_m=float(length),R_Pa_s_m3=float(R),source=str(src),method='constant-section swept port volume divided by measured cap area; approximate collar contribution retained')
R={};rows=[]
for i,j in parents.items():
 length=np.linalg.norm(xyz[i]-xyz[j])*1e-3
 radii=np.linspace(swc[j,5],swc[i,5],101)*1e-3
 resistance=8*mu/np.pi*length*np.trapezoid(1/radii**4,dx=.01)
 R[i]=resistance
 rows.append(dict(child=int(swc[i,0]),parent=int(swc[j,0]),length_m=length,r_start_m=radii[0],r_end_m=radii[-1],R_Pa_s_m3=resistance,Q_target_m3_s=flows[i],drop_Pa=resistance*flows[i]))
drops={k:sum(R[i]*flows[i] for i in path)+tubes[k]['R_Pa_s_m3']*1e-7 for k,path in paths.items()}
pressure={k:max(drops.values())-v for k,v in drops.items()}
result=dict(method='Poiseuille resistance along fitted centerline plus constant-area engineered tube estimate; preliminary design, not a 3D flow result',
 Q_in_m3_s=3e-7,Q_out_target_m3_s=1e-7,outlet_pressures_Pa=pressure,inlet_estimate_Pa=max(drops.values())+tubes['I1']['R_Pa_s_m3']*3e-7,
 ports=tips,tubes=tubes,root_to_outlet_drop_Pa=drops,uniform_target_fractions={k:1/3 for k in paths},
 limitations=['Inertia and three-dimensional junction/curvature effects not captured by this resistance estimate.','Final boundary pressures and achieved splits require actual CFD verification.'])
(ROOT/'inputs/boundary_design.json').write_text(json.dumps(result,indent=2)+'\n')
with (ROOT/'reports/network_edges.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps(result,indent=2))
