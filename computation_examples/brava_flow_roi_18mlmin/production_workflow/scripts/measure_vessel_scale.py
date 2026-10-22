"""Measure transformed centerline sections and local cells, never infer diameter from wall distance."""
from pathlib import Path
import json,hashlib,csv,collections,itertools
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as PolygonPath
R=Path(__file__).resolve().parents[1];B=R.parent;OUT=R/'vessel_scale';OUT.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 swc=np.loadtxt(B/'inputs/fitted_centerline.swc');spec=json.loads((B/'inputs/geometry.json').read_text());T=np.array(json.loads((B/'inputs/print_transform.json').read_text())['transform_4x4']);xyz=(swc[:,2:5]@T[:3,:3].T+T[:3,3]-spec['origin_source_mm'])*1e-3
 ids={int(r[0]):i for i,r in enumerate(swc)};parents={i:ids[int(r[6])] for i,r in enumerate(swc) if r[6]>=0};adj=collections.defaultdict(list)
 for i,j in parents.items():adj[i].append(j);adj[j].append(i)
 branch=np.array([len(adj[i])>2 for i in range(len(swc))]);root=next(i for i in range(len(swc)) if i not in parents);s=np.zeros(len(swc));weights=np.zeros(len(swc))
 for i in range(len(swc)):
  k=i
  while k in parents:s[i]+=np.linalg.norm(xyz[k]-xyz[parents[k]]);k=parents[k]
 for i,j in parents.items():L=np.linalg.norm(xyz[i]-xyz[j]);weights[i]+=L/2;weights[j]+=L/2
 mesh=pv.read(B/'mesh/SV_MESH/mesh-complete.mesh.vtu');surface=pv.read(B/'inputs/labelled_core_m.vtp');tet=mesh.cells.reshape(-1,5)[:,1:];mp=np.asarray(mesh.points);edge=np.array(list(itertools.combinations(range(4),2)));h=np.linalg.norm(mp[tet[:,edge[:,0]]]-mp[tet[:,edge[:,1]]],axis=2).mean(1)
 cell=mesh.find_containing_cell(xyz);rows=[]
 for i,x in enumerate(xyz):
  children=[j for j in adj[i] if j!=parents.get(i)]
  if i in parents:
   tangent=(xyz[children[0]]-xyz[parents[i]]) if len(children)==1 else x-xyz[parents[i]]
  else:tangent=xyz[children[0]]-x
  tangent/=np.linalg.norm(tangent);ref=np.eye(3)[np.argmin(abs(tangent))];a=np.cross(tangent,ref);a/=np.linalg.norm(a);b=np.cross(tangent,a)
  cut=surface.slice(normal=tangent,origin=x);pts=np.asarray(cut.points);links=collections.defaultdict(set);arr=cut.lines;k=0
  while k<len(arr):
   n=int(arr[k]);vs=arr[k+1:k+1+n];k+=n+1
   for u,v in zip(vs[:-1],vs[1:]):links[int(u)].add(int(v));links[int(v)].add(int(u))
  seen=set();loops=[]
  for start in links:
   if start in seen:continue
   component=set();todo=[start]
   while todo:
    v=todo.pop()
    if v in component:continue
    component.add(v);todo.extend(links[v]-component)
   seen|=component
   if not all(len(links[v])==2 for v in component):continue
   ordered=[start];prev=None;cur=start
   while True:
    nxt=next(v for v in links[cur] if v!=prev)
    if nxt==start:break
    ordered.append(nxt);prev,cur=cur,nxt
    if len(ordered)>len(component):break
   if len(ordered)!=len(component):continue
   uv=np.column_stack([(pts[ordered]-x)@a,(pts[ordered]-x)@b]);area=abs(np.sum(uv[:,0]*np.roll(uv[:,1],-1)-uv[:,1]*np.roll(uv[:,0],-1)))/2
   if area>0 and PolygonPath(uv).contains_point([0,0]):loops.append((area,uv))
  # Choose the enclosing local lumen, not remote vessel/tube intersections.
  area=min((q[0] for q in loops),default=np.nan);D=2*np.sqrt(area/np.pi);local_h=float(h[cell[i]]) if cell[i]>=0 else np.nan
  rows.append(dict(node_id=int(swc[i,0]),parent_id=int(swc[i,6]),distance_from_root_m=float(s[i]),x_m=float(x[0]),y_m=float(x[1]),z_m=float(x[2]),length_weight_m=float(weights[i]),model_radius_m=float(swc[i,5]*1e-3),diameter_area_m=float(D),model_diameter_m=float(2*swc[i,5]*1e-3),section_valid=bool(np.isfinite(D)),branch_node=bool(branch[i]),cell_id=int(cell[i]),local_mean_tet_edge_m=local_h,N_across=float(D/local_h),section_vs_model_relative=float(D/(2*swc[i,5]*1e-3)-1)))
 vals=np.array([r['diameter_area_m'] for r in rows]);valid=np.isfinite(vals);assert valid.sum()>.85*len(rows),(valid.sum(),len(rows))
 def quantiles(v,w):
  order=np.argsort(v);v=v[order];w=w[order];x=(np.cumsum(w)-w/2)/w.sum();return np.interp([0,.05,.25,.5,.75,.95,1],x,v)
 keys=['minimum','P5','P25','median','P75','P95','maximum'];distribution=dict(zip(keys,map(float,quantiles(vals[valid],weights[valid]))));radii=dict(R1_small_m=distribution['P5']/2,R2_medium_m=distribution['median']/2)
 # Artificial tubes are reported separately; their known port radii are not anatomical samples.
 ports={k:2*v['radius_m'] for k,v in spec['ports'].items()};na=np.array([r['N_across'] for r in rows]);good=np.isfinite(na)
 summary=dict(status='PASS_MEASURED_GEOMETRY',method='Closed cross-section area containing transformed fitted centerline; D_eq=2*sqrt(A/pi). Length-weighted distribution over existing centerline nodes. No nearest-wall distance or bounding-box diameter.',samples=len(rows),valid_sections=int(valid.sum()),missing_node_ids=[r['node_id'] for r in rows if not r['section_valid']],diameter_distribution_m=distribution,model_radius_source='inputs/fitted_centerline.swc',radius_for_pipe=radii,artificial_port_diameters_m=ports,current_N_across_distribution=dict(zip(keys,map(float,quantiles(na[good],weights[good])))),current_h_definition='Mean of six edges in the tetra containing the centerline sample; D/h is a resolution proxy, not an actual wall-normal layer count.',branch_node_ids=[r['node_id'] for r in rows if r['branch_node']],small_vessel_node_ids=[r['node_id'] for r in rows if r['diameter_area_m']<=distribution['P25']],source_sha256={rel:sha(B/rel) for rel in ['inputs/fitted_centerline.swc','inputs/print_transform.json','inputs/geometry.json','inputs/labelled_core_m.vtp','mesh/SV_MESH/mesh-complete.mesh.vtu']},interpretation='Manufactured compensated human ROI lumen, not unmodified anatomical physiology. Junction sections may be noncircular; area-equivalent D is reported explicitly. Engineered ports are separated.',nearwall_hotspots='Historical cells 388769 and 138584 and frozen slow-fast cohort remain required downstream audit locations; D/h does not estimate particle occupancy.')
 (OUT/'SCALE.json').write_text(json.dumps(summary,indent=2)+'\n')
 with (OUT/'centerline_sections.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 np.savez_compressed(OUT/'centerline_scale.npz',points_m=xyz,diameter_m=vals,model_radius_m=swc[:,5]*1e-3,parents=np.array(list(parents.items())),length_weights_m=weights,branch=branch,N_across=na)
 plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'figure.constrained_layout.use':True,'pdf.fonttype':42})
 fig,ax=plt.subplots(1,2,figsize=(12,4.5));ax[0].hist(vals[valid]*1e3,bins=18,weights=weights[valid]*1e3,color='#327a9f');ax[0].set(xlabel='Area-equivalent lumen diameter (mm)',ylabel='Represented centerline length (mm)',title='Manufactured BraVa ROI | tubes separate')
 colors=ax[1].scatter(s[valid]*1e3,vals[valid]*1e3,c=na[valid],s=22,cmap='viridis');fig.colorbar(colors,ax=ax[1],label='Current D / h (proxy)');ax[1].scatter(s[branch]*1e3,vals[branch]*1e3,marker='x',c='red',label='Branch nodes');ax[1].set(xlabel='Centerline distance from root (mm)',ylabel='Diameter (mm)',title='All branches share root distance');ax[1].legend()
 fig.savefig(R/'figures/vessel_diameter_distribution.png',dpi=300);fig.savefig(R/'figures/vessel_diameter_distribution.pdf');plt.close(fig)
 fig=plt.figure(figsize=(10,7));ax=fig.add_subplot(projection='3d');c=ax.scatter(*xyz.T*1e3,c=vals*1e3,s=20,cmap='viridis')
 for i,j in parents.items():ax.plot(*xyz[[i,j]].T*1e3,c='.5',lw=.7)
 ax.set(xlabel='x (mm)',ylabel='y (mm)',zlabel='z (mm)',title='Measured lumen diameter along fitted centerline');fig.colorbar(c,ax=ax,label='Area-equivalent diameter (mm)',shrink=.7);fig.savefig(R/'figures/centerline_diameter.png',dpi=300);fig.savefig(R/'figures/centerline_diameter.pdf');plt.close(fig)
 lines=['# BraVa 血管尺度测量','', '状态：**PASS_MEASURED_GEOMETRY**。使用现有拟合中心线的刚体变换，逐节点切实际 lumen 表面，选择包含中心点的闭合截面；报告面积等效直径，不用最近 WALL 距离充当直径。','', '| 统计量（按中心线长度加权） | 直径 mm | 当前 D/h |','|---|---:|---:|']
 for k in keys:lines.append(f"| {k} | {distribution[k]*1e3:.6f} | {summary['current_N_across_distribution'][k]:.3f} |")
 lines+=['',f'截面有效 {valid.sum()}/{len(rows)}；未闭合节点 {summary["missing_node_ids"]} 明确保留。h 为球心所在 tetra 六边平均长度，D/h 是近似分辨率，不等同真实层数。分叉用面积等效直径，不能假设圆形。', '',f'直管代表半径：R1=P5/2={radii["R1_small_m"]*1e6:.3f} µm；R2=median/2={radii["R2_medium_m"]*1e6:.3f} µm。', '', '这是工程补偿与打印几何对应的 human ROI，不是未经修改的生理解剖直径。人工接管单独记录在 SCALE.json，未混入原 ROI 分布。小血管与分叉节点标记见 CSV；当前 D/h 不能证明近壁精度或粒子停留概率。', '', '[分布图](figures/vessel_diameter_distribution.png) · [沿中心线直径](figures/centerline_diameter.png) · [全部截面数据](vessel_scale/centerline_sections.csv) · [测量身份](vessel_scale/SCALE.json)']
 (R/'BRAVA_VESSEL_SCALE_REPORT_ZH.md').write_text('\n'.join(lines)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
