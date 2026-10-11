"""Actual-native-field measurements with explicit coordinate-format identity.
Native SV writes vtkPoints float32 and velocity/pressure float64. Bijection to
frozen input doubles is accepted only on exact float32 roundtrip equality.
All integrals use the actual double input mesh. No interpolation/extrapolation.
"""
from analysis.dual_track_v1.common import *
from vascular_validation.postprocess import SolutionMeasurements
from vascular_validation.validation import triangle_flux,parse_result_vtu
from scipy.spatial import cKDTree
class Measurements(SolutionMeasurements):
 def __init__(self,*args):
  super().__init__(*args);self.wall_nodes=np.unique(self.boundary[np.isin(self.tags,[1,5])]);self.identity=[]
 def read(self,path):
  grid,u,p=parse_result_vtu(path);assert grid.n_points==len(self.points) and grid.n_cells==len(self.tetra)
  dist,idx=cKDTree(np.asarray(grid.points,float)).query(self.points);assert len(np.unique(idx))==len(idx)
  assert np.array_equal(np.asarray(grid.points)[idx].astype(np.float32),self.points.astype(np.float32)),'Native coordinate identity not an exact float32 roundtrip'
  assert np.asarray(u).dtype==np.float64 and np.asarray(p).dtype==np.float64,'Native velocity/pressure precision changed'
  self.identity.append(dict(snapshot=str(path),max_coordinate_roundoff_m=float(dist.max()),bijection=True,exact_float32_roundtrip=True,velocity_dtype=str(np.asarray(u).dtype),pressure_dtype=str(np.asarray(p).dtype)))
  return np.asarray(u,float)[idx],np.asarray(p,float).reshape(-1)[idx]
 def measure(self,u,p):
  flows={};pressure={};backflow={}
  for tag,name in [(2,'INLET'),(3,'OUTLET_A'),(4,'OUTLET_B')]:
   tri=self.boundary[self.tags==tag];xyz=self.points[tri];av=.5*np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]);area=np.linalg.norm(av,axis=1);flow=np.einsum('ij,ij->i',u[tri].mean(1),av);flows[name]=float(flow.sum());pressure[name]=float(np.dot(area,p[tri].mean(1))/area.sum());backflow[name]=float(area[flow<0].sum()/area.sum())
  incoming=-flows['INLET'];outlets={n:flows[n] for n in ['OUTLET_A','OUTLET_B']};total=sum(outlets.values());wall=float(np.linalg.norm(u[self.wall_nodes],axis=1).max());clot=np.unique(self.boundary[self.tags==5]);clotmax=float(np.linalg.norm(u[clot],axis=1).max()) if len(clot) else None
  return dict(Q_target_m3_s=self.Q,Q_in_m3_s=incoming,Q_out_total_m3_s=total,outlet_flows_m3_s=outlets,signed_outward_boundary_flows_m3_s=flows,epsilon_Q=abs(incoming-self.Q)/self.Q,epsilon_mass=abs(total-incoming)/self.Q,outlet_fractions={n:v/total for n,v in outlets.items()},area_average_pressure_pa=pressure,negative_outward_flux_area_fraction=backflow,velocity_L2=self.velocity_l2(u),velocity_max_m_s=float(np.linalg.norm(u,axis=1).max()),velocity_max_location_m=self.points[np.linalg.norm(u,axis=1).argmax()].tolist(),pressure_range_pa=[float(p.min()),float(p.max())],velocity_finite=bool(np.isfinite(u).all()),pressure_finite=bool(np.isfinite(p).all()),wall_velocity_max_m_s=wall,clot_velocity_max_m_s=clotmax,wall_noslip_pass=wall<=5e-13)
