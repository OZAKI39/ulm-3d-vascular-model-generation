"""Independent native-output checks and existing P1 wall-traction recovery in SI."""
from pathlib import Path
import argparse,json,hashlib,csv,shutil,sys
import numpy as np
import pyvista as pv
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'vendor'))
from vascular_validation.postprocess import SolutionMeasurements
from flow_solver_support.wss_case import recover,material,coordinate_identity
from flow_solver_support.flow_parser import parse_solver_log
from flow_solver_support.solver_checks import linear_gate,nonlinear_gate

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('case',type=Path);a=ap.parse_args();case=a.case.resolve()
 e=json.loads((case/'reports/execution.json').read_text());assert e['status']=='PASS','Never use an incomplete/failed CFD result'
 policy=json.loads((case/'policy.json').read_text());mu=material(case)['mu_Pa_s']
 for path,h in json.loads((case/'input_hashes.json').read_text()).items():assert sha(case/path)==h,path
 history=parse_solver_log((case/'run/solver.log').read_text(errors='replace'),policy['dt_s'])
 linear_gate(dict(exit_code=e['exit_code'],history=history));nonlinear_gate(history)
 final=case/e['final_vtu'];mesh=pv.read(case/'SV_MESH/mesh-complete.mesh.vtu');flow=pv.read(final)
 coordinate_identity(flow.points,mesh.points)
 canonical_tet=mesh.cells.reshape(-1,5)[:,1:];native_tet=flow.cells.reshape(-1,5)[:,1:]
 assert native_tet.shape==canonical_tet.shape and np.array_equal(np.sort(native_tet,axis=1),np.sort(canonical_tet,axis=1)), 'Native tetra rows or node membership mismatch'
 correspondence=dict(point_order_identical=True,tetra_rows_identical_ignoring_local_vertex_order=True,local_order_changed_cells=int(np.count_nonzero(np.any(native_tet!=canonical_tet,axis=1))),gradients_use_canonical_connectivity=True)
 measure=SolutionMeasurements(case/'SV_MESH/mesh_arrays.npz',policy['Q_target_m3_s'],policy['Umean_m_s'])
 u,p=measure.read(final);metrics=measure.measure(u,p)
 assert max(metrics['epsilon_Q'],metrics['epsilon_mass'])<=policy['mass_limit'] and metrics['wall_noslip_pass']
 arrays=np.load(case/'SV_MESH/mesh_arrays.npz');x=arrays['points_m'];t=arrays['tetra'];b=arrays['boundary_triangles'];tags=arrays['facet_tags']
 surface,raw=recover(x,t,b,tags,u,p,mu)
 area=raw['area'];w=raw['magnitude']
 # Last saved intervals check pressure and raw wall stress as well as monitored velocity/flux.
 results=sorted((case/'run'/f"{policy.get('MPI_ranks',1)}-procs").glob('result_*.vtu'),key=lambda q:int(q.stem.rsplit('_',1)[1]))
 prev_u,prev_p=measure.read(results[-2]);prev_s,prev_raw=recover(x,t,b,tags,prev_u,prev_p,mu)
 delta=dict(velocity_relative_L2=measure.velocity_l2(u-prev_u)/measure.velocity_l2(u),
 pressure_relative_L2=measure.velocity_l2((p-prev_p)[:,None])/max(measure.velocity_l2(p[:,None]),1e-30),
 wss_area_relative_L2=float(np.sqrt(np.sum(area*(w-prev_raw['magnitude'])**2)/np.sum(area*w*w))))
 assert max(delta.values())<1e-4,delta
 def weighted_q(q):
  idx=np.argsort(w);return float(np.interp(q,np.cumsum(area[idx])/area.sum(),w[idx]))
 wss=dict(area_mean_Pa=float(np.average(w,weights=area)),min_Pa=float(w.min()),P5_Pa=weighted_q(.05),P50_Pa=weighted_q(.5),P95_Pa=weighted_q(.95),max_Pa=float(w.max()),
 maximum_abs_normal_traction_Pa=float(np.max(np.abs(np.einsum('ij,ij->i',raw['traction'],raw['normal'])))),
 raw_cells=len(w),method='mu*(G+G.T)*n followed by tangential projection; raw P1 adjacent-tetra value per triangle; area-weighted nodal average only for display')
 out=case/'postprocess';out.mkdir(exist_ok=True);surface.save(out/'wall_wss_si.vtp')
 exterior=pv.read(case/'SV_MESH/mesh-complete.exterior.vtp');ids=exterior['GlobalNodeID'].astype(int)-1
 exterior.point_data['Pressure_Pa']=p[ids];exterior.point_data['Velocity']=u[ids];exterior.save(out/'pressure_surface_si.vtp')
 frozen=case/'frozen_flow';frozen.mkdir(exist_ok=True);shutil.copy2(final,frozen/'steady_flow.vtu')
 np.savez_compressed(frozen/'flow_arrays_si.npz',points_m=x,tetra=t,boundary_triangles=b,facet_tags=tags,velocity_m_s=u,pressure_pa=p)
 rows=[dict(boundary='INLET',Q_m3_s=metrics['Q_in_m3_s'],Q_uL_min=metrics['Q_in_m3_s']*6e10,fraction_percent=100.,mean_pressure_Pa=metrics['area_average_pressure_pa']['INLET'])]
 rows.extend(dict(boundary=r,Q_m3_s=q,Q_uL_min=q*6e10,fraction_percent=metrics['outlet_fractions'][r]*100,mean_pressure_Pa=metrics['area_average_pressure_pa'][r]) for r,q in metrics['outlet_flows_m3_s'].items())
 with (out/'boundary_summary.csv').open('w') as f:
  writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
 record=dict(status='PASS',mesh_flow_correspondence=correspondence,case=case.name,mode=policy['mode'],flow_sha256=sha(frozen/'steady_flow.vtu'),mesh_sha256=sha(case/'SV_MESH/mesh-complete.mesh.vtu'),mu_Pa_s=mu,
 measurements=metrics,raw_WSS=wss,last_interval_changes=delta,final_step=e['final_step'],physical_time_s=e['final_step']*policy['dt_s'],
 max_split_deviation_percentage_points=max(abs(q-1/3)*100 for q in metrics['outlet_fractions'].values()),
 split_is_prescribed_for_calibration_only=policy['mode']=='calibration',single_mesh_not_mesh_independence=True)
 dump(out/'INDEPENDENT_CHECK.json',record)
 if policy['mode']=='calibration':
  pressures={r:metrics['area_average_pressure_pa'][r] for r in metrics['outlet_flows_m3_s']};base=min(pressures.values());pressures={r:v-base for r,v in pressures.items()}
  dump(out/'designed_outlet_pressures_Pa.json',pressures)
  dump(out/'pressure_design_note.json',dict(method='Measure area-average pressure at fully extended equal-flow calibration outlets; subtract common minimum; re-solve with natural pressure traction at every outlet',reference_shift_Pa=base,final_3D_split_not_yet_verified=True))
 else:
  assert record['max_split_deviation_percentage_points']<=.5,'Final achieved split needs a pressure correction, not relabeling as balanced'
  dump(ROOT/'reports/ACTIVE_FLOW.json',{**record,'case':str(case.relative_to(ROOT))})
 print(json.dumps(record,indent=2))
if __name__=='__main__':main()
