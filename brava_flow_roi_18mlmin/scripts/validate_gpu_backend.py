from pathlib import Path
import sys,json,hashlib
import numpy as np,pyvista as pv
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'vendor'))
from vascular_validation.postprocess import SolutionMeasurements
from flow_solver_support.flow_parser import parse_solver_log
from flow_solver_support.solver_checks import linear_gate,nonlinear_gate
from flow_solver_support.wss_case import recover
cpu=ROOT/'cases/diagnostic_cpu8_dt0p0002';gpu=ROOT/'cases/diagnostic_gpu8_host_sync'
records={}
for c in [cpu,gpu]:
 e=json.loads((c/'reports/execution.json').read_text());assert e['exit_code']==0 and e['final_step']==2 and not e['health_failures']
 h=parse_solver_log((c/'run/solver.log').read_text(),.0002);linear_gate(dict(exit_code=0,history=h));nonlinear_gate(h)
 records[c.name]={'exit_code':e['exit_code'],'native_steps':e['final_step'],'linear_solves':len(h['linear_solves']),'solver_sha256':e['solver_sha256']}
a=np.load(ROOT/'mesh/SV_MESH/mesh_arrays.npz');measure=SolutionMeasurements(ROOT/'mesh/SV_MESH/mesh_arrays.npz',3e-7,.2696);rows=[]
for step in [1,2]:
 paths=[c/'run/8-procs'/f'result_{step:03d}.vtu' for c in [cpu,gpu]]
 uc,pc=measure.read(paths[0]);ug,pg=measure.read(paths[1])
 wc,rc=recover(a['points_m'],a['tetra'],a['boundary_triangles'],a['facet_tags'],uc,pc,.00345312)
 wg,rg=recover(a['points_m'],a['tetra'],a['boundary_triangles'],a['facet_tags'],ug,pg,.00345312)
 errors={'velocity_relative_L2':measure.velocity_l2(ug-uc)/measure.velocity_l2(uc),
 'pressure_relative_L2':measure.velocity_l2((pg-pc)[:,None])/measure.velocity_l2(pc[:,None]),
 'WSS_area_relative_L2':float(np.sqrt(np.sum(rc['area']*(rg['magnitude']-rc['magnitude'])**2)/np.sum(rc['area']*rc['magnitude']**2)))}
 assert max(errors.values())<1e-8,errors
 rows.append(dict(step=step,errors=errors,CPU=measure.measure(uc,pc),GPU=measure.measure(ug,pg),sha256=[hashlib.sha256(p.read_bytes()).hexdigest() for p in paths]))
out=dict(passed=True,method='Identical real mesh, material, boundary conditions, dt0.0002s and two native steps; CPU8 original executable vs GPU8 isolated host-synchronization fix',criterion_relative_L2=1e-8,records=records,comparisons=rows,scope='Backend equivalence only; two transient steps are not steady CFD validation')
(ROOT/'gpu_solver_fix/backend_equivalence.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='comparisons'},indent=2));print([r['errors'] for r in rows])
