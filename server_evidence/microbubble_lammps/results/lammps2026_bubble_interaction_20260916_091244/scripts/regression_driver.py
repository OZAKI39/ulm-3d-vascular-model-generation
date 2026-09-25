from pathlib import Path
import sys,json,importlib.util,hashlib
import numpy as np
R=Path(__file__).resolve().parents[1];W=Path('/workspace/lammps_migration/new_20260916_091244');stage=sys.argv[1];root=R/'migration'/stage
spec=importlib.util.spec_from_file_location('baseline_runner',root/'scripts/run_validation.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
m.ROOT=root;m.WORK=W/stage;m.ENGINE=W
with (root/'MIGRATION_VALIDATION_STARTED.json').open('x') as f:json.dump({'no_retry':True,'actual_new_root':str(W)},f)
assert m.integrity(root)
if stage=='coupling':
 out=root/'validation/cpp';out.mkdir()
 jobs=[('uniform','UNIFORM_FIELD.h5'),('linear','LINEAR_FIELD.h5'),('real_nodes','FROZEN_FLOW_FIELD_V0.h5'),('real_interior','FROZEN_FLOW_FIELD_V0.h5'),('invalid_domain','UNIFORM_FIELD.h5'),('invalid_missing','MISSING_CORNER_FIELD.h5'),('invalid_solid','FROZEN_FLOW_FIELD_V0.h5'),('force','UNIFORM_FIELD.h5')]
 for label,field in jobs:m.command([str(m.WORK/'build/flow_audit_cli'),'force' if label=='force' else 'query',str(root/'fields'/field),str(root/'validation/queries'/(label+'.csv')),str(out/(label+'.csv'))])
 r=m.unit_audit(root);m.save(root/'validation/UNIT_AUDIT.json',r);assert r['status']=='PASS',r
 for sp in json.loads((root/'contracts/PALABOS_LAMMPS_COUPLING_V0_CONTRACT.json').read_text())['cases']:m.run_case(sp)
else:
 m.cmd([str(m.WORK/'build/wall_audit'),str(root/'provenance/closed_geometry_m.stl'),str(root/'validation/WALL_QUERY_REFERENCE.csv'),str(root/'validation/WALL_QUERY_CPP.csv')])
 a=m.table(root/'validation/WALL_QUERY_REFERENCE.csv');b=m.table(root/'validation/WALL_QUERY_CPP.csv');err=float(np.max(np.abs(a['wall_distance_m']-b['distance_m'])));assert err<=1e-12;m.save(root/'validation/WALL_DISTANCE_AUDIT.json',{'status':'PASS','maximum_error_m':err})
 lookup={s['name']:s for s in json.loads((root/'contracts/CASE_INDEX.json').read_text())}
 names=[f'case{i}_{s}' for i in range(3) for s in ['c050','c025','c0125']]+['case3_'+q+'_c025' for q in ['d10','d50','d90']]+['case4_mpi1_c025','case4_mpi4_c025','case5_c025','kokkos_c025']
 for name in names:m.run_case(lookup[name])
 from finalize_passive_transport_v0 import trajectory,xyz,vel
 old=Path('/workspace/microbubble_lammps/results/passive_transport_v0_20260915_233516');a=trajectory(old/'cases/case2_c025');b=trajectory(root/'cases/case2_c025');assert np.array_equal(a['time_s'],b['time_s']) and np.array_equal(a['particle_id'],b['particle_id'])
 dif=np.linalg.norm(xyz(a)-xyz(b),axis=1);r={'final_position_difference_m':float(dif[-1]),'max_trajectory_difference_m':float(dif.max()),'max_velocity_difference_m_s':float(np.linalg.norm(vel(a)-vel(b),axis=1).max()),'same_record_times':True,'status':'PASS' if dif.max()<=1e-10 else 'FAIL'};m.save(root/'NEW_VS_OLD_CASE2.json',r);assert r['status']=='PASS'
r=m.audit_all(root);m.save(root/'REMOTE_INDEPENDENT_FINALIZER.json',r);assert r['status']=='PASS',r
print(json.dumps({'stage':stage,'status':'PASS'}),flush=True)
