from launch import *
import numpy as np,pyvista as pv
run('p1_baseline','initial')
cwd=HERE/'runtime/p1_baseline/run';data=np.vstack([np.loadtxt(p,delimiter=',') for p in cwd.glob('initial_rank_*.csv')]);v=pv.read(P1/'run/8-procs/result_010.vtu');ids=data[:,0].astype(int)
# Native global IDs are zero-based; independently compare positions and both saved fields.
check=dict(position_max_abs_m=float(np.max(abs(data[:,1:4]-v.points[ids]))),velocity_max_abs_m_s=float(np.max(abs(data[:,4:7]-v['Velocity'][ids]))),pressure_max_abs_Pa=float(np.max(abs(data[:,7]-v['Pressure'][ids]))))
mesh=pv.read(P1/'SV_MESH/mesh-complete.mesh.vtu')
check['position_original_mesh_max_abs_m']=float(np.max(abs(data[:,1:4]-mesh.points[ids])))
check['output_coordinate_dtype']=str(v.points.dtype)
check['output_coordinates_equal_at_serialized_precision']=bool(np.array_equal(data[:,1:4].astype(v.points.dtype),v.points[ids]))
check['PASS']=check['position_original_mesh_max_abs_m']==0 and check['output_coordinates_equal_at_serialized_precision'] and check['velocity_max_abs_m_s']<1e-12 and check['pressure_max_abs_Pa']<1e-9
check['checkpoint_sha256']=sha(P1/'run/8-procs/stFile_010.bin');check['VTU_sha256']=sha(P1/'run/8-procs/result_010.vtu');put('P1_RESTART_CORRESPONDENCE.json',check);assert check['PASS'],check
run('p1_baseline','parity')
