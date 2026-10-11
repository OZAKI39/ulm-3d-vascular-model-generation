"""Audit the already-generated mesh, restore frozen PLC doubles, then first CFD.
No additional volume meshing and no previous CFD solve to restart.
"""
from analysis.dual_track_v1.common import *
from analysis.dual_track_v1.mesh_b import audit
from analysis.dual_track_v1.cfd_b import main as cfd
import time,traceback,pyvista as pv
from scipy.spatial import cKDTree
out=B/'baseline';start=time.monotonic()
try:
 assert not (out/'FORMAL_ATTEMPT.json').exists();(out/'MESH_QUALITY_AUDIT.json').rename(out/'MESH_QUALITY_NATIVE_FLOAT32.json');(out/'SV_MESH').rename(out/'SV_MESH_NATIVE_FLOAT32')
 for name in ['baseline_STOP.json','baseline_EXECUTION.json']:
  (B/name).rename(B/name.replace('baseline_','baseline_native_float32_'))
 d=json.loads((B/'DESIGN_FREEZE.json').read_text());geom=json.loads((B/'GEOMETRY_AUDIT.json').read_text());r=audit(out/'raw',out/'SV_MESH',Path('/workspace/y_clot/deps/y_cfd_vendor'));grid=pv.read(out/'SV_MESH/mesh-complete.mesh.vtu');surface=pv.read(out/'SV_MESH/mesh-complete.exterior.vtp');plc=pv.read(B/'clot_present_lumen.vtp');dist,_=cKDTree(plc.points).query(surface.points);relative=abs(grid.volume/geom['fluid_volume_m3']-1)
 checks=dict(positive_quality=r['quality']['q_min']>=d['mesh']['quality']['min_SICN'],volume=relative<=d['mesh']['quality']['relative_volume_to_PLC'],points=r['points']<=d['mesh']['budgets']['max_points'],elements=r['tetrahedra']<=d['mesh']['budgets']['max_tetrahedra'],closed=surface.n_open_edges==0,unchanged_boundary=dist.max()<1e-10)
 r.update(checks=checks,status='PASS' if all(checks.values()) else 'FAIL',volume_m3=float(grid.volume),relative_volume_PLC=relative,PLC_sha256=d['geometry_sha256'],mesh_sha256=sha(out/'SV_MESH/mesh-complete.mesh.vtu'),max_boundary_node_distance_from_PLC_m=float(dist.max()),local_target_m=.0001,postprocess_wall_s=time.monotonic()-start);save(out/'MESH_QUALITY_AUDIT.json',r);save(B/'MESH_QUALITY_AUDIT.json',dict(level='baseline',**r));assert all(checks.values())
 save(out/'FORMAL_ATTEMPT.json',dict(formal_CFD_attempt=1,level='baseline',started_unix=time.time()));cfd('baseline')
except BaseException as e:save(B/'baseline_STOP.json',dict(stage='precision_restore_or_CFD',error=repr(e),traceback=traceback.format_exc(),wall_s=time.monotonic()-start));raise
finally:save(B/'baseline_EXECUTION.json',dict(wall_s=time.monotonic()-start,original_mesh_execution=str(out/'MESH_EXECUTION.json'),no_remesh=True))
