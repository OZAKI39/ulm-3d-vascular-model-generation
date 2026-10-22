"""Read-only interface preflight on an explicitly labelled backend diagnostic."""
from pathlib import Path
import sys,json,hashlib
import numpy as np,pyvista as pv
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'vendor'))
sys.path.insert(0,'/workspace/particle9a5_formal_trajectories_20260925T080115Z_dt1ms/particle_3d/src')
from flow_solver_support.wss_case import coordinate_identity,material
from particle_3d.field import FrozenFEMField
meshpath=ROOT/'mesh/SV_MESH/mesh-complete.mesh.vtu'
flowpath=ROOT/'cases/diagnostic_gpu8_host_sync/run/8-procs/result_002.vtu'
mesh=pv.read(meshpath);flow=pv.read(flowpath)
coordinates=coordinate_identity(flow.points,mesh.points)
field=FrozenFEMField.from_grids(mesh,flow)
canonical=mesh.cells.reshape(-1,5)[:,1:];native=flow.cells.reshape(-1,5)[:,1:]
assert np.array_equal(np.sort(canonical,axis=1),np.sort(native,axis=1))
permutation=np.argmax(native[:,:,None]==canonical[:,None,:],axis=2)
result=dict(PASS=True,scope='Interface only: two-step transient diagnostic is NOT the active steady flow; no microbubble integration performed',
            coordinates=coordinates,node_count=len(field.points_m),tetra_count=len(field.tetra),
            local_vertex_permutations=np.unique(permutation,axis=0).tolist(),
            original_particle_pair_loader='PASS',material=material(ROOT/'cases/calibration_gpu8_production'),
            files={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [meshpath,flowpath]})
(ROOT/'reports/field_interface_preflight.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
