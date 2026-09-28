"""New geometry adapter; existing finite-bubble physics stays byte-identical."""
from pathlib import Path
from types import SimpleNamespace
import hashlib,json,sys
import numpy as np
import pyvista as pv
ROOT=Path(__file__).resolve().parents[1]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def environment(source_root,flow_path):
 source_root=Path(source_root);sys.path.insert(0,str(source_root/'particle_3d/src'))
 from particle_3d.field import FrozenFEMField
 from particle_3d.wall_geometry import WallGeometry
 from particle_3d.validation_boundary import ValidationBoundaryClassifier
 from particle_3d.inlet_flux import frozen_boundary_flux
 from particle_3d.injection_method_c import TruncatedSonoVue
 from particle_3d.injection_admission import FiniteSizeAdmission
 from particle_3d.open_boundary_rim import build_rim_topology
 spec=json.loads((ROOT/'reports/ACTIVE_FLOW.json').read_text());case=ROOT/spec['case']
 assert Path(flow_path).resolve()==(case/'frozen_flow/steady_flow.vtu').resolve()
 assert sha(flow_path)==spec['flow_sha256']
 for relative,expected in json.loads((case/'input_hashes.json').read_text()).items():
  assert sha(case/relative)==expected, 'Changed CFD input: '+relative
 mesh=pv.read(ROOT/'mesh/SV_MESH/mesh-complete.mesh.vtu');flow=pv.read(flow_path)
 field=FrozenFEMField.from_grids(mesh,flow)
 boundaries={r:pv.read(ROOT/'mesh/SV_MESH/mesh-surfaces'/f'{r}.vtp') for r in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03','WALL']}
 w=boundaries['WALL'];faces=w.faces.reshape(-1,4)[:,1:]
 wall=WallGeometry(np.asarray(w.points)[faces],provenance={'role':'REAL_BRAVA_FEM_WALL','source':str(ROOT/'mesh/SV_MESH/mesh-surfaces/WALL.vtp'),'physical_coating_thickness_m':0.,'wall_rigid':True,'wall_stationary':True})
 wall.global_node_ids=np.asarray(w.point_data['GlobalNodeID'])[faces]-1
 caps={r:np.asarray(s.point_data['GlobalNodeID'])[s.faces.reshape(-1,4)[:,1:]]-1 for r,s in boundaries.items() if r!='WALL'}
 hashes={r:sha(ROOT/'mesh/SV_MESH/mesh-surfaces'/f'{r}.vtp') for r in boundaries}
 wall.open_boundary_topology=build_rim_topology(wall.global_node_ids,caps,hashes)
 classifier=ValidationBoundaryClassifier({r:s for r,s in boundaries.items() if r.startswith('OUTLET_')})
 audit,samplers=frozen_boundary_flux(mesh,flow,boundaries)
 distribution=TruncatedSonoVue(source_root/'sonovue_size_distribution_v0')
 env=SimpleNamespace(mesh=mesh,field=field,wall=wall,boundaries=boundaries,classifier=classifier,
  mu=spec['mu_Pa_s'],viscosity={'mu_Pa_s':spec['mu_Pa_s'],'source':'active CFD XML'},
  provenance=spec,audit=audit,flux_audit=audit,sampler=samplers['INLET'],distribution=distribution,
  repo_root=source_root,new_flow_path=Path(flow_path),new_flow_sha256=sha(flow_path))
 env.checker=FiniteSizeAdmission(None,None,wall=wall,field=field)
 return env
