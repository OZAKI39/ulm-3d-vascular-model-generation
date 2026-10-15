from common import *
from pd_clot.lammps_fragmentation.y_bifurcation_sources import YSources
from pd_clot.portable_paths import relocate
s=YSources();register=R/'verification/y_bifurcation_clot_placement_visual_integration_v1/SOURCE_ASSET_REGISTER.json';assets=[]
for name,info in s.files.items():
 p=relocate(info['local_readonly_copy']);assert sha(p)==info['sha256'];assets.append(dict(relative_name=name,path=str(p),bytes=p.stat().st_size,sha256=sha(p)))
dump(O/'REUSED_CFD_ASSETS.json',dict(status='PASS',source_register=str(register),register_sha256=sha(register),files=assets,geometry_id=s.geometry_id,mesh_id=s.mesh_id,flow_id=s.flow_id,CFD_recomputed=False,used='Historical unobstructed Y background and native WALL surfaces, no new clot-CFD nodal traction map'))
# Preserve small imported source modules, with original absolute path and digest.
folder=O/'source_identity/reused';folder.mkdir(exist_ok=True)
for rel in ['pd_clot/lammps_fragmentation/native_y_wall_query.py','pd_clot/lammps_fragmentation/native_y_wall_contact.py','pd_clot/lammps_fragmentation/background_flow_sampler.py','pd_clot/lammps_fragmentation/y_bifurcation_sources.py','pd_clot/lammps_fragmentation/y_mural_clot_geometry.py','pd_clot/lammps_fragmentation/early_flow_coupling.py','pd_clot/lammps_fragmentation/detached_fragment_metrics.py','pd_clot/portable_paths.py']:
 p=R/rel;(folder/p.name).write_bytes(p.read_bytes())
print('DEPENDENCIES_VERIFIED',len(assets))
