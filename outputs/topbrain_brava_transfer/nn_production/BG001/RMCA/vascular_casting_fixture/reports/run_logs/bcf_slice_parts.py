from pathlib import Path
import sys,json
import numpy as np,trimesh
from scipy.spatial.transform import Rotation
sys.path.insert(0,'/home/lzy/projects/vascular_printing')
from vascular_processing.casting_fixture_geometry import port_matrix
from vascular_processing.casting_fixture_pipeline import default_output
from vascular_processing.casting_fixture_profiles import slice_part
p=default_output();engine=json.loads((p/'config/profiles/bambu_engine.json').read_text());profiles={k:json.loads((p/'bambu/frame/profiles'/f'{k}.json').read_text()) for k in ['machine','process','filament']};layout=json.loads((p/'assembly/layout.json').read_text());results={}
for port in layout['ports']:
 for kind in ['closure','locking_pin']:
  ident=port['port_id']+'_'+kind;mesh=trimesh.load_mesh(p/'closures'/f'port_{ident}.stl');mesh.apply_transform(np.linalg.inv(port_matrix(port)))
  specs=[('y',90),('y',-90),('x',-90),('x',90),('x',0),('x',180)] if kind=='closure' and port['role']=='MASTER_DATUM' else [('x',-90),('x',90),('y',90),('y',-90),('x',0),('x',180)]
  for axis,angle in specs:
   r=slice_part(mesh,engine,profiles,p/'bambu/orientation_checks'/f'{ident}_{axis}_{angle}',Rotation.from_euler(axis,angle,degrees=True).as_matrix());print(ident,axis,angle,r['slicer_success'],r['support_free_slice_passed'],r['warning_audit'],flush=True)
   if r['support_free_slice_passed']:results[ident]=dict(axis=axis,angle=angle);break
(p/'bambu/part_print_orientations.json').write_text(json.dumps(results,indent=2))
