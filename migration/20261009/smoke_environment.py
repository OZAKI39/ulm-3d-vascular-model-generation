"""Real native imports, immutable mesh, Boolean/collision and off-screen frame."""
from pathlib import Path
import sys,json,importlib,importlib.metadata,hashlib
root=Path('/workspace/vascular_printing');sys.path.insert(0,str(root));out=root/'migration/20261009';checks={}
checks['python']=dict(version=sys.version,prefix=sys.prefix,base_prefix=sys.base_prefix,isolated=sys.prefix!=sys.base_prefix and Path(sys.prefix)==root/'.venv')
assert checks['python']['isolated']
for name in ['numpy','scipy','vtk','pyvista','trimesh','manifold3d','fcl','rtree','geomdl','open3d','cadquery','nibabel','skimage','h5py','yaml','ot','embreex']:
 m=importlib.import_module(name);checks[name]=dict(imported=True,version=getattr(m,'__version__',None))
import trimesh,numpy as np,pyvista as pv
from vascular_processing.core_support_inputs import discover_source,branches
from vascular_processing.core_support_p1s import default_output,profiles,evaluated
source=discover_source();mesh=trimesh.load_mesh(source['core']);assert mesh.is_volume
bs,top=branches(source,.21)
checks['frozen_core']=dict(path=str(source['core']),sha256=hashlib.sha256(source['core'].read_bytes()).hexdigest(),watertight=mesh.is_watertight,faces=len(mesh.faces),branches=len(bs),bounds_mm=mesh.bounds.tolist())
assert checks['frozen_core']['sha256']=='84f1989392fb8d210432750a8022a830943b9bff890b044aa9d310822b619ce0'
a=trimesh.creation.box();b=a.copy();b.apply_translation([.5,0,0]);c=trimesh.boolean.intersection([a,b],engine='manifold');assert abs(c.volume-.5)<1e-6
manager=trimesh.collision.CollisionManager();manager.add_object('box',a);assert manager.in_collision_single(b)
checks['boolean_and_collision']=True
archived=default_output();s=json.loads((archived/'regeneration_summary.json').read_text());bundle=profiles(archived/'bambu/profiles',review_only=True)
checks['archived_P1S_evaluation']=evaluated(archived/'candidates'/s['selected_candidate'],bundle,readonly=True)['slice']['slicer_success']
checks['manufacturing_status_unchanged']=s['status']
plot=pv.Plotter(off_screen=True,window_size=(1000,750));plot.set_background('white');plot.add_mesh(pv.wrap(mesh),color='tomato');plot.view_isometric();plot.reset_camera();pixels=plot.screenshot(out/'frozen_core_offscreen.png');plot.close();assert pixels.std()>5
checks['offscreen']=dict(shape=list(pixels.shape),std=float(pixels.std()),path=str(out/'frozen_core_offscreen.png'))
(out/'environment_smoke.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
