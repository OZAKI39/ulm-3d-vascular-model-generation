"""Run actual surface preparation/import in isolation after retiring the smoke report."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys

import numpy as np
import pyvista as pv

C=Path('/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular')
R=Path(__file__).resolve().parent
W=R/'isolated_surface_check'
W.mkdir(exist_ok=False)
(W/'mesh_generate/scripts').mkdir(parents=True)
shutil.copytree(C/'mesh_generate/src',W/'mesh_generate/src',ignore=shutil.ignore_patterns('__pycache__'))
for name in ['prepare_surface.py','import_surface_model.py']:
    shutil.copy2(C/'mesh_generate/scripts'/name,W/'mesh_generate/scripts'/name)
(W/'inputs').symlink_to(C/'inputs',target_is_directory=True)
(W/'configs').mkdir()
(W/'reports/mesh_and_flow').mkdir(parents=True)
api_path=W/'reports/mesh_and_flow/simvascular_api.json'
original_api=json.loads((C/'reports/mesh_and_flow/simvascular_api.json').read_text())
api_path.write_text(json.dumps(original_api))
distribution=json.loads((C/'reports/mesh_and_flow/simvascular_distribution.json').read_text())
launcher=C/distribution['path']/'simvascular'
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONNOUSERSITE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',
         QT_QPA_PLATFORM='offscreen',PYVISTA_OFF_SCREEN='true',XDG_CONFIG_HOME=str(W/'app_config'),XDG_CACHE_HOME=str(W/'app_cache'))
env.pop('PYTHONPATH',None)
commands={
    'prepare_surface':[sys.executable,str(W/'mesh_generate/scripts/prepare_surface.py')],
    'import_surface_model':[str(launcher),'-python','--',str(W/'mesh_generate/scripts/import_surface_model.py')],
}

def run(name,command,success=True):
    result=subprocess.run(command,cwd=W,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=90)
    (R/(name+'.log')).write_text(result.stdout)
    if success:assert result.returncode==0,(name,result.stdout[-1800:])
    else:assert result.returncode!=0 and 'Official meshing API check must pass' in result.stdout,(name,result.stdout[-1800:])

assert not (W/'reports/mesh_and_flow/official_smoke.json').exists()
run('prepare_surface',commands['prepare_surface'])
run('import_surface_model',commands['import_surface_model'])
for filename in ['source.vtp','imported.vtp']:
    current=pv.read(C/'outputs/mesh_and_flow/model'/filename)
    actual=pv.read(W/'outputs/mesh_and_flow/model'/filename)
    assert np.array_equal(current.points,actual.points),filename
    assert np.array_equal(current.faces,actual.faces),filename
    for key in current.cell_data:assert np.array_equal(current.cell_data[key],actual.cell_data[key]),(filename,key)
assert json.loads((W/'configs/face_map.json').read_text())==json.loads((C/'configs/face_map.json').read_text())
assert json.loads((W/'reports/mesh_and_flow/geometry_import.json').read_text())==json.loads((C/'reports/mesh_and_flow/geometry_import.json').read_text())

# Reject bad prerequisites through the real scripts, before any model is rewritten.
model_dir=W/'outputs/mesh_and_flow/model'
before={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in model_dir.iterdir()}
for failure in ['failed_status','missing_meshing_api']:
    api=json.loads(json.dumps(original_api))
    if failure=='failed_status':api['status']='FAIL'
    else:api['result']['meshing_api_present']=False
    api_path.write_text(json.dumps(api))
    for name,command in commands.items():run(name+'_'+failure,command,success=False)
    assert {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in model_dir.iterdir()}==before

result=dict(all_pass=True,old_smoke_report_required=False,actual_surface_preparation_pass=True,
            official_surface_import_pass=True,source_and_imported_geometry_arrays_identical=True,
            face_map_identical=True,geometry_import_report_identical=True,
            invalid_API_status_rejected=True,unavailable_meshing_API_rejected=True,
            production_data_overwritten=False)
(R/'surface_entrypoints_verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
