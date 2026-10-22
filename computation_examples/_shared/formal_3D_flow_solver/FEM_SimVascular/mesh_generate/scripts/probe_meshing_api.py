# Executed by the official distribution's embedded Python (3.5-compatible syntax).
import json
import os
import platform
import sys
import sv
import vtk
result={'status':'PASS','python':sys.version,'architecture':platform.machine(),
        'sv_module_file':getattr(sv,'__file__',None),'sv_version':str(getattr(sv,'__version__','not exported')),
        'vtk_version':vtk.vtkVersion.GetVTKVersion(),'actual_mesher_type':str(type(sv.meshing.TetGen())),
        'model_type':str(type(sv.modeling.PolyData())),
        'meshing_api_present':all(hasattr(sv.meshing.TetGen(),n) for n in ['load_model','set_walls','generate_mesh','write_mesh']),
        'geometry_loaded':False,'mesh_generated':False}
options=sv.meshing.TetGenOptions(global_edge_size=1.0,surface_mesh_flag=True,volume_mesh_flag=True)
options.use_mmg=False
result['options_surface_mesh_flag']=options.surface_mesh_flag
result['options_volume_mesh_flag']=options.volume_mesh_flag
result['options_use_mmg']=options.use_mmg
assert result['meshing_api_present']
print('MESH_API_PROBE_JSON='+json.dumps(result,sort_keys=True))
