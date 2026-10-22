"""Executed only by SimVascular embedded Python 3.5; no mesh generation."""
import json
import sys
import sv
import vtk

print(json.dumps(dict(status='PASS', python=sys.version, vtk=vtk.vtkVersion.GetVTKVersion(),
                     sv_module=sv.__file__, tetgen_class=str(sv.meshing.TetGen),
                     new_meshes=0, CFD_runs=0)))
