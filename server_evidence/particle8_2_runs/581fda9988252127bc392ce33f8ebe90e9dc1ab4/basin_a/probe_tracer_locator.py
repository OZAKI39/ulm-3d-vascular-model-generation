from pathlib import Path
import sys,json,time
import numpy as np
sys.path.insert(0,str(Path.cwd()/'particle_3d/src'))
from particle_3d.particle82_tracers import init_environment,trace_positions
import vtk
base=Path.cwd().parent
positions=np.load(base/'point_smoke/birth_positions.npz')['positions']
env=init_environment()
original=vtk.vtkStreamTracer
class Factory:
 LENGTH_UNIT=original.LENGTH_UNIT
 def __new__(cls):
  obj=original();obj.SetInterpolatorTypeToCellLocator();obj.SetCellLocatorToStaticCellLocator();return obj
vtk.vtkStreamTracer=Factory
results=[]
for step,error in [(2e-7,1e-11),(5e-8,1e-12),(2.5e-8,1e-12)]:
 start=time.time();rows=trace_positions(positions,step_m=step,error=error)
 counts={o:sum(r['outlet']==o for r in rows) for o in ['OUTLET_01','OUTLET_02','OUTLET_03']}
 bad=[{'id':r['seed_index'],'reason':r['end_reason'],'vtk':r['vtk_reason'],'n':len(r['path'])} for r in rows if not r['outlet']]
 result=dict(step_m=step,error=error,counts=counts,bad=bad,wall_s=time.time()-start)
 results.append(result);print(json.dumps(result),flush=True)
(base/'TRACER_LOCATOR_PROBE.json').write_text(json.dumps(results,indent=2))
