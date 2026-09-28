from pathlib import Path
import ast,json,sys,hashlib,shutil,xml.etree.ElementTree as ET
import numpy as np
import pyvista as pv
V=Path(__file__).resolve().parents[1];C=V.parent
sys.path.insert(0,str(C/'solver_support/src'))
from flow_solver_support import wss as p
from flow_solver_support.wss_case import material,coordinate_identity
A=V/'stage1';(A/'data').mkdir(exist_ok=True)
MU=.00345312
# Reuse the previous audit's discriminating checks, changing only imported production module/output.
old=(C/'wss_audit/scripts/analytical_validation.py').read_text()
body=old.split('def main():\n',1)[1].split('    R=4e-6;',1)[0]
exec('def run_tensor_checks():\n'+body+'\nrun_tensor_checks()')
baseline=pv.read(C/'rotate_visualization/input_data/field_diagnostics/data/wall_wss_si.vtp')
new=pv.read(A/'recomputed/data/wall_wss_si.vtp')
assert np.array_equal(baseline.points,new.points) and np.array_equal(baseline.faces,new.faces)
checks={}
for name in baseline.array_names:
 if name in new.array_names:
  checks[name]=float(np.max(np.abs(baseline[name]-new[name])))
assert checks['WSS_raw_Pa']<=1e-9 and checks['WSS_display_Pa']<=1e-9
src=ast.parse((V/'evidence/before/compute_field_diagnostics.py').read_text());dst=ast.parse(Path(p.__file__).read_text());identity={}
for node in dst.body:
 if isinstance(node,ast.FunctionDef):
  original=next(n for n in src.body if isinstance(n,ast.FunctionDef) and n.name==node.name)
  if node.name=='surface':node.body=node.body[1:]
  identity[node.name]=ast.dump(node,include_attributes=False)==ast.dump(original,include_attributes=False)
assert all(identity.values())
case=Path('/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1')
props=material(case)
example=np.array([[1.23456789e-6,2.34567891e-6,3.45678912e-6]])
serialization=coordinate_identity(example.astype(np.float32),example)
try:coordinate_identity(example+1e-9,example)
except ValueError:serialization['arbitrary_coordinate_change_rejected']=True
else:raise AssertionError('arbitrary geometry mismatch accepted')
negative=A/'material_mismatch';negative.mkdir(exist_ok=True);(negative/'run').mkdir(exist_ok=True)
shutil.copy2(case/'policy.json',negative/'policy.json')
tree=ET.parse(case/'run/solver.xml');tree.find('.//Viscosity/Value').text=str(MU*2);tree.write(negative/'run/solver.xml')
try:material(negative)
except ValueError as e:rejected=str(e)
else:raise AssertionError('mismatched material was not rejected')
result=dict(status='PASS',baseline_max_abs_difference=checks,regression_tolerance_Pa=1e-9,
 core_function_AST_unchanged=identity,source_recovered_sha256=hashlib.sha256((V/'evidence/before/compute_field_diagnostics.py').read_bytes()).hexdigest(),
 material=props,mismatched_case_material_rejected=rejected,coordinate_serialization_test=serialization,production_entry='rotate_visualization/prepare_surface_data.py')
(A/'regression.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
