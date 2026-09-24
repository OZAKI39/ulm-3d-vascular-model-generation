import ast
from particle_3d.injection_population import OLD_RBC_DRIFT

def test_no_double_count(repo,contract):
 assert OLD_RBC_DRIFT==contract['OLD_RBC_DRIFT']=='DISABLED'
 for name in ['injection_population.py','injection_admission.py','particle7_lifecycle.py']:
  tree=ast.parse((repo/'particle_3d/src/particle_3d'/name).read_text())
  imports=[ast.dump(x) for x in ast.walk(tree) if isinstance(x,(ast.Import,ast.ImportFrom))]
  assert not any('red_blood_cell_transport' in x or 'hematocrit_drift' in x for x in imports)
