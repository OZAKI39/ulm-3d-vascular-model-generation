"""Run migrated mesh audit against actual data, exporting to an isolated folder."""
from pathlib import Path
import importlib.util
import json
import shutil
import sys
import numpy as np

C = Path('/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular')
R = Path(__file__).resolve().parent
W = R / 'isolated_mesh_check'
W.mkdir(exist_ok=False)
(W / 'inputs').symlink_to(C / 'inputs', target_is_directory=True)
(W / 'configs').symlink_to(C / 'configs', target_is_directory=True)
(W / 'outputs/sv1').mkdir(parents=True)
(W / 'outputs/sv1/mesh_generation').symlink_to(C / 'outputs/sv1/mesh_generation', target_is_directory=True)
(W / 'reports/sv1').mkdir(parents=True)
shutil.copy2(C / 'reports/sv1/mesh_policy_provenance.json', W / 'reports/sv1/mesh_policy_provenance.json')
spec = importlib.util.spec_from_file_location('migrated_mesh_audit', C / 'mesh_generate/scripts/audit_mesh.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
assert module.ROOT == C
module.ROOT = W  # Only this test process redirects audit outputs away from production.
sys.argv = [str(C / 'mesh_generate/scripts/audit_mesh.py'), '--attempt', 'primary']
module.main()
with np.load(C / 'outputs/sv1/SV_MESH/mesh_arrays.npz') as old, np.load(W / 'outputs/sv1/SV_MESH/mesh_arrays.npz') as new:
    assert set(old.files) == set(new.files)
    for key in old.files:
        assert np.array_equal(old[key], new[key]), key
    result = dict(all_pass=True, arrays_identical=old.files,
                  vertices=len(new['points_m']), tetrahedra=len(new['tetra']),
                  audit_output=str(W), production_mesh_overwritten=False)
(R / 'mesh_verification.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
