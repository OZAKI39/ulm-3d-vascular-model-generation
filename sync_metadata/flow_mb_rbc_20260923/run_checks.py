"""Run the migration checks against this checkout, without a full production solve."""
from pathlib import Path
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
FLOW = 'formal_3D_flow_solver/FEM_SimVascular/tests/flow_2mmps/'
selectors = [
    'particle_3d/tests/particle0',
    'particle_3d/tests/rbc_mb_flow_rotation',
    'particle_3d/tests/rbc_mb_coflow',
    'particle_3d/tests/rbc_mb_normal',
    'particle_3d/tests/rbc_mb_interaction',
    FLOW + 'test_physics_guards.py',
    FLOW + 'test_field_diagnostics.py',
    FLOW + 'test_residual_semantics.py',
    FLOW + 'test_deliverables.py::test_native_solve_and_steady_evidence',
    FLOW + 'test_deliverables.py::test_new_export_is_actual_native_solution_and_original_mesh',
    FLOW + 'test_deliverables.py::test_old_frozen_artifacts_unchanged',
]
env = dict(os.environ)
env.update(PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           SONOVUE_ROOT=str(ROOT / 'sonovue_size_distribution_v0'),
           PYTHONPATH=os.pathsep.join(str(ROOT / p) for p in
                                     ['particle_3d/src', 'formal_3D_flow_solver/FEM_SimVascular/src']))
command = [sys.executable, '-m', 'pytest', '-q', '-p', 'no:cacheprovider']
subprocess.run(command + selectors, cwd=ROOT, env=env, check=True)
subprocess.run(command + [
    'particle_3d/shear_lift_audit/tests',
    '--ignore=particle_3d/shear_lift_audit/tests/test_no_solver_modification.py',
    '--ignore=particle_3d/shear_lift_audit/tests/test_saved_trajectory_readonly.py',
], cwd=ROOT, env=env, check=True)
# Those three historical read-only checks contain original-machine absolute paths.
# verify_snapshot.py covers the retained source/trajectory hashes in this checkout.
