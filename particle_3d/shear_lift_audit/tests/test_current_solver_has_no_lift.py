import ast
from pathlib import Path
import numpy as np
from particle_3d.particle65_motion import assemble_v1
from particle_3d.particle_shapes import Sphere
from particle_3d.hydrodynamic_resistance import sphere_self_diagonal

def test_current_rhs_is_only_stokes_background():
    shape=Sphere(np.array([0.,0.,0.]),1e-6);u=np.array([1e-3,2e-3,0,0,0,10.]);mu=.00345312
    s=assemble_v1({1:shape},{1:u},mu)
    np.testing.assert_array_equal(s.rhs,sphere_self_diagonal(shape,mu)*u)
    np.testing.assert_array_equal(s.matrix.toarray(),np.diag(sphere_self_diagonal(shape,mu)))

def test_audit_has_no_formal_import_path_and_no_integrator():
    src=Path(__file__).resolve().parents[2]/'src/particle_3d'
    assert src.is_dir()
    for p in src.glob('*.py'):
        tree=ast.parse(p.read_text())
        for n in ast.walk(tree):
            if isinstance(n,ast.ImportFrom):assert 'shear_lift_audit' not in (n.module or '')
    audit=Path(__file__).resolve().parents[1]/'code/run_audit.py'
    names={n.func.attr for n in ast.walk(ast.parse(audit.read_text())) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)}
    assert not names.intersection({'step_to','integrate_one','v1_trial','advance_orientation'})
