import ast
import inspect
from dataclasses import fields
from particle_3d.microbubble import MicrobubbleState
from particle_3d.integrator import advance_single_microbubble


def test_one_sphere_state_and_caller_required_dt():
    assert {f.name for f in fields(MicrobubbleState)} == {'particle_id','position_m','radius_m','velocity_m_s','angular_velocity_s_inv'}
    assert inspect.signature(advance_single_microbubble).parameters['dt_s'].default is inspect.Parameter.empty


def test_runtime_modules_do_not_depend_on_future_physics_or_solver(repo):
    for name in ['microbubble.py','integrator.py','validation_boundary.py','sonovue_adapter.py']:
        tree = ast.parse((repo/'particle_3d/src/particle_3d'/name).read_text())
        imports = [n.module or '' for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)]
        imports += [a.name for n in ast.walk(tree) if isinstance(n,ast.Import) for a in n.names]
        assert not any(any(word in name.lower() for word in ['lammps','red_blood','petsc','solver','collision','adhesion']) for name in imports)
