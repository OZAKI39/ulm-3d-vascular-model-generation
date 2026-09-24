import ast
from particle_3d.particle2_audit import check_dependencies


def test_p0_p1_accepted_dependencies_are_byte_identical(p2_repo):
    p0,p1=check_dependencies(p2_repo)
    assert p0["tests"]["passed"]==46 and p1["particle1_tests"]["passed"]==67


def test_p2_runtime_has_no_future_physics_or_solver_imports(p2_repo,contract):
    root=p2_repo/"particle_3d/src/particle_3d"
    for name in ["rbc.py","rbc_distribution.py","rbc_orientation.py","rbc_integrator.py"]:
        tree=ast.parse((root/name).read_text())
        imports=[n.module or "" for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)]
        imports += [a.name for n in ast.walk(tree) if isinstance(n,ast.Import) for a in n.names]
        assert not any(any(bad in n.lower() for bad in ["lammps","solver","collision","deform","membrane","adhesion"]) for n in imports)
    assert contract["finite_size_wall_clearance"]=="NOT_VALIDATED_PARTICLE3"
    assert contract["no_rbc_deformation"] and not contract["particle3_started"]
