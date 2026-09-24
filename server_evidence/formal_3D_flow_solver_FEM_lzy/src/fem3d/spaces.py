"""Native [P2]^3 / P1 / Real spaces, with one owned global Real DOF."""


def create_spaces(domain):
    import basix.ufl
    from dolfinx import fem
    import ufl
    velocity = fem.functionspace(domain, basix.ufl.element("Lagrange", domain.basix_cell(), 2, shape=(3,)))
    pressure = fem.functionspace(domain, basix.ufl.element("Lagrange", domain.basix_cell(), 1))
    real = fem.functionspace(domain, basix.ufl.real_element(domain.basix_cell()))
    mixed = ufl.MixedFunctionSpace(velocity, pressure, real)
    assert real.dofmap.index_map.size_global * real.dofmap.index_map_bs == 1
    return velocity, pressure, real, mixed


def real_owned_value(function):
    """Reduce the single owner's value; ghosts are never independent unknowns."""
    from mpi4py import MPI
    space = function.function_space
    count = space.dofmap.index_map.size_local * space.dofmap.index_map_bs
    assert count in (0, 1)
    return float(space.mesh.comm.allreduce(float(function.x.array[:count].sum()), op=MPI.SUM))
