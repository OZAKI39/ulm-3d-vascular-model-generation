"""Only wall velocity is essential; inlet/outlet velocities remain free."""
import numpy as np


def wall_condition(velocity_space,facet_tags,wall_tag=1):
    from dolfinx import fem
    from petsc4py import PETSc
    from mpi4py import MPI
    dofs=fem.locate_dofs_topological(velocity_space,2,facet_tags.find(wall_tag))
    condition=fem.dirichletbc(np.zeros(3,dtype=PETSc.ScalarType),dofs,velocity_space)
    owned=velocity_space.dofmap.index_map.size_local
    count=velocity_space.mesh.comm.allreduce(int(np.count_nonzero(dofs<owned))*3,op=MPI.SUM)
    return condition,{"essential_boundary":"WALL only (includes the wall/cap intersection rim)",
        "velocity_dirichlet_dofs":int(count),"pressure_dirichlet_dofs":0,
        "inlet_velocity_profile":None,"outlet_velocity_profile":None,"pressure_nullspace_registered":False}
