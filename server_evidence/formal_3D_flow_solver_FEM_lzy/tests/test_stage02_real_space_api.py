"""Run the native API probe before any Stokes solver is implemented."""
import json
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def run_probe(output):
    import basix
    import dolfinx
    import numpy as np
    import ufl
    from mpi4py import MPI
    from petsc4py import PETSc
    from dolfinx import fem, mesh
    from dolfinx.fem import petsc
    from fem3d.audit import timestamp, write_json
    from fem3d.spaces import create_spaces, real_owned_value

    comm = MPI.COMM_WORLD
    domain = mesh.create_unit_cube(comm, 2, 2, 2, cell_type=mesh.CellType.tetrahedron)
    V, P, R, W = create_spaces(domain)
    sizes = [s.dofmap.index_map.size_global * s.dofmap.index_map_bs for s in (V, P, R)]
    owned = comm.allgather(R.dofmap.index_map.size_local * R.dofmap.index_map_bs)
    assert sizes == [375, 27, 1] and sum(owned) == 1
    inlet = mesh.locate_entities_boundary(domain, 2, lambda x: np.isclose(x[0], 0))
    facet_tags = mesh.meshtags(domain, 2, np.sort(inlet), np.full(len(inlet), 4, dtype=np.int32))
    ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)
    dx = ufl.Measure("dx", domain=domain)
    normal = ufl.FacetNormal(domain)
    velocity, pressure, multiplier = ufl.TrialFunctions(W)
    test_velocity, test_pressure, test_real = ufl.TestFunctions(W)
    # This is an assembly probe with mass diagonals, NOT the Stokes PDE.
    block = [[ufl.inner(velocity, test_velocity)*dx, -pressure*ufl.div(test_velocity)*dx,
              multiplier*ufl.dot(test_velocity, normal)*ds(4)],
             [test_pressure*ufl.div(velocity)*dx, pressure*test_pressure*dx, None],
             [test_real*ufl.dot(velocity, normal)*ds(4), None, multiplier*test_real*dx]]
    cache = str(ROOT / "outputs/stage02/jit_cache")
    compiled = fem.form(block, jit_options={"cache_dir": cache})
    matrix = petsc.assemble_matrix(compiled, kind="nest")
    matrix.assemble()
    dims = [[list(matrix.getNestSubMatrix(i,j).getSize()) if block[i][j] is not None else None for j in range(3)] for i in range(3)]
    for i in range(3):
        for j in range(3):
            if dims[i][j] is not None:
                assert dims[i][j] == [sizes[i], sizes[j]]
    column = matrix.getNestSubMatrix(0, 2)
    row = matrix.getNestSubMatrix(2, 0)
    one = column.createVecRight()
    one.set(1)
    boundary_load = column.createVecLeft()
    column.mult(one, boundary_load)
    assert abs(boundary_load.sum() + 1) < 1e-13
    # PETSc transpose() without an output argument mutates the matrix itself.
    # Copy the block first so this check cannot zero the original C block.
    transpose = column.copy()
    transpose.transpose()
    transpose.axpy(-1, row, structure=PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN)
    assert transpose.norm() < 1e-13
    assert column.norm() > 0
    real_mass = matrix.getNestSubMatrix(2, 2)
    integrated_one = real_mass.createVecLeft()
    real_mass.mult(one, integrated_one)
    assert abs(integrated_one.sum() - 1) < 1e-13
    monolithic = petsc.assemble_matrix(compiled, kind="mpi")
    monolithic.assemble()
    assert monolithic.getSize() == (sum(sizes), sum(sizes))
    scalar = fem.Function(R)
    count = R.dofmap.index_map.size_local
    scalar.x.array[:count] = 7.25
    scalar.x.scatter_forward()
    assert real_owned_value(scalar) == 7.25
    result = {"status": "PASS", "timestamp": timestamp(), "mpi_ranks": comm.size,
        "dolfinx_version": dolfinx.__version__, "basix_version": basix.__version__,
        "petsc_version": list(PETSc.Sys.getVersion()), "gpu_used": False,
        "spaces": {"velocity_degree": 2, "velocity_components": 3, "pressure_degree": 1,
                   "velocity_global_dofs": sizes[0], "pressure_global_dofs": sizes[1], "real_global_dofs": sizes[2]},
        "real_owned_dofs_per_rank": owned,
        "real_ghost_dofs_per_rank": comm.allgather(R.dofmap.index_map.num_ghosts),
        "block_dimensions": dims, "monolithic_dimensions": list(monolithic.getSize()),
        "C_nonzero": column.norm() > 0, "C_transpose_defect": transpose.norm(),
        "unit_velocity_inlet_flux": float(boundary_load.sum()), "real_mass_integral": float(integrated_one.sum()),
        "shared_real_value": real_owned_value(scalar), "pde_solved": False,
        "api": "basix.ufl.real_element + dolfinx.fem.functionspace; ufl.MixedFunctionSpace; dolfinx.fem.petsc.assemble_matrix(kind=nest/mpi)",
        "ownership": "Native DOLFINx index map owns the scalar once globally; other ranks carry ghost copies only"}
    if comm.rank == 0:
        if Path(output).exists():
            import shutil
            import os
            previous = Path(output).with_name(Path(output).stem + "_before_" + os.environ.get("FEM3D_RUN_ID", "rerun") + ".json")
            shutil.copy2(output, previous)
        write_json(output, result)
        print(json.dumps(result, indent=2))
    for obj in (integrated_one, transpose, boundary_load, one, monolithic, matrix):
        obj.destroy()
    return result


def test_native_real_api_one_and_two_ranks():
    for ranks in (1, 2):
        result = json.loads((ROOT / f"outputs/stage02/api/real_r{ranks}.json").read_text())
        assert result["status"] == "PASS" and result["mpi_ranks"] == ranks
        assert result["spaces"]["real_global_dofs"] == sum(result["real_owned_dofs_per_rank"]) == 1
        assert result["spaces"]["velocity_degree"] == 2 and result["spaces"]["pressure_degree"] == 1
        assert result["C_nonzero"] and abs(result["unit_velocity_inlet_flux"] + 1) < 1e-13
        assert not result["pde_solved"]


if __name__ == "__main__":
    from mpi4py import MPI
    run_probe(ROOT / f"outputs/stage02/api/real_r{MPI.COMM_WORLD.size}.json")
