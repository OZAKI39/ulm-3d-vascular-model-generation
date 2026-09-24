"""Tiny 3D scalar Poisson environment test. This is NOT a blood-flow test.

Uses the official DOLFINx 0.11 LinearProblem API; authored independently:
https://docs.fenicsproject.org/dolfinx/v0.11.0.post0/python/demos/demo_poisson.html
"""
import json
import os
import resource
import time
from pathlib import Path

ENTRY_TIME = time.perf_counter()
ENTRY_EPOCH = time.time()


def run_smoke(output):
    from mpi4py import MPI
    from petsc4py import PETSc
    import dolfinx
    from dolfinx import fem, mesh
    from dolfinx.fem.petsc import LinearProblem
    import numpy as np
    import ufl

    comm = MPI.COMM_WORLD
    expected = int(os.environ.get("FEM3D_EXPECTED_RANKS", comm.size))
    assert comm.size == expected, f"MPI launcher and mpi4py disagree: {comm.size} != {expected}"
    assert dolfinx.__version__.startswith("0.11."), dolfinx.__version__
    comm.Barrier()
    imported_time = time.perf_counter()
    domain = mesh.create_unit_cube(comm, 3, 3, 3, cell_type=mesh.CellType.tetrahedron)
    V = fem.functionspace(domain, ("Lagrange", 1))
    tdim = domain.topology.dim
    domain.topology.create_connectivity(tdim-1,tdim)
    boundary = mesh.exterior_facet_indices(domain.topology)
    dofs = fem.locate_dofs_topological(V, tdim-1, boundary)
    exact = fem.Function(V)
    exact.interpolate(lambda x: x[0]+2*x[1]+3*x[2])
    bc = fem.dirichletbc(exact,dofs)
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    a = ufl.inner(ufl.grad(u),ufl.grad(v))*ufl.dx
    rhs = ufl.inner(fem.Constant(domain,PETSc.ScalarType(0)),v)*ufl.dx
    volume = comm.allreduce(fem.assemble_scalar(fem.form(fem.Constant(domain,PETSc.ScalarType(1))*ufl.dx)),op=MPI.SUM)
    problem = LinearProblem(a,rhs,bcs=[bc],petsc_options_prefix="stage00_poisson_",
        petsc_options={"ksp_type":"cg","pc_type":"jacobi","ksp_rtol":1e-11,"ksp_atol":1e-13,
                       "ksp_error_if_not_converged":True},
        jit_options={"cache_dir":str(Path("outputs/stage00/jit_cache").resolve())})
    uh = problem.solve()
    uh.x.scatter_forward()
    owned = V.dofmap.index_map.size_local*V.dofmap.index_map_bs
    finite = comm.allreduce(bool(np.isfinite(uh.x.array[:owned]).all()),op=MPI.LAND)
    norm = comm.allreduce(float(np.abs(uh.x.array[:owned]).sum()),op=MPI.SUM)
    error = float(np.sqrt(comm.allreduce(fem.assemble_scalar(fem.form(ufl.inner(uh-exact,uh-exact)*ufl.dx)),op=MPI.SUM)))
    reason = int(problem.solver.getConvergedReason())
    # An explicit algebraic residual guards against a false solver-success flag.
    residual = problem.b.duplicate()
    problem.A.mult(problem.x,residual)
    residual.axpy(-1,problem.b)
    relres = float(residual.norm()/max(problem.b.norm(),1e-300))
    residual.destroy()
    assert domain.topology.dim == 3 and abs(float(volume)-1)<1e-12
    assert finite and norm>0 and reason>0 and error<1e-10 and relres<1e-10
    comm.Barrier()
    result = {"status":"PASS","kind":"3D scalar affine Poisson infrastructure smoke, not blood flow",
              "dolfinx_version":dolfinx.__version__,"petsc_version":list(PETSc.Sys.getVersion()),
              "mpi_ranks":comm.size,"mpi_library":MPI.Get_library_version(),
              "mesh_cells":domain.topology.index_map(tdim).size_global,
              "global_dofs":V.dofmap.index_map.size_global,"assembled_volume":float(volume),
              "all_owned_solution_values_finite":finite,"solution_l2_error":error,
              "ksp_converged_reason":reason,"ksp_iterations":int(problem.solver.getIterationNumber()),
              "relative_algebraic_residual":relres,
              "max_import_and_mpi_init_s":comm.allreduce(imported_time-ENTRY_TIME,op=MPI.MAX),
              "max_assembly_and_solve_s":comm.allreduce(time.perf_counter()-imported_time,op=MPI.MAX),
              "max_rank_rss_kib":comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,op=MPI.MAX),
              "launch_to_script_entry_max_s":comm.allreduce(ENTRY_EPOCH-float(os.environ.get("FEM3D_LAUNCH_EPOCH",ENTRY_EPOCH)),op=MPI.MAX),
              "remote_work_dir":str(Path.cwd()),"input_geometry":"Analytic unit cube; vascular input not used"}
    if comm.rank==0:
        output = Path(output)
        output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(json.dumps(result,indent=2)+"\n")
        print(json.dumps(result))
    return result


def test_dolfinx_remote_smoke(tmp_path):
    import pytest
    pytest.importorskip("dolfinx", reason="Executed separately on remote 1/2/4-rank DOLFINx environment")
    assert run_smoke(tmp_path/"smoke.json")["status"] == "PASS"


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    run_smoke(args.output)
