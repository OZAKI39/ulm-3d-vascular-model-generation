"""General SI Stokes core. Density is deliberately absent from this API."""
import time
import resource
import numpy as np
from .spaces import create_spaces,real_owned_value
from .boundary import wall_condition
from .formulation import block_forms


def solve_stokes(domain,facet_tags,mu_pa_s,Q_target_m3_s,length_scale_m,cache_dir):
    """Production-facing API: positive inflow magnitude, natural end tractions."""
    if not np.isfinite(Q_target_m3_s) or Q_target_m3_s<=0:
        raise ValueError("Q_target must be a positive inflow magnitude in m^3/s")
    return _solve(domain,facet_tags,mu_pa_s,Q_target_m3_s,length_scale_m,cache_dir)


def _solve(domain,facet_tags,mu,signed_flow,length_scale,cache_dir,reference_geometry=None):
    """Internal regression interface allows reversed RHS and analytic traction."""
    import dolfinx
    import basix
    import ufl
    from dolfinx import fem
    from dolfinx.fem.petsc import LinearProblem
    from mpi4py import MPI
    from petsc4py import PETSc
    if not np.isfinite([mu,signed_flow,length_scale]).all() or mu<=0 or signed_flow==0 or length_scale<=0:
        raise ValueError("Invalid SI viscosity, nonzero regression flow or length scale")
    start=time.perf_counter()
    spaces=create_spaces(domain)
    V,P,R,_=spaces
    sizes=[s.dofmap.index_map.size_global*s.dofmap.index_map_bs for s in (V,P,R)]
    if sum(sizes)>1500000:
        raise RuntimeError("Development DOF limit exceeded; no automatic resource escalation")
    bc,boundary_metadata=wall_condition(V,facet_tags)
    options={"cache_dir":str(cache_dir)}
    volume=domain.comm.allreduce(fem.assemble_scalar(fem.form(1*ufl.dx(domain=domain,metadata={"quadrature_degree":6}),jit_options=options)),op=MPI.SUM)
    a,rhs,scales=block_forms(domain,facet_tags,spaces,mu,signed_flow,length_scale,volume,reference_geometry=reference_geometry)
    fields=[fem.Function(V,name="velocity"),fem.Function(P,name="pressure_gauge_pa"),fem.Function(R,name="inlet_normal_traction_multiplier")]
    problem=LinearProblem(a,rhs,u=fields,bcs=[bc],kind="mpi",petsc_options_prefix="stage02_stokes_",
        petsc_options={"ksp_type":"preonly","pc_type":"lu","pc_factor_mat_solver_type":"mumps","ksp_error_if_not_converged":True},jit_options=options)
    assembled_start=time.perf_counter()
    problem.solve()
    solve_seconds=time.perf_counter()-assembled_start
    reason=int(problem.solver.getConvergedReason())
    residual=problem.b.duplicate()
    problem.A.mult(problem.x,residual)
    residual.axpy(-1,problem.b)
    residual_norm=float(residual.norm())
    relative_residual=residual_norm/max(float(problem.b.norm()),1e-300)
    residual.destroy()
    for field,scale in zip(fields,(scales["velocity_m_s"],scales["pressure_pa"],scales["pressure_pa"])):
        field.x.array[:]*=scale
        field.x.scatter_forward()
    finite=domain.comm.allreduce(all(np.isfinite(f.x.array).all() for f in fields),op=MPI.LAND)
    multiplier=real_owned_value(fields[2])
    assert finite and reason>0
    assert multiplier*signed_flow>0, "Multiplier sign must follow the prescribed signed flow"
    result={"status":"PASS","converged_reason":reason,"ksp_type":problem.solver.getType(),
        "pc_type":problem.solver.getPC().getType(),"factor_solver":problem.solver.getPC().getFactorSolverType(),
        "relative_algebraic_residual":relative_residual,"residual_norm_scaled":residual_norm,
        "assembly_and_factor_solve_s":solve_seconds,"wall_time_s":time.perf_counter()-start,
        "max_rank_peak_rss_kib":int(domain.comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,op=MPI.MAX)),
        "mpi_ranks":domain.comm.size,"gpu_used":False,"dolfinx_version":dolfinx.__version__,"basix_version":basix.__version__,
        "petsc_version":list(PETSc.Sys.getVersion()),"velocity_dofs":sizes[0],"pressure_dofs":sizes[1],"real_global_dofs":sizes[2],
        "real_owned_dofs_per_rank":domain.comm.allgather(R.dofmap.index_map.size_local),
        "boundary_conditions":boundary_metadata,"physical_fields_finite":bool(finite),"lambda_pa":multiplier,
        "algebraic_scaling":scales,"dynamic_viscosity_pa_s":mu,
        "boundary_model":"natural_traction" if reference_geometry is None else "poiseuille_traction_reference",
        "pressure_gauge":"outlet traction; no essential pressure condition or registered nullspace"}
    # Keep the high-level problem alive while returning fields; it owns PETSc objects.
    return {"velocity":fields[0],"pressure":fields[1],"multiplier":fields[2],"lambda_pa":multiplier,
            "solver":result,"problem":problem,"volume_m3":float(volume)}
