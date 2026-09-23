"""Integrate computed fields and sample disk quadrature in the pipe's local frame."""
import numpy as np
from .analytic import pipe_frame,pipe_theory,poiseuille_at


def integrator(domain,cache_dir):
    from mpi4py import MPI
    from dolfinx import fem
    def integrate(expression):
        value=fem.assemble_scalar(fem.form(expression,jit_options={"cache_dir":str(cache_dir)}))
        return float(domain.comm.allreduce(value,op=MPI.SUM))
    return integrate


def flux_diagnostics(domain,facet_tags,velocity,signed_flow,cache_dir):
    import ufl
    integrate=integrator(domain,cache_dir)
    ds=ufl.Measure("ds",domain=domain,subdomain_data=facet_tags,metadata={"quadrature_degree":6})
    normal=ufl.FacetNormal(domain)
    inflow=integrate(ufl.dot(velocity,normal)*ds(4))
    outflow=integrate(ufl.dot(velocity,normal)*ds(2))
    wall=integrate(ufl.dot(velocity,normal)*ds(1))
    constraint=inflow+signed_flow
    closure=inflow+outflow
    constraint_rel=abs(constraint)/abs(signed_flow)
    closure_rel=abs(closure)/abs(signed_flow)
    return {"target_Q_m3_s":signed_flow,"Q_in_signed_m3_s":inflow,"Q_out_signed_m3_s":outflow,
        "actual_Q_in_m3_s":-inflow,"actual_Q_out_m3_s":outflow,"wall_flux_m3_s":wall,
        "constraint_residual_m3_s":constraint,"mass_closure_residual_m3_s":closure,
        "relative_inlet_constraint_error":constraint_rel,"relative_mass_closure":closure_rel,
        "gate":1e-10,"status":"PASS" if max(constraint_rel,closure_rel)<=1e-10 else "FAIL",
        "method":"FEM solution boundary integral using FacetNormal, tagged ds and one MPI SUM reduction"}


def evaluate_points(function,points):
    from dolfinx import geometry
    from mpi4py import MPI
    domain=function.function_space.mesh
    comm=domain.comm
    points=np.asarray(points,dtype=domain.geometry.x.dtype)
    tree=geometry.bb_tree(domain,domain.topology.dim,padding=1e-13)
    candidates=geometry.compute_collisions_points(tree,points)
    # A convex-hull collision alone is insufficient for curved P2 tetrahedra.
    # Confirm membership in the actual reference tetra using the supported
    # nonlinear coordinate pull-back API, then evaluate in that same cell.
    cells=np.full(len(points),-1,dtype=np.int32)
    for i,point in enumerate(points):
        for cell in candidates.links(i):
            coordinates=domain.geometry.x[domain.geometry.dofmaps[0][cell]]
            reference=domain.geometry.cmap.pull_back(point.reshape(1,3),coordinates,tol=1e-12,maxit=30)[0]
            if reference.min()>=-1e-10 and reference.sum()<=1+1e-10:
                cells[i]=cell
                break
    owners=np.where(cells>=0,comm.rank,comm.size).astype(np.int32)
    final_owners=np.empty_like(owners)
    comm.Allreduce(owners,final_owners,op=MPI.MIN)
    if np.any(final_owners==comm.size):
        bad=points[final_owners==comm.size]
        raise ValueError(f"FEM sampling points not found in mesh: {len(bad)}, first={bad[0]}")
    local=np.flatnonzero(final_owners==comm.rank)
    shape=function.function_space.element.value_shape
    width=int(np.prod(shape)) if len(shape) else 1
    values=np.zeros((len(points),width),dtype=float)
    if len(local):
        values[local]=function.eval(points[local],cells[local],tol=1e-12,maxit=30).reshape(-1,width)
    result=np.empty_like(values)
    comm.Allreduce(values,result,op=MPI.SUM)
    return result


def analytic_diagnostics(domain,velocity,pressure,multiplier,geometry,mu,signed_flow,rho,cache_dir):
    import ufl
    origin,axis,first,second=pipe_frame(geometry["origin_m"],geometry["axis"])
    radius,length=geometry["radius_m"],geometry["length_m"]
    theory=pipe_theory(radius,length,mu,signed_flow,rho)
    radial_nodes,radial_weights=np.polynomial.legendre.leggauss(12)
    radial=radius*np.sqrt((radial_nodes+1)/2)
    theta=np.arange(64)*2*np.pi/64
    offsets=(radial[:,None,None]*(np.cos(theta)[None,:,None]*first+np.sin(theta)[None,:,None]*second)).reshape(-1,3)
    weights=np.repeat(radial_weights/2/len(theta),len(theta))
    sections=[.25,.5,.75]
    velocity_rows=[]
    profile_rows=[]
    error_sum=0.; exact_sum=0.
    for fraction in sections:
        points=origin+fraction*length*axis+offsets
        values=evaluate_points(velocity,points)
        exact,_=poiseuille_at(points,origin,axis,radius,length,mu,signed_flow)
        error=float(np.sum(weights*np.sum((values-exact)**2,axis=1)))
        norm=float(np.sum(weights*np.sum(exact**2,axis=1)))
        axial=values@axis
        transverse=values-axial[:,None]*axis
        velocity_rows.append({"s_over_L":fraction,"relative_L2_error":float(np.sqrt(error/norm)),
            "relative_max_error":float(np.linalg.norm(values-exact,axis=1).max()/abs(theory["centerline_velocity_m_s"])),
            "maximum_transverse_velocity_over_Umean":float(np.linalg.norm(transverse,axis=1).max()/abs(theory["mean_velocity_m_s"]))})
        error_sum+=error; exact_sum+=norm
        diameter=np.linspace(-.98,.98,121)
        sample_points=origin+fraction*length*axis+diameter[:,None]*radius*first
        sampled=evaluate_points(velocity,sample_points)
        analytic,_=poiseuille_at(sample_points,origin,axis,radius,length,mu,signed_flow)
        profile_rows.append({"s_over_L":fraction,"r_over_R":diameter.tolist(),
            "axial_velocity_m_s":(sampled@axis).tolist(),"analytic_axial_velocity_m_s":(analytic@axis).tolist(),
            "transverse_velocity_m_s":np.linalg.norm(sampled-(sampled@axis)[:,None]*axis,axis=1).tolist()})
    pressure_rows=[]
    for fraction in np.linspace(.1,.9,9):
        points=origin+fraction*length*axis+offsets
        values=evaluate_points(pressure,points)[:,0]
        expected=theory["delta_pressure_pa"]*(1-fraction)
        pressure_rows.append({"s_over_L":float(fraction),"s_m":float(fraction*length),
            "mean_pressure_pa":float(weights@values),"analytic_pressure_pa":float(expected)})
    actual_pressure=np.array([r["mean_pressure_pa"] for r in pressure_rows])
    exact_pressure=np.array([r["analytic_pressure_pa"] for r in pressure_rows])
    slope,intercept=np.polyfit(np.array([r["s_m"] for r in pressure_rows]),actual_pressure,1)
    center_velocity=float(evaluate_points(velocity,np.array([origin+.5*length*axis]))[0]@axis)
    dx=ufl.Measure("dx",domain=domain,metadata={"quadrature_degree":6})
    integrate=integrator(domain,cache_dir)
    x=ufl.SpatialCoordinate(domain)
    offset=x-ufl.as_vector(origin)
    s=ufl.dot(offset,ufl.as_vector(axis))
    radial2=ufl.dot(offset,offset)-s*s
    exact_u=2*theory["mean_velocity_m_s"]*(1-radial2/radius**2)*ufl.as_vector(axis)
    exact_p=theory["delta_pressure_pa"]*(1-s/length)
    global_u_error=np.sqrt(integrate(ufl.inner(velocity-exact_u,velocity-exact_u)*dx)/integrate(ufl.inner(exact_u,exact_u)*dx))
    global_p_error=np.sqrt(integrate((pressure-exact_p)**2*dx)/integrate(exact_p**2*dx))
    return {"lambda_pa":multiplier,"analytic_delta_p_pa":theory["delta_pressure_pa"],
        "lambda_relative_error":abs(multiplier-theory["delta_pressure_pa"])/abs(theory["delta_pressure_pa"]),
        "velocity_L2_relative_error":float(np.sqrt(error_sum/exact_sum)),
        "velocity_global_L2_relative_error":float(global_u_error),"pressure_global_L2_relative_error":float(global_p_error),
        "pressure_profile_error":float(np.linalg.norm(actual_pressure-exact_pressure)/np.linalg.norm(exact_pressure)),
        "pressure_drop_from_internal_section_fit_pa":float(-slope*length),"pressure_fit_outlet_extrapolation_pa":float(intercept+slope*length),
        "centerline_velocity_m_s":center_velocity,"analytic_centerline_velocity_m_s":theory["centerline_velocity_m_s"],
        "reynolds_number":theory["reynolds_number"],"mean_velocity_m_s":theory["mean_velocity_m_s"],
        "velocity_L2_norm_m_pow_2p5_s":float(np.sqrt(integrate(ufl.inner(velocity,velocity)*dx))),
        "pressure_integral_pa_m3":integrate(pressure*dx),"divergence_L2":float(np.sqrt(integrate(ufl.div(velocity)**2*dx))),
        "velocity_sections":velocity_rows,"pressure_sections":pressure_rows,"velocity_profiles":profile_rows,
        "point_location_method":"Native nonlinear pull-back verifies reference-tetra membership; eval tol=1e-12",
        "section_quadrature":{"radial_gauss_points":12,"azimuthal_points":64,"disk_radius_m":radius,
            "description":"Area-weighted Gauss quadrature on nominal analytic disks; every quadrature point located in actual FEM mesh. Mean over many points, never a single pressure node."}}
