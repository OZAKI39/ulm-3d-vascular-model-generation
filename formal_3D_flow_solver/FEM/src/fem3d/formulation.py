"""Symmetric full-stress Stokes + one native global flux multiplier."""


def block_forms(domain,facet_tags,spaces,mu,signed_flow,length_scale,volume,quadrature_degree=6,reference_geometry=None):
    import numpy as np
    import ufl
    V,P,R,W=spaces
    trial_velocity,trial_pressure,trial_real=ufl.TrialFunctions(W)
    test_velocity,test_pressure,test_real=ufl.TestFunctions(W)
    # Algebraic congruence scaling only: mesh coordinates and all physics stay SI.
    # Physical u=U*uh, p=P*ph, lambda=P*lh; test scales match trial scales.
    velocity_scale=abs(signed_flow)/(np.pi*length_scale**2)
    pressure_scale=mu*velocity_scale/length_scale
    power_scale=mu*velocity_scale**2*length_scale
    velocity=velocity_scale*trial_velocity
    pressure=pressure_scale*trial_pressure
    multiplier=pressure_scale*trial_real
    test_v=velocity_scale*test_velocity
    test_q=pressure_scale*test_pressure
    test_eta=pressure_scale*test_real
    dx=ufl.Measure("dx",domain=domain,metadata={"quadrature_degree":quadrature_degree})
    ds=ufl.Measure("ds",domain=domain,subdomain_data=facet_tags,metadata={"quadrature_degree":quadrature_degree})
    normal=ufl.FacetNormal(domain)
    # The continuity equation is multiplied by -1 to obtain B=-D.
    # This does not change the physical pressure sign or divergence constraint.
    a=[[2*mu/power_scale*ufl.inner(ufl.sym(ufl.grad(velocity)),ufl.sym(ufl.grad(test_v)))*dx,
        -pressure/power_scale*ufl.div(test_v)*dx,multiplier/power_scale*ufl.dot(test_v,normal)*ds(4)],
       [-test_q/power_scale*ufl.div(velocity)*dx,None,None],
       [test_eta/power_scale*ufl.dot(velocity,normal)*ds(4),None,None]]
    from dolfinx import fem
    from petsc4py import PETSc
    zero=fem.Constant(domain,PETSc.ScalarType(0))
    load=zero*ufl.dot(test_v,ufl.as_vector((1.,1.,1.)))*dx
    if reference_geometry is not None:
        # Authorized analytic-only case. No velocity data are prescribed.
        from .analytic import pipe_frame
        origin,axis,_,_=pipe_frame(reference_geometry["origin_m"],reference_geometry["axis"])
        axis=ufl.as_vector(axis)
        offset=ufl.SpatialCoordinate(domain)-ufl.as_vector(origin)
        radial=offset-ufl.dot(offset,axis)*axis
        radius=reference_geometry["radius_m"]
        mean=signed_flow/(np.pi*radius**2)
        tangential_traction=ufl.dot(normal,axis)*(-4*mu*mean/radius**2)*radial
        load=ufl.dot(tangential_traction,test_v)*(ds(4)+ds(2))
    rhs=[(1/power_scale)*load,zero*test_q/power_scale*dx,
         -signed_flow/(volume*power_scale)*test_eta*dx]
    return a,rhs,{"velocity_m_s":velocity_scale,"pressure_pa":pressure_scale,"power_w":power_scale,
        "description":"Symmetric algebraic congruence scaling; physical SI fields restored after solve"}
