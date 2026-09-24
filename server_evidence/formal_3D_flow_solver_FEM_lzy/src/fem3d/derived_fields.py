"""Derived SI fields evaluated at reference-tetra centroids (DG0 samples)."""


def derive(velocity,cache_dir):
    import basix.ufl
    import ufl
    from dolfinx import fem
    domain=velocity.function_space.mesh
    tensor=fem.functionspace(domain,basix.ufl.element("DG",domain.basix_cell(),0,shape=(3,3)))
    vector=fem.functionspace(domain,basix.ufl.element("DG",domain.basix_cell(),0,shape=(3,)))
    fields={}
    for name,expression,space in (("velocity_gradient_s_inv",ufl.grad(velocity),tensor),
        ("strain_rate_tensor_s_inv",ufl.sym(ufl.grad(velocity)),tensor),
        ("vorticity_s_inv",ufl.curl(velocity),vector)):
        field=fem.Function(space,name=name)
        compiled=fem.Expression(expression,space.element.interpolation_points,jit_options={"cache_dir":str(cache_dir)})
        field.interpolate(compiled)
        field.x.scatter_forward()
        fields[name]=field
    metadata={"gradient_definition":"gradient[i,j] = d u_i / d x_j",
        "tensor_component_order":["xx","xy","xz","yx","yy","yz","zx","zy","zz"],
        "strain_definition":"0.5 * (gradient + transpose(gradient))",
        "vorticity_definition":"curl(u); Cartesian components x,y,z", "units":"s^-1",
        "representation":"DG0 values sampled at mapped reference-tetra centroids; not claimed to represent the varying derivative everywhere in a curved cell",
        "wss_computed":False}
    return fields,metadata
