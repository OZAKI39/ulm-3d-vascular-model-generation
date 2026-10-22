"""Original per-facet P1 wall traction; nodal averaging is display-only."""
from flow_solver_support.wss_case import recover

def calculate(field):
    a=field.arrays
    return recover(a['points_m'],a['tetra'],a['boundary_triangles'],a['facet_tags'],
                   field.velocity,field.pressure,field.policy['mu_Pa_s'])

