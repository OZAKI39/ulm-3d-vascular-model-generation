"""Original P1 wall shear stress: tau=(I-nn^T) mu (grad U+grad U^T)n, in Pa."""
from .solver_support.wss import (p1_gradients, boundary_owners, wall_geometry,
    tangential_traction, nodal_average, surface)
