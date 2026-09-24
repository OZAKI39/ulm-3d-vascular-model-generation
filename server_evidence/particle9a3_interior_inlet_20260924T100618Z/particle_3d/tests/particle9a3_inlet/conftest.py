import numpy as np
import pyvista as pv
import pytest

from particle_3d.interior_section import root_topology, locator_for


@pytest.fixture
def tube():
    grid = pv.ImageData(dimensions=(3, 3, 5), spacing=(.5e-6, .5e-6, .25e-6),
                        origin=(-.5e-6, -.5e-6, 0)).cast_to_unstructured_grid().triangulate()
    grid.points = np.asarray(grid.points, float)
    grid['Velocity'] = np.tile([0., 0., .002], (grid.n_points, 1))
    exterior = grid.extract_surface(algorithm='dataset_surface')
    z = exterior.cell_centers().points[:, 2]
    inlet_mask, outlet_mask = z == 0, z == 1e-6
    surface = lambda mask: exterior.extract_cells(mask).extract_surface(algorithm='dataset_surface')
    caps = {'INLET': surface(inlet_mask), 'OUTLET': surface(outlet_mask)}
    wall = surface(~(inlet_mask | outlet_mask))
    swc = np.array([[1, 3, 0, 0, .25, .5, -1], [2, 3, 0, 0, .8, .5, 1],
                    [3, 3, .25, 0, 1, .5, 2], [4, 3, -.25, 0, 1, .5, 2]])
    topo = root_topology(swc, np.zeros(3))
    candidate = dict(candidate_id=0, refinement=1, center_m=np.array([0, 0, .5e-6]),
                     normal=np.array([0., 0., 1.]), arclength_m=.5e-6,
                     distance_to_first_junction_m=.3e-6)
    return grid, candidate, topo, locator_for(wall), {k: locator_for(v) for k, v in caps.items()}
