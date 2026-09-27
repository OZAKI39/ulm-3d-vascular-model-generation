"""Independent integrals on the actual native solution, in SI units."""
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from .geometry import ROLE_IDS
from .validation import parse_result_vtu, require, triangle_flux


class SolutionMeasurements:
    def __init__(self, arrays_path, Q_target, Umean):
        with np.load(arrays_path) as data:
            self.points = data['points_m'].copy()
            self.tetra = data['tetra'].copy()
            self.boundary = data['boundary_triangles'].copy()
            self.tags = data['facet_tags'].copy()
        self.Q = Q_target
        self.Umean = Umean
        xyz = self.points[self.tetra]
        self.volumes = np.linalg.det(xyz[:, 1:]-xyz[:, :1])/6.
        self.tree = cKDTree(self.points)
        self.wall_nodes = np.unique(self.boundary[self.tags == 1])

    def read(self, path):
        grid, velocity, pressure = parse_result_vtu(path)
        require(grid.n_points == len(self.points), 'Solution and actual mesh point counts differ')
        require(grid.n_cells == len(self.tetra), 'Solution and actual mesh cell counts differ')
        if not np.array_equal(np.asarray(grid.points), self.points):
            distances, indices = cKDTree(np.asarray(grid.points, float)).query(self.points)
            tolerance = 1e-10 * np.ptp(self.points, axis=0).max()
            require(distances.max() <= tolerance and len(np.unique(indices)) == len(indices), 'Solution geometry cannot be matched uniquely')
            velocity, pressure = velocity[indices], pressure[indices]
        return np.asarray(velocity, float), np.asarray(pressure, float).reshape(-1)

    def velocity_l2(self, velocity):
        values = velocity[self.tetra]
        squared = (np.sum(values*values, axis=(1, 2)) + np.sum(values.sum(axis=1)**2, axis=1)) / 20.
        return float(np.sqrt(np.dot(self.volumes, squared)))

    def measure(self, velocity, pressure):
        flows, average_pressure = {}, {}
        for face_id, role in ROLE_IDS.items():
            if role == 'WALL':
                continue
            triangles = self.boundary[self.tags == face_id]
            flows[role] = triangle_flux(self.points, triangles, velocity)
            xyz = self.points[triangles]
            areas = .5*np.linalg.norm(np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0]), axis=1)
            average_pressure[role] = float(np.dot(areas, pressure[triangles].mean(axis=1))/areas.sum())
        incoming = -flows['INLET']
        outlets = {name: flows[name] for name in ('OUTLET_01', 'OUTLET_02', 'OUTLET_03')}
        total = sum(outlets.values())
        wall_speed = np.linalg.norm(velocity[self.wall_nodes], axis=1)
        return {'Q_target_m3_s': self.Q, 'Q_in_m3_s': incoming,
                'outlet_flows_m3_s': outlets, 'Q_out_total_m3_s': total,
                'signed_outward_boundary_flows_m3_s': flows,
                'epsilon_Q': abs(incoming-self.Q)/self.Q,
                'epsilon_mass': abs(total-incoming)/self.Q,
                'outlet_fractions': {name: value/total for name, value in outlets.items()} if total != 0 else None,
                'area_average_pressure_pa': average_pressure,
                'velocity_L2': self.velocity_l2(velocity),
                'velocity_max_m_s': float(np.linalg.norm(velocity, axis=1).max()),
                'pressure_range_pa': [float(pressure.min()), float(pressure.max())],
                'velocity_finite': bool(np.isfinite(velocity).all()),
                'pressure_finite': bool(np.isfinite(pressure).all()),
                'wall_velocity_max_m_s': float(wall_speed.max()),
                'wall_velocity_P95_m_s': float(np.quantile(wall_speed, .95)),
                'wall_noslip_tolerance_m_s': 1e-10*self.Umean,
                'wall_noslip_pass': bool(wall_speed.max() <= 1e-10*self.Umean)}
