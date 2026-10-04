"""Single-inlet/outlet integrals for the unchanged production steady monitor."""
import numpy as np
from scipy.spatial import cKDTree
from vascular_validation.postprocess import SolutionMeasurements
from vascular_validation.validation import triangle_flux, parse_result_vtu


def serialized_node_mapping(points, tetra, grid):
    """Prove a bijection at the writer's coordinate precision and identical cells.

    svMultiPhysics VTK coordinates may be float32 even for a float64 mesh. This
    compares exact representations, not a relaxed nearest-neighbor tolerance.
    The original physical mesh remains the integration/interpolation geometry.
    """
    xyz=np.asarray(grid.points)
    if len(xyz)!=len(points) or grid.n_cells!=len(tetra) or not np.all(grid.celltypes==10):
        raise ValueError('Solution mesh dimensions or cell types differ')
    encoded=np.asarray(points,dtype=xyz.dtype)
    if not np.isfinite(xyz).all() or len(np.unique(encoded,axis=0))!=len(points):
        raise ValueError('Nonfinite coordinates or serialization coordinate collision')
    distance,index=cKDTree(xyz.astype(float)).query(encoded.astype(float))
    if not np.all(distance==0) or len(np.unique(index))!=len(index):
        raise ValueError('Solution coordinates differ from exact mesh serialization')
    inverse=np.empty(len(index),dtype=int);inverse[index]=np.arange(len(index))
    actual=inverse[grid.cells.reshape(-1,5)[:,1:]]
    def canonical(cells):
        rows=np.sort(cells,axis=1)
        return rows[np.lexsort(rows.T[::-1])]
    if not np.array_equal(canonical(actual),canonical(tetra)):
        raise ValueError('Solution tetra connectivity differs from physical mesh')
    return index


class PipeMeasurements(SolutionMeasurements):
    def read(self,path):
        grid,velocity,pressure=parse_result_vtu(path)
        index=serialized_node_mapping(self.points,self.tetra,grid)
        self.last_output_identity=dict(coordinate_dtype=str(grid.points.dtype),
            maximum_coordinate_serialization_error_m=float(np.linalg.norm(grid.points[index]-self.points,axis=1).max()),
            matched_nodes=len(index),matched_tetra=len(self.tetra),bijective=True,connectivity_equal=True)
        return np.asarray(velocity[index],float),np.asarray(pressure[index],float).reshape(-1)

    def measure(self, velocity, pressure):
        flows, average_pressure = {}, {}
        for tag, role in [(2, 'INLET'), (3, 'OUTLET_01')]:
            face = self.boundary[self.tags == tag]
            flows[role] = triangle_flux(self.points, face, velocity)
            xyz = self.points[face]
            area = .5 * np.linalg.norm(np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]), axis=1)
            average_pressure[role] = float(np.dot(area, pressure[face].mean(axis=1)) / area.sum())
        incoming, outgoing = -flows['INLET'], flows['OUTLET_01']
        wall_speed = np.linalg.norm(velocity[self.wall_nodes], axis=1)
        return dict(Q_target_m3_s=self.Q, Q_in_m3_s=incoming, Q_out_total_m3_s=outgoing,
                    outlet_flows_m3_s={'OUTLET_01': outgoing}, signed_outward_boundary_flows_m3_s=flows,
                    epsilon_Q=abs(incoming - self.Q) / self.Q, epsilon_mass=abs(outgoing - incoming) / self.Q,
                    outlet_fractions={'OUTLET_01': 1.} if outgoing != 0 else None,
                    area_average_pressure_pa=average_pressure, velocity_L2=self.velocity_l2(velocity),
                    velocity_max_m_s=float(np.linalg.norm(velocity, axis=1).max()),
                    pressure_range_pa=[float(pressure.min()), float(pressure.max())],
                    velocity_finite=bool(np.isfinite(velocity).all()), pressure_finite=bool(np.isfinite(pressure).all()),
                    wall_velocity_max_m_s=float(wall_speed.max()), wall_velocity_P95_m_s=float(np.quantile(wall_speed, .95)),
                    wall_noslip_tolerance_m_s=1e-10 * self.Umean,
                    wall_noslip_pass=bool(wall_speed.max() <= 1e-10 * self.Umean))
