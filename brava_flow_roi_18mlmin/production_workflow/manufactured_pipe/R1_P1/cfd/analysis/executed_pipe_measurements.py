"""Single-inlet/outlet integrals for the unchanged production steady monitor."""
import numpy as np
from vascular_validation.postprocess import SolutionMeasurements
from vascular_validation.validation import triangle_flux


class PipeMeasurements(SolutionMeasurements):
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
