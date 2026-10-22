"""SOFTWARE_TEST_ONLY: smooth bounded incompressible jet-splitting test field."""
import numpy as np
from scipy.special import expit
from mevo_mb.flow.fields import FlowField


class SyntheticBifurcationFlowField(FlowField):
    status = "SOFTWARE_TEST_ONLY"

    def __init__(self, geometry, config):
        self.geometry = geometry
        self.config = config
        self.maximum_velocity = float(config["centerline_velocity_m_s"])

    def velocity(self, positions_m, time_s):
        x = np.asarray(positions_m)
        length = self.config["transition_length_m"]
        width = self.config.get("jet_width_m", 0.0008)
        slope = np.tan(np.deg2rad(self.geometry.config["branch_angle_deg"]))
        center = slope * length * np.logaddexp(0, x[:, 0]/length)
        derivative = slope * expit(x[:, 0]/length)
        upper = np.exp(-((x[:, 1]-center)/width)**2)
        lower = np.exp(-((x[:, 1]+center)/width)**2)
        z_factor = np.exp(-(x[:, 2]/width)**2)
        velocity = np.zeros_like(x, dtype=float)
        velocity[:, 0] = self.maximum_velocity/2 * (upper+lower) * z_factor
        velocity[:, 1] = self.maximum_velocity/2 * derivative * (upper-lower) * z_factor
        return velocity
