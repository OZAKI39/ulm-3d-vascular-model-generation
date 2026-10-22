"""Wall interfaces. Geometric projection is not hydrodynamic wall interaction."""

from abc import ABC, abstractmethod
from dataclasses import dataclass

from mevo_mb.coupling.geometry import Cylinder


class WallInteractionModel(ABC):
    @abstractmethod
    def apply(self, positions_m, radii_m, velocity_m_s, time_s):
        """Return feasible positions and velocities, plus model status."""


@dataclass
class GeometricNoPenetration(WallInteractionModel):
    cylinder: Cylinder

    def apply(self, positions_m, radii_m, velocity_m_s, time_s):
        x, hit = self.cylinder.project(positions_m, radii_m)
        v = self.cylinder.constrain_velocity(x, radii_m, velocity_m_s)
        return x, v, {"physics_status": "demo only", "wall_hit": hit, "model": "geometric no penetration"}


class DoinikovWallReference(WallInteractionModel):
    """TODO: PMB60,7909–7925; DOI10.1088/0031-9155/60/20/7909.

    Requires Marmottant shell, planar-wall boundary conditions, separation
    convergence and a reference curve. No mirror-bubble approximation is
    silently substituted. Deliberately unavailable in production.
    """

    def apply(self, positions_m, radii_m, velocity_m_s, time_s):
        raise NotImplementedError("Doinikov 2015 wall reference has not been implemented/validated")
