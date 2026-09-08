"""Read-only vessel import, explicit boundary identity, and lossless packages.

No legacy project, reconstruction code, or solver is imported here.
"""

from .model import BoundaryPatch, GeometryError, VesselGeometry

__all__ = ["BoundaryPatch", "GeometryError", "VesselGeometry"]
