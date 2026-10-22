"""Finite-wall/contact entry points; zero-thickness rigid walls, no new adhesion model."""
from .wall_geometry import WallGeometry
from .nearfield_regularization import NearFieldRegularizationV1
from .nearfield_handoff import exact_interactions
from .kinematic_contact import project_contacts
