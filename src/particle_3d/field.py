"""Canonical tetra sampler. All core quantities use SI and float64."""
from dataclasses import dataclass, fields
import numpy as np
from .geometry import TetraGeometry, readonly


@dataclass(frozen=True)
class FlowSample:
    velocity_m_s: np.ndarray
    pressure_pa: float
    velocity_gradient_s_inv: np.ndarray
    vorticity_s_inv: np.ndarray
    strain_rate_s_inv: np.ndarray
    inside_lumen: bool
    tetra_id: int


@dataclass(frozen=True)
class FlowBatch:
    """Same fields as FlowSample, with a leading query dimension."""
    velocity_m_s: np.ndarray
    pressure_pa: np.ndarray
    velocity_gradient_s_inv: np.ndarray
    vorticity_s_inv: np.ndarray
    strain_rate_s_inv: np.ndarray
    inside_lumen: np.ndarray
    tetra_id: np.ndarray


class FrozenFEMField:
    """Read-only nodal linear field, including the closed domain boundary.

    Multiple containing cells: smallest canonical zero-based cell ID wins.
    Outside: inside_lumen=False, tetra_id=-1, every physical quantity NaN.
    Nonfinite/malformed coordinates: ValueError (including any batch row).
    Boundary membership does not establish finite-radius particle clearance.
    No smoothing, particle state, dynamics or flow solve exists in this class.
    """

    def __init__(self, points_m, tetra, velocity_m_s, pressure_pa):
        self.geometry = TetraGeometry(points_m, tetra)
        self.points_m = self.geometry.points
        self.tetra = self.geometry.tetra
        self.velocity_nodes_m_s = readonly(velocity_m_s)
        self.pressure_nodes_pa = readonly(pressure_pa)
        n = len(self.points_m)
        if self.velocity_nodes_m_s.shape != (n, 3) or self.pressure_nodes_pa.shape != (n,):
            raise ValueError("nodal field shapes must be (N,3) and (N,)")
        if not np.isfinite(self.velocity_nodes_m_s).all() or not np.isfinite(self.pressure_nodes_pa).all():
            raise ValueError("nodal physical fields must be finite")
        # Difference form avoids cancellation from a large constant velocity.
        nodal = self.velocity_nodes_m_s[self.tetra]
        du = nodal[:, 1:] - nodal[:, :1]
        g = np.einsum("nai,naj->nij", du, self.geometry.inverse)
        self.gradients_s_inv = readonly(g)  # G[i,j] = du_i/dx_j
        self.strain_s_inv = readonly(0.5 * (g + np.swapaxes(g, 1, 2)))
        self.curl_s_inv = readonly(np.stack((g[:, 2, 1] - g[:, 1, 2],
                                             g[:, 0, 2] - g[:, 2, 0],
                                             g[:, 1, 0] - g[:, 0, 1]), axis=1))
        # VTK only enumerates intersecting cell bounding boxes.
        import pyvista as pv
        import vtk
        cells = np.column_stack((np.full(len(self.tetra), 4), self.tetra))
        self._grid = pv.UnstructuredGrid(cells.ravel(), np.full(len(self.tetra), 10, np.uint8), self.points_m.copy())
        self._locator = vtk.vtkStaticCellLocator()
        self._locator.SetDataSet(self._grid)
        self._locator.BuildLocator()
        self.provenance = None

    @classmethod
    def from_frozen(cls, fem_root):
        from .audit import read_frozen
        summary, mesh, flow, _ = read_frozen(fem_root)
        obj = cls.from_grids(mesh, flow)
        obj.provenance = summary
        return obj

    @classmethod
    def from_grids(cls, canonical_volume, flow):
        """Low-level pair loader (no hash provenance); useful for ordering tests."""
        from .audit import validate_pair
        tetra, _ = validate_pair(canonical_volume, flow)
        return cls(canonical_volume.points, tetra, flow.point_data["Velocity"], flow.point_data["Pressure"])

    @staticmethod
    def _position(position_m):
        position = np.asarray(position_m, dtype=np.float64)
        if position.shape != (3,) or not np.isfinite(position).all():
            raise ValueError("position_m must be finite with shape (3,) in m")
        return position

    def locate(self, position_m):
        """Return (canonical ID, own barycentric weights), or (-1, four NaNs)."""
        import vtk
        position = self._position(position_m)
        pad = self.geometry.candidate_padding_m
        bounds = np.column_stack((position - pad, position + pad)).ravel()
        candidates = vtk.vtkIdList()
        self._locator.FindCellsWithinBounds(bounds, candidates)
        ids = np.array([candidates.GetId(i) for i in range(candidates.GetNumberOfIds())], dtype=np.int64)
        if len(ids):
            ids.sort()
            weights = self.geometry.weights(position, ids)
            valid = np.flatnonzero(self.geometry.contains(weights, ids))
            if len(valid):
                k = valid[0]
                return int(ids[k]), weights[k]
        return -1, np.full(4, np.nan)

    def sample(self, position_m):
        cell, weights = self.locate(position_m)
        if cell < 0:
            return FlowSample(readonly(np.full(3, np.nan)), float("nan"), readonly(np.full((3, 3), np.nan)),
                              readonly(np.full(3, np.nan)), readonly(np.full((3, 3), np.nan)), False, -1)
        nodes = self.tetra[cell]
        return FlowSample(readonly(weights @ self.velocity_nodes_m_s[nodes]),
                          float(weights @ self.pressure_nodes_pa[nodes]),
                          readonly(self.gradients_s_inv[cell]), readonly(self.curl_s_inv[cell]),
                          readonly(self.strain_s_inv[cell]), True, cell)

    def sample_many(self, positions_m):
        positions = np.asarray(positions_m, dtype=np.float64)
        if positions.ndim != 2 or positions.shape[1] != 3 or not np.isfinite(positions).all():
            raise ValueError("positions_m must be finite with shape (N,3) in m")
        samples = [self.sample(p) for p in positions]
        shapes = {"velocity_m_s": (3,), "pressure_pa": (), "velocity_gradient_s_inv": (3, 3),
                  "vorticity_s_inv": (3,), "strain_rate_s_inv": (3, 3), "inside_lumen": (), "tetra_id": ()}
        data = {}
        for f in fields(FlowBatch):
            dtype = bool if f.name == "inside_lumen" else np.int64 if f.name == "tetra_id" else np.float64
            values = [getattr(s, f.name) for s in samples]
            data[f.name] = readonly(np.asarray(values, dtype=dtype).reshape((len(samples), *shapes[f.name])), dtype)
        return FlowBatch(**data)
