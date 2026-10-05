"""Complex SI pressure fields with bounded trilinear interpolation and HDF5 storage."""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import RegularGridInterpolator


@dataclass
class ComplexAcousticField:
    axes_m: tuple[np.ndarray, ...]
    pressure_Pa: np.ndarray
    gradient_Pa_m: np.ndarray
    frequency_Hz: float

    def __post_init__(self):
        self._p = RegularGridInterpolator(self.axes_m, self.pressure_Pa, bounds_error=True)
        self._g = RegularGridInterpolator(self.axes_m, self.gradient_Pa_m, bounds_error=True)

    def interpolate_complex_pressure(self, position_m: np.ndarray) -> np.ndarray:
        return self._p(position_m)

    def interpolate_complex_gradient(self, position_m: np.ndarray) -> np.ndarray:
        return self._g(position_m)


class AcousticBasisDatabase:
    """H is pressure per applied source Pa; control entries carry Pa."""

    def __init__(self, axes_m, pressure_basis, gradient_basis, metadata):
        if metadata.get("BonA") is not None:
            raise ValueError("Nonlinear propagation cannot be combined as a linear basis")
        self.axes_m = tuple(np.asarray(a) for a in axes_m)
        self.pressure_basis = np.asarray(pressure_basis)
        self.gradient_basis = np.asarray(gradient_basis)
        self.metadata = metadata
        self.field = None

    def combine(self, control_vector: np.ndarray) -> ComplexAcousticField:
        u = np.asarray(control_vector)
        if u.shape != (len(self.pressure_basis),):
            raise ValueError("Control vector must match source count")
        p = np.einsum("i,ixyz->xyz", u, self.pressure_basis)
        g = np.einsum("i,ixyzj->xyzj", u, self.gradient_basis)
        self.field = ComplexAcousticField(self.axes_m, p, g, self.metadata["frequency_Hz"])
        return self.field

    def interpolate_complex_pressure(self, position_m):
        if self.field is None:
            raise RuntimeError("Call combine first")
        return self.field.interpolate_complex_pressure(position_m)

    def interpolate_complex_gradient(self, position_m):
        if self.field is None:
            raise RuntimeError("Call combine first")
        return self.field.interpolate_complex_gradient(position_m)

    def save_basis(self, path: str | Path) -> None:
        with h5py.File(path, "x") as f:
            f.attrs["metadata_json"] = json.dumps(self.metadata)
            f.attrs["harmonic_convention"] = "Re[P exp(+i omega t)]"
            f.attrs["basis_units"] = "Pa / source_Pa; gradient 1/m"
            for axis, values in zip("xyz", self.axes_m, strict=True):
                f.create_dataset(f"axes/{axis}_m", data=values)
            f.create_dataset("P_complex", data=self.pressure_basis.astype(np.complex64), compression="gzip", compression_opts=1)
            for j, axis in enumerate("xyz"):
                f.create_dataset(f"gradP_complex_{axis}", data=self.gradient_basis[..., j].astype(np.complex64), compression="gzip", compression_opts=1)

    @classmethod
    def load_basis(cls, path: str | Path):
        with h5py.File(path) as f:
            return cls(tuple(f[f"axes/{a}_m"][:] for a in "xyz"), f["P_complex"][:],
                       np.stack([f[f"gradP_complex_{a}"][:] for a in "xyz"], axis=-1), json.loads(f.attrs["metadata_json"]))
