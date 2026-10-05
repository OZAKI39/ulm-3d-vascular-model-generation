"""Vectorized, bounded complex response lookup (NumPy / CuPy).

Interpolation is in SI radius, frequency and amplitude; never in wrapped phase.
Cell validation is empirical midpoint evidence, not a certified error bound.
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np

from .response import BubbleResponseModel


class ResponseTableBoundsError(ValueError):
    pass


class ResponseTableBubbleModel(BubbleResponseModel):
    def __init__(self, path: str | Path, *, backend="numpy", require_validated_cells=True):
        if backend not in ("numpy", "cupy"):
            raise ValueError("Backend must be numpy or cupy")
        if backend == "cupy":
            import cupy as xp
        else:
            xp = np
        self.xp = xp
        self.path = Path(path)
        self.require_validated_cells = require_validated_cells
        with h5py.File(path) as h5:
            if h5.attrs.get("schema_version") != 1 or not h5.attrs.get("complete", False):
                raise ValueError("Incomplete or unsupported response table")
            self.metadata = json.loads(h5.attrs["metadata_json"])
            axes = [h5[k][:] for k in ("radius_m", "frequency_Hz", "pressure_Pa")]
            for a in axes:
                if len(a) < 2 or not np.isfinite(a).all() or np.any(np.diff(a) <= 0):
                    raise ValueError("Finite strictly increasing table axes required")
            shape = tuple(len(a) for a in axes)
            values = h5["harmonics_m3"][:]
            if values.shape[:3] != shape or values.ndim != 4 or not np.isfinite(values).all():
                raise ValueError("Invalid table shape or nonfinite responses")
            self.axes = tuple(xp.asarray(a) for a in axes)
            self.harmonics = xp.asarray(values)
            self.node_valid = xp.asarray(h5["node_valid"][:], dtype=bool)
            self.cell_valid = xp.asarray(h5["cell_valid"][:], dtype=bool)
            self.max_R = xp.asarray(h5["max_R_over_R0"][:])
            self.min_R = xp.asarray(h5["min_R_over_R0"][:])
            self.harmonic_resolved = xp.asarray(h5["harmonic_resolved"][:], dtype=bool)
            if (
                self.node_valid.shape != shape
                or self.harmonic_resolved.shape != shape
                or self.cell_valid.shape != tuple(n - 1 for n in shape)
            ):
                raise ValueError("Invalid response validity flag shape")
            if (
                self.max_R.shape != shape
                or self.min_R.shape != shape
                or bool(
                    xp.any(
                        ~xp.isfinite(self.max_R)
                        | ~xp.isfinite(self.min_R)
                        | (self.min_R <= 0)
                        | (self.max_R < self.min_R)
                    )
                )
            ):
                raise ValueError("Invalid radial response extrema")

    def _locate(self, R0, frequency, pressure):
        xp = self.xp
        q = xp.broadcast_arrays(
            xp.asarray(R0, dtype=float), xp.asarray(frequency, dtype=float), xp.asarray(pressure, dtype=float)
        )
        indices, weights = [], []
        for x, axis in zip(q, self.axes):
            if bool(xp.any(~xp.isfinite(x) | (x < axis[0]) | (x > axis[-1]))):
                raise ResponseTableBoundsError("Outside response table: extrapolation is forbidden")
            i = xp.clip(xp.searchsorted(axis, x, side="right") - 1, 0, len(axis) - 2)
            indices.append(i)
            weights.append((x - axis[i]) / (axis[i + 1] - axis[i]))
        exact = xp.ones(q[0].shape, dtype=bool)
        node_indices = []
        for i, w in zip(indices, weights):
            # Only actual axis values qualify; near a node is still interpolation.
            exact &= (w == 0) | (w == 1)
            node_indices.append(i + (w == 1).astype(int))
        return q, indices, weights, exact, node_indices

    def query(self, R0_m, frequency_Hz, pressure_amplitude_Pa):
        xp = self.xp
        q, idx, weights, exact, node_idx = self._locate(R0_m, frequency_Hz, pressure_amplitude_Pa)
        value = xp.zeros(q[0].shape + (self.harmonics.shape[-1],), dtype=complex)
        rmax, rmin = xp.zeros_like(q[0]), xp.zeros_like(q[0])
        envelope_max = xp.full_like(q[0], -xp.inf)
        envelope_min = xp.full_like(q[0], xp.inf)
        node_ok = xp.ones(q[0].shape, dtype=bool)
        resolved = node_ok.copy()
        for a in (0, 1):
            for b in (0, 1):
                for c in (0, 1):
                    corner = (idx[0] + a, idx[1] + b, idx[2] + c)
                    w = (
                        (weights[0] if a else 1 - weights[0])
                        * (weights[1] if b else 1 - weights[1])
                        * (weights[2] if c else 1 - weights[2])
                    )
                    value += w[..., None] * self.harmonics[corner]
                    rmax += w * self.max_R[corner]
                    rmin += w * self.min_R[corner]
                    envelope_max = xp.maximum(envelope_max, xp.where(w > 0, self.max_R[corner], -xp.inf))
                    envelope_min = xp.minimum(envelope_min, xp.where(w > 0, self.min_R[corner], xp.inf))
                    node_ok &= self.node_valid[corner] | (w == 0)
                    resolved &= self.harmonic_resolved[corner] | (w == 0)
        interpolation_ok = self.cell_valid[tuple(idx)] | exact
        if bool(xp.any(~node_ok)):
            raise ResponseTableBoundsError("Nonperiodic or compressibility-invalid radial response")
        if self.require_validated_cells and bool(xp.any(~interpolation_ok)):
            raise ResponseTableBoundsError("Unvalidated interpolation cell: refine table before production")
        return dict(
            harmonics_m3=value,
            max_R_over_R0=rmax,
            min_R_over_R0=rmin,
            corner_envelope_max_R_over_R0=envelope_max,
            corner_envelope_min_R_over_R0=envelope_min,
            harmonic_resolved=resolved,
            within_response_table=xp.ones(q[0].shape, bool),
            interpolation_validated=interpolation_ok,
            parameter_calibration_validated=False,
        )

    def evaluate_batch(self, R0_m, f, pressure_amplitude):
        return self.query(R0_m, f, pressure_amplitude)["harmonics_m3"][..., 0]

    def response(self, R0_m, frequency_Hz, local_complex_pressure_Pa, medium_properties):
        xp = self.xp
        self.validate_medium(medium_properties)
        p = xp.asarray(local_complex_pressure_Pa)
        return self.evaluate_batch(R0_m, frequency_Hz, xp.abs(p)) * xp.exp(1j * xp.angle(p))

    def validate_medium(self, medium_properties):
        # A table for water must not be silently used in blood.
        for key in ("density_kg_m3", "sound_speed_m_s", "viscosity_Pa_s"):
            value = medium_properties.get(key)
            if key == "viscosity_Pa_s":
                # Effective radial viscosity may include thermal losses; this
                # optional explicit value must not change translation drag.
                value = medium_properties.get("radial_effective_viscosity_Pa_s", value)
            if value is not None and not np.isclose(value, self.metadata["medium"][key], rtol=1e-10, atol=0):
                raise ValueError(f"Response table medium mismatch: {key}")

    def harmonic_response(self, R0_m, frequency_Hz, local_complex_pressure_Pa):
        xp = self.xp
        p = xp.asarray(local_complex_pressure_Pa)
        q = self.query(R0_m, frequency_Hz, xp.abs(p))
        h = xp.arange(1, self.harmonics.shape[-1] + 1)
        q["harmonics_m3"] *= xp.exp(1j * xp.angle(p)[..., None] * h)
        return q
