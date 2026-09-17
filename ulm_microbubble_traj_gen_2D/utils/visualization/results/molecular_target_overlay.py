"""Load and prepare molecular-target walls for trajectory visualization."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import ndimage


@dataclass(frozen=True)
class MolecularTargetOverlay:
    """Static target-wall information aligned with the CFD visualization grid."""

    source_path: Path | None
    enabled: bool
    target_wall_mask: np.ndarray
    target_density_field_molecules_per_m2: np.ndarray
    target_density_molecules_per_m2: float

    @property
    def target_wall_site_count(self) -> int:
        """Return the number of target-positive wall cells."""

        return int(np.count_nonzero(self.target_wall_mask))

    @property
    def visible(self) -> bool:
        """Return whether the renderer has at least one target wall cell to draw."""

        return bool(self.enabled and self.target_wall_site_count > 0)


def load_molecular_target_overlay(
    result_dir: Path,
    x_coordinates_um: np.ndarray,
    z_coordinates_um: np.ndarray,
) -> MolecularTargetOverlay:
    """Load the optional target field and verify that it matches the CFD grid."""

    result_path = Path(result_dir)
    target_path = result_path / "molecular_target_field.npz"
    expected_x = np.asarray(x_coordinates_um, dtype=np.float64)
    expected_z = np.asarray(z_coordinates_um, dtype=np.float64)
    expected_shape = (expected_x.size, expected_z.size)

    if not target_path.is_file():
        return _empty_overlay(expected_shape)

    with np.load(target_path, allow_pickle=False) as target:
        stored_x = np.asarray(target["x_coordinates_um"], dtype=np.float64)
        stored_z = np.asarray(target["z_coordinates_um"], dtype=np.float64)
        if not np.array_equal(stored_x, expected_x) or not np.array_equal(stored_z, expected_z):
            raise ValueError(
                "The molecular target coordinates do not match the visualized CFD grid: "
                f"{target_path}"
            )

        target_wall_mask = np.asarray(target["target_wall_mask"], dtype=bool)
        target_density_field = np.asarray(
            target["target_density_field_molecules_per_m2"],
            dtype=np.float64,
        )
        if target_wall_mask.shape != expected_shape or target_density_field.shape != expected_shape:
            raise ValueError(
                "The molecular target arrays must match the visualized CFD grid shape "
                f"{expected_shape}: {target_path}"
            )

        enabled = _scalar_bool(target, "enabled", fallback=bool(np.any(target_wall_mask)))
        density = _scalar_float(target, "target_density_molecules_per_m2", fallback=0.0)

    return MolecularTargetOverlay(
        source_path=target_path,
        enabled=enabled,
        target_wall_mask=np.ascontiguousarray(target_wall_mask),
        target_density_field_molecules_per_m2=np.ascontiguousarray(target_density_field),
        target_density_molecules_per_m2=density,
    )


def target_display_masks(
    target_wall_mask: np.ndarray,
    *,
    halo_cells: int = 2,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the exact target mask and a separate thin visibility halo."""

    exact = np.asarray(target_wall_mask, dtype=bool)
    if not np.any(exact) or int(halo_cells) <= 0:
        return exact, np.zeros(exact.shape, dtype=bool)
    expanded = ndimage.binary_dilation(exact, iterations=int(halo_cells))
    return exact, np.asarray(expanded & ~exact, dtype=bool)


def _empty_overlay(shape: tuple[int, int]) -> MolecularTargetOverlay:
    """Create a harmless no-target overlay for legacy result directories."""

    return MolecularTargetOverlay(
        source_path=None,
        enabled=False,
        target_wall_mask=np.zeros(shape, dtype=bool),
        target_density_field_molecules_per_m2=np.zeros(shape, dtype=np.float64),
        target_density_molecules_per_m2=0.0,
    )


def _scalar_bool(target: np.lib.npyio.NpzFile, key: str, *, fallback: bool) -> bool:
    """Read one optional Boolean value from an NPZ file."""

    if key not in target:
        return bool(fallback)
    values = np.asarray(target[key]).reshape(-1)
    return bool(values[0]) if values.size else bool(fallback)


def _scalar_float(target: np.lib.npyio.NpzFile, key: str, *, fallback: float) -> float:
    """Read one optional floating-point value from an NPZ file."""

    if key not in target:
        return float(fallback)
    values = np.asarray(target[key], dtype=np.float64).reshape(-1)
    return float(values[0]) if values.size else float(fallback)
