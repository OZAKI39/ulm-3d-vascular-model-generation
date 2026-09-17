"""Strict entry point for volumetric vascular models.

The current production transport data model is intentionally two-dimensional:
its grid, finite-element cache, continuous wall, molecular target, particle
contact, and renderers all live in X-Z.  A volumetric model must therefore
never be projected into that workflow implicitly.
"""

from __future__ import annotations

from ..io.vascular_io import load_physics_input, validate_physics_input_geometry


class VolumetricPipelineUnavailableError(RuntimeError):
    """Raised when a real 3-D model is selected without a complete 3-D backend."""


def run_volumetric_generation(
    cfg,
    *,
    render_artifacts: bool = True,
    reuse_field_from=None,
):
    """Validate a 3-D input and stop before any planar operation is performed."""

    physics_input = load_physics_input(cfg.model_dir)
    validate_physics_input_geometry(
        cfg.model_dir,
        physics_input,
        expected_mode="volumetric_3d",
    )
    raise VolumetricPipelineUnavailableError(
        "The selected vascular result is volumetric_3d, but this repository "
        "does not yet contain a complete three-dimensional CFD and particle "
        "transport backend. A valid 3-D path requires a volumetric lumen mesh, "
        "a three-component 3-D velocity/pressure solve, 3-D wall contact and "
        "molecular-target geometry, 3-D particle advection, field reuse, and "
        "3-D output/rendering contracts. The model was rejected before X-Z "
        "domain construction; it was not projected into the planar solver. "
        f"Selected model: {cfg.model_dir}"
    )
