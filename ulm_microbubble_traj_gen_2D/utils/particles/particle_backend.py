"""Acceleration-backend selection for particle transport."""

from __future__ import annotations

try:
    import numba as _numba
    from numba import njit as _njit
except ImportError:  # pragma: no cover - depends on the active environment
    _numba = None
    _njit = None

from ..geometry.continuous_vessel_numba import _EXACT_STATE_PARALLEL_MIN_BATCH


def resolve_particle_backend(requested: str) -> str:
    """Resolve the configured backend to the implementation used at runtime."""

    value = str(requested).lower()
    if value == "auto":
        return "numba_cpu" if _njit is not None else "python"
    if value == "numba_cpu":
        if _njit is None:
            raise ImportError(
                "Numba was requested but is not installed in the active Python environment."
            )
        return value
    if value == "python":
        return value
    raise ValueError(
        "particles.acceleration_backend must be one of: auto, numba_cpu, python."
    )


def configure_particle_numba_worker_threads(backend: str) -> int:
    """Use the measured small-batch worker count for parallel swept paths."""

    if str(backend).lower() != "numba_cpu" or _numba is None:
        return 1
    current = int(_numba.get_num_threads())
    workers = max(1, min(current, 4))
    if workers != current:
        _numba.set_num_threads(workers)
    return int(_numba.get_num_threads())


def particle_backend_details(backend: str) -> dict[str, object]:
    """Return additive reproducibility metadata for the resolved backend."""

    value = str(backend).lower()
    if value == "numba_cpu" and _numba is not None:
        return {
            "particle_numeric_kernel_family": "numba_batched_component_kernels_v18",
            "particle_numba_version": str(_numba.__version__),
            "particle_numba_thread_capacity": int(
                _numba.config.NUMBA_NUM_THREADS
            ),
            "particle_numba_worker_threads": int(_numba.get_num_threads()),
            "particle_numba_cache_enabled": True,
            "particle_numba_parallel_particle_loops": True,
            "particle_numba_parallel_swept_path_queries": True,
            "particle_numba_parallel_exact_wall_state_queries": True,
            "particle_numba_exact_wall_state_parallel_min_batch": int(
                _EXACT_STATE_PARALLEL_MIN_BATCH
            ),
            "particle_numba_fused_pulsatile_sampling": True,
            "particle_numba_predictive_contact_geometry": True,
            "particle_numba_directed_inlet_crossing_guard": True,
            "particle_continuous_wall_inward_normals_cached": True,
            "particle_numba_cached_predictive_wall_endpoints": True,
            "particle_numba_frame_transaction_batching": True,
            "particle_numba_outlet_spatial_index": "uniform_grid_csr",
            "particle_numba_exact_solid_face_queries": True,
            "particle_numba_continuous_wall_segment_queries": True,
            "particle_numba_swept_disc_audit": True,
            "particle_numba_scalar_diagnostic_reduction": True,
            "particle_taichi_backend_used": False,
        }
    return {
        "particle_numeric_kernel_family": "python_reference",
        "particle_numba_version": "not_used",
        "particle_numba_thread_capacity": 1,
        "particle_numba_worker_threads": 1,
        "particle_numba_cache_enabled": False,
        "particle_numba_parallel_particle_loops": False,
        "particle_numba_parallel_swept_path_queries": False,
        "particle_numba_parallel_exact_wall_state_queries": False,
        "particle_numba_exact_wall_state_parallel_min_batch": 0,
        "particle_numba_fused_pulsatile_sampling": False,
        "particle_numba_predictive_contact_geometry": False,
        "particle_numba_directed_inlet_crossing_guard": False,
        "particle_continuous_wall_inward_normals_cached": False,
        "particle_numba_cached_predictive_wall_endpoints": False,
        "particle_numba_frame_transaction_batching": False,
        "particle_numba_outlet_spatial_index": "not_used",
        "particle_numba_exact_solid_face_queries": False,
        "particle_numba_continuous_wall_segment_queries": False,
        "particle_numba_swept_disc_audit": False,
        "particle_numba_scalar_diagnostic_reduction": False,
        "particle_taichi_backend_used": False,
    }
