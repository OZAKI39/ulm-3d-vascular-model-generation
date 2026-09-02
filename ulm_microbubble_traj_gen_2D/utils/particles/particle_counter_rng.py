"""Order-independent counter-based random numbers for particle transport."""

from __future__ import annotations

import math

import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover - production installs include Numba
    njit = None


COUNTER_RNG_ALGORITHM = "splitmix64_box_muller_v1"
_MASK64 = (1 << 64) - 1
_STREAM_1 = 0xD2B74407B1CE6E93
_STREAM_2 = 0xCA5A826395121157
_TWO_POW_MINUS_53 = 1.0 / float(1 << 53)


def _splitmix64_python(value: int) -> int:
    value = (int(value) + 0x9E3779B97F4A7C15) & _MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & _MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & _MASK64
    return (value ^ (value >> 31)) & _MASK64


def counter_hashes(
    seed: int,
    permanent_particle_id: int,
    global_internal_step: int,
    component: int = 0,
) -> tuple[int, int]:
    """Return the two fixed stream hashes for one immutable particle key."""

    values = (seed, permanent_particle_id, global_internal_step, component)
    if any(isinstance(value, bool) or int(value) != value or value < 0 for value in values):
        raise ValueError("Counter-RNG key values must be non-negative integers.")
    state = _splitmix64_python(int(seed))
    state = _splitmix64_python(state ^ int(permanent_particle_id))
    state = _splitmix64_python(state ^ int(global_internal_step))
    state = _splitmix64_python(state ^ int(component))
    return (
        _splitmix64_python(state ^ _STREAM_1),
        _splitmix64_python(state ^ _STREAM_2),
    )


def counter_normal(
    seed: int,
    permanent_particle_id: int,
    global_internal_step: int,
    component: int = 0,
) -> float:
    """Return one standard normal variate without mutable RNG state."""

    first, second = counter_hashes(
        seed, permanent_particle_id, global_internal_step, component
    )
    uniform_1 = ((first >> 11) + 0.5) * _TWO_POW_MINUS_53
    uniform_2 = ((second >> 11) + 0.5) * _TWO_POW_MINUS_53
    return math.sqrt(-2.0 * math.log(uniform_1)) * math.cos(
        2.0 * math.pi * uniform_2
    )


def counter_normal_batch(
    seed: int,
    permanent_particle_id: np.ndarray,
    global_internal_step: int,
    component: int = 0,
    *,
    use_numba: bool,
) -> np.ndarray:
    """Vector form whose values depend on permanent IDs, never lane order."""

    ids = np.ascontiguousarray(permanent_particle_id, dtype=np.int64)
    if int(seed) < 0 or int(global_internal_step) < 0 or int(component) < 0:
        raise ValueError("Counter-RNG key values must be non-negative integers.")
    if np.any(ids < 0):
        raise ValueError("Permanent particle IDs must be non-negative.")
    kernel = (
        _counter_normal_batch_numba
        if use_numba and njit is not None
        else None
    )
    if kernel is None:
        return _counter_normal_batch_numpy(
            int(seed), ids, int(global_internal_step), int(component)
        )
    return kernel(
        np.uint64(seed),
        ids,
        np.uint64(global_internal_step),
        np.uint64(component),
    )


def _splitmix64_numpy(value: np.ndarray) -> np.ndarray:
    with np.errstate(over="ignore"):
        value = value + np.uint64(0x9E3779B97F4A7C15)
        value = (value ^ (value >> np.uint64(30))) * np.uint64(
            0xBF58476D1CE4E5B9
        )
        value = (value ^ (value >> np.uint64(27))) * np.uint64(
            0x94D049BB133111EB
        )
    return value ^ (value >> np.uint64(31))


def _counter_normal_batch_numpy(seed, permanent_ids, global_step, component):
    state = np.full(permanent_ids.size, np.uint64(seed), dtype=np.uint64)
    state = _splitmix64_numpy(state)
    state = _splitmix64_numpy(state ^ permanent_ids.astype(np.uint64))
    state = _splitmix64_numpy(state ^ np.uint64(global_step))
    state = _splitmix64_numpy(state ^ np.uint64(component))
    first = _splitmix64_numpy(state ^ np.uint64(_STREAM_1))
    second = _splitmix64_numpy(state ^ np.uint64(_STREAM_2))
    uniform_1 = ((first >> np.uint64(11)).astype(np.float64) + 0.5) * (
        _TWO_POW_MINUS_53
    )
    uniform_2 = ((second >> np.uint64(11)).astype(np.float64) + 0.5) * (
        _TWO_POW_MINUS_53
    )
    return np.sqrt(-2.0 * np.log(uniform_1)) * np.cos(2.0 * np.pi * uniform_2)


def _splitmix64_kernel(value):
    value = value + np.uint64(0x9E3779B97F4A7C15)
    value = (value ^ (value >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    value = (value ^ (value >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return value ^ (value >> np.uint64(31))


def _counter_normal_batch_kernel(seed, permanent_ids, global_step, component):
    count = permanent_ids.size
    result = np.empty(count, dtype=np.float64)
    for lane in range(count):
        state = _splitmix64_kernel(seed)
        state = _splitmix64_kernel(state ^ np.uint64(permanent_ids[lane]))
        state = _splitmix64_kernel(state ^ global_step)
        state = _splitmix64_kernel(state ^ component)
        first = _splitmix64_kernel(state ^ np.uint64(_STREAM_1))
        second = _splitmix64_kernel(state ^ np.uint64(_STREAM_2))
        uniform_1 = (
            float(first >> np.uint64(11)) + 0.5
        ) * _TWO_POW_MINUS_53
        uniform_2 = (
            float(second >> np.uint64(11)) + 0.5
        ) * _TWO_POW_MINUS_53
        result[lane] = math.sqrt(-2.0 * math.log(uniform_1)) * math.cos(
            2.0 * math.pi * uniform_2
        )
    return result


if njit is not None:
    _splitmix64_kernel = njit(cache=True, nogil=True)(_splitmix64_kernel)
    _counter_normal_batch_numba = njit(cache=True, nogil=True)(
        _counter_normal_batch_kernel
    )
else:  # pragma: no cover
    _counter_normal_batch_numba = None
