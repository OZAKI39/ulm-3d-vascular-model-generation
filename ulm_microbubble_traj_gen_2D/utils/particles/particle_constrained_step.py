"""Side-effect-free state and physical-time helpers for constrained stepping.

The continuous-perfusion driver owns the full Euler/Heun algorithm, particle
lifecycle, molecular-state updates, and output recording.  This module keeps
the transactional batch state and physical-time refinement rules independent of
that driver. Revised-v15 predictive contact is derived afresh from the exact
true wall gap, and its single-wall solve lives in a dedicated module.

None of the helpers shortens a requested displacement while pretending that
the original time has passed.  A failed constrained trial must instead be
retried over the returned physical half intervals.  Reaching the configured
refinement limit raises an error; silently holding particles in place would
reintroduce the numerical wall lock that this module is intended to remove.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


class PhysicalTimeRefinementError(RuntimeError):
    """Raised when a failed physical step cannot be divided safely again.

    ``reason`` is intentionally machine-readable so callers can distinguish a
    configured depth limit from the independent floating-point time floor.
    The optional numeric fields preserve backwards compatibility with callers
    that construct this exception from only a message while allowing production
    diagnostics to report the rejected interval precisely.
    """

    def __init__(
        self,
        message: str,
        *,
        reason: str = "unknown",
        duration_s: float | None = None,
        local_time_ulp_s: float | None = None,
        refinement_depth: int | None = None,
        maximum_refinement_depth: int | None = None,
    ) -> None:
        super().__init__(message)
        self.reason = str(reason)
        self.duration_s = None if duration_s is None else float(duration_s)
        self.local_time_ulp_s = (
            None if local_time_ulp_s is None else float(local_time_ulp_s)
        )
        self.refinement_depth = (
            None if refinement_depth is None else int(refinement_depth)
        )
        self.maximum_refinement_depth = (
            None if maximum_refinement_depth is None else int(maximum_refinement_depth)
        )


@dataclass(frozen=True)
class BatchLocalState:
    """A temporary active-particle state used by one trial time interval.

    The arrays intentionally mirror the four values changed by particle
    integration.  Permanent IDs and lifecycle flags remain owned by the
    perfusion driver.  Keeping the temporary state separate lets the driver
    discard a failed trial without accidentally committing position, rotation,
    or molecular-bond changes.

    The time integrator treats these arrays as immutable inputs and creates new
    arrays for every accepted result.  This keeps rejected trials from changing
    permanent particle state without paying for an unused preliminary copy.
    """

    position_grid: np.ndarray
    rotation_angle_rad: np.ndarray
    bond_count_expected: np.ndarray | None = None
    bond_total_tangential_extension_um: np.ndarray | None = None


@dataclass(frozen=True)
class PhysicalTimeInterval:
    """One half-open physical interval and its current refinement depth."""

    start_time_s: float
    end_time_s: float
    refinement_depth: int = 0

    def __post_init__(self) -> None:
        """Reject intervals that cannot represent forward physical time."""

        start = float(self.start_time_s)
        end = float(self.end_time_s)
        depth = int(self.refinement_depth)
        if not math.isfinite(start) or not math.isfinite(end):
            raise ValueError("Physical interval endpoints must be finite.")
        if end <= start:
            raise ValueError("A physical time interval must have positive duration.")
        if depth != self.refinement_depth or depth < 0:
            raise ValueError("refinement_depth must be a non-negative integer.")

    @property
    def duration_s(self) -> float:
        """Return the physical duration covered by this interval."""

        return float(self.end_time_s) - float(self.start_time_s)


def split_physical_time_interval(
    interval: PhysicalTimeInterval,
    maximum_refinement_depth: int,
) -> tuple[PhysicalTimeInterval, PhysicalTimeInterval]:
    """Bisect one failed interval without losing, duplicating, or faking time.

    Both returned intervals have depth ``parent.depth + 1`` and share the same
    midpoint object value.  Their outer endpoints are copied directly from the
    parent, which gives exact endpoint coverage with no gap or overlap.  The
    maximum depth is checked *before* constructing children so exhaustion is a
    clear numerical failure rather than an implicit zero-displacement step.
    """

    if not isinstance(interval, PhysicalTimeInterval):
        raise TypeError("interval must be a PhysicalTimeInterval.")
    maximum_depth = int(maximum_refinement_depth)
    if maximum_depth != maximum_refinement_depth or maximum_depth < 0:
        raise ValueError("maximum_refinement_depth must be a non-negative integer.")
    duration = float(interval.duration_s)
    local_time_ulp = max(
        math.ulp(float(interval.start_time_s)),
        math.ulp(float(interval.end_time_s)),
    )
    if interval.refinement_depth >= maximum_depth:
        raise PhysicalTimeRefinementError(
            "Wall-contact physical-time refinement reached "
            f"maximum_refinement_depth={maximum_depth}; the failed step was "
            "not accepted and particles were not silently held in place. "
            f"reason=depth_limit; duration_s={duration:.17g}; "
            f"local_time_ulp_s={local_time_ulp:.17g}; "
            f"refinement_depth={interval.refinement_depth}.",
            reason="depth_limit",
            duration_s=duration,
            local_time_ulp_s=local_time_ulp,
            refinement_depth=int(interval.refinement_depth),
            maximum_refinement_depth=maximum_depth,
        )

    # This midpoint form avoids overflow from (start + end) / 2.  Extremely
    # short intervals can still run out of representable floating-point times;
    # that condition must fail explicitly rather than create a zero-time child.
    midpoint = float(interval.start_time_s) + 0.5 * duration
    if not (float(interval.start_time_s) < midpoint < float(interval.end_time_s)):
        raise PhysicalTimeRefinementError(
            "The physical interval is too small to bisect in floating-point time. "
            f"reason=unrepresentable_midpoint; duration_s={duration:.17g}; "
            f"local_time_ulp_s={local_time_ulp:.17g}; "
            f"refinement_depth={interval.refinement_depth}; "
            f"maximum_refinement_depth={maximum_depth}.",
            reason="unrepresentable_midpoint",
            duration_s=duration,
            local_time_ulp_s=local_time_ulp,
            refinement_depth=int(interval.refinement_depth),
            maximum_refinement_depth=maximum_depth,
        )

    child_depth = int(interval.refinement_depth) + 1
    first = PhysicalTimeInterval(
        start_time_s=float(interval.start_time_s),
        end_time_s=midpoint,
        refinement_depth=child_depth,
    )
    second = PhysicalTimeInterval(
        start_time_s=midpoint,
        end_time_s=float(interval.end_time_s),
        refinement_depth=child_depth,
    )
    validate_physical_time_partition(interval, first, second)
    return first, second


def validate_physical_time_partition(
    parent: PhysicalTimeInterval,
    first: PhysicalTimeInterval,
    second: PhysicalTimeInterval,
) -> None:
    """Verify that two children exactly tile their parent in physical time."""

    expected_depth = int(parent.refinement_depth) + 1
    if (
        first.start_time_s != parent.start_time_s
        or first.end_time_s != second.start_time_s
        or second.end_time_s != parent.end_time_s
        or first.refinement_depth != expected_depth
        or second.refinement_depth != expected_depth
    ):
        raise ValueError(
            "Refined physical intervals must exactly cover the parent without "
            "a time gap, overlap, or depth mismatch."
        )
