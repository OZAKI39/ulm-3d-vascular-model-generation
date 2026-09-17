"""Shared runtime progress helpers for flow and particle calculations."""

from __future__ import annotations

import sys
import threading
from typing import Any

from .console_output import PROGRESS_BAR_FORMAT


def create_progress_bar(
    total: int,
    *,
    description: str,
    unit: str,
    position: int = 0,
    leave: bool = True,
) -> Any:
    """Create a consistently styled ``tqdm`` bar or its no-op fallback."""

    try:
        from tqdm.auto import tqdm

        return tqdm(
            total=max(0, int(total)),
            desc=description,
            unit=unit,
            file=sys.stdout,
            ascii=True,
            dynamic_ncols=True,
            bar_format=PROGRESS_BAR_FORMAT,
            position=int(position),
            leave=bool(leave),
        )
    except ImportError:  # pragma: no cover - only used in minimal environments
        return _NoProgress()


def create_particle_progress_bar(total_internal_steps: int) -> Any:
    """Create the sole Stage 06 bar using physical integration steps."""

    return create_progress_bar(
        total_internal_steps,
        description="Stage 06 simulation",
        unit="substep",
        position=0,
        leave=True,
    )


def create_stage_progress_bar(stage: int, total_modules: int) -> Any:
    """Create a stage bar with a live bar for its currently running module."""

    return _StageProgress(
        create_progress_bar(
            total_modules,
            description=f"Stage {int(stage):02d} modules",
            unit="module",
            position=0,
            leave=True,
        )
    )


class _StageProgress:
    """Pair a stage-wide bar with a live, elapsed-time current-module bar."""

    def __init__(self, stage_bar: Any) -> None:
        self._stage_bar = stage_bar
        self._current_module: _LiveModuleProgress | None = None
        self._closed = False

    def __enter__(self):
        return self

    def __exit__(self, _exc_type, _exc_value, _traceback) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def update(self, amount: int = 1) -> None:
        self._finish_current_module(completed=True)
        self._stage_bar.update(amount)

    def set_postfix(self, **kwargs: object) -> None:
        module = kwargs.get("module")
        self._finish_current_module(completed=False)
        self._stage_bar.set_postfix(**kwargs)
        if module is not None:
            self._current_module = _LiveModuleProgress(str(module))

    def set_submodule_progress(
        self,
        *,
        completed: int,
        total: int,
        submodule: str,
    ) -> None:
        """Show quantitative checkpoints inside the active module."""

        if self._current_module is not None:
            self._current_module.set_progress(
                completed=completed,
                total=total,
                submodule=submodule,
            )

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        self._finish_current_module(completed=False)
        self._stage_bar.close()

    def _finish_current_module(self, *, completed: bool) -> None:
        if self._current_module is None:
            return
        self._current_module.close(completed=completed)
        self._current_module = None


class _LiveModuleProgress:
    """Refresh a one-task bar so elapsed time changes during blocking work."""

    _REFRESH_INTERVAL_S = 1.0

    def __init__(self, module: str) -> None:
        self._total = 1
        self._completed = 0
        self._bar = create_progress_bar(
            self._total,
            description=f"  Current: {module}",
            unit="step",
            position=1,
            leave=False,
        )
        self._stop = threading.Event()
        self._closed = False
        self._thread = threading.Thread(
            target=_refresh_progress_until_stopped,
            args=(self._bar, self._stop, self._REFRESH_INTERVAL_S),
            name="ulm-current-module-progress",
            daemon=True,
        )
        self._thread.start()

    def set_progress(
        self,
        *,
        completed: int,
        total: int,
        submodule: str,
    ) -> None:
        """Update planned checkpoints and identify the active inner operation."""

        requested_total = max(1, int(total))
        requested_completed = min(
            max(0, int(completed)),
            requested_total,
        )
        if requested_total != self._total:
            self._total = requested_total
            self._bar.total = requested_total
        delta = requested_completed - self._completed
        if delta > 0:
            self._bar.update(delta)
        self._completed = max(self._completed, requested_completed)
        self._bar.set_postfix(submodule=str(submodule), refresh=True)

    def close(self, *, completed: bool) -> None:
        if self._closed:
            return
        self._closed = True
        self._stop.set()
        self._thread.join(timeout=self._REFRESH_INTERVAL_S)
        if completed:
            remaining = max(0, self._total - self._completed)
            if remaining:
                self._bar.update(remaining)
            self._completed = self._total
        self._bar.close()

    def __del__(self) -> None:
        try:
            self.close(completed=False)
        except Exception:
            pass


def _refresh_progress_until_stopped(
    progress: Any,
    stop: threading.Event,
    interval_s: float,
) -> None:
    """Refresh elapsed time while a long-running module owns the child bar."""

    while not stop.wait(float(interval_s)):
        try:
            progress.refresh()
        except Exception:
            return


class _NoProgress:
    """Subset of the tqdm API used by both numerical drivers."""

    def __enter__(self):
        return self

    def __exit__(self, _exc_type, _exc_value, _traceback) -> None:
        self.close()

    def update(self, _amount: int = 1) -> None:
        return None

    def set_postfix(self, **_kwargs: object) -> None:
        return None

    def refresh(self) -> None:
        return None

    def close(self) -> None:
        return None
