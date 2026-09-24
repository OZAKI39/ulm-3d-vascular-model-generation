"""Adapters for the existing representative-ROI sampling output."""

from __future__ import annotations

from pathlib import Path

from utils.sampling.sampling_io import load_sampling_display_rois
from utils.sampling.sampling_types import ROIRecord


def _is_sampling_run(path: Path) -> bool:
    return (
        (path / "manifests" / "candidate_rois.csv").is_file()
        and (path / "manifests" / "selected_rois.csv").is_file()
        and (path / "roi_library").is_dir()
    )


def resolve_sampling_run(path: Path | None, *, project_root: Path) -> Path:
    """Resolve an explicit run or the latest complete run under outputs/sampling."""

    if path is not None:
        resolved = Path(path).resolve()
        if not _is_sampling_run(resolved):
            raise FileNotFoundError(f"Not a complete sampling run: {resolved}")
        return resolved
    base = Path(project_root).resolve() / "outputs" / "sampling"
    candidates = sorted(
        (candidate for candidate in base.iterdir() if candidate.is_dir() and _is_sampling_run(candidate)),
        key=lambda candidate: (candidate.stat().st_mtime_ns, candidate.name),
        reverse=True,
    ) if base.is_dir() else []
    if not candidates:
        raise FileNotFoundError(
            f"No complete sampling run found under {base}; provide --sampling-run explicitly"
        )
    return candidates[0].resolve()


def load_sampling_rois(
    run_root: Path,
    *,
    roi_id: str | None = None,
    selected_only: bool = True,
) -> list[ROIRecord]:
    """Load geometry through the sampling project's existing NPZ/manifest adapter."""

    root = Path(run_root).resolve()
    if roi_id is not None:
        rois = load_sampling_display_rois(root, selected_only=False)
        matches = [roi for roi in rois if roi.roi_id == roi_id]
        if not matches:
            available = ", ".join(roi.roi_id for roi in rois[:5])
            raise KeyError(f"ROI {roi_id!r} is absent from {root}; first available IDs: {available}")
        return matches
    rois = load_sampling_display_rois(root, selected_only=selected_only)
    if not rois:
        kind = "selected representative" if selected_only else "candidate"
        raise ValueError(f"Sampling run contains no loadable {kind} ROIs: {root}")
    return sorted(rois, key=lambda roi: (roi.selection_rank if roi.selection_rank >= 0 else 10**9, roi.roi_id))
