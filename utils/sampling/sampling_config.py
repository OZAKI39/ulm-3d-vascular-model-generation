"""Configuration for reproducible real-ROI sampling."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True, slots=True)
class SamplingConfig:
    output_root: Path
    seed: int = 42
    anchor_mode: str = "farthest_point"
    min_anchor_distance_um: float = 45.0
    max_candidate_anchors: int = 80
    roi_size_um: tuple[float, float, float] = (80.0, 80.0, 120.0)
    min_branch_count: int = 2
    max_cut_ports: int | None = None
    feature_mode: str = "radius_plus_structure"
    radius_quantiles: tuple[float, ...] = (0.10, 0.25, 0.50, 0.75, 0.90)
    scaler: str = "robust"
    radius_feature_weight: float = 1.0
    structure_feature_weight: float = 1.0
    clustering_method: str = "kmeans"
    n_clusters: int = 5
    exploratory_k: tuple[int, ...] = (2, 3, 4, 5, 8, 10)
    kmeans_max_iter: int = 200
    selection_mode: str = "coverage_balanced"
    target_selected_count: int = 10
    representatives_per_cluster: int = 2
    max_selected_overlap: float = 0.25
    min_representative_distance_um: float = 20.0
    compare_feature_modes: bool = True
    max_roi_previews: int = 12

    def validate(self) -> None:
        if self.anchor_mode not in {"random", "farthest_point", "poisson_disk"}:
            raise ValueError(f"Unsupported anchor_mode: {self.anchor_mode}")
        if self.feature_mode not in {
            "radius_only",
            "radius_plus_structure",
            "extended_morphology",
        }:
            raise ValueError(f"Unsupported feature_mode: {self.feature_mode}")
        if self.scaler != "robust":
            raise ValueError("Only the reproducible robust scaler is supported")
        if self.clustering_method != "kmeans":
            raise ValueError("Only dependency-free deterministic kmeans is currently supported")
        if self.selection_mode not in {"distribution_preserving", "coverage_balanced"}:
            raise ValueError(f"Unsupported selection_mode: {self.selection_mode}")
        if any(value <= 0 for value in self.roi_size_um):
            raise ValueError("roi_size_um values must be positive")
        if self.min_anchor_distance_um <= 0:
            raise ValueError("min_anchor_distance_um must be positive")
        if self.max_candidate_anchors < 1:
            raise ValueError("max_candidate_anchors must be positive")
        if self.min_branch_count < 1:
            raise ValueError("min_branch_count must be positive")
        if self.n_clusters < 1:
            raise ValueError("n_clusters must be positive")
        if self.target_selected_count < 1 or self.representatives_per_cluster < 1:
            raise ValueError("selection counts must be positive")
        if not 0.0 <= self.max_selected_overlap <= 1.0:
            raise ValueError("max_selected_overlap must be in [0, 1]")
        if self.min_representative_distance_um < 0:
            raise ValueError("min_representative_distance_um cannot be negative")
        if self.max_cut_ports is not None and self.max_cut_ports < 0:
            raise ValueError("max_cut_ports cannot be negative")
        if any(not 0.0 <= value <= 1.0 for value in self.radius_quantiles):
            raise ValueError("radius_quantiles must be in [0, 1]")
        if self.radius_feature_weight <= 0 or self.structure_feature_weight <= 0:
            raise ValueError("feature group weights must be positive")
        if self.kmeans_max_iter < 1:
            raise ValueError("kmeans_max_iter must be positive")

    def report(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["output_root"] = str(self.output_root.resolve())
        payload["scientific_scope"] = (
            "full real vascular network -> representative real ROI library"
        )
        payload["physiological_flow_inference"] = False
        return payload
