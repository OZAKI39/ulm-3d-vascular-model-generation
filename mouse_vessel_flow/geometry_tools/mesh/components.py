"""Transparent selection and classification of disconnected mesh components."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True, slots=True)
class MainNetworkSelection:
    policy: str
    selected_component_id: int
    selection_method: str
    requested_component_id: int | None
    largest_surface_area_component_id: int
    largest_triangle_count_component_id: int
    largest_diagonal_component_id: int
    ranking_agrees: bool
    top_candidates: list[dict[str, Any]]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def select_main_network(
    *,
    policy: str,
    requested_component_id: int | None,
    triangle_counts: np.ndarray,
    surface_areas: np.ndarray,
    diagonals: np.ndarray,
    candidate_limit: int = 10,
) -> MainNetworkSelection:
    component_count = len(surface_areas)
    if component_count == 0:
        raise ValueError("No connected mesh component is available")
    if requested_component_id is not None and not 0 <= requested_component_id < component_count:
        raise ValueError(
            f"main_component_id={requested_component_id} is outside [0, {component_count - 1}]"
        )

    largest_area_id = int(np.argmax(surface_areas))
    largest_faces_id = int(np.argmax(triangle_counts))
    largest_diagonal_id = int(np.argmax(diagonals))
    if requested_component_id is None:
        selected_id = largest_area_id
        method = "largest surface area"
    else:
        selected_id = requested_component_id
        method = "explicit component ID"

    order = np.argsort(surface_areas)[::-1][:candidate_limit]
    top_candidates = [
        {
            "component_id": int(component_id),
            "surface_area_um2": float(surface_areas[component_id]),
            "triangle_count": int(triangle_counts[component_id]),
            "bbox_diagonal_um": float(diagonals[component_id]),
            "selected": int(component_id) == selected_id,
        }
        for component_id in order
    ]
    return MainNetworkSelection(
        policy=policy,
        selected_component_id=selected_id,
        selection_method=method,
        requested_component_id=requested_component_id,
        largest_surface_area_component_id=largest_area_id,
        largest_triangle_count_component_id=largest_faces_id,
        largest_diagonal_component_id=largest_diagonal_id,
        ranking_agrees=len({largest_area_id, largest_faces_id, largest_diagonal_id}) == 1,
        top_candidates=top_candidates,
    )


def is_small_fragment(
    *,
    triangle_count: int,
    surface_area_um2: float,
    diagonal_um: float,
    min_faces: int,
    min_area_um2: float,
    min_diagonal_um: float,
) -> bool:
    return (
        triangle_count < min_faces
        and surface_area_um2 < min_area_um2
        and diagonal_um < min_diagonal_um
    )

