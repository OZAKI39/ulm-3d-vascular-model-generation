"""Deterministic synthetic molecular-target selection from vessel-bed candidates."""

from __future__ import annotations

from dataclasses import dataclass
import math

from .molecular_target_candidates import (
    MolecularTargetCandidate,
    MolecularTargetCandidateCatalog,
)


MINIMUM_EXPECTED_BUBBLE_VISITS = 1.0
AUTOMATIC_CANDIDATE_KIND = "downstream_subtree"


@dataclass(frozen=True)
class AutomaticTargetCandidateEvaluation:
    """Eligibility and ranking information for one catalog candidate."""

    candidate_id: str
    eligible: bool
    rejection_reason: str | None
    endothelial_wall_area_fraction: float
    expected_bubble_visits: float
    area_log_error: float


@dataclass(frozen=True)
class AutomaticInfluenceAnchorResult:
    """One accessible downstream bed used only to locate a tissue-space region."""

    requested_influence_wall_area_fraction: float
    anchor_candidate_id: str
    anchor_x_um: float
    anchor_z_um: float
    eligible_candidate_count: int
    evaluations: tuple[AutomaticTargetCandidateEvaluation, ...]


def select_automatic_influence_anchor(
    catalog: MolecularTargetCandidateCatalog,
    influence_wall_area_fraction: float,
) -> AutomaticInfluenceAnchorResult:
    """Choose a reproducible coarse anchor without making its subtree target-positive."""

    target_fraction = float(influence_wall_area_fraction)
    if not math.isfinite(target_fraction) or not 0.0 < target_fraction < 1.0:
        raise ValueError(
            "influence_wall_area_fraction must be finite and strictly between 0 and 1."
        )
    if not catalog.automatic_metrics_available:
        raise ValueError(
            "Revised-v10 automatic target generation requires a v3 candidate catalog with "
            "physical wall-area weights and accessibility. Rebuild the candidate catalog."
        )
    evaluations = tuple(
        _evaluate_candidate(candidate, target_fraction)
        for candidate in catalog.candidates
    )
    evaluation_by_id = {
        evaluation.candidate_id: evaluation for evaluation in evaluations
    }
    eligible = [
        candidate
        for candidate in catalog.candidates
        if evaluation_by_id[candidate.candidate_id].eligible
    ]
    if not eligible:
        raise ValueError(
            "No accessible downstream-subtree candidate can anchor the influence region."
        )
    anchor = min(
        eligible,
        key=lambda candidate: _individual_rank(candidate, target_fraction),
    )
    return AutomaticInfluenceAnchorResult(
        requested_influence_wall_area_fraction=target_fraction,
        anchor_candidate_id=anchor.candidate_id,
        anchor_x_um=anchor.wall_area_centroid_x_um,
        anchor_z_um=anchor.wall_area_centroid_z_um,
        eligible_candidate_count=len(eligible),
        evaluations=evaluations,
    )


def _evaluate_candidate(
    candidate: MolecularTargetCandidate,
    target_fraction: float,
) -> AutomaticTargetCandidateEvaluation:
    rejection_reason: str | None = None
    if candidate.kind != AUTOMATIC_CANDIDATE_KIND:
        rejection_reason = "not_a_downstream_subtree"
    elif not math.isfinite(candidate.endothelial_wall_area_fraction) or (
        candidate.endothelial_wall_area_fraction <= 0.0
    ):
        rejection_reason = "invalid_endothelial_wall_area"
    elif candidate.topology_depth < 0 or not math.isfinite(
        candidate.radius_of_gyration_um
    ) or candidate.radius_of_gyration_um < 0.0:
        rejection_reason = "invalid_structural_metric"
    elif not math.isfinite(candidate.network_flow_fraction) or (
        candidate.network_flow_fraction <= 0.0
    ):
        rejection_reason = "not_perfused"
    elif not math.isfinite(candidate.expected_bubble_visits):
        rejection_reason = "missing_accessibility_estimate"
    elif candidate.expected_bubble_visits < MINIMUM_EXPECTED_BUBBLE_VISITS:
        rejection_reason = "expected_bubble_visits_below_one"

    area_error = (
        _area_log_error(candidate.endothelial_wall_area_fraction, target_fraction)
        if math.isfinite(candidate.endothelial_wall_area_fraction)
        and candidate.endothelial_wall_area_fraction > 0.0
        else math.inf
    )
    return AutomaticTargetCandidateEvaluation(
        candidate_id=candidate.candidate_id,
        eligible=rejection_reason is None,
        rejection_reason=rejection_reason,
        endothelial_wall_area_fraction=float(candidate.endothelial_wall_area_fraction),
        expected_bubble_visits=float(candidate.expected_bubble_visits),
        area_log_error=float(area_error),
    )


def _individual_rank(
    candidate: MolecularTargetCandidate,
    target_fraction: float,
) -> tuple[float, int, float, str]:
    return (
        _area_log_error(candidate.endothelial_wall_area_fraction, target_fraction),
        -candidate.topology_depth,
        candidate.radius_of_gyration_um,
        candidate.candidate_id,
    )


def _area_log_error(candidate_fraction: float, target_fraction: float) -> float:
    return abs(math.log(float(candidate_fraction) / float(target_fraction)))
