from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from utils.cfd_flow.adaptive_target_population import (
    AdaptiveTargetPopulationError,
    audit_adaptive_target_population,
)


MASS_FLOW = 2.8901803804796421e-12


@pytest.mark.parametrize("count", [1, 2, 375, 376, 1000])
def test_constant_target_is_count_invariant(count: int) -> None:
    result = audit_adaptive_target_population(
        [[MASS_FLOW] * count],
        [list(range(1, count + 1))],
        [[1] * count],
        stale_mesh_population=count,
    )
    assert result["active_count_global"] == count
    assert result["target_mass_flow"] == pytest.approx(MASS_FLOW, rel=2.0e-14)


def test_mpi_partition_reduces_same_active_population() -> None:
    result = audit_adaptive_target_population(
        [[], [MASS_FLOW] * 125, [MASS_FLOW] * 250, []],
        [[], list(range(1, 126)), list(range(126, 376)), []],
        [[], [4] * 125, [3] * 250, []],
        stale_mesh_population=376,
    )
    assert [item["active_count"] for item in result["per_rank"]] == [0, 125, 250, 0]
    assert result["active_count_global"] == 375
    assert result["sample_count_global"] == 375
    assert result["target_mass_flow"] == pytest.approx(MASS_FLOW, rel=2.0e-14)


def test_duplicate_positive_point_index_is_legal_and_counted() -> None:
    result = audit_adaptive_target_population(
        [[MASS_FLOW, MASS_FLOW]], [[7, 7]], [[1, 2]]
    )
    assert result["per_rank"][0]["point_duplicate_count"] == 1
    assert result["active_count_global"] == 2


@pytest.mark.parametrize("bad_index", [0, -1])
def test_invalid_point_index_fails_fast(bad_index: int) -> None:
    with pytest.raises(AdaptiveTargetPopulationError, match="invalid point index"):
        audit_adaptive_target_population([[MASS_FLOW]], [[bad_index]], [[1]])


def test_nonfinite_sample_fails_fast() -> None:
    with pytest.raises(AdaptiveTargetPopulationError, match="non-finite"):
        audit_adaptive_target_population([[np.nan]], [[1]], [[1]])


def test_inactive_element_in_active_list_fails_fast() -> None:
    with pytest.raises(AdaptiveTargetPopulationError, match="inactive element"):
        audit_adaptive_target_population([[MASS_FLOW]], [[1]], [[0]])


def test_population_length_mismatch_fails_fast() -> None:
    with pytest.raises(AdaptiveTargetPopulationError, match="population mismatch"):
        audit_adaptive_target_population([[MASS_FLOW]], [[1, 2]], [[1]])


def test_375_over_376_regression_is_removed() -> None:
    result = audit_adaptive_target_population(
        [[MASS_FLOW] * 375, []],
        [list(range(1, 376)), []],
        [[1] * 375, []],
        stale_mesh_population=376,
    )
    old_target = result["mass_global"] / 376
    assert old_target / MASS_FLOW == pytest.approx(375 / 376, rel=2.0e-14)
    assert result["target_mass_flow"] == pytest.approx(MASS_FLOW, rel=2.0e-14)


def test_musubi_patch_uses_active_denominator_and_fail_fast_invariants() -> None:
    patch = (Path(__file__).parents[1] / "patches/musubi/adaptive_flux_pressure.patch").read_text(
        encoding="utf-8"
    )
    assert "targetMassFlow = massGlobal / real(activeCountGlobal, kind=rk)" in patch
    assert "targetMassFlow = massGlobal / real(globBC%nElems_total, kind=rk)" not in patch
    assert "point population mismatch" in patch
    assert "invalid point index" in patch
    assert "non-finite mass-flow sample" in patch
    assert "inactive element in active list" in patch
    assert "sampleCountGlobal /= activeCountGlobal" in patch
