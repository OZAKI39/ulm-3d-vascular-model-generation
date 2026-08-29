from __future__ import annotations

import numpy as np
import pytest

from utils.cfd_flow.exact_link_flux import (
    INVERSE_DIRECTIONS,
    MFR_TARGET_MISMATCH,
    classify_exact_flux,
    closest_discrete_direction,
    equilibrium_pdf,
    mfr_eq_area_proxy,
    pull_fetch_pdfs,
    reconstruct_boundary,
    signed_mass_balance,
)
from utils.cfd_flow.restart_decode import D3Q19_DIRECTIONS, D3Q19_WEIGHTS


def test_mfr_eq_area_proxy_uses_boundary_element_count_times_dx_squared() -> None:
    assert mfr_eq_area_proxy(311, 2.0e-7) == pytest.approx(1.244e-11)
    assert mfr_eq_area_proxy(334, 2.0e-7) == pytest.approx(1.336e-11)


def test_axis_aligned_normal_ind_mapping() -> None:
    assert closest_discrete_direction(np.asarray((8, 0, 0))) == 3


def test_oblique_normal_ind_mapping() -> None:
    index = closest_discrete_direction(np.asarray((1, -7, -9)))
    assert D3Q19_DIRECTIONS[index].tolist() == [0, -1, -1]


def test_boundary_bitmask_uses_inverse_direction_and_weighted_normal() -> None:
    boundary_ids = np.zeros((1, 26), dtype=np.int64)
    boundary_ids[0, 4] = 5  # outward +y -> incoming -y
    result = reconstruct_boundary(
        boundary_ids, np.asarray([7]), label="inlet", boundary_id=5
    )
    assert result.cell_indices.tolist() == [7]
    assert result.raw_normals.tolist() == [[0, -4, 0]]
    assert np.flatnonzero(result.incoming_masks[0]).tolist() == [1]
    assert int(INVERSE_DIRECTIONS[4]) == 1


def test_strict_d3q19_equilibrium_recovers_density_and_velocity() -> None:
    rho = 1.02
    velocity = np.asarray((0.012, -0.006, 0.003))
    pdf = equilibrium_pdf(rho, velocity)
    assert np.sum(pdf) == pytest.approx(rho, abs=1.0e-15)
    momentum = pdf @ D3Q19_DIRECTIONS.astype(np.float64)
    assert momentum == pytest.approx(rho * velocity, abs=1.0e-15)
    assert equilibrium_pdf(1.0, np.zeros(3)) == pytest.approx(D3Q19_WEIGHTS)


def test_pull_fetch_reads_source_neighbor_and_bounceback_inverse() -> None:
    coordinates = np.asarray(((0, 0, 0), (1, 0, 0)), dtype=np.int64)
    pdf = np.arange(38, dtype=np.float64).reshape(2, 19)
    fetched = pull_fetch_pdfs(pdf, coordinates, np.asarray([1]))
    assert fetched[0, 3] == pdf[0, 3]  # c=+x fetches source at target-c
    assert (
        fetched[0, 0] == pdf[1, INVERSE_DIRECTIONS[0]]
    )  # missing +x source bounces back


def test_incoming_outgoing_replacement_sign() -> None:
    old_outgoing = 0.03
    new_incoming = 0.05
    inlet_domain_gain = new_incoming - old_outgoing
    outlet_signed_outward = old_outgoing - new_incoming
    assert inlet_domain_gain > 0.0
    assert outlet_signed_outward < 0.0


def test_signed_mass_balance_preserves_negative_outlet() -> None:
    result = signed_mass_balance(10.0, (6.0, -1.0, 5.0))
    assert result["outlet_signed_sum_kg_s"] == 10.0
    assert result["relative_error"] == 0.0


def test_negative_outlet_is_unresolved_when_signed_balance_fails() -> None:
    result = classify_exact_flux(0.2, 0.08, -1.0)
    assert result["status"] == MFR_TARGET_MISMATCH
    assert result["outlet_02_backflow_confirmed"] == "UNRESOLVED"


def test_negative_outlet_is_confirmed_only_when_balance_passes() -> None:
    result = classify_exact_flux(0.2, 0.005, -1.0)
    assert result["outlet_02_backflow_confirmed"] == "YES"
