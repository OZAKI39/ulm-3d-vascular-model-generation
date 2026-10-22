"""OSQP runs the same physical-limit and independent reference checks as V1."""

from functools import partial
from pathlib import Path
import numpy as np
import pytest
from scipy import sparse
import test_stage5a_motion as reference
from mevo_mb.bubbles.constrained_motion import ConstrainedLubrication
from mevo_mb.bubbles.osqp_motion import OSQPMotionQP

pytest.importorskip("osqp")


@pytest.fixture
def osqp_reference(monkeypatch):
    monkeypatch.setattr(
        reference, "ConstrainedLubrication", partial(ConstrainedLubrication, qp_solver="osqp")
    )
    monkeypatch.setattr(reference, "SparseActiveSetQP", OSQPMotionQP)
    return reference


@pytest.mark.parametrize(
    "name",
    [
        "test_T1_stokes_force_and_flow",
        "test_T2_common_translation",
        "test_T3_squeeze_large_dt",
        "test_T4_separating_inactive",
        "test_T5_tangential_not_locked",
        "test_T6_wall_large_dt",
        "test_T10_components_global",
        "test_T11_warm_cold_and_cache_ids",
        "test_T12_predictor_corrector_refinement",
        "test_rejection_does_not_commit_cache",
    ],
)
def test_osqp_physical_limits(osqp_reference, name):
    getattr(osqp_reference, name)()


@pytest.mark.parametrize("count", [3, 30])
def test_osqp_chain(osqp_reference, count):
    osqp_reference.test_contact_chain(count)


@pytest.mark.parametrize("seed", range(6))
def test_osqp_independent_dense_reference(osqp_reference, seed):
    osqp_reference.test_T9_random_sparse_dense(seed)


@pytest.mark.parametrize("seed", range(4))
def test_osqp_frozen_stage4l_reference(osqp_reference, seed):
    osqp_reference.test_T9_frozen_Stage4L_contact_reference(seed)


def test_osqp_observed_failure_cold_and_warm():
    folder = Path("tests/data/stage5a_failed_qp")
    d = dict(np.load(folder / "FAILED_QP_STATE.npz"))
    a, g = sparse.load_npz(folder / "A.npz"), sparse.load_npz(folder / "G.npz")
    results = []
    for warm in [False, True]:
        results.append(
            OSQPMotionQP().solve(
                a,
                d["b"],
                g,
                d["lower"],
                warm_primal=d["warm_primal"] if warm else None,
                warm_dual=d["warm_dual"] if warm else None,
                feasible=d["feasible"],
            )
        )
        assert results[-1][2]["kkt"] <= 1e-9
        attempt = results[-1][2]["attempts"][-1]
        assert attempt["status"] == "solved" or attempt["external_sparse_polish"]["accepted"]
    np.testing.assert_allclose(results[0][0], results[1][0], atol=1e-9, rtol=1e-8)
    assert abs(results[0][2]["objective"] - results[1][2]["objective"]) < 1e-9


def test_osqp_workspace_update_and_rejected_trial_reset():
    qp = OSQPMotionQP()
    a, g = sparse.eye(2, format="csc"), sparse.eye(2, format="csc")
    qp.solve(a, np.array([-2.0, 1.0]), g, np.zeros(2))
    u, _, diag = qp.solve(2 * a, np.array([-1.0, 3.0]), g, np.zeros(2))
    assert diag["workspace_reused"]
    assert diag["cold_start"]
    np.testing.assert_allclose(u, [0.0, 1.5], atol=1e-10)


def test_unknown_qp_rejected():
    with pytest.raises(ValueError, match="Unknown constrained QP"):
        ConstrainedLubrication(0.001, qp_solver="guess")
