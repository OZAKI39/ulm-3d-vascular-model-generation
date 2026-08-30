from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from utils.cfd_flow import grid_convergence as gc
from utils.cfd_flow.io import FlowError


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_base_observables_are_recovered_from_archive_and_referee_evidence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    archive = tmp_path / "outputs" / "cfd_flow" / "archive"
    header = archive / "restart" / "accepted_header.lua"
    binary = archive / "restart" / "accepted.lsb"
    header.parent.mkdir(parents=True)
    header.write_text("header", encoding="utf-8")
    binary.write_bytes(b"pdf")
    _write_json(
        archive / "qc" / "accepted_steady_manifest.json",
        {
            "status": "PASS",
            "iteration": 100,
            "archived_header": str(header),
            "archived_binary": str(binary),
            "binary_sha256": "accepted-sha",
        },
    )
    referee = tmp_path / "outputs" / "cfd_flow" / "referee" / "qc" / "referee"
    residual_fields = [
        "iteration", "restart_sha256", "mean_outlet_01", "mean_outlet_02",
        "mean_outlet_03",
    ]
    _write_csv(
        referee / "corrected_residual_history.csv",
        residual_fields,
        [{
            "iteration": 100,
            "restart_sha256": "accepted-sha",
            "mean_outlet_01": 2.0,
            "mean_outlet_02": 2.0,
            "mean_outlet_03": 2.0,
        }],
    )
    flux_fields = [
        "iteration", "inlet_into_domain_kg_s", "outlet_01_outward_kg_s",
        "outlet_02_outward_kg_s", "outlet_03_outward_kg_s",
    ]
    _write_csv(
        referee / "corrected_boundary_flux_history.csv",
        flux_fields,
        [
            {"iteration": 80, "inlet_into_domain_kg_s": 10, "outlet_01_outward_kg_s": 1,
             "outlet_02_outward_kg_s": 2, "outlet_03_outward_kg_s": 2},
            {"iteration": 90, "inlet_into_domain_kg_s": 10, "outlet_01_outward_kg_s": 2,
             "outlet_02_outward_kg_s": 2, "outlet_03_outward_kg_s": 2},
            {"iteration": 100, "inlet_into_domain_kg_s": 10, "outlet_01_outward_kg_s": 3,
             "outlet_02_outward_kg_s": 2, "outlet_03_outward_kg_s": 2},
        ],
    )
    monkeypatch.setattr(gc, "ARCHIVE_RUN", "archive")
    monkeypatch.setattr(gc, "RUN_NAME", "referee")
    monkeypatch.setattr(gc, "BASE_MESH_RUN", "mesh")
    monkeypatch.setattr(gc, "ACCEPTED_ITERATION", 100)
    monkeypatch.setattr(gc, "ACCEPTED_SHA256", "accepted-sha")
    monkeypatch.setattr(gc, "EXPECTED_CELLS", 2)
    monkeypatch.setattr(gc, "AUDIT_WINDOW_ITERATIONS", 20)
    monkeypatch.setattr(gc, "sha256_file", lambda _: "accepted-sha")
    monkeypatch.setattr(
        gc,
        "parse_restart_header",
        lambda _: SimpleNamespace(
            iteration=100, n_elems=2, n_components=19, binary_path=binary
        ),
    )
    state = {
        "global_mean_velocity_m_s": 1.0,
        "global_mean_pressure_gauge_pa": 2.0,
        "maximum_physical_velocity_m_s": 3.0,
    }
    monkeypatch.setattr(gc, "_grid_pdf_state", lambda *_: (state, None, np.zeros((2, 19))))
    monkeypatch.setattr(gc, "load_mesh_contract", lambda *_, **__: object())
    monkeypatch.setattr(
        gc,
        "replay_boundary_step",
        lambda *_, **__: {"details": {"inlet": {"rho": 0.75}}},
    )

    recovered = gc._base_grid_observables(tmp_path)

    assert recovered["outlet_flow_fractions"] == pytest.approx(
        {"outlet_01": 0.2, "outlet_02": 0.2, "outlet_03": 0.2}
    )
    assert sum(recovered["outlet_flow_fractions"].values()) == pytest.approx(0.6)
    assert recovered["provenance"]["accepted_restart_sha256"] == "accepted-sha"
    assert recovered["provenance"]["long_window_sample_iterations"] == [80, 90, 100]
    assert recovered["provenance"]["fractions_not_renormalized"] is True


def test_grid_specs_use_constant_r_and_diffusive_time_scaling() -> None:
    specs = gc.grid_specs()
    assert specs["coarse"].dx_m / specs["base"].dx_m == pytest.approx(1.3)
    assert specs["base"].dx_m / specs["fine"].dx_m == pytest.approx(1.3)
    for spec in specs.values():
        assert spec.dt_s == pytest.approx(
            gc.DT_S * (spec.dx_m / gc.DX_M) ** 2, rel=1.0e-14
        )


def test_monotonic_sequence_reports_order_richardson_and_both_gcis() -> None:
    r = 1.3
    fine = 1.01
    base = 1.0 + 0.01 * r**2
    coarse = 1.0 + 0.01 * r**4
    result = gc.three_grid_scalar_analysis(coarse, base, fine, refinement_ratio=r)
    assert result["classification"] == "ASYMPTOTIC_MONOTONIC"
    assert result["observed_order_p"] == pytest.approx(2.0)
    assert result["richardson_extrapolation"] == pytest.approx(1.0)
    assert result["gci_coarse_base"] is not None
    assert result["gci_base_fine"] is not None


@pytest.mark.parametrize(
    ("values", "classification"),
    [
        ((1.1, 0.9, 1.0), "OSCILLATORY"),
        ((1.0, 1.0, 1.0), "STALLED"),
        ((3.0, 2.0, 0.0), "NON_ASYMPTOTIC"),
    ],
)
def test_invalid_three_grid_trends_do_not_fabricate_gci(
    values: tuple[float, float, float], classification: str
) -> None:
    result = gc.three_grid_scalar_analysis(*values)
    assert result["status"] == "UNAVAILABLE"
    assert result["classification"] == classification
    assert result["observed_order_p"] is None
    assert result["richardson_extrapolation"] is None
    assert result["gci_coarse_base"] is None
    assert result["gci_base_fine"] is None


def _analysis(relative: float, *, available: bool = True) -> dict[str, object]:
    return {
        "status": "AVAILABLE" if available else "UNAVAILABLE",
        "classification": "ASYMPTOTIC_MONOTONIC" if available else "OSCILLATORY",
        "base_to_fine_relative_difference": relative,
    }


def test_five_percent_is_hard_gate_and_two_percent_is_preferred_gate() -> None:
    primary = ("a", "b")
    engineering = gc.evaluate_grid_convergence_gate(
        {"a": _analysis(0.049), "b": _analysis(0.03)}, primary
    )
    assert engineering["status"] == "PASS"
    assert engineering["base_fine_primary_within_preferred_2_percent"] is False
    preferred = gc.evaluate_grid_convergence_gate(
        {"a": _analysis(0.019), "b": _analysis(0.02)}, primary
    )
    assert preferred["status"] == "PASS"
    assert preferred["base_fine_primary_within_preferred_2_percent"] is True
    too_large = gc.evaluate_grid_convergence_gate(
        {"a": _analysis(0.051), "b": _analysis(0.01)}, primary
    )
    assert too_large["status"] == "FAIL"
    assert too_large["next"] == "DESIGN ONE FINER GRID FOR CONVERGENCE CONFIRMATION"


def test_oscillatory_primary_metric_selects_geometry_review() -> None:
    gate = gc.evaluate_grid_convergence_gate(
        {"a": _analysis(0.01, available=False)}, ("a",)
    )
    assert gate["status"] == "FAIL"
    assert gate["next"] == "REVIEW VOXELIZED GEOMETRY GRID-SENSITIVITY"


def test_deferred_attach_passes_550k_without_launching_musubi(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    captured: dict[str, object] = {}

    def fake_monitor(*args: object, **kwargs: object) -> dict[str, str]:
        captured["args"] = args
        captured["kwargs"] = kwargs
        return {"status": "IN_PROGRESS"}

    monkeypatch.setattr(gc, "monitor_existing_grid_steady", fake_monitor)
    result = gc.monitor_existing_grid_steady_deferred(
        tmp_path, "fine", "/home/lzy/u3da/gc_fine_continue_test"
    )
    assert result["status"] == "IN_PROGRESS"
    assert captured["kwargs"] == {"defer_full_audit_until_iteration": 550_000}


def test_duplicate_fine_musubi_launch_is_blocked(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        gc,
        "_active_grid_musubi_processes",
        lambda: [{"runtime_root_wsl": "/home/lzy/u3da/gc_fine_continue_live"}],
    )
    with pytest.raises(FlowError, match="already active"):
        gc._assert_no_duplicate_grid_musubi("fine")


def test_grid_process_history_counts_real_records(tmp_path: Path) -> None:
    steady = tmp_path / "steady"
    _write_json(steady / "history" / "initial_process_provenance.json", {"segment": "initial"})
    _write_json(
        steady / "continuations" / "from_10" / "process_provenance.json",
        {"segment": "from_10"},
    )
    assert [item["segment"] for item in gc._grid_process_history(steady)] == [
        "initial", "from_10"
    ]


def test_solved_grid_observables_carry_referee_v2_long_window_provenance(
    tmp_path: Path,
) -> None:
    summary = {
        "status": "PASS",
        "steady_iteration": 50,
        "accepted_restart_header": "accepted.lua",
        "accepted_restart_binary": "accepted.lsb",
        "accepted_restart_sha256": "sha",
        "referee_v2": {
            "inlet_gauge_pressure_pa": 10.0,
            "pressure_drops_pa": {"outlet_01": 1.0, "outlet_02": 2.0, "outlet_03": 3.0},
            "outlet_flow_fractions": {"outlet_01": 0.1, "outlet_02": 0.2, "outlet_03": 0.3},
            "global_mean_velocity_m_s": 0.4,
            "global_mean_pressure_gauge_pa": 0.5,
            "maximum_physical_velocity_m_s": 0.6,
            "referee_revision": gc.REFEREE_REVISION_NEW,
            "long_boundary_window": {
                "sample_iterations": [30, 35, 40, 45, 50],
                "mean_port_flows_kg_s": {
                    "inlet": 1.0, "outlet_01": 0.1,
                    "outlet_02": 0.2, "outlet_03": 0.3,
                },
            },
        },
    }
    _write_json(
        tmp_path / "outputs" / "cfd_flow" / gc.GRID_RUN / "grids" / "coarse"
        / "steady" / "steady_summary.json",
        summary,
    )
    result = gc._solved_grid_observables(tmp_path, "coarse")
    provenance = result["provenance"]
    assert provenance["observable_semantics"] == (
        "REFEREE_V2_ACCEPTED_RESTART_AND_20K_LONG_WINDOW"
    )
    assert provenance["long_window_sample_iterations"] == [30, 35, 40, 45, 50]
    assert provenance["fractions_not_renormalized"] is True


def test_clean_nonsteady_segment_is_archived_as_history_not_numerical_failure(
    tmp_path: Path,
) -> None:
    steady = tmp_path / "steady"
    _write_json(
        steady / "steady_summary.json",
        {"status": "FAIL", "returncode": 0, "latest_iteration": 9431, "wall_time_s": 2.0},
    )
    gc._archive_intermediate_summary(steady)
    historical = json.loads(
        (steady / "history" / "initial_segment_summary.json").read_text(encoding="utf-8")
    )
    process = json.loads(
        (steady / "history" / "initial_process_provenance.json").read_text(encoding="utf-8")
    )
    assert historical["historical_status"] == "COMPLETED_NOT_STEADY"
    assert historical["status_was_not_a_numerical_failure"] is True
    assert process["status"] == "COMPLETED_NOT_STEADY"


def test_process_reconciliation_records_shutdown_iteration_after_gate_pass(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    steady = (
        tmp_path
        / "outputs"
        / "cfd_flow"
        / gc.GRID_RUN
        / "grids"
        / "fine"
        / "steady"
    )
    _write_json(
        steady / "steady_summary.json",
        {
            "status": "PASS",
            "steady_iteration": 554431,
            "runtime_root_wsl": "/runtime/fine",
        },
    )
    _write_json(
        steady / "history" / "initial_process_provenance.json",
        {"segment": "initial", "end_iteration": 9431, "status": "COMPLETED"},
    )
    continuation = steady / "continuations" / "from_9431"
    _write_json(
        continuation / "progress.json",
        {
            "runtime_root_wsl": "/runtime/fine",
            "latest_iteration": 556148,
            "elapsed_s": 23296.2,
        },
    )
    monkeypatch.setattr(gc, "_grid_runtime_is_running", lambda *_: False)

    result = gc.reconcile_grid_process_provenance(tmp_path, "fine")

    assert result["musubi_process_count"] == 2
    assert result["processes"][-1]["end_iteration"] == 556148
    assert result["processes"][-1]["returncode"] == 0
    summary = json.loads((steady / "steady_summary.json").read_text(encoding="utf-8"))
    assert summary["gate_pass_iteration"] == 554431
    assert summary["stop_request_iteration"] == 554431
    assert summary["actual_final_iteration"] == 556148


def test_reconciliation_replaces_stale_active_status_but_preserves_history(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run = tmp_path / "outputs" / "cfd_flow" / gc.GRID_RUN
    _write_json(
        run / "qc" / "grid_convergence_status.json",
        {"final_status": "CFD_FLOW_GRID_SEEDER_PATH_FIX_PASS"},
    )
    for label in ("coarse", "fine"):
        _write_json(
            run / "grids" / label / "steady" / "steady_summary.json",
            {"status": "PASS", "steady_iteration": 10, "selected_ranks": 4},
        )
        _write_json(run / "grids" / label / "qc" / "mesh_qc.json", {"status": "PASS"})
        _write_json(
            run / "grids" / label / "benchmark" / "benchmark_summary.json",
            {"status": "PASS"},
        )
    _write_json(
        tmp_path / "outputs" / "cfd_flow" / gc.ARCHIVE_RUN / "qc"
        / "accepted_steady_manifest.json",
        {
            "status": "PASS",
            "iteration": gc.ACCEPTED_ITERATION,
            "binary_sha256": gc.ACCEPTED_SHA256,
        },
    )
    monkeypatch.setattr(gc, "_git", lambda *_: "head")
    monkeypatch.setattr(gc, "_production_paths", lambda _: tuple())
    monkeypatch.setattr(
        gc,
        "reconcile_grid_process_provenance",
        lambda *_: {"status": "PASS", "musubi_process_count": 2, "processes": []},
    )
    reconciled = gc.reconcile_grid_convergence_evidence(tmp_path)
    assert reconciled["evidence_reconciliation"] == "PASS"
    assert reconciled["grid_convergence"] == "READY_FOR_GCI"
    assert reconciled["historical_status"]["final_status"] == (
        "CFD_FLOW_GRID_SEEDER_PATH_FIX_PASS"
    )


def test_runtime_path_is_recovered_from_unc_restart_evidence() -> None:
    value = r"\\wsl.localhost\Ubuntu\home\lzy\u3da\gc_fine_steady_x\restart\a.lsb"
    assert gc._runtime_from_restart_path(value) == "/home/lzy/u3da/gc_fine_steady_x"
